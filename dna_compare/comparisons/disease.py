from __future__ import annotations

import gzip
import math
from functools import lru_cache
from pathlib import Path

from dna_compare._native import native_available, score_pgs_batch_native
from dna_compare.config import (
    DISEASE_COVERAGE_PERCENTILE_PCT,
    DISEASE_COVERAGE_READY_PCT,
    DISEASE_ELEVATED_PERCENTILE,
    PGS_CATALOG_SCORES,
    PGS_REFERENCE_EUR,
    Settings,
)
from dna_compare.genotype_index import build_rsid_index
from dna_compare.liftover import ensure_hg38_to_hg19_chain
from dna_compare.models import ComparisonBlock, DiseaseEstimate
from dna_compare.vcf_parser import normalize_chrom

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}

DISEASE_NOTES = [
    "Published polygenic scores from the PGS Catalog. A research overlay — not a diagnosis, not medical advice, and not a lifetime-risk percentage.",
    "Common-disease rows are published SNP scores: schizophrenia, depression, bipolar disorder, ADHD, anxiety, fluid intelligence, educational attainment, coronary artery disease, type 2 diabetes, stroke, Alzheimer's, asthma, plus breast / prostate / colorectal cancer. There is no single “all disease” or “all cancer” number.",
    "Percentiles and row highlighting require ≥50% SNP overlap. Highlighting is reserved for disease traits at or above the 90th percentile of a European reference (published cohort stats or GWAS effect-allele frequencies). Low coverage is shown as a badge, not row color.",
    "Fluid intelligence and educational attainment PGS are research polygenic sums — not an IQ test, not a g-factor estimate, and not framed as disease risk.",
    "These GWAS were mostly European. A South Asian kit is scored on overlapping SNPs only; coverage will be incomplete and uncalibrated scores are not percentiles.",
    "SNP arrays miss most rare pathogenic variants (BRCA1/2, Lynch genes, high-impact SCZ CNVs). A low score does not rule those out.",
]


def _normal_cdf(z: float) -> float:
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))


@lru_cache(maxsize=1)
def _load_pgs_reference(path: str) -> dict[str, dict[str, float | str]]:
    refs: dict[str, dict[str, float | str]] = {}
    ref_path = Path(path)
    if not ref_path.exists():
        return refs
    with ref_path.open("rt") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 3:
                continue
            pgs_id = cols[0].strip()
            try:
                mean = float(cols[1])
                sd = float(cols[2])
            except ValueError:
                continue
            if sd <= 0:
                continue
            refs[pgs_id] = {
                "mean": mean,
                "sd": sd,
                "ref_pop": cols[3].strip() if len(cols) > 3 else "European",
                "source": "published",
            }
    return refs


@lru_cache(maxsize=64)
def _reference_from_scoring_file(path: str) -> tuple[float, float, int, int] | None:
    mean = 0.0
    var = 0.0
    n_af = 0
    n_total = 0
    ref_path = Path(path)
    if not ref_path.exists():
        return None
    with _open_text(ref_path) as handle:
        header: list[str] | None = None
        for raw in handle:
            if raw.startswith("#"):
                continue
            cols = raw.rstrip("\n").split("\t")
            if header is None:
                header = cols
                continue
            if not header or len(cols) < len(header):
                continue
            rec = {key: cols[i] if i < len(cols) else "" for i, key in enumerate(header)}
            try:
                weight = float(rec.get("effect_weight") or "")
            except ValueError:
                continue
            n_total += 1
            try:
                freq = float(rec.get("allelefrequency_effect") or "")
            except ValueError:
                continue
            if not 0.0 <= freq <= 1.0:
                continue
            dosage = 2.0 * freq
            mean += weight * dosage
            var += (weight**2) * 2.0 * freq * (1.0 - freq)
            n_af += 1
    if n_total == 0 or n_af == 0:
        return None
    if n_af / n_total < 0.5:
        return None
    sd = math.sqrt(var)
    if sd <= 0:
        return None
    return mean, sd, n_af, n_total


def _resolve_pgs_reference(
    pgs_id: str,
    scoring_path: Path | None,
    tsv_refs: dict[str, dict[str, float | str]],
) -> dict[str, float | str] | None:
    if pgs_id in tsv_refs:
        return dict(tsv_refs[pgs_id])
    if scoring_path is None:
        return None
    computed = _reference_from_scoring_file(str(scoring_path))
    if computed is None:
        return None
    mean, sd, _n_af, _n_total = computed
    return {
        "mean": mean,
        "sd": sd,
        "ref_pop": "European (GWAS effect-allele AF)",
        "source": "gwas_af",
    }


def _coverage_status(*, used: int, score: float | None, coverage_pct: float) -> str:
    if used == 0 or score is None:
        return "missing"
    if coverage_pct < 10:
        return "low"
    if coverage_pct < DISEASE_COVERAGE_READY_PCT:
        return "partial"
    return "ready"


def _percentile_from_reference(
    score: float | None,
    *,
    coverage_pct: float,
    ref: dict[str, float | str] | None,
) -> float | None:
    if score is None or coverage_pct < DISEASE_COVERAGE_PERCENTILE_PCT or ref is None:
        return None
    mean = float(ref["mean"])
    sd = float(ref["sd"])
    z = (score - mean) / sd
    return round(100.0 * _normal_cdf(z), 1)


def _relative_level(
    *,
    trait_category: str,
    coverage_pct: float,
    percentile: float | None,
) -> str:
    if trait_category == "cognitive":
        return ""
    if coverage_pct < DISEASE_COVERAGE_PERCENTILE_PCT:
        return ""
    if percentile is None:
        return "uncalibrated"
    if percentile >= DISEASE_ELEVATED_PERCENTILE:
        return "elevated"
    return "typical"


def _open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open("rt")


def _is_snp_allele(allele: str) -> bool:
    return len(allele) == 1 and allele in COMPLEMENT


def _match_dosage(row: dict, effect: str, other: str) -> float | None:
    dosage = row.get("dosage_alt")
    if dosage is None:
        return None
    ref = str(row.get("ref") or "").upper()
    alt = str(row.get("alt") or "").upper()
    effect = effect.upper()
    other = (other or "").upper()
    if not _is_snp_allele(effect) or not _is_snp_allele(alt) or not _is_snp_allele(ref):
        return None
    if alt == effect and (not other or ref == other):
        return float(dosage)
    if ref == effect and (not other or alt == other):
        return 2.0 - float(dosage)
    if other and COMPLEMENT.get(alt) == effect and COMPLEMENT.get(ref) == other:
        return float(dosage)
    if other and COMPLEMENT.get(ref) == effect and COMPLEMENT.get(alt) == other:
        return 2.0 - float(dosage)
    return None


def _read_pgs_variant(rec: dict) -> dict | None:
    rsid = (rec.get("rsID") or rec.get("hm_rsID") or "").strip()
    chrom = normalize_chrom(rec.get("hm_chr") or rec.get("chr_name") or "")
    pos_raw = rec.get("hm_pos") or rec.get("chr_position") or ""
    try:
        pos = int(pos_raw)
    except ValueError:
        pos = 0
    effect = (rec.get("effect_allele") or "").strip().upper()
    other = (rec.get("other_allele") or "").strip().upper()
    try:
        weight = float(rec.get("effect_weight") or "")
    except ValueError:
        return None
    if not _is_snp_allele(effect):
        return None
    if other and not _is_snp_allele(other):
        return None
    if not rsid and not (chrom and pos):
        return None
    return {
        "rsid": rsid if rsid not in {".", ""} else "",
        "chrom": chrom,
        "pos": pos,
        "effect": effect,
        "other": other,
        "weight": weight,
    }


@lru_cache(maxsize=32)
def _pgs_lookup_cached(path: str) -> tuple[dict[str, dict], dict[tuple[str, int], dict], int]:
    by_rsid: dict[str, dict] = {}
    by_pos: dict[tuple[str, int], dict] = {}
    n_total = 0
    with _open_text(Path(path)) as handle:
        header: list[str] | None = None
        for raw in handle:
            if raw.startswith("#"):
                continue
            cols = raw.rstrip("\n").split("\t")
            if header is None:
                header = cols
                continue
            if not header or len(cols) < len(header):
                continue
            rec = {key: cols[i] if i < len(cols) else "" for i, key in enumerate(header)}
            var = _read_pgs_variant(rec)
            if var is None:
                continue
            n_total += 1
            rsid = str(var.get("rsid") or "").lower()
            if rsid:
                by_rsid[rsid] = var
            chrom = var.get("chrom")
            pos = var.get("pos")
            if chrom and pos:
                by_pos[(chrom, pos)] = var
    return by_rsid, by_pos, n_total


def parse_pgs_variants(path: Path | str) -> list[dict]:
    by_rsid, by_pos, _n_total = _pgs_lookup_cached(str(path))
    variants: list[dict] = []
    seen: set[int] = set()
    for var in by_rsid.values():
        key = id(var)
        if key in seen:
            continue
        seen.add(key)
        variants.append(var)
    for var in by_pos.values():
        key = id(var)
        if key in seen:
            continue
        seen.add(key)
        variants.append(var)
    return variants


def overlap_pgs_variants_for_query(
    path: Path | str,
    *,
    by_rsid: dict[str, dict],
    query_index: dict[tuple[str, int], dict],
) -> tuple[list[dict], int]:
    pgs_by_rsid, pgs_by_pos, n_total = _pgs_lookup_cached(str(path))
    overlap: list[dict] = []
    seen: set[int] = set()
    for rsid in by_rsid:
        var = pgs_by_rsid.get(rsid.lower())
        if var is None:
            continue
        key = id(var)
        if key in seen:
            continue
        seen.add(key)
        overlap.append(var)
    for chrom_pos in query_index:
        var = pgs_by_pos.get(chrom_pos)
        if var is None:
            continue
        key = id(var)
        if key in seen:
            continue
        seen.add(key)
        overlap.append(var)
    return overlap, n_total


def ensure_pgs_file(spec: dict, settings: Settings) -> Path | None:
    dest = Path(settings.disease_dir) / str(spec["filename"])
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    if not spec.get("auto_download") or not settings.auto_download_pgs:
        return dest if dest.exists() else None
    dest.parent.mkdir(parents=True, exist_ok=True)
    return ensure_hg38_to_hg19_chain(dest, url=str(spec["url"]), download=True)


def score_pgs(
    index: dict[tuple[str, int], dict],
    variants: list[dict] | tuple[dict, ...],
    *,
    by_rsid: dict[str, dict] | None = None,
    n_score: int | None = None,
) -> tuple[float, int, int]:
    if by_rsid is None:
        by_rsid = build_rsid_index(index)
    by_pos = index
    total = 0.0
    used = 0
    for var in variants:
        row = by_rsid.get(str(var.get("rsid") or "").lower()) if var.get("rsid") else None
        if row is None and var.get("chrom") and var.get("pos"):
            row = by_pos.get((var["chrom"], var["pos"]))
        if row is None:
            continue
        dosage = _match_dosage(row, var["effect"], var["other"])
        if dosage is None:
            continue
        total += dosage * var["weight"]
        used += 1
    total_variants = n_score if n_score is not None else len(variants)
    return total, used, total_variants


def compare_disease(
    query_index: dict[tuple[str, int], dict],
    *,
    settings: Settings,
    by_rsid: dict[str, dict] | None = None,
) -> ComparisonBlock:
    notes = list(DISEASE_NOTES)
    estimates: list[DiseaseEstimate] = []
    Path(settings.disease_dir).mkdir(parents=True, exist_ok=True)
    references = _load_pgs_reference(str(PGS_REFERENCE_EUR))
    if by_rsid is None:
        by_rsid = build_rsid_index(query_index)

    specs_with_paths: list[tuple[dict, Path | None]] = []
    for spec in PGS_CATALOG_SCORES:
        path = ensure_pgs_file(spec, settings)
        specs_with_paths.append((spec, path))

    variant_batches: list[list[dict]] = []
    score_totals: list[int] = []
    ready_specs: list[tuple[dict, Path]] = []
    for spec, path in specs_with_paths:
        if path is None or not path.exists():
            estimates.append(
                DiseaseEstimate(
                    trait=str(spec["trait"]),
                    pgs_id=str(spec["pgs_id"]),
                    score=None,
                    n_snps=0,
                    n_score=int(spec.get("n_variants") or 0),
                    coverage_pct=None,
                    coverage_status="missing",
                    trait_category=str(spec.get("trait_category") or "disease"),
                    citation=str(spec["citation"]),
                    available=False,
                    message=f"Scoring file missing: {spec['filename']}. Run with network to download, or place it under data/references/disease/.",
                )
            )
            continue
        ready_specs.append((spec, path))
        overlap, n_total = overlap_pgs_variants_for_query(
            path,
            by_rsid=by_rsid,
            query_index=query_index,
        )
        variant_batches.append(overlap)
        score_totals.append(int(spec.get("n_variants") or 0) or n_total)

    batch_scores: list[tuple[float, int, int]] | None = None
    if native_available() and variant_batches:
        batch_scores = score_pgs_batch_native(by_rsid, query_index, variant_batches)

    for i, (spec, path) in enumerate(ready_specs):
        n_score = score_totals[i]
        if batch_scores is not None:
            raw, used, _ = batch_scores[i]
        else:
            raw, used, _ = score_pgs(
                query_index,
                variant_batches[i],
                by_rsid=by_rsid,
                n_score=n_score,
            )
        coverage = (100.0 * used / n_score) if n_score else 0.0
        score = round(raw, 4) if used else None
        coverage_status = _coverage_status(used=used, score=score, coverage_pct=coverage)
        trait_category = str(spec.get("trait_category") or "disease")
        pgs_id = str(spec["pgs_id"])
        ref = _resolve_pgs_reference(pgs_id, path, references)
        percentile = _percentile_from_reference(score, coverage_pct=coverage, ref=ref)
        relative_level = _relative_level(
            trait_category=trait_category,
            coverage_pct=coverage,
            percentile=percentile,
        )
        message = None
        if used == 0:
            message = "No overlapping SNPs in this VCF."
        elif coverage < 10:
            message = "Very low overlap; treat the raw score as incomplete."
        elif coverage < DISEASE_COVERAGE_PERCENTILE_PCT and trait_category == "disease":
            message = "Overlap below 50%; percentile not shown until more score SNPs are present."
        elif relative_level == "uncalibrated" and trait_category == "disease":
            message = "Raw score only — no European reference distribution for this PGS yet."
        elif ref is not None and ref.get("source") == "gwas_af" and percentile is not None:
            message = "Percentile vs European GWAS effect-allele frequencies (research context only)."
        estimates.append(
            DiseaseEstimate(
                trait=str(spec["trait"]),
                pgs_id=str(spec["pgs_id"]),
                score=score,
                n_snps=used,
                n_score=n_score,
                coverage_pct=round(coverage, 1),
                coverage_status=coverage_status,
                trait_category=trait_category,
                percentile=percentile,
                relative_level=relative_level,
                citation=str(spec["citation"]),
                available=used > 0,
                message=message,
            )
        )
    available = any(item.available and item.score is not None for item in estimates)
    if not available:
        notes.append("Need published PGS Catalog files under data/references/disease/ and overlapping SNPs in the VCF.")
    return ComparisonBlock(kind="disease", available=available, estimates=estimates, notes=notes)
