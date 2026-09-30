from __future__ import annotations

import re
from collections import defaultdict
from functools import lru_cache
from pathlib import Path
from typing import Callable

from dna_compare.comparisons.additional import _coverage_status, _find_row, _lookups, _norm_rsid
from dna_compare.config import Settings
from dna_compare.models import DrugBlock, DrugResponse, PharmacogeneCall
from dna_compare.vcf_parser import normalize_chrom

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
SUITABILITY_RANK = {"suitable": 0, "partial": 1, "not_suitable": 2}
CONFIDENCE_RANK = {"none": 0, "low": 1, "medium": 2, "high": 3}

METABOLIZER_LABELS = {
    "UM": "ultrarapid metabolizer",
    "RM": "rapid metabolizer",
    "NM": "normal metabolizer",
    "IM": "intermediate metabolizer",
    "IM1.5": "intermediate metabolizer (activity score 1.5)",
    "IM1": "intermediate metabolizer (activity score 1)",
    "PM": "poor metabolizer",
}

DRUG_NOTES = [
    "Research overlay from a few tag SNPs, simplified from CPIC guidelines — not a clinical pharmacogenomic test. Do not start, stop, or change any medicine without your doctor or pharmacist.",
    "Suitable = standard use expected for this genotype; partially suitable = dose change, closer monitoring, or a preferred alternative; not suitable = guideline advises avoiding the drug for this genotype; unknown = the needed SNPs were not called.",
    "Genotypes are unphased: two different variants in one gene are assumed to sit on different copies of the gene.",
    "Coding alleles in the pharmacogene table are scored when the exome calls them. Alleles missing from the file are assumed absent, which is the right reading of a variant-only exome and an incomplete reading of a capture that never covered the site.",
    "Not assessable here: HLA-B*15:02, HLA-A*31:01, and HLA-B*58:01 (they need an HLA type), CYP2D6 deletions and duplications, and UGT1A1*28 (a promoter repeat).",
    "A 'suitable' row only means the typed variants were absent. Other genes, other medicines, kidney/liver function, and age still matter.",
]

DRUG_MESSAGES = {
    "Phenytoin / fosphenytoin": "HLA-B*15:02 (severe skin-reaction risk) is not typed here; clinical testing is advised before starting, especially with South or Southeast Asian ancestry.",
}


def observed_alleles(row: dict | None) -> list[str] | None:
    if row is None:
        return None
    ref = str(row.get("ref") or "").upper()
    alts = [a.strip().upper() for a in str(row.get("alt_all") or row.get("alt") or "").split(",") if a.strip()]
    choices = [ref] + alts
    tokens = re.split(r"[/|]", str(row.get("genotype") or ""))
    try:
        idx = [int(token) for token in tokens]
    except ValueError:
        return None
    if not idx or any(i < 0 or i >= len(choices) for i in idx):
        return None
    return [choices[i] for i in idx]


def effect_copies(row: dict | None, effect: str, other: str) -> tuple[float, str] | None:
    alleles = observed_alleles(row)
    if not alleles:
        return None
    allowed = {effect, other}
    if all(a in allowed for a in alleles):
        return float(alleles.count(effect)), "/".join(alleles)
    if len(effect) == 1 and len(other) == 1 and COMPLEMENT.get(effect) != other:
        flipped = [COMPLEMENT.get(a, a) for a in alleles]
        if all(a in allowed for a in flipped):
            return float(flipped.count(effect)), "/".join(flipped)
    return None


@lru_cache(maxsize=4)
def load_pgx_markers(path: Path) -> dict[str, list[dict]]:
    grouped: dict[str, list[dict]] = defaultdict(list)
    for rec in _read_tsv(path):
        gene = rec.get("gene", "").strip()
        rsid = rec.get("rsid", "").strip()
        if not gene or not rsid:
            continue
        grouped[gene].append(
            {
                "gene": gene,
                "allele": rec.get("allele", "").strip(),
                "rsid": rsid,
                "chrom": normalize_chrom(rec.get("chrom") or ""),
                "pos_hg19": _int(rec.get("pos_hg19")),
                "pos_hg38": _int(rec.get("pos_hg38")),
                "effect_allele": rec.get("effect_allele", "").strip().upper(),
                "other_allele": rec.get("other_allele", "").strip().upper(),
                "function": rec.get("function", "").strip(),
                "activity": rec.get("activity", "").strip(),
                "tier": (rec.get("tier") or "core").strip().lower() or "core",
                "source": rec.get("source", "").strip(),
            }
        )
    return dict(grouped)


@lru_cache(maxsize=4)
def load_drug_guidance(path: Path) -> dict[str, dict]:
    drugs: dict[str, dict] = {}
    for rec in _read_tsv(path):
        drug = rec.get("drug", "").strip()
        gene = rec.get("gene", "").strip()
        phenotype = rec.get("phenotype", "").strip()
        suitability = rec.get("suitability", "").strip()
        if not drug or not gene or not phenotype or suitability not in SUITABILITY_RANK:
            continue
        entry = drugs.setdefault(drug, {"drug_class": rec.get("drug_class", "").strip(), "genes": {}})
        entry["genes"].setdefault(gene, {})[phenotype] = {
            "suitability": suitability,
            "recommendation": rec.get("recommendation", "").strip(),
            "source": rec.get("source", "").strip(),
        }
    return drugs


def _read_tsv(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    rows: list[dict] = []
    header: list[str] | None = None
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            if raw.startswith("#") or not raw.strip():
                continue
            cols = raw.rstrip("\n").split("\t")
            if header is None:
                header = cols
                continue
            rows.append({key: cols[i] if i < len(cols) else "" for i, key in enumerate(header)})
    return rows


def _int(value: str | None) -> int:
    try:
        return int(value or 0)
    except ValueError:
        return 0


def _type_markers(markers: list[dict], lookups) -> tuple[dict[str, float], list[str], list[str], list[str]]:
    by_rsid, by_pos = lookups
    copies: dict[str, float] = {}
    evidence: list[str] = []
    missing_core: list[str] = []
    missing_extra: list[str] = []
    for marker in markers:
        got = effect_copies(_find_row(marker, by_rsid, by_pos), marker["effect_allele"], marker["other_allele"])
        bucket = missing_extra if marker.get("tier") == "extra" else missing_core
        if got is None:
            bucket.append(f"{marker['allele']} ({marker['rsid']})")
            continue
        n, genotype = got
        key = _norm_rsid(marker["rsid"])
        if key in copies:
            key = f"{key}#{marker['effect_allele']}"
        copies[key] = n
        evidence.append(f"{marker['rsid']} {marker['gene']} {marker['allele']}: {genotype} ({n:.0f}× {marker['effect_allele']})")
    return copies, evidence, missing_core, missing_extra


def _spec_from_markers(markers: list[dict]) -> list[tuple[str, str, float]]:
    spec: list[tuple[str, str, float]] = []
    for marker in markers:
        raw = str(marker.get("activity") or "").strip()
        if not raw:
            continue
        spec.append((_norm_rsid(marker["rsid"]), marker["allele"], float(raw)))
    return spec


def _allele_copies(copies: dict[str, float], marker: dict) -> float:
    rsid = _norm_rsid(marker["rsid"])
    if rsid in copies:
        return copies[rsid]
    return copies.get(f"{rsid}#{marker['effect_allele']}", 0.0)


def _stars(copies: dict[str, float], spec: list[tuple[str, str, float]]) -> list[tuple[str, float]]:
    """Variant star alleles (worst function first) implied by unphased copy counts, capped at two."""
    found: list[tuple[str, float]] = []
    for rsid, star, value in sorted(spec, key=lambda item: item[2]):
        found.extend([(star, value)] * int(round(copies.get(rsid, 0.0))))
    return found[:2]


def _activity_call(
    copies: dict[str, float],
    spec: list[tuple[str, str, float]],
    *,
    normal: tuple[str, float] = ("*1", 1.0),
) -> tuple[float, str]:
    alleles = _stars(copies, spec)
    while len(alleles) < 2:
        alleles.insert(0, normal)
    score = sum(value for _star, value in alleles)
    shown = sorted(alleles, key=lambda item: -item[1])
    return score, "/".join(star for star, _value in shown)


def _call_cyp2c19(copies, ctx):
    adjusted = dict(copies)
    # *2 haplotypes carry rs12769205, so that SNP is *35 only when *2 is absent.
    n2 = adjusted.get("rs4244285", 0.0)
    adjusted["rs12769205"] = max(0.0, adjusted.get("rs12769205", 0.0) - n2)
    spec = []
    for marker in ctx["markers"]:
        rsid = _norm_rsid(marker["rsid"])
        if marker.get("function") == "no":
            spec.append((rsid, marker["allele"], 0.0))
        elif marker.get("function") == "increased":
            spec.append((rsid, marker["allele"], 2.0))
    alleles = _stars(adjusted, spec)
    while len(alleles) < 2:
        alleles.insert(0, ("*1", 1))
    n_no = sum(1 for _s, v in alleles if v == 0)
    n_inc = sum(1 for _s, v in alleles if v == 2)
    if n_no == 2:
        pheno = "PM"
    elif n_no == 1:
        pheno = "IM"
    elif n_inc == 2:
        pheno = "UM"
    elif n_inc == 1:
        pheno = "RM"
    else:
        pheno = "NM"
    return pheno, "/".join(s for s, _v in alleles), METABOLIZER_LABELS[pheno], ""


def _call_cyp2c9(copies, ctx):
    score, diplo = _activity_call(copies, _spec_from_markers(ctx["markers"]))
    pheno = "NM" if score >= 2 else "IM1.5" if score >= 1.5 else "IM1" if score >= 1 else "PM"
    return pheno, diplo, METABOLIZER_LABELS[pheno], ""


def _call_cyp2d6(copies, ctx):
    adjusted = dict(copies)
    if "rs1065852" in adjusted:
        adjusted["rs1065852"] = max(0.0, adjusted["rs1065852"] - adjusted.get("rs3892097", 0.0))
    score, diplo = _activity_call(adjusted, _spec_from_markers(ctx["markers"]))
    pheno = "PM" if score == 0 else "IM" if score < 1.25 else "NM"
    message = (
        "Deletions (*5) and gene duplications are invisible to array SNPs, so ultrarapid metabolizers cannot be "
        "identified and a normal call may hide a *5 deletion. *10 is read from rs1065852 after removing copies explained by *4."
    )
    return pheno, diplo, f"{METABOLIZER_LABELS[pheno]} (activity score {score:g})", message


def _call_dpyd(copies, ctx):
    score, diplo = _activity_call(copies, _spec_from_markers(ctx["markers"]), normal=("normal", 1.0))
    pheno = "NM" if score >= 2 else "IM" if score >= 1 else "PM"
    return pheno, diplo, f"{METABOLIZER_LABELS[pheno]} (activity score {score:g})", ""


def _call_tpmt(copies, ctx):
    c719 = int(round(copies.get("rs1142345", 0.0)))
    c460 = int(round(copies.get("rs1800460", 0.0)))
    both = min(c719, c460)
    variants = ["*3A"] * both + ["*3C"] * (c719 - both) + ["*3B"] * (c460 - both)
    for marker in ctx["markers"]:
        rsid = _norm_rsid(marker["rsid"])
        if rsid in {"rs1142345", "rs1800460"} or marker.get("function") != "no":
            continue
        variants.extend([marker["allele"]] * int(round(_allele_copies(copies, marker))))
    alleles = (["*1", "*1"] + variants[:2])[-2:]
    n = min(2, len(variants))
    pheno = "PM" if n >= 2 else "IM" if n == 1 else "NM"
    message = ""
    if c719 == 1 and c460 == 1:
        message = "Heterozygous at both c.460 and c.719 is read as *1/*3A (usual); the rare *3B/*3C (poor metabolizer) needs phased or clinical testing to exclude."
    return pheno, "/".join(alleles), METABOLIZER_LABELS[pheno], message


def _single_star(rsid: str, star: str, labels: dict[str, str] | None = None):
    def call(copies, _ctx):
        n = int(round(copies.get(rsid, 0.0)))
        pheno = {0: "NM", 1: "IM"}.get(n, "PM")
        diplo = {0: "*1/*1", 1: f"*1/{star}"}.get(n, f"{star}/{star}")
        return pheno, diplo, (labels or METABOLIZER_LABELS)[pheno], ""

    return call


def _call_vkorc1(copies, _ctx):
    n = int(round(copies.get("rs9923231", 0.0)))
    pheno = {0: "GG", 1: "GA"}.get(n, "AA")
    label = {
        "GG": "typical warfarin sensitivity",
        "GA": "moderately increased warfarin sensitivity",
        "AA": "high warfarin sensitivity",
    }[pheno]
    return pheno, f"-1639 {pheno[0]}/{pheno[1]}", label, ""


def _call_slco1b1(copies, _ctx):
    n = int(round(copies.get("rs4149056", 0.0)))
    pheno = {0: "normal", 1: "decreased"}.get(n, "poor")
    diplo = {0: "c.521 T/T", 1: "c.521 T/C"}.get(n, "c.521 C/C")
    return pheno, diplo, f"{pheno} transporter function", ""


def _call_hla_b5701(copies, _ctx):
    n = int(round(copies.get("rs2395029", 0.0)))
    pheno = "positive" if n >= 1 else "negative"
    message = "rs2395029 (HCP5) tags HLA-B*57:01 well in European ancestry but less reliably in others; it is not an HLA type."
    label = "HLA-B*57:01 tag present" if n else "HLA-B*57:01 tag absent"
    return pheno, f"rs2395029 {n}× G", label, message


def _call_cyp4f2(copies, ctx):
    n = int(round(sum(_allele_copies(copies, marker) for marker in ctx["markers"])))
    if n >= 2:
        return "hom", "*3/*3", "CYP4F2*3 homozygous — higher warfarin dose expected", ""
    if n == 1:
        return "carrier", "*1/*3", "one CYP4F2*3 allele — small warfarin dose increase", ""
    return "typical", "*1/*1", "no CYP4F2*3 allele", ""


def _call_ugt1a1(copies, ctx):
    score, diplo = _activity_call(copies, _spec_from_markers(ctx["markers"]))
    pheno = "NM" if score >= 2 else "IM" if score >= 1.5 else "PM"
    return pheno, diplo, METABOLIZER_LABELS[pheno], "UGT1A1*28, the TA-repeat promoter allele, is not read from an exome."


def _call_nat2(copies, ctx):
    n = min(2, int(round(sum(_allele_copies(copies, marker) for marker in ctx["markers"]))))
    if n >= 2:
        return "slow", "slow/slow", "slow acetylator", ""
    if n == 1:
        return "intermediate", "rapid/slow", "intermediate acetylator", ""
    return "rapid", "rapid/rapid", "rapid acetylator", ""


def _call_bche(copies, ctx):
    n = int(round(sum(_allele_copies(copies, marker) for marker in ctx["markers"])))
    message = "Only the atypical BCHE allele is typed."
    if n >= 2:
        return "PM", "atypical/atypical", "atypical cholinesterase homozygous", message
    if n == 1:
        return "IM", "*1/atypical", "atypical cholinesterase carrier", message
    return "NM", "*1/*1", "no atypical BCHE allele", message


def _call_mh(copies, ctx):
    hits = []
    for marker in ctx["markers"]:
        n = int(round(_allele_copies(copies, marker)))
        if n:
            hits.append(f"{marker['allele']}×{n}")
    if hits:
        shown = ", ".join(hits[:4]) + (f" (+{len(hits) - 4} more)" if len(hits) > 4 else "")
        return "positive", shown, "malignant-hyperthermia allele present", ""
    return (
        "negative",
        "none",
        "no curated malignant-hyperthermia allele at typed sites",
        "A negative call covers this curated list only.",
    )


def _call_g6pd(copies, ctx):
    n = int(round(sum(copies.values())))
    message = "Curated G6PD deficiency alleles only, so a normal call does not rule out deficiency."
    if ctx["male"]:
        pheno = "deficient" if n >= 1 else "normal"
        diplo = "hemizygous variant" if n >= 1 else "hemizygous reference"
    else:
        pheno = "deficient" if n >= 2 else "variable" if n == 1 else "normal"
        diplo = {0: "reference/reference", 1: "heterozygous"}.get(n, "two variant copies")
        message += " No Y-chromosome calls, so scored as female; in a male any variant copy means deficient."
    return pheno, diplo, f"G6PD {pheno}", message


def _call_nudt15(copies, ctx):
    score, diplo = _activity_call(copies, _spec_from_markers(ctx["markers"]))
    pheno = "NM" if score >= 2 else "IM" if score >= 1 else "PM"
    return pheno, diplo, METABOLIZER_LABELS[pheno], ""


GENE_CALLERS: dict[str, Callable] = {
    "CYP2C19": _call_cyp2c19,
    "CYP2C9": _call_cyp2c9,
    "VKORC1": _call_vkorc1,
    "SLCO1B1": _call_slco1b1,
    "TPMT": _call_tpmt,
    "NUDT15": _call_nudt15,
    "DPYD": _call_dpyd,
    "CYP3A5": _single_star(
        "rs776746",
        "*3",
        {
            "NM": "normal metabolizer (CYP3A5 expresser)",
            "IM": "intermediate metabolizer (CYP3A5 expresser)",
            "PM": "poor metabolizer (CYP3A5 non-expresser)",
        },
    ),
    "CYP2D6": _call_cyp2d6,
    "HLA-B": _call_hla_b5701,
    "G6PD": _call_g6pd,
    "CYP4F2": _call_cyp4f2,
    "UGT1A1": _call_ugt1a1,
    "NAT2": _call_nat2,
    "BCHE": _call_bche,
    "RYR1": _call_mh,
    "CACNA1S": _call_mh,
}

LOW_CONFIDENCE_GENES = {"CYP2D6", "HLA-B"}


def call_gene(gene: str, markers: list[dict], lookups, ctx: dict) -> PharmacogeneCall:
    copies, evidence, missing_core, missing_extra = _type_markers(markers, lookups)
    core = [marker for marker in markers if marker.get("tier") != "extra"]
    denom = core or markers
    n_found = sum(1 for marker in denom if _norm_rsid(marker["rsid"]) in copies or f"{_norm_rsid(marker['rsid'])}#{marker['effect_allele']}" in copies)
    caller = GENE_CALLERS.get(gene)
    extra_only = bool(markers) and all(marker.get("tier") == "extra" for marker in markers)
    if caller is None or (not copies and not extra_only):
        missing = missing_core or missing_extra
        return PharmacogeneCall(
            gene=gene,
            phenotype=None,
            label=f"{gene} not called on this file",
            n_markers=len(denom),
            message=f"Needed: {', '.join(missing)}." if missing else "",
        )
    if not copies and extra_only:
        return PharmacogeneCall(
            gene=gene,
            phenotype=None,
            label=f"{gene} not in this file",
            n_markers=len(markers),
            coverage_status="optional",
            message=f"{gene} variants were not in this file.",
        )
    local = dict(ctx)
    local["markers"] = markers
    phenotype, diplotype, label, message = caller(copies, local)
    notes = [message] if message else []
    if missing_core:
        notes.insert(0, f"Not typed on this file: {', '.join(missing_core)} — those alleles were assumed absent.")
    if missing_extra:
        notes.append("Other curated alleles were not in this file and were assumed absent.")
    if gene in LOW_CONFIDENCE_GENES:
        confidence = "low"
    elif missing_core:
        confidence = "medium" if n_found > len(missing_core) else "low"
    else:
        confidence = "high"
    return PharmacogeneCall(
        gene=gene,
        phenotype=phenotype,
        diplotype=diplotype,
        label=f"{gene} {diplotype} — {label}",
        evidence=evidence,
        confidence=confidence,
        n_snps=n_found,
        n_markers=len(denom),
        coverage_status=_coverage_status(n_found, len(denom)),
        message=" ".join(notes),
    )


def _unique(items) -> list[str]:
    seen: list[str] = []
    for item in items:
        if item and item not in seen:
            seen.append(item)
    return seen


def evaluate_drug(drug: str, entry: dict, calls: dict[str, PharmacogeneCall]) -> DrugResponse:
    genes = list(entry["genes"])
    verdicts: list[tuple[str, dict]] = []
    unknown: list[str] = []
    for gene in genes:
        call = calls.get(gene)
        if call and call.coverage_status == "optional":
            continue
        rule = entry["genes"][gene].get(call.phenotype) if call and call.phenotype else None
        if rule is None:
            unknown.append(gene)
        else:
            verdicts.append((gene, rule))

    gene_calls = [calls[g] for g in genes if g in calls and calls[g].coverage_status != "optional"]
    messages = [c.message for c in gene_calls] + [DRUG_MESSAGES.get(drug, "")]
    if verdicts:
        worst = max(SUITABILITY_RANK[rule["suitability"]] for _g, rule in verdicts)
        suitability = next(k for k, v in SUITABILITY_RANK.items() if v == worst)
        if unknown and suitability == "suitable":
            suitability = "unknown"
        flagged = [rule["recommendation"] for _g, rule in verdicts if SUITABILITY_RANK[rule["suitability"]] == worst]
        recommendation = " ".join(_unique(flagged))
        if unknown:
            messages.insert(0, f"{', '.join(unknown)} not called; verdict uses {', '.join(g for g, _r in verdicts)} only.")
    else:
        suitability = "unknown"
        recommendation = f"Needed SNPs for {', '.join(genes)} were not called on this file."

    confidences = [c.confidence for c in gene_calls] or ["none"]
    n_snps = sum(c.n_snps for c in gene_calls)
    n_markers = sum(c.n_markers for c in gene_calls)
    sources = _unique(
        [rule["source"] for _g, rule in verdicts]
        or [r["source"] for g in genes for r in entry["genes"][g].values()][:1]
    )
    return DrugResponse(
        drug=drug,
        drug_class=entry["drug_class"],
        genes=", ".join(genes),
        suitability=suitability,
        recommendation=recommendation,
        genotype="; ".join(c.label for c in gene_calls),
        evidence="; ".join(e for c in gene_calls for e in c.evidence),
        source="; ".join(sources),
        available=bool(verdicts),
        confidence=min(confidences, key=lambda c: CONFIDENCE_RANK.get(c, 0)),
        n_snps=n_snps,
        n_markers=n_markers,
        coverage_status=_coverage_status(n_snps, n_markers),
        message=" ".join(_unique(messages)),
    )


def compare_drugs(query_index: dict[tuple[str, int], dict], *, settings: Settings) -> DrugBlock:
    markers = load_pgx_markers(Path(settings.pgx_markers))
    guidance = load_drug_guidance(Path(settings.drug_guidance))
    notes = list(DRUG_NOTES)
    if not markers or not guidance:
        missing = [str(p) for p in (settings.pgx_markers, settings.drug_guidance) if not Path(p).is_file()]
        return DrugBlock(kind="drugs", available=False, notes=notes + [f"Missing drug tables: {', '.join(missing)}."])

    lookups = _lookups(query_index)
    ctx = {"male": any(chrom == "Y" for chrom, _pos in query_index)}
    calls = {gene: call_gene(gene, gene_markers, lookups, ctx) for gene, gene_markers in markers.items()}
    estimates = [evaluate_drug(drug, entry, calls) for drug, entry in guidance.items()]

    counts = {key: sum(1 for e in estimates if e.suitability == key) for key in ("suitable", "partial", "not_suitable", "unknown")}
    notes.append(
        f"{len(estimates)} drugs: {counts['suitable']} suitable, {counts['partial']} partially suitable, "
        f"{counts['not_suitable']} not suitable, {counts['unknown']} unknown on this file."
    )
    return DrugBlock(
        kind="drugs",
        available=any(e.available for e in estimates),
        estimates=estimates,
        notes=notes,
        genes=list(calls.values()),
    )
