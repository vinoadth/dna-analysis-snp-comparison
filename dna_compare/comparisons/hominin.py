from __future__ import annotations

from pathlib import Path

import numpy as np

from dna_compare.align import PanelAlignment, align_query_to_panel, query_allele1_freq
from dna_compare.config import ARCHAIC_AADR_SAMPLES, HOMININ_DOWNLOADS, Settings
from dna_compare.models import ComparisonBlock, HomininEstimate, ReferenceStatus
from dna_compare.vcf_parser import normalize_chrom


def _missing_hominin_vcfs(settings: Settings) -> list[ReferenceStatus]:
    missing = []
    for key, spec in HOMININ_DOWNLOADS.items():
        filename = spec["filename"]
        path = settings.hominin_dir / filename
        missing.append(
            ReferenceStatus(
                key=key,
                filename=filename,
                present=path.exists(),
                path=str(path) if path.exists() else None,
                note=spec.get("label"),
            )
        )
    return missing


def _f4(p_a: np.ndarray, p_b: np.ndarray, p_c: np.ndarray, p_d: np.ndarray) -> tuple[float, int]:
    mask = np.isfinite(p_a) & np.isfinite(p_b) & np.isfinite(p_c) & np.isfinite(p_d)
    if not np.any(mask):
        return float("nan"), 0
    value = float(np.mean((p_a[mask] - p_b[mask]) * (p_c[mask] - p_d[mask])))
    return value, int(mask.sum())


def _sharing(query: np.ndarray, archaic: np.ndarray) -> tuple[float | None, int]:
    mask = np.isfinite(query) & np.isfinite(archaic)
    if mask.sum() < 50:
        return None, int(mask.sum())
    # Probability both carry the same allele1 frequency product proxy via IBS on dosage.
    ibs = 1.0 - np.mean(np.abs(query[mask] - archaic[mask]))
    return float(100.0 * ibs), int(mask.sum())


def compare_hominin(
    query_index: dict[tuple[str, int], dict],
    *,
    panel,
    settings: Settings,
    alignment: PanelAlignment | None = None,
) -> ComparisonBlock:
    missing = _missing_hominin_vcfs(settings)
    notes = [
        "Optional EVA high-coverage VCFs are not required when AADR HO archaic samples are present.",
        "f4-ratio values are statistical ancestry estimates, not exact genealogical percentages.",
    ]
    present_optional = [item.filename for item in missing if item.present]
    if present_optional:
        notes.append("Using extra VCF/TSV files: " + ", ".join(present_optional))

    extra_estimates = _estimates_from_optional_files(query_index, settings)

    if panel is None or not panel.available:
        if extra_estimates:
            return ComparisonBlock(
                kind="hominin",
                available=True,
                estimates=extra_estimates,
                missing_files=missing,
                notes=notes,
            )
        return ComparisonBlock(
            kind="hominin",
            available=False,
            missing_files=missing,
            notes=notes + ["AADR panel missing; download the suggested hominin VCFs to enable comparison."],
        )

    if alignment is None:
        alignment = align_query_to_panel(panel, query_index)
    query = alignment.query_freqs if alignment is not None else np.array([], dtype=np.float32)
    n_overlap = alignment.n_overlap if alignment is not None else 0

    if n_overlap < 50:
        block = ComparisonBlock(
            kind="hominin",
            available=bool(extra_estimates),
            estimates=extra_estimates,
            missing_files=missing,
            notes=notes + [f"Query overlaps only {n_overlap} HO SNPs."],
        )
        return block

    freqs = {}
    for label, sample_id in ARCHAIC_AADR_SAMPLES.items():
        try:
            freqs[label] = panel.dosage_allele1(sample_id) / 2.0
        except KeyError:
            notes.append(f"Missing AADR sample {sample_id} ({label})")

    mbuti = panel.population_allele1_freq("Mbuti")
    yoruba = panel.population_allele1_freq("Yoruba")

    estimates: list[HomininEstimate] = []
    for label in ("Altai_Neanderthal", "Vindija_Neanderthal", "Chagyrskaya_Neanderthal", "Denisova"):
        if label not in freqs:
            continue
        pct, n = _sharing(query, freqs[label])
        estimates.append(
            HomininEstimate(
                label=f"{label} allele sharing",
                method="ibs_allele1",
                percent=None if pct is None else round(pct, 3),
                n_snps=n,
                available=pct is not None,
                message="IBS-like sharing vs AADR genotype (not an admixture fraction).",
            )
        )

    if "Altai_Neanderthal" in freqs and "Chimp" in freqs:
        num, n_num = _f4(mbuti, query, freqs["Altai_Neanderthal"], freqs["Chimp"])
        den, n_den = _f4(mbuti, yoruba, freqs["Altai_Neanderthal"], freqs["Chimp"])
        # Standard-style ratio using Yoruba as a near-zero archaic baseline in the denominator
        # is weak; use Vindija/Altai when both present.
        percent = None
        method = "f4_mbuti_test_altai_chimp"
        if "Vindija_Neanderthal" in freqs:
            den, n_den = _f4(mbuti, freqs["Vindija_Neanderthal"], freqs["Altai_Neanderthal"], freqs["Chimp"])
            method = "f4(Mbuti,Test;Altai,Chimp)/f4(Mbuti,Vindija;Altai,Chimp)"
        if np.isfinite(num) and np.isfinite(den) and abs(den) > 1e-8:
            percent = 100.0 * num / den
        estimates.append(
            HomininEstimate(
                label="Neanderthal ancestry",
                method=method,
                percent=None if percent is None else round(float(np.clip(percent, 0, 100)), 3),
                n_snps=min(n_num, n_den),
                detail={"f4_num": num, "f4_den": den},
                available=percent is not None,
                message=None if percent is not None else "f4 denominator was zero or undefined.",
            )
        )

    if "Denisova" in freqs and "Chimp" in freqs:
        num, n_num = _f4(mbuti, query, freqs["Denisova"], freqs["Chimp"])
        den, n_den = _f4(mbuti, freqs["Denisova"], freqs["Denisova"], freqs["Chimp"])
        method = "f4(Mbuti,Test;Denisova,Chimp)/f4(Mbuti,Denisova;Denisova,Chimp)"
        if "Altai_Neanderthal" in freqs:
            den, n_den = _f4(mbuti, freqs["Denisova"], freqs["Altai_Neanderthal"], freqs["Chimp"])
            method = "f4(Mbuti,Test;Denisova,Chimp)/f4(Mbuti,Denisova;Altai,Chimp)"
        percent = None
        if np.isfinite(num) and np.isfinite(den) and abs(den) > 1e-8:
            percent = 100.0 * num / den
        estimates.append(
            HomininEstimate(
                label="Denisovan ancestry",
                method=method,
                percent=None if percent is None else round(float(np.clip(percent, 0, 100)), 3),
                n_snps=min(n_num, n_den),
                detail={"f4_num": num, "f4_den": den},
                available=percent is not None,
                message=None if percent is not None else "f4 denominator was zero or undefined.",
            )
        )

    estimates.extend(extra_estimates)
    return ComparisonBlock(
        kind="hominin",
        available=any(item.available and item.percent is not None for item in estimates),
        estimates=estimates,
        missing_files=missing,
        notes=notes,
    )


def _load_informative_tsv(path: Path) -> list[tuple[str, int, str]]:
    rows = []
    with path.open() as handle:
        header = handle.readline()
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 5:
                continue
            chrom = normalize_chrom(cols[0])
            pos = int(cols[1])
            archaic_allele = cols[4].upper()
            rows.append((chrom, pos, archaic_allele))
    return rows


def _estimates_from_optional_files(query_index, settings: Settings) -> list[HomininEstimate]:
    estimates = []
    mapping = {
        "neanderthal_informative_snps": "Neanderthal informative-site carrying rate",
        "denisovan_informative_snps": "Denisovan informative-site carrying rate",
    }
    for key, label in mapping.items():
        spec = HOMININ_DOWNLOADS[key]
        path = settings.hominin_dir / spec["filename"]
        if not path.exists():
            continue
        sites = _load_informative_tsv(path)
        hits = 0
        used = 0
        for chrom, pos, archaic_allele in sites:
            row = query_index.get((chrom, pos))
            if row is None:
                continue
            used += 1
            if row["alt"].upper() == archaic_allele and (row.get("dosage_alt") or 0) > 0:
                hits += 1
            elif row["ref"].upper() == archaic_allele and row.get("dosage_alt") == 0:
                hits += 1
        percent = None if used == 0 else 100.0 * hits / used
        estimates.append(
            HomininEstimate(
                label=label,
                method="informative_site_carrier_rate",
                percent=None if percent is None else round(percent, 3),
                n_snps=used,
                available=percent is not None,
                detail={"file": path.name},
            )
        )
    return estimates
