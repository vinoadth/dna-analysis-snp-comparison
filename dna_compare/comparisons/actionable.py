from __future__ import annotations

from functools import lru_cache
from pathlib import Path

from dna_compare.comparisons.additional import _find_row, _lookups
from dna_compare.comparisons.drugs import _read_tsv, effect_copies
from dna_compare.config import Settings
from dna_compare.models import AdditionalDetail, ComparisonBlock
from dna_compare.vcf_parser import normalize_chrom

ACTIONABLE_NOTES = [
    "Research list of curated coding variants from an exome. It is not a clinical report and it is not a search of every pathogenic variant.",
    "Rows appear only when the variant allele is called. A site that is absent from the file is not shown and is not a negative result.",
    "Recessive diseases: one copy is a carrier result, two copies is the disease genotype. Dominant rows are risk alleles. This is not a diagnosis.",
]


@lru_cache(maxsize=4)
def load_actionable_variants(path: str) -> tuple[dict, ...]:
    rows = []
    for rec in _read_tsv(Path(path)):
        rsid = (rec.get("rsid") or "").strip()
        if not rsid:
            continue
        try:
            pos_hg19 = int(rec.get("pos_hg19") or 0)
        except ValueError:
            pos_hg19 = 0
        try:
            pos_hg38 = int(rec.get("pos_hg38") or 0)
        except ValueError:
            pos_hg38 = 0
        rows.append(
            {
                "gene": (rec.get("gene") or "").strip(),
                "allele": (rec.get("allele") or "").strip(),
                "rsid": rsid,
                "chrom": normalize_chrom(rec.get("chrom") or ""),
                "pos_hg19": pos_hg19,
                "pos_hg38": pos_hg38,
                "effect_allele": (rec.get("effect_allele") or "").strip().upper(),
                "other_allele": (rec.get("other_allele") or "").strip().upper(),
                "condition": (rec.get("condition") or "").strip(),
                "source": (rec.get("source") or "").strip(),
                "inheritance": (rec.get("inheritance") or "dominant").strip().lower() or "dominant",
            }
        )
    return tuple(rows)


def _zygosity(inheritance: str, n: float, *, male: bool) -> str:
    if inheritance == "recessive":
        return "homozygous — disease genotype" if n >= 2 else "carrier (heterozygous)"
    if inheritance == "x-linked":
        if male and n >= 1:
            return "hemizygous"
        return "homozygous" if n >= 2 else "heterozygous"
    return "homozygous" if n >= 2 else "heterozygous"


def compare_actionable(query_index: dict[tuple[str, int], dict], *, settings: Settings) -> ComparisonBlock:
    path = Path(settings.actionable_variants)
    notes = list(ACTIONABLE_NOTES)
    if not path.is_file():
        return ComparisonBlock(
            kind="actionable",
            available=False,
            notes=notes + [f"Missing variant table: {path}."],
        )
    markers = load_actionable_variants(str(path))
    lookups = _lookups(query_index)
    by_rsid, by_pos = lookups
    male = any(chrom == "Y" for chrom, _pos in query_index)
    hits: list[AdditionalDetail] = []
    typed = 0
    for marker in markers:
        row = _find_row(marker, by_rsid, by_pos)
        if row is None:
            continue
        typed += 1
        got = effect_copies(row, marker["effect_allele"], marker["other_allele"])
        if got is None:
            continue
        n, genotype = got
        if n < 1:
            continue
        zygosity = _zygosity(marker.get("inheritance") or "dominant", n, male=male)
        hits.append(
            AdditionalDetail(
                key=marker["rsid"],
                topic=marker["condition"] or marker["gene"],
                finding=f"{marker['gene']} {marker['allele']} {zygosity}.",
                evidence=f"{marker['rsid']} {genotype} ({n:.0f}× {marker['effect_allele']})",
                source=marker["source"],
                available=True,
                confidence="high" if n >= 2 else "medium",
                n_snps=1,
                n_markers=1,
                coverage_status="ready",
                parent=zygosity,
            )
        )
    notes.append(
        f"{len(hits)} curated variant allele(s) called; {typed} of {len(markers)} listed sites were present in this file."
    )
    return ComparisonBlock(kind="actionable", available=True, estimates=hits, notes=notes)
