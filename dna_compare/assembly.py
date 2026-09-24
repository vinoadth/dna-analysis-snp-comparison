from __future__ import annotations

import re

from dna_compare.models import VcfSummary

# Distinctive assembled lengths. MT is 16569 on both and is not used.
_HG19_LEN = {
    "1": 249250621,
    "2": 243199373,
    "X": 155270560,
    "Y": 59373566,
}
_HG38_LEN = {
    "1": 248956422,
    "2": 242193529,
    "X": 156040895,
    "Y": 57227415,
}

_GRCH38 = re.compile(r"\b(grch38|hg38|b38)\b", re.I)
_GRCH37 = re.compile(r"\b(grch37|hg19|b37)\b", re.I)
_GSA_V3 = re.compile(r"gsa-24v3|gtc2vcf", re.I)


def _header_blob(summary: VcfSummary) -> str:
    parts = [summary.reference or "", summary.source or ""]
    return " ".join(part for part in parts if part)


def _contig_vote(lengths: dict[str, int]) -> tuple[str | None, str]:
    hg19 = hg38 = 0
    matched: list[str] = []
    for chrom, size in lengths.items():
        if _HG19_LEN.get(chrom) == size:
            hg19 += 1
            matched.append(f"chr{chrom}={size}")
        if _HG38_LEN.get(chrom) == size:
            hg38 += 1
            matched.append(f"chr{chrom}={size}")
    if hg38 > hg19 and hg38 > 0:
        return "GRCh38", "contig lengths " + ", ".join(matched[:3])
    if hg19 > hg38 and hg19 > 0:
        return "GRCh37", "contig lengths " + ", ".join(matched[:3])
    return None, ""


def detect_assembly(summary: VcfSummary, assume: str | None = None) -> VcfSummary:
    """Fill assembly / evidence on a parsed VCF summary. Does not rewrite positions."""
    if assume in {"GRCh37", "GRCh38"}:
        summary.assembly = assume
        summary.assembly_evidence = f"explicit --assembly {assume}"
        return summary

    blob = _header_blob(summary)
    contig_assembly, contig_why = _contig_vote(summary.contig_lengths)
    if contig_assembly:
        summary.assembly = contig_assembly
        summary.assembly_evidence = contig_why
        return summary
    if _GRCH38.search(blob) and not _GRCH37.search(blob):
        summary.assembly = "GRCh38"
        summary.assembly_evidence = "VCF ##reference/##source names GRCh38/hg38"
        return summary
    if _GRCH37.search(blob) and not _GRCH38.search(blob):
        summary.assembly = "GRCh37"
        summary.assembly_evidence = "VCF ##reference/##source names GRCh37/hg19"
        return summary
    if _GSA_V3.search(blob):
        # Default Illumina GSA v3 BPM/manifest is GRCh38; gtc2vcf does not lift.
        summary.assembly = "GRCh38"
        summary.assembly_evidence = "GSA-24v3 / bcftools_gtc2vcf header (default manifest is GRCh38)"
        return summary
    summary.assembly = "unknown"
    summary.assembly_evidence = "no contig lengths, reference, or GSA/gtc2vcf source in the header"
    return summary


def assembly_note(summary: VcfSummary) -> str | None:
    if summary.lifted_to:
        return (
            f"Coordinates lifted {summary.assembly} → {summary.lifted_to} "
            f"({summary.n_lifted} sites; {summary.n_unmapped} unmapped) so they match AADR v66.1 HO / hg19. "
            f"Detected from {summary.assembly_evidence}."
        )
    if summary.assembly == "GRCh38":
        return (
            "This VCF looks like GRCh38 (hg38), usually a GSA v3 / gtc2vcf export, while AADR HO is GRCh37 (hg19). "
            "Y haplogroups still match by rsID and published hg38 positions. Autosomal AADR overlap needs "
            "data/references/liftover/hg38ToHg19.over.chain.gz."
        )
    if summary.assembly == "GRCh37":
        extra = ""
        blob = f"{summary.source or ''} {summary.reference or ''}".lower()
        if "gsa-24v3" in blob or "gtc2vcf" in blob:
            extra = " GSA/gtc2vcf is present, but the FASTA/contigs are GRCh37, so coordinates were not lifted."
        return f"VCF treated as GRCh37/hg19 ({summary.assembly_evidence}), same assembly as AADR HO.{extra}"
    return None
