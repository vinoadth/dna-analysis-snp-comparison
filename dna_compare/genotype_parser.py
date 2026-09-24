from __future__ import annotations

from pathlib import Path
from typing import BinaryIO, TextIO

from dna_compare.config import VARIANT_PREVIEW_LIMIT
from dna_compare.models import VcfSummary
from dna_compare.raw23_parser import Raw23ParseError, parse_raw23, validate_raw23
from dna_compare.vcf_parser import parse_vcf


def is_raw23_filename(filename: str | None) -> bool:
    if not filename:
        return False
    name = Path(filename).name.lower()
    return name.endswith(".txt") or name.endswith(".23andme")


def parse_genotype_file(
    source: Path | str | TextIO | BinaryIO,
    *,
    filename: str | None = None,
    preview_limit: int = VARIANT_PREVIEW_LIMIT,
    min_raw23_rows: int = 1000,
) -> tuple[VcfSummary, dict[tuple[str, int], dict]]:
    """Parse a SNP VCF or 23andMe-compatible raw text export."""
    name = filename
    if name is None and isinstance(source, (str, Path)):
        name = Path(source).name
    if name is None:
        name = getattr(source, "name", None)
        if isinstance(name, str):
            name = Path(name).name

    if is_raw23_filename(name):
        return parse_raw23(source, preview_limit=preview_limit, min_snp_rows=min_raw23_rows)

    if isinstance(source, (str, Path)):
        path = Path(source)
        if path.suffix.lower() == ".txt":
            return parse_raw23(source, preview_limit=preview_limit, min_snp_rows=min_raw23_rows)

    return parse_vcf(source, preview_limit=preview_limit)


def validate_genotype_upload(
    source: Path | str | TextIO | BinaryIO,
    *,
    filename: str | None = None,
    min_raw23_rows: int = 1000,
) -> None:
    """Validate an upload before analysis; raises Raw23ParseError for bad raw files."""
    name = filename
    if name is None and isinstance(source, (str, Path)):
        name = Path(source).name
    if is_raw23_filename(name) or (isinstance(source, (str, Path)) and Path(source).suffix.lower() == ".txt"):
        validate_raw23(source, min_snp_rows=min_raw23_rows)
