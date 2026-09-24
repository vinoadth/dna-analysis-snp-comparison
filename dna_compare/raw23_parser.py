from __future__ import annotations

import io
import re
from collections import Counter
from pathlib import Path
from typing import BinaryIO, Iterator, TextIO

from dna_compare.config import VARIANT_PREVIEW_LIMIT
from dna_compare.models import VariantRow, VcfSummary
from dna_compare.vcf_parser import chrom_sort_key, normalize_chrom, preview_rows

_BASES = frozenset("ACGT")
_HEADER_COLS = ("rsid", "chromosome", "position", "genotype")
_MIN_SNP_ROWS = 1000
_CHROM_RE = re.compile(r"^(?:[1-9]|1[0-9]|2[0-2]|X|Y|MT|M)$", re.I)
# Autosomal indels kept for blood-type inference (23andMe reports these as I/D).
_BLOOD_TYPE_INDEL_IDS = frozenset({"rs8176719", "i4001527"})


def _valid_rsid(rsid: str) -> bool:
    token = (rsid or "").strip()
    return bool(token) and not any(ch in token for ch in " \t")


class Raw23ParseError(ValueError):
    """Raised when a file is not valid 23andMe-compatible raw genotype data."""


def _pick_alt(ref: str) -> str:
    for base in "ACGT":
        if base != ref:
            return base
    return "N"


def _decode_genotype(raw: str) -> tuple[str, str, str, float | None, bool]:
    """Map a raw genotype token to ref, alt, VCF-style GT, dosage_alt, is_snp."""
    gt = (raw or "").strip().upper()
    if gt in {"", ".", "--"}:
        return (".", ".", ".", None, False)
    if gt in {"II", "I"}:
        return ("I", "D", "0/0", 0.0, False)
    if gt in {"DD", "D"}:
        return ("I", "D", "1/1", 2.0, False)
    if gt in {"DI", "ID"}:
        return ("I", "D", "0/1", 1.0, False)
    if len(gt) == 1 and gt in _BASES:
        ref = gt
        alt = _pick_alt(ref)
        return (ref, alt, "0/0", 0.0, True)
    if len(gt) == 2 and gt[0] in _BASES and gt[1] in _BASES:
        a, b = gt[0], gt[1]
        if a == b:
            ref, alt = a, _pick_alt(a)
            return (ref, alt, "0/0", 0.0, True)
        return (a, b, "0/1", 1.0, True)
    return (".", ".", ".", None, False)


def _open_text(source: Path | str | TextIO | BinaryIO) -> tuple[TextIO, str, bool]:
    filename = "upload.txt"
    close = False
    if isinstance(source, (str, Path)):
        path = Path(source)
        filename = path.name
        handle = path.open("rt", encoding="utf-8", errors="replace")
        close = True
        return handle, filename, close
    name = getattr(source, "name", filename)
    if isinstance(name, str):
        filename = Path(name).name
    raw = source.read()
    if isinstance(raw, bytes):
        handle = io.StringIO(raw.decode("utf-8", errors="replace"))
        close = True
    else:
        handle = io.StringIO(raw)
        close = True
    return handle, filename, close


def _peek_lines(handle: TextIO, limit: int = 50) -> list[str]:
    pos = handle.tell()
    lines = []
    for _ in range(limit):
        line = handle.readline()
        if not line:
            break
        lines.append(line.rstrip("\n\r"))
    handle.seek(pos)
    return lines


def _looks_like_vcf(lines: list[str]) -> bool:
    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        return stripped.startswith("##fileformat=") or stripped.startswith("#CHROM")
    return False


def _header_matches(line: str) -> bool:
    if not line.startswith("#"):
        return False
    body = line.lstrip("#").strip().lower()
    if "rsid" not in body or "genotype" not in body:
        return False
    cols = [part.strip() for part in re.split(r"[\t,|]", body) if part.strip()]
    if len(cols) < 4:
        return False
    normalized = [cols[0], cols[1], cols[2], cols[3]]
    return (
        normalized[0] == "rsid"
        and normalized[1] in {"chromosome", "chrom"}
        and normalized[2] in {"position", "pos"}
        and normalized[3] == "genotype"
    )


def _split_fields(line: str, *, require_tabs: bool = False) -> list[str]:
    if "\t" in line:
        return line.split("\t")
    if require_tabs:
        return [line]
    if "," in line and line.count(",") >= 3:
        return line.split(",")
    return [line]


def _has_raw23_banner(comments: list[str]) -> bool:
    blob = " ".join(comments).lower()
    return any(token in blob for token in ("23andme", "gedmatch", "raw genotype", "raw data"))


def validate_raw23_lines(
    lines: Iterator[str],
    *,
    min_snp_rows: int = _MIN_SNP_ROWS,
) -> tuple[list[str], int, int]:
    """Validate 23andMe-compatible raw data. Returns (header_comments, n_valid, first_bad_line)."""
    header_comments: list[str] = []
    saw_header = False
    n_valid = 0
    n_checked = 0
    first_bad: tuple[int, str] | None = None

    for line_no, raw in enumerate(lines, start=1):
        line = raw.rstrip("\n\r")
        if not line.strip():
            continue
        if line.startswith("#"):
            header_comments.append(line)
            if _header_matches(line):
                saw_header = True
            continue

        if "\t" not in line:
            if first_bad is None:
                first_bad = (
                    line_no,
                    f"line {line_no} is not tab-separated (23andMe raw exports use tabs between rsid, chromosome, position, and genotype)",
                )
            continue

        fields = _split_fields(line.strip(), require_tabs=True)
        if len(fields) != 4:
            if first_bad is None:
                first_bad = (
                    line_no,
                    f"line {line_no} has {len(fields)} tab-separated fields (expected 4: rsid, chromosome, position, genotype)",
                )
            continue

        rsid, chrom, pos_raw, genotype = (field.strip() for field in fields)
        if not _valid_rsid(rsid):
            if first_bad is None:
                first_bad = (line_no, f"line {line_no} has an empty or malformed rsid")
            continue
        chrom_norm = normalize_chrom(chrom)
        if not _CHROM_RE.match(chrom_norm):
            if first_bad is None:
                first_bad = (line_no, f"line {line_no} has invalid chromosome {chrom!r}")
            continue
        try:
            pos = int(pos_raw)
        except ValueError:
            if first_bad is None:
                first_bad = (line_no, f"line {line_no} has non-numeric position {pos_raw!r}")
            continue
        if pos <= 0:
            if first_bad is None:
                first_bad = (line_no, f"line {line_no} has non-positive position {pos}")
            continue

        ref, alt, _gt, _dosage, is_snp = _decode_genotype(genotype)
        if not is_snp and ref == ".":
            continue

        n_valid += 1
        n_checked += 1
        if n_checked >= 200 and not saw_header and n_valid < 5:
            break

    if not saw_header and not _has_raw23_banner(header_comments):
        raise Raw23ParseError(
            "Not a 23andMe-compatible raw file: missing header "
            "'# rsid\\tchromosome\\tposition\\tgenotype' (or a 23andMe/GEDmatch comment block)."
        )
    if not saw_header and first_bad is not None and n_valid == 0:
        raise Raw23ParseError(
            "Not a 23andMe-compatible raw file: expected header "
            "'# rsid\\tchromosome\\tposition\\tgenotype'. "
            f"{first_bad[1]}."
        )

    if first_bad is not None:
        raise Raw23ParseError(f"Not a 23andMe-compatible raw file: {first_bad[1]}.")

    if n_valid < min_snp_rows:
        raise Raw23ParseError(
            f"Not a 23andMe-compatible raw file: too few SNP rows ({n_valid}); "
            f"expected at least {min_snp_rows}."
        )

    return header_comments, n_valid, first_bad[0] if first_bad else 0


def validate_raw23(
    source: Path | str | TextIO | BinaryIO,
    *,
    min_snp_rows: int = _MIN_SNP_ROWS,
) -> None:
    handle, _filename, close = _open_text(source)
    try:
        peek = _peek_lines(handle)
        if _looks_like_vcf(peek):
            raise Raw23ParseError(
                "This file looks like a VCF (starts with ##fileformat=). Upload it as a VCF, not raw 23andMe text."
            )
        if not any(_header_matches(line) for line in peek):
            non_comment = [line for line in peek if line.strip() and not line.startswith("#")]
            if non_comment:
                fields = _split_fields(non_comment[0].strip())
                if len(fields) != 4:
                    raise Raw23ParseError(
                        "Not a 23andMe-compatible raw file: expected header "
                        "'# rsid\\tchromosome\\tposition\\tgenotype' and tab-separated data rows."
                    )
        validate_raw23_lines(handle, min_snp_rows=min_snp_rows)
    finally:
        if close:
            handle.close()


def _assembly_from_comments(comments: list[str]) -> tuple[str | None, str | None]:
    blob = " ".join(comments).lower()
    if "grch38" in blob or "hg38" in blob or "build 38" in blob:
        return "GRCh38", "23andMe raw header mentions GRCh38/hg38"
    if "grch37" in blob or "hg19" in blob or "build 37" in blob:
        return "GRCh37", "23andMe raw header mentions GRCh37/hg19"
    return None, None


def parse_raw23(
    source: Path | str | TextIO | BinaryIO,
    *,
    preview_limit: int = VARIANT_PREVIEW_LIMIT,
    min_snp_rows: int = _MIN_SNP_ROWS,
) -> tuple[VcfSummary, dict[tuple[str, int], dict]]:
    handle, filename, close = _open_text(source)
    try:
        peek = _peek_lines(handle)
        if _looks_like_vcf(peek):
            raise Raw23ParseError(
                "This file looks like a VCF (starts with ##fileformat=). Upload it as a VCF, not raw 23andMe text."
            )
        handle.seek(0)
        header_comments, _n_valid, _ = validate_raw23_lines(handle, min_snp_rows=min_snp_rows)
        handle.seek(0)

        preview: list[VariantRow] = []
        index: dict[tuple[str, int], dict] = {}
        chrom_counts: Counter[str] = Counter()
        n_records = 0
        n_snps = 0
        n_skip = 0
        sample_id = Path(filename).stem

        for raw in handle:
            line = raw.rstrip("\n\r")
            if not line.strip() or line.startswith("#"):
                continue
            if "\t" not in line:
                continue
            fields = _split_fields(line.strip(), require_tabs=True)
            if len(fields) != 4:
                continue
            rsid, chrom_raw, pos_raw, genotype_raw = (field.strip() for field in fields)
            if not _valid_rsid(rsid):
                continue
            chrom = normalize_chrom(chrom_raw)
            if not _CHROM_RE.match(chrom):
                continue
            try:
                pos = int(pos_raw)
            except ValueError:
                continue
            if pos <= 0:
                continue

            ref, alt, gt, dosage, is_snp = _decode_genotype(genotype_raw)
            n_records += 1
            row = {
                "chrom": chrom,
                "pos": pos,
                "rsid": rsid,
                "ref": ref,
                "alt": alt,
                "alt_all": alt,
                "is_snp": is_snp,
                "genotype": gt,
                "dosage_alt": dosage,
                "vcf_filter": ".",
                "n_samples": 1,
                "qual": None,
                "gq": None,
                "dp": None,
                "igc": None,
                "ad": None,
            }
            if not is_snp:
                n_skip += 1
                keep = chrom in {"Y", "MT"} or rsid.strip().lower() in _BLOOD_TYPE_INDEL_IDS
                if keep:
                    index[(chrom, pos)] = row
                continue

            n_snps += 1
            chrom_counts[chrom] += 1
            index[(chrom, pos)] = row
            if len(preview) < preview_limit:
                preview.append(
                    VariantRow(
                        chrom=chrom,
                        pos=pos,
                        rsid=rsid,
                        ref=ref,
                        alt=alt,
                        genotype=gt,
                        dosage_alt=dosage,
                    )
                )

        assembly_hint, assembly_evidence = _assembly_from_comments(header_comments)
        summary = VcfSummary(
            sample_id=sample_id,
            n_records=n_records,
            n_snps=n_snps,
            n_non_snp_skipped=n_skip,
            n_samples_in_file=1,
            chrom_counts=dict(sorted(chrom_counts.items(), key=lambda kv: chrom_sort_key(kv[0]))),
            preview=preview,
            reference=assembly_hint,
            source="23andMe-compatible raw data",
            contig_lengths={},
        )
        if assembly_hint:
            summary.assembly = assembly_hint
            summary.assembly_evidence = assembly_evidence or ""
        if len(preview) < preview_limit and index:
            summary.preview = preview_rows(index, preview_limit, stratified=True)
        return summary, index
    finally:
        if close:
            handle.close()
