from __future__ import annotations

import gzip
import re
from collections import Counter
from pathlib import Path
from typing import BinaryIO, Iterator, TextIO

from dna_compare.config import VARIANT_PREVIEW_LIMIT
from dna_compare.models import VariantRow, VcfSummary

_CONTIG_RE = re.compile(r"ID=([^,>]+).*length=(\d+)", re.I)
_ASM_HINTS = ("grch38", "hg38", "grch37", "hg19", "gsa-24v3", "gtc2vcf", "human_g1k_v37")


def _open_text(path: Path) -> TextIO:
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open("rt")


def normalize_chrom(chrom: str) -> str:
    chrom = chrom.strip()
    if chrom.startswith("chr"):
        chrom = chrom[3:]
    if chrom == "23":
        return "X"
    if chrom == "24":
        return "Y"
    if chrom in {"25", "26", "M", "MT", "chrM"}:
        return "MT"
    return chrom


def _is_snp(ref: str, alt: str) -> bool:
    if ref in {".", ""} or alt in {".", ""}:
        return False
    if "," in alt:
        alt = alt.split(",", 1)[0]
    return len(ref) == 1 and len(alt) == 1 and ref != "*" and alt != "*"


def _parse_gt(gt: str) -> tuple[str, float | None]:
    if not gt or gt == "." or gt.startswith("."):
        return (".", None)
    token = gt.split(":", 1)[0]
    seps = "/" if "/" in token else "|" if "|" in token else None
    if seps is None:
        if token in {"0", "1"}:
            dosage = float(token)
            return (token, dosage)
        return (token, None)
    alleles = token.split(seps)
    if any(a == "." for a in alleles):
        return (token, None)
    try:
        dosage = sum(int(a) for a in alleles)
    except ValueError:
        return (token, None)
    return (token.replace("|", "/"), float(dosage))


def _maybe_int(value: str | None) -> int | None:
    if value in {None, "", "."}:
        return None
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def _maybe_float(value: str | None) -> float | None:
    if value in {None, "", "."}:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _format_fields(fmt: str, sample: str) -> dict[str, str]:
    keys = fmt.split(":") if fmt and fmt != "." else []
    vals = sample.split(":") if sample else []
    return {key: vals[i] if i < len(vals) else "" for i, key in enumerate(keys)}


def _sample_metrics(qual: str, fields: dict[str, str]) -> dict:
    ad = fields.get("AD")
    return {
        "qual": _maybe_float(qual),
        "gq": _maybe_int(fields.get("GQ")),
        "dp": _maybe_int(fields.get("DP")),
        "igc": _maybe_float(fields.get("IGC")),
        "ad": None if ad in {None, "", "."} else ad,
    }


def _ingest_header(line: str, meta: dict) -> None:
    if line.startswith("##reference="):
        meta["reference"] = line[len("##reference=") :].strip().strip('"')
        return
    if line.startswith("##source="):
        sources = meta.setdefault("sources", [])
        sources.append(line[len("##source=") :].strip().strip('"'))
        return
    if line.startswith("##contig="):
        match = _CONTIG_RE.search(line)
        if match:
            meta.setdefault("contig_lengths", {})[normalize_chrom(match.group(1))] = int(match.group(2))
        return
    low = line.lower()
    if any(hint in low for hint in _ASM_HINTS):
        meta.setdefault("sources", []).append(line[2:].strip()[:240])


def iter_vcf_snps(handle: TextIO) -> tuple[list[str], dict, Iterator[dict]]:
    samples: list[str] = []
    meta: dict = {"reference": None, "sources": [], "contig_lengths": {}}

    def _rows() -> Iterator[dict]:
        for raw in handle:
            if raw.startswith("##"):
                _ingest_header(raw.rstrip("\n"), meta)
                continue
            if raw.startswith("#CHROM"):
                header = raw.rstrip("\n").split("\t")
                samples.clear()
                samples.extend(header[9:])
                continue
            if not raw.strip() or raw.startswith("#"):
                continue
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            ref = cols[3]
            alt = cols[4]
            snp = _is_snp(ref, alt)
            alt_primary = alt.split(",", 1)[0]
            gt_raw = cols[9] if len(cols) > 9 else "."
            gt, dosage = _parse_gt(gt_raw)
            fields = _format_fields(cols[8] if len(cols) > 8 else "GT", gt_raw)
            metrics = _sample_metrics(cols[5], fields)
            yield {
                "chrom": normalize_chrom(cols[0]),
                "pos": int(cols[1]),
                "rsid": cols[2],
                "ref": ref,
                "alt": alt_primary,
                "alt_all": alt,
                "is_snp": snp,
                "genotype": gt,
                "dosage_alt": dosage,
                "vcf_filter": cols[6] if len(cols) > 6 else ".",
                "n_samples": max(0, len(cols) - 9),
                **metrics,
            }

    return samples, meta, _rows()


def chrom_sort_key(chrom: str) -> tuple[int, str]:
    key = str(chrom).replace("chr", "").replace("CHR", "").upper()
    if key == "X":
        return (23, key)
    if key == "Y":
        return (24, key)
    if key in {"MT", "M"}:
        return (25, key)
    if key.isdigit():
        return (int(key), key)
    return (100, key)


def _to_variant_row(row: dict) -> VariantRow:
    return VariantRow(
        chrom=str(row.get("chrom") or ""),
        pos=int(row.get("pos") or 0),
        rsid=str(row.get("rsid") or ""),
        ref=str(row.get("ref") or ""),
        alt=str(row.get("alt") or ""),
        genotype=str(row.get("genotype") or ""),
        dosage_alt=row.get("dosage_alt"),
    )


def preview_rows(
    index: dict[tuple[str, int], dict],
    limit: int | None,
    *,
    stratified: bool = False,
) -> list[VariantRow]:
    eligible = [
        row
        for row in index.values()
        if row.get("is_snp") or row.get("chrom") in {"Y", "MT"}
    ]
    eligible.sort(key=lambda row: (chrom_sort_key(str(row["chrom"])), int(row["pos"])))
    if limit is None or limit <= 0 or len(eligible) <= limit:
        return [_to_variant_row(row) for row in eligible]
    if not stratified:
        return [_to_variant_row(row) for row in eligible[:limit]]
    by_chrom: dict[str, list[dict]] = {}
    for row in eligible:
        by_chrom.setdefault(str(row["chrom"]), []).append(row)
    chroms = sorted(by_chrom, key=chrom_sort_key)
    per = max(1, limit // max(1, len(chroms)))
    picked: list[dict] = []
    for chrom in chroms:
        picked.extend(by_chrom[chrom][:per])
        if len(picked) >= limit:
            break
    picked.sort(key=lambda row: (chrom_sort_key(str(row["chrom"])), int(row["pos"])))
    return [_to_variant_row(row) for row in picked[:limit]]


def parse_vcf(
    source: Path | str | TextIO | BinaryIO,
    *,
    preview_limit: int = VARIANT_PREVIEW_LIMIT,
) -> tuple[VcfSummary, dict[tuple[str, int], dict]]:
    """Parse a SNP VCF into a UI summary plus (chrom, pos) lookup.

    `source` may be a path, an open text handle, or a binary buffer (API upload).
    """
    import io

    close = False
    handle: TextIO
    filename = "upload.vcf"
    if isinstance(source, (str, Path)):
        path = Path(source)
        filename = path.name
        handle = _open_text(path)
        close = True
    else:
        name = getattr(source, "name", filename)
        if isinstance(name, str):
            filename = Path(name).name
        raw = source.read()
        if isinstance(raw, bytes):
            if raw[:2] == b"\x1f\x8b":
                handle = gzip.open(io.BytesIO(raw), "rt")
            else:
                handle = io.StringIO(raw.decode("utf-8", errors="replace"))
            close = True
        else:
            handle = io.StringIO(raw)

    try:
        return _parse_open(handle, filename=filename, preview_limit=preview_limit)
    finally:
        if close:
            handle.close()


def _parse_open(handle: TextIO, *, filename: str, preview_limit: int) -> tuple[VcfSummary, dict[tuple[str, int], dict]]:
    samples, meta, rows = iter_vcf_snps(handle)
    preview: list[VariantRow] = []
    index: dict[tuple[str, int], dict] = {}
    chrom_counts: Counter[str] = Counter()
    n_records = 0
    n_snps = 0
    n_skip = 0
    sample_id = Path(filename).stem
    n_file_samples = 0

    for row in rows:
        n_records += 1
        n_file_samples = max(n_file_samples, row["n_samples"])
        if not row["is_snp"]:
            n_skip += 1
            # Keep Y/MT indels (M17 is an insertion) for haplogroup scoring.
            if row["chrom"] in {"Y", "MT"}:
                index[(row["chrom"], row["pos"])] = row
            continue
        n_snps += 1
        chrom_counts[row["chrom"]] += 1
        key = (row["chrom"], row["pos"])
        index[key] = row
        if len(preview) < preview_limit:
            preview.append(
                VariantRow(
                    chrom=row["chrom"],
                    pos=row["pos"],
                    rsid=row["rsid"],
                    ref=row["ref"],
                    alt=row["alt"],
                    genotype=row["genotype"],
                    dosage_alt=row["dosage_alt"],
                )
            )

    if samples:
        sample_id = samples[0]
    sources = meta.get("sources") or []
    summary = VcfSummary(
        sample_id=sample_id,
        n_records=n_records,
        n_snps=n_snps,
        n_non_snp_skipped=n_skip,
        n_samples_in_file=n_file_samples or len(samples),
        chrom_counts=dict(sorted(chrom_counts.items(), key=lambda kv: kv[0])),
        preview=preview,
        reference=meta.get("reference"),
        source="; ".join(sources) if sources else None,
        contig_lengths=dict(meta.get("contig_lengths") or {}),
    )
    return summary, index
