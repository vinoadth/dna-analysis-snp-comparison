#!/usr/bin/env python3
"""Merge a WES VCF with a SNP-array VCF into one single-sample VCF.

The WES VCF is the base: every WES record is kept unchanged. SNP-array records
are added only where WES has nothing to say:
  - the (chrom, pos) is absent from the WES VCF, or
  - the WES record at that (chrom, pos) is a no-call (`./.`) and the array
    has a called genotype (the WES row is then replaced by the array row).

SNP-array no-calls are never added. Chromosome names follow the WES naming
style (`chr1` vs `1`, `chrM` vs `MT`), records are sorted by chrom/pos, and
every record gets `INFO/SRC=WES` or `INFO/SRC=SNP`. Both inputs must be on
the same genome build; the merge refuses to mix GRCh37 and GRCh38.

The output has one sample column, so it can go straight into
`python main.py analyze`.

Output next to the WES input (or `--outdir`):
  <wes-name>_wes_snp_merged.vcf

Example:
  python merge-wes-snp-vcf.py sample.wes.vcf.gz sample.snps.vcf
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import TextIO

from dna_compare.assembly import detect_assembly
from dna_compare.models import VcfSummary
from dna_compare.vcf_parser import _ingest_header, chrom_sort_key

TOOL = "merge-wes-snp-vcf.py"
SRC_INFO = '##INFO=<ID=SRC,Number=1,Type=String,Description="Record source: WES (base) or SNP (filled from SNP-array VCF)">\n'
_META_ID_RE = re.compile(r"^##(INFO|FORMAT|FILTER|ALT|contig)=<ID=([^,>]+)", re.I)


@dataclass
class LoadedVcf:
    path: Path
    meta: list[str]
    sample_name: str
    records: dict[tuple[str, int], list[list[str]]] = field(default_factory=dict)
    raw_chrom_names: dict[str, str] = field(default_factory=dict)
    assembly: str = "unknown"
    assembly_evidence: str = ""


def input_name(path: Path) -> str:
    name = path.name
    lower = name.lower()
    for suffix in (".vcf.gz", ".vcf"):
        if lower.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def output_path(wes_path: Path, outdir: Path) -> Path:
    return outdir / f"{input_name(wes_path)}_wes_snp_merged.vcf"


def _open_text(path: Path) -> TextIO:
    if str(path).lower().endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open("rt")


def normalize_chrom_token(chrom: str) -> str:
    token = chrom.strip()
    if token.upper().startswith("CHR"):
        token = token[3:]
    token = token.upper()
    if token == "23":
        return "X"
    if token == "24":
        return "Y"
    if token in {"M", "MT", "25", "26"}:
        return "MT"
    return token


def is_no_call(cols: list[str]) -> bool:
    if len(cols) < 10:
        return True
    gt = cols[9].split(":", 1)[0]
    return not gt or all(allele == "." for allele in re.split(r"[/|]", gt))


def _load_vcf(path: Path, sample_index: int) -> LoadedVcf:
    if not path.exists():
        raise FileNotFoundError(path)
    meta: list[str] = []
    header_meta: dict = {"reference": None, "sources": [], "contig_lengths": {}}
    header: list[str] | None = None
    loaded: LoadedVcf | None = None
    with _open_text(path) as handle:
        for raw in handle:
            line = raw if raw.endswith("\n") else raw + "\n"
            if line.startswith("##"):
                meta.append(line)
                _ingest_header(line.rstrip("\n"), header_meta)
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                samples = header[9:]
                if samples and not 0 <= sample_index < len(samples):
                    raise ValueError(f"sample index {sample_index} out of range for {samples} in {path}")
                sample_name = samples[sample_index] if samples else input_name(path)
                loaded = LoadedVcf(path=path, meta=meta, sample_name=sample_name)
                continue
            if not line.strip() or line.startswith("#"):
                continue
            if loaded is None:
                raise ValueError(f"No #CHROM header in {path}")
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            sample = cols[9 + sample_index] if len(cols) > 9 + sample_index else "."
            fmt = cols[8] if len(cols) > 8 else "GT"
            cols = cols[:8] + [fmt, sample]
            chrom = normalize_chrom_token(cols[0])
            loaded.raw_chrom_names.setdefault(chrom, cols[0])
            loaded.records.setdefault((chrom, int(cols[1])), []).append(cols)
    if loaded is None:
        raise ValueError(f"No #CHROM header in {path}")

    summary = VcfSummary(
        sample_id=loaded.sample_name,
        n_records=0,
        n_snps=0,
        n_non_snp_skipped=0,
        n_samples_in_file=0,
        reference=header_meta.get("reference"),
        source="; ".join(header_meta.get("sources") or []) or None,
        contig_lengths=dict(header_meta.get("contig_lengths") or {}),
    )
    detect_assembly(summary)
    loaded.assembly = summary.assembly
    loaded.assembly_evidence = summary.assembly_evidence
    return loaded


class ChromNamer:
    """Render normalized chrom tokens in the WES file's naming style."""

    def __init__(self, wes: LoadedVcf):
        self.known = dict(wes.raw_chrom_names)
        for line in wes.meta:
            match = _META_ID_RE.match(line)
            if match and match.group(1).lower() == "contig":
                self.known.setdefault(normalize_chrom_token(match.group(2)), match.group(2))
        names = list(self.known.values())
        self.use_chr = bool(names) and sum(n.lower().startswith("chr") for n in names) * 2 >= len(names)

    def __call__(self, chrom: str) -> str:
        if chrom in self.known:
            return self.known[chrom]
        if not self.use_chr:
            return chrom
        return "chrM" if chrom == "MT" else f"chr{chrom}"


def _meta_key(line: str) -> tuple[str, str] | None:
    match = _META_ID_RE.match(line)
    if not match:
        return None
    kind = match.group(1).lower()
    ident = normalize_chrom_token(match.group(2)) if kind == "contig" else match.group(2)
    return kind, ident


def merged_meta(wes: LoadedVcf, snp: LoadedVcf, namer: ChromNamer) -> list[str]:
    wes_meta = [line for line in wes.meta if not line.startswith(("##INFO=<ID=SRC,", f"##{TOOL}"))]
    if not wes_meta or not wes_meta[0].startswith("##fileformat"):
        wes_meta.insert(0, "##fileformat=VCFv4.2\n")
    seen = {key for line in wes_meta if (key := _meta_key(line))}
    extra: list[str] = []
    for line in snp.meta:
        key = _meta_key(line)
        if key is None or key in seen or key == ("info", "SRC"):
            continue
        seen.add(key)
        if key[0] == "contig":
            line = re.sub(r"ID=[^,>]+", f"ID={namer(key[1])}", line, count=1)
        extra.append(line)
    tool_line = (
        f'##{TOOL}=<ID=merge,Description="Base={wes.path.name}; filled from {snp.path.name} '
        f'at positions absent from WES or WES no-calls">\n'
    )
    return [*wes_meta, *extra, SRC_INFO, tool_line]


def _with_src(cols: list[str], src: str, chrom_name: str | None = None) -> list[str]:
    info = cols[7]
    new_info = f"SRC={src}" if info in {".", ""} else f"{info};SRC={src}"
    return [chrom_name or cols[0], *cols[1:7], new_info, *cols[8:]]


def merge_vcfs(
    wes_path: Path,
    snp_path: Path,
    *,
    outdir: Path | None = None,
    wes_sample_index: int = 0,
    snp_sample_index: int = 0,
    allow_unknown_assembly: bool = True,
) -> tuple[Path, dict[str, int | str]]:
    wes = _load_vcf(wes_path, wes_sample_index)
    snp = _load_vcf(snp_path, snp_sample_index)

    known = {wes.assembly, snp.assembly} - {"unknown"}
    if len(known) > 1:
        raise ValueError(
            f"assembly mismatch: WES is {wes.assembly} ({wes.assembly_evidence}), "
            f"SNP is {snp.assembly} ({snp.assembly_evidence}); lift one file over first"
        )
    if not allow_unknown_assembly and "unknown" in {wes.assembly, snp.assembly}:
        raise ValueError(
            f"could not detect assembly (WES={wes.assembly}, SNP={snp.assembly}); "
            "pass --allow-unknown-assembly if both files are known to be on the same build"
        )

    outdir = outdir or wes_path.parent
    outdir.mkdir(parents=True, exist_ok=True)
    out_path = output_path(wes_path, outdir)
    namer = ChromNamer(wes)

    stats: dict[str, int | str] = {
        "wes_records": sum(len(rows) for rows in wes.records.values()),
        "snp_records": sum(len(rows) for rows in snp.records.values()),
        "kept_wes": 0,
        "added_snp": 0,
        "replaced_wes_nocall": 0,
        "snp_overlap_skipped": 0,
        "snp_nocall_skipped": 0,
        "ref_mismatch": 0,
        "wes_assembly": wes.assembly,
        "snp_assembly": snp.assembly,
    }

    merged: dict[tuple[str, int], list[list[str]]] = {}
    for key, rows in wes.records.items():
        merged[key] = [_with_src(cols, "WES") for cols in rows]
        stats["kept_wes"] += len(rows)

    for key, rows in snp.records.items():
        called = [cols for cols in rows if not is_no_call(cols)]
        stats["snp_nocall_skipped"] += len(rows) - len(called)
        if not called:
            continue
        wes_rows = wes.records.get(key)
        if wes_rows:
            if wes_rows[0][3].upper() != called[0][3].upper():
                stats["ref_mismatch"] += 1
            if not all(is_no_call(cols) for cols in wes_rows):
                stats["snp_overlap_skipped"] += len(called)
                continue
            stats["kept_wes"] -= len(wes_rows)
            stats["replaced_wes_nocall"] += len(wes_rows)
        merged[key] = [_with_src(cols, "SNP", namer(key[0])) for cols in called]
        stats["added_snp"] += len(called)

    header_cols = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", wes.sample_name]
    with out_path.open("w") as out:
        out.writelines(merged_meta(wes, snp, namer))
        out.write("\t".join(header_cols) + "\n")
        for key in sorted(merged, key=lambda k: (chrom_sort_key(k[0]), k[1])):
            for cols in merged[key]:
                out.write("\t".join(cols) + "\n")

    stats["output_records"] = sum(len(rows) for rows in merged.values())
    return out_path, stats


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Merge a WES VCF (base) with a SNP-array VCF, filling positions WES does not cover."
    )
    parser.add_argument("wes_vcf", type=Path, help="WES VCF (.vcf or .vcf.gz); used as the base")
    parser.add_argument("snp_vcf", type=Path, help="SNP-array VCF (.vcf or .vcf.gz); fills the gaps")
    parser.add_argument("-o", "--outdir", type=Path, default=None, help="Output directory (default: next to WES input)")
    parser.add_argument("--wes-sample-index", type=int, default=0, help="Sample column in the WES VCF (0-based)")
    parser.add_argument("--snp-sample-index", type=int, default=0, help="Sample column in the SNP VCF (0-based)")
    parser.add_argument(
        "--require-assembly",
        action="store_true",
        help="Fail unless both files' genome build can be detected from the header",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        out_path, stats = merge_vcfs(
            args.wes_vcf,
            args.snp_vcf,
            outdir=args.outdir,
            wes_sample_index=args.wes_sample_index,
            snp_sample_index=args.snp_sample_index,
            allow_unknown_assembly=not args.require_assembly,
        )
    except (FileNotFoundError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print(f"wes:    {args.wes_vcf}  records={stats['wes_records']}  assembly={stats['wes_assembly']}")
    print(f"snp:    {args.snp_vcf}  records={stats['snp_records']}  assembly={stats['snp_assembly']}")
    print(f"merged: {out_path}  records={stats['output_records']}")
    print(
        f"kept_wes={stats['kept_wes']}  added_snp={stats['added_snp']}  "
        f"replaced_wes_nocall={stats['replaced_wes_nocall']}  "
        f"snp_overlap_skipped={stats['snp_overlap_skipped']}  "
        f"snp_nocall_skipped={stats['snp_nocall_skipped']}"
    )
    if "unknown" in {stats["wes_assembly"], stats["snp_assembly"]}:
        print("note: genome build could not be detected for one input; make sure both are on the same build.")
    if stats["ref_mismatch"]:
        print(
            f"warning: {stats['ref_mismatch']} shared positions have different REF alleles; "
            "check that both files use the same reference build."
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
