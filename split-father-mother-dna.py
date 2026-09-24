#!/usr/bin/env python3
"""Split a diploid VCF into paternal and maternal haplotype VCFs.

VCF convention: for a genotype `A|B` (or unphased `A/B`), allele A is written
to the father file and allele B to the mother file. Uniparental chromosomes
are assigned wholly: Y → father, MT → mother. A male haploid X goes to mother.

Each haplotype is encoded as a homozygous diploid genotype (`0/0` or `1/1`)
so the existing `python main.py analyze` pipeline can consume the outputs.

Outputs next to the input (or `--outdir`):
  <input-name>_father.vcf
  <input-name>_mother.vcf

Example:
  python split-father-mother-dna.py data/samples/demo.snps.vcf
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path
from typing import TextIO

EXTRA_META = [
    '##INFO=<ID=POO,Number=1,Type=String,Description="Parent-of-origin split: paternal or maternal haplotype">\n',
    '##split-father-mother-dna.py=<ID=convention,Description="Allele before | or / is paternal; allele after is maternal. Y=father, MT=mother, haploid X=mother">\n',
]


def input_name(path: Path) -> str:
    name = path.name
    lower = name.lower()
    for suffix in (".vcf.gz", ".vcf"):
        if lower.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def output_paths(vcf_path: Path, outdir: Path) -> tuple[Path, Path]:
    stem = input_name(vcf_path)
    return outdir / f"{stem}_father.vcf", outdir / f"{stem}_mother.vcf"


def _open_text(path: Path) -> TextIO:
    if str(path).lower().endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open("rt")


def normalize_chrom_token(chrom: str) -> str:
    token = chrom.strip()
    if token.upper().startswith("CHR"):
        token = token[3:]
    return token.upper()


def chrom_bucket(chrom: str) -> str:
    token = normalize_chrom_token(chrom)
    if token in {"M", "MT", "26"}:
        return "mother"
    if token in {"Y", "24"}:
        return "father"
    if token in {"X", "23"}:
        return "x"
    return "auto"


def parse_gt_alleles(sample_field: str) -> tuple[str | None, str | None, str]:
    """Return (paternal_allele, maternal_allele, gt_style).

    gt_style is 'phased', 'unphased', 'haploid', or 'missing'.
    """
    if not sample_field or sample_field == ".":
        return (None, None, "missing")
    token = sample_field.split(":", 1)[0]
    if not token or token == ".":
        return (None, None, "missing")
    if "|" in token:
        parts = token.split("|")
        style = "phased"
    elif "/" in token:
        parts = token.split("/")
        style = "unphased"
    else:
        parts = [token]
        style = "haploid"
    alleles: list[str | None] = [None if p == "." else p for p in parts]
    if len(alleles) == 1:
        return (alleles[0], alleles[0], style)
    return (alleles[0], alleles[1], style)


def homozygous_gt(allele: str | None) -> str:
    if allele is None:
        return "./."
    return f"{allele}/{allele}"


def rewrite_sample(sample_field: str, allele: str | None, n_fmt: int) -> str:
    gt = homozygous_gt(allele)
    if n_fmt <= 1:
        return gt
    rest = sample_field.split(":")[1:]
    return ":".join([gt, *rest])


def allele_for_side(
    bucket: str,
    side: str,
    father_a: str | None,
    mother_a: str | None,
    style: str,
) -> str | None:
    """Which allele belongs in this output file, or None to omit the site."""
    if bucket == "father":
        if side != "father":
            return None
        return father_a if father_a is not None else mother_a
    if bucket == "mother":
        if side != "mother":
            return None
        return mother_a if mother_a is not None else father_a
    if bucket == "x" and style == "haploid":
        if side != "mother":
            return None
        return father_a if father_a is not None else mother_a
    return father_a if side == "father" else mother_a


def _load_vcf(vcf_path: Path) -> tuple[list[str], str, list[str]]:
    meta: list[str] = []
    header: str | None = None
    body: list[str] = []
    with _open_text(vcf_path) as handle:
        for raw in handle:
            line = raw if raw.endswith("\n") else raw + "\n"
            if line.startswith("##"):
                meta.append(line)
            elif line.startswith("#CHROM"):
                header = line
            elif line.strip():
                body.append(line)
    if header is None:
        raise ValueError(f"No #CHROM header in {vcf_path}")
    return meta, header, body


def _write_header(out, meta: list[str], header_cols: list[str]) -> None:
    if meta:
        out.write(meta[0])
        for line in EXTRA_META:
            out.write(line)
        for line in meta[1:]:
            out.write(line)
    else:
        out.write("##fileformat=VCFv4.2\n")
        for line in EXTRA_META:
            out.write(line)
    out.write("\t".join(header_cols) + "\n")


def split_vcf(
    vcf_path: Path,
    *,
    outdir: Path | None = None,
    sample_index: int = 0,
) -> tuple[Path, Path, dict[str, int]]:
    if not vcf_path.exists():
        raise FileNotFoundError(vcf_path)
    outdir = outdir or vcf_path.parent
    outdir.mkdir(parents=True, exist_ok=True)
    father_path, mother_path = output_paths(vcf_path, outdir)

    meta, header, body = _load_vcf(vcf_path)
    cols_header = header.rstrip("\n").split("\t")
    samples = cols_header[9:]
    if samples:
        if sample_index < 0 or sample_index >= len(samples):
            raise ValueError(f"sample_index {sample_index} out of range for {samples}")
        sample_name = samples[sample_index]
    else:
        sample_name = input_name(vcf_path)

    stats = {
        "records": 0,
        "father": 0,
        "mother": 0,
        "unphased": 0,
        "phased": 0,
        "haploid": 0,
        "skipped_no_gt": 0,
    }

    handles = {
        "father": father_path.open("w"),
        "mother": mother_path.open("w"),
    }
    try:
        _write_header(handles["father"], meta, cols_header[:9] + [f"{sample_name}_father"])
        _write_header(handles["mother"], meta, cols_header[:9] + [f"{sample_name}_mother"])

        for line in body:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            stats["records"] += 1
            bucket = chrom_bucket(cols[0])
            sample_field = cols[9 + sample_index] if len(cols) > 9 + sample_index else "."
            father_a, mother_a, style = parse_gt_alleles(sample_field)
            if style in stats:
                stats[style] += 1
            n_fmt = len(cols[8].split(":")) if len(cols) > 8 else 1
            info = cols[7] if len(cols) > 7 else "."
            fmt = cols[8] if len(cols) > 8 else "GT"

            for side in ("father", "mother"):
                allele = allele_for_side(bucket, side, father_a, mother_a, style)
                if allele is None:
                    continue
                poo = "paternal" if side == "father" else "maternal"
                new_info = f"POO={poo}" if info in {".", ""} else f"{info};POO={poo}"
                new_sample = rewrite_sample(sample_field, allele, n_fmt)
                out_cols = cols[:7] + [new_info, fmt, new_sample]
                handles[side].write("\t".join(out_cols) + "\n")
                stats[side] += 1
            if father_a is None and mother_a is None:
                stats["skipped_no_gt"] += 1
    finally:
        for handle in handles.values():
            handle.close()

    return father_path, mother_path, stats


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Split a VCF into <name>_father.vcf and <name>_mother.vcf haplotype files."
    )
    parser.add_argument("vcf", type=Path, help="Input VCF (.vcf or .vcf.gz)")
    parser.add_argument("-o", "--outdir", type=Path, default=None, help="Output directory (default: next to input)")
    parser.add_argument(
        "--sample-index",
        type=int,
        default=0,
        help="Which sample column to split when the VCF has multiple samples (0-based)",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        father_path, mother_path, stats = split_vcf(
            args.vcf,
            outdir=args.outdir,
            sample_index=args.sample_index,
        )
    except (FileNotFoundError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print(f"input:  {args.vcf}")
    print(f"father: {father_path}  variants={stats['father']}")
    print(f"mother: {mother_path}  variants={stats['mother']}")
    print(
        f"records={stats['records']}  phased={stats['phased']}  "
        f"unphased={stats['unphased']}  haploid={stats['haploid']}  "
        f"skipped_no_gt={stats['skipped_no_gt']}"
    )
    if stats["unphased"] and not stats["phased"]:
        print(
            "note: genotypes are unphased; alleles are split in file order "
            "(first=father, second=mother), not proven parent-of-origin."
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
