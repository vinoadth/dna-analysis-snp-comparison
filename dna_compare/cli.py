from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from webbrowser import open as open_browser

from dna_compare.config import (
    CASTE_DOWNLOADS,
    HOMININ_DOWNLOADS,
    POPULATION_PACKS,
    REFERENCE_DIR,
    UPLOAD_DIR,
    default_settings,
)
from dna_compare.eigenstrat import AadrPanel
from dna_compare.report import write_report
from dna_compare.service import AnalysisService


def _print_downloads() -> None:
    print("AADR v66.1 Human Origins lives in data/references/aadr/ and is enough for")
    print("Greek, Chinese, Persian/Iranian, Indian-caste, and hominin comparisons.\n")
    print("Optional extra files (already extracted where possible):")
    print("Hominin (data/references/hominin/)")
    for spec in HOMININ_DOWNLOADS.values():
        print(f"  - {spec['filename']}")
    print("Indian extras (data/references/caste/) — not required if AADR is present")
    for spec in CASTE_DOWNLOADS.values():
        print(f"  - {spec['filename']}")
    print(f"\nPopulation packs: {', '.join(POPULATION_PACKS)}")
    print(f"Root: {REFERENCE_DIR}")


def _list_pops(query: str | None) -> None:
    panel = AadrPanel(default_settings())
    if not panel.available:
        print("AADR .ind not found. Expected data/references/aadr/v66.p1_HO.aadr.patch.PUB.ind")
        return
    counts: dict[str, int] = {}
    for rec in panel.inds():
        counts[rec.population] = counts.get(rec.population, 0) + 1
    q = (query or "").lower()
    rows = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    if q:
        rows = [(p, n) for p, n in rows if q in p.lower()]
    print(f"{len(rows)} AADR groups" + (f" matching {query!r}" if query else ""))
    for pop, n in rows[:200]:
        print(f"  {n:4d}  {pop}")
    if len(rows) > 200:
        print(f"  ... {len(rows) - 200} more")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compare a SNP VCF to AADR HO populations (Greek, Chinese, Persian, caste) and hominins."
    )
    sub = parser.add_subparsers(dest="command", required=False)
    analyze = sub.add_parser("analyze", help="Parse a VCF and open an HTML dashboard of the results")
    analyze.add_argument("vcf", type=Path)
    analyze.add_argument(
        "--other",
        type=Path,
        default=None,
        help="Optional second SNP VCF (parent, cousin, friend, …) for a kinship / closeness estimate",
    )
    analyze.add_argument("--json", action="store_true", help="Print the API-shaped JSON payload")
    analyze.add_argument("--text", action="store_true", help="Print a text summary instead of HTML")
    analyze.add_argument("--html", type=Path, default=None, help="Write the HTML report to this path")
    analyze.add_argument("--no-open", action="store_true", help="Do not open the HTML report in a browser")
    analyze.add_argument("--no-hominin", action="store_true")
    analyze.add_argument("--no-caste", action="store_true")
    analyze.add_argument("--no-populations", action="store_true")
    analyze.add_argument("--no-ancestry", action="store_true")
    analyze.add_argument("--no-haplogroups", action="store_true")
    analyze.add_argument("--no-disease", action="store_true")
    analyze.add_argument(
        "--groups",
        default="greek,chinese,persian,caste",
        help="Comma-separated packs: greek,chinese,persian,caste",
    )
    analyze.add_argument(
        "--pops",
        default="",
        help="Extra exact AADR group IDs, comma-separated (see list-pops)",
    )
    analyze.add_argument("--preview", type=int, default=None, help="How many SNPs to include for the HTML table")
    analyze.add_argument(
        "--assembly",
        choices=("auto", "GRCh37", "GRCh38"),
        default="auto",
        help="Query VCF assembly. auto reads the header (GSA-24v3 / gtc2vcf defaults to GRCh38). "
        "GRCh38 files are lifted to hg19 before AADR overlap.",
    )
    sub.add_parser("list-references", help="Print reference layout")
    pops = sub.add_parser("list-pops", help="Search AADR population labels")
    pops.add_argument("query", nargs="?", default="")
    sub.add_parser("fetch-references", help="Rebuild compact optional extracts")
    cache_cmd = sub.add_parser(
        "build-freq-cache",
        help="Pre-build AADR population frequency caches (one TGENO pass; run after placing AADR files)",
    )
    cache_cmd.add_argument(
        "--force",
        action="store_true",
        help="Delete existing aadr_freq_*.npy caches and rebuild",
    )
    cache_cmd.add_argument(
        "--groups",
        default="greek,chinese,persian,caste",
        help="Comma-separated population packs (same as analyze --groups)",
    )
    serve_cmd = sub.add_parser("serve", help="Open the HTML dashboard (upload or bundled samples)")
    serve_cmd.add_argument("--host", default="127.0.0.1")
    serve_cmd.add_argument("--port", type=int, default=8765)
    serve_cmd.add_argument("--no-open", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv if argv is not None else (sys.argv[1:] or ["serve"]))
    if args.command in (None, "serve"):
        from dna_compare.web import serve

        host = getattr(args, "host", "127.0.0.1")
        port = getattr(args, "port", 8765)
        open_tab = not getattr(args, "no_open", False)
        serve(host=host, port=port, open_browser_tab=open_tab)
        return 0
    if args.command == "list-references":
        _print_downloads()
        return 0
    if args.command == "list-pops":
        _list_pops(args.query or None)
        return 0
    if args.command == "fetch-references":
        from dna_compare.fetch_references import fetch_references

        for line in fetch_references():
            print(line)
        return 0
    if args.command == "build-freq-cache":
        settings = default_settings()
        packs = tuple(p.strip() for p in args.groups.split(",") if p.strip())
        settings.population_packs = packs
        panel = AadrPanel(settings)
        if not panel.available:
            print("AADR HO files not found under data/references/aadr/")
            return 1
        built, skipped = panel.build_all_population_freq_caches(force=args.force)
        print(
            f"Frequency caches: built {built}, already present {skipped} "
            f"(dir: {settings.cache_dir})"
        )
        return 0

    settings = default_settings()
    settings.auto_download_pgs = True
    if args.preview is not None:
        settings.variant_preview_limit = args.preview
    if args.assembly != "auto":
        settings.assume_assembly = args.assembly
    packs = tuple(p.strip() for p in args.groups.split(",") if p.strip())
    if args.no_caste:
        packs = tuple(p for p in packs if p != "caste")
    settings.population_packs = packs
    extra = {}
    if args.pops:
        for name in args.pops.split(","):
            name = name.strip()
            if name:
                extra[name] = (name,)
        settings.extra_pops = extra
    service = AnalysisService(settings)
    result = service.analyze(
        args.vcf,
        compare_hominin_flag=not args.no_hominin,
        compare_caste_flag=not args.no_caste,
        compare_populations_flag=not args.no_populations,
        compare_ancestry_flag=not args.no_ancestry,
        compare_haplogroups_flag=not args.no_haplogroups,
        compare_disease_flag=not args.no_disease,
        other_source=args.other,
        other_filename=None if args.other is None else args.other.name,
    )
    payload = result.to_dict()
    if args.json:
        print(json.dumps(payload, indent=2))
        return 0 if result.ok else 1
    if args.text:
        _print_human(payload)
        return 0 if result.ok else 1
    dest = args.html or (UPLOAD_DIR / f"{args.vcf.stem}.report.html")
    write_report(payload, dest)
    print(dest)
    if not args.no_open:
        open_browser(dest.resolve().as_uri())
    return 0 if result.ok else 1


def _print_human(payload: dict) -> None:
    print(f"File: {payload['source_filename']}  ok={payload['ok']}")
    if payload["errors"]:
        for err in payload["errors"]:
            print(f"  error: {err}")
    vcf = payload.get("vcf") or {}
    if vcf:
        print(
            f"Sample {vcf.get('sample_id')}  SNPs={vcf.get('n_snps')}  "
            f"skipped_non_snp={vcf.get('n_non_snp_skipped')}  "
            f"assembly={vcf.get('assembly')}"
            + (f"→{vcf.get('lifted_to')} lifted={vcf.get('n_lifted')}" if vcf.get("lifted_to") else "")
        )
        preview = vcf.get("preview") or []
        if preview:
            print(f"Preview ({len(preview)} rows):")
            for row in preview[:10]:
                print(
                    f"  {row['chrom']}:{row['pos']}  {row['rsid']}  "
                    f"{row['ref']}>{row['alt']}  GT={row['genotype']}"
                )
    for kind in ("hominin", "ancestry", "ancestry_5", "community_ref", "disease", "additional", "populations", "caste"):
        block = payload.get(kind) or {}
        if block.get("hidden"):
            continue
        print(f"\n{kind} comparison  available={block['available']}")
        for note in block.get("notes") or []:
            print(f"  note: {note}")
        for est in block.get("estimates") or []:
            if "finding" in est:
                print(
                    f"  {est.get('topic')}: {est.get('finding')}  "
                    f"snps={est.get('n_snps')}/{est.get('n_markers')}"
                )
            elif "trait" in est:
                print(
                    f"  {est['trait']} ({est.get('pgs_id')}): "
                    f"score={est.get('score')}  "
                    f"snps={est.get('n_snps')}/{est.get('n_score')}  "
                    f"coverage={est.get('coverage_pct')}%"
                )
            elif "population" in est:
                print(
                    f"  {est['population']}: {est['percent']}%  "
                    f"snps={est['n_snps']}  ibs={est.get('mean_ibs')}"
                )
            else:
                print(
                    f"  {est['label']}: {est['percent']}%  "
                    f"method={est['method']}  snps={est['n_snps']}"
                )
    related = payload.get("relatedness") or {}
    print(
        f"\nrelatedness  available={related.get('available')}  "
        f"other={related.get('other_filename')}  "
        f"kinship={related.get('kinship')}  "
        f"call={related.get('relationship')}  "
        f"reliability={related.get('reliability')}"
    )
    if related.get("available"):
        print(
            f"  shared={related.get('shared_pct')}%  "
            f"IBS0={related.get('ibs0_pct')}%  "
            f"IBS1={related.get('ibs1_pct')}%  "
            f"IBS2={related.get('ibs2_pct')}%"
        )
    for note in related.get("notes") or []:
        print(f"  note: {note}")
    haplo = payload.get("haplogroups") or {}
    print(f"\nhaplogroups  available={haplo.get('available')}  best={haplo.get('sample_best')}")
    for note in haplo.get("notes") or []:
        print(f"  note: {note}")
    if haplo.get("available"):
        for row in haplo.get("rows") or []:
            bits = []
            for name, cell in (row.get("groups") or {}).items():
                n_called = cell.get("n_called") or 0
                if not n_called:
                    continue
                pct = cell.get("percent")
                extra = f" {pct:.0f}%" if pct is not None else ""
                bits.append(f"{name} {cell.get('n')}/{n_called}{extra}")
            print(f"  {row.get('haplogroup')} [{row.get('sample_status')}]: " + "; ".join(bits[:8]))
    print(f"\nmtDNA haplogroups  available={haplo.get('mt_available')}  best={haplo.get('mt_sample_best')}")
    for note in haplo.get("mt_notes") or []:
        print(f"  note: {note}")
    if haplo.get("mt_available"):
        for row in haplo.get("mt_rows") or []:
            bits = []
            for name, cell in (row.get("groups") or {}).items():
                n_called = cell.get("n_called") or 0
                if not n_called:
                    continue
                pct = cell.get("percent")
                extra = f" {pct:.0f}%" if pct is not None else ""
                bits.append(f"{name} {cell.get('n')}/{n_called}{extra}")
            print(f"  {row.get('haplogroup')} [{row.get('sample_status')}]: " + "; ".join(bits[:8]))
    if haplo.get("status_notes"):
        print("\nhaplogroup status notes")
        for note in haplo.get("status_notes") or []:
            print(f"  note: {note}")


if __name__ == "__main__":
    raise SystemExit(main())
