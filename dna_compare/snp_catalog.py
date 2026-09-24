from __future__ import annotations

import re
from pathlib import Path

from dna_compare.config import Settings
from dna_compare.vcf_parser import chrom_sort_key, preview_rows

CATALOG_SUFFIX = ".snp_catalog.tsv"
_SAFE_NAME = re.compile(r"[^A-Za-z0-9._-]+")


def catalog_id_for(filename: str) -> str:
    safe = _SAFE_NAME.sub("_", Path(filename).name).strip("._") or "upload.vcf"
    return safe + CATALOG_SUFFIX


def catalog_path(filename: str, settings: Settings) -> Path:
    return settings.cache_dir / "snp_lists" / catalog_id_for(filename)


def _include_row(row: dict) -> bool:
    return bool(row.get("is_snp")) or row.get("chrom") in {"Y", "MT"}


def gt_carries_alt(gt: str, dosage=None) -> bool:
    """True when GT has a non-reference allele (0/1, 1/0, 1/1, haploid 1). 0/0 is false."""
    token = (gt or "").split(":", 1)[0].strip()
    if not token or token in {".", "./.", ".|."}:
        if dosage is None or dosage == "":
            return False
        try:
            return float(dosage) > 0
        except (TypeError, ValueError):
            return False
    for part in token.replace("|", "/").split("/"):
        if part.isdigit() and int(part) > 0:
            return True
    if dosage not in {None, ""}:
        try:
            return float(dosage) > 0
        except (TypeError, ValueError):
            return False
    return False


def write_snp_catalog(index: dict[tuple[str, int], dict], path: Path) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = [row for row in index.values() if _include_row(row)]
    rows.sort(key=lambda row: (chrom_sort_key(str(row["chrom"])), int(row["pos"])))
    with path.open("w", encoding="utf-8") as handle:
        handle.write("chrom\tpos\trsid\tref\talt\tgt\tdosage\n")
        for row in rows:
            dosage = row.get("dosage_alt")
            dosage_s = "" if dosage is None else str(dosage)
            handle.write(
                "\t".join(
                    [
                        str(row.get("chrom") or ""),
                        str(row.get("pos") or ""),
                        str(row.get("rsid") or ""),
                        str(row.get("ref") or ""),
                        str(row.get("alt") or ""),
                        str(row.get("genotype") or ""),
                        dosage_s,
                    ]
                )
                + "\n"
            )
    return len(rows)


def query_snp_catalog(
    path: Path,
    *,
    chrom: str | None = None,
    q: str | None = None,
    called_only: bool = False,
    offset: int = 0,
    limit: int = 200,
) -> dict:
    if not path.is_file():
        return {"ok": False, "error": "SNP list not found. Analyze a VCF first.", "rows": [], "total": 0}
    want = (chrom or "").strip()
    if want.lower() in {"", "all"}:
        want = ""
    needle = (q or "").strip().lower()
    offset = max(0, int(offset))
    limit = max(1, min(int(limit), 1000))
    matched: list[dict] = []
    total = 0
    with path.open(encoding="utf-8") as handle:
        header = handle.readline()
        if not header:
            return {"ok": True, "rows": [], "total": 0, "offset": offset, "limit": limit}
        for line in handle:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6:
                continue
            if want and cols[0] != want:
                continue
            if needle and needle not in line.lower():
                continue
            if called_only and not gt_carries_alt(cols[5], cols[6] if len(cols) > 6 else None):
                continue
            if offset <= total < offset + limit:
                matched.append(
                    {
                        "chrom": cols[0],
                        "pos": int(cols[1]) if cols[1].isdigit() else cols[1],
                        "rsid": cols[2],
                        "ref": cols[3],
                        "alt": cols[4],
                        "genotype": cols[5],
                        "dosage_alt": None if len(cols) < 7 or cols[6] == "" else float(cols[6]),
                    }
                )
            total += 1
    return {
        "ok": True,
        "rows": matched,
        "total": total,
        "offset": offset,
        "limit": limit,
        "chrom": want or None,
        "q": needle or None,
        "called_only": bool(called_only),
    }


def filter_catalog_tsv(path: Path, *, called_only: bool = False, chrom: str | None = None) -> bytes:
    if not path.is_file():
        return b""
    want = (chrom or "").strip()
    if want.lower() in {"", "all"}:
        want = ""
    out: list[str] = []
    with path.open(encoding="utf-8") as handle:
        header = handle.readline()
        if not header:
            return b""
        out.append(header if header.endswith("\n") else header + "\n")
        for line in handle:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6:
                continue
            if want and cols[0] != want:
                continue
            if called_only and not gt_carries_alt(cols[5], cols[6] if len(cols) > 6 else None):
                continue
            out.append(line if line.endswith("\n") else line + "\n")
    return "".join(out).encode("utf-8")


def publish_snp_list(
    summary,
    index: dict[tuple[str, int], dict],
    *,
    filename: str,
    settings: Settings,
) -> Path:
    path = catalog_path(filename, settings)
    n = write_snp_catalog(index, path)
    summary.snp_catalog_id = path.name
    summary.n_catalog = n
    embed_limit = getattr(settings, "variant_embed_limit", 5000)
    if n <= embed_limit:
        summary.preview = preview_rows(index, limit=None)
    else:
        summary.preview = preview_rows(index, limit=settings.variant_preview_limit, stratified=True)
    return path
