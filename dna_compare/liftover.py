from __future__ import annotations

import gzip
import subprocess
import urllib.request
from bisect import bisect_right
from dataclasses import dataclass
from pathlib import Path

from dna_compare.config import HG38_TO_HG19_CHAIN, HG38_TO_HG19_CHAIN_URL
from dna_compare.vcf_parser import normalize_chrom

_COMP = str.maketrans("ACGTacgt", "TGCAtgca")
_SKIP = {"MT", "M"}


@dataclass(frozen=True)
class _Block:
    t_start: int
    t_end: int
    q_chrom: str
    q_start: int
    q_strand: str
    q_size: int


class ChainMap:
    """Minimal UCSC chain reader for SNP coordinate conversion (1-based VCF positions)."""

    def __init__(self, blocks: dict[str, list[_Block]]):
        self._blocks = {chrom: tuple(sorted(items, key=lambda b: b.t_start)) for chrom, items in blocks.items()}
        self._starts = {chrom: [item.t_start for item in items] for chrom, items in self._blocks.items()}

    def convert(self, chrom: str, pos_1based: int) -> tuple[str, int, str] | None:
        chrom = normalize_chrom(chrom)
        if chrom in _SKIP:
            return chrom, pos_1based, "+"
        blocks = self._blocks.get(chrom)
        if not blocks:
            return None
        pos0 = pos_1based - 1
        idx = bisect_right(self._starts[chrom], pos0) - 1
        if idx < 0:
            return None
        block = blocks[idx]
        if not (block.t_start <= pos0 < block.t_end):
            return None
        offset = pos0 - block.t_start
        if block.q_strand == "+":
            q0 = block.q_start + offset
            strand = "+"
        else:
            q0 = block.q_size - (block.q_start + offset) - 1
            strand = "-"
        return block.q_chrom, q0 + 1, strand


def _open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open("rt")


def load_chain(path: Path) -> ChainMap:
    by_chrom: dict[str, list[_Block]] = {}
    with _open_text(path) as handle:
        header: list[str] | None = None
        align_rows: list[str] = []
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                if header and align_rows:
                    _flush_chain(header, align_rows, by_chrom)
                header = None
                align_rows = []
                continue
            if line.startswith("chain"):
                if header and align_rows:
                    _flush_chain(header, align_rows, by_chrom)
                header = line.split()
                align_rows = []
                continue
            if header:
                align_rows.append(line)
        if header and align_rows:
            _flush_chain(header, align_rows, by_chrom)
    if not by_chrom:
        raise ValueError(f"No chain intervals in {path}")
    return ChainMap(by_chrom)


def _flush_chain(fields: list[str], rows: list[str], by_chrom: dict[str, list[_Block]]) -> None:
    if len(fields) < 12:
        return
    t_chrom = normalize_chrom(fields[2])
    q_chrom = normalize_chrom(fields[7])
    q_size = int(fields[8])
    q_strand = fields[9]
    t_pos = int(fields[5])
    q_pos = int(fields[10])
    for row in rows:
        parts = row.split()
        size = int(parts[0])
        by_chrom.setdefault(t_chrom, []).append(
            _Block(
                t_start=t_pos,
                t_end=t_pos + size,
                q_chrom=q_chrom,
                q_start=q_pos,
                q_strand=q_strand,
                q_size=q_size,
            )
        )
        if len(parts) >= 3:
            t_pos += size + int(parts[1])
            q_pos += size + int(parts[2])


def _revcomp(seq: str) -> str:
    if not seq or seq in {".", "*", "D", "I"}:
        return seq
    return seq.translate(_COMP)[::-1]


def lift_query_index(
    index: dict[tuple[str, int], dict],
    chain: ChainMap,
) -> tuple[dict[tuple[str, int], dict], dict[str, int]]:
    lifted: dict[tuple[str, int], dict] = {}
    leftover: list[tuple[tuple[str, int], dict]] = []
    n_lifted = n_unmapped = n_kept = n_flipped = 0
    for (chrom, pos), row in index.items():
        mapped = chain.convert(chrom, pos)
        if mapped is None:
            leftover.append(((chrom, pos), row))
            n_unmapped += 1
            continue
        new_chrom, new_pos, strand = mapped
        if new_chrom == chrom and new_pos == pos and chrom in _SKIP:
            lifted[(chrom, pos)] = row
            n_kept += 1
            continue
        item = dict(row)
        item["chrom"] = new_chrom
        item["pos"] = new_pos
        if strand == "-":
            item["ref"] = _revcomp(str(item.get("ref") or ""))
            item["alt"] = _revcomp(str(item.get("alt") or ""))
            n_flipped += 1
        lifted[(new_chrom, new_pos)] = item
        n_lifted += 1
    for key, row in leftover:
        lifted.setdefault(key, row)
    return lifted, {
        "n_lifted": n_lifted,
        "n_unmapped": n_unmapped,
        "n_kept": n_kept,
        "n_flipped": n_flipped,
    }


def ensure_hg38_to_hg19_chain(
    path: Path | None = None,
    *,
    url: str = HG38_TO_HG19_CHAIN_URL,
    download: bool = True,
) -> Path | None:
    dest = Path(path) if path is not None else HG38_TO_HG19_CHAIN
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    if not download:
        return None
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    try:
        urllib.request.urlretrieve(url, tmp)
        tmp.replace(dest)
    except OSError:
        if tmp.exists():
            tmp.unlink()
        if not _curl_download(url, dest):
            return None
    return dest if dest.exists() and dest.stat().st_size > 0 else None


def _curl_download(url: str, dest: Path) -> bool:
    tmp = dest.with_suffix(dest.suffix + ".part")
    try:
        completed = subprocess.run(
            ["curl", "-fsSL", "--retry", "2", "-o", str(tmp), url],
            check=False,
            capture_output=True,
        )
    except OSError:
        return False
    if completed.returncode != 0 or not tmp.exists() or tmp.stat().st_size == 0:
        if tmp.exists():
            tmp.unlink()
        return False
    tmp.replace(dest)
    return True
