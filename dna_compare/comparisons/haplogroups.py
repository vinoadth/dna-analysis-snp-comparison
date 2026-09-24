from __future__ import annotations

import re
from collections import Counter
from pathlib import Path

from dna_compare.config import Settings
from dna_compare.models import (
    HaplogroupGroupCount,
    HaplogroupMarkerCall,
    HaplogroupResult,
    HaplogroupRow,
)

# Published backbone Y markers. `positions` are hg19 / GRCh37; `positions_hg38` are GRCh38.
# Alleles from ISOGG / FTDNA / Karafet. M17 is an indel and is ignored if it contradicts M173.
Y_MARKERS: tuple[dict, ...] = (
    {"haplogroup": "R (M207)", "marker": "M207", "rsids": ("rs2032658",), "positions": (15581983,), "positions_hg38": (13470103,), "ancestral": "A", "derived": "G", "backbone": True},
    {"haplogroup": "R1 (M173)", "marker": "M173", "rsids": ("rs2032624",), "positions": (15026424,), "positions_hg38": (12914512,), "ancestral": "A", "derived": "C", "backbone": True},
    {"haplogroup": "R1a1 (M17)", "marker": "M17", "rsids": ("rs3908",), "positions": (21733165, 21733168), "positions_hg38": (19571279, 19571280, 19571281, 19571282), "ancestral": "D", "derived": "I", "backbone": True},
    {"haplogroup": "R1a1a (M198)", "marker": "M198", "rsids": ("rs2020857",), "positions": (15030752,), "positions_hg38": (12918840,), "ancestral": "C", "derived": "T", "backbone": True},
    {"haplogroup": "R1a-Z93", "marker": "Z93", "rsids": ("Z93",), "positions": (7552356,), "positions_hg38": (7684315,), "ancestral": "G", "derived": "A", "backbone": True},
    {"haplogroup": "R1b (M343)", "marker": "M343", "rsids": ("rs9786184",), "positions": (2887824,), "positions_hg38": (3019783,), "ancestral": "C", "derived": "A", "backbone": True},
    {"haplogroup": "R2 (M124)", "marker": "M124", "rsids": ("M124",), "positions": (21764501,), "positions_hg38": (19602615,), "ancestral": "G", "derived": "A", "backbone": True},
    {"haplogroup": "H (M69)", "marker": "M69", "rsids": ("rs2032673",), "positions": (21894058,), "positions_hg38": (19732172,), "ancestral": "T", "derived": "C", "backbone": True},
    {"haplogroup": "H1 (M52)", "marker": "M52", "rsids": ("rs376769460",), "positions": (21753199,), "positions_hg38": (19591313,), "ancestral": "A", "derived": "C", "backbone": True},
    {"haplogroup": "H3 (Z5857)", "marker": "Z5857", "rsids": ("rs569006329",), "positions": (2759285,), "positions_hg38": (2891244,), "ancestral": "C", "derived": "G", "backbone": True},
    {"haplogroup": "H3b (Z13871)", "marker": "Z13871", "rsids": ("Z13871",), "positions": (2878605,), "positions_hg38": (3010564,), "ancestral": "T", "derived": "G", "backbone": True},
    {"haplogroup": "L (M20)", "marker": "M20", "rsids": ("rs3911",), "positions": (21733454,), "positions_hg38": (19571568,), "ancestral": "A", "derived": "G", "backbone": True},
    {"haplogroup": "L1 (M27)", "marker": "M27", "rsids": ("M27",), "positions": (21739646,), "positions_hg38": (19577760,), "ancestral": "C", "derived": "G", "backbone": True},
    {"haplogroup": "L3 (M357)", "marker": "M357", "rsids": ("M357",), "positions": (2888252,), "positions_hg38": (3020211,), "ancestral": "C", "derived": "A", "backbone": True},
    {"haplogroup": "J2 (M172)", "marker": "M172", "rsids": ("rs2032604",), "positions": (14969634,), "positions_hg38": (12857709,), "ancestral": "T", "derived": "G", "backbone": True},
    {"haplogroup": "J2a (M410)", "marker": "M410", "rsids": ("M410",), "positions": (2751678,), "positions_hg38": (2883637,), "ancestral": "A", "derived": "G", "backbone": True},
    {"haplogroup": "J2b (M241)", "marker": "M241", "rsids": ("M241",), "positions": (15018459,), "positions_hg38": (12906548,), "ancestral": "G", "derived": "A", "backbone": True},
    {"haplogroup": "J1 (M267)", "marker": "M267", "rsids": ("rs9341313",), "positions": (22741818,), "positions_hg38": (20579932,), "ancestral": "T", "derived": "G", "backbone": True},
    {"haplogroup": "I (M170)", "marker": "M170", "rsids": ("rs2032597",), "positions": (14847792,), "positions_hg38": (12735858,), "ancestral": "A", "derived": "C", "backbone": True},
    {"haplogroup": "G (M201)", "marker": "G-M201", "rsids": ("rs2032636",), "positions": (15027529,), "positions_hg38": (12915617,), "ancestral": "G", "derived": "T", "backbone": True},
    {"haplogroup": "C (M130)", "marker": "M130", "rsids": ("rs35284970",), "positions": (2734854,), "positions_hg38": (2866813,), "ancestral": "C", "derived": "T", "backbone": True},
    {"haplogroup": "Q (M242)", "marker": "M242", "rsids": ("rs8179021",), "positions": (15018582,), "positions_hg38": (12906671,), "ancestral": "C", "derived": "T", "backbone": True},
    {"haplogroup": "T (M70)", "marker": "M70", "rsids": ("M70",), "positions": (21893881,), "positions_hg38": (19731995,), "ancestral": "A", "derived": "C", "backbone": True},
    {"haplogroup": "E (M96)", "marker": "M96", "rsids": ("M96",), "positions": (21778998,), "positions_hg38": (19617112,), "ancestral": "C", "derived": "G", "backbone": True},
)

# Deepest-first labels used for the sample call.
_SPECIFICITY = (
    "R1a-Z93",
    "R1a1a (M198)",
    "R1a1 (M17)",
    "R1b (M343)",
    "R2 (M124)",
    "R1 (M173)",
    "R (M207)",
    "J2b (M241)",
    "J2a (M410)",
    "J2 (M172)",
    "J1 (M267)",
    "H3b (Z13871)",
    "H3 (Z5857)",
    "H1 (M52)",
    "H (M69)",
    "L3 (M357)",
    "L1 (M27)",
    "L (M20)",
    "G (M201)",
    "C (M130)",
    "Q (M242)",
    "T (M70)",
    "E (M96)",
    "I (M170)",
)

# Rows in the AADR frequency table (ISOGG prefix → display).
# Third field is the sample-marker haplogroup key when it differs from the row label.
_FREQ_ROWS: tuple[tuple[str, str, str], ...] = (
    ("R1a1 (M17)", "M17", "R1a1 (M17)"),
    ("R1a-Z93", "Z93", "R1a-Z93"),
    ("R1b (M343)", "M343", "R1b (M343)"),
    ("R2 (M124)", "M124", "R2 (M124)"),
    ("H (M69)", "M69", "H (M69)"),
    ("H1 (M52)", "M52", "H1 (M52)"),
    ("H3 (Z5857)", "Z5857", "H3 (Z5857)"),
    ("H3b (Z13871)", "Z13871", "H3b (Z13871)"),
    ("L (M20)", "M20", "L (M20)"),
    ("L1 (M27)", "M27", "L1 (M27)"),
    ("L3 (M357)", "M357", "L3 (M357)"),
    ("J2 (M172)", "M172", "J2 (M172)"),
    ("J2a (M410)", "M410", "J2a (M410)"),
    ("J2b (M241)", "M241", "J2b (M241)"),
    ("J1 (M267)", "M267", "J1 (M267)"),
    ("C", "M130", "C (M130)"),
    ("O", "M175", ""),
    ("G", "M201", "G (M201)"),
    ("Q", "M242", "Q (M242)"),
    ("T (M70)", "M70", "T (M70)"),
    ("E (M96)", "M96", "E (M96)"),
    ("I (M170)", "M170", "I (M170)"),
    ("other", "", ""),
)

def _mt_digit_excludes(prefix: str) -> tuple[str, ...]:
    """M4 should not count M40; U1 should not count U10."""
    return tuple(f"{prefix}{d}" for d in range(10))


# rCRS / PhyloTree Build 17 positions (rCRS is H2a2a1, so REF=rCRS is often already derived).
# Coding SNPs only; skip recurrent sites (G1888A = M5 and T; T16362C; T13500C = U7 and R2).
MT_MARKERS: tuple[dict, ...] = (
    {"haplogroup": "M (10400T)", "marker": "10400T", "chrom": "MT", "rsids": ("rs28358278",), "positions": (10400,), "ancestral": "C", "derived": "T", "backbone": True},
    {"haplogroup": "M2 (1780C)", "marker": "1780C", "chrom": "MT", "rsids": ("rs2854127",), "positions": (1780,), "ancestral": "T", "derived": "C", "backbone": True},
    {"haplogroup": "M3 (4580C)", "marker": "4580C", "chrom": "MT", "rsids": (), "positions": (4580,), "ancestral": "T", "derived": "C", "backbone": True},
    {"haplogroup": "M4 (6620C)", "marker": "6620C", "chrom": "MT", "rsids": (), "positions": (6620,), "ancestral": "T", "derived": "C", "backbone": True},
    {"haplogroup": "M5 (12477C)", "marker": "12477C", "chrom": "MT", "rsids": (), "positions": (12477,), "ancestral": "T", "derived": "C", "backbone": True},
    {"haplogroup": "M6 (10640C)", "marker": "10640C", "chrom": "MT", "rsids": (), "positions": (10640,), "ancestral": "T", "derived": "C", "backbone": True},
    {"haplogroup": "M30 (15431A)", "marker": "15431A", "chrom": "MT", "rsids": (), "positions": (15431,), "ancestral": "G", "derived": "A", "backbone": True},
    {"haplogroup": "M33 (2361A)", "marker": "2361A", "chrom": "MT", "rsids": (), "positions": (2361,), "ancestral": "G", "derived": "A", "backbone": True},
    {"haplogroup": "M36 (7271G)", "marker": "7271G", "chrom": "MT", "rsids": (), "positions": (7271,), "ancestral": "A", "derived": "G", "backbone": True},
    {"haplogroup": "R (12705C)", "marker": "12705C", "chrom": "MT", "rsids": (), "positions": (12705,), "ancestral": "T", "derived": "C", "backbone": True},
    {"haplogroup": "R5 (8594C)", "marker": "8594C", "chrom": "MT", "rsids": ("rs372365688",), "positions": (8594,), "ancestral": "T", "derived": "C", "backbone": True},
    {"haplogroup": "R6 (12285C)", "marker": "12285C", "chrom": "MT", "rsids": (), "positions": (12285,), "ancestral": "T", "derived": "C", "backbone": True},
    {"haplogroup": "R7 (13105G)", "marker": "13105G", "chrom": "MT", "rsids": (), "positions": (13105,), "ancestral": "A", "derived": "G", "backbone": True},
    {"haplogroup": "R8 (13215C)", "marker": "13215C", "chrom": "MT", "rsids": (), "positions": (13215,), "ancestral": "T", "derived": "C", "backbone": True},
    {"haplogroup": "R30 (8584A)", "marker": "8584A", "chrom": "MT", "rsids": ("rs3135028",), "positions": (8584,), "ancestral": "G", "derived": "A", "backbone": True},
    {"haplogroup": "U (12308G)", "marker": "12308G", "chrom": "MT", "rsids": (), "positions": (12308,), "ancestral": "A", "derived": "G", "backbone": True},
    {"haplogroup": "U1 (12879C)", "marker": "12879C", "chrom": "MT", "rsids": (), "positions": (12879,), "ancestral": "T", "derived": "C", "backbone": True},
    {"haplogroup": "U2 (16051G)", "marker": "16051G", "chrom": "MT", "rsids": ("rs117565943",), "positions": (16051,), "ancestral": "A", "derived": "G", "backbone": True},
    {"haplogroup": "U7 (5360T)", "marker": "5360T", "chrom": "MT", "rsids": (), "positions": (5360,), "ancestral": "C", "derived": "T", "backbone": True},
    {"haplogroup": "HV (14766C)", "marker": "14766C", "chrom": "MT", "rsids": ("rs3135031",), "positions": (14766,), "ancestral": "T", "derived": "C", "backbone": True},
    {"haplogroup": "H (7028C)", "marker": "7028C", "chrom": "MT", "rsids": (), "positions": (7028,), "ancestral": "T", "derived": "C", "backbone": True},
    {"haplogroup": "H3 (6776C)", "marker": "6776C", "chrom": "MT", "rsids": (), "positions": (6776,), "ancestral": "T", "derived": "C", "backbone": True},
    {"haplogroup": "H3b (2581G)", "marker": "2581G", "chrom": "MT", "rsids": (), "positions": (2581,), "ancestral": "A", "derived": "G", "backbone": True},
    {"haplogroup": "JT (16126C)", "marker": "16126C", "chrom": "MT", "rsids": ("rs147029798",), "positions": (16126,), "ancestral": "T", "derived": "C", "backbone": True},
    {"haplogroup": "J (16069T)", "marker": "16069T", "chrom": "MT", "rsids": (), "positions": (16069,), "ancestral": "C", "derived": "T", "backbone": True},
    {"haplogroup": "T (14905A)", "marker": "14905A", "chrom": "MT", "rsids": (), "positions": (14905,), "ancestral": "G", "derived": "A", "backbone": True},
    {"haplogroup": "N1 (10238C)", "marker": "10238C", "chrom": "MT", "rsids": (), "positions": (10238,), "ancestral": "T", "derived": "C", "backbone": True},
    {"haplogroup": "W (3505G)", "marker": "3505G", "chrom": "MT", "rsids": (), "positions": (3505,), "ancestral": "A", "derived": "G", "backbone": True},
)

_MT_SPECIFICITY = (
    "H3b (2581G)",
    "H3 (6776C)",
    "H (7028C)",
    "HV (14766C)",
    "U7 (5360T)",
    "U2 (16051G)",
    "U1 (12879C)",
    "U (12308G)",
    "J (16069T)",
    "T (14905A)",
    "JT (16126C)",
    "R8 (13215C)",
    "R7 (13105G)",
    "R6 (12285C)",
    "R5 (8594C)",
    "R30 (8584A)",
    "W (3505G)",
    "N1 (10238C)",
    "M36 (7271G)",
    "M33 (2361A)",
    "M30 (15431A)",
    "M6 (10640C)",
    "M5 (12477C)",
    "M4 (6620C)",
    "M3 (4580C)",
    "M2 (1780C)",
    "M (10400T)",
    "R (12705C)",
)

# display, marker, sample key, AADR prefix, prefixes to exclude
_MT_FREQ_ROWS: tuple[tuple[str, str, str, str, tuple[str, ...]], ...] = (
    ("M (10400T)", "10400T", "M (10400T)", "M", ()),
    ("M2 (1780C)", "1780C", "M2 (1780C)", "M2", _mt_digit_excludes("M2")),
    ("M3", "4580C", "M3 (4580C)", "M3", ("M30", "M33", "M36")),
    ("M4 (6620C)", "6620C", "M4 (6620C)", "M4", _mt_digit_excludes("M4")),
    ("M5 (12477C)", "12477C", "M5 (12477C)", "M5", _mt_digit_excludes("M5")),
    ("M6 (10640C)", "10640C", "M6 (10640C)", "M6", _mt_digit_excludes("M6")),
    ("M30 (15431A)", "15431A", "M30 (15431A)", "M30", ()),
    ("M33 (2361A)", "2361A", "M33 (2361A)", "M33", ()),
    ("M36 (7271G)", "7271G", "M36 (7271G)", "M36", ()),
    ("R (12705C)", "12705C", "R (12705C)", "R", ()),
    ("R5 (8594C)", "8594C", "R5 (8594C)", "R5", ()),
    ("R6 (12285C)", "12285C", "R6 (12285C)", "R6", _mt_digit_excludes("R6")),
    ("R7 (13105G)", "13105G", "R7 (13105G)", "R7", ()),
    ("R8 (13215C)", "13215C", "R8 (13215C)", "R8", ()),
    ("R30 (8584A)", "8584A", "R30 (8584A)", "R30", ()),
    ("U (12308G)", "12308G", "U (12308G)", "U", ()),
    ("U1 (12879C)", "12879C", "U1 (12879C)", "U1", _mt_digit_excludes("U1")),
    ("U2 (16051G)", "16051G", "U2 (16051G)", "U2", ()),
    ("U7 (5360T)", "5360T", "U7 (5360T)", "U7", ()),
    ("HV (14766C)", "14766C", "HV (14766C)", "HV", ()),
    ("H (7028C)", "7028C", "H (7028C)", "H", ("HV",)),
    ("H3 (6776C)", "6776C", "H3 (6776C)", "H3", _mt_digit_excludes("H3")),
    ("H3b (2581G)", "2581G", "H3b (2581G)", "H3b", ()),
    ("JT (16126C)", "16126C", "JT (16126C)", "JT", ()),
    ("J (16069T)", "16069T", "J (16069T)", "J", ("JT",)),
    ("T (14905A)", "14905A", "T (14905A)", "T", ()),
    ("N1 (10238C)", "10238C", "N1 (10238C)", "N1", _mt_digit_excludes("N1") + ("N21", "N22")),
    ("W (3505G)", "3505G", "W (3505G)", "W", ()),
)

_COMP = str.maketrans("ACGT", "TGCA")
_ANNO_GROUP = 14
_ANNO_SEX = 30
_ANNO_Y_ISOGG = 35
_ANNO_Y_MANUAL = 36
_ANNO_MT = 38
_RSID_RE = re.compile(r"(rs\d+)", re.I)


def _norm_rsid(rsid: str) -> str:
    text = (rsid or "").split(",")[0].strip()
    match = _RSID_RE.search(text)
    if match:
        return match.group(1).lower()
    if text.lower().startswith("exm-"):
        text = text[4:]
    return text


def _observed_allele(row: dict) -> str | None:
    gt = row.get("genotype") or ""
    ref = str(row.get("ref") or "").upper()
    alt = str(row.get("alt") or "").upper()
    if not gt or gt == "." or gt.startswith("."):
        return None
    token = gt.split(":", 1)[0].replace("|", "/")
    parts = token.split("/") if "/" in token else [token]
    alleles: list[str] = []
    for part in parts:
        if part == "0":
            alleles.append(ref)
        elif part == "1":
            alleles.append(alt)
        elif part.upper() in {ref, alt, "D", "I", "A", "C", "G", "T"}:
            alleles.append(part.upper())
        elif part == ".":
            return None
    uniq = {a for a in alleles if a}
    if len(uniq) == 1:
        return next(iter(uniq))
    if len(uniq) > 1:
        return "het"
    return None


def _rsid_index(query_index: dict[tuple[str, int], dict]) -> dict[tuple[str, str], dict]:
    out: dict[tuple[str, str], dict] = {}
    for row in query_index.values():
        rsid = _norm_rsid(str(row.get("rsid") or ""))
        if rsid and rsid not in {".", ""}:
            out.setdefault((str(row.get("chrom") or ""), rsid), row)
    return out


def _allele_compatible(row: dict, ancestral: str, derived: str) -> bool:
    observed = _observed_allele(row)
    if observed is None or observed == "het":
        return False
    obs = observed.upper()
    anc, der = ancestral.upper(), derived.upper()
    if obs in {anc, der}:
        return True
    if anc not in {"D", "I"} and der not in {"D", "I"}:
        return obs in {anc.translate(_COMP), der.translate(_COMP)}
    return False


def _lookup_marker(
    index: dict[tuple[str, int], dict],
    spec: dict,
    rsid_index: dict[tuple[str, str], dict] | None = None,
) -> dict | None:
    chrom = spec.get("chrom") or "Y"
    want = {_norm_rsid(r) for r in spec.get("rsids") or () if r}
    # rsID first: a GSA/hg38 site must not be read as a different hg19 SNP.
    if want:
        lookup = rsid_index if rsid_index is not None else _rsid_index(index)
        for rsid in want:
            row = lookup.get((chrom, rsid))
            if row:
                return row
    positions = list(spec.get("positions") or ()) + list(spec.get("positions_hg38") or ())
    hits = []
    seen: set[int] = set()
    for pos in positions:
        row = index.get((chrom, int(pos)))
        if row and id(row) not in seen:
            seen.add(id(row))
            hits.append(row)
    if not hits:
        return None
    if len(hits) == 1:
        return hits[0]
    for row in hits:
        if _allele_compatible(row, spec["ancestral"], spec["derived"]):
            return row
    return hits[0]


def _status(observed: str | None, ancestral: str, derived: str) -> str:
    if observed is None:
        return "no-call"
    if observed == "het":
        return "het"
    obs = observed.upper()
    anc, der = ancestral.upper(), derived.upper()
    if obs == der:
        return "derived"
    if obs == anc:
        return "ancestral"
    if anc not in {"D", "I"} and der not in {"D", "I"}:
        if obs == der.translate(_COMP) and anc == anc:
            return "derived"
        if obs == anc.translate(_COMP):
            return "ancestral"
    return "mismatch"


def _quality_from_row(row: dict | None) -> dict:
    if not row:
        return {}
    return {
        "qual": row.get("qual"),
        "gq": row.get("gq"),
        "dp": row.get("dp"),
        "igc": row.get("igc"),
        "ad": row.get("ad"),
    }


def _median(values: list[float]) -> float | None:
    if not values:
        return None
    xs = sorted(values)
    mid = len(xs) // 2
    if len(xs) % 2:
        return xs[mid]
    return (xs[mid - 1] + xs[mid]) / 2.0


def marker_quality_note(calls: list[HaplogroupMarkerCall], *, label: str) -> str:
    used = [call for call in calls if call.status != "no-call"]
    gqs = [float(call.gq) for call in used if call.gq is not None]
    dps = [float(call.dp) for call in used if call.dp is not None]
    igcs = [float(call.igc) for call in used if call.igc is not None]
    if not used:
        return f"No {label} marker was called, so there is no GQ or DP to show."
    if not (gqs or dps or igcs):
        return f"This VCF has genotype only on {label} markers (no GQ, DP, or IGC)."
    parts: list[str] = []
    if gqs:
        parts.append(f"median GQ {_median(gqs):.0f}")
    if dps:
        parts.append(f"median DP {_median(dps):.0f}")
    if igcs:
        parts.append(f"median IGC {_median(igcs):.2f}")
    extra = ""
    if not dps:
        extra = " No sequencing DP — typical of SNP-array files (GQ / IGC instead)."
    return f"{label} call support: {', '.join(parts)}.{extra}"


def score_markers(query_index: dict[tuple[str, int], dict], specs: tuple[dict, ...]) -> list[HaplogroupMarkerCall]:
    calls: list[HaplogroupMarkerCall] = []
    rsid_index = _rsid_index(query_index)
    for spec in specs:
        row = _lookup_marker(query_index, spec, rsid_index)
        observed = _observed_allele(row) if row else None
        status = _status(observed, spec["ancestral"], spec["derived"])
        rsids = spec.get("rsids") or ()
        calls.append(
            HaplogroupMarkerCall(
                haplogroup=spec["haplogroup"],
                marker=spec["marker"],
                rsid=rsids[0] if rsids else spec["marker"],
                chrom=spec.get("chrom") or "Y",
                pos=int(row["pos"]) if row else int(spec["positions"][0]),
                ancestral=spec["ancestral"],
                derived=spec["derived"],
                observed=observed,
                status=status,
                genotype=None if row is None else str(row.get("genotype") or ""),
                backbone=bool(spec.get("backbone", True)),
                **_quality_from_row(row),
            )
        )
    return calls


def _child_needs_parent(by_label: dict[str, HaplogroupMarkerCall], child: str, parent: str) -> None:
    c_call, p_call = by_label.get(child), by_label.get(parent)
    if c_call and p_call and c_call.status == "derived" and p_call.status in {"ancestral", "conflict"}:
        c_call.status = "conflict"


def _incompatible_if_other_derived(by_label: dict[str, HaplogroupMarkerCall], label: str, other: str) -> None:
    call, other_call = by_label.get(label), by_label.get(other)
    if call and other_call and call.status == "derived" and other_call.status == "derived":
        call.status = "conflict"


def score_y_markers(query_index: dict[tuple[str, int], dict]) -> list[HaplogroupMarkerCall]:
    calls = score_markers(query_index, Y_MARKERS)
    by_label = {call.haplogroup: call for call in calls}

    # M17 indel on chips often contradicts M173; do not let it win.
    m17 = by_label.get("R1a1 (M17)")
    m173 = by_label.get("R1 (M173)")
    if m17 and m173 and m17.status == "derived" and m173.status == "ancestral":
        m17.status = "conflict"

    _child_needs_parent(by_label, "R1a1a (M198)", "R1a1 (M17)")
    _child_needs_parent(by_label, "R1a-Z93", "R1a1a (M198)")
    _child_needs_parent(by_label, "R1a-Z93", "R1a1 (M17)")
    _child_needs_parent(by_label, "R1b (M343)", "R1 (M173)")
    _child_needs_parent(by_label, "R1 (M173)", "R (M207)")
    _child_needs_parent(by_label, "R2 (M124)", "R (M207)")
    _child_needs_parent(by_label, "H1 (M52)", "H (M69)")
    _child_needs_parent(by_label, "H3 (Z5857)", "H (M69)")
    _child_needs_parent(by_label, "H3b (Z13871)", "H3 (Z5857)")
    _child_needs_parent(by_label, "H3b (Z13871)", "H (M69)")
    _child_needs_parent(by_label, "L1 (M27)", "L (M20)")
    _child_needs_parent(by_label, "L3 (M357)", "L (M20)")
    _child_needs_parent(by_label, "J2a (M410)", "J2 (M172)")
    _child_needs_parent(by_label, "J2b (M241)", "J2 (M172)")
    return calls


def best_sample_haplogroup(
    calls: list[HaplogroupMarkerCall],
    specificity: tuple[str, ...] = _SPECIFICITY,
) -> str | None:
    derived = {c.haplogroup for c in calls if c.backbone and c.status == "derived"}
    if not derived:
        return None
    for label in specificity:
        if label in derived:
            return label
    return sorted(derived)[0]


def collapse_y_haplogroup(raw: str) -> str | None:
    text = (raw or "").strip().rstrip("~").replace(" ", "")
    if not text or text in {".", ".."} or text.lower().startswith("n/a"):
        return None
    if text.startswith("R1a"):
        if "Z93" in text or text.startswith("R1a1a1b2"):
            return "R1a-Z93"
        return "R1a1 (M17)"
    if text.startswith("R1b"):
        return "R1b (M343)"
    if text.startswith("R2"):
        return "R2 (M124)"
    if text.startswith("H"):
        if "Z13871" in text or text.startswith("H3b"):
            return "H3b (Z13871)"
        if "Z5857" in text or (text.startswith("H3") and not text.startswith("H3b")):
            return "H3 (Z5857)"
        if "M52" in text or "M82" in text or text.startswith("H1"):
            return "H1 (M52)"
        return "H (M69)"
    if text.startswith("LT"):
        return "other"
    if text.startswith("L"):
        if text.startswith("L1") or "M27" in text:
            return "L1 (M27)"
        if text.startswith("L3") or "M357" in text:
            return "L3 (M357)"
        return "L (M20)"
    if text.startswith("J2"):
        if text.startswith("J2b") or "M241" in text:
            return "J2b (M241)"
        if text.startswith("J2a") or "M410" in text:
            return "J2a (M410)"
        return "J2 (M172)"
    if text.startswith("J"):
        return "J1 (M267)" if text.startswith("J1") else "J"
    if text.startswith("C"):
        return "C"
    if text.startswith("O"):
        return "O"
    if text.startswith("G"):
        return "G"
    if text.startswith("Q"):
        return "Q"
    if text.startswith("T"):
        return "T (M70)"
    if text.startswith("E"):
        return "E (M96)"
    if text.startswith("I"):
        return "I (M170)"
    return "other"


def _y_call_from_anno(cols: list[str]) -> str | None:
    if len(cols) <= _ANNO_Y_MANUAL:
        return None
    manual = cols[_ANNO_Y_MANUAL].strip()
    auto = cols[_ANNO_Y_ISOGG].strip() if len(cols) > _ANNO_Y_ISOGG else ""
    for raw in (manual, auto):
        collapsed = collapse_y_haplogroup(raw)
        if collapsed is not None:
            return collapsed
    return None


def collapse_mt_haplogroup(raw: str) -> str | None:
    """Keep the first AADR token (M36d;M → M36d) for prefix counts."""
    text = (raw or "").split(";")[0].strip().replace(" ", "")
    if not text or text in {".", ".."} or text.lower().startswith("n/a"):
        return None
    return text


def score_mt_markers(query_index: dict[tuple[str, int], dict]) -> list[HaplogroupMarkerCall]:
    calls = score_markers(query_index, MT_MARKERS)
    by_label = {call.haplogroup: call for call in calls}
    m_call = by_label.get("M (10400T)")
    r_call = by_label.get("R (12705C)")
    if m_call and r_call and m_call.status == "derived" and r_call.status == "derived":
        m_call.status = "conflict"
        r_call.status = "conflict"
    j_call = by_label.get("J (16069T)")
    t_call = by_label.get("T (14905A)")
    if j_call and t_call and j_call.status == "derived" and t_call.status == "derived":
        j_call.status = "conflict"
        t_call.status = "conflict"
    for child, parent in (
        ("M2 (1780C)", "M (10400T)"),
        ("M3 (4580C)", "M (10400T)"),
        ("M4 (6620C)", "M (10400T)"),
        ("M5 (12477C)", "M (10400T)"),
        ("M6 (10640C)", "M (10400T)"),
        ("M30 (15431A)", "M (10400T)"),
        ("M33 (2361A)", "M (10400T)"),
        ("M36 (7271G)", "M (10400T)"),
        ("R5 (8594C)", "R (12705C)"),
        ("R6 (12285C)", "R (12705C)"),
        ("R7 (13105G)", "R (12705C)"),
        ("R8 (13215C)", "R (12705C)"),
        ("R30 (8584A)", "R (12705C)"),
        ("U (12308G)", "R (12705C)"),
        ("U1 (12879C)", "U (12308G)"),
        ("U2 (16051G)", "U (12308G)"),
        ("U7 (5360T)", "U (12308G)"),
        ("HV (14766C)", "R (12705C)"),
        ("H (7028C)", "HV (14766C)"),
        ("H3 (6776C)", "H (7028C)"),
        ("H3b (2581G)", "H3 (6776C)"),
        ("H3b (2581G)", "H (7028C)"),
        ("JT (16126C)", "R (12705C)"),
        ("J (16069T)", "JT (16126C)"),
        ("T (14905A)", "JT (16126C)"),
    ):
        _child_needs_parent(by_label, child, parent)
    for label in ("N1 (10238C)", "W (3505G)"):
        _incompatible_if_other_derived(by_label, label, "M (10400T)")
        _incompatible_if_other_derived(by_label, label, "R (12705C)")
    _incompatible_if_other_derived(by_label, "N1 (10238C)", "W (3505G)")
    _incompatible_if_other_derived(by_label, "W (3505G)", "N1 (10238C)")
    return calls


def load_anno_mt_by_group(anno_path: Path) -> dict[str, list[str | None]]:
    """AADR group ID → mt haplogroup tokens (all sexes; None = uncalled)."""
    out: dict[str, list[str | None]] = {}
    if not anno_path.exists():
        return out
    with anno_path.open(encoding="utf-8", errors="replace") as handle:
        header = handle.readline()
        if not header:
            return out
        for line in handle:
            cols = line.rstrip("\n").split("\t")
            if len(cols) <= _ANNO_MT:
                continue
            group = cols[_ANNO_GROUP].strip()
            out.setdefault(group, []).append(collapse_mt_haplogroup(cols[_ANNO_MT]))
    return out


def _prefix_count(tokens: list[str], prefix: str, exclude: tuple[str, ...]) -> int:
    n = 0
    for token in tokens:
        if not token.startswith(prefix):
            continue
        if any(token.startswith(ex) for ex in exclude):
            continue
        n += 1
    return n


def aadr_mt_frequency_rows(
    anno_by_pop: dict[str, list[str | None]],
    groups: dict[str, tuple[str, ...]],
    sample_status: dict[str, str],
) -> list[HaplogroupRow]:
    rows: list[HaplogroupRow] = []
    for label, marker, sample_key, prefix, exclude in _MT_FREQ_ROWS:
        group_counts: dict[str, HaplogroupGroupCount] = {}
        for display, pops in groups.items():
            tokens = [tok for pop in pops for tok in anno_by_pop.get(pop, [])]
            called = [tok for tok in tokens if tok is not None]
            n_called = len(called)
            n = _prefix_count(called, prefix, exclude) if prefix else 0
            percent = (100.0 * n / n_called) if n_called else None
            group_counts[display] = HaplogroupGroupCount(n=n, n_called=n_called, percent=percent)
        if not prefix and not any(item.n for item in group_counts.values()):
            continue
        rows.append(
            HaplogroupRow(
                haplogroup=label,
                marker=marker,
                sample_status=sample_status.get(sample_key or label, "no-call"),
                groups=group_counts,
            )
        )
    other_groups: dict[str, HaplogroupGroupCount] = {}
    named_prefixes = tuple(prefix for _l, _m, _s, prefix, _e in _MT_FREQ_ROWS if prefix)
    for display, pops in groups.items():
        called = [tok for pop in pops for tok in anno_by_pop.get(pop, []) if tok is not None]
        n_called = len(called)
        n = sum(1 for tok in called if not any(tok.startswith(p) for p in named_prefixes))
        other_groups[display] = HaplogroupGroupCount(
            n=n,
            n_called=n_called,
            percent=(100.0 * n / n_called) if n_called else None,
        )
    if any(item.n for item in other_groups.values()):
        rows.append(
            HaplogroupRow(
                haplogroup="other",
                marker="",
                sample_status="no-call",
                groups=other_groups,
            )
        )
    return rows


def load_anno_y_by_group(anno_path: Path) -> dict[str, list[str | None]]:
    """AADR group ID → Y haplogroup labels for males (None = uncalled)."""
    out: dict[str, list[str | None]] = {}
    if not anno_path.exists():
        return out
    with anno_path.open(encoding="utf-8", errors="replace") as handle:
        header = handle.readline()
        if not header:
            return out
        for line in handle:
            cols = line.rstrip("\n").split("\t")
            if len(cols) <= _ANNO_Y_ISOGG:
                continue
            group = cols[_ANNO_GROUP].strip()
            sex = cols[_ANNO_SEX].strip().upper() if len(cols) > _ANNO_SEX else ""
            if sex.startswith("F"):
                continue
            out.setdefault(group, []).append(_y_call_from_anno(cols))
    return out


def _group_counts(calls: list[str | None]) -> dict[str, HaplogroupGroupCount]:
    called = [c for c in calls if c is not None]
    n_called = len(called)
    counts = Counter(called)
    rows: dict[str, HaplogroupGroupCount] = {}
    for label, _marker, _prefix in _FREQ_ROWS:
        n = counts.get(label, 0)
        if label == "other":
            named = {name for name, _, _ in _FREQ_ROWS if name != "other"}
            n = sum(v for k, v in counts.items() if k not in named)
        percent = (100.0 * n / n_called) if n_called else None
        rows[label] = HaplogroupGroupCount(n=n, n_called=n_called, percent=percent)
    return rows


def aadr_frequency_rows(
    anno_by_pop: dict[str, list[str | None]],
    groups: dict[str, tuple[str, ...]],
    sample_status: dict[str, str],
) -> list[HaplogroupRow]:
    per_group = {
        display: _group_counts([call for pop in pops for call in anno_by_pop.get(pop, [])])
        for display, pops in groups.items()
    }
    rows: list[HaplogroupRow] = []
    for label, marker, sample_key in _FREQ_ROWS:
        if label == "other":
            continue
        rows.append(
            HaplogroupRow(
                haplogroup=label,
                marker=marker,
                sample_status=sample_status.get(sample_key or label, "no-call"),
                groups={name: per_group[name][label] for name in groups},
            )
        )
    other_groups = {name: per_group[name]["other"] for name in groups}
    if any(item.n for item in other_groups.values()):
        rows.append(
            HaplogroupRow(
                haplogroup="other",
                marker="",
                sample_status="no-call",
                groups=other_groups,
            )
        )
    return rows


SAMPLE_STATUS_NOTES: tuple[str, ...] = (
    "This sample column: derived = yes, this file has the mutation that defines that haplogroup.",
    "Ancestral = no, this file has the older allele, so that haplogroup is ruled out.",
    "No-call = that defining SNP is missing or unreadable in this VCF, so we cannot say yes or no.",
    "Conflict = markers disagree (a downstream SNP looks like yes while a parent lineage is no). Do not treat the downstream yes as a call.",
)


def compare_haplogroups(
    query_index: dict[tuple[str, int], dict],
    *,
    settings: Settings,
) -> HaplogroupResult:
    notes = [
        "One person has one Y haplogroup. Percentages are AADR male counts in each scored group, not mixture weights.",
        "R1a1 here is the M17 / M198 lineage (also written R1a1a). M17 on SNP arrays is an indel and is marked conflict if M173 is ancestral.",
        "H3b is Z13871 (ISOGG / YFull), under H3-Z5857 and H-M69. GSA files often lack Z13871, so H3b is no-call even when a sequencing test reports H3b; that is missing chip content, not a rival haplogroup.",
        "Same for other South Asian backbones added from ISOGG/Yleaf: R2-M124, R1a-Z93, J2a-M410, J2b-M241, L1-M27, L3-M357, T-M70, E-M96. No-call means the SNP is not on this chip.",
        "Y markers match by rsID first, then published hg19 and hg38 positions, so a GSA/gtc2vcf file is not scored at the wrong physical site.",
        "AADR Y labels are the ISOGG column in v66.1 HO .anno (YFull-based automatic calls). "
        "Several HO community labels have only a few samples (often 2–9).",
    ]
    mt_notes = [
        "One person has one mtDNA haplogroup. Percentages are AADR published mt calls, not mixture weights.",
        "Markers use rCRS / PhyloTree Build 17 positions. rCRS is haplogroup H2a2a1, so REF=rCRS can already be the derived allele (for example 12705C = R).",
        "M4, M36, and other South Asian / West Eurasian subclades use PhyloTree coding SNPs (same parent/child conflict rule as Y). Recurrent sites such as 1888 (M5 and T) are skipped; M5 is the M5a T12477C SNP.",
        "GSA files often carry these mt sites. No-call still means missing chip content, not a rival haplogroup.",
        "A dash means that label has no mt haplogroup in the AADR .anno (common for 1000 Genomes HO groups such as STU, ITU, GIH, PJL, BEB).",
        "Nested rows (M3 inside M, H3 inside H) can add up to more than 100%.",
    ]
    n_y = sum(1 for (chrom, _pos) in query_index if chrom == "Y")
    n_mt = sum(1 for (chrom, _pos) in query_index if chrom == "MT")
    y_present = n_y > 0
    mt_present = n_mt > 0

    calls = score_y_markers(query_index) if y_present else []
    best = best_sample_haplogroup(calls) if y_present else None
    sample_status = {c.haplogroup: c.status for c in calls}
    mt_calls = score_mt_markers(query_index) if mt_present else []
    mt_best = best_sample_haplogroup(mt_calls, _MT_SPECIFICITY) if mt_present else None
    mt_status = {c.haplogroup: c.status for c in mt_calls}

    if y_present:
        if best:
            notes.append(f"Deepest consistent derived marker in this file: {best}.")
        else:
            notes.append("No backbone Y marker is derived in this file (or derived calls conflict).")
        notes.append(marker_quality_note(calls, label="Y"))
    else:
        notes = [
            "No chrY SNPs in this VCF (typical of a mother / filtered extract), so the Y haplogroup table is omitted."
        ]
    if mt_present:
        if mt_best:
            mt_notes.append(f"Deepest consistent derived mt marker in this file: {mt_best}.")
        else:
            mt_notes.append("No backbone mt marker is derived in this file (or calls conflict).")
        mt_notes.append(marker_quality_note(mt_calls, label="mtDNA"))
    else:
        mt_notes = [
            "No chrMT SNPs in this VCF, so the mtDNA haplogroup table is omitted."
        ]

    anno_path = settings.aadr_anno
    groups = settings.caste_groups()
    rows: list[HaplogroupRow] = []
    mt_rows: list[HaplogroupRow] = []
    if not anno_path.exists():
        if y_present:
            notes.append("AADR .anno not found; group percentages need data/references/aadr/*.anno.")
        if mt_present:
            mt_notes.append("AADR .anno not found; group percentages need data/references/aadr/*.anno.")
        return HaplogroupResult(
            available=y_present,
            sample_best=best,
            markers=calls,
            rows=[],
            notes=notes,
            mt_available=mt_present,
            mt_sample_best=mt_best,
            mt_markers=mt_calls,
            mt_rows=[],
            mt_notes=mt_notes,
            status_notes=list(SAMPLE_STATUS_NOTES) if y_present or mt_present else [],
        )
    if y_present:
        anno_by_pop = load_anno_y_by_group(anno_path)
        rows = aadr_frequency_rows(anno_by_pop, groups, sample_status)
        n_males = sum(len(anno_by_pop.get(pop, [])) for pops in groups.values() for pop in pops)
        notes.append(f"AADR males in scored caste/community groups: {n_males}.")
    if mt_present:
        mt_by_pop = load_anno_mt_by_group(anno_path)
        mt_rows = aadr_mt_frequency_rows(mt_by_pop, groups, mt_status)
        n_mt_called = sum(
            1
            for pops in groups.values()
            for pop in pops
            for tok in mt_by_pop.get(pop, [])
            if tok is not None
        )
        mt_with_calls = [
            name
            for name, pops in groups.items()
            if any(tok is not None for pop in pops for tok in mt_by_pop.get(pop, []))
        ]
        if mt_with_calls:
            mt_notes.append(
                "Labels with at least one AADR mt call: "
                + ", ".join(mt_with_calls)
                + f" ({n_mt_called} calls). Example: Vellalar/VLR is a small HO set that happens to have mt filled in."
            )
        else:
            mt_notes.append("No scored label has an AADR mt haplogroup in this .anno.")
    return HaplogroupResult(
        available=y_present,
        sample_best=best,
        markers=calls,
        rows=rows,
        notes=notes,
        mt_available=mt_present,
        mt_sample_best=mt_best,
        mt_markers=mt_calls,
        mt_rows=mt_rows,
        mt_notes=mt_notes,
        status_notes=list(SAMPLE_STATUS_NOTES) if y_present or mt_present else [],
    )
