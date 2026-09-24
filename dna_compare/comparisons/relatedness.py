from __future__ import annotations

from dna_compare._native import king_kinship_native
from dna_compare.genotype_index import build_rsid_index, normalize_rsid
from dna_compare.models import RelatednessResult

_AUTO = {str(i) for i in range(1, 23)}
_COMP = str.maketrans("ACGT", "TGCA")
_MIN_OVERLAP = 200
_THIN_BP = 100_000
_THIN_IF_AT_LEAST = 2500
_MIN_IGC = 0.15
_MIN_GQ = 7

RELATEDNESS_NOTES = (
    "Optional second VCF: parent, relative, friend, or anyone else. They do not have to be related.",
    "Kinship is the KING-robust estimator on overlapping autosomal SNPs (not a legal relationship test).",
    "Same-population strangers can look slightly closer than people from different continents.",
)


def _diploid_alt_dosage(row: dict) -> float | None:
    """Alt-allele count 0/1/2. Haploid 0 or 1 on autosomes is treated as 0/0 or 1/1."""
    token = str(row.get("genotype") or "").split(":", 1)[0].replace("|", "/")
    if "/" in token:
        parts = token.split("/")
        if len(parts) != 2 or any(part == "." for part in parts):
            return None
        try:
            return float(int(parts[0]) + int(parts[1]))
        except ValueError:
            return None
    if token in {"0", "1"}:
        return 0.0 if token == "0" else 2.0
    return None


def _aligned_dosages(left: dict, right: dict) -> tuple[float, float] | None:
    a_ref, a_alt = str(left.get("ref") or "").upper(), str(left.get("alt") or "").upper()
    b_ref, b_alt = str(right.get("ref") or "").upper(), str(right.get("alt") or "").upper()
    da, db = _diploid_alt_dosage(left), _diploid_alt_dosage(right)
    if da is None or db is None or len(a_ref) != 1 or len(a_alt) != 1 or len(b_ref) != 1 or len(b_alt) != 1:
        return None
    if a_ref == b_ref and a_alt == b_alt:
        return da, db
    if a_ref == b_alt and a_alt == b_ref:
        return da, 2.0 - db
    a_ref_c, a_alt_c = a_ref.translate(_COMP), a_alt.translate(_COMP)
    if a_ref_c == b_ref and a_alt_c == b_alt:
        return da, db
    if a_ref_c == b_alt and a_alt_c == b_ref:
        return da, 2.0 - db
    return None


def _usable_rsid(value: object) -> str | None:
    return normalize_rsid(value)


def _qc_ok(row: dict) -> bool:
    filt = str(row.get("vcf_filter") or ".").strip()
    if filt not in {".", "", "PASS", "pass"} and not filt.upper().startswith("PASS"):
        return False
    igc = row.get("igc")
    if igc is not None:
        try:
            if float(igc) < _MIN_IGC:
                return False
        except (TypeError, ValueError):
            return False
    gq = row.get("gq")
    if gq is not None:
        try:
            if int(gq) < _MIN_GQ:
                return False
        except (TypeError, ValueError):
            return False
    return True


def _rsid_index(index: dict[tuple[str, int], dict]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for key, row in index.items():
        if key[0] not in _AUTO:
            continue
        rsid = _usable_rsid(row.get("rsid"))
        if rsid and rsid not in out:
            out[rsid] = row
    return out


def _n_auto(index: dict[tuple[str, int], dict]) -> int:
    return sum(1 for key, row in index.items() if key[0] in _AUTO and _diploid_alt_dosage(row) is not None)


def _pair_sites(
    query_index: dict[tuple[str, int], dict],
    other_index: dict[tuple[str, int], dict],
) -> tuple[list[tuple[tuple[str, int], tuple[float, float], str]], int]:
    other_by_rsid = _rsid_index(other_index)
    used_other: set[tuple[str, int]] = set()
    pairs: list[tuple[tuple[str, int], tuple[float, float], str]] = []
    n_qc_dropped = 0
    for key, left in query_index.items():
        if key[0] not in _AUTO:
            continue
        how = "pos"
        right = other_index.get(key)
        if right is None:
            rsid = _usable_rsid(left.get("rsid"))
            right = other_by_rsid.get(rsid) if rsid else None
            how = "rsid"
        if right is None:
            continue
        rkey = (str(right.get("chrom") or key[0]), int(right.get("pos") or key[1]))
        if rkey in used_other:
            continue
        if not _qc_ok(left) or not _qc_ok(right):
            n_qc_dropped += 1
            continue
        aligned = _aligned_dosages(left, right)
        if aligned is None:
            continue
        used_other.add(rkey)
        pairs.append((key, aligned, how))
    return pairs, n_qc_dropped


def _thin_pairs(
    pairs: list[tuple[tuple[str, int], tuple[float, float], str]],
) -> tuple[list[tuple[tuple[str, int], tuple[float, float], str]], int]:
    if len(pairs) < _THIN_IF_AT_LEAST:
        return pairs, 0
    ordered = sorted(pairs, key=lambda item: (int(item[0][0]) if item[0][0].isdigit() else 99, item[0][1]))
    kept: list[tuple[tuple[str, int], tuple[float, float], str]] = []
    last_chrom = None
    last_pos = -10**18
    for item in ordered:
        chrom, pos = item[0]
        if chrom != last_chrom or pos - last_pos >= _THIN_BP:
            kept.append(item)
            last_chrom, last_pos = chrom, pos
    if len(kept) < _MIN_OVERLAP:
        return pairs, 0
    return kept, len(pairs) - len(kept)


def classify_relationship(kinship: float | None, ibs0: int, n_snps: int) -> str:
    if kinship is None or n_snps < _MIN_OVERLAP:
        return "not enough overlapping autosomal SNPs"
    ibs0_rate = ibs0 / n_snps
    if kinship >= 0.354:
        return "same person or identical twin (duplicate / very high sharing)"
    if kinship >= 0.177:
        # Real arrays have a few opposite-homozygote errors; synthetic parent–child is ~0.
        if ibs0_rate < 0.005:
            return "parent–child (first-degree; almost no opposite homozygotes)"
        return "first-degree (full siblings, or parent–child with some array noise)"
    if kinship >= 0.088:
        return "second-degree (half-sibling, uncle/aunt–niece/nephew, grandparent)"
    if kinship >= 0.044:
        return "third-degree (first cousin or similar)"
    if kinship >= 0.022:
        return "fourth-degree or distant relative"
    return "unrelated or very distant"


def _reliability(n_snps: int) -> str:
    if n_snps >= 8000:
        return "high"
    if n_snps >= 1000:
        return "ok"
    return "low"


def compare_relatedness(
    query_index: dict[tuple[str, int], dict],
    other_index: dict[tuple[str, int], dict],
    *,
    other_filename: str,
    other_sample_id: str | None = None,
    query_sample_id: str | None = None,
    query_filename: str | None = None,
    query_assembly: str | None = None,
    other_assembly: str | None = None,
    query_lifted_to: str | None = None,
    other_lifted_to: str | None = None,
    query_by_rsid: dict[str, dict] | None = None,
    other_by_rsid: dict[str, dict] | None = None,
) -> RelatednessResult:
    if query_by_rsid is None:
        query_by_rsid = build_rsid_index(query_index)
    if other_by_rsid is None:
        other_by_rsid = build_rsid_index(other_index)

    native = king_kinship_native(query_index, other_index, query_by_rsid, other_by_rsid)
    if native is not None:
        return _relatedness_from_native(
            native,
            other_filename=other_filename,
            other_sample_id=other_sample_id,
            query_sample_id=query_sample_id,
            query_filename=query_filename,
            query_assembly=query_assembly,
            other_assembly=other_assembly,
            query_lifted_to=query_lifted_to,
            other_lifted_to=other_lifted_to,
            query_index=query_index,
            other_index=other_index,
        )

    pairs, n_qc_dropped = _pair_sites(query_index, other_index)
    pairs, n_pruned = _thin_pairs(pairs)
    n_matched_pos = sum(1 for _key, _aligned, how in pairs if how == "pos")
    n_matched_rsid = sum(1 for _key, _aligned, how in pairs if how == "rsid")
    ibs0 = ibs1 = ibs2 = 0
    shared_sum = 0.0
    het_a = het_b = both_het = 0
    for _key, (da, db), _how in pairs:
        shared = 2 - abs(da - db)
        if da in {0.0, 2.0} and db in {0.0, 2.0} and da != db:
            ibs0 += 1
        elif da == 1.0 and db == 1.0:
            ibs2 += 1
            both_het += 1
        elif {da, db} <= {0.0, 1.0, 2.0} and abs(da - db) == 1:
            ibs1 += 1
        elif da == db:
            ibs2 += 1
        else:
            ibs1 += 1
        if da == 1.0:
            het_a += 1
        if db == 1.0:
            het_b += 1
        shared_sum += shared
    n = ibs0 + ibs1 + ibs2
    notes = list(RELATEDNESS_NOTES)
    n_query = _n_auto(query_index)
    n_other = _n_auto(other_index)
    empty = RelatednessResult(
        available=False,
        other_filename=other_filename,
        other_sample_id=other_sample_id,
        query_sample_id=query_sample_id,
        n_snps=n,
        n_matched_pos=n_matched_pos,
        n_matched_rsid=n_matched_rsid,
        n_qc_dropped=n_qc_dropped,
        n_pruned=n_pruned,
        notes=notes,
    )
    if n < _MIN_OVERLAP:
        notes.append(
            f"Only {n} overlapping autosomal SNPs after QC"
            + (f" ({n_matched_rsid} matched by rsID)" if n_matched_rsid else "")
            + "; need the same chip or a denser overlap. "
            f"Query has {n_query} autosomal SNPs, second file has {n_other}."
        )
        if n_qc_dropped:
            notes.append(f"Dropped {n_qc_dropped} overlapping sites with low IGC/GQ or a non-PASS filter.")
        return empty
    mean_ibs = shared_sum / (2.0 * n)
    het_sum = het_a + het_b
    if het_sum == 0:
        kinship = 0.5 if ibs0 == 0 else 0.0
    else:
        # KING within-family / robust (Manichaikul 2010 eq. 9):
        # identical ≈ 0.5, parent–child ≈ 0.25, first cousin ≈ 0.0625, unrelated ≈ 0.
        kinship = (both_het - 2 * ibs0) / het_sum
    relationship = classify_relationship(kinship, ibs0, n)
    reliability = _reliability(n)
    het_rate_a = round(het_a / n, 4)
    het_rate_b = round(het_b / n, 4)
    notes.append(
        f"Compared {n} autosomal SNPs to {other_sample_id or other_filename}"
        + (f" ({n_matched_pos} by position, {n_matched_rsid} by rsID)" if n_matched_rsid else "")
        + ". IBS2 = both alleles match, IBS0 = opposite homozygotes."
    )
    if n_qc_dropped:
        notes.append(f"Dropped {n_qc_dropped} overlapping sites with low IGC/GQ or a non-PASS filter.")
    if n_pruned:
        notes.append(
            f"Thinned {n_pruned} nearby SNPs (kept one per {_THIN_BP // 1000} kb) so dense-array LD does not inflate kinship."
        )
    smaller = min(n_query, n_other) or 1
    overlap_frac = n / smaller
    if smaller >= 1000 and overlap_frac < 0.08:
        notes.append(
            f"Only {overlap_frac:.1%} of the smaller file overlapped "
            f"({n_query} vs {n_other} autosomal SNPs). Different chip or assembly; treat distant calls cautiously."
        )
    query_build = query_lifted_to or query_assembly
    other_build = other_lifted_to or other_assembly
    if query_assembly and other_assembly and query_build != other_build:
        notes.append(
            f"Coordinates may still disagree (query {query_assembly}"
            + (f"→{query_lifted_to}" if query_lifted_to else "")
            + f" vs other {other_assembly}"
            + (f"→{other_lifted_to}" if other_lifted_to else "")
            + "). rsID matching is used when positions do not line up."
        )
    if n_matched_rsid and not (query_assembly and other_assembly and query_build != other_build):
        notes.append(
            "Some SNPs matched by rsID rather than chromosome:position. Typical when one VCF is hg38 and the other is hg19."
        )
    if het_rate_a < 0.08 or het_rate_a > 0.50 or het_rate_b < 0.08 or het_rate_b > 0.50:
        notes.append(
            f"Heterozygosity on overlapping SNPs looks unusual for a SNP array "
            f"(query {het_rate_a:.2f}, other {het_rate_b:.2f}). Check no-calls, haploid encoding, or a failed lift."
        )
    elif abs(het_rate_a - het_rate_b) > 0.12:
        notes.append(
            f"Heterozygosity differs a lot between the two files ({het_rate_a:.2f} vs {het_rate_b:.2f}); "
            "one may have more no-calls or a different calling pipeline."
        )
    if query_sample_id and other_sample_id and query_sample_id == other_sample_id:
        notes.append(f"Both files use sample ID {query_sample_id}; check that this is not the same person twice.")
    if query_filename and other_filename and query_filename == other_filename:
        notes.append("Both files have the same filename.")
    if reliability == "low":
        notes.append(
            f"Only {n} SNPs after QC/thinning — first-degree calls can still work, but cousin vs unrelated is noisy."
        )
    elif reliability == "high":
        notes.append(f"{n} well-spaced autosomal SNPs is enough for a stable 3rd-degree call on real array data.")
    kinship_r = None if kinship is None else round(float(kinship), 4)
    return _build_relatedness_result(
        available=True,
        other_filename=other_filename,
        other_sample_id=other_sample_id,
        query_sample_id=query_sample_id,
        n=n,
        n_matched_pos=n_matched_pos,
        n_matched_rsid=n_matched_rsid,
        n_qc_dropped=n_qc_dropped,
        n_pruned=n_pruned,
        het_rate_a=het_rate_a,
        het_rate_b=het_rate_b,
        reliability=reliability,
        mean_ibs=mean_ibs,
        kinship_r=kinship_r,
        ibs0=ibs0,
        ibs1=ibs1,
        ibs2=ibs2,
        relationship=relationship,
        notes=notes,
    )


def _build_relatedness_result(
    *,
    available: bool,
    other_filename: str,
    other_sample_id: str | None,
    query_sample_id: str | None,
    n: int,
    n_matched_pos: int,
    n_matched_rsid: int,
    n_qc_dropped: int,
    n_pruned: int,
    het_rate_a: float,
    het_rate_b: float,
    reliability: str | None,
    mean_ibs: float,
    kinship_r: float | None,
    ibs0: int,
    ibs1: int,
    ibs2: int,
    relationship: str | None,
    notes: list[str],
) -> RelatednessResult:
    return RelatednessResult(
        available=available,
        other_filename=other_filename,
        other_sample_id=other_sample_id,
        query_sample_id=query_sample_id,
        n_snps=n,
        n_matched_pos=n_matched_pos,
        n_matched_rsid=n_matched_rsid,
        n_qc_dropped=n_qc_dropped,
        n_pruned=n_pruned,
        het_rate_query=het_rate_a,
        het_rate_other=het_rate_b,
        reliability=reliability,
        mean_ibs=round(mean_ibs, 4),
        kinship=kinship_r,
        ibs0=ibs0,
        ibs1=ibs1,
        ibs2=ibs2,
        ibs0_pct=round(100.0 * ibs0 / n, 2) if n else None,
        ibs1_pct=round(100.0 * ibs1 / n, 2) if n else None,
        ibs2_pct=round(100.0 * ibs2 / n, 2) if n else None,
        shared_pct=None if kinship_r is None else round(200.0 * kinship_r, 2),
        relationship=relationship,
        notes=notes,
    )


def _relatedness_from_native(
    native: dict,
    *,
    other_filename: str,
    other_sample_id: str | None,
    query_sample_id: str | None,
    query_filename: str | None,
    query_assembly: str | None,
    other_assembly: str | None,
    query_lifted_to: str | None,
    other_lifted_to: str | None,
    query_index: dict[tuple[str, int], dict],
    other_index: dict[tuple[str, int], dict],
) -> RelatednessResult:
    n = int(native.get("n_snps") or 0)
    notes = list(RELATEDNESS_NOTES)
    n_query = _n_auto(query_index)
    n_other = _n_auto(other_index)
    if n < _MIN_OVERLAP:
        notes.append(
            f"Only {n} overlapping autosomal SNPs after QC"
            + (f" ({int(native.get('n_matched_rsid') or 0)} matched by rsID)" if native.get("n_matched_rsid") else "")
            + "; need the same chip or a denser overlap. "
            f"Query has {n_query} autosomal SNPs, second file has {n_other}."
        )
        return RelatednessResult(
            available=False,
            other_filename=other_filename,
            other_sample_id=other_sample_id,
            query_sample_id=query_sample_id,
            n_snps=n,
            n_matched_pos=int(native.get("n_matched_pos") or 0),
            n_matched_rsid=int(native.get("n_matched_rsid") or 0),
            n_qc_dropped=int(native.get("n_qc_dropped") or 0),
            n_pruned=int(native.get("n_pruned") or 0),
            notes=notes,
        )
    kinship = float(native.get("kinship") or 0.0)
    ibs0 = int(native.get("ibs0") or 0)
    relationship = classify_relationship(kinship, ibs0, n)
    reliability = _reliability(n)
    n_matched_pos = int(native.get("n_matched_pos") or 0)
    n_matched_rsid = int(native.get("n_matched_rsid") or 0)
    n_qc_dropped = int(native.get("n_qc_dropped") or 0)
    n_pruned = int(native.get("n_pruned") or 0)
    het_rate_a = float(native.get("het_rate_query") or 0.0)
    het_rate_b = float(native.get("het_rate_other") or 0.0)
    mean_ibs = float(native.get("mean_ibs") or 0.0)
    notes.append(
        f"Compared {n} autosomal SNPs to {other_sample_id or other_filename}"
        + (f" ({n_matched_pos} by position, {n_matched_rsid} by rsID)" if n_matched_rsid else "")
        + ". IBS2 = both alleles match, IBS0 = opposite homozygotes."
    )
    if n_qc_dropped:
        notes.append(f"Dropped {n_qc_dropped} overlapping sites with low IGC/GQ or a non-PASS filter.")
    if n_pruned:
        notes.append(
            f"Thinned {n_pruned} nearby SNPs (kept one per {_THIN_BP // 1000} kb) so dense-array LD does not inflate kinship."
        )
    return _build_relatedness_result(
        available=True,
        other_filename=other_filename,
        other_sample_id=other_sample_id,
        query_sample_id=query_sample_id,
        n=n,
        n_matched_pos=n_matched_pos,
        n_matched_rsid=n_matched_rsid,
        n_qc_dropped=n_qc_dropped,
        n_pruned=n_pruned,
        het_rate_a=het_rate_a,
        het_rate_b=het_rate_b,
        reliability=reliability,
        mean_ibs=mean_ibs,
        kinship_r=round(kinship, 4),
        ibs0=ibs0,
        ibs1=int(native.get("ibs1") or 0),
        ibs2=int(native.get("ibs2") or 0),
        relationship=relationship,
        notes=notes,
    )
