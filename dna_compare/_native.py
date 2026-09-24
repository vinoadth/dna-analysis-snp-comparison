from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from dna_compare.eigenstrat import AadrPanel

try:
    import dna_compare_rs as _rs
except ImportError:
    _rs = None

_AUTO = {str(i) for i in range(1, 23)}


def native_available() -> bool:
    return _rs is not None


def population_allele1_freq_native(
    geno_path: str,
    sample_indices: list[int],
    nsnp: int,
) -> np.ndarray | None:
    if _rs is None:
        return None
    return np.asarray(
        _rs.population_allele1_freq(geno_path, sample_indices, nsnp),
        dtype=np.float32,
    )


def build_population_freqs_single_pass_native(
    geno_path: str,
    nsnp: int,
    ind_indices: list[int],
    ind_group_offsets: list[int],
    ind_group_ids: list[int],
    n_groups: int,
) -> list[np.ndarray] | None:
    if _rs is None:
        return None
    arrays = _rs.build_population_freqs_single_pass(
        geno_path,
        nsnp,
        ind_indices,
        ind_group_offsets,
        ind_group_ids,
        n_groups,
    )
    return [np.asarray(arr, dtype=np.float32) for arr in arrays]


def align_query_to_panel_native(panel: AadrPanel, query_index: dict[tuple[str, int], dict]):
    if _rs is None:
        return None
    snps = panel.snps()
    chroms = [snp.chrom for snp in snps]
    positions = [int(snp.pos) for snp in snps]
    allele1 = [snp.allele1 for snp in snps]
    allele2 = [snp.allele2 for snp in snps]
    q_chroms: list[str] = []
    q_positions: list[int] = []
    q_ref: list[str] = []
    q_alt: list[str] = []
    q_dosage: list[float] = []
    for (chrom, pos), row in query_index.items():
        dosage = row.get("dosage_alt")
        if dosage is None:
            continue
        ref = str(row.get("ref") or "")
        alt = str(row.get("alt") or "")
        if len(ref) != 1 or len(alt) != 1:
            continue
        q_chroms.append(str(chrom))
        q_positions.append(int(pos))
        q_ref.append(ref.upper())
        q_alt.append(alt.upper())
        q_dosage.append(float(dosage))
    query_freqs, n_overlap = _rs.align_query_to_panel(
        chroms,
        positions,
        allele1,
        allele2,
        q_chroms,
        q_positions,
        q_ref,
        q_alt,
        q_dosage,
    )
    return np.asarray(query_freqs, dtype=np.float32), int(n_overlap)


def score_pgs_batch_native(
    by_rsid: dict[str, dict],
    by_pos: dict[tuple[str, int], dict],
    variant_batches: list[list[dict]],
) -> list[tuple[float, int, int]] | None:
    if _rs is None:
        return None
    flat: list[dict] = []
    offsets: list[int] = [0]
    for batch in variant_batches:
        flat.extend(batch)
        offsets.append(len(flat))
    if not flat:
        return [(0.0, 0, 0) for _ in variant_batches]
    rsids = [str(v.get("rsid") or "") for v in flat]
    chroms = [str(v.get("chrom") or "") for v in flat]
    positions = [int(v.get("pos") or 0) for v in flat]
    effects = [str(v.get("effect") or "") for v in flat]
    others = [str(v.get("other") or "") for v in flat]
    weights = [float(v.get("weight") or 0.0) for v in flat]
    q_keys = list(by_pos.keys())
    q_chroms = [str(k[0]) for k in q_keys]
    q_positions = [int(k[1]) for k in q_keys]
    q_ref = [str(by_pos[k].get("ref") or "").upper() for k in q_keys]
    q_alt = [str(by_pos[k].get("alt") or "").upper() for k in q_keys]
    q_dosage = [float(by_pos[k].get("dosage_alt") or 0.0) for k in q_keys]
    q_rsids = [str(by_pos[k].get("rsid") or "") for k in q_keys]
    totals, used, counts = _rs.score_pgs_batch(
        rsids,
        chroms,
        positions,
        effects,
        others,
        weights,
        offsets,
        q_rsids,
        q_chroms,
        q_positions,
        q_ref,
        q_alt,
        q_dosage,
    )
    out: list[tuple[float, int, int]] = []
    for i in range(len(offsets) - 1):
        n_score = offsets[i + 1] - offsets[i]
        u = int(used[i])
        out.append((float(totals[i]), u, n_score))
    return out


def king_kinship_native(
    query_index: dict[tuple[str, int], dict],
    other_index: dict[tuple[str, int], dict],
    query_by_rsid: dict[str, dict],
    other_by_rsid: dict[str, dict],
) -> dict | None:
    if _rs is None:
        return None
    q = _pack_index(query_index)
    o = _pack_index(other_index)
    return dict(
        _rs.king_kinship(
            q["chroms"],
            q["positions"],
            q["ref"],
            q["alt"],
            q["dosage"],
            q["genotype"],
            q["rsid"],
            q["filter"],
            q["igc"],
            q["gq"],
            o["chroms"],
            o["positions"],
            o["ref"],
            o["alt"],
            o["dosage"],
            o["genotype"],
            o["rsid"],
            o["filter"],
            o["igc"],
            o["gq"],
        )
    )


def _pack_index(index: dict[tuple[str, int], dict]) -> dict:
    keys = [(c, p) for (c, p) in index if c in _AUTO]
    return {
        "chroms": [str(c) for c, _ in keys],
        "positions": [int(p) for _, p in keys],
        "ref": [str(index[k].get("ref") or "").upper() for k in keys],
        "alt": [str(index[k].get("alt") or "").upper() for k in keys],
        "dosage": [float(index[k].get("dosage_alt") or -1.0) for k in keys],
        "genotype": [str(index[k].get("genotype") or "") for k in keys],
        "rsid": [str(index[k].get("rsid") or "") for k in keys],
        "filter": [str(index[k].get("vcf_filter") or ".") for k in keys],
        "igc": [_float_or_neg(index[k].get("igc")) for k in keys],
        "gq": [_int_or_neg(index[k].get("gq")) for k in keys],
    }


def _float_or_neg(value: object) -> float:
    try:
        return float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return -1.0


def _int_or_neg(value: object) -> int:
    try:
        return int(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return -1
