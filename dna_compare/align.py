from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from dna_compare._native import align_query_to_panel_native


def query_allele1_freq(row: dict, allele1: str, allele2: str) -> float | None:
    dosage = row.get("dosage_alt")
    if dosage is None:
        return None
    q_ref, q_alt = row["ref"].upper(), row["alt"].upper()
    allele1, allele2 = allele1.upper(), allele2.upper()
    copies = 2.0
    gt = row.get("genotype") or ""
    if gt and "/" not in gt and "|" not in gt:
        copies = 1.0
    q_alt_freq = float(dosage) / copies
    if q_ref == allele2 and q_alt == allele1:
        return q_alt_freq
    if q_ref == allele1 and q_alt == allele2:
        return 1.0 - q_alt_freq
    return None


@dataclass
class PanelAlignment:
    query_freqs: np.ndarray
    n_overlap: int


def _align_python(panel, query_index: dict[tuple[str, int], dict]) -> PanelAlignment:
    snps = panel.snps()
    query = np.full(len(snps), np.nan, dtype=np.float32)
    n_overlap = 0
    for snp in snps:
        row = query_index.get((snp.chrom, snp.pos))
        if row is None:
            continue
        freq = query_allele1_freq(row, snp.allele1, snp.allele2)
        if freq is None:
            continue
        query[snp.index] = freq
        n_overlap += 1
    return PanelAlignment(query_freqs=query, n_overlap=n_overlap)


def align_query_to_panel(panel, query_index: dict[tuple[str, int], dict]) -> PanelAlignment | None:
    if panel is None or not panel.available:
        return None
    native = align_query_to_panel_native(panel, query_index)
    if native is not None:
        query_freqs, n_overlap = native
        return PanelAlignment(query_freqs=query_freqs, n_overlap=n_overlap)
    return _align_python(panel, query_index)


def sparse_overlap_columns(
    alignment: PanelAlignment,
    freq_arrays: list[np.ndarray],
) -> tuple[np.ndarray, list[np.ndarray], int]:
    q_vals: list[float] = []
    columns: list[list[float]] = [[] for _ in freq_arrays]
    for idx, qv in enumerate(alignment.query_freqs):
        if not np.isfinite(qv):
            continue
        freqs = [float(freq[idx]) for freq in freq_arrays]
        if any(not np.isfinite(val) for val in freqs):
            continue
        q_vals.append(float(qv))
        for i, val in enumerate(freqs):
            columns[i].append(val)
    if not q_vals:
        return np.asarray([], dtype=float), [np.asarray([], dtype=float) for _ in freq_arrays], 0
    return (
        np.asarray(q_vals, dtype=float),
        [np.asarray(col, dtype=float) for col in columns],
        len(q_vals),
    )


def aligned_columns_from_freqs(
    query_index: dict[tuple[str, int], dict],
    panel,
    freq_arrays: list[np.ndarray],
    alignment: PanelAlignment | None = None,
) -> tuple[np.ndarray, list[np.ndarray]]:
    if alignment is not None:
        target, cols, _n = sparse_overlap_columns(alignment, freq_arrays)
        return target, cols
    target: list[float] = []
    columns: list[list[float]] = [[] for _ in freq_arrays]
    for snp in panel.snps():
        row = query_index.get((snp.chrom, snp.pos))
        if row is None:
            continue
        qv = query_allele1_freq(row, snp.allele1, snp.allele2)
        if qv is None:
            continue
        vals = [float(freq[snp.index]) for freq in freq_arrays]
        if any(not np.isfinite(val) for val in vals):
            continue
        target.append(qv)
        for i, val in enumerate(vals):
            columns[i].append(val)
    return np.asarray(target, dtype=float), [np.asarray(col, dtype=float) for col in columns]
