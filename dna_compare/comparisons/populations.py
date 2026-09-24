from __future__ import annotations

import numpy as np

from dna_compare.align import PanelAlignment, aligned_columns_from_freqs, sparse_overlap_columns
from dna_compare.config import Settings
from dna_compare.models import ComparisonBlock, PopulationEstimate


def _nnls_percentages(F: np.ndarray, q: np.ndarray) -> np.ndarray:
    weights, *_ = np.linalg.lstsq(F, q, rcond=None)
    weights = np.clip(weights, 0.0, None)
    total = float(weights.sum())
    if total <= 0:
        return np.zeros(F.shape[1])
    return 100.0 * weights / total


def compare_populations(
    query_index: dict[tuple[str, int], dict],
    *,
    panel,
    settings: Settings,
    groups: dict[str, tuple[str, ...]] | None = None,
    kind: str = "populations",
    alignment: PanelAlignment | None = None,
) -> ComparisonBlock:
    notes = [
        "Percentages are mixture weights on overlapping Human Origins SNPs from AADR v66.1, "
        "not ethnic, national, or caste identity."
    ]
    if panel is None or not panel.available:
        return ComparisonBlock(
            kind=kind,
            available=False,
            notes=["AADR HO panel not found under data/references/aadr/."],
        )
    groups = groups if groups is not None else settings.selected_population_groups()
    if not groups:
        return ComparisonBlock(kind=kind, available=False, notes=["No population packs selected."])

    labels: list[str] = []
    pop_freqs: list[np.ndarray] = []
    for label, aadr_pops in groups.items():
        try:
            pop_freqs.append(panel.population_allele1_freq(aadr_pops))
            labels.append(label)
        except KeyError:
            notes.append(f"No AADR samples for {label} ({', '.join(aadr_pops)})")

    if not labels:
        return ComparisonBlock(kind=kind, available=False, notes=notes)

    if alignment is not None:
        q_vec, columns, n_used = sparse_overlap_columns(alignment, pop_freqs)
    else:
        q_vec, columns = aligned_columns_from_freqs(query_index, panel, pop_freqs)
        n_used = int(q_vec.size)
    if n_used < 50:
        return ComparisonBlock(
            kind=kind,
            available=False,
            notes=notes + [f"Only {n_used} overlapping SNPs; need a denser SNP VCF."],
        )

    F = np.column_stack(columns)
    q_vec = np.asarray(q_vec, dtype=float)
    percents = _nnls_percentages(F, q_vec)
    estimates = []
    for label, pct, col in zip(labels, percents, columns):
        ibs = 1.0 - float(np.mean(np.abs(q_vec - np.asarray(col))))
        estimates.append(
            PopulationEstimate(
                population=label,
                percent=round(float(pct), 3),
                n_snps=n_used,
                mean_ibs=round(ibs, 4),
            )
        )
    estimates.sort(key=lambda item: item.percent, reverse=True)
    return ComparisonBlock(kind=kind, available=True, estimates=estimates, notes=notes)
