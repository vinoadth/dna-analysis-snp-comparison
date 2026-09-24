from __future__ import annotations

from dna_compare.comparisons.populations import compare_populations
from dna_compare.config import (
    CASTE_AADR_POPS,
    CASTE_GENERAL_NOTES,
    CASTE_LABEL_NOTES,
    CASTE_MISSING_PANEL_NOTES,
    Settings,
)
from dna_compare.models import ComparisonBlock


def caste_alias_notes(groups: dict[str, tuple[str, ...]] | None = None) -> list[str]:
    """Notes for scored labels plus jatis that have no public HO bar."""
    selected = groups if groups is not None else CASTE_AADR_POPS
    notes = list(CASTE_GENERAL_NOTES)
    notes.extend(CASTE_LABEL_NOTES[name] for name in selected if name in CASTE_LABEL_NOTES)
    notes.extend(CASTE_MISSING_PANEL_NOTES)
    return notes


def compare_caste(
    query_index: dict[tuple[str, int], dict],
    *,
    panel,
    settings: Settings,
    alignment=None,
) -> ComparisonBlock:
    block = compare_populations(
        query_index,
        panel=panel,
        settings=settings,
        groups=settings.caste_groups(),
        kind="caste",
        alignment=alignment,
    )
    block.notes = list(block.notes) + caste_alias_notes(settings.caste_groups())
    return block
