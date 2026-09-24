from __future__ import annotations

from dna_compare.config import Settings, default_settings


def collect_population_groups(settings: Settings | None = None) -> list[tuple[str, ...]]:
    """Population tuples whose AADR allele-frequency caches analysis may need."""
    settings = settings or default_settings()
    groups: list[tuple[str, ...]] = []
    for mapping in (
        settings.selected_population_groups(),
        settings.caste_groups(),
        settings.ancestry_groups(),
        settings.ancestry_5_groups(),
    ):
        groups.extend(tuple(pops) for pops in mapping.values())
    groups.extend((pop,) for pop in ("Mbuti", "Yoruba"))
    return groups
