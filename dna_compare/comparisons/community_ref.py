from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path

from dna_compare.config import (
    BENGALI_COMMUNITY_REF,
    COMMUNITY_PANEL_ALIASES,
    GUJARATI_COMMUNITY_REF,
    KERALA_COMMUNITY_REF,
    MARATHI_COMMUNITY_REF,
    PUNJABI_COMMUNITY_REF,
    Settings,
    TAMIL_COMMUNITY_REF,
)
from dna_compare.models import CommunityRefMatch, ComparisonBlock, HaplogroupResult

TAMIL_FILENAME_HINTS = (
    "tamil",
    "vellalar",
    "iyer",
    "iyengar",
    "nadar",
    "thevar",
    "mukkulathor",
    "vanniyar",
    "gounder",
    "paraiyar",
    "parayar",
    "irula",
    "pillai",
    "chennai",
    "madurai",
    "coimbatore",
    "stu",
)

PUNJABI_FILENAME_HINTS = (
    "punjabi",
    "pjl",
    "lahore",
    "amritsar",
    "khatri",
    "arora",
    "jat",
    "jatt",
    "sikh",
    "panjab",
    "punjab",
)

BENGALI_FILENAME_HINTS = (
    "bengali",
    "beb",
    "bangladesh",
    "dhaka",
    "kolkata",
    "calcutta",
    "kayastha",
    "baidya",
    "namasudra",
    "rajbanshi",
    "santhal",
)

GUJARATI_FILENAME_HINTS = (
    "gujarati",
    "gujarat",
    "gih",
    "ahmedabad",
    "surat",
    "vadodara",
    "rajkot",
    "anavil",
    "nagar",
    "patel",
    "kadva",
    "leva",
    "bania",
    "oswal",
    "kapol",
    "koli",
    "parsi",
)

TAMIL_NOTES = [
    "Tamil fit uses qpAdm-style AASI_Onge and Steppe_MLBA against published Tamil-community ranges.",
    "The table's AASI/ASI column is not the same quantity as AASI_Onge (Onge is a hunter-gatherer proxy; ASI often includes Iran/IVC).",
    "ANI/Steppe is compared to Steppe_MLBA (Sintashta). Inside the published range scores higher.",
    "Y percentages in the table are community frequencies. One person has one Y haplogroup; an uncalled or conflicting Y is left out of the fit.",
    "This is a reference-table overlay, not a caste assignment and not an AADR HO bar.",
]

PUNJABI_NOTES = [
    "Punjabi fit compares AASI_Onge, Steppe_MLBA, and Indus_Periphery to published community ranges.",
    "Table ASI/AASI and Iranian-farmer columns map to AASI_Onge and Indus_Periphery in the 3-source model — not identical labels.",
    "ANI (total) from the source table is not scored directly; Steppe and Indus ranges carry most of the north/west Asian signal.",
    "This is a reference-table overlay, not a caste assignment and not an AADR HO bar.",
]

BENGALI_NOTES = [
    "Bengali fit compares ASI to AASI_Onge, ANI to Steppe_MLBA, and East Asian % to East_Asian (Dai) from the 5-source model.",
    "These published labels are not identical to qpAdm source names; treat as a research overlay only.",
    "Y percentages are community frequencies from the source table; one person has one Y haplogroup.",
    "This is a reference-table overlay, not a caste assignment and not an AADR HO bar.",
]

GUJARATI_NOTES = [
    "Gujarati fit compares ASI to AASI_Onge and ANI to Steppe_MLBA from the 3-source model.",
    "Table ANI/ASI columns are not identical to qpAdm source names; treat as a research overlay only.",
    "Parsi ANI totals include substantial Iranian ancestry (~50%) within the published ANI range.",
    "This is a reference-table overlay, not a caste assignment and not an AADR HO bar.",
]

MARATHI_NOTES = [
    "Marathi fit compares ASI to AASI_Onge and ANI to Steppe_MLBA from the 3-source model.",
    "There is no Maratha / Marathi HO bar; panels trigger from filename hints, aliases, or HO Brahmin plus a Maharashtra filename hint.",
    "HO Brahmin (two samples) is not Deshastha, Chitpavan, or Karhade — use this table for Maharashtra Brahmin ranges.",
    "Y percentages are community frequencies from the source table; one person has one Y haplogroup.",
    "This is a reference-table overlay, not a caste assignment and not an AADR HO bar.",
]

KERALA_NOTES = [
    "Kerala fit compares ASI/AASI to AASI_Onge and ANI to Steppe_MLBA from the 3-source model.",
    "There is no Nair / Ezhava / Namboothiri HO bar; panels trigger from filename hints, aliases, or HO Brahmin plus a Kerala filename hint.",
    "Cochin Jewish (Jew_Cochin) is a separate small HO set — not Syrian Christian, Nair, or general Malayali identity.",
    "Y percentages are community frequencies from the source table; one person has one Y haplogroup.",
    "This is a reference-table overlay, not a caste assignment and not an AADR HO bar.",
]

_PANEL_ALIASES_CACHE: dict[str, tuple[str, ...]] | None = None

_Y_ALIASES = (
    ("R1a", "R1a"),
    ("H-M69", "H-M69"),
    ("H (M69)", "H-M69"),
    ("H1", "H-M69"),
    ("H3", "H-M69"),
    ("L-M20", "L-M20"),
    ("L (M20)", "L-M20"),
    ("L1", "L-M20"),
    ("L3", "L-M20"),
    ("J2", "J2"),
    ("R2", "R2"),
    ("C-M130", "C-M130"),
    ("C (M130)", "C-M130"),
    ("O-M175", "O-M175"),
    ("O2a", "O-M175"),
    ("O2", "O-M175"),
    ("O", "O-M175"),
)


@dataclass(frozen=True)
class CommunityRefRow:
    community_id: str
    display_name: str
    aasi_min: float
    aasi_max: float
    steppe_min: float
    steppe_max: float
    y_haplogroups: dict[str, float]
    note: str
    indus_min: float | None = None
    indus_max: float | None = None
    east_asian_min: float | None = None
    east_asian_max: float | None = None


@dataclass(frozen=True)
class CommunityRefPanelSpec:
    panel_id: str
    title: str
    path: Path
    notes: tuple[str, ...]
    caste_label: str
    filename_hints: tuple[str, ...]
    caste_labels: tuple[str, ...] = ()
    proxy_caste_labels: tuple[str, ...] = ()
    proxy_region_hints: tuple[str, ...] = ()
    ref_aasi_label: str = "Ref AASI"
    ref_steppe_label: str = "Ref Steppe"
    ref_indus_label: str = "Ref Indus"
    ref_east_asian_label: str = ""
    sample_aasi_label: str = "This AASI"
    sample_steppe_label: str = "This Steppe"
    sample_indus_label: str = "This Indus"
    sample_east_asian_label: str = ""


COMMUNITY_REF_PANELS: tuple[CommunityRefPanelSpec, ...] = (
    CommunityRefPanelSpec(
        panel_id="tamil",
        title="Tamil community reference ranges",
        path=TAMIL_COMMUNITY_REF,
        notes=tuple(TAMIL_NOTES),
        caste_label="Tamil",
        filename_hints=TAMIL_FILENAME_HINTS,
        caste_labels=("Tamil", "Vellalar", "Irula"),
    ),
    CommunityRefPanelSpec(
        panel_id="punjabi",
        title="Punjabi community reference ranges",
        path=PUNJABI_COMMUNITY_REF,
        notes=tuple(PUNJABI_NOTES),
        caste_label="Punjabi",
        filename_hints=PUNJABI_FILENAME_HINTS,
    ),
    CommunityRefPanelSpec(
        panel_id="bengali",
        title="Bengali community reference ranges",
        path=BENGALI_COMMUNITY_REF,
        notes=tuple(BENGALI_NOTES),
        caste_label="Bengali",
        filename_hints=BENGALI_FILENAME_HINTS,
        ref_aasi_label="Ref ASI",
        ref_steppe_label="Ref ANI",
        sample_aasi_label="This ASI",
        sample_steppe_label="This ANI",
        ref_east_asian_label="Ref E Asian",
        sample_east_asian_label="This E Asian",
    ),
    CommunityRefPanelSpec(
        panel_id="gujarati",
        title="Gujarati community reference ranges",
        path=GUJARATI_COMMUNITY_REF,
        notes=tuple(GUJARATI_NOTES),
        caste_label="Gujarati",
        filename_hints=GUJARATI_FILENAME_HINTS,
        ref_aasi_label="Ref ASI",
        ref_steppe_label="Ref ANI",
        sample_aasi_label="This ASI",
        sample_steppe_label="This ANI",
    ),
    CommunityRefPanelSpec(
        panel_id="marathi",
        title="Marathi community reference ranges",
        path=MARATHI_COMMUNITY_REF,
        notes=tuple(MARATHI_NOTES),
        caste_label="",
        filename_hints=(),
        proxy_caste_labels=("Brahmin",),
        proxy_region_hints=(
            "maratha",
            "marathi",
            "maharashtra",
            "pune",
            "mumbai",
            "nagpur",
            "kolhapur",
            "satara",
            "chitpavan",
            "ckp",
            "deshastha",
            "karhade",
            "kunbi",
            "dhangar",
            "mahar",
        ),
        ref_aasi_label="Ref ASI",
        ref_steppe_label="Ref ANI",
        sample_aasi_label="This ASI",
        sample_steppe_label="This ANI",
    ),
    CommunityRefPanelSpec(
        panel_id="kerala",
        title="Kerala community reference ranges",
        path=KERALA_COMMUNITY_REF,
        notes=tuple(KERALA_NOTES),
        caste_label="",
        filename_hints=(),
        proxy_caste_labels=("Brahmin",),
        proxy_region_hints=(
            "kerala",
            "malayalam",
            "malayali",
            "malabar",
            "kochi",
            "cochin",
            "trivandrum",
            "thiruvananthapuram",
            "kozhikode",
            "calicut",
            "nair",
            "nayar",
            "ezhava",
            "thiyya",
            "namboothiri",
            "nambudiri",
            "pulaya",
            "nasrani",
            "thrissur",
            "kannur",
        ),
        ref_aasi_label="Ref ASI",
        ref_steppe_label="Ref ANI",
        sample_aasi_label="This ASI",
        sample_steppe_label="This ANI",
    ),
)


def _panel_by_id(panel_id: str) -> CommunityRefPanelSpec:
    for panel in COMMUNITY_REF_PANELS:
        if panel.panel_id == panel_id:
            return panel
    raise KeyError(panel_id)


def load_community_panel_aliases(path: Path | None = None) -> dict[str, tuple[str, ...]]:
    global _PANEL_ALIASES_CACHE
    src = path or COMMUNITY_PANEL_ALIASES
    if path is None and _PANEL_ALIASES_CACHE is not None:
        return _PANEL_ALIASES_CACHE
    out: dict[str, list[str]] = {}
    if src.exists():
        with src.open(encoding="utf-8") as handle:
            header_line = handle.readline()
            if header_line:
                header = [col.strip() for col in header_line.rstrip("\n").split("\t")]
                for line in handle:
                    if not line.strip() or line.startswith("#"):
                        continue
                    cols = line.rstrip("\n").split("\t")
                    rec = {header[i]: cols[i].strip() if i < len(cols) else "" for i in range(len(header))}
                    panel_id = rec.get("panel_id") or ""
                    alias = (rec.get("alias") or "").lower()
                    if panel_id and alias:
                        out.setdefault(panel_id, []).append(alias)
    merged = {panel_id: tuple(sorted(set(aliases))) for panel_id, aliases in out.items()}
    if path is None:
        _PANEL_ALIASES_CACHE = merged
    return merged


def _panel_filename_hints(panel: CommunityRefPanelSpec, aliases: dict[str, tuple[str, ...]] | None = None) -> tuple[str, ...]:
    alias_map = aliases if aliases is not None else load_community_panel_aliases()
    extra = alias_map.get(panel.panel_id, ())
    return panel.filename_hints + extra


def _panel_caste_labels(panel: CommunityRefPanelSpec) -> tuple[str, ...]:
    if panel.caste_labels:
        return panel.caste_labels
    if panel.caste_label:
        return (panel.caste_label,)
    return ()


def collapse_ref_y(label: str | None) -> str | None:
    text = (label or "").strip()
    if not text or text.lower() == "other":
        return None
    for prefix, key in _Y_ALIASES:
        if text == prefix or text.startswith(prefix):
            return key
    if text.startswith("H"):
        return "H-M69"
    if text.startswith("L"):
        return "L-M20"
    return None


def interval_score(value: float | None, lo: float, hi: float) -> float | None:
    if value is None:
        return None
    if lo <= value <= hi:
        return 1.0
    dist = lo - value if value < lo else value - hi
    return max(0.0, 1.0 - dist / 15.0)


def parse_y_haplogroups(raw: str) -> dict[str, float]:
    out: dict[str, float] = {}
    for part in (raw or "").split(";"):
        token = part.strip()
        if not token or ":" not in token:
            continue
        name, pct = token.split(":", 1)
        key = collapse_ref_y(name.strip()) or name.strip()
        out[key] = float(pct)
    return out


def _collapse_uniform_notes(rows: list[CommunityRefRow]) -> list[CommunityRefRow]:
    """Drop per-row notes when every populated row repeats the same text."""
    populated = {row.note.strip() for row in rows if row.note.strip()}
    if len(populated) != 1:
        return rows
    return [replace(row, note="") for row in rows]


def _float_col(rec: dict[str, str], key: str) -> float | None:
    raw = (rec.get(key) or "").strip()
    if not raw:
        return None
    try:
        return float(raw)
    except ValueError:
        return None


def load_community_reference(path: Path | None = None) -> list[CommunityRefRow]:
    src = path or TAMIL_COMMUNITY_REF
    if not src.exists():
        return []
    rows: list[CommunityRefRow] = []
    with src.open(encoding="utf-8") as handle:
        header_line = handle.readline()
        if not header_line:
            return []
        header = [col.strip() for col in header_line.rstrip("\n").split("\t")]
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            rec = {header[i]: cols[i].strip() if i < len(cols) else "" for i in range(len(header))}
            if not rec.get("id") or not rec.get("display_name"):
                continue
            rows.append(
                CommunityRefRow(
                    community_id=rec["id"],
                    display_name=rec["display_name"],
                    aasi_min=float(rec["aasi_min"]),
                    aasi_max=float(rec["aasi_max"]),
                    steppe_min=float(rec["steppe_min"]),
                    steppe_max=float(rec["steppe_max"]),
                    indus_min=_float_col(rec, "indus_min"),
                    indus_max=_float_col(rec, "indus_max"),
                    east_asian_min=_float_col(rec, "east_asian_min"),
                    east_asian_max=_float_col(rec, "east_asian_max"),
                    y_haplogroups=parse_y_haplogroups(rec.get("y_haplogroups") or ""),
                    note=rec.get("note") or "",
                )
            )
    return _collapse_uniform_notes(rows)


def load_tamil_community_reference(path: Path | None = None) -> list[CommunityRefRow]:
    return load_community_reference(path or TAMIL_COMMUNITY_REF)


def load_punjabi_community_reference(path: Path | None = None) -> list[CommunityRefRow]:
    return load_community_reference(path or PUNJABI_COMMUNITY_REF)


def load_bengali_community_reference(path: Path | None = None) -> list[CommunityRefRow]:
    return load_community_reference(path or BENGALI_COMMUNITY_REF)


def load_gujarati_community_reference(path: Path | None = None) -> list[CommunityRefRow]:
    return load_community_reference(path or GUJARATI_COMMUNITY_REF)


def load_marathi_community_reference(path: Path | None = None) -> list[CommunityRefRow]:
    return load_community_reference(path or MARATHI_COMMUNITY_REF)


def load_kerala_community_reference(path: Path | None = None) -> list[CommunityRefRow]:
    return load_community_reference(path or KERALA_COMMUNITY_REF)


def _ancestry_pct(ancestry: ComparisonBlock | None, name: str) -> float | None:
    if ancestry is None:
        return None
    for est in ancestry.estimates:
        label = getattr(est, "population", None) or (est.get("population") if isinstance(est, dict) else None)
        if label == name:
            return float(getattr(est, "percent", est.get("percent") if isinstance(est, dict) else 0))
    return None


def _estimate_label(est) -> str:
    if isinstance(est, dict):
        return str(est.get("population") or est.get("label") or "")
    return str(getattr(est, "population", None) or getattr(est, "label", "") or "")


def _estimate_percent(est) -> float:
    if isinstance(est, dict):
        return float(est.get("percent") or 0)
    return float(getattr(est, "percent", 0) or 0)


def panel_reference_applicable(
    panel: CommunityRefPanelSpec,
    caste: ComparisonBlock | None,
    filename: str | None = None,
    *,
    aliases: dict[str, tuple[str, ...]] | None = None,
) -> bool:
    name = (filename or "").lower()
    if any(hint in name for hint in _panel_filename_hints(panel, aliases=aliases)):
        return True
    if caste is None or not caste.available or not caste.estimates:
        return False
    ranked = sorted(caste.estimates, key=_estimate_percent, reverse=True)
    top3 = {_estimate_label(est) for est in ranked[:3]}
    direct = set(_panel_caste_labels(panel))
    if direct and top3 & direct:
        return True
    if panel.proxy_caste_labels and top3 & set(panel.proxy_caste_labels):
        region_hints = panel.proxy_region_hints or _panel_filename_hints(panel, aliases=aliases)
        return any(hint in name for hint in region_hints)
    return False


def tamil_reference_applicable(
    caste: ComparisonBlock | None,
    filename: str | None = None,
) -> bool:
    return panel_reference_applicable(_panel_by_id("tamil"), caste, filename=filename)


def punjabi_reference_applicable(
    caste: ComparisonBlock | None,
    filename: str | None = None,
) -> bool:
    return panel_reference_applicable(_panel_by_id("punjabi"), caste, filename=filename)


def bengali_reference_applicable(
    caste: ComparisonBlock | None,
    filename: str | None = None,
) -> bool:
    return panel_reference_applicable(_panel_by_id("bengali"), caste, filename=filename)


def gujarati_reference_applicable(
    caste: ComparisonBlock | None,
    filename: str | None = None,
) -> bool:
    return panel_reference_applicable(_panel_by_id("gujarati"), caste, filename=filename)


def marathi_reference_applicable(
    caste: ComparisonBlock | None,
    filename: str | None = None,
) -> bool:
    return panel_reference_applicable(_panel_by_id("marathi"), caste, filename=filename)


def kerala_reference_applicable(
    caste: ComparisonBlock | None,
    filename: str | None = None,
) -> bool:
    return panel_reference_applicable(_panel_by_id("kerala"), caste, filename=filename)


def _panel_path(panel: CommunityRefPanelSpec, settings: Settings | None) -> Path:
    if settings is None:
        return panel.path
    if panel.panel_id == "tamil":
        return settings.tamil_community_ref
    if panel.panel_id == "punjabi":
        return settings.punjabi_community_ref
    if panel.panel_id == "bengali":
        return settings.bengali_community_ref
    if panel.panel_id == "gujarati":
        return settings.gujarati_community_ref
    if panel.panel_id == "marathi":
        return settings.marathi_community_ref
    if panel.panel_id == "kerala":
        return settings.kerala_community_ref
    return panel.path


def _score_panel(
    panel: CommunityRefPanelSpec,
    ancestry: ComparisonBlock | None,
    ancestry_5: ComparisonBlock | None,
    haplogroups: HaplogroupResult | None,
    *,
    path: Path,
) -> ComparisonBlock:
    refs = load_community_reference(path)
    notes = list(panel.notes)
    if not refs:
        return ComparisonBlock(
            kind="community_ref",
            available=False,
            notes=notes + [f"{panel.title}: reference table not found at {path}."],
        )
    sample_aasi = _ancestry_pct(ancestry, "AASI_Onge")
    sample_steppe = _ancestry_pct(ancestry, "Steppe_MLBA")
    sample_indus = _ancestry_pct(ancestry, "Indus_Periphery")
    sample_east_asian = _ancestry_pct(ancestry_5, "East_Asian")
    uses_indus = any(row.indus_min is not None and row.indus_max is not None for row in refs)
    uses_east_asian = any(
        row.east_asian_min is not None and row.east_asian_max is not None for row in refs
    )
    if (
        sample_aasi is None
        and sample_steppe is None
        and (not uses_indus or sample_indus is None)
        and (not uses_east_asian or sample_east_asian is None)
    ):
        return ComparisonBlock(
            kind="community_ref",
            available=False,
            notes=notes + [f"{panel.title}: need the ancestry comparison to score this table."],
        )
    sample_y = haplogroups.sample_best if haplogroups is not None else None
    sample_y_key = collapse_ref_y(sample_y)
    matches: list[CommunityRefMatch] = []
    n_snps = 0
    if ancestry is not None:
        for est in ancestry.estimates:
            n_snps = int(getattr(est, "n_snps", 0) or 0)
            if n_snps:
                break
    for row in refs:
        aasi_score = interval_score(sample_aasi, row.aasi_min, row.aasi_max)
        steppe_score = interval_score(sample_steppe, row.steppe_min, row.steppe_max)
        indus_score = None
        if row.indus_min is not None and row.indus_max is not None:
            indus_score = interval_score(sample_indus, row.indus_min, row.indus_max)
        east_asian_score = None
        if row.east_asian_min is not None and row.east_asian_max is not None:
            east_asian_score = interval_score(sample_east_asian, row.east_asian_min, row.east_asian_max)
        y_score = None
        y_note = "Y not used (no derived backbone call in this VCF)."
        if sample_y_key and row.y_haplogroups:
            typical = row.y_haplogroups.get(sample_y_key)
            if typical is not None:
                y_score = min(1.0, typical / 25.0)
                y_note = f"Sample {sample_y} maps to {sample_y_key} (table ~{typical:.0f}%)."
            else:
                y_score = 0.0
                y_note = f"Sample {sample_y} is not among this row's listed Y haplogroups."
        parts: list[tuple[float, float]] = []
        if steppe_score is not None:
            if uses_indus:
                steppe_w = 0.45
            elif uses_east_asian:
                steppe_w = 0.3
            else:
                steppe_w = 0.6 if y_score is None else 0.5
            parts.append((steppe_score, steppe_w))
        if aasi_score is not None:
            if uses_indus:
                aasi_w = 0.35
            elif uses_east_asian:
                aasi_w = 0.3
            else:
                aasi_w = 0.4 if y_score is None else 0.35
            parts.append((aasi_score, aasi_w))
        if indus_score is not None:
            parts.append((indus_score, 0.25 if y_score is None else 0.2))
        if east_asian_score is not None:
            parts.append((east_asian_score, 0.25 if y_score is None else 0.2))
        if y_score is not None:
            parts.append((y_score, 0.15))
        weight = sum(w for _s, w in parts) or 1.0
        fit = 100.0 * sum(score * w for score, w in parts) / weight
        indus_in_range = bool(
            sample_indus is not None
            and row.indus_min is not None
            and row.indus_max is not None
            and row.indus_min <= sample_indus <= row.indus_max
        )
        east_asian_in_range = bool(
            sample_east_asian is not None
            and row.east_asian_min is not None
            and row.east_asian_max is not None
            and row.east_asian_min <= sample_east_asian <= row.east_asian_max
        )
        matches.append(
            CommunityRefMatch(
                population=row.display_name,
                percent=round(fit, 1),
                n_snps=n_snps,
                panel=panel.title,
                ref_aasi_label=panel.ref_aasi_label,
                ref_steppe_label=panel.ref_steppe_label,
                ref_indus_label=panel.ref_indus_label,
                ref_east_asian_label=panel.ref_east_asian_label,
                sample_aasi_label=panel.sample_aasi_label,
                sample_steppe_label=panel.sample_steppe_label,
                sample_indus_label=panel.sample_indus_label,
                sample_east_asian_label=panel.sample_east_asian_label,
                aasi_in_range=bool(sample_aasi is not None and row.aasi_min <= sample_aasi <= row.aasi_max),
                steppe_in_range=bool(
                    sample_steppe is not None and row.steppe_min <= sample_steppe <= row.steppe_max
                ),
                indus_in_range=indus_in_range,
                east_asian_in_range=east_asian_in_range,
                aasi_score=None if aasi_score is None else round(aasi_score, 3),
                steppe_score=None if steppe_score is None else round(steppe_score, 3),
                indus_score=None if indus_score is None else round(indus_score, 3),
                east_asian_score=None if east_asian_score is None else round(east_asian_score, 3),
                sample_aasi=sample_aasi,
                sample_steppe=sample_steppe,
                sample_indus=sample_indus,
                sample_east_asian=sample_east_asian,
                ref_aasi=f"{row.aasi_min:.0f}–{row.aasi_max:.0f}%",
                ref_steppe=f"{row.steppe_min:.0f}–{row.steppe_max:.0f}%",
                ref_indus=(
                    f"{row.indus_min:.0f}–{row.indus_max:.0f}%"
                    if row.indus_min is not None and row.indus_max is not None
                    else ""
                ),
                ref_east_asian=(
                    f"{row.east_asian_min:.0f}–{row.east_asian_max:.0f}%"
                    if row.east_asian_min is not None and row.east_asian_max is not None
                    else ""
                ),
                ref_y="; ".join(f"{name} ~{pct:.0f}%" for name, pct in row.y_haplogroups.items()),
                sample_y=sample_y,
                y_note=y_note,
                note=row.note,
            )
        )
    matches.sort(key=lambda item: item.percent, reverse=True)
    if matches:
        notes.append(
            f"{panel.title}: closest range {matches[0].population} (fit {matches[0].percent:.0f}%). "
            "Fit is distance to the table, not a probability you belong to that community."
        )
    return ComparisonBlock(kind="community_ref", available=True, estimates=matches, notes=notes)


def score_community_reference(
    ancestry: ComparisonBlock | None,
    haplogroups: HaplogroupResult | None,
    *,
    ancestry_5: ComparisonBlock | None = None,
    caste: ComparisonBlock | None = None,
    filename: str | None = None,
    settings: Settings | None = None,
) -> ComparisonBlock:
    applicable = [
        panel
        for panel in COMMUNITY_REF_PANELS
        if panel_reference_applicable(panel, caste, filename=filename)
    ]
    if not applicable:
        return ComparisonBlock(
            kind="community_ref",
            available=False,
            hidden=True,
            notes=[
                "Community reference ranges hidden: no matching panel "
                f"({', '.join(panel.title for panel in COMMUNITY_REF_PANELS)} — "
                "triggered by top community weights, filename hints, or aliases under data/references/caste/community_panel_aliases.tsv)."
            ],
        )
    combined_estimates: list[CommunityRefMatch] = []
    combined_notes: list[str] = []
    any_available = False
    for panel in applicable:
        path = _panel_path(panel, settings)
        result = _score_panel(panel, ancestry, ancestry_5, haplogroups, path=path)
        combined_notes.extend(result.notes)
        if result.available:
            any_available = True
            combined_estimates.extend(result.estimates)
    if not any_available:
        return ComparisonBlock(
            kind="community_ref",
            available=False,
            notes=combined_notes,
        )
    return ComparisonBlock(
        kind="community_ref",
        available=True,
        estimates=combined_estimates,
        notes=combined_notes,
    )
