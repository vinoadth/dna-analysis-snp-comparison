from __future__ import annotations

import re
from collections import defaultdict
from functools import lru_cache
from pathlib import Path

from dna_compare.comparisons.disease import _match_dosage
from dna_compare.config import Settings
from dna_compare.models import AdditionalDetail, ComparisonBlock
from dna_compare.vcf_parser import normalize_chrom

_RSID_RE = re.compile(r"(rs\d+)", re.I)
_PARENT_MOTHER = ("mother", "mom", "maternal")
_PARENT_FATHER = ("father", "dad", "paternal")

ADDITIONAL_NOTES = [
    "Research overlay from a few published markers — not a diagnosis, not diet advice, not a blood-bank type, and not a caste or community ID.",
    "Coverage: ready = all listed markers called on this file; partial = some; missing = none (common on sparse or imputed-only files).",
    "ABO and Rh rows are research guesses, not blood-bank typing. Rh prefers 23andMe i4001527 (RHD deletion) when present, else rs590787 proxy — less reliable outside European-ancestry deletion backgrounds.",
    "ABO, HFE, MTHFR, clotting, and warfarin rows are short genotype tables. They are not clinical tests.",
    "Vegetarian-community rows use chr 11 FADS1/FADS2 SNPs that tag the plant-diet / South Asian insertion haplotype (Kothapalli 2016).",
    "Those FADS SNPs are autosomal: a heterozygous call means one parent contributed the haplotype. Mother vs father is named only when a parental VCF is attached.",
    "Add later traits as extra keys in data/references/additional/details.tsv; they appear as new rows in this table.",
]

DISPLAY_ORDER = (
    "abo_blood_type",
    "rh_factor",
    "blood_type",
    "vegetarian_community",
    "heart_disease",
    "celiac",
    "ibd",
    "autoimmune_ptpn22",
    "hfe",
    "clotting",
    "mthfr",
    "warfarin",
    "cyp2c19",
    "tpmt",
    "slco1b1",
    "nat2",
    "alcohol",
    "lactase",
    "pigmentation",
    "eye_color",
    "mc1r",
    "hair_color",
    "nicotine",
    "obesity",
    "t2d_single",
    "apoe",
    "cilantro",
    "photic_sneeze",
    "motion_sickness",
    "bitter_taste",
    "actn3",
    "baldness",
    "vitamin_d",
    "sleep_preferences",
    "earwax",
    "asparagus",
    "caffeine",
)


def _norm_rsid(rsid: str) -> str:
    text = (rsid or "").split(",")[0].strip()
    text = re.sub(r"^(gsa-|exm-|ilm-)", "", text, flags=re.I)
    match = _RSID_RE.search(text)
    return match.group(1).lower() if match else text.lower()


def parent_role(filename: str | None) -> str | None:
    name = Path(filename or "").name.lower()
    if any(token in name for token in _PARENT_MOTHER):
        return "mother"
    if any(token in name for token in _PARENT_FATHER):
        return "father"
    return None


@lru_cache(maxsize=4)
def load_additional_markers(path: Path) -> dict[str, list[dict]]:
    grouped: dict[str, list[dict]] = defaultdict(list)
    if not path.is_file():
        return {}
    header: list[str] | None = None
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            if raw.startswith("#") or not raw.strip():
                continue
            cols = raw.rstrip("\n").split("\t")
            if header is None:
                header = cols
                continue
            rec = {key: cols[i] if i < len(cols) else "" for i, key in enumerate(header)}
            key = (rec.get("key") or "").strip()
            rsid = (rec.get("rsid") or "").strip()
            if not key or not rsid:
                continue
            try:
                pos_hg19 = int(rec.get("pos_hg19") or 0)
            except ValueError:
                pos_hg19 = 0
            try:
                pos_hg38 = int(rec.get("pos_hg38") or 0)
            except ValueError:
                pos_hg38 = 0
            grouped[key].append(
                {
                    "key": key,
                    "topic": (rec.get("topic") or key).strip(),
                    "chrom": normalize_chrom(rec.get("chrom") or "11"),
                    "rsid": rsid,
                    "pos_hg19": pos_hg19,
                    "pos_hg38": pos_hg38,
                    "effect_allele": (rec.get("effect_allele") or "").strip().upper(),
                    "other_allele": (rec.get("other_allele") or "").strip().upper(),
                    "gene": (rec.get("gene") or "").strip(),
                    "source": (rec.get("source") or "").strip(),
                    "hom_finding": (rec.get("hom_finding") or "").strip(),
                    "het_finding": (rec.get("het_finding") or "").strip(),
                    "ref_finding": (rec.get("ref_finding") or "").strip(),
                }
            )
    return dict(grouped)


def _lookups(index: dict[tuple[str, int], dict]) -> tuple[dict[str, dict], dict[tuple[str, int], dict]]:
    by_rsid: dict[str, dict] = {}
    by_pos: dict[tuple[str, int], dict] = {}
    for key, row in index.items():
        by_pos[key] = row
        rsid = _norm_rsid(str(row.get("rsid") or ""))
        if rsid and rsid not in {".", ""}:
            by_rsid[rsid] = row
    return by_rsid, by_pos


def _find_row(marker: dict, by_rsid: dict[str, dict], by_pos: dict[tuple[str, int], dict]) -> dict | None:
    row = by_rsid.get(_norm_rsid(marker["rsid"]))
    if row is not None:
        return row
    chrom = marker["chrom"]
    for pos in (marker.get("pos_hg19"), marker.get("pos_hg38")):
        if pos:
            row = by_pos.get((chrom, int(pos)))
            if row is not None:
                return row
    return None


def allele_dosage(row: dict | None, effect: str, other: str = "") -> float | None:
    if row is None:
        return None
    got = _match_dosage(row, effect, other)
    if got is not None:
        return float(got)
    ref = str(row.get("ref") or "").upper()
    alt = str(row.get("alt") or "").upper()
    dosage = row.get("dosage_alt")
    if dosage is None:
        return None
    effect = (effect or "").upper()
    other = (other or "").upper()
    aliases = {
        "I": {"I", "INS", "G"},
        "D": {"D", "DEL", "-"},
    }

    def same(token: str, want: str) -> bool:
        return bool(want) and (token == want or token in aliases.get(want, set()))

    if same(alt, effect) and (not other or same(ref, other)):
        return float(dosage)
    if same(ref, effect) and (not other or same(alt, other)):
        return 2.0 - float(dosage)
    return None


def score_marker_group(
    index: dict[tuple[str, int], dict],
    markers: list[dict],
    *,
    lookups: tuple[dict[str, dict], dict[tuple[str, int], dict]] | None = None,
) -> dict:
    if lookups is None:
        by_rsid, by_pos = _lookups(index)
    else:
        by_rsid, by_pos = lookups
    dosages: list[float] = []
    evidence: list[str] = []
    for marker in markers:
        row = _find_row(marker, by_rsid, by_pos)
        if row is None:
            continue
        copies = allele_dosage(row, marker["effect_allele"], marker["other_allele"])
        if copies is None:
            continue
        dosages.append(float(copies))
        gene = f" {marker['gene']}" if marker.get("gene") else ""
        evidence.append(
            f"{marker['rsid']}{gene} {row.get('genotype') or '?'} "
            f"({copies:.0f}× {marker['effect_allele']})"
        )
    mean = sum(dosages) / len(dosages) if dosages else None
    if mean is None:
        status = "missing"
    elif mean >= 1.5:
        status = "both"
    elif mean >= 0.5:
        status = "one"
    else:
        status = "none"
    spread = max(dosages) - min(dosages) if dosages else 0
    if not dosages:
        confidence = "none"
    elif spread <= 1 and len(dosages) >= 3:
        confidence = "high"
    elif spread <= 1:
        confidence = "medium"
    else:
        confidence = "low"
    return {
        "status": status,
        "mean": mean,
        "n_snps": len(dosages),
        "n_markers": len(markers),
        "evidence": evidence,
        "confidence": confidence,
    }


def _opposite_parent(role: str) -> str:
    return "father" if role == "mother" else "mother"


def _vegetarian_finding(child: dict, other: dict | None, other_filename: str | None, markers: list[dict]) -> tuple[str, str | None, str]:
    first = markers[0]
    status = child["status"]
    if status == "missing":
        return (
            "chr 11 FADS rsIDs were not called on this array.",
            None,
            "Need rs174547, rs174570, or linked FADS1/FADS2 SNPs (GSA usually has them; a chr 1–22-only imputed file will not).",
        )
    if status == "both":
        return first["hom_finding"], "both", ""
    if status == "none":
        return first["ref_finding"], "none", ""
    role = parent_role(other_filename)
    other_status = (other or {}).get("status")
    if role and other_status == "both":
        return f"{role.capitalize()} belongs to a vegetarian community.", role, ""
    if role and other_status == "none":
        other_role = _opposite_parent(role)
        return f"{other_role.capitalize()} belongs to a vegetarian community.", other_role, ""
    if role and other_status == "one":
        return (
            f"{role.capitalize()} is heterozygous at these FADS SNPs; one grandparental lineage on that side looks vegetarian-adapted.",
            role,
            "",
        )
    return (
        first["het_finding"],
        "one",
        "Attach a mother or father VCF to name which parent contributed the haplotype.",
    )


def _abo_copies(
    index: dict[tuple[str, int], dict],
    markers: list[dict],
    *,
    lookups: tuple[dict[str, dict], dict[tuple[str, int], dict]] | None = None,
) -> tuple[float | None, float | None, list[str]]:
    if lookups is None:
        by_rsid, by_pos = _lookups(index)
    else:
        by_rsid, by_pos = lookups
    o_copies = None
    a_copies = None
    evidence: list[str] = []
    for marker in markers:
        row = _find_row(marker, by_rsid, by_pos)
        rsid = _norm_rsid(marker["rsid"])
        if rsid == "rs8176719":
            o_copies = allele_dosage(row, "D", "I")
            if o_copies is None:
                o_copies = allele_dosage(row, "D", "G")
            if row is not None and o_copies is not None:
                evidence.append(f"{marker['rsid']} ABO {row.get('genotype') or '?'} ({o_copies:.0f}× O)")
        elif rsid == "rs8176746":
            a_copies = allele_dosage(row, marker["effect_allele"], marker["other_allele"])
            if a_copies is None:
                a_copies = allele_dosage(row, "C", "A")
            if row is not None and a_copies is not None:
                evidence.append(f"{marker['rsid']} ABO {row.get('genotype') or '?'} ({a_copies:.0f}× A)")
    return o_copies, a_copies, evidence


def _abo_finding(
    index: dict[tuple[str, int], dict],
    markers: list[dict],
    *,
    lookups: tuple[dict[str, dict], dict[tuple[str, int], dict]] | None = None,
) -> tuple[str, str | None, str, dict]:
    o_copies, a_copies, evidence = _abo_copies(index, markers, lookups=lookups)
    n = len(evidence)
    score = {
        "status": "missing" if n == 0 else "both",
        "n_snps": n,
        "n_markers": len(markers),
        "evidence": evidence,
        "confidence": "high" if n == 2 else ("medium" if n == 1 else "none"),
    }
    if n == 0:
        return (
            "ABO markers were not called on this file.",
            None,
            "Need rs8176719 (O indel) and rs8176746 (A vs B). GSA usually has both.",
            score,
        )
    if o_copies == 2:
        return "Blood type O.", "O", "", score
    if o_copies == 0 and a_copies == 2:
        return "Blood type A.", "A", "", score
    if o_copies == 0 and a_copies == 0:
        return "Blood type B.", "B", "", score
    if o_copies == 0 and a_copies == 1:
        return "Blood type AB.", "AB", "", score
    if o_copies == 1 and a_copies == 2:
        return "Blood type A.", "A", "", score
    if o_copies == 1 and a_copies == 0:
        return "Blood type B.", "B", "", score
    if o_copies == 1 and a_copies == 1:
        return (
            "Blood type A or B.",
            "A/B",
            "One O allele and a heterozygous A/B SNP — the array is unphased, so A versus B cannot be named.",
            score,
        )
    if o_copies == 2:
        return "Blood type O.", "O", "", score
    if a_copies == 2:
        return "Blood type A or AB.", "A", "O indel was not called.", score
    if a_copies == 0:
        return "Blood type B or O.", "B", "O indel was not called.", score
    return "ABO genotype is incomplete on this file.", None, "", score


_RH_PROXY_CAVEAT = (
    "rs590787 is a proxy near RHD, not direct deletion typing — often ~90% concordant in "
    "European-ancestry samples but less reliable elsewhere. Weak D, partial D, and DEL variants "
    "are not resolved."
)
_RH_DELETION_CAVEAT = (
    "i4001527 assays the common RHD whole-gene deletion on 23andMe arrays. Non-deletion Rh− "
    "mechanisms and weak/partial D are not resolved."
)
_BLOOD_TYPE_CAVEAT = "Research guess from published markers — not for transfusion or clinical use."


def _marker_id(marker: dict) -> str:
    return _norm_rsid(marker["rsid"])


def _rh_deletion_copies(row: dict | None) -> float | None:
    if row is None:
        return None
    copies = allele_dosage(row, "D", "I")
    if copies is not None:
        return copies
    return allele_dosage(row, "D", "G")


def _rh_neg_tag_copies(row: dict | None, marker: dict) -> float | None:
    if row is None:
        return None
    for effect, other in (
        ("C", "T"),  # SNPedia / minus-strand reporting
        ("C", "A"),  # dbSNP plus-strand (ref A, alt C)
        ("G", "A"),  # alternate plus-strand encoding
    ):
        copies = allele_dosage(row, effect, other)
        if copies is not None:
            return copies
    return allele_dosage(row, marker["effect_allele"], marker["other_allele"])


def _rh_call_from_copies(neg_copies: float) -> tuple[str, str, str]:
    if neg_copies >= 1.5:
        return "Likely Rh negative (−).", "−", "both"
    if neg_copies >= 0.5:
        return "Likely Rh positive (+).", "+", "one"
    return "Likely Rh positive (+).", "+", "none"


def _rh_finding(
    index: dict[tuple[str, int], dict],
    markers: list[dict],
    *,
    lookups: tuple[dict[str, dict], dict[tuple[str, int], dict]] | None = None,
) -> tuple[str, str | None, str, dict]:
    if lookups is None:
        by_rsid, by_pos = _lookups(index)
    else:
        by_rsid, by_pos = lookups
    deletion_copies: float | None = None
    proxy_neg_copies: float | None = None
    evidence: list[str] = []
    for marker in markers:
        row = _find_row(marker, by_rsid, by_pos)
        mid = _marker_id(marker)
        if mid == "i4001527":
            deletion_copies = _rh_deletion_copies(row)
            if row is not None and deletion_copies is not None:
                evidence.append(
                    f"{marker['rsid']} RHD {row.get('genotype') or '?'} ({deletion_copies:.0f}× deletion)"
                )
        elif mid == "rs590787":
            proxy_neg_copies = _rh_neg_tag_copies(row, marker)
            if row is not None and proxy_neg_copies is not None:
                evidence.append(
                    f"{marker['rsid']} RHD {row.get('genotype') or '?'} ({proxy_neg_copies:.0f}× Rh− tag)"
                )
    n = len(evidence)
    if n == 0:
        score = {
            "status": "missing",
            "n_snps": 0,
            "n_markers": len(markers),
            "evidence": evidence,
            "confidence": "none",
        }
        return (
            "Rh factor markers were not called on this file.",
            None,
            "Need i4001527 (23andMe RHD deletion) and/or rs590787 (RHD proxy). Often missing from imputed-only files.",
            score,
        )
    del_neg = deletion_copies is not None and deletion_copies >= 1.5
    del_pos = deletion_copies is not None and deletion_copies < 1.5
    proxy_neg = proxy_neg_copies is not None and proxy_neg_copies >= 1.5
    proxy_pos = proxy_neg_copies is not None and proxy_neg_copies < 1.5
    if deletion_copies is not None and proxy_neg_copies is not None and del_neg != proxy_neg:
        score = {
            "status": "one",
            "n_snps": n,
            "n_markers": len(markers),
            "evidence": evidence,
            "confidence": "low",
        }
        return (
            "Rh status is unclear (i4001527 and rs590787 disagree).",
            None,
            "Treat as inconclusive; serologic typing or a clinical RHD assay is needed.",
            score,
        )
    if deletion_copies is not None:
        finding, parent, status = _rh_call_from_copies(deletion_copies)
        message = _RH_DELETION_CAVEAT
        if proxy_neg_copies is not None:
            message += " rs590787 agrees."
        confidence = "high" if proxy_neg_copies is not None else "medium"
    else:
        finding, parent, status = _rh_call_from_copies(proxy_neg_copies or 0.0)
        message = _RH_PROXY_CAVEAT
        confidence = "medium"
    score = {
        "status": status,
        "n_snps": n,
        "n_markers": len(markers),
        "evidence": evidence,
        "confidence": confidence,
    }
    return finding, parent, message, score


def _combined_blood_type(
    abo: AdditionalDetail | None,
    rh: AdditionalDetail | None,
) -> AdditionalDetail | None:
    if abo is None or rh is None or not abo.parent or not rh.parent:
        return None
    if not abo.available or not rh.available:
        return None
    if abo.parent == "A/B":
        finding = f"Blood type A or B{rh.parent}."
        message = abo.message or "ABO A versus B is unphased on this file."
    else:
        finding = f"Blood type {abo.parent}{rh.parent}."
        message = ""
    parts = [part for part in (message, rh.message, _BLOOD_TYPE_CAVEAT) if part]
    confidence = "medium"
    if abo.confidence == "high" and rh.confidence == "high":
        confidence = "high"
    elif abo.confidence == "none" or rh.confidence in {"none", "low"}:
        confidence = "low"
    evidence = "; ".join(filter(None, [abo.evidence, rh.evidence]))
    n_snps = abo.n_snps + rh.n_snps
    n_markers = abo.n_markers + rh.n_markers
    return AdditionalDetail(
        key="blood_type",
        topic="Blood type (ABO + Rh)",
        finding=finding,
        evidence=evidence,
        source=f"{abo.source}; {rh.source}",
        available=True,
        confidence=confidence,
        n_snps=n_snps,
        n_markers=n_markers,
        coverage_status=_coverage_status(n_snps, n_markers),
        parent=f"{abo.parent}{rh.parent}",
        message=" ".join(parts),
    )


def _coverage_status(n_found: int, n_total: int) -> str:
    if n_total == 0 or n_found == 0:
        return "missing"
    if n_found >= n_total:
        return "ready"
    return "partial"


def _generic_finding(score: dict, markers: list[dict]) -> tuple[str, str | None, str]:
    first = markers[0]
    status = score["status"]
    if status == "missing":
        return f"{first['topic']} markers were not called on this file.", None, ""
    if status == "both":
        return first["hom_finding"] or f"{first['topic']}: homozygous for the listed alleles.", "both", ""
    if status == "one":
        return first["het_finding"] or f"{first['topic']}: heterozygous.", "one", ""
    return first["ref_finding"] or f"{first['topic']}: homozygous for the other alleles.", "none", ""


def compare_additional(
    query_index: dict[tuple[str, int], dict],
    *,
    settings: Settings,
    other_index: dict[tuple[str, int], dict] | None = None,
    other_filename: str | None = None,
) -> ComparisonBlock:
    path = Path(settings.additional_details)
    groups = load_additional_markers(path)
    notes = list(ADDITIONAL_NOTES)
    if not groups:
        return ComparisonBlock(
            kind="additional",
            available=False,
            notes=notes + [f"Missing marker table: {path}."],
        )
    estimates: list[AdditionalDetail] = []
    abo_detail: AdditionalDetail | None = None
    rh_detail: AdditionalDetail | None = None
    lookups = _lookups(query_index)
    other_lookups = _lookups(other_index) if other_index else None
    keys = [key for key in DISPLAY_ORDER if key in groups]
    keys.extend(key for key in groups if key not in DISPLAY_ORDER)
    for key in keys:
        if key == "blood_type":
            continue
        markers = groups[key]
        child = score_marker_group(query_index, markers, lookups=lookups)
        other = score_marker_group(other_index, markers, lookups=other_lookups) if other_index else None
        if key == "abo_blood_type":
            finding, parent, message, child = _abo_finding(query_index, markers, lookups=lookups)
        elif key == "rh_factor":
            finding, parent, message, child = _rh_finding(query_index, markers, lookups=lookups)
        elif key == "vegetarian_community":
            finding, parent, message = _vegetarian_finding(child, other, other_filename, markers)
        else:
            from dna_compare.comparisons.additional_combos import COMBO_HANDLERS

            handler = COMBO_HANDLERS.get(key)
            if handler is not None:
                finding, parent, message, child = handler(query_index, markers, lookups=lookups)
            else:
                finding, parent, message = _generic_finding(child, markers)
        n_snps = int(child["n_snps"])
        n_markers = int(child["n_markers"])
        detail = AdditionalDetail(
            key=key,
            topic=markers[0]["topic"],
            finding=finding,
            evidence="; ".join(child["evidence"][:8]),
            source=markers[0]["source"],
            available=child["status"] != "missing",
            confidence=child["confidence"],
            n_snps=n_snps,
            n_markers=n_markers,
            coverage_status=_coverage_status(n_snps, n_markers),
            parent=parent,
            message=message,
        )
        estimates.append(detail)
        if key == "abo_blood_type":
            abo_detail = detail
        elif key == "rh_factor":
            rh_detail = detail
            combined = _combined_blood_type(abo_detail, rh_detail)
            if combined is not None:
                estimates.append(combined)
    n_callable = sum(1 for item in estimates if item.available)
    n_ready = sum(1 for item in estimates if item.coverage_status == "ready")
    notes.append(
        f"Marker coverage on this file: {n_ready} ready, {n_callable} callable, "
        f"{len(estimates) - n_callable} missing."
    )
    return ComparisonBlock(
        kind="additional",
        available=any(item.available for item in estimates),
        estimates=estimates,
        notes=notes,
    )
