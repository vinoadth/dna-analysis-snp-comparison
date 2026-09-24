from __future__ import annotations

from dna_compare.comparisons.additional import allele_dosage, _find_row, _lookups, _norm_rsid


def _pack(n: int, n_markers: int, evidence: list[str], *, missing: bool = False) -> dict:
    return {
        "status": "missing" if missing or n == 0 else "both",
        "n_snps": n,
        "n_markers": n_markers,
        "evidence": evidence,
        "confidence": "none" if n == 0 else ("high" if n >= 2 else "medium"),
    }


def rsid_copies(
    index: dict[tuple[str, int], dict],
    markers: list[dict],
    *,
    lookups: tuple[dict[str, dict], dict[tuple[str, int], dict]] | None = None,
) -> tuple[dict[str, float], list[str]]:
    if lookups is None:
        by_rsid, by_pos = _lookups(index)
    else:
        by_rsid, by_pos = lookups
    copies: dict[str, float] = {}
    evidence: list[str] = []
    for marker in markers:
        row = _find_row(marker, by_rsid, by_pos)
        dosage = allele_dosage(row, marker["effect_allele"], marker["other_allele"])
        if row is None or dosage is None:
            continue
        rsid = _norm_rsid(marker["rsid"])
        copies[rsid] = float(dosage)
        gene = f" {marker['gene']}" if marker.get("gene") else ""
        evidence.append(
            f"{marker['rsid']}{gene} {row.get('genotype') or '?'} "
            f"({dosage:.0f}× {marker['effect_allele']})"
        )
    return copies, evidence


def _has_y(index: dict[tuple[str, int], dict]) -> bool:
    return any(chrom == "Y" for chrom, _pos in index)


def _get(copies: dict[str, float], rsid: str) -> float:
    return copies.get(rsid, 0.0)


def hfe_finding(index, markers, *, lookups=None):
    copies, evidence = rsid_copies(index, markers, lookups=lookups)
    score = _pack(len(copies), len(markers), evidence, missing=not copies)
    if not copies:
        return "HFE markers were not called on this file.", None, "", score
    y, d, s = _get(copies, "rs1800562"), _get(copies, "rs1799945"), _get(copies, "rs1800730")
    if y >= 2:
        return "C282Y homozygous — HFE haemochromatosis-risk genotype (research only, not a diagnosis).", "C282Y/C282Y", "", score
    if y >= 1 and d >= 1:
        return "C282Y/H63D compound heterozygote — milder HFE risk genotype.", "C282Y/H63D", "", score
    if y >= 1 and s >= 1:
        return "C282Y/S65C compound heterozygote — milder HFE risk genotype.", "C282Y/S65C", "", score
    if d >= 2:
        return "H63D homozygous — low-penetrance HFE genotype.", "H63D/H63D", "", score
    if y >= 1:
        return "C282Y carrier.", "C282Y carrier", "", score
    if d >= 1 or s >= 1:
        return "H63D or S65C carrier.", "carrier", "", score
    return "No common HFE iron-overload variants.", "wild-type", "", score


def clotting_finding(index, markers, *, lookups=None):
    copies, evidence = rsid_copies(index, markers, lookups=lookups)
    score = _pack(len(copies), len(markers), evidence, missing=not copies)
    if not copies:
        return "Clotting markers were not called on this file.", None, "", score
    fvl, pt = _get(copies, "rs6025"), _get(copies, "rs1799963")
    if fvl >= 2 or pt >= 2 or (fvl >= 1 and pt >= 1):
        return "Higher-risk clotting-variant combination (research only, not a diagnosis).", "high", "", score
    if fvl >= 1:
        return "Factor V Leiden carrier.", "FVL", "", score
    if pt >= 1:
        return "Prothrombin G20210A carrier.", "F2", "", score
    return "No Factor V Leiden or prothrombin G20210A variant.", "wild-type", "", score


def mthfr_finding(index, markers, *, lookups=None):
    copies, evidence = rsid_copies(index, markers, lookups=lookups)
    score = _pack(len(copies), len(markers), evidence, missing=not copies)
    if not copies:
        return "MTHFR markers were not called on this file.", None, "", score
    t677, c1298 = _get(copies, "rs1801133"), _get(copies, "rs1801131")
    if t677 >= 2:
        return "MTHFR 677 TT — reduced-activity genotype (research only, not a diagnosis).", "677TT", "", score
    if t677 >= 1 and c1298 >= 1:
        return "MTHFR 677 CT and 1298 AC — compound variant genotype.", "compound", "", score
    if c1298 >= 2:
        return "MTHFR 1298 CC — reduced-activity genotype.", "1298CC", "", score
    if t677 >= 1 or c1298 >= 1:
        return "One MTHFR variant — typical or mildly reduced activity.", "carrier", "", score
    return "Typical MTHFR genotypes at C677T and A1298C.", "wild-type", "", score


def warfarin_finding(index, markers, *, lookups=None):
    copies, evidence = rsid_copies(index, markers, lookups=lookups)
    score = _pack(len(copies), len(markers), evidence, missing=not copies)
    if not copies:
        return "Warfarin-sensitivity markers were not called on this file.", None, "", score
    sens = _get(copies, "rs9923231") + _get(copies, "rs1799853") + _get(copies, "rs1057910")
    if sens >= 2:
        return "Increased warfarin-sensitivity (lower-dose) genotype at VKORC1/CYP2C9.", "sensitive", "", score
    if sens >= 1:
        return "Intermediate warfarin-sensitivity genotype at VKORC1/CYP2C9.", "intermediate", "", score
    return "Typical warfarin-sensitivity genotype at VKORC1/CYP2C9.", "typical", "", score


def alcohol_finding(index, markers, *, lookups=None):
    copies, evidence = rsid_copies(index, markers, lookups=lookups)
    score = _pack(len(copies), len(markers), evidence, missing=not copies)
    if not copies:
        return "Alcohol-metabolism markers were not called on this file.", None, "", score
    aldh = _get(copies, "rs671")
    adh = _get(copies, "rs1229984")
    if aldh >= 1:
        return "ALDH2 deficient — alcohol flush / acetaldehyde buildup is more likely.", "ALDH2", "", score
    if adh >= 1:
        return "Faster ADH1B (*2) without ALDH2 flush — ethanol is converted faster to acetaldehyde.", "ADH1B", "", score
    return "Typical ADH1B / ALDH2 alcohol-metabolism genotype.", "typical", "", score


def pigmentation_finding(index, markers, *, lookups=None):
    copies, evidence = rsid_copies(index, markers, lookups=lookups)
    score = _pack(len(copies), len(markers), evidence, missing=not copies)
    if not copies:
        return "Pigmentation markers were not called on this file.", None, "", score
    light = sum(copies.values())
    frac = light / (2 * len(copies))
    if frac >= 0.6:
        return "Lighter / more European-like pigmentation alleles.", "light", "", score
    if frac >= 0.25:
        return "Mixed pigmentation alleles.", "mixed", "", score
    return "Ancestral pigmentation haplotype (common in South Asia).", "ancestral", "", score


def eye_color_finding(index, markers, *, lookups=None):
    copies, evidence = rsid_copies(index, markers, lookups=lookups)
    score = _pack(len(copies), len(markers), evidence, missing=not copies)
    if not copies:
        return "Eye-color markers were not called on this file.", None, "", score
    blue = _get(copies, "rs12913832")
    if blue >= 2:
        return "Likely blue or green eyes.", "blue", "", score
    if blue >= 1:
        return "Intermediate eye-color genotype; brown or hazel is more likely.", "intermediate", "", score
    return "Likely brown eyes.", "brown", "", score


def mc1r_finding(index, markers, *, lookups=None):
    copies, evidence = rsid_copies(index, markers, lookups=lookups)
    score = _pack(len(copies), len(markers), evidence, missing=not copies)
    if not copies:
        return "MC1R markers were not called on this file.", None, "", score
    n = sum(copies.values())
    if n >= 2:
        return "Two or more MC1R red-hair variants — fair skin / red hair is more likely.", "red", "", score
    if n >= 1:
        return "One MC1R variant — typical hair color; slightly fairer skin is possible.", "carrier", "", score
    return "No common MC1R red-hair variants.", "none", "", score


def bitter_diplotype_finding(index, markers, *, lookups=None):
    copies, evidence = rsid_copies(index, markers, lookups=lookups)
    score = _pack(len(copies), len(markers), evidence, missing=not copies)
    if not copies:
        return "TAS2R38 markers were not called on this file.", None, "", score
    taster = [_get(copies, r) for r in ("rs713598", "rs1726866", "rs10246939") if r in copies]
    if taster and all(v >= 1.5 for v in taster):
        return "TAS2R38 PAV/PAV — likely a bitter (PTC) taster.", "PAV/PAV", "", score
    if taster and all(v <= 0.5 for v in taster):
        return "TAS2R38 AVI/AVI — likely a bitter (PTC) non-taster.", "AVI/AVI", "", score
    return "TAS2R38 PAV/AVI or mixed sites — intermediate bitter taste.", "PAV/AVI", "", score


def actn3_ace_finding(index, markers, *, lookups=None):
    copies, evidence = rsid_copies(index, markers, lookups=lookups)
    score = _pack(len(copies), len(markers), evidence, missing=not copies)
    if not copies:
        return "ACTN3/ACE markers were not called on this file.", None, "", score
    r_copies = _get(copies, "rs1815739")
    d_copies = _get(copies, "rs4343")
    if "rs4343" not in copies:
        d_copies = _get(copies, "rs4341")
    actn = "RR" if r_copies >= 1.5 else ("XX" if r_copies <= 0.5 and "rs1815739" in copies else "RX")
    ace = "DD" if d_copies >= 1.5 else ("II" if d_copies <= 0.5 else "ID")
    if actn == "RR" and ace == "DD":
        return "Sprint/power lean (ACTN3 RR + ACE DD).", "power", "", score
    if actn == "XX" and ace == "II":
        return "Endurance lean (ACTN3 XX + ACE II).", "endurance", "", score
    return f"Mixed sprint and endurance (ACTN3 {actn} + ACE {ace}).", "mixed", "", score


def baldness_finding(index, markers, *, lookups=None):
    copies, evidence = rsid_copies(index, markers, lookups=lookups)
    score = _pack(len(copies), len(markers), evidence, missing=not copies)
    if not copies:
        return "Baldness markers were not called on this file.", None, "", score
    male = _has_y(index)
    auto = _get(copies, "rs2180439") + _get(copies, "rs1160312")
    x_raw = _get(copies, "rs2497938") + _get(copies, "rs1385699")
    x_score = (1.0 if _get(copies, "rs2497938") >= 1 else 0.0) + (
        1.0 if _get(copies, "rs1385699") >= 1 else 0.0
    )
    total = auto + (x_score if male else x_raw)
    note = "" if male else "Y was not called; X SNPs are scored as diploid. This row is more relevant on a male kit."
    if total >= 3:
        finding = "Higher male-pattern baldness tendency at these SNPs."
        label = "higher"
    elif total >= 1:
        finding = "Intermediate male-pattern baldness tendency at these SNPs."
        label = "intermediate"
    else:
        finding = "Lower male-pattern baldness tendency at these SNPs."
        label = "lower"
    return finding, label, note, score


def vitamin_d_finding(index, markers, *, lookups=None):
    copies, evidence = rsid_copies(index, markers, lookups=lookups)
    score = _pack(len(copies), len(markers), evidence, missing=not copies)
    if not copies:
        return "Vitamin D markers were not called on this file.", None, "", score
    low = sum(copies.values())
    if low >= 3:
        return "More alleles linked to lower vitamin D status.", "lower", "", score
    if low >= 1:
        return "Intermediate vitamin D–related genotype.", "intermediate", "", score
    return "Fewer alleles linked to lower vitamin D status.", "higher", "", score


def celiac_finding(index, markers, *, lookups=None):
    copies, evidence = rsid_copies(index, markers, lookups=lookups)
    score = _pack(len(copies), len(markers), evidence, missing=not copies)
    if not copies:
        return "Celiac-tag markers were not called on this file.", None, "", score
    dq25 = _get(copies, "rs2187668")
    dq8 = _get(copies, "rs7454108")
    extra = _get(copies, "rs7775228")
    risk = dq25 + dq8 + extra
    if dq25 >= 1.5 or (dq25 >= 0.5 and extra >= 0.5):
        return "HLA-DQ2.5 celiac-susceptibility tag present (research only, not diagnostic).", "DQ2.5", "", score
    if dq8 >= 1:
        return "HLA-DQ8 celiac-susceptibility tag present (research only, not diagnostic).", "DQ8", "", score
    if risk >= 1:
        return "One celiac-susceptibility tag allele present.", "tag", "", score
    return "No common celiac-susceptibility tags at these SNPs.", "none", "", score


def ibd_finding(index, markers, *, lookups=None):
    copies, evidence = rsid_copies(index, markers, lookups=lookups)
    score = _pack(len(copies), len(markers), evidence, missing=not copies)
    if not copies:
        return "IBD-tag markers were not called on this file.", None, "", score
    nod2 = _get(copies, "rs2066844")
    il23 = _get(copies, "rs11209026")
    if nod2 >= 1.5:
        return "NOD2 Crohn-risk allele homozygous — higher IBD-risk tag.", "NOD2", "", score
    if nod2 >= 1 and il23 <= 0:
        return "NOD2 risk allele without IL23R protective allele — Crohn-risk tag.", "risk", "", score
    if il23 >= 1:
        return "IL23R protective allele present — lower Crohn-risk tag.", "protective", "", score
    if nod2 >= 1:
        return "One NOD2 risk allele present.", "NOD2-het", "", score
    return "No common NOD2/IL23R IBD tags at these SNPs.", "none", "", score


def autoimmune_ptpn22_finding(index, markers, *, lookups=None):
    copies, evidence = rsid_copies(index, markers, lookups=lookups)
    score = _pack(len(copies), len(markers), evidence, missing=not copies)
    if not copies:
        return "PTPN22 markers were not called on this file.", None, "", score
    risk = _get(copies, "rs2476601")
    if risk >= 2:
        return "PTPN22 R620W homozygous — autoimmune-risk tag (research only).", "hom", "", score
    if risk >= 1:
        return "PTPN22 R620W carrier — autoimmune-risk tag (research only).", "het", "", score
    return "No PTPN22 R620W risk allele at rs2476601.", "none", "", score


def cyp2c19_finding(index, markers, *, lookups=None):
    copies, evidence = rsid_copies(index, markers, lookups=lookups)
    score = _pack(len(copies), len(markers), evidence, missing=not copies)
    if not copies:
        return "CYP2C19 markers were not called on this file.", None, "", score
    star2, star3, star4 = _get(copies, "rs4244285"), _get(copies, "rs4986893"), _get(copies, "rs28399504")
    loss = star2 + star3 + star4
    if star2 >= 2 or star3 >= 2 or (star2 >= 1 and star3 >= 1):
        return "Likely CYP2C19 poor metabolizer — clopidogrel/PPI effect tags (research only).", "poor", "", score
    if loss >= 1:
        return "One CYP2C19 reduced-function allele — intermediate metabolizer lean.", "intermediate", "", score
    return "Typical CYP2C19 metabolizer at these tag SNPs.", "typical", "", score


def tpmt_finding(index, markers, *, lookups=None):
    copies, evidence = rsid_copies(index, markers, lookups=lookups)
    score = _pack(len(copies), len(markers), evidence, missing=not copies)
    if not copies:
        return "TPMT markers were not called on this file.", None, "", score
    risk = _get(copies, "rs1142345")
    if risk >= 2:
        return "TPMT reduced-function homozygous — thiopurine toxicity risk tag (research only).", "high", "", score
    if risk >= 1:
        return "TPMT reduced-function carrier — thiopurine toxicity risk tag.", "carrier", "", score
    return "No common TPMT *3A allele at rs1142345.", "typical", "", score


def slco1b1_finding(index, markers, *, lookups=None):
    copies, evidence = rsid_copies(index, markers, lookups=lookups)
    score = _pack(len(copies), len(markers), evidence, missing=not copies)
    if not copies:
        return "SLCO1B1 markers were not called on this file.", None, "", score
    risk = _get(copies, "rs4149056")
    if risk >= 2:
        return "SLCO1B1 reduced-function homozygous — statin myopathy risk tag.", "high", "", score
    if risk >= 1:
        return "SLCO1B1 reduced-function carrier — statin myopathy risk tag.", "carrier", "", score
    return "Typical SLCO1B1 genotype at rs4149056.", "typical", "", score


def nat2_finding(index, markers, *, lookups=None):
    copies, evidence = rsid_copies(index, markers, lookups=lookups)
    score = _pack(len(copies), len(markers), evidence, missing=not copies)
    if not copies:
        return "NAT2 markers were not called on this file.", None, "", score
    slow = sum(_get(copies, r) for r in ("rs1800566", "rs1799930", "rs1208") if r in copies)
    if slow >= 3:
        return "Slow NAT2 acetylator tag — slower clearance of some drugs (research only).", "slow", "", score
    if slow >= 1:
        return "Intermediate NAT2 acetylator tag.", "intermediate", "", score
    return "Fast NAT2 acetylator lean at these tag SNPs.", "fast", "", score


def hair_color_finding(index, markers, *, lookups=None):
    copies, evidence = rsid_copies(index, markers, lookups=lookups)
    score = _pack(len(copies), len(markers), evidence, missing=not copies)
    if not copies:
        return "Hair-color markers were not called on this file.", None, "", score
    light = sum(copies.values())
    if light >= 2:
        return "Lighter hair-color alleles (blond/light brown lean).", "light", "", score
    if light >= 1:
        return "Mixed hair-color alleles.", "mixed", "", score
    return "Darker hair-color alleles lean.", "dark", "", score


def nicotine_finding(index, markers, *, lookups=None):
    copies, evidence = rsid_copies(index, markers, lookups=lookups)
    score = _pack(len(copies), len(markers), evidence, missing=not copies)
    if not copies:
        return "Nicotine-dependence markers were not called on this file.", None, "", score
    risk = _get(copies, "rs16969968") + _get(copies, "rs1051730")
    if risk >= 2:
        return "Higher nicotine-dependence risk alleles at CHRNA5/CHRNA3.", "higher", "", score
    if risk >= 1:
        return "One nicotine-dependence risk allele.", "intermediate", "", score
    return "Lower nicotine-dependence risk alleles at these SNPs.", "lower", "", score


def obesity_finding(index, markers, *, lookups=None):
    copies, evidence = rsid_copies(index, markers, lookups=lookups)
    score = _pack(len(copies), len(markers), evidence, missing=not copies)
    if not copies:
        return "FTO obesity markers were not called on this file.", None, "", score
    risk = _get(copies, "rs9939609") + _get(copies, "rs1421085")
    if risk >= 3:
        return "More FTO obesity-risk alleles.", "higher", "", score
    if risk >= 1:
        return "One or two FTO obesity-risk alleles.", "intermediate", "", score
    return "Fewer FTO obesity-risk alleles.", "lower", "", score


def t2d_single_finding(index, markers, *, lookups=None):
    copies, evidence = rsid_copies(index, markers, lookups=lookups)
    score = _pack(len(copies), len(markers), evidence, missing=not copies)
    if not copies:
        return "T2D single-SNP markers were not called on this file.", None, "", score
    risk = _get(copies, "rs7903146") + _get(copies, "rs5219")
    if risk >= 2:
        return "Multiple T2D susceptibility-tag alleles (TCF7L2/KCNJ11).", "higher", "", score
    if risk >= 1:
        return "One T2D susceptibility-tag allele.", "tag", "", score
    return "No T2D risk alleles at these tag SNPs.", "none", "", score


def _apoe_allele(c429: float, c7412: float) -> str:
    e4 = c429 >= 1
    e2 = c7412 >= 1
    if e4 and e2:
        return "e2/e4"
    if c429 >= 2:
        return "e4/e4"
    if e4:
        return "e3/e4"
    if c7412 >= 2:
        return "e2/e2"
    if e2:
        return "e2/e3"
    return "e3/e3"


def apoe_finding(index, markers, *, lookups=None):
    copies, evidence = rsid_copies(index, markers, lookups=lookups)
    score = _pack(len(copies), len(markers), evidence, missing=not copies)
    if not copies:
        return "APOE markers were not called on this file.", None, "", score
    e4 = _get(copies, "rs429358")
    e2 = _get(copies, "rs7412")
    if "rs429358" not in copies or "rs7412" not in copies:
        if e4 >= 1:
            return "APOE e4 allele present (partial call — research only).", "e4", "Only one APOE tag SNP was called.", score
        if e2 >= 1:
            return "APOE e2 allele present (partial call — research only).", "e2", "Only one APOE tag SNP was called.", score
        return "APOE tag SNPs incomplete on this file.", None, "", score
    hap = _apoe_allele(e4, e2)
    note = "Research tag only — not diagnostic for Alzheimer's disease."
    if "e4/e4" in hap or hap == "e3/e4":
        return f"Likely APOE {hap} — e4 Alzheimer-risk tag present.", hap, note, score
    if "e2" in hap:
        return f"Likely APOE {hap}.", hap, note, score
    return f"Likely APOE {hap} — no e4 allele at these tag SNPs.", hap, note, score


def cilantro_finding(index, markers, *, lookups=None):
    copies, evidence = rsid_copies(index, markers, lookups=lookups)
    score = _pack(len(copies), len(markers), evidence, missing=not copies)
    if not copies:
        return "Cilantro-taste markers were not called on this file.", None, "", score
    soapy = _get(copies, "rs8021378")
    if "rs2736100" in copies:
        soapy += _get(copies, "rs2736100")
    if soapy >= 2:
        return "Likely to taste cilantro as soapy.", "soapy", "", score
    if soapy >= 1:
        return "Intermediate cilantro/soapy taste.", "mixed", "", score
    return "Unlikely to taste cilantro as soapy.", "typical", "", score


def photic_sneeze_finding(index, markers, *, lookups=None):
    copies, evidence = rsid_copies(index, markers, lookups=lookups)
    score = _pack(len(copies), len(markers), evidence, missing=not copies)
    if not copies:
        return "Photic-sneeze markers were not called on this file.", None, "", score
    risk = _get(copies, "rs6151429")
    if risk >= 2:
        return "Likely photic sneeze reflex (ACHOO).", "likely", "", score
    if risk >= 1:
        return "Intermediate photic sneeze reflex.", "intermediate", "", score
    return "Unlikely photic sneeze reflex.", "unlikely", "", score


def motion_sickness_finding(index, markers, *, lookups=None):
    copies, evidence = rsid_copies(index, markers, lookups=lookups)
    score = _pack(len(copies), len(markers), evidence, missing=not copies)
    if not copies:
        return "Motion-sickness markers were not called on this file.", None, "", score
    risk = _get(copies, "rs10512472")
    if risk >= 2:
        return "Higher motion-sickness susceptibility alleles.", "higher", "", score
    if risk >= 1:
        return "One motion-sickness susceptibility allele.", "intermediate", "", score
    return "Lower motion-sickness susceptibility alleles.", "lower", "", score


COMBO_HANDLERS = {
    "hfe": hfe_finding,
    "clotting": clotting_finding,
    "mthfr": mthfr_finding,
    "warfarin": warfarin_finding,
    "alcohol": alcohol_finding,
    "pigmentation": pigmentation_finding,
    "eye_color": eye_color_finding,
    "mc1r": mc1r_finding,
    "bitter_taste": bitter_diplotype_finding,
    "actn3": actn3_ace_finding,
    "baldness": baldness_finding,
    "vitamin_d": vitamin_d_finding,
    "celiac": celiac_finding,
    "ibd": ibd_finding,
    "autoimmune_ptpn22": autoimmune_ptpn22_finding,
    "cyp2c19": cyp2c19_finding,
    "tpmt": tpmt_finding,
    "slco1b1": slco1b1_finding,
    "nat2": nat2_finding,
    "hair_color": hair_color_finding,
    "nicotine": nicotine_finding,
    "obesity": obesity_finding,
    "t2d_single": t2d_single_finding,
    "apoe": apoe_finding,
    "cilantro": cilantro_finding,
    "photic_sneeze": photic_sneeze_finding,
    "motion_sickness": motion_sickness_finding,
}
