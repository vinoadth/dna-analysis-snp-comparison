import tempfile
import unittest
from pathlib import Path

from dna_compare.comparisons.additional import compare_additional, parent_role, score_marker_group
from dna_compare.config import default_settings
from dna_compare.service import AnalysisService
from dna_compare.vcf_parser import parse_vcf

CHILD_VCF = """\
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	CHILD
11	61552680	rs174537	G	T	.	PASS	.	GT	0/1
11	61569830	rs174546	C	T	.	PASS	.	GT	0/1
11	61570783	rs174547	T	C	.	PASS	.	GT	0/1
11	61571348	rs174548	C	G	.	PASS	.	GT	0/1
11	61571478	rs174550	T	C	.	PASS	.	GT	0/1
11	61597212	rs174570	C	T	.	PASS	.	GT	0/1
"""

FATHER_VCF = """\
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	FATHER
11	61552680	rs174537	G	T	.	PASS	.	GT	0/0
11	61569830	rs174546	C	T	.	PASS	.	GT	0/0
11	61570783	rs174547	T	C	.	PASS	.	GT	0/0
11	61571348	rs174548	C	G	.	PASS	.	GT	0/0
11	61571478	rs174550	T	C	.	PASS	.	GT	0/0
11	61597212	rs174570	C	T	.	PASS	.	GT	0/0
"""

MOTHER_VCF = """\
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	MOTHER
11	61552680	rs174537	G	T	.	PASS	.	GT	1/1
11	61569830	rs174546	C	T	.	PASS	.	GT	1/1
11	61570783	rs174547	T	C	.	PASS	.	GT	1/1
11	61571348	rs174548	C	G	.	PASS	.	GT	1/1
11	61571478	rs174550	T	C	.	PASS	.	GT	1/1
11	61597212	rs174570	C	T	.	PASS	.	GT	1/1
"""

ABO_O_VCF = """\
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S
9	136132908	rs8176719	D	I	.	PASS	.	GT	0/0
9	136131322	rs8176746	G	T	.	PASS	.	GT	0/0
"""

ABO_AB_VCF = """\
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S
9	136132908	rs8176719	D	I	.	PASS	.	GT	1/1
9	136131322	rs8176746	G	T	.	PASS	.	GT	0/1
"""

ABO_A_VCF = """\
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S
9	136132908	rs8176719	D	I	.	PASS	.	GT	0/1
9	136131322	rs8176746	G	T	.	PASS	.	GT	0/0
"""

RH_NEG_VCF = """\
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S
1	25629943	rs590787	A	C	.	PASS	.	GT	1/1
"""

RH_POS_VCF = """\
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S
1	25629943	rs590787	A	C	.	PASS	.	GT	0/0
"""

RH_POS_HET_VCF = """\
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S
1	25629943	rs590787	A	C	.	PASS	.	GT	0/1
"""

BLOOD_O_POS_VCF = """\
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S
9	136132908	rs8176719	D	I	.	PASS	.	GT	0/0
9	136131322	rs8176746	G	T	.	PASS	.	GT	0/0
1	25629943	rs590787	A	C	.	PASS	.	GT	0/0
"""

RH_DEL_VCF = """\
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S
1	25600000	i4001527	I	D	.	PASS	.	GT	1/1
"""

RH_DEL_CONFLICT_VCF = """\
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S
1	25600000	i4001527	I	D	.	PASS	.	GT	1/1
1	25629943	rs590787	A	C	.	PASS	.	GT	0/0
"""

FIXTURE = Path(__file__).resolve().parents[1] / "data" / "samples" / "demo.snps.vcf"


def _index(text: str):
    with tempfile.NamedTemporaryFile("w", suffix=".vcf", delete=False) as handle:
        handle.write(text)
        path = Path(handle.name)
    try:
        _summary, index = parse_vcf(path)
    finally:
        path.unlink(missing_ok=True)
    return index


def _topic(block, key: str):
    return next(row for row in block.estimates if row.key == key)


class AdditionalDetailsTests(unittest.TestCase):
    def test_parent_role_from_filename(self):
        self.assertEqual(parent_role("C8XY_processed_father.vcf"), "father")
        self.assertEqual(parent_role("C8XY_processed_mother.vcf"), "mother")
        self.assertIsNone(parent_role("demo.snps.vcf"))

    def test_homozygous_veg_says_community(self):
        block = compare_additional(_index(FATHER_VCF), settings=default_settings())
        row = _topic(block, "vegetarian_community")
        self.assertTrue(row.available)
        self.assertEqual(row.parent, "both")
        self.assertIn("vegetarian community", row.finding.lower())
        self.assertNotIn("father", row.finding.lower())
        self.assertNotIn("mother", row.finding.lower())

    def test_heterozygous_without_parent_is_one_parent(self):
        block = compare_additional(_index(CHILD_VCF), settings=default_settings())
        row = _topic(block, "vegetarian_community")
        self.assertEqual(row.parent, "one")
        self.assertEqual(row.finding, "One parent belongs to a vegetarian community.")
        self.assertIn("mother or father VCF", row.message)

    def test_homozygous_other_is_not_vegetarian(self):
        block = compare_additional(_index(MOTHER_VCF), settings=default_settings())
        row = _topic(block, "vegetarian_community")
        self.assertEqual(row.parent, "none")
        self.assertIn("do not show", row.finding.lower())

    def test_father_vcf_names_father(self):
        block = compare_additional(
            _index(CHILD_VCF),
            settings=default_settings(),
            other_index=_index(FATHER_VCF),
            other_filename="C8XY_processed_father.vcf",
        )
        row = _topic(block, "vegetarian_community")
        self.assertEqual(row.finding, "Father belongs to a vegetarian community.")
        self.assertEqual(row.parent, "father")

    def test_mother_vcf_names_father_as_the_other_parent(self):
        block = compare_additional(
            _index(CHILD_VCF),
            settings=default_settings(),
            other_index=_index(MOTHER_VCF),
            other_filename="C8XY_processed_mother.vcf",
        )
        row = _topic(block, "vegetarian_community")
        self.assertEqual(row.finding, "Father belongs to a vegetarian community.")
        self.assertEqual(row.parent, "father")

    def test_abo_blood_type_is_listed_first(self):
        block = compare_additional(_index(ABO_O_VCF), settings=default_settings())
        self.assertEqual(block.estimates[0].key, "abo_blood_type")
        self.assertEqual(block.estimates[0].finding, "Blood type O.")
        self.assertEqual(_topic(block, "abo_blood_type").parent, "O")

    def test_abo_ab_and_a(self):
        ab = _topic(compare_additional(_index(ABO_AB_VCF), settings=default_settings()), "abo_blood_type")
        self.assertEqual(ab.finding, "Blood type AB.")
        a = _topic(compare_additional(_index(ABO_A_VCF), settings=default_settings()), "abo_blood_type")
        self.assertEqual(a.finding, "Blood type A.")

    def test_rh_factor_follows_abo(self):
        block = compare_additional(_index(ABO_O_VCF), settings=default_settings())
        self.assertEqual(block.estimates[0].key, "abo_blood_type")
        self.assertEqual(block.estimates[1].key, "rh_factor")
        self.assertNotIn("blood_type", {row.key for row in block.estimates})

    def test_combined_blood_type_when_abo_and_rh_available(self):
        block = compare_additional(_index(BLOOD_O_POS_VCF), settings=default_settings())
        combined = _topic(block, "blood_type")
        self.assertEqual(combined.finding, "Blood type O+.")
        self.assertEqual(combined.parent, "O+")
        self.assertEqual(block.estimates[2].key, "blood_type")

    def test_rh_negative_and_positive(self):
        neg = _topic(compare_additional(_index(RH_NEG_VCF), settings=default_settings()), "rh_factor")
        self.assertEqual(neg.finding, "Likely Rh negative (−).")
        self.assertEqual(neg.parent, "−")
        self.assertIn("proxy", neg.message.lower())
        pos = _topic(compare_additional(_index(RH_POS_VCF), settings=default_settings()), "rh_factor")
        self.assertEqual(pos.finding, "Likely Rh positive (+).")
        self.assertEqual(pos.parent, "+")
        het = _topic(compare_additional(_index(RH_POS_HET_VCF), settings=default_settings()), "rh_factor")
        self.assertEqual(het.finding, "Likely Rh positive (+).")
        self.assertEqual(het.parent, "+")

    def test_rh_deletion_marker_i4001527(self):
        neg = _topic(compare_additional(_index(RH_DEL_VCF), settings=default_settings()), "rh_factor")
        self.assertEqual(neg.finding, "Likely Rh negative (−).")
        self.assertEqual(neg.parent, "−")
        self.assertIn("i4001527", neg.evidence)

    def test_rh_markers_disagree(self):
        row = _topic(compare_additional(_index(RH_DEL_CONFLICT_VCF), settings=default_settings()), "rh_factor")
        self.assertIn("unclear", row.finding.lower())
        self.assertIsNone(row.parent)

    def test_coverage_status_on_rows(self):
        block = compare_additional(_index(ABO_O_VCF), settings=default_settings())
        abo = _topic(block, "abo_blood_type")
        self.assertEqual(abo.coverage_status, "ready")
        rh = _topic(block, "rh_factor")
        self.assertEqual(rh.coverage_status, "missing")
        self.assertIn("Marker coverage on this file", block.notes[-1])

    def test_analyze_includes_additional_table_row(self):
        result = AnalysisService().analyze(
            FIXTURE,
            filename="demo.snps.vcf",
            compare_hominin_flag=False,
            compare_caste_flag=False,
            compare_populations_flag=False,
            compare_ancestry_flag=False,
            compare_haplogroups_flag=False,
            compare_disease_flag=False,
        )
        payload = result.to_dict()
        self.assertEqual(payload["additional"]["kind"], "additional")
        rows = payload["additional"]["estimates"]
        topics = {row["topic"] for row in rows}
        self.assertEqual(rows[0]["topic"], "ABO blood type")
        self.assertEqual(rows[0]["finding"], "Blood type O.")
        self.assertIn("Vegetarian community", topics)
        self.assertIn("HFE iron overload", topics)
        self.assertIn("Clotting variants", topics)
        self.assertIn("MTHFR", topics)
        self.assertIn("Warfarin sensitivity", topics)
        self.assertIn("Alcohol metabolism", topics)
        self.assertIn("Milk / lactase", topics)
        self.assertIn("Pigmentation haplotype", topics)
        self.assertIn("Bitter taste (PTC)", topics)
        self.assertIn("Sprint vs endurance", topics)
        self.assertIn("Male-pattern baldness", topics)
        self.assertIn("Vitamin D related", topics)
        self.assertIn("Coronary artery disease (9p21)", topics)
        self.assertIn("Sleep preferences", topics)
        self.assertIn("Red hair / fair skin", topics)
        self.assertEqual(_topic(result.additional, "hfe").finding, "No common HFE iron-overload variants.")
        self.assertEqual(_topic(result.additional, "clotting").finding, "No Factor V Leiden or prothrombin G20210A variant.")
        self.assertIn("Typical MTHFR", _topic(result.additional, "mthfr").finding)
        self.assertIn("Typical warfarin", _topic(result.additional, "warfarin").finding)
        self.assertIn("Typical ADH1B", _topic(result.additional, "alcohol").finding)
        self.assertIn("PAV/AVI", _topic(result.additional, "bitter_taste").finding)
        self.assertIn("ACTN3 RX", _topic(result.additional, "actn3").finding)

    def test_hfe_c282y_homozygous(self):
        vcf = """\
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S
6	26093141	rs1800562	G	A	.	PASS	.	GT	1/1
6	26091179	rs1799945	C	G	.	PASS	.	GT	0/0
6	26091185	rs1800730	A	T	.	PASS	.	GT	0/0
"""
        row = _topic(compare_additional(_index(vcf), settings=default_settings()), "hfe")
        self.assertIn("C282Y homozygous", row.finding)

    def test_factor_v_leiden_carrier(self):
        vcf = """\
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S
1	169519049	rs6025	T	C	.	PASS	.	GT	0/1
11	46761055	rs1799963	G	A	.	PASS	.	GT	0/0
"""
        row = _topic(compare_additional(_index(vcf), settings=default_settings()), "clotting")
        self.assertEqual(row.finding, "Factor V Leiden carrier.")

    def test_generic_future_key_uses_tsv_findings(self):
        settings = default_settings()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "details.tsv"
            path.write_text(
                "key\ttopic\tchrom\trsid\tpos_hg19\tpos_hg38\teffect_allele\tother_allele\tgene\tsource\thom_finding\thet_finding\tref_finding\n"
                "height\tAdult height\t11\trs174547\t61570783\t61803311\tT\tC\tFADS1\tdemo\tLikely taller than the reference.\tHeight is intermediate.\tLikely shorter than the reference.\n",
                encoding="utf-8",
            )
            settings.additional_details = path
            block = compare_additional(_index(FATHER_VCF), settings=settings)
            self.assertEqual(block.estimates[0].finding, "Likely taller than the reference.")

    def test_score_uses_ref_copies_when_effect_is_ref(self):
        markers = [
            {
                "chrom": "11",
                "rsid": "rs174547",
                "pos_hg19": 61570783,
                "pos_hg38": 0,
                "effect_allele": "T",
                "other_allele": "C",
                "gene": "FADS1",
            }
        ]
        scored = score_marker_group(_index(FATHER_VCF), markers)
        self.assertEqual(scored["status"], "both")
        self.assertAlmostEqual(scored["mean"], 2.0)
        scored_child = score_marker_group(_index(CHILD_VCF), markers)
        self.assertEqual(scored_child["status"], "one")


if __name__ == "__main__":
    unittest.main()
