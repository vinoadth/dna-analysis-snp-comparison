import unittest
from io import StringIO

from dna_compare.comparisons.haplogroups import (
    _MT_SPECIFICITY,
    best_sample_haplogroup,
    collapse_mt_haplogroup,
    collapse_y_haplogroup,
    compare_haplogroups,
    marker_quality_note,
    score_mt_markers,
    score_y_markers,
)
from dna_compare.config import default_settings
from dna_compare.vcf_parser import parse_vcf

M17_CONFLICT_VCF = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TEST1
Y	15026424	rs2032624	A	C	.	PASS	.	GT	0
Y	21733165	rs3908	D	I	.	PASS	.	GT	1
Y	14969634	rs2032604	T	G	.	PASS	.	GT	0
"""

QUALITY_VCF = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TEST1
Y	15026424	rs2032624	A	C	.	PASS	.	GT:GQ:IGC	0:5:0.65
Y	21733165	rs3908	D	I	12	PASS	.	GT:GQ:DP	1:40:18
Y	15581983	rs2032658	A	G	.	PASS	.	GT:GQ:IGC	1:9:0.86
"""

R1A1_VCF = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TEST1
Y	15581983	rs2032658	A	G	.	PASS	.	GT	1
Y	15026424	rs2032624	A	C	.	PASS	.	GT	1
Y	21733165	rs3908	D	I	.	PASS	.	GT	1
"""


class HaplogroupTests(unittest.TestCase):
    def test_collapse_r1a1_and_vellalar_style_labels(self):
        self.assertEqual(collapse_y_haplogroup("R1a1a1"), "R1a1 (M17)")
        self.assertEqual(collapse_y_haplogroup("R1a1a1b2a2~"), "R1a-Z93")
        self.assertEqual(collapse_y_haplogroup("L1a1b"), "L1 (M27)")
        self.assertEqual(collapse_y_haplogroup("J2b2a2b1a~"), "J2b (M241)")
        self.assertEqual(collapse_y_haplogroup("R2a"), "R2 (M124)")
        self.assertEqual(collapse_y_haplogroup("R1a-Z93"), "R1a-Z93")
        self.assertEqual(collapse_y_haplogroup("R1a1a1b2a"), "R1a-Z93")
        self.assertEqual(collapse_y_haplogroup("LT"), "other")
        self.assertEqual(collapse_y_haplogroup("T1a"), "T (M70)")
        self.assertEqual(collapse_y_haplogroup("E1b1b"), "E (M96)")
        self.assertEqual(collapse_y_haplogroup("H3b1a"), "H3b (Z13871)")
        self.assertEqual(collapse_y_haplogroup("H-Z13871"), "H3b (Z13871)")
        self.assertEqual(collapse_y_haplogroup("H3-Z5857"), "H3 (Z5857)")
        self.assertEqual(collapse_y_haplogroup("H1a1a"), "H1 (M52)")
        self.assertEqual(collapse_y_haplogroup("H-M69"), "H (M69)")
        self.assertIsNone(collapse_y_haplogroup("n/a (female)"))
        self.assertIsNone(collapse_y_haplogroup(".."))

    def test_m17_indel_conflicts_when_m173_ancestral(self):
        _summary, index = parse_vcf(StringIO(M17_CONFLICT_VCF))
        self.assertIn(("Y", 21733165), index)
        calls = {c.marker: c for c in score_y_markers(index)}
        self.assertEqual(calls["M173"].status, "ancestral")
        self.assertEqual(calls["M17"].status, "conflict")
        self.assertIsNone(best_sample_haplogroup(list(calls.values())))

    def test_derived_m17_with_r_and_r1(self):
        _summary, index = parse_vcf(StringIO(R1A1_VCF))
        calls = score_y_markers(index)
        self.assertEqual(best_sample_haplogroup(calls), "R1a1 (M17)")

    def test_marker_calls_include_rsid_and_allele_fields(self):
        _summary, index = parse_vcf(StringIO(R1A1_VCF))
        calls = score_y_markers(index)
        m173 = next(c for c in calls if c.marker == "M173")
        self.assertEqual(m173.rsid, "rs2032624")
        self.assertEqual(m173.chrom, "Y")
        self.assertEqual(m173.pos, 15026424)
        self.assertEqual(m173.ancestral, "A")
        self.assertEqual(m173.derived, "C")
        self.assertEqual(m173.observed, "C")
        self.assertEqual(m173.status, "derived")
        payload = m173.to_dict()
        for key in ("rsid", "chrom", "pos", "ancestral", "derived", "observed", "status", "genotype"):
            self.assertIn(key, payload)

    def test_compare_includes_aadr_percent_rows_when_anno_present(self):
        settings = default_settings()
        _summary, index = parse_vcf(StringIO(M17_CONFLICT_VCF))
        result = compare_haplogroups(index, settings=settings)
        if not settings.aadr_anno.exists():
            self.skipTest("AADR .anno not present")
        labels = [row.haplogroup for row in result.rows]
        self.assertIn("R1a1 (M17)", labels)
        vellalar = next(row for row in result.rows if row.haplogroup == "R1a1 (M17)").groups["Vellalar"]
        self.assertEqual(vellalar.n, 0)
        self.assertGreater(vellalar.n_called, 0)
        tamil = next(row for row in result.rows if row.haplogroup == "R1a1 (M17)").groups["Tamil"]
        self.assertGreater(tamil.n, 0)

    def test_collapse_mt_takes_first_aadr_token(self):
        self.assertEqual(collapse_mt_haplogroup("M36d;M "), "M36d")
        self.assertEqual(collapse_mt_haplogroup("R8a2"), "R8a2")
        self.assertIsNone(collapse_mt_haplogroup("n/a (<2x)"))
        self.assertIsNone(collapse_mt_haplogroup(".."))

    def test_mt_r_derived_when_m_ancestral(self):
        vcf = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TEST1
MT	10400	rs28358278	C	T	.	PASS	.	GT	0
MT	12705	R	C	T	.	PASS	.	GT	0
MT	12308	U	A	G	.	PASS	.	GT	0
MT	14766	rs3135031	C	T	.	PASS	.	GT	1
"""
        _summary, index = parse_vcf(StringIO(vcf))
        calls = {c.haplogroup: c for c in score_mt_markers(index)}
        self.assertEqual(calls["M (10400T)"].status, "ancestral")
        self.assertEqual(calls["R (12705C)"].status, "derived")
        self.assertEqual(calls["U (12308G)"].status, "ancestral")
        self.assertEqual(calls["HV (14766C)"].status, "ancestral")
        self.assertEqual(best_sample_haplogroup(list(calls.values()), (
            "R8 (13215C)", "U (12308G)", "M (10400T)", "R (12705C)",
        )), "R (12705C)")

    def test_mt_m4_derived_when_m_agrees(self):
        vcf = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TEST1
MT	10400	rs28358278	C	T	.	PASS	.	GT	1
MT	12705	R	C	T	.	PASS	.	GT	1
MT	6620	M4	T	C	.	PASS	.	GT	1
MT	7271	M36	A	G	.	PASS	.	GT	0
"""
        _summary, index = parse_vcf(StringIO(vcf))
        calls = {c.haplogroup: c for c in score_mt_markers(index)}
        self.assertEqual(calls["M (10400T)"].status, "derived")
        self.assertEqual(calls["R (12705C)"].status, "ancestral")
        self.assertEqual(calls["M4 (6620C)"].status, "derived")
        self.assertEqual(calls["M36 (7271G)"].status, "ancestral")
        self.assertEqual(best_sample_haplogroup(list(calls.values()), _MT_SPECIFICITY), "M4 (6620C)")

    def test_mt_m4_conflicts_when_m_ancestral(self):
        vcf = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TEST1
MT	10400	rs28358278	C	T	.	PASS	.	GT	0
MT	6620	M4	T	C	.	PASS	.	GT	1
"""
        _summary, index = parse_vcf(StringIO(vcf))
        calls = {c.haplogroup: c for c in score_mt_markers(index)}
        self.assertEqual(calls["M (10400T)"].status, "ancestral")
        self.assertEqual(calls["M4 (6620C)"].status, "conflict")
        self.assertIsNone(best_sample_haplogroup(list(calls.values()), _MT_SPECIFICITY))

    def test_mt_h3b_and_j_backbone_calls(self):
        vcf = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TEST1
MT	12705	R	C	T	.	PASS	.	GT	0
MT	14766	rs3135031	C	T	.	PASS	.	GT	0
MT	7028	H	C	T	.	PASS	.	GT	0
MT	6776	H3	T	C	.	PASS	.	GT	1
MT	2581	H3b	A	G	.	PASS	.	GT	1
MT	16126	rs147029798	T	C	.	PASS	.	GT	0
MT	16069	J	C	T	.	PASS	.	GT	0
MT	14905	T	G	A	.	PASS	.	GT	0
"""
        _summary, index = parse_vcf(StringIO(vcf))
        calls = {c.haplogroup: c for c in score_mt_markers(index)}
        self.assertEqual(calls["H (7028C)"].status, "derived")
        self.assertEqual(calls["H3 (6776C)"].status, "derived")
        self.assertEqual(calls["H3b (2581G)"].status, "derived")
        self.assertEqual(calls["JT (16126C)"].status, "ancestral")
        self.assertEqual(calls["J (16069T)"].status, "ancestral")
        self.assertEqual(best_sample_haplogroup(list(calls.values()), _MT_SPECIFICITY), "H3b (2581G)")

    def test_mt_n1_conflicts_when_m_derived(self):
        vcf = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TEST1
MT	10400	rs28358278	C	T	.	PASS	.	GT	1
MT	12705	R	C	T	.	PASS	.	GT	1
MT	10238	N1	T	C	.	PASS	.	GT	1
"""
        _summary, index = parse_vcf(StringIO(vcf))
        calls = {c.haplogroup: c for c in score_mt_markers(index)}
        self.assertEqual(calls["M (10400T)"].status, "derived")
        self.assertEqual(calls["N1 (10238C)"].status, "conflict")
        self.assertEqual(best_sample_haplogroup(list(calls.values()), _MT_SPECIFICITY), "M (10400T)")

    def test_omit_y_table_when_vcf_has_no_chrY(self):
        vcf = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	MOM1
MT	10400	rs28358278	C	T	.	PASS	.	GT	0
MT	12705	R	C	T	.	PASS	.	GT	0
1	752566	rs3094315	G	A	.	PASS	.	GT	0/1
"""
        _summary, index = parse_vcf(StringIO(vcf))
        result = compare_haplogroups(index, settings=default_settings())
        self.assertFalse(result.available)
        self.assertEqual(result.rows, [])
        self.assertEqual(result.markers, [])
        self.assertTrue(result.mt_available)
        self.assertTrue(any("Y haplogroup table is omitted" in note for note in result.notes))

    def test_omit_mt_table_when_vcf_has_no_chrMT(self):
        _summary, index = parse_vcf(StringIO(M17_CONFLICT_VCF))
        result = compare_haplogroups(index, settings=default_settings())
        self.assertTrue(result.available)
        self.assertFalse(result.mt_available)
        self.assertEqual(result.mt_rows, [])
        self.assertEqual(result.mt_markers, [])
        self.assertTrue(any("mtDNA haplogroup table is omitted" in note for note in result.mt_notes))

    def test_compare_mt_counts_only_groups_with_anno_calls(self):
        settings = default_settings()
        mt_vcf = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TEST1
MT	10400	rs28358278	C	T	.	PASS	.	GT	0
MT	12705	R	C	T	.	PASS	.	GT	0
"""
        _summary, index = parse_vcf(StringIO(mt_vcf))
        result = compare_haplogroups(index, settings=settings)
        if not settings.aadr_anno.exists():
            self.skipTest("AADR .anno not present")
        labels = [row.haplogroup for row in result.mt_rows]
        self.assertIn("R (12705C)", labels)
        self.assertIn("M (10400T)", labels)
        self.assertIn("M4 (6620C)", labels)
        self.assertIn("M36 (7271G)", labels)
        self.assertIn("H3 (6776C)", labels)
        self.assertIn("U7 (5360T)", labels)
        r_row = next(row for row in result.mt_rows if row.haplogroup == "R (12705C)")
        m_row = next(row for row in result.mt_rows if row.haplogroup == "M (10400T)")
        called = [name for name, cell in r_row.groups.items() if cell.n_called]
        self.assertTrue(called)
        self.assertTrue(any("dash means" in note.lower() or "no mt haplogroup" in note.lower() for note in result.mt_notes))
        joined = " ".join(result.status_notes).lower()
        self.assertIn("derived = yes", joined)
        self.assertIn("ancestral = no", joined)
        self.assertIn("no-call", joined)
        self.assertIn("conflict", joined)
        self.assertFalse(any("derived = yes" in (n or "").lower() for n in result.notes + result.mt_notes))
        # Vellalar/VLR is one HO label that happens to have mt filled in.
        if "Vellalar" in r_row.groups and r_row.groups["Vellalar"].n_called:
            self.assertEqual(r_row.groups["Vellalar"].n, 5)
            self.assertEqual(m_row.groups["Vellalar"].n, 4)

    def test_hg38_positions_call_r1a1_without_chain(self):
        vcf = """##fileformat=VCFv4.2
##source=bcftools_gtc2vcf
##reference=GSA-24v3-0_A1
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TEST1
Y	13470103	rs2032658	A	G	.	PASS	.	GT	1
Y	12914512	rs2032624	A	C	.	PASS	.	GT	1
Y	19571279	rs3908	D	I	.	PASS	.	GT	1
"""
        _summary, index = parse_vcf(StringIO(vcf))
        calls = score_y_markers(index)
        self.assertEqual(best_sample_haplogroup(calls), "R1a1 (M17)")

    def test_rsid_wins_over_decoy_at_hg19_coordinate(self):
        vcf = """##fileformat=VCFv4.2
##source=bcftools_gtc2vcf
##reference=GSA-24v3-0_A1
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TEST1
Y	15026424	not_m173	T	A	.	PASS	.	GT	0
Y	13470103	exm-rs2032658	A	G	.	PASS	.	GT	1
Y	12914512	rs2032624	A	C	.	PASS	.	GT	1
Y	19571279	rs3908	D	I	.	PASS	.	GT	1
"""
        _summary, index = parse_vcf(StringIO(vcf))
        calls = {c.marker: c for c in score_y_markers(index)}
        self.assertEqual(calls["M173"].status, "derived")
        self.assertEqual(calls["M17"].status, "derived")
        self.assertEqual(best_sample_haplogroup(list(calls.values())), "R1a1 (M17)")

    def test_h3b_derived_when_h_and_h3_agree(self):
        vcf = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TEST1
Y	21894058	rs2032673	T	C	.	PASS	.	GT	1
Y	21753199	rs376769460	A	C	.	PASS	.	GT	0
Y	2759285	rs569006329	C	G	.	PASS	.	GT	1
Y	2878605	Z13871	T	G	.	PASS	.	GT	1
"""
        _summary, index = parse_vcf(StringIO(vcf))
        calls = {c.marker: c for c in score_y_markers(index)}
        self.assertEqual(calls["M69"].status, "derived")
        self.assertEqual(calls["M52"].status, "ancestral")
        self.assertEqual(calls["Z5857"].status, "derived")
        self.assertEqual(calls["Z13871"].status, "derived")
        self.assertEqual(best_sample_haplogroup(list(calls.values())), "H3b (Z13871)")

    def test_h3b_conflicts_when_m69_ancestral(self):
        vcf = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TEST1
Y	21894058	rs2032673	T	C	.	PASS	.	GT	0
Y	2878605	Z13871	T	G	.	PASS	.	GT	1
"""
        _summary, index = parse_vcf(StringIO(vcf))
        calls = {c.marker: c for c in score_y_markers(index)}
        self.assertEqual(calls["M69"].status, "ancestral")
        self.assertEqual(calls["Z13871"].status, "conflict")
        self.assertIsNone(best_sample_haplogroup(list(calls.values())))

    def test_r2_and_j2b_backbone_calls(self):
        vcf = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TEST1
Y	15581983	rs2032658	A	G	.	PASS	.	GT	1
Y	21764501	M124	G	A	.	PASS	.	GT	1
Y	14969634	rs2032604	T	G	.	PASS	.	GT	1
Y	15018459	M241	G	A	.	PASS	.	GT	1
"""
        _summary, index = parse_vcf(StringIO(vcf))
        calls = {c.marker: c for c in score_y_markers(index)}
        self.assertEqual(calls["M124"].status, "derived")
        self.assertEqual(calls["M241"].status, "derived")
        self.assertEqual(best_sample_haplogroup(list(calls.values())), "R2 (M124)")

    def test_marker_quality_from_gq_dp_igc(self):
        _summary, index = parse_vcf(StringIO(QUALITY_VCF))
        self.assertEqual(index[("Y", 15026424)]["gq"], 5)
        self.assertAlmostEqual(index[("Y", 15026424)]["igc"], 0.65)
        self.assertIsNone(index[("Y", 15026424)]["dp"])
        self.assertEqual(index[("Y", 21733165)]["dp"], 18)
        self.assertEqual(index[("Y", 21733165)]["gq"], 40)
        self.assertEqual(index[("Y", 21733165)]["qual"], 12.0)
        calls = {c.marker: c for c in score_y_markers(index)}
        self.assertEqual(calls["M173"].gq, 5)
        self.assertAlmostEqual(calls["M173"].igc, 0.65)
        self.assertIsNone(calls["M173"].dp)
        self.assertEqual(calls["M17"].dp, 18)
        self.assertEqual(calls["M17"].gq, 40)
        self.assertEqual(calls["M17"].qual, 12.0)
        note = marker_quality_note(list(calls.values()), label="Y")
        self.assertIn("median GQ", note)
        self.assertIn("median DP", note)
        self.assertIn("median IGC", note)
        result = compare_haplogroups(index, settings=default_settings())
        self.assertTrue(any("call support" in n for n in result.notes))


if __name__ == "__main__":
    unittest.main()
