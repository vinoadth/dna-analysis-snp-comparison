import gzip
import tempfile
import unittest
from pathlib import Path

from dna_compare.comparisons.disease import (
    _coverage_status,
    _percentile_from_reference,
    _reference_from_scoring_file,
    _relative_level,
    overlap_pgs_variants_for_query,
    parse_pgs_variants,
    score_pgs,
)
from dna_compare.genotype_index import build_rsid_index
from dna_compare.config import DISEASE_DIR
from dna_compare.models import AnalysisResult
from dna_compare.service import AnalysisService
from dna_compare.vcf_parser import parse_vcf

FIXTURE = Path(__file__).resolve().parents[1] / "data" / "samples" / "demo.snps.vcf"


class DiseasePgsTests(unittest.TestCase):
    def test_parse_bundled_cancer_scores(self):
        breast = parse_pgs_variants(DISEASE_DIR / "PGS000004.txt.gz")
        prostate = parse_pgs_variants(DISEASE_DIR / "PGS000030.txt.gz")
        crc = parse_pgs_variants(DISEASE_DIR / "PGS000765.txt.gz")
        cad = parse_pgs_variants(DISEASE_DIR / "PGS000011.txt.gz")
        t2d = parse_pgs_variants(DISEASE_DIR / "PGS004226.txt.gz")
        depression = parse_pgs_variants(DISEASE_DIR / "PGS000145.txt.gz")
        self.assertGreaterEqual(len(breast), 250)
        self.assertGreaterEqual(len(prostate), 100)
        self.assertGreaterEqual(len(crc), 80)
        self.assertGreaterEqual(len(cad), 40)
        self.assertGreaterEqual(len(t2d), 40)
        self.assertGreaterEqual(len(depression), 20000)
        self.assertTrue(all(len(row["effect"]) == 1 for row in breast))

    def test_scores_known_demo_snp(self):
        _summary, index = parse_vcf(FIXTURE)
        variants = [
            {
                "rsid": "rs3094315",
                "chrom": "1",
                "pos": 752566,
                "effect": "A",
                "other": "G",
                "weight": 0.5,
            }
        ]
        raw, used, n_score = score_pgs(index, variants)
        self.assertEqual(n_score, 1)
        self.assertEqual(used, 1)
        self.assertAlmostEqual(raw, 0.5)

    def test_flipped_effect_allele_uses_ref_copies(self):
        _summary, index = parse_vcf(FIXTURE)
        variants = [
            {
                "rsid": "rs3094315",
                "chrom": "1",
                "pos": 752566,
                "effect": "G",
                "other": "A",
                "weight": 1.0,
            }
        ]
        raw, used, _n = score_pgs(index, variants)
        self.assertEqual(used, 1)
        self.assertAlmostEqual(raw, 1.0)

    def test_analyze_includes_disease_block(self):
        result = AnalysisService().analyze(
            FIXTURE,
            filename="demo.snps.vcf",
            compare_hominin_flag=False,
            compare_caste_flag=False,
            compare_populations_flag=False,
            compare_ancestry_flag=False,
            compare_haplogroups_flag=False,
        )
        self.assertIsInstance(result, AnalysisResult)
        payload = result.to_dict()
        self.assertIn("disease", payload)
        self.assertEqual(payload["disease"]["kind"], "disease")
        traits = {row["trait"] for row in payload["disease"]["estimates"]}
        self.assertIn("Schizophrenia", traits)
        self.assertIn("Major depressive disorder", traits)
        self.assertIn("Bipolar disorder", traits)
        self.assertIn("ADHD", traits)
        self.assertIn("Anxiety disorders", traits)
        self.assertIn("Fluid intelligence", traits)
        self.assertIn("Educational attainment", traits)
        self.assertIn("Breast cancer", traits)
        self.assertIn("Prostate cancer", traits)
        self.assertIn("Colorectal cancer", traits)
        self.assertIn("Coronary artery disease", traits)
        self.assertIn("Type 2 diabetes", traits)
        self.assertIn("Stroke", traits)
        self.assertIn("Alzheimer's disease", traits)
        self.assertIn("Asthma", traits)

    def test_overlap_scoring_matches_full_scan(self):
        _summary, index = parse_vcf(FIXTURE)
        by_rsid = build_rsid_index(index)
        path = DISEASE_DIR / "PGS000011.txt.gz"
        overlap, n_total = overlap_pgs_variants_for_query(
            path,
            by_rsid=by_rsid,
            query_index=index,
        )
        self.assertGreater(n_total, 0)
        self.assertLessEqual(len(overlap), n_total)
        raw_overlap, used_overlap, _ = score_pgs(index, overlap, by_rsid=by_rsid, n_score=n_total)
        raw_full, used_full, _ = score_pgs(
            index,
            parse_pgs_variants(path),
            by_rsid=by_rsid,
            n_score=n_total,
        )
        self.assertEqual(used_overlap, used_full)
        self.assertAlmostEqual(raw_overlap, raw_full)

    def test_coverage_status_and_elevated_percentile(self):
        self.assertEqual(_coverage_status(used=0, score=None, coverage_pct=0.0), "missing")
        self.assertEqual(_coverage_status(used=5, score=1.0, coverage_pct=5.0), "low")
        self.assertEqual(_coverage_status(used=5, score=1.0, coverage_pct=15.0), "partial")
        self.assertEqual(_coverage_status(used=40, score=3.4, coverage_pct=82.0), "ready")
        refs = {"PGS000011": {"mean": 3.82, "sd": 0.43, "ref_pop": "European"}}
        low_pct = _percentile_from_reference(3.402, coverage_pct=82.0, ref=refs["PGS000011"])
        self.assertIsNotNone(low_pct)
        self.assertLess(low_pct, 50.0)
        high_pct = _percentile_from_reference(4.5, coverage_pct=82.0, ref=refs["PGS000011"])
        self.assertIsNotNone(high_pct)
        self.assertGreaterEqual(high_pct, 90.0)
        self.assertIsNone(_percentile_from_reference(4.5, coverage_pct=30.0, ref=refs["PGS000011"]))
        self.assertEqual(
            _relative_level(trait_category="disease", coverage_pct=82.0, percentile=high_pct),
            "elevated",
        )
        self.assertEqual(
            _relative_level(trait_category="disease", coverage_pct=82.0, percentile=low_pct),
            "typical",
        )
        self.assertEqual(
            _relative_level(trait_category="cognitive", coverage_pct=82.0, percentile=95.0),
            "",
        )
        cad_ref = _reference_from_scoring_file(str(DISEASE_DIR / "PGS000011.txt.gz"))
        self.assertIsNotNone(cad_ref)
        mean, sd, n_af, n_total = cad_ref
        self.assertGreaterEqual(n_af, 40)
        self.assertAlmostEqual(mean, 3.82, delta=0.05)
        self.assertAlmostEqual(sd, 0.43, delta=0.08)

    def test_analyze_includes_disease_presentation_fields(self):
        result = AnalysisService().analyze(
            FIXTURE,
            filename="demo.snps.vcf",
            compare_hominin_flag=False,
            compare_caste_flag=False,
            compare_populations_flag=False,
            compare_ancestry_flag=False,
            compare_haplogroups_flag=False,
        )
        row = payload = result.to_dict()["disease"]["estimates"][0]
        self.assertIn("coverage_status", row)
        self.assertIn("relative_level", row)
        self.assertIn("trait_category", row)

    def test_parse_custom_scoring_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "tiny.txt.gz"
            with gzip.open(path, "wt") as handle:
                handle.write("#pgs_id=TEST\n")
                handle.write("rsID\tchr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\n")
                handle.write("rs3094315\t1\t752566\tA\tG\t0.25\n")
                handle.write("indel\t1\t1\tAT\tA\t9.0\n")
            rows = parse_pgs_variants(path)
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["rsid"], "rs3094315")


if __name__ == "__main__":
    unittest.main()
