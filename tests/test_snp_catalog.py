import tempfile
import unittest
from pathlib import Path

from dna_compare.config import default_settings
from dna_compare.service import AnalysisService
from dna_compare.snp_catalog import catalog_path, gt_carries_alt, query_snp_catalog, write_snp_catalog
from dna_compare.vcf_parser import chrom_sort_key, parse_vcf, preview_rows

FIXTURE = Path(__file__).resolve().parents[1] / "data" / "samples" / "demo.snps.vcf"
ROOT = Path(__file__).resolve().parents[1]


class SnpCatalogTests(unittest.TestCase):
    def test_chrom_sort_puts_x_after_autosomes(self):
        ordered = sorted(["X", "2", "10", "1", "MT", "Y"], key=chrom_sort_key)
        self.assertEqual(ordered, ["1", "2", "10", "X", "Y", "MT"])

    def test_stratified_preview_covers_multiple_chromosomes(self):
        _summary, index = parse_vcf(FIXTURE)
        rows = preview_rows(index, limit=24, stratified=True)
        chroms = {row.chrom for row in rows}
        self.assertIn("1", chroms)
        self.assertTrue(chroms - {"1"})

    def test_catalog_lists_every_chromosome(self):
        _summary, index = parse_vcf(FIXTURE)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "demo.snp_catalog.tsv"
            n = write_snp_catalog(index, path)
            self.assertGreater(n, 200)
            all_rows = query_snp_catalog(path, limit=1000)
            self.assertTrue(all_rows["ok"])
            self.assertEqual(all_rows["total"], n)
            chroms = {row["chrom"] for row in all_rows["rows"]}
            self.assertIn("1", chroms)
            self.assertIn("22", chroms)
            self.assertIn("X", chroms)
            only_x = query_snp_catalog(path, chrom="X", limit=50)
            self.assertGreaterEqual(only_x["total"], 1)
            self.assertTrue(all(row["chrom"] == "X" for row in only_x["rows"]))
            page = query_snp_catalog(path, offset=20, limit=10)
            self.assertEqual(len(page["rows"]), 10)
            self.assertEqual(page["offset"], 20)
            called = query_snp_catalog(path, called_only=True, limit=1000)
            self.assertLess(called["total"], all_rows["total"])
            self.assertTrue(called["called_only"])
            self.assertTrue(all(gt_carries_alt(row["genotype"], row.get("dosage_alt")) for row in called["rows"]))

    def test_analyze_writes_catalog_for_api(self):
        settings = default_settings()
        result = AnalysisService(settings).analyze(
            FIXTURE,
            filename="demo.snps.vcf",
            compare_hominin_flag=False,
            compare_caste_flag=False,
            compare_populations_flag=False,
            compare_ancestry_flag=False,
            compare_haplogroups_flag=False,
        )
        self.assertTrue(result.vcf.snp_catalog_id)
        path = catalog_path("demo.snps.vcf", settings)
        self.assertTrue(path.is_file())
        payload = query_snp_catalog(path, chrom="22", limit=50)
        self.assertGreater(payload["total"], 0)
        self.assertTrue(all(row["chrom"] == "22" for row in payload["rows"]))

    def test_gt_carries_alt(self):
        self.assertFalse(gt_carries_alt("0/0"))
        self.assertFalse(gt_carries_alt("0|0"))
        self.assertFalse(gt_carries_alt("0"))
        self.assertTrue(gt_carries_alt("0/1"))
        self.assertTrue(gt_carries_alt("1/0"))
        self.assertTrue(gt_carries_alt("1/1"))
        self.assertTrue(gt_carries_alt("1|0"))
        self.assertTrue(gt_carries_alt("1"))
        self.assertFalse(gt_carries_alt("./."))

    def test_makefile_and_dockerfile(self):
        makefile = (ROOT / "Makefile").read_text(encoding="utf-8")
        dockerfile = (ROOT / "Dockerfile").read_text(encoding="utf-8")
        self.assertRegex(makefile, r"(?m)^build:")
        self.assertRegex(makefile, r"(?m)^serve:")
        self.assertIn("docker build", makefile)
        self.assertIn("docker run", makefile)
        self.assertIn("0.0.0.0", dockerfile)
        self.assertIn("8765", dockerfile)


if __name__ == "__main__":
    unittest.main()
