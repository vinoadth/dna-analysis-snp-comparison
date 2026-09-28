import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

from dna_compare.vcf_parser import parse_vcf

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "merge-wes-snp-vcf.py"


def load_merger():
    spec = importlib.util.spec_from_file_location("merge_wes_snp_vcf", MODULE_PATH)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


WES = """##fileformat=VCFv4.2
##reference=GRCh38
##contig=<ID=chr1,length=248956422>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	WES_S1	OTHER
chr1	100	rsW1	A	G	50	PASS	DP=30	GT:DP	0/1:30	1/1:10
chr1	300	rsW3	C	T	50	PASS	.	GT:DP	./.:0	0/1:5
chr2	50	.	G	GA	50	PASS	.	GT	0/1	0/0
"""

SNP = """##fileformat=VCFv4.2
##source=bcftools_gtc2vcf GSA-24v3
##contig=<ID=1,length=248956422>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	ARRAY_S1
1	100	rsA1	A	G	.	PASS	.	GT:GQ	1/1:99
1	200	rsA2	T	C	.	PASS	.	GT:GQ	0/1:99
1	300	rsA3	C	T	.	PASS	.	GT:GQ	1/1:99
1	400	rsA4	G	A	.	PASS	.	GT:GQ	./.:0
MT	500	rsM	A	G	.	PASS	.	GT	1
"""


class MergeWesSnpTests(unittest.TestCase):
    def _merge(self, wes_text=WES, snp_text=SNP):
        mod = load_merger()
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        wes = Path(tmp.name) / "sample.wes.vcf"
        snp = Path(tmp.name) / "sample.snps.vcf"
        wes.write_text(wes_text)
        snp.write_text(snp_text)
        return mod.merge_vcfs(wes, snp)

    def test_wes_base_with_snp_fill(self):
        out, stats = self._merge()
        self.assertEqual(out.name, "sample.wes_wes_snp_merged.vcf")
        body = [ln.split("\t") for ln in out.read_text().splitlines() if not ln.startswith("#")]
        by_pos = {(c[0], int(c[1])): c for c in body}

        self.assertEqual(by_pos[("chr1", 100)][2], "rsW1")
        self.assertEqual(by_pos[("chr1", 100)][9], "0/1:30")
        self.assertIn("SRC=WES", by_pos[("chr1", 100)][7])
        self.assertEqual(by_pos[("chr1", 200)][2], "rsA2")
        self.assertIn("SRC=SNP", by_pos[("chr1", 200)][7])
        self.assertEqual(by_pos[("chr1", 300)][2], "rsA3")
        self.assertNotIn(("chr1", 400), by_pos)
        self.assertIn(("chrM", 500), by_pos)
        self.assertIn(("chr2", 50), by_pos)
        self.assertTrue(all(len(c) == 10 for c in body))
        self.assertEqual([c[0] for c in body], ["chr1", "chr1", "chr1", "chr2", "chrM"])

        self.assertEqual(stats["added_snp"], 3)
        self.assertEqual(stats["replaced_wes_nocall"], 1)
        self.assertEqual(stats["snp_overlap_skipped"], 1)
        self.assertEqual(stats["snp_nocall_skipped"], 1)
        self.assertEqual(stats["output_records"], 5)

        text = out.read_text()
        self.assertIn("##INFO=<ID=SRC,", text)
        self.assertIn("##FORMAT=<ID=GQ,", text)
        self.assertIn("\tWES_S1\n", text)

        summary, index = parse_vcf(out)
        self.assertEqual(summary.sample_id, "WES_S1")
        self.assertEqual(index[("1", 100)]["genotype"], "0/1")
        self.assertEqual(index[("1", 200)]["genotype"], "0/1")
        self.assertEqual(index[("1", 300)]["genotype"], "1/1")

    def test_assembly_mismatch_is_rejected(self):
        snp37 = SNP.replace("##source=bcftools_gtc2vcf GSA-24v3", "##reference=GRCh37").replace(
            "length=248956422", "length=249250621"
        )
        with self.assertRaises(ValueError):
            self._merge(snp_text=snp37)


if __name__ == "__main__":
    unittest.main()
