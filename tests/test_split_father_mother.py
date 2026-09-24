import importlib.util
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "split-father-mother-dna.py"


def load_splitter():
    spec = importlib.util.spec_from_file_location("split_father_mother_dna", MODULE_PATH)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


class SplitFatherMotherTests(unittest.TestCase):
    def test_phased_and_uniparental(self):
        mod = load_splitter()
        vcf = """##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	DEMO
1	100	rsA	A	G	.	PASS	.	GT	0|1
Y	200	rsY	C	T	.	PASS	.	GT	1
MT	300	rsM	A	G	.	PASS	.	GT	1
X	400	rsX	A	C	.	PASS	.	GT	0
"""
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "sample.snps.vcf"
            src.write_text(vcf)
            father, mother, stats = mod.split_vcf(src)
            self.assertEqual(father.name, "sample.snps_father.vcf")
            self.assertEqual(mother.name, "sample.snps_mother.vcf")
            ftxt = father.read_text()
            mtxt = mother.read_text()
            self.assertIn("rsY", ftxt)
            self.assertNotIn("rsM", ftxt)
            self.assertIn("rsM", mtxt)
            self.assertNotIn("rsY", mtxt)
            self.assertIn("rsX", mtxt)
            self.assertNotIn("rsX", ftxt)
            self.assertIn("DEMO_father", ftxt)
            self.assertIn("DEMO_mother", mtxt)
            father_gt = [ln.split("\t")[-1] for ln in ftxt.splitlines() if ln.startswith("1\t")][0]
            mother_gt = [ln.split("\t")[-1] for ln in mtxt.splitlines() if ln.startswith("1\t")][0]
            self.assertEqual(father_gt, "0/0")
            self.assertEqual(mother_gt, "1/1")
            self.assertEqual(stats["phased"], 1)
            self.assertEqual(stats["haploid"], 3)

    def test_unphased_file_order(self):
        mod = load_splitter()
        vcf = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
2	10	rsHet	G	A	.	PASS	.	GT	0/1
"""
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "demo.vcf"
            src.write_text(vcf)
            father, mother, _stats = mod.split_vcf(src)
            self.assertEqual(father.name, "demo_father.vcf")
            self.assertEqual(mother.name, "demo_mother.vcf")
            self.assertIn("0/0", father.read_text())
            self.assertIn("1/1", mother.read_text())


if __name__ == "__main__":
    unittest.main()
