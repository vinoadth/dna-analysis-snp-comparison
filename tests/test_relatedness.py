import unittest
from pathlib import Path
from io import StringIO

from dna_compare.comparisons.relatedness import classify_relationship, compare_relatedness
from dna_compare.models import AnalysisResult
from dna_compare.service import AnalysisService
from dna_compare.vcf_parser import parse_vcf

SAMPLES = Path(__file__).resolve().parents[1] / "data" / "samples"


def _row(chrom: str, pos: int, ref: str, alt: str, genotype: str, dosage: float) -> dict:
    return {
        "chrom": chrom,
        "pos": pos,
        "rsid": f"rs{pos}",
        "ref": ref,
        "alt": alt,
        "genotype": genotype,
        "dosage_alt": dosage,
        "is_snp": True,
    }


def _gt(dosage: float) -> tuple[str, float]:
    return {0.0: ("0/0", 0.0), 1.0: ("0/1", 1.0), 2.0: ("1/1", 2.0)}[dosage]


def _index(dosages: list[float], *, chrom: str = "1", start: int = 1000) -> dict[tuple[str, int], dict]:
    index = {}
    for i, dosage in enumerate(dosages):
        pos = start + i
        gt, dose = _gt(dosage)
        index[(chrom, pos)] = _row(chrom, pos, "A", "G", gt, dose)
    return index


def _vcf(sample: str, dosages: list[float]) -> StringIO:
    lines = [
        "##fileformat=VCFv4.2",
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}",
    ]
    for i, dosage in enumerate(dosages):
        gt, _dose = _gt(dosage)
        lines.append(f"1\t{1000 + i}\trs{1000 + i}\tA\tG\t.\tPASS\t.\tGT\t{gt}")
    return StringIO("\n".join(lines) + "\n")


class RelatednessTests(unittest.TestCase):
    def test_identical_files_are_duplicate_kinship(self):
        dosages = ([0.0, 1.0, 2.0] * 90)[:250]
        left = _index(dosages)
        result = compare_relatedness(left, dict(left), other_filename="same.vcf", other_sample_id="TWIN")
        self.assertTrue(result.available)
        self.assertEqual(result.n_snps, 250)
        self.assertEqual(result.ibs0, 0)
        self.assertAlmostEqual(result.kinship, 0.5, places=3)
        self.assertAlmostEqual(result.shared_pct, 100.0, places=1)
        self.assertAlmostEqual(result.ibs0_pct, 0.0, places=1)
        self.assertGreater(result.ibs2_pct, 50.0)
        self.assertIn("same person", result.relationship)

    def test_parent_child_has_no_ibs0_and_first_degree_kinship(self):
        rng = __import__("random").Random(7)
        parent = []
        child = []
        p = 0.3
        for _ in range(400):
            p1 = (int(rng.random() < p), int(rng.random() < p))
            p2 = (int(rng.random() < p), int(rng.random() < p))
            transmitted = p1[rng.randrange(2)]
            other = p2[rng.randrange(2)]
            parent.append(float(sum(p1)))
            child.append(float(transmitted + other))
        result = compare_relatedness(
            _index(parent),
            _index(child),
            other_filename="parent.vcf",
            other_sample_id="PARENT",
        )
        self.assertTrue(result.available)
        self.assertEqual(result.ibs0, 0)
        self.assertGreater(result.kinship, 0.177)
        self.assertLess(result.kinship, 0.354)
        self.assertIn("parent–child", result.relationship)

    def test_unrelated_pair_is_not_close(self):
        left = ([0.0, 0.0, 2.0, 2.0] * 70)[:250]
        right = ([2.0, 2.0, 0.0, 0.0] * 70)[:250]
        result = compare_relatedness(_index(left), _index(right), other_filename="friend.vcf")
        self.assertTrue(result.available)
        self.assertGreater(result.ibs0, 200)
        self.assertLess(result.kinship, 0.022)
        self.assertIn("unrelated", result.relationship)

    def test_strand_swap_still_matches(self):
        left = _index([0.0, 1.0, 2.0] * 80)
        right = {}
        for key, row in left.items():
            flipped = dict(row)
            flipped["ref"] = row["alt"]
            flipped["alt"] = row["ref"]
            flipped["dosage_alt"] = 2.0 - row["dosage_alt"]
            flipped["genotype"] = _gt(flipped["dosage_alt"])[0]
            right[key] = flipped
        result = compare_relatedness(left, right, other_filename="swap.vcf")
        self.assertTrue(result.available)
        self.assertEqual(result.ibs0, 0)
        self.assertAlmostEqual(result.kinship, 0.5, places=3)

    def test_too_few_overlapping_snps_is_unavailable(self):
        result = compare_relatedness(_index([0.0, 1.0, 2.0]), _index([0.0, 1.0, 2.0]), other_filename="tiny.vcf")
        self.assertFalse(result.available)
        self.assertEqual(result.n_snps, 3)
        self.assertIn("overlapping autosomal SNPs", result.notes[-1])

    def test_sex_chromosomes_are_ignored(self):
        auto = _index([1.0] * 250)
        left = dict(auto)
        right = dict(auto)
        left[("Y", 1)] = _row("Y", 1, "A", "G", "1", 1.0)
        right[("Y", 1)] = _row("Y", 1, "A", "G", "0", 0.0)
        result = compare_relatedness(left, right, other_filename="y.vcf")
        self.assertEqual(result.n_snps, 250)
        self.assertEqual(result.ibs0, 0)

    def test_classify_relationship_bins(self):
        self.assertIn("same person", classify_relationship(0.5, 0, 1000))
        self.assertIn("parent–child", classify_relationship(0.25, 0, 1000))
        self.assertIn("first-degree", classify_relationship(0.25, 20, 1000))
        self.assertIn("second-degree", classify_relationship(0.12, 10, 1000))
        self.assertIn("third-degree", classify_relationship(0.06, 10, 1000))
        self.assertIn("fourth-degree", classify_relationship(0.03, 10, 1000))
        self.assertIn("unrelated", classify_relationship(0.0, 100, 1000))
        self.assertIn("not enough", classify_relationship(0.5, 0, 10))
        self.assertIn("parent–child", classify_relationship(0.25, 4, 1000))

    def test_rsid_matches_when_coordinates_differ(self):
        dosages = ([0.0, 1.0, 2.0] * 90)[:250]
        left = _index(dosages)
        right = {}
        for (chrom, pos), row in left.items():
            moved = dict(row)
            moved["pos"] = pos + 50_000_000
            right[(chrom, pos + 50_000_000)] = moved
        result = compare_relatedness(left, right, other_filename="hg38.vcf")
        self.assertTrue(result.available)
        self.assertEqual(result.n_matched_rsid, 250)
        self.assertEqual(result.n_matched_pos, 0)
        self.assertAlmostEqual(result.kinship, 0.5, places=3)
        self.assertTrue(any("rsID" in note for note in result.notes))

    def test_haploid_calls_are_treated_as_homozygotes(self):
        left = {}
        for i in range(250):
            token = "0" if i % 2 == 0 else "1"
            left[("1", 1000 + i)] = _row("1", 1000 + i, "A", "G", token, 0.0 if token == "0" else 1.0)
        result = compare_relatedness(left, dict(left), other_filename="haploid.vcf")
        self.assertTrue(result.available)
        self.assertEqual(result.ibs0, 0)
        self.assertAlmostEqual(result.kinship, 0.5, places=3)

    def test_low_quality_sites_are_dropped(self):
        dosages = ([0.0, 1.0, 2.0] * 90)[:250]
        left = _index(dosages)
        right = dict(left)
        for i in range(80):
            pos = 500_000 + i
            bad_left = _row("1", pos, "A", "G", "0/0", 0.0)
            bad_right = _row("1", pos, "A", "G", "1/1", 2.0)
            bad_left["igc"] = 0.05
            bad_right["igc"] = 0.05
            left[("1", pos)] = bad_left
            right[("1", pos)] = bad_right
        result = compare_relatedness(left, right, other_filename="noisy.vcf")
        self.assertTrue(result.available)
        self.assertEqual(result.n_qc_dropped, 80)
        self.assertEqual(result.ibs0, 0)
        self.assertAlmostEqual(result.kinship, 0.5, places=3)

    def test_analyze_without_other_file_marks_unavailable(self):
        dosages = [0.0, 1.0, 2.0] * 80
        service = AnalysisService()
        result = service.analyze(
            _vcf("DEMO1", dosages),
            filename="demo.vcf",
            compare_hominin_flag=False,
            compare_caste_flag=False,
            compare_populations_flag=False,
            compare_ancestry_flag=False,
            compare_haplogroups_flag=False,
        )
        self.assertIsInstance(result, AnalysisResult)
        self.assertFalse(result.relatedness.available)
        payload = result.to_dict()
        self.assertIn("relatedness", payload)
        self.assertFalse(payload["relatedness"]["available"])
        self.assertIn("second VCF", payload["relatedness"]["notes"][0])

    def test_analyze_with_other_file_reports_kinship(self):
        dosages = ([0.0, 1.0, 2.0] * 90)[:250]
        service = AnalysisService()
        result = service.analyze(
            _vcf("CHILD", dosages),
            filename="child.vcf",
            other_source=_vcf("COUSIN", dosages),
            other_filename="cousin.vcf",
            compare_hominin_flag=False,
            compare_caste_flag=False,
            compare_populations_flag=False,
            compare_ancestry_flag=False,
            compare_haplogroups_flag=False,
        )
        self.assertTrue(result.relatedness.available)
        self.assertEqual(result.relatedness.other_filename, "cousin.vcf")
        self.assertEqual(result.relatedness.other_sample_id, "COUSIN")
        self.assertAlmostEqual(result.relatedness.kinship, 0.5, places=3)
        self.assertIn("same person", result.relatedness.relationship)
        payload = result.to_dict()
        self.assertIn("relatedness", payload)
        self.assertAlmostEqual(payload["relatedness"]["kinship"], 0.5, places=3)

    def test_bundled_demo_and_demo2_are_first_cousins(self):
        demo = SAMPLES / "demo.snps.vcf"
        demo2 = SAMPLES / "demo2.snps.vcf"
        self.assertTrue(demo.is_file())
        self.assertTrue(demo2.is_file())
        _, left = parse_vcf(demo)
        _, right = parse_vcf(demo2)
        result = compare_relatedness(left, right, other_filename=demo2.name, other_sample_id="DEMO2")
        self.assertTrue(result.available)
        self.assertGreaterEqual(result.n_snps, 200)
        self.assertGreaterEqual(result.kinship, 0.044)
        self.assertLess(result.kinship, 0.088)
        self.assertIn("cousin", result.relationship)
        service = AnalysisService()
        payload = service.analyze(
            demo,
            filename=demo.name,
            other_source=demo2,
            other_filename=demo2.name,
            compare_hominin_flag=False,
            compare_caste_flag=False,
            compare_populations_flag=False,
            compare_ancestry_flag=False,
            compare_haplogroups_flag=False,
        ).to_dict()
        self.assertTrue(payload["relatedness"]["available"])
        self.assertIn("cousin", payload["relatedness"]["relationship"])


if __name__ == "__main__":
    unittest.main()
