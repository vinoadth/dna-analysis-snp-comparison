import unittest
from io import BytesIO
from pathlib import Path

from dna_compare.genotype_parser import parse_genotype_file
from dna_compare.raw23_parser import Raw23ParseError, parse_raw23, validate_raw23
from dna_compare.service import AnalysisService

FIXTURE = Path(__file__).resolve().parents[1] / "data" / "samples" / "demo.raw23.txt"
VCF_FIXTURE = Path(__file__).resolve().parents[1] / "data" / "samples" / "demo.snps.vcf"

VALID_HEADER = "# rsid\tchromosome\tposition\tgenotype\n"


def _mini_raw23(n_rows: int = 12, *, include_header: bool = True) -> str:
    lines = []
    if include_header:
        lines.append(VALID_HEADER)
    for i in range(n_rows):
        pos = 1_000_000 + i
        gt = "AG" if i % 3 == 0 else "TT"
        lines.append(f"rs{i:06d}\t1\t{pos}\t{gt}\n")
    return "".join(lines)


class Raw23ParserTests(unittest.TestCase):
    def test_parse_demo_fixture(self):
        summary, index = parse_raw23(FIXTURE, min_snp_rows=1000)
        self.assertGreaterEqual(summary.n_snps, 1000)
        self.assertEqual(summary.source, "23andMe-compatible raw data")
        self.assertEqual(summary.assembly, "GRCh37")
        self.assertIn(("1", 567667), index)
        self.assertEqual(index[("1", 567667)]["genotype"], "0/0")
        self.assertEqual(index[("1", 567667)]["ref"], "C")
        self.assertIn(("1", 838555), index)
        self.assertEqual(index[("1", 838555)]["genotype"], "0/1")

    def test_genotype_file_routes_txt_by_filename(self):
        payload = FIXTURE.read_bytes()
        summary, index = parse_genotype_file(BytesIO(payload), filename="kit.txt", min_raw23_rows=1000)
        self.assertGreaterEqual(summary.n_snps, 1000)
        self.assertIn(("1", 567667), index)

    def test_rejects_vcf_as_raw23(self):
        payload = VCF_FIXTURE.read_bytes()
        with self.assertRaises(Raw23ParseError) as ctx:
            parse_genotype_file(BytesIO(payload), filename="kit.txt", min_raw23_rows=5)
        self.assertIn("looks like a VCF", str(ctx.exception))

    def test_rejects_missing_header_and_bad_columns(self):
        bad = "rs1 1 100 AA\n" * 5
        with self.assertRaises(Raw23ParseError) as ctx:
            validate_raw23(BytesIO(bad.encode()), min_snp_rows=3)
        self.assertIn("header", str(ctx.exception).lower())

    def test_rejects_space_separated_rows(self):
        body = VALID_HEADER + "rs123 1 1000 AG\n" + _mini_raw23(20, include_header=False)
        with self.assertRaises(Raw23ParseError) as ctx:
            validate_raw23(BytesIO(body.encode()), min_snp_rows=10)
        self.assertIn("tab-separated", str(ctx.exception))

    def test_rejects_too_few_rows(self):
        tiny = _mini_raw23(12)
        with self.assertRaises(Raw23ParseError) as ctx:
            validate_raw23(BytesIO(tiny.encode()), min_snp_rows=1000)
        self.assertIn("too few SNP rows", str(ctx.exception))

    def test_skips_unsupported_genotype_rows(self):
        body = VALID_HEADER + "rs123\t1\t1000\tXY\n" + _mini_raw23(5, include_header=False)
        with self.assertRaises(Raw23ParseError) as ctx:
            validate_raw23(BytesIO(body.encode()), min_snp_rows=10)
        self.assertIn("too few SNP rows", str(ctx.exception))

    def test_keeps_blood_type_indels_on_autosomes(self):
        body = (
            VALID_HEADER
            + "rs8176719\t9\t136132908\tDD\n"
            + "i4001527\t1\t25600000\tDI\n"
            + _mini_raw23(12, include_header=False)
        )
        _summary, index = parse_raw23(BytesIO(body.encode()), min_snp_rows=10)
        self.assertIn(("9", 136132908), index)
        self.assertEqual(index[("9", 136132908)]["genotype"], "1/1")
        self.assertIn(("1", 25600000), index)
        self.assertEqual(index[("1", 25600000)]["genotype"], "0/1")

    def test_service_returns_clear_error_for_invalid_raw(self):
        result = AnalysisService().analyze(
            BytesIO(b"not a genotype file\n"),
            filename="broken.txt",
        )
        self.assertFalse(result.ok)
        self.assertTrue(result.errors)
        self.assertIn("23andMe-compatible", result.errors[0])

    def test_service_analyzes_valid_raw_fixture(self):
        result = AnalysisService().analyze(
            FIXTURE,
            filename="demo.raw23.txt",
            compare_hominin_flag=False,
            compare_caste_flag=False,
            compare_populations_flag=False,
            compare_ancestry_flag=False,
        )
        self.assertTrue(result.ok, result.errors)
        self.assertGreaterEqual(result.vcf.n_snps, 1000)


if __name__ == "__main__":
    unittest.main()
