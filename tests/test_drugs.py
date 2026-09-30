import tempfile
import unittest
from pathlib import Path

from dna_compare.comparisons.actionable import compare_actionable
from dna_compare.comparisons.additional import compare_additional
from dna_compare.comparisons.drugs import compare_drugs
from dna_compare.config import default_settings
from dna_compare.genotype_parser import parse_genotype_file
from dna_compare.service import AnalysisService
from dna_compare.vcf_parser import parse_vcf

HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
FIXTURE = Path(__file__).resolve().parents[1] / "data" / "samples" / "demo.snps.vcf"


def _vcf(*rows: str):
    with tempfile.NamedTemporaryFile("w", suffix=".vcf", delete=False) as handle:
        handle.write(HEADER + "".join(row.replace(" ", "\t") + "\n" for row in rows))
        path = Path(handle.name)
    try:
        return parse_vcf(path)[1]
    finally:
        path.unlink(missing_ok=True)


def _raw23(*rows: str):
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as handle:
        handle.write("# rsid\tchromosome\tposition\tgenotype\n" + "".join(row.replace(" ", "\t") + "\n" for row in rows))
        path = Path(handle.name)
    try:
        return parse_genotype_file(path, filename=path.name, min_raw23_rows=1)[1]
    finally:
        path.unlink(missing_ok=True)


def _drug(block, name: str):
    return next(row for row in block.estimates if row.drug.startswith(name))


def _gene(block, name: str):
    return next(call for call in block.genes if call.gene == name)


def _run(index):
    return compare_drugs(index, settings=default_settings())


class DrugResponseTests(unittest.TestCase):
    def test_cyp2c19_poor_metabolizer(self):
        block = _run(_vcf(
            "10 96541616 rs4244285 G A . PASS . GT 1/1",
            "10 96540410 rs4986893 G A . PASS . GT 0/0",
            "10 96522463 rs28399504 A G . PASS . GT 0/0",
            "10 96521657 rs12248560 C T . PASS . GT 0/0",
        ))
        self.assertEqual(_gene(block, "CYP2C19").diplotype, "*2/*2")
        self.assertEqual(_drug(block, "Clopidogrel").suitability, "not_suitable")
        self.assertEqual(_drug(block, "Omeprazole").suitability, "partial")
        self.assertEqual(_drug(block, "Sertraline").suitability, "partial")
        self.assertEqual(_drug(block, "Clopidogrel").coverage_status, "ready")

    def test_cyp2c19_rapid_metabolizer(self):
        block = _run(_vcf(
            "10 96541616 rs4244285 G A . PASS . GT 0/0",
            "10 96521657 rs12248560 C T . PASS . GT 0/1",
        ))
        call = _gene(block, "CYP2C19")
        self.assertEqual((call.phenotype, call.diplotype), ("RM", "*1/*17"))
        self.assertEqual(_drug(block, "Clopidogrel").suitability, "suitable")
        self.assertEqual(_drug(block, "Voriconazole").suitability, "not_suitable")
        self.assertIn("Not typed", call.message)
        self.assertEqual(call.coverage_status, "partial")

    def test_reference_homozygous_star4_is_not_a_carrier(self):
        block = _run(_vcf("10 96522463 rs28399504 A G . PASS . GT 0/0"))
        self.assertEqual(_gene(block, "CYP2C19").phenotype, "NM")

    def test_opposite_strand_calls_are_counted(self):
        block = _run(_raw23("rs4244285 10 96541616 TT"))
        self.assertEqual(_gene(block, "CYP2C19").phenotype, "PM")

    def test_raw23_homozygous_effect_allele(self):
        block = _run(_raw23("rs4149056 12 21331549 CC"))
        self.assertEqual(_gene(block, "SLCO1B1").phenotype, "poor")
        self.assertEqual(_drug(block, "Simvastatin").suitability, "not_suitable")
        self.assertEqual(_drug(block, "Rosuvastatin").suitability, "partial")

    def test_warfarin_needs_both_genes_to_be_suitable(self):
        typical = _run(_vcf("16 31107689 rs9923231 C T . PASS . GT 0/0"))
        warfarin = _drug(typical, "Warfarin")
        self.assertEqual(warfarin.suitability, "unknown")
        self.assertIn("CYP2C9 not called", warfarin.message)
        sensitive = _run(_vcf("16 31107689 rs9923231 C T . PASS . GT 1/1"))
        self.assertEqual(_drug(sensitive, "Warfarin").suitability, "partial")

    def test_cyp2c9_activity_score(self):
        block = _run(_vcf(
            "10 96702047 rs1799853 C T . PASS . GT 0/1",
            "10 96741053 rs1057910 A C . PASS . GT 0/1",
        ))
        call = _gene(block, "CYP2C9")
        self.assertEqual((call.phenotype, call.diplotype), ("PM", "*2/*3"))
        self.assertEqual(_drug(block, "Meloxicam").suitability, "not_suitable")

    def test_thiopurines_take_worst_gene(self):
        block = _run(_vcf(
            "6 18130918 rs1142345 T C . PASS . GT 0/1",
            "6 18139228 rs1800460 C T . PASS . GT 0/1",
            "13 48619855 rs116855232 C T . PASS . GT 1/1",
        ))
        self.assertEqual(_gene(block, "TPMT").diplotype, "*1/*3A")
        self.assertEqual(_gene(block, "NUDT15").phenotype, "PM")
        self.assertEqual(_drug(block, "Azathioprine").suitability, "not_suitable")

    def test_dpyd_intermediate(self):
        block = _run(_vcf(
            "1 97915614 rs3918290 C T . PASS . GT 0/1",
            "1 97981343 rs55886062 A C . PASS . GT 0/0",
            "1 97547947 rs67376798 T A . PASS . GT 0/0",
            "1 98045449 rs75017182 G C . PASS . GT 0/0",
        ))
        self.assertEqual(_gene(block, "DPYD").phenotype, "IM")
        self.assertEqual(_drug(block, "Fluorouracil").suitability, "partial")

    def test_hla_b5701_tag(self):
        block = _run(_vcf("6 31431780 rs2395029 T G . PASS . GT 0/1"))
        self.assertEqual(_drug(block, "Abacavir").suitability, "not_suitable")

    def test_g6pd_uses_sex(self):
        male = _run(_vcf("X 153764217 rs1050828 C T . PASS . GT 1", "Y 2655180 rs11575897 G A . PASS . GT 0"))
        self.assertEqual(_gene(male, "G6PD").phenotype, "deficient")
        self.assertEqual(_drug(male, "Rasburicase").suitability, "not_suitable")
        female = _run(_vcf("X 153764217 rs1050828 C T . PASS . GT 0/1"))
        self.assertEqual(_gene(female, "G6PD").phenotype, "variable")
        self.assertEqual(_drug(female, "Rasburicase").suitability, "partial")

    def test_cyp2d6_is_low_confidence(self):
        block = _run(_vcf(
            "22 42524947 rs3892097 C T . PASS . GT 0/1",
            "22 42526694 rs1065852 G A . PASS . GT 0/1",
            "22 42523805 rs28371725 C T . PASS . GT 0/0",
        ))
        call = _gene(block, "CYP2D6")
        self.assertEqual((call.phenotype, call.diplotype, call.confidence), ("IM", "*1/*4", "low"))
        self.assertEqual(_drug(block, "Codeine").suitability, "partial")

    def test_empty_file_is_unknown(self):
        block = _run(_vcf("1 1000 rs1 A G . PASS . GT 0/1"))
        self.assertFalse(block.available)
        self.assertTrue(all(row.suitability == "unknown" for row in block.estimates))

    def test_missing_tables(self):
        settings = default_settings()
        settings.pgx_markers = Path("/nonexistent/pgx.tsv")
        block = compare_drugs({}, settings=settings)
        self.assertFalse(block.available)
        self.assertIn("Missing drug tables", block.notes[-1])

    def test_cyp2c19_star35_is_not_counted_on_top_of_star2(self):
        both = _run(_vcf(
            "10 96541616 rs4244285 G A . PASS . GT 0/1",
            "10 96535124 rs12769205 A G . PASS . GT 0/1",
            "10 96540410 rs4986893 G A . PASS . GT 0/0",
            "10 96522463 rs28399504 A G . PASS . GT 0/0",
            "10 96521657 rs12248560 C T . PASS . GT 0/0",
        ))
        self.assertEqual(_gene(both, "CYP2C19").diplotype, "*1/*2")
        star35 = _run(_vcf(
            "10 96541616 rs4244285 G A . PASS . GT 0/0",
            "10 96535124 rs12769205 A G . PASS . GT 1/1",
            "10 96540410 rs4986893 G A . PASS . GT 0/0",
            "10 96522463 rs28399504 A G . PASS . GT 0/0",
            "10 96521657 rs12248560 C T . PASS . GT 0/0",
        ))
        self.assertEqual(_gene(star35, "CYP2C19").diplotype, "*35/*35")
        self.assertEqual(_gene(star35, "CYP2C19").phenotype, "PM")

    def test_cyp2c9_star8_lowers_activity(self):
        block = _run(_vcf(
            "10 96702047 rs1799853 C T . PASS . GT 0/0",
            "10 96741053 rs1057910 A C . PASS . GT 0/0",
            "10 96702066 rs7900194 G A . PASS . GT 1/1",
        ))
        self.assertEqual(_gene(block, "CYP2C9").diplotype, "*8/*8")
        self.assertEqual(_gene(block, "CYP2C9").phenotype, "IM1")

    def test_g6pd_orissa_is_deficient_in_a_male(self):
        block = _run(_vcf(
            "X 153764383 rs78478128 G C . PASS . GT 1",
            "Y 2655180 rs11575897 G A . PASS . GT 0",
        ))
        self.assertEqual(_gene(block, "G6PD").phenotype, "deficient")
        self.assertEqual(_drug(block, "Rasburicase").suitability, "not_suitable")

    def test_cyp2d6_star3_indel_is_a_poor_metabolizer_when_homozygous(self):
        block = _run(_vcf(
            "22 42524947 rs3892097 C T . PASS . GT 0/0",
            "22 42526694 rs1065852 G A . PASS . GT 0/0",
            "22 42523805 rs28371725 C T . PASS . GT 0/0",
            "22 42524243 rs35742686 CT C . PASS . GT 1/1",
        ))
        call = _gene(block, "CYP2D6")
        self.assertEqual(call.phenotype, "PM")
        self.assertIn("*3", call.diplotype)
        self.assertEqual(_drug(block, "Codeine").suitability, "not_suitable")

    def test_ryr1_variant_makes_volatiles_not_suitable(self):
        block = _run(_vcf("19 38948185 rs118192172 C T . PASS . GT 0/1"))
        self.assertEqual(_drug(block, "Volatile").suitability, "not_suitable")

    def test_untyped_exome_only_drugs_stay_unknown(self):
        block = _run(_vcf("1 1000 rs1 A G . PASS . GT 0/1"))
        self.assertEqual(_drug(block, "Irinotecan").suitability, "unknown")
        self.assertEqual(_drug(block, "Volatile").suitability, "unknown")
        self.assertEqual(_drug(block, "Isoniazid").suitability, "unknown")

    def test_analyze_payload_has_drugs(self):
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
        self.assertEqual(payload["drugs"]["kind"], "drugs")
        self.assertIn("genes", payload["drugs"])
        self.assertIn("suitability", payload["drugs"]["estimates"][0])
        self.assertEqual(payload["actionable"]["kind"], "actionable")

    def test_hemoglobin_and_actionable_report_hbs(self):
        index = _vcf("11 5248232 rs334 T A . PASS . GT 0/1")
        additional = compare_additional(index, settings=default_settings())
        row = next(item for item in additional.estimates if item.key == "hemoglobin")
        self.assertIn("HbS carrier", row.finding)
        actionable = compare_actionable(index, settings=default_settings())
        self.assertEqual(len(actionable.estimates), 1)
        self.assertIn("HbS", actionable.estimates[0].finding)
        self.assertEqual(compare_actionable(_vcf("11 5248232 rs334 T A . PASS . GT 0/0"), settings=default_settings()).estimates, [])
        self.assertIn("carrier", actionable.estimates[0].finding)

    def test_cftr_carrier_is_not_called_affected(self):
        het = compare_actionable(
            _vcf("7 117199645 rs113993960 TCTT T . PASS . GT 0/1"),
            settings=default_settings(),
        )
        self.assertEqual(len(het.estimates), 1)
        self.assertIn("carrier", het.estimates[0].finding)
        self.assertIn("F508del", het.estimates[0].finding)
        hom = compare_actionable(
            _vcf("7 117199645 rs113993960 TCTT T . PASS . GT 1/1"),
            settings=default_settings(),
        )
        self.assertIn("disease genotype", hom.estimates[0].finding)

    def test_mybpc3_deletion_and_g6pd_hemizygous(self):
        hcm = compare_actionable(
            _vcf("11 47353826 rs36212066 GAGAGGGAGGGAAGCCATCCAGGCTGAGAGGG GAGAGGG . PASS . GT 0/1"),
            settings=default_settings(),
        )
        self.assertEqual(len(hcm.estimates), 1)
        self.assertIn("heterozygous", hcm.estimates[0].finding)
        self.assertIn("MYBPC3", hcm.estimates[0].finding)
        male = compare_actionable(
            _vcf(
                "Y 1000 rsY A G . PASS . GT 0/0",
                "X 153764383 rs78478128 G C . PASS . GT 0/1",
            ),
            settings=default_settings(),
        )
        self.assertEqual(len(male.estimates), 1)
        self.assertIn("hemizygous", male.estimates[0].finding)

    def test_brca_founder_indels_are_reported(self):
        brca2 = compare_actionable(
            _vcf("13 32914437 rs80359550 GT G . PASS . GT 0/1"),
            settings=default_settings(),
        )
        self.assertEqual(len(brca2.estimates), 1)
        self.assertIn("6174delT", brca2.estimates[0].finding)
        self.assertIn("heterozygous", brca2.estimates[0].finding)
        brca1 = compare_actionable(
            _vcf("17 41276045 rs80357914 CTCT CT . PASS . GT 0/1"),
            settings=default_settings(),
        )
        self.assertIn("185delAG", brca1.estimates[0].finding)


if __name__ == "__main__":
    unittest.main()
