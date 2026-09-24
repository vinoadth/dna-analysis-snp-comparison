import unittest

from dna_compare.comparisons.community_ref import (
    bengali_reference_applicable,
    collapse_ref_y,
    gujarati_reference_applicable,
    interval_score,
    load_bengali_community_reference,
    load_community_panel_aliases,
    load_gujarati_community_reference,
    load_kerala_community_reference,
    load_marathi_community_reference,
    load_punjabi_community_reference,
    load_tamil_community_reference,
    kerala_reference_applicable,
    marathi_reference_applicable,
    parse_y_haplogroups,
    punjabi_reference_applicable,
    score_community_reference,
    tamil_reference_applicable,
)
from dna_compare.models import ComparisonBlock, HaplogroupResult, PopulationEstimate


class CommunityRefTests(unittest.TestCase):
    def test_loads_published_tamil_rows(self):
        rows = load_tamil_community_reference()
        names = [row.display_name for row in rows]
        self.assertIn("Tamil Brahmin (Iyer)", names)
        self.assertIn("Vanniyar", names)
        self.assertIn("Irula (Tribal)", names)
        vanniyar = next(row for row in rows if row.community_id == "vanniyar")
        self.assertEqual(vanniyar.steppe_min, 6)
        self.assertEqual(vanniyar.steppe_max, 12)
        self.assertEqual(vanniyar.y_haplogroups["H-M69"], 28)

    def test_parse_and_collapse_y(self):
        parsed = parse_y_haplogroups("R1a:25;H-M69:18;L-M20:15;J2:20")
        self.assertEqual(parsed["R1a"], 25)
        self.assertEqual(collapse_ref_y("R1a1 (M17)"), "R1a")
        self.assertEqual(collapse_ref_y("H3b (Z13871)"), "H-M69")
        self.assertEqual(collapse_ref_y("J2b (M241)"), "J2")
        self.assertEqual(collapse_ref_y("O2a (M175)"), "O-M175")
        self.assertIsNone(collapse_ref_y("other"))

    def test_interval_score_inside_and_outside(self):
        self.assertEqual(interval_score(10.9, 6, 12), 1.0)
        self.assertLess(interval_score(10.9, 18, 25), 0.6)
        self.assertEqual(interval_score(1.0, 0, 2), 1.0)
        self.assertEqual(interval_score(20.0, 0, 2), 0.0)

    def test_steppe_11_ranks_vanniyar_above_iyer_and_irula(self):
        ancestry = ComparisonBlock(
            kind="ancestry",
            available=True,
            estimates=[
                PopulationEstimate(population="AASI_Onge", percent=49.0, n_snps=1000),
                PopulationEstimate(population="Steppe_MLBA", percent=10.9, n_snps=1000),
                PopulationEstimate(population="Indus_Periphery", percent=40.1, n_snps=1000),
            ],
        )
        haplo = HaplogroupResult(available=True, sample_best=None)
        caste = ComparisonBlock(
            kind="caste",
            available=True,
            estimates=[
                PopulationEstimate(population="Tamil", percent=22.0, n_snps=1000),
                PopulationEstimate(population="Bengali", percent=21.0, n_snps=1000),
            ],
        )
        result = score_community_reference(ancestry, haplo, caste=caste)
        self.assertTrue(result.available)
        ranked = [est.population for est in result.estimates]
        self.assertLess(ranked.index("Vanniyar"), ranked.index("Tamil Brahmin (Iyer)"))
        self.assertLess(ranked.index("Vanniyar"), ranked.index("Irula (Tribal)"))
        vanniyar = next(est for est in result.estimates if est.population == "Vanniyar")
        self.assertTrue(vanniyar.steppe_in_range)
        self.assertFalse(vanniyar.aasi_in_range)
        iyer = next(est for est in result.estimates if "Iyer" in est.population)
        self.assertFalse(iyer.steppe_in_range)

    def test_loads_punjabi_rows(self):
        rows = load_punjabi_community_reference()
        names = [row.display_name for row in rows]
        self.assertIn("Khatri", names)
        self.assertIn("Jat", names)
        khatri = next(row for row in rows if row.community_id == "khatri")
        self.assertEqual(khatri.steppe_min, 25)
        self.assertEqual(khatri.indus_min, 22)
        self.assertEqual(khatri.aasi_max, 38)
        self.assertEqual(khatri.note, "")

    def test_collapse_uniform_notes(self):
        rows = load_punjabi_community_reference()
        self.assertTrue(all(row.note == "" for row in rows))
        tamil = load_tamil_community_reference()
        self.assertTrue(any(row.note for row in tamil))
        self.assertGreater(len({row.note for row in tamil if row.note}), 1)

    def test_loads_bengali_rows(self):
        rows = load_bengali_community_reference()
        names = [row.display_name for row in rows]
        self.assertIn("Kayastha", names)
        self.assertIn("Santhal (tribal)", names)
        kayastha = next(row for row in rows if row.community_id == "kayastha")
        self.assertEqual(kayastha.steppe_min, 48)
        self.assertEqual(kayastha.east_asian_min, 7)
        self.assertEqual(kayastha.y_haplogroups["R1a"], 40)

    def test_bengali_panel_scores_with_east_asian(self):
        ancestry = ComparisonBlock(
            kind="ancestry",
            available=True,
            estimates=[
                PopulationEstimate(population="AASI_Onge", percent=36.0, n_snps=1000),
                PopulationEstimate(population="Steppe_MLBA", percent=52.0, n_snps=1000),
                PopulationEstimate(population="Indus_Periphery", percent=12.0, n_snps=1000),
            ],
        )
        ancestry_5 = ComparisonBlock(
            kind="ancestry_5",
            available=True,
            estimates=[
                PopulationEstimate(population="East_Asian", percent=9.0, n_snps=1000),
            ],
        )
        caste = ComparisonBlock(
            kind="caste",
            available=True,
            estimates=[
                PopulationEstimate(population="Bengali", percent=45.0, n_snps=1000),
            ],
        )
        result = score_community_reference(
            ancestry,
            HaplogroupResult(available=True, sample_best="R1a1"),
            ancestry_5=ancestry_5,
            caste=caste,
        )
        self.assertTrue(result.available)
        panels = {row.panel for row in result.estimates}
        self.assertIn("Bengali community reference ranges", panels)
        kayastha = next(row for row in result.estimates if row.population == "Kayastha")
        self.assertTrue(kayastha.aasi_in_range)
        self.assertTrue(kayastha.steppe_in_range)
        self.assertTrue(kayastha.east_asian_in_range)
        self.assertEqual(kayastha.ref_aasi_label, "Ref ASI")

    def test_loads_gujarati_rows(self):
        rows = load_gujarati_community_reference()
        names = [row.display_name for row in rows]
        self.assertIn("Anavil Brahmin (south Gujarat)", names)
        self.assertIn("Leva / Kadva Patel", names)
        self.assertIn("Parsi", names)
        anavil = next(row for row in rows if row.community_id == "anavil_brahmin")
        self.assertEqual(anavil.steppe_min, 65)
        self.assertEqual(anavil.steppe_max, 72)
        self.assertEqual(anavil.aasi_min, 28)
        self.assertEqual(anavil.aasi_max, 35)
        parsi = next(row for row in rows if row.community_id == "parsi")
        self.assertIn("Iranian", parsi.note)

    def test_gujarati_panel_scores_ani_asi(self):
        ancestry = ComparisonBlock(
            kind="ancestry",
            available=True,
            estimates=[
                PopulationEstimate(population="AASI_Onge", percent=32.0, n_snps=1000),
                PopulationEstimate(population="Steppe_MLBA", percent=68.0, n_snps=1000),
                PopulationEstimate(population="Indus_Periphery", percent=0.0, n_snps=1000),
            ],
        )
        caste = ComparisonBlock(
            kind="caste",
            available=True,
            estimates=[
                PopulationEstimate(population="Gujarati", percent=45.0, n_snps=1000),
            ],
        )
        result = score_community_reference(
            ancestry,
            HaplogroupResult(available=True, sample_best=None),
            caste=caste,
        )
        self.assertTrue(result.available)
        panels = {row.panel for row in result.estimates}
        self.assertIn("Gujarati community reference ranges", panels)
        anavil = next(row for row in result.estimates if row.population == "Anavil Brahmin (south Gujarat)")
        self.assertTrue(anavil.aasi_in_range)
        self.assertTrue(anavil.steppe_in_range)
        self.assertEqual(anavil.ref_aasi_label, "Ref ASI")

    def test_loads_marathi_rows(self):
        rows = load_marathi_community_reference()
        names = [row.display_name for row in rows]
        self.assertIn("Maratha", names)
        self.assertIn("Deshastha Brahmin", names)
        self.assertIn("Mahar", names)
        maratha = next(row for row in rows if row.community_id == "maratha")
        self.assertEqual(maratha.steppe_min, 40)
        self.assertEqual(maratha.aasi_max, 55)
        self.assertEqual(maratha.y_haplogroups["R1a"], 12)

    def test_panel_aliases_include_marathi_and_tamil_groups(self):
        aliases = load_community_panel_aliases()
        self.assertIn("marathi", aliases["marathi"])
        self.assertIn("malayalam", aliases["kerala"])
        self.assertIn("chettiyar", aliases["tamil"])
        self.assertTrue(marathi_reference_applicable(None, "pune_kit.vcf"))
        self.assertTrue(tamil_reference_applicable(None, "nattukottai_chettiyar.vcf"))

    def test_marathi_panel_via_brahmin_bar_and_region_filename(self):
        caste = ComparisonBlock(
            kind="caste",
            available=True,
            estimates=[
                PopulationEstimate(population="Brahmin", percent=35.0, n_snps=100),
                PopulationEstimate(population="Gujarati", percent=25.0, n_snps=100),
            ],
        )
        self.assertFalse(marathi_reference_applicable(caste, "kit.vcf"))
        self.assertTrue(marathi_reference_applicable(caste, "mumbai_sample.vcf"))

    def test_tamil_panel_via_vellalar_ho_bar(self):
        caste = ComparisonBlock(
            kind="caste",
            available=True,
            estimates=[
                PopulationEstimate(population="Vellalar", percent=35.0, n_snps=100),
                PopulationEstimate(population="Gujarati", percent=25.0, n_snps=100),
            ],
        )
        self.assertTrue(tamil_reference_applicable(caste, "kit.vcf"))

    def test_loads_kerala_rows(self):
        rows = load_kerala_community_reference()
        names = [row.display_name for row in rows]
        self.assertIn("Namboothiri Brahmin", names)
        self.assertIn("Nair", names)
        self.assertIn("Ezhava / Thiyya", names)
        nair = next(row for row in rows if row.community_id == "nair")
        self.assertEqual(nair.steppe_min, 40)
        self.assertEqual(nair.aasi_max, 55)
        self.assertEqual(nair.y_haplogroups["H-M69"], 30)
        syrian = next(row for row in rows if row.community_id == "syrian_christian")
        self.assertIn("Cochin Jewish", syrian.note)

    def test_kerala_panel_via_filename_and_brahmin_proxy(self):
        self.assertTrue(kerala_reference_applicable(None, "malayalam_nair.vcf"))
        caste = ComparisonBlock(
            kind="caste",
            available=True,
            estimates=[
                PopulationEstimate(population="Brahmin", percent=35.0, n_snps=100),
                PopulationEstimate(population="Gujarati", percent=25.0, n_snps=100),
            ],
        )
        self.assertFalse(kerala_reference_applicable(caste, "kit.vcf"))
        self.assertTrue(kerala_reference_applicable(caste, "kochi_sample.vcf"))

    def test_kerala_panel_scores_ani_asi(self):
        ancestry = ComparisonBlock(
            kind="ancestry",
            available=True,
            estimates=[
                PopulationEstimate(population="AASI_Onge", percent=48.0, n_snps=1000),
                PopulationEstimate(population="Steppe_MLBA", percent=45.0, n_snps=1000),
                PopulationEstimate(population="Indus_Periphery", percent=7.0, n_snps=1000),
            ],
        )
        result = score_community_reference(
            ancestry,
            HaplogroupResult(available=True, sample_best="H-M69"),
            caste=None,
            filename="kerala_ezhava.vcf",
        )
        self.assertTrue(result.available)
        panels = {row.panel for row in result.estimates}
        self.assertIn("Kerala community reference ranges", panels)
        nair = next(row for row in result.estimates if row.population == "Nair")
        self.assertTrue(nair.aasi_in_range)
        self.assertTrue(nair.steppe_in_range)

    def test_marathi_panel_scores_ani_asi(self):
        ancestry = ComparisonBlock(
            kind="ancestry",
            available=True,
            estimates=[
                PopulationEstimate(population="AASI_Onge", percent=52.0, n_snps=1000),
                PopulationEstimate(population="Steppe_MLBA", percent=42.0, n_snps=1000),
                PopulationEstimate(population="Indus_Periphery", percent=6.0, n_snps=1000),
            ],
        )
        result = score_community_reference(
            ancestry,
            HaplogroupResult(available=True, sample_best="R1a1"),
            caste=None,
            filename="maratha_pune.vcf",
        )
        self.assertTrue(result.available)
        panels = {row.panel for row in result.estimates}
        self.assertIn("Marathi community reference ranges", panels)
        maratha = next(row for row in result.estimates if row.population == "Maratha")
        self.assertTrue(maratha.aasi_in_range)
        self.assertTrue(maratha.steppe_in_range)

    def test_punjabi_panel_scores_with_indus(self):
        ancestry = ComparisonBlock(
            kind="ancestry",
            available=True,
            estimates=[
                PopulationEstimate(population="AASI_Onge", percent=32.0, n_snps=1000),
                PopulationEstimate(population="Steppe_MLBA", percent=28.0, n_snps=1000),
                PopulationEstimate(population="Indus_Periphery", percent=26.0, n_snps=1000),
            ],
        )
        caste = ComparisonBlock(
            kind="caste",
            available=True,
            estimates=[
                PopulationEstimate(population="Punjabi", percent=40.0, n_snps=1000),
                PopulationEstimate(population="Gujarati", percent=20.0, n_snps=1000),
            ],
        )
        result = score_community_reference(
            ancestry,
            HaplogroupResult(available=True, sample_best=None),
            caste=caste,
        )
        self.assertTrue(result.available)
        panels = {row.panel for row in result.estimates}
        self.assertIn("Punjabi community reference ranges", panels)
        khatri = next(row for row in result.estimates if row.population == "Khatri")
        self.assertTrue(khatri.aasi_in_range)
        self.assertTrue(khatri.steppe_in_range)
        self.assertTrue(khatri.indus_in_range)

    def test_hidden_when_file_is_not_tamil(self):
        self.assertTrue(tamil_reference_applicable(None, "C8XY_iyer.vcf"))
        self.assertFalse(tamil_reference_applicable(None, "demo.snps.vcf"))
        north = ComparisonBlock(
            kind="caste",
            available=True,
            estimates=[
                PopulationEstimate(population="Punjabi", percent=40.0, n_snps=100),
                PopulationEstimate(population="Gujarati", percent=30.0, n_snps=100),
                PopulationEstimate(population="Marathi", percent=15.0, n_snps=100),
                PopulationEstimate(population="Tamil", percent=5.0, n_snps=100),
                PopulationEstimate(population="Bengali", percent=5.0, n_snps=100),
            ],
        )
        self.assertFalse(tamil_reference_applicable(north, "kit.vcf"))
        self.assertTrue(punjabi_reference_applicable(north, "kit.vcf"))
        self.assertFalse(bengali_reference_applicable(north, "kit.vcf"))
        bengali = ComparisonBlock(
            kind="caste",
            available=True,
            estimates=[
                PopulationEstimate(population="Bengali", percent=40.0, n_snps=100),
                PopulationEstimate(population="Gujarati", percent=30.0, n_snps=100),
            ],
        )
        self.assertTrue(bengali_reference_applicable(bengali, "kit.vcf"))
        gujarati = ComparisonBlock(
            kind="caste",
            available=True,
            estimates=[
                PopulationEstimate(population="Gujarati", percent=40.0, n_snps=100),
                PopulationEstimate(population="Marathi", percent=30.0, n_snps=100),
            ],
        )
        self.assertTrue(gujarati_reference_applicable(gujarati, "kit.vcf"))
        third = ComparisonBlock(
            kind="caste",
            available=True,
            estimates=[
                PopulationEstimate(population="Bengali", percent=22.0, n_snps=100),
                PopulationEstimate(population="Telugu", percent=20.0, n_snps=100),
                PopulationEstimate(population="Tamil", percent=18.0, n_snps=100),
                PopulationEstimate(population="Punjabi", percent=10.0, n_snps=100),
            ],
        )
        self.assertTrue(tamil_reference_applicable(third, "kit.vcf"))
        self.assertFalse(punjabi_reference_applicable(third, "kit.vcf"))
        hidden = score_community_reference(
            ComparisonBlock(kind="ancestry", available=True, estimates=[]),
            HaplogroupResult(available=False),
            caste=ComparisonBlock(
                kind="caste",
                available=True,
                estimates=[
                    PopulationEstimate(population="Marathi", percent=40.0, n_snps=100),
                    PopulationEstimate(population="Rajasthani", percent=30.0, n_snps=100),
                ],
            ),
            filename="kit.vcf",
        )
        self.assertTrue(hidden.hidden)
        self.assertFalse(hidden.available)
        punjabi_only = score_community_reference(
            ComparisonBlock(kind="ancestry", available=True, estimates=[]),
            HaplogroupResult(available=False),
            caste=north,
            filename="kit.vcf",
        )
        self.assertFalse(punjabi_only.hidden)


if __name__ == "__main__":
    unittest.main()
