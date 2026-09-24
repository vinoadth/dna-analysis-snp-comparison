import unittest
from io import BytesIO
from pathlib import Path

from dna_compare.comparisons.caste import caste_alias_notes
from dna_compare.models import AnalysisResult
from dna_compare.service import AnalysisService
from dna_compare.vcf_parser import parse_vcf

FIXTURE = Path(__file__).resolve().parents[1] / "data" / "samples" / "demo.snps.vcf"


class VcfAndApiTests(unittest.TestCase):
    def test_parse_skips_indels_and_normalizes_chrom(self):
        summary, index = parse_vcf(FIXTURE)
        self.assertEqual(summary.sample_id, "DEMO1")
        self.assertGreaterEqual(summary.n_snps, 200)
        self.assertEqual(summary.n_non_snp_skipped, 1)
        self.assertIn(("1", 752566), index)
        self.assertIn(("3", 12345), index)
        self.assertEqual(index[("1", 752566)]["dosage_alt"], 1.0)
        self.assertEqual(index[("1", 891021)]["dosage_alt"], 2.0)
        self.assertIsNone(index[("1", 752566)]["gq"])
        self.assertIsNone(index[("1", 752566)]["dp"])
        self.assertEqual(summary.preview[0].rsid, "rs3094315")

    def test_parse_from_bytes_for_future_api_upload(self):
        payload = FIXTURE.read_bytes()
        summary, index = parse_vcf(BytesIO(payload), preview_limit=3)
        self.assertGreaterEqual(summary.n_snps, 200)
        self.assertEqual(len(summary.preview), 3)
        self.assertIn(("X", 2699624), index)

    def test_analyze_json_shape(self):
        service = AnalysisService()
        result = service.analyze(
            FIXTURE,
            filename="demo.snps.vcf",
            compare_hominin_flag=False,
            compare_caste_flag=False,
            compare_populations_flag=False,
            compare_ancestry_flag=False,
        )
        self.assertIsInstance(result, AnalysisResult)
        payload = result.to_dict()
        self.assertTrue(payload["ok"])
        self.assertEqual(payload["source_filename"], "demo.snps.vcf")
        self.assertIn("preview", payload["vcf"])
        self.assertTrue(payload["vcf"].get("snp_catalog_id"))
        self.assertGreater(payload["vcf"]["n_catalog"], 200)
        preview_chroms = {row["chrom"] for row in payload["vcf"]["preview"]}
        self.assertIn("1", preview_chroms)
        self.assertIn("X", preview_chroms)
        self.assertGreaterEqual(len(preview_chroms), 20)
        self.assertEqual(payload["hominin"]["kind"], "hominin")
        self.assertEqual(payload["caste"]["kind"], "caste")
        self.assertEqual(payload["populations"]["kind"], "populations")
        self.assertEqual(payload["ancestry"]["kind"], "ancestry")
        self.assertIn("ancestry_5", payload)
        self.assertEqual(payload["ancestry_5"]["kind"], "ancestry_5")
        self.assertFalse(payload["ancestry_5"]["available"])
        self.assertIn("community_ref", payload)
        self.assertTrue(payload["community_ref"].get("hidden"))
        self.assertIn("haplogroups", payload)
        self.assertIn("markers", payload["haplogroups"])
        self.assertIn("rows", payload["haplogroups"])
        self.assertIn("mt_rows", payload["haplogroups"])
        self.assertIn("mt_notes", payload["haplogroups"])
        self.assertIn("status_notes", payload["haplogroups"])
        self.assertIn("assembly", payload["vcf"])
        self.assertIn(payload["vcf"]["assembly"], {"unknown", "GRCh37", "GRCh38"})
        self.assertIn("disease", payload)
        self.assertEqual(payload["disease"]["kind"], "disease")
        self.assertIn("additional", payload)
        self.assertEqual(payload["additional"]["kind"], "additional")
        self.assertIn("relatedness", payload)
        self.assertIn("coverage_status", payload["additional"]["estimates"][0])
        self.assertFalse(payload["relatedness"]["available"])
        self.assertIn("second VCF", payload["relatedness"]["notes"][0])

    def test_ancestry_aadr_labels_exist(self):
        from dna_compare.config import (
            ANCESTRY_5_AADR_POPS,
            ANCESTRY_5_RIGHT_POPS,
            ANCESTRY_AADR_POPS,
            ANCESTRY_RIGHT_POPS,
            default_settings,
        )
        from dna_compare.eigenstrat import AadrPanel

        panel = AadrPanel(default_settings())
        if not panel.available:
            self.skipTest("AADR HO panel not present")
        present = {rec.population for rec in panel.inds()}
        for label, pops in ANCESTRY_AADR_POPS.items():
            self.assertTrue(any(pop in present for pop in pops), f"{label} missing {pops}")
        for pop in ANCESTRY_RIGHT_POPS:
            self.assertIn(pop, present, f"outgroup {pop} missing")
        for label, pops in ANCESTRY_5_AADR_POPS.items():
            self.assertTrue(any(pop in present for pop in pops), f"5-source {label} missing {pops}")
        for pop in ANCESTRY_5_RIGHT_POPS:
            self.assertIn(pop, present, f"5-source outgroup {pop} missing")
        self.assertNotIn("Han", ANCESTRY_5_RIGHT_POPS)
        self.assertNotIn("French", ANCESTRY_5_RIGHT_POPS)

    def test_caste_alias_notes_include_pillai_and_missing_panels(self):
        notes = caste_alias_notes()
        joined = " ".join(notes)
        self.assertTrue(any("Pillai" in note for note in notes))
        self.assertIn("Vellalar", joined)
        self.assertIn("not a Tamil Nadu jati", joined)
        self.assertIn("Chettiyar", joined)
        self.assertIn("Vanniyar", joined)
        self.assertIn("Parayar", joined)
        self.assertIn("Cochin Jewish", joined)
        self.assertIn("not linked to Y-DNA J haplogroup", joined)
        self.assertNotIn("includes all", joined.lower())

    def test_multipart_form_parses_text_and_file_fields(self):
        from email.message import EmailMessage

        from dna_compare.web import parse_multipart

        body = (
            b"--xyz\r\n"
            b'Content-Disposition: form-data; name="sample"\r\n\r\n'
            b"demo.snps.vcf\r\n"
            b"--xyz\r\n"
            b'Content-Disposition: form-data; name="vcf"; filename="kit.vcf"\r\n'
            b"Content-Type: text/plain\r\n\r\n"
            b"##fileformat=VCFv4.2\n\r\n"
            b"--xyz--\r\n"
        )
        headers = EmailMessage()
        headers["Content-Type"] = "multipart/form-data; boundary=xyz"
        form = parse_multipart(headers, body)
        self.assertEqual(form.getfirst("sample"), "demo.snps.vcf")
        self.assertTrue("vcf" in form)
        self.assertEqual(form["vcf"].filename, "kit.vcf")
        self.assertEqual(form["vcf"].file.read(), b"##fileformat=VCFv4.2\n")

    def test_bootstrap_is_vendored_for_offline_ui(self):
        from dna_compare.report import STATIC_DIR, render_report_html
        from dna_compare.web import _static_file

        css = STATIC_DIR / "vendor" / "cerulean.min.css"
        js = STATIC_DIR / "vendor" / "bootstrap.bundle.min.js"
        self.assertTrue(css.is_file())
        self.assertTrue(js.is_file())
        head = css.read_bytes()[:120]
        self.assertIn(b"Bootswatch", head)
        self.assertIn(b"cerulean", head)
        page = (STATIC_DIR / "index.html").read_text(encoding="utf-8")
        self.assertIn("/static/vendor/cerulean.min.css", page)
        self.assertIn("/static/vendor/bootstrap.bundle.min.js", page)
        self.assertIn("other-select", page)
        self.assertIn("other-file", page)
        self.assertNotIn("cdn.", page.lower())
        self.assertIsNotNone(_static_file("/static/vendor/cerulean.min.css"))
        self.assertIsNone(_static_file("/static/../web.py"))
        report = render_report_html({"ok": True, "source_filename": "demo.snps.vcf", "vcf": {}})
        self.assertIn("cerulean", report)
        self.assertIn("navbar", report)
        self.assertIn("card", report)
        js = (STATIC_DIR / "app.js").read_text(encoding="utf-8")
        self.assertIn("communityRefCard", js)
        self.assertIn("relatednessTables", js)
        self.assertIn("snpListCard", js)
        self.assertIn("SNP list (all chromosomes)", js)
        self.assertIn("resultTabs", js)
        self.assertIn("diseaseCard", js)
        self.assertIn("Disease (research)", js)
        self.assertIn("disease-elevated", js)
        self.assertIn("<th>Coverage</th><th>Relative level</th><th>Percentile</th>", js)
        self.assertIn("additionalCard", js)
        self.assertIn("Additional details", js)
        self.assertIn("<th>Topic</th><th>Coverage</th><th>Finding</th><th>SNPs used / listed</th>", js)
        self.assertIn("<th>Evidence</th><th>Source</th><th>Note</th>", js)
        self.assertIn("cmp-disease", page)
        self.assertIn("nav-tabs", js)
        self.assertIn("tab-pane", js)
        self.assertIn("ancestry5Card", js)
        self.assertIn("qpAdm-style 5-source", js)
        self.assertIn("if (!block || !block.other_filename)", js)
        self.assertIn("hgPctClass", js)
        self.assertIn("table-warning", js)
        self.assertIn("table-success", js)
        self.assertIn("text-bg-warning", js)
        self.assertNotIn("haploShareTable", js)


if __name__ == "__main__":
    unittest.main()
