import tempfile
import unittest
from io import StringIO
from pathlib import Path

from dna_compare.assembly import detect_assembly
from dna_compare.config import Settings, default_settings
from dna_compare.liftover import lift_query_index, load_chain
from dna_compare.service import AnalysisService
from dna_compare.vcf_parser import parse_vcf

GSA_HEADER_VCF = """##fileformat=VCFv4.2
##source=bcftools_gtc2vcf
##reference=GSA-24v3-0_A1.bpm
##contig=<ID=1,length=248956422>
##contig=<ID=Y,length=57227415>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TEST1
Y	12914512	rs2032624	A	C	.	PASS	.	GT	1
Y	13470103	rs2032658	A	G	.	PASS	.	GT	1
1	100	rs1	A	G	.	PASS	.	GT	0/0
"""

HG19_CONTIG_VCF = """##fileformat=VCFv4.2
##reference=GRCh37
##contig=<ID=1,length=249250621>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TEST1
1	752566	rs3094315	G	A	.	PASS	.	GT	0/1
"""

TINY_CHAIN = """chain 1000 chrY 57227415 + 12914511 12914512 chrY 59373566 + 15026423 15026424 1
1
"""

GTC2VCF_GRCH37_VCF = """##fileformat=VCFv4.2
##BPM=GSA-24v3-0_A1.bpm
##EGT=GSA-24v3-0_A1_ClusterFile.egt
##bcftools_gtc2vcfVersion=1.21+htslib-1.21
##bcftools_gtc2vcfCommand=gtc2vcf -b GSA-24v3-0_A1.bpm -f /refs/genome/grch37/human_g1k_v37.fasta
##contig=<ID=1,length=249250621>
##contig=<ID=Y,length=59373566>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TEST1
Y	15026424	rs2032624	A	C	.	PASS	.	GT	0
Y	21733165	rs3908	D	I	.	PASS	.	GT	1
"""


class AssemblyLiftoverTests(unittest.TestCase):
    def test_gsa_gtc2vcf_header_is_grch38(self):
        summary, _index = parse_vcf(StringIO(GSA_HEADER_VCF))
        detect_assembly(summary)
        self.assertEqual(summary.assembly, "GRCh38")
        self.assertIn("248956422", summary.assembly_evidence)
        self.assertEqual(summary.reference, "GSA-24v3-0_A1.bpm")
        self.assertIn("gtc2vcf", summary.source or "")

    def test_hg19_contig_is_grch37(self):
        summary, _index = parse_vcf(StringIO(HG19_CONTIG_VCF))
        detect_assembly(summary)
        self.assertEqual(summary.assembly, "GRCh37")

    def test_assume_assembly_overrides_header(self):
        summary, _index = parse_vcf(StringIO(GSA_HEADER_VCF))
        detect_assembly(summary, assume="GRCh37")
        self.assertEqual(summary.assembly, "GRCh37")
        self.assertIn("explicit", summary.assembly_evidence)

    def test_gtc2vcf_with_grch37_fasta_is_not_lifted(self):
        summary, _index = parse_vcf(StringIO(GTC2VCF_GRCH37_VCF))
        detect_assembly(summary)
        self.assertEqual(summary.assembly, "GRCh37")
        self.assertIn("gtc2vcf", (summary.source or "").lower())
        settings = default_settings()
        settings.auto_download_chain = False
        settings.compare_hominin = False
        settings.compare_caste = False
        settings.compare_populations = False
        settings.compare_ancestry = False
        result = AnalysisService(settings).analyze(StringIO(GTC2VCF_GRCH37_VCF), filename="c8xy.vcf")
        self.assertIsNone(result.vcf.lifted_to)
        calls = {c.marker: c for c in result.haplogroups.markers}
        self.assertEqual(calls["M173"].status, "ancestral")
        self.assertEqual(calls["M17"].status, "conflict")

    def test_tiny_chain_lifts_m173(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "tiny.over.chain"
            path.write_text(TINY_CHAIN, encoding="utf-8")
            chain = load_chain(path)
            self.assertEqual(chain.convert("Y", 12914512), ("Y", 15026424, "+"))
            _summary, index = parse_vcf(StringIO(GSA_HEADER_VCF))
            lifted, stats = lift_query_index(index, chain)
            self.assertIn(("Y", 15026424), lifted)
            self.assertNotIn(("Y", 12914512), lifted)
            self.assertGreaterEqual(stats["n_lifted"], 1)
            self.assertEqual(lifted[("Y", 15026424)]["rsid"], "rs2032624")

    def test_analyze_lifts_when_chain_present(self):
        with tempfile.TemporaryDirectory() as tmp:
            chain_path = Path(tmp) / "tiny.over.chain"
            chain_path.write_text(TINY_CHAIN, encoding="utf-8")
            settings = default_settings()
            settings.liftover_chain = chain_path
            settings.auto_download_chain = False
            settings.compare_hominin = False
            settings.compare_caste = False
            settings.compare_populations = False
            settings.compare_ancestry = False
            result = AnalysisService(settings).analyze(StringIO(GSA_HEADER_VCF), filename="gsa.vcf")
            self.assertTrue(result.ok)
            self.assertEqual(result.vcf.assembly, "GRCh38")
            self.assertEqual(result.vcf.lifted_to, "GRCh37")
            self.assertGreaterEqual(result.vcf.n_lifted, 1)
            self.assertTrue(any("lifted" in note.lower() for note in result.haplogroups.notes))
            m173 = next(call for call in result.haplogroups.markers if call.marker == "M173")
            self.assertEqual(m173.status, "derived")
            self.assertEqual(m173.pos, 15026424)

    def test_analyze_without_chain_still_scores_hg38_markers(self):
        settings = Settings(
            auto_download_chain=False,
            liftover_chain=Path("/tmp/missing-hg38ToHg19.over.chain.gz"),
            compare_hominin=False,
            compare_caste=False,
            compare_populations=False,
            compare_ancestry=False,
        )
        result = AnalysisService(settings).analyze(StringIO(GSA_HEADER_VCF), filename="gsa.vcf")
        self.assertTrue(result.ok)
        self.assertEqual(result.vcf.assembly, "GRCh38")
        self.assertIsNone(result.vcf.lifted_to)
        self.assertTrue(any("GRCh38" in note for note in result.haplogroups.notes))
        m173 = next(call for call in result.haplogroups.markers if call.marker == "M173")
        self.assertEqual(m173.status, "derived")


if __name__ == "__main__":
    unittest.main()
