from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "data"
REFERENCE_DIR = DATA_DIR / "references"
HOMININ_DIR = REFERENCE_DIR / "hominin"
CASTE_DIR = REFERENCE_DIR / "caste"
TAMIL_COMMUNITY_REF = CASTE_DIR / "tamil_community_reference.tsv"
PUNJABI_COMMUNITY_REF = CASTE_DIR / "punjabi_community_reference.tsv"
BENGALI_COMMUNITY_REF = CASTE_DIR / "bengali_community_reference.tsv"
GUJARATI_COMMUNITY_REF = CASTE_DIR / "gujarati_community_reference.tsv"
MARATHI_COMMUNITY_REF = CASTE_DIR / "marathi_community_reference.tsv"
KERALA_COMMUNITY_REF = CASTE_DIR / "kerala_community_reference.tsv"
COMMUNITY_PANEL_ALIASES = CASTE_DIR / "community_panel_aliases.tsv"
DISEASE_DIR = REFERENCE_DIR / "disease"
PGS_REFERENCE_EUR = DISEASE_DIR / "pgs_reference_eur.tsv"
DISEASE_ELEVATED_PERCENTILE = 90.0
DISEASE_COVERAGE_READY_PCT = 20.0
DISEASE_COVERAGE_PERCENTILE_PCT = 50.0
ADDITIONAL_DIR = REFERENCE_DIR / "additional"
ADDITIONAL_DETAILS = ADDITIONAL_DIR / "details.tsv"
CACHE_DIR = REFERENCE_DIR / "cache"
SAMPLE_DIR = DATA_DIR / "samples"
UPLOAD_DIR = DATA_DIR / "uploads"

AADR_DIR = REFERENCE_DIR / "aadr"
AADR_STEM = AADR_DIR / "v66.p1_HO.aadr.patch.PUB"
AADR_GENO = Path(str(AADR_STEM) + ".geno")
AADR_IND = Path(str(AADR_STEM) + ".ind")
AADR_SNP = Path(str(AADR_STEM) + ".snp")
AADR_ANNO = AADR_DIR / "v66.p1_HO.aadr.PUB.anno"

LIFTOVER_DIR = REFERENCE_DIR / "liftover"
HG38_TO_HG19_CHAIN = LIFTOVER_DIR / "hg38ToHg19.over.chain.gz"
HG38_TO_HG19_CHAIN_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"

# Legacy locations (pre-move). default_settings() prefers these if still in the project root.
_LEGACY_STEM = PROJECT_ROOT / "v66.p1_HO.aadr.patch.PUB"

# High-coverage archaic VCFs (hg19/GRCh37 to match HO coordinates).
# Download later into data/references/hominin/ using these exact names.
HOMININ_DOWNLOADS = {
    "altai_neanderthal": {
        "filename": "AltaiNea.hg19_1000g.vcf.gz",
        "index": "AltaiNea.hg19_1000g.vcf.gz.tbi",
        "source": "http://cdna.eva.mpg.de/neandertal/altai/AltaiNeandertal/VCF/",
        "label": "Altai Neanderthal (high-coverage VCF)",
    },
    "vindija_neanderthal": {
        "filename": "Vindija33.19.hg19_1000g.vcf.gz",
        "index": "Vindija33.19.hg19_1000g.vcf.gz.tbi",
        "source": "http://cdna.eva.mpg.de/neandertal/Vindija/VCF/",
        "label": "Vindija 33.19 Neanderthal (high-coverage VCF)",
    },
    "chagyrskaya_neanderthal": {
        "filename": "Chagyrskaya-Phalanx.hg19.vcf.gz",
        "index": "Chagyrskaya-Phalanx.hg19.vcf.gz.tbi",
        "source": "http://ftp.eva.mpg.de/neandertal/Chagyrskaya/VCF/",
        "label": "Chagyrskaya Neanderthal VCF",
    },
    "denisovan": {
        "filename": "DenisovaPinky.hg19_1000g.vcf.gz",
        "index": "DenisovaPinky.hg19_1000g.vcf.gz.tbi",
        "source": "http://cdna.eva.mpg.de/denisova/",
        "label": "Denisova 3 (Pinky) high-coverage VCF",
    },
    "neanderthal_informative_snps": {
        "filename": "neanderthal_informative_snps.tsv",
        "source": "Sankararaman / Prüfer archaic-informative sites (build yourself or from paper supplements)",
        "label": "Optional SNP list: chrom pos ancestral derived source=Neanderthal",
        "columns": "chrom\tpos\tref\talt\tarchaic_allele\tsource",
    },
    "denisovan_informative_snps": {
        "filename": "denisovan_informative_snps.tsv",
        "source": "Browning Sprime / Sankararaman Denisovan-informative sites",
        "label": "Optional SNP list: chrom pos ancestral derived source=Denisovan",
        "columns": "chrom\tpos\tref\talt\tarchaic_allele\tsource",
    },
}

# Extra Indian-caste frequency / genotype files (not in AADR HO, or denser).
# Place into data/references/caste/ using these names.
CASTE_DOWNLOADS = {
    "caste_allele_frequencies": {
        "filename": "indian_caste_allele_frequencies.tsv",
        "label": "Preferred compact table you can build from any genotype dataset",
        "columns": "chrom\tpos\trsid\tref\talt\tpopulation\talt_freq\tn_haplotypes",
        "source": "Derived from 1000G / GenomeAsia / published Indian caste arrays",
    },
    "1000g_gih": {
        "filename": "1000G_GIH.snps.vcf.gz",
        "index": "1000G_GIH.snps.vcf.gz.tbi",
        "label": "1000 Genomes Gujarati Indian in Houston (GIH)",
        "source": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/",
    },
    "1000g_pjl": {
        "filename": "1000G_PJL.snps.vcf.gz",
        "index": "1000G_PJL.snps.vcf.gz.tbi",
        "label": "1000 Genomes Punjabi in Lahore (PJL)",
        "source": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/",
    },
    "1000g_beb": {
        "filename": "1000G_BEB.snps.vcf.gz",
        "index": "1000G_BEB.snps.vcf.gz.tbi",
        "label": "1000 Genomes Bengali in Bangladesh (BEB)",
        "source": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/",
    },
    "1000g_stu": {
        "filename": "1000G_STU.snps.vcf.gz",
        "index": "1000G_STU.snps.vcf.gz.tbi",
        "label": "1000 Genomes Sri Lankan Tamil in the UK (STU)",
        "source": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/",
    },
    "1000g_itu": {
        "filename": "1000G_ITU.snps.vcf.gz",
        "index": "1000G_ITU.snps.vcf.gz.tbi",
        "label": "1000 Genomes Indian Telugu in the UK (ITU)",
        "source": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/",
    },
    "genomeasia_india": {
        "filename": "GenomeAsia100K_India.snps.vcf.gz",
        "index": "GenomeAsia100K_India.snps.vcf.gz.tbi",
        "label": "GenomeAsia 100K Indian-caste subset (request access)",
        "source": "https://browser.genomeasia100k.org/",
    },
    "estonian_biocentre_india": {
        "filename": "EstonianBiocentre_IndianCastes.eigenstrat.geno",
        "companion": [
            "EstonianBiocentre_IndianCastes.eigenstrat.snp",
            "EstonianBiocentre_IndianCastes.eigenstrat.ind",
        ],
        "label": "Estonian Biocentre / published Indian caste Human Origins-style genotypes",
        "source": "Metspalu / Reich / Moorjani Indian-caste supplements",
    },
    "nakatsuka2017_india": {
        "filename": "nakatsuka2017_india.snp",
        "companion": ["nakatsuka2017_india.ind", "nakatsuka2017_india.geno"],
        "label": "Nakatsuka et al. 2017 Cell India genotype pack (many jati/caste labels)",
        "source": "https://doi.org/10.1016/j.cell.2017.09.019",
    },
}

# AADR group IDs already on disk in v66.p1_HO (used immediately).
ARCHAIC_AADR_SAMPLES = {
    "Altai_Neanderthal": "AltaiNeanderthal.DG",
    "Vindija_Neanderthal": "Vindija.DG",
    "Chagyrskaya_Neanderthal": "Chagyrskaya8.DG",
    "Denisova": "Denisova3.DG",
    "Chimp": "Chimp.REF",
}

OUTGROUP_AADR_POPS = ("Mbuti", "Yoruba")

# Right / outgroup pops for qpAdm-style ancestry (must not overlap sources).
ANCESTRY_RIGHT_POPS: tuple[str, ...] = (
    "Mbuti",
    "Yoruba",
    "Ju_hoan_North",
    "Mandenka",
    "French",
    "Han",
    "Papuan",
    "Karitiana",
    "Ulchi",
)

# Map display caste/community → AADR .ind population labels.
CASTE_AADR_POPS: dict[str, tuple[str, ...]] = {
    "Brahmin": ("Brahmin",),
    "Yadava": ("Yadava",),
    "Kapu": ("Kapu",),
    "Mala": ("Mala",),
    "Madiga": ("Madiga",),
    "Irula": ("Irula",),
    "Relli": ("Relli_1", "Relli_2"),
    "Gujarati": ("GujaratiA", "GujaratiB", "GujaratiC", "GujaratiD", "GIH"),
    "Punjabi": ("Punjabi", "PJL"),
    "Bengali": ("BEB",),
    "Telugu": ("ITU",),
    "Tamil": ("STU", "STU-1", "STU-2"),
    "Vellalar": ("VLR",),
    "Cochin_Jew": ("Jew_Cochin",),
}

# Generic first; per-label notes are coverage caveats, not a preferred caste list.
CASTE_GENERAL_NOTES: tuple[str, ...] = (
    "Bars are the AADR HO community labels that exist in this file — not a ranked or complete caste list.",
    "A title or subcaste may sit under a scored label. Example: Pillai is often a Vellalar title; "
    "the HO VLR set is 9 samples, not every Vellalar subdivision.",
)

# Short UI notes: titles/subcastes often associated with a scored label.
CASTE_LABEL_NOTES: dict[str, str] = {
    "Tamil": (
        "Tamil (STU): Sri Lankan Tamil in 1000 Genomes — a language/region label, not a Tamil Nadu jati."
    ),
    "Brahmin": (
        "Brahmin: two HO samples only. Iyer / Iyengar are not this file."
    ),
    "Telugu": (
        "Telugu (ITU): Indian Telugu in the UK (1000 Genomes), not a single caste."
    ),
    "Vellalar": (
        "Vellalar (VLR): 9 Mondal samples. Not Gounder, Mudaliar, or every Pillai-using family."
    ),
    "Kapu": (
        "Kapu: Andhra community in AADR. Not Chettiyar / Nagarathar."
    ),
    "Mala": (
        "Mala: Andhra group in AADR. Not Pallar or Parayar."
    ),
    "Madiga": (
        "Madiga: Andhra group in AADR. Not Pallar or Parayar."
    ),
    "Irula": (
        "Irula: AADR tribal/community samples from South India (also GenomeAsia IRU, restricted)."
    ),
    "Gujarati": "Gujarati: 1000 Genomes GIH plus HO Gujarati A–D, not a single jati.",
    "Punjabi": "Punjabi: 1000 Genomes PJL plus HO Punjabi.",
    "Bengali": "Bengali: 1000 Genomes BEB.",
    "Cochin_Jew": (
        "Cochin Jewish (Jew_Cochin): five AADR Kerala Jewish reference samples. "
        "The bar is a small SNP mixture weight — not Jewish identity, not Middle Eastern ancestry, "
        "and not linked to Y-DNA J haplogroup (South Asian J2/J2b is a separate lineage from Jewish J)."
    ),
}

# Mentioned often, but no public HO bar.
CASTE_MISSING_PANEL_NOTES: tuple[str, ...] = (
    "No public HO panel for Chettiyar / Nattukottai Chettiar, Vanniyar, or Parayar. "
    "Do not read a missing name from a nearby bar (Kapu is not Chettiyar; Mala/Madiga are not Parayar).",
    "No public HO panel for Maratha, Kunbi, Dhangar, Mahar, or Maharashtra Brahmin subdivisions. "
    "HO Brahmin is two generic samples — not Deshastha / Chitpavan / Karhade. "
    "Use Community reference ranges when the filename suggests Maharashtra or a listed jati.",
    "No public HO panel for Nair, Ezhava, Namboothiri, Syrian Christian, or Pulaya. "
    "Cochin Jewish (Jew_Cochin) is a small Kerala Jewish reference set — not Malayali jati identity. "
    "Use Community reference ranges when the filename suggests Kerala / Malayalam or a listed community.",
)

# Named packs of AADR HO groups already in v66.p1 (no extra download).
POPULATION_PACKS: dict[str, dict[str, tuple[str, ...]]] = {
    "greek": {
        "Greek_modern": ("Greek", "Greek_WGA"),
        "Greece_Crete_EMBA": ("Greece_Crete_HgCharalambos_EMBA",),
        "Greece_Crete_LBA": ("Greece_Crete_LBA",),
        "Greece_Peloponnese_LBA": ("Greece_Peloponnese_LBA",),
        "Greece_PalaceofNestor_BA": ("Greece_PalaceofNestor_BA",),
        "Greece_Kastrouli_BA": ("Greece_Kastrouli_BA",),
        "Greece_LBA": ("Greece_LBA",),
        "Greece_Crete_N": ("Greece_Crete_Aposelemis_N",),
        "Cycladic_EBA": ("Greece_EpanoKoufonisi_EBA_Cycladic",),
    },
    "chinese": {
        "Han": ("Han",),
        "Dai": ("Dai",),
        "Naxi": ("Naxi",),
        "China_Yangshao_LN": ("China_Baligang_LN_Yangshao",),
        "China_Xiaoheyan_LN": ("China_LN_Xiaoheyan",),
        "China_Shang": ("China_Jinan_LiuJiaZhuang_Shang", "China_Henan_Xisima_LShang"),
        "China_Zhou": ("China_Weifang_XinZhi_Zhou", "China_Baligang_BA_EasternZhou"),
        "China_Tibet_IA": ("China_Tibet_Gebusailu_IA",),
    },
    "persian": {
        "Iranian": ("Iranian",),
        "Iranian_Zoroastrian": ("Iranian_Zoroastrian",),
        "Iran_Hasanlu_IA": ("Iran_Hasanlu_IA",),
        "Iran_GanjDareh_N": ("Iran_GanjDareh_N",),
        "Iran_TepeHissar_C": ("Iran_TepeHissar_C",),
        "Iran_Parthian": ("Iran_LiarsangBon_Parthian",),
        "Turkey_Ancient_Persian": ("Turkey_Ancient_Persian",),
    },
    "caste": CASTE_AADR_POPS,
}

# 3-source South Asia model (Narasimhan-style left pops). One Steppe number.
ANCESTRY_AADR_POPS: dict[str, tuple[str, ...]] = {
    "Steppe_MLBA": ("Russia_Chelyabinsk_MLBA_Sintashta",),
    "Indus_Periphery": ("Iran_ShahriSokhta_BA1-1", "Iran_ShahriSokhta_BA2-2"),
    "AASI_Onge": ("ONG",),
}

# 5-source: same three plus Anatolian farmer and East Asian. Han/French stay off the right set.
ANCESTRY_5_AADR_POPS: dict[str, tuple[str, ...]] = {
    "AASI_Onge": ("ONG",),
    "Indus_Periphery": ("Iran_ShahriSokhta_BA1-1", "Iran_ShahriSokhta_BA2-2"),
    "Steppe_MLBA": ("Russia_Chelyabinsk_MLBA_Sintashta",),
    "Anatolia_N": ("Turkey_N",),
    "East_Asian": ("Dai",),
}

ANCESTRY_5_RIGHT_POPS: tuple[str, ...] = (
    "Mbuti",
    "Yoruba",
    "Ju_hoan_North",
    "Mandenka",
    "Papuan",
    "Karitiana",
    "Ulchi",
)

DEFAULT_POPULATION_PACKS = ("greek", "chinese", "persian", "caste")

# PGS Catalog scoring files (https://www.pgscatalog.org/). Not a diagnosis.
PGS_CATALOG_SCORES: tuple[dict[str, object], ...] = (
    {
        "pgs_id": "PGS000133",
        "trait": "Schizophrenia",
        "filename": "PGS000133.txt.gz",
        "url": "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS000133/ScoringFiles/PGS000133.txt.gz",
        "citation": "Zheutlin AB et al. Am J Psychiatry (2019)",
        "n_variants": 604645,
        "auto_download": True,
    },
    {
        "pgs_id": "PGS000145",
        "trait": "Major depressive disorder",
        "filename": "PGS000145.txt.gz",
        "url": "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS000145/ScoringFiles/PGS000145.txt.gz",
        "citation": "Cai N et al. Nat Genet (2020)",
        "n_variants": 21510,
        "auto_download": False,
    },
    {
        "pgs_id": "PGS002786",
        "trait": "Bipolar disorder",
        "filename": "PGS002786.txt.gz",
        "url": "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS002786/ScoringFiles/PGS002786.txt.gz",
        "citation": "Gui Y et al. Transl Psychiatry (2022)",
        "n_variants": 948996,
        "auto_download": True,
    },
    {
        "pgs_id": "PGS002746",
        "trait": "ADHD",
        "filename": "PGS002746.txt.gz",
        "url": "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS002746/ScoringFiles/PGS002746.txt.gz",
        "citation": "Lahey BB et al. J Psychiatr Res (2022)",
        "n_variants": 513659,
        "auto_download": True,
    },
    {
        "pgs_id": "PGS004451",
        "trait": "Anxiety disorders",
        "filename": "PGS004451.txt.gz",
        "url": "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS004451/ScoringFiles/PGS004451.txt.gz",
        "citation": "Jung H et al. Commun Biol (2024)",
        "n_variants": 1059939,
        "auto_download": True,
    },
    {
        "pgs_id": "PGS004427",
        "trait": "Fluid intelligence",
        "filename": "PGS004427.txt.gz",
        "url": "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS004427/ScoringFiles/PGS004427.txt.gz",
        "citation": "Jung H et al. Commun Biol (2024)",
        "n_variants": 1059939,
        "auto_download": True,
        "trait_category": "cognitive",
    },
    {
        "pgs_id": "PGS002012",
        "trait": "Educational attainment",
        "filename": "PGS002012.txt.gz",
        "url": "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS002012/ScoringFiles/PGS002012.txt.gz",
        "citation": "Privé F et al. Am J Hum Genet (2022)",
        "n_variants": 50413,
        "auto_download": True,
        "trait_category": "cognitive",
    },
    {
        "pgs_id": "PGS000011",
        "trait": "Coronary artery disease",
        "filename": "PGS000011.txt.gz",
        "url": "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS000011/ScoringFiles/PGS000011.txt.gz",
        "citation": "Tada H et al. Eur Heart J (2015) GRS50",
        "n_variants": 50,
        "auto_download": False,
    },
    {
        "pgs_id": "PGS004226",
        "trait": "Type 2 diabetes",
        "filename": "PGS004226.txt.gz",
        "url": "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS004226/ScoringFiles/PGS004226.txt.gz",
        "citation": "Liu J et al. Nutrients (2023) PRS50_T2DEur",
        "n_variants": 50,
        "auto_download": False,
    },
    {
        "pgs_id": "PGS000038",
        "trait": "Stroke",
        "filename": "PGS000038.txt.gz",
        "url": "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS000038/ScoringFiles/PGS000038.txt.gz",
        "citation": "Rutten-Jacobs LC et al. BMJ (2018) PRS90",
        "n_variants": 90,
        "auto_download": False,
    },
    {
        "pgs_id": "PGS000025",
        "trait": "Alzheimer's disease",
        "filename": "PGS000025.txt.gz",
        "url": "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS000025/ScoringFiles/PGS000025.txt.gz",
        "citation": "Chouraki V et al. J Alzheimers Dis (2016)",
        "n_variants": 19,
        "auto_download": False,
    },
    {
        "pgs_id": "PGS004252",
        "trait": "Asthma",
        "filename": "PGS004252.txt.gz",
        "url": "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS004252/ScoringFiles/PGS004252.txt.gz",
        "citation": "Zhu Y et al. Ecotoxicol Environ Saf (2023) PRS212",
        "n_variants": 212,
        "auto_download": False,
    },
    {
        "pgs_id": "PGS000004",
        "trait": "Breast cancer",
        "filename": "PGS000004.txt.gz",
        "url": "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS000004/ScoringFiles/PGS000004.txt.gz",
        "citation": "Mavaddat N et al. Am J Hum Genet (2018) PRS313",
        "n_variants": 313,
        "auto_download": False,
    },
    {
        "pgs_id": "PGS000030",
        "trait": "Prostate cancer",
        "filename": "PGS000030.txt.gz",
        "url": "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS000030/ScoringFiles/PGS000030.txt.gz",
        "citation": "Schumacher FR et al. Nat Genet (2018)",
        "n_variants": 147,
        "auto_download": False,
    },
    {
        "pgs_id": "PGS000765",
        "trait": "Colorectal cancer",
        "filename": "PGS000765.txt.gz",
        "url": "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS000765/ScoringFiles/PGS000765.txt.gz",
        "citation": "Huyghe JR et al. Nat Genet (2018) PRS_CRC95",
        "n_variants": 95,
        "auto_download": False,
    },
)

MAX_SAMPLES_PER_POP = 24
VARIANT_PREVIEW_LIMIT = 200
VARIANT_EMBED_LIMIT = 5000
SNP_PAGE_SIZE = 200

# Chromosome codes used in EIGENSTRAT/HO .snp files.
CHROM_TO_CODE = {str(i): i for i in range(1, 23)}
CHROM_TO_CODE.update(
    {
        "X": 23,
        "Y": 24,
        "XY": 25,
        "MT": 26,
        "M": 26,
    }
)
CODE_TO_CHROM = {v: k for k, v in CHROM_TO_CODE.items() if k not in {"M", "XY"}}
CODE_TO_CHROM[25] = "XY"
CODE_TO_CHROM[26] = "MT"


@dataclass
class Settings:
    aadr_geno: Path = AADR_GENO
    aadr_ind: Path = AADR_IND
    aadr_snp: Path = AADR_SNP
    aadr_anno: Path = AADR_ANNO
    hominin_dir: Path = HOMININ_DIR
    caste_dir: Path = CASTE_DIR
    tamil_community_ref: Path = TAMIL_COMMUNITY_REF
    punjabi_community_ref: Path = PUNJABI_COMMUNITY_REF
    bengali_community_ref: Path = BENGALI_COMMUNITY_REF
    gujarati_community_ref: Path = GUJARATI_COMMUNITY_REF
    marathi_community_ref: Path = MARATHI_COMMUNITY_REF
    kerala_community_ref: Path = KERALA_COMMUNITY_REF
    community_panel_aliases: Path = COMMUNITY_PANEL_ALIASES
    disease_dir: Path = DISEASE_DIR
    additional_details: Path = ADDITIONAL_DETAILS
    cache_dir: Path = CACHE_DIR
    variant_preview_limit: int = VARIANT_PREVIEW_LIMIT
    variant_embed_limit: int = VARIANT_EMBED_LIMIT
    max_samples_per_pop: int = MAX_SAMPLES_PER_POP
    compare_hominin: bool = True
    compare_caste: bool = True
    compare_populations: bool = True
    compare_ancestry: bool = True
    compare_haplogroups: bool = True
    compare_disease: bool = True
    analyze_timings: bool = False
    auto_download_pgs: bool = False
    assume_assembly: str | None = None
    auto_download_chain: bool = True
    liftover_chain: Path = HG38_TO_HG19_CHAIN
    population_packs: tuple[str, ...] = DEFAULT_POPULATION_PACKS
    extra_pops: dict[str, tuple[str, ...]] = field(default_factory=dict)
    extra_caste_groups: dict[str, tuple[str, ...]] = field(default_factory=dict)
    ancestry_right_pops: tuple[str, ...] = ANCESTRY_RIGHT_POPS
    ancestry_5_right_pops: tuple[str, ...] = ANCESTRY_5_RIGHT_POPS

    def caste_groups(self) -> dict[str, tuple[str, ...]]:
        merged = dict(CASTE_AADR_POPS)
        merged.update(self.extra_caste_groups)
        return merged

    def ancestry_groups(self) -> dict[str, tuple[str, ...]]:
        return dict(ANCESTRY_AADR_POPS)

    def ancestry_5_groups(self) -> dict[str, tuple[str, ...]]:
        return dict(ANCESTRY_5_AADR_POPS)

    def selected_population_groups(self) -> dict[str, tuple[str, ...]]:
        groups: dict[str, tuple[str, ...]] = {}
        for pack in self.population_packs:
            groups.update(POPULATION_PACKS.get(pack, {}))
        groups.update(self.extra_pops)
        groups.update(self.extra_caste_groups)
        return groups


def _existing_aadr_stem() -> Path:
    for stem in (AADR_STEM, _LEGACY_STEM):
        if Path(str(stem) + ".geno").exists() and Path(str(stem) + ".ind").exists():
            return stem
    return AADR_STEM


def default_settings() -> Settings:
    stem = _existing_aadr_stem()
    return Settings(
        aadr_geno=Path(str(stem) + ".geno"),
        aadr_ind=Path(str(stem) + ".ind"),
        aadr_snp=Path(str(stem) + ".snp"),
        aadr_anno=AADR_ANNO,
    )
