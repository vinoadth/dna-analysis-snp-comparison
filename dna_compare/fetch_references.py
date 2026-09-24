from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np

from dna_compare.config import (
    ARCHAIC_AADR_SAMPLES,
    HG38_TO_HG19_CHAIN,
    HOMININ_DIR,
    HOMININ_DOWNLOADS,
    Settings,
    default_settings,
)
from dna_compare.liftover import ensure_hg38_to_hg19_chain
from dna_compare.eigenstrat import AadrPanel, PACKED_MISSING
from dna_compare.vcf_parser import normalize_chrom

GT_MAP = {0: "0/0", 1: "0/1", 2: "1/1", 3: "./."}

HOMININ_VCF_SAMPLES = {
    "AltaiNea.hg19_1000g.vcf.gz": "AltaiNeanderthal.DG",
    "Vindija33.19.hg19_1000g.vcf.gz": "Vindija.DG",
    "Chagyrskaya-Phalanx.hg19.vcf.gz": "Chagyrskaya8.DG",
    "DenisovaPinky.hg19_1000g.vcf.gz": "Denisova3.DG",
}


def _write_vcf(path: Path, sample_ids: list[str], genotypes: np.ndarray, snps) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    header = [
        "##fileformat=VCFv4.2",
        "##source=AADR v66.p1_HO TGENO subset (full remote VCFs are 50-70GB/genome; disk-limited extract)",
        "##INFO=<ID=HO,Number=0,Type=Flag,Description=\"Human Origins SNP\">",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(sample_ids),
    ]
    with gzip.open(path, "wt") as handle:
        handle.write("\n".join(header) + "\n")
        n_samples = genotypes.shape[0]
        for snp in snps:
            gts = []
            any_called = False
            for i in range(n_samples):
                code = int(genotypes[i, snp.index])
                if code != PACKED_MISSING:
                    any_called = True
                gts.append(GT_MAP.get(code, "./."))
            if not any_called:
                continue
            handle.write(
                f"{snp.chrom}\t{snp.pos}\t{snp.rsid}\t{snp.allele2}\t{snp.allele1}\t"
                f".\tPASS\tHO\tGT\t" + "\t".join(gts) + "\n"
            )


def _convert_sprime_bed(bed_path: Path, out_path: Path) -> int:
    n = 0
    with bed_path.open() as src, out_path.open("w") as dst:
        dst.write("chrom\tpos\tref\talt\tarchaic_allele\tsource\n")
        for line in src:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6:
                continue
            chrom = normalize_chrom(cols[0])
            pos = int(cols[2]) if cols[2].isdigit() else int(cols[1]) + 1
            ref, alt = cols[4].upper(), cols[5].upper()
            if len(ref) != 1 or len(alt) != 1:
                continue
            dst.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{alt}\tBrowning2018_Sprime_NeanderthalMatch\n")
            n += 1
    return n


def _write_denisovan_informative(panel: AadrPanel, out_path: Path) -> int:
    snps = panel.snps()
    den = panel.dosage_allele1(ARCHAIC_AADR_SAMPLES["Denisova"]) / 2.0
    altai = panel.dosage_allele1(ARCHAIC_AADR_SAMPLES["Altai_Neanderthal"]) / 2.0
    mbuti = panel.population_allele1_freq("Mbuti")
    n = 0
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as dst:
        dst.write("chrom\tpos\tref\talt\tarchaic_allele\tsource\n")
        for snp in snps:
            d, a, m = den[snp.index], altai[snp.index], mbuti[snp.index]
            if not (np.isfinite(d) and np.isfinite(a) and np.isfinite(m)):
                continue
            if d >= 0.9 and m <= 0.1 and a <= 0.2:
                dst.write(
                    f"{snp.chrom}\t{snp.pos}\t{snp.allele2}\t{snp.allele1}\t"
                    f"{snp.allele1}\tAADR_HO_Denisova3_vs_Mbuti\n"
                )
                n += 1
    return n


def fetch_references(settings: Settings | None = None) -> list[str]:
    settings = settings or default_settings()
    panel = AadrPanel(settings)
    if not panel.available:
        raise FileNotFoundError("AADR HO files are required to build compact reference extracts.")
    notes = []
    HOMININ_DIR.mkdir(parents=True, exist_ok=True)
    snps = panel.snps()

    sprime = HOMININ_DIR / "ALL1KG_NMATCH_sprime_results.bed"
    neo_tsv = HOMININ_DIR / HOMININ_DOWNLOADS["neanderthal_informative_snps"]["filename"]
    if sprime.exists():
        n = _convert_sprime_bed(sprime, neo_tsv)
        notes.append(f"Wrote {neo_tsv.name} ({n} SNPs) from Browning/Sprime BED")
    else:
        notes.append(f"Missing {sprime.name}; skip neanderthal informative TSV from Sprime")

    den_tsv = HOMININ_DIR / HOMININ_DOWNLOADS["denisovan_informative_snps"]["filename"]
    n = _write_denisovan_informative(panel, den_tsv)
    notes.append(f"Wrote {den_tsv.name} ({n} SNPs) from AADR Denisova3 vs Mbuti")

    for filename, sample_id in HOMININ_VCF_SAMPLES.items():
        path = HOMININ_DIR / filename
        rec = panel.sample_by_id(sample_id)
        if rec is None:
            notes.append(f"Missing AADR sample {sample_id}; skip {filename}")
            continue
        geno = panel.geno().read_individual(rec.index)
        _write_vcf(path, [sample_id], geno[None, :], snps)
        notes.append(f"Wrote {filename} from AADR {sample_id}")

    notes.append("Greek/Chinese/Persian/caste comparisons use AADR HO directly; extra caste VCFs are not written.")
    chain = ensure_hg38_to_hg19_chain(HG38_TO_HG19_CHAIN, download=True)
    if chain is not None:
        notes.append(f"Ready {chain} for GSA/gtc2vcf GRCh38 → GRCh37 liftover")
    else:
        notes.append(f"Could not download {HG38_TO_HG19_CHAIN.name}; place the UCSC hg38ToHg19 chain there.")
    return notes
