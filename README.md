# DNA Compare

Compare your SNP data against the [Allen Ancient DNA Resource (AADR)](https://reich.hms.harvard.edu/allen-ancient-dna-resource-aadr-downloads) v66.1 **Human Origins** panel. Upload a VCF or 23andMe-style raw export and get mixture-style population weights, South Asian community reference bars, ancestry models, Y/mt haplogroup context, archaic hominin overlap, polygenic score coverage, and optional pairwise relatedness.

The web dashboard and CLI share the same analysis engine. Mixture percentages are computed from overlapping Human Origins SNPs on **GRCh37 / hg19** — they describe genetic similarity to reference groups, not ethnic identity or caste assignment.

## Features

- **Population packs** — Greek, Chinese, Persian/Iranian, and Indian community labels from AADR HO
- **Ancestry models** — 3-source and 5-source South Asia qpAdm-style estimates
- **Haplogroups** — Y and mtDNA backbone markers scored from your file, with AADR group frequency context (requires the `.anno` file)
- **Hominin** — Neanderthal/Denisova informative sites (compact extracts ship with the repo; full VCFs optional)
- **Disease (PGS)** — overlap with selected polygenic scores from the PGS Catalog
- **Relatedness** — optional second sample for kinship / IBS estimates
- **Input formats** — `.vcf`, `.vcf.gz`, 23andMe / GEDmatch-style raw `.txt`, `.23andme`
- **Assembly handling** — GRCh38 exports (e.g. Illumina GSA v3 / `gtc2vcf`) are auto-lifted to hg19 when a UCSC chain file is available

Bundled demo files under `data/samples/` (`demo.snps.vcf`, `demo.raw23.txt`, …) let you try the UI before adding your own data.

## Required reference data (AADR)

The AADR Human Origins Eigenstrat panel is **not** included in this repository (it is several gigabytes). Download it from Harvard Dataverse:

**https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/FFIDCW**

From that dataset, place these four files in **`data/references/aadr/`** using **exact** filenames:

| File | Role |
| --- | --- |
| `v66.p1_HO.aadr.patch.PUB.geno` | Packed genotype matrix (~4 GB) |
| `v66.p1_HO.aadr.patch.PUB.ind` | Sample metadata (population labels) |
| `v66.p1_HO.aadr.patch.PUB.snp` | SNP list (positions, alleles) |
| `v66.p1_HO.aadr.PUB.anno` | Annotations (Y ISOGG, mtDNA calls for haplogroup tables) |

Expected layout:

```
data/references/aadr/
├── v66.p1_HO.aadr.patch.PUB.geno
├── v66.p1_HO.aadr.patch.PUB.ind
├── v66.p1_HO.aadr.patch.PUB.snp
└── v66.p1_HO.aadr.PUB.anno
```

The `.geno`, `.ind`, and `.snp` files share the `v66.p1_HO.aadr.patch.PUB` stem; the annotation file uses a slightly different name (`v66.p1_HO.aadr.PUB.anno`, without `patch`).

Without the first three files, population and ancestry comparisons are unavailable. Without `.anno`, haplogroup **group percentage** rows are omitted (your sample haplogroup calls still work from VCF markers).

Verify the install:

```bash
python main.py list-pops
python main.py list-references
```

## Optional reference files

These are **not** required for core AADR comparisons (Greek, Chinese, Persian, Indian HO labels, and built-in hominin extracts already work once AADR is present):

- **Liftover chain** — `data/references/liftover/hg38ToHg19.over.chain.gz` for GRCh38 VCFs. The app can download this automatically on first use, or you can fetch it from [UCSC](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz).
- **Full hominin VCFs** — `data/references/hominin/` (multi-GB; compact `.vcf` subsets are already in the repo)
- **Extra caste frequency tables** — `data/references/caste/`

See [`data/references/DOWNLOADS.md`](data/references/DOWNLOADS.md) for filenames, sources, and notes on South Asian community labels.

## Requirements

- **Python 3.12+** (stdlib HTTP server; no Flask/Django)
- **numpy** ≥ 1.24 (see `requirements.txt`)

## Run with Docker

Build and start the web UI on port **8765**, mounting your local `data/` directory so AADR files and uploads persist:

```bash
make serve
```

Equivalent manual commands:

```bash
docker build -t dna-compare .
docker run --rm --name dna-compare -p 8765:8765 \
  -v "$(pwd)/data:/app/data" \
  dna-compare
```

Open **http://127.0.0.1:8765/** in a browser. Override the host port with `make serve PORT=9000`.

## Run with Python directly

```bash
python3 -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt

# Web dashboard (default command)
python main.py serve

# Custom bind address
python main.py serve --host 0.0.0.0 --port 8765 --no-open
```

### CLI analysis (no server)

Analyze a file and open an HTML report:

```bash
python main.py analyze data/samples/demo.snps.vcf
python main.py analyze my_genome.vcf --other relative.vcf   # relatedness
python main.py analyze my_genome.vcf --json                   # API-shaped JSON
python main.py analyze my_genome.vcf --text                   # terminal summary
```

Useful flags:

| Flag | Purpose |
| --- | --- |
| `--groups greek,chinese,persian,caste` | Population packs to include |
| `--pops Brahmin,Yadava` | Extra exact AADR group IDs |
| `--assembly auto\|GRCh37\|GRCh38` | Force or auto-detect VCF assembly |
| `--no-hominin`, `--no-caste`, … | Disable comparison sections |
| `--html path/to/report.html` | Write report without opening browser |

Other commands:

```bash
python main.py list-pops tamil      # search AADR population labels
python main.py list-references      # print reference layout
python main.py fetch-references     # rebuild compact hominin extracts from AADR
```

## Project layout

```
.
├── main.py                 # CLI entry point
├── dna_compare/            # parsers, comparisons, web UI, report generator
├── data/
│   ├── references/
│   │   ├── aadr/           # ← place downloaded AADR files here
│   │   ├── hominin/        # archaic reference VCFs / TSVs
│   │   ├── caste/          # community metadata and optional frequency tables
│   │   └── liftover/       # hg38→hg19 chain (optional; auto-downloadable)
│   ├── samples/            # bundled demo VCFs / raw files
│   └── uploads/            # generated reports and uploaded files (gitignored)
├── Dockerfile
├── Makefile
└── requirements.txt
```

## Tests

```bash
pip install pytest
pytest
```

Some tests skip automatically when AADR files are not installed.

## Rust acceleration (optional)

CPU-heavy steps (AADR TGENO reads, HO SNP alignment, batch PGS scoring, kinship) can use a native extension when built:

```bash
pip install maturin
maturin develop --release
```

Without the extension, pure-Python fallbacks still work. Docker images build the extension automatically. Set `DNA_COMPARE_TIMINGS=1` to include per-stage timing in analyze JSON output.

After placing AADR files, pre-build population frequency caches in **one TGENO scan** (recommended before first analyze):

```bash
python main.py build-freq-cache
python main.py build-freq-cache --force   # rebuild existing aadr_freq_*.npy files
```

**Python 3.14:** This repo uses PyO3 0.25+, which supports 3.14. If you still see `newer than PyO3's maximum supported version`, upgrade dependencies (`git pull`) or use a 3.12/3.13 venv. As a temporary workaround on an older PyO3 pin:

```bash
PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1 maturin develop --release
```

## Notes

- The **first** comparison after install can take several minutes while SNP frequency caches are built under `data/references/cache/` unless you run `python main.py build-freq-cache` first (single pass over AADR individuals, much faster than rebuilding group-by-group).
- AADR HO is **GRCh37 / hg19**. If your VCF header indicates GRCh38, coordinates are lifted before autosomal overlap; haplogroup markers also match by rsID and published hg38 positions.
- Community/caste bars show AADR HO labels that exist in the panel — not a complete or ranked list of social groups.
- For licensing and citation of AADR, follow the terms on the [Harvard Dataverse dataset page](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/FFIDCW) and the [AADR documentation](https://reich.hms.harvard.edu/allen-ancient-dna-resource-aadr-downloads).
