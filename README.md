# RNA_seq_pip

A comprehensive toolchain for analyzing next-generation sequencing (NGS) data, covering RNA-seq, deep sequencing (whole-genome copy number / replication fork analysis), and metagenome mapping for bacterial genomes and synthetic genetic circuits.

## Supported Analysis Types

- **RNA-seq**: Gene expression quantification (FPKM, TPM, RPM), rRNA/tRNA ratio analysis
- **Deep-seq**: Whole-genome coverage analysis for replication fork dynamics and copy number variation
- **Metagenome mapping**: Single-strain mapping against reference genomes with quality statistics
- **Coverage visualization**: Strand-specific read coverage plots for plasmids and genetic circuits
- **Genome annotation**: GenBank/GFF conversion, BLAST-based re-annotation, NCBI/UniProt integration

## Repository Structure

### Core Pipeline

| Script | Purpose |
|---|---|
| `miniBac-seq_pipeline.sh` | Master shell orchestrator: runs cleaning then alignment/quantification |
| `miniBac-seq_raw_data_process_clean_demix_20250401.py` | Raw data QC: UMI removal, adapter trimming via fastp |
| `miniBac-seq_pip_20250105.py` | Main analysis pipeline: alignment (Bowtie2) and gene expression quantification |

### Core Libraries

| Script | Purpose |
|---|---|
| `RNA_seq_analyzer.py` | `DNASeqAnalyzer` and `RNASeqAnalyzer` classes for alignment, counting, and normalization |
| `seq_utility.py` | `BAMile` and `GeneFeature` classes; GFF/FASTA parsing; strand-specific coverage; BAM manipulation |

### Alternative Pipelines

| Script | Purpose |
|---|---|
| `RNA_seq_pip_CP.py` | RNA-seq pipeline for plasmid/circuit samples, driven by `RNASeqSampleInfo.json` |
| `RNA-seq_pip.py` | Legacy RNA-seq pipeline for standard E. coli samples |

### Deep Sequencing

| Script | Purpose |
|---|---|
| `Deep_seq_for_C.py` | Whole-genome coverage analysis for replication fork detection; ori/ter ratio and slope metrics |

### Metagenome

| Script | Purpose |
|---|---|
| `meta_genome_map.py` | Maps metagenome reads against a reference; per-position coverage, mapping quality, and depth stats |

### Annotation Tools

| Script | Purpose |
|---|---|
| `genome_annotation_utilities.py` | GenBank to GFF3/GFF2/FASTA conversion; NCBI Entrez gene-to-UniProt lookup |
| `reannotate_gb.py` | BLAST-based comparative genome re-annotation; transfers annotations between genomes |
| `Script_find_uniprot_and_annotate.py` | Batch UniProt annotation retrieval for gene lists |
| `uniport_api.py` | UniProt REST API helper |

### Visualization

| Script | Purpose |
|---|---|
| `plot_circuits_coverage.py` | Strand-specific read coverage plots for synthetic genetic circuits and plasmids |
| `plot_mRNA_coverage.py` | Origin-of-replication RNA coverage metrics (RNAI/RNAII quantification) |
| `sciplot.py` | Publication-quality matplotlib styling (ggplot2 theme, error bars, scientific notation) |

---

## Environment Setup

A conda environment file is provided. Create and activate the environment:

```shell
conda env create -f env.yml
conda activate bioinfo
```

**Key dependencies:** Python 3.12, Bowtie2, samtools, HTSeq, fastp (v1.1.0), cutadapt, pysam, Biopython, pandas, numpy, scipy, matplotlib, seaborn, scikit-learn, joblib, tqdm, pysamstats.

---

## Pipeline Workflows

### 1. miniBac-seq Pipeline (Full Workflow)

The master pipeline performs two steps: raw data cleaning followed by alignment and quantification.

```shell
bash miniBac-seq_pipeline.sh -i <raw_data_dir> -g <gff_file> -f <fasta_file> -o <clean_data_dir> -p <output_dir> [options]
```

**Required arguments:**

| Flag | Description |
|---|---|
| `-i`, `--input-dir` | Directory containing raw `*.fastq.gz` files (paired-end: `sample_R1.fastq.gz` / `sample_R2.fastq.gz`) |
| `-g`, `--gff-file` | GFF annotation file |
| `-f`, `--fasta-file` | Reference genome FASTA file |

**Optional arguments:**

| Flag | Description | Default |
|---|---|---|
| `-c`, `--cpu` | Number of CPUs for fastp/alignment | `8` |
| `-d`, `--dedup` | Enable deduplication during fastp processing | `false` |
| `-r`, `--run-pipeline` | Run without interactive prompts | — |
| `-t`, `--threading-max` | Max parallel sample processing threads | `8` |
| `-o`, `--clean-output` | Output directory for cleaned FASTQ files | `./cleanData` |
| `-p`, `--ret-output` | Output directory for alignment results | `./` |

**Step 1 — Data Cleaning** (`miniBac-seq_raw_data_process_clean_demix_20250401.py`):

- Removes UMI barcode (7 bp) from the 5' end of read 1
- Trims adapters and low-quality bases using fastp
- Minimum read length after trimming: 10 bp
- Optionally performs deduplication
- Outputs cleaned FASTQ files and fastp QC reports (HTML + JSON)

**Step 2 — Alignment and Quantification** (`miniBac-seq_pip_20250105.py`):

- Builds Bowtie2 index from the reference FASTA
- Aligns paired-end reads with Bowtie2 using the following parameters:
  - `-X 1000`: Maximum fragment length for valid paired-end alignments
  - `-I 18`: Minimum fragment length for valid paired-end alignments
  - `--no-mixed`: Suppress unpaired alignments when paired alignment fails
  - `--no-discordant`: Suppress discordant alignments
- Converts SAM to sorted BAM with samtools
- Counts reads per gene feature using HTSeq-count
- Calculates FPKM, TPM, and RPM normalization
- Outputs per-sample expression statistics CSV files, sorted BAM files, and alignment logs

**Example:**

```shell
bash miniBac-seq_pipeline.sh \
    -i ./raw_data \
    -g ./annotation_file/NH3.23.gff \
    -f ./annotation_file/NH3.23.fasta \
    -o ./cleanData \
    -p ./results \
    -c 16 -r
```

### 2. Plasmid/Circuit RNA-seq Pipeline

For samples involving synthetic genetic circuits or plasmids, use `RNA_seq_pip_CP.py`, which reads sample configuration from `RNASeqSampleInfo.json`.

**JSON configuration schema** (`RNASeqSampleInfo.json`):

```json
{
    "genome_annotations": {
        "<sample_type>": {
            "fasta_ps": "<path_to_fasta>",
            "gff_ps": "<path_to_gff>"
        }
    },
    "samples": {
        "<sample_id>": ["<sample_type>", "<data_directory>"]
    },
    "save_dir": "<output_directory>"
}
```

Each sample can have multiple replicates; the script expects replicate FASTQ files organized by sample folder.

### 3. Deep Sequencing — Replication Fork Analysis

`Deep_seq_for_C.py` analyzes whole-genome sequencing coverage to detect replication fork movement and chromosome organization.

**Key features:**

- Calculates binned coverage depth across the genome
- Maps absolute genomic positions to relative positions anchored at the origin of replication (oriC) and terminus (ter)
- Computes log2-normalized coverage
- Performs linear regression on leading and lagging strand arms to detect replication bias
- Outputs depth statistics CSV and coverage plots (SVG)

**Usage:**

```shell
python Deep_seq_for_C.py \
    --parent_dir <clean_reads_dir> \
    --exp_dir <output_dir> \
    [--ref_ps <reference_fasta>] \
    [--ori_site <oriC_position>] \
    [--ter_site <ter_position>] \
    [--bin_length <bin_size>]
```

**Arguments:**

| Flag | Required | Default | Description |
|---|---|---|---|
| `--parent_dir` | ✅ yes | — | Directory containing sample subdirectories with cleaned FASTQ (`.gz`) files |
| `--exp_dir` | ✅ yes | — | Output directory for alignment results and coverage plots |
| `--ref_ps` | no | `annotation_file/Genome_ref/1655_genome_Liu_lab_20220322.fa` | Reference genome FASTA file |
| `--ori_site` | no | `1` | Genomic position of oriC |
| `--ter_site` | no | `2305111` | Genomic position of the terminus site |
| `--bin_length` | no | `5000` | Bin size (bp) for coverage binning |

**Example:**

```shell
python Deep_seq_for_C.py \
    --parent_dir /data/clean \
    --exp_dir /data/deep_seq_results \
    --ori_site 3925859 \
    --ter_site 2305111 \
    --bin_length 5000
```

### 4. Coverage Visualization

**Strand-specific circuit coverage** (`plot_circuits_coverage.py`):

- Generates bidirectional coverage plots (forward/reverse strands) for specified genomic regions
- Overlays gene annotations with directional arrows
- Color-coded by strand: pink for forward (+), cyan for reverse (-)
- Input: BAM file (binary format), GFF annotation, genomic range (start-end)

**Origin-of-replication RNA coverage** (`plot_mRNA_coverage.py`):

- Calculates RNAI and RNAII coverage from forward/reverse strand reads
- Computes counts and CPKM (counts per kilobase per million) for replication control RNA elements
- Outputs per-sample CSV with ori coverage metrics

### 5. Metagenome Mapping

`meta_genome_map.py` maps metagenome reads to a single reference genome with comprehensive quality statistics:

- Per-position coverage depth (raw and binned)
- Mapping quality distribution along the genome
- High/low quality read ratio per genomic region (MAPQ thresholds: 30 and 20)

### 6. Genome Annotation Utilities

**Format conversion** (`genome_annotation_utilities.py`):

```shell
python genome_annotation_utilities.py <input.gb>
```

- Converts GenBank (`.gb`) files to GFF3, GFF2/GTF, and FASTA
- Queries NCBI Entrez to find UniProtKB IDs from NCBI Gene IDs

**Comparative re-annotation** (`reannotate_gb.py`):

- Uses BLAST to compare a reference genome against a target genome
- Identifies unique and shared genes
- Transfers annotations (db_xref, protein_id, gene_synonym) from reference to target
- Accounts for indels and strand orientation
- Outputs an updated GenBank file and an Excel summary report

---

## Manual Adapter Trimming and Demultiplexing

If the automated fastp-based cleaning is insufficient, manual cutadapt commands are available:

**Quality filtering and adapter trimming:**

```shell
cutadapt -q 20,20 --minimum-length 100:100 --max-n 3 --pair-filter=any \
    -o "${sample}_QF_R1.fastq.gz" -p "${sample}_QF_R2.fastq.gz" \
    "${sample}_R1.fastq.gz" "${sample}_R2.fastq.gz" > "${sample}_logs/QF.log"
```

**Demultiplexing** (when adapters are strictly anchored at the 5' end of R1):

```shell
cutadapt --pair-filter=any --no-indels --minimum-length 20 --times 1 -e 0.15 --overlap 7 \
    -g file:adaptor_R1.fasta -A file:adaptor_R2.fasta \
    -o "end5-{name}.R1.fastq.gz" -p "end5-{name}.R2.fastq.gz" \
    "${sample}_QF_R1.fastq.gz" "${sample}_QF_R2.fastq.gz" > "${sample}_logs/demultiplexing.log"
```

---

## Configuration Files

| File | Purpose |
|---|---|
| `env.yml` | Conda environment specification with all dependencies |
| `RNASeqSampleInfo.json` | Sample metadata and genome annotation paths for circuit/plasmid pipelines |
| `entrez.json` | NCBI Entrez API credentials (email and API key) |

### `entrez.json` format

```json
{
    "ENTREZ_EMAIL": "your_email@example.com",
    "ENTREZ_API_KEY": "your_api_key"
}
```

---

## Input/Output Summary

### Input Formats

| Format | Used By |
|---|---|
| Paired-end FASTQ (`.fastq.gz`, `.fq.gz`) | All pipelines |
| FASTA (`.fasta`, `.fa`, `.fna`) | Reference genomes |
| GFF (`.gff`, `.gff3`) | Gene annotations |
| GenBank (`.gb`, `.gbff`) | Annotation conversion, re-annotation |
| JSON (`.json`) | Sample configuration, Entrez credentials |

### Output Formats

| Format | Produced By |
|---|---|
| Cleaned FASTQ (`.fastq.gz`) | `miniBac-seq_raw_data_process_clean_demix_20250401.py` |
| fastp QC reports (`.html`, `.json`) | `miniBac-seq_raw_data_process_clean_demix_20250401.py` |
| Sorted BAM + BAI index | `RNA_seq_analyzer.py` |
| Expression statistics CSV (FPKM, TPM, RPM) | `RNA_seq_analyzer.py` |
| Coverage plots (SVG, PNG) | `plot_circuits_coverage.py`, `plot_mRNA_coverage.py`, `Deep_seq_for_C.py` |
| Re-annotated GenBank (`.gb`) | `reannotate_gb.py` |
| GFF3 / GFF2 / FASTA | `genome_annotation_utilities.py` |

---

## Directory Layout

```
RNA_seq_pip/
    miniBac-seq_pipeline.sh          # Master shell pipeline
    miniBac-seq_pip_20250105.py      # Alignment and quantification
    miniBac-seq_raw_data_process_clean_demix_20250401.py  # Raw data QC
    RNA_seq_analyzer.py              # Core analyzer classes
    seq_utility.py                   # Sequence utility library
    RNA_seq_pip_CP.py                # Circuit/plasmid pipeline
    RNA-seq_pip.py                   # Legacy RNA-seq pipeline
    Deep_seq_for_C.py                # Replication fork analysis
    meta_genome_map.py               # Metagenome mapping
    plot_circuits_coverage.py        # Circuit coverage plots
    plot_mRNA_coverage.py            # ori RNA coverage plots
    sciplot.py                       # Plot styling
    genome_annotation_utilities.py   # Annotation conversion
    reannotate_gb.py                 # BLAST re-annotation
    Script_find_uniprot_and_annotate.py  # UniProt annotation
    uniport_api.py                   # UniProt API helper
    env.yml                          # Conda environment
    RNASeqSampleInfo.json            # Sample configuration
    entrez.json                      # NCBI API credentials
    annotation_file/                 # Reference genomes and annotations
    example_data/                    # Example sequencing data
    deprecated/                      # Deprecated scripts
    Bin/                             # Binary/helper scripts
```

