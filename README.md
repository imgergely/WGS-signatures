# WGS mutational signature analysis


## Repository contents

| File | Description |
|---|---|
| `run_samtools_qc.sh` | Runs per-clone `samtools flagstat` and `samtools coverage`, and calculates the reference-span-weighted mean depth used in Table S2 |
| `isomut_config.py` | Python 2.7 configuration for the joint IsoMut call across the 42 analyzed clones |
| `wgs_signature_pipeline.R` | Implements:  cohort validation, blacklist filtering, SBS96 analysis, reference comparison, and output generation |
| `Table_S1_clone_key.csv` | Mapping between protocol names, original clone IDs, ENA BioSample accessions, and data availability |
| `Table_S2_perclone_QC.csv` | Per-clone sequencing QC metrics |
| `Table_S3_mutation_counts.csv` | Per-clone SNV counts before and after ENCODE blacklist filtering |
| `figures/fig2A_burden.pdf` | Filtered SNV burden by treatment |
| `figures/fig2B_spectra.pdf` | Mean treatment-level SBS96 spectra |
| `figures/fig2C_cosine.pdf` | Per-clone cosine similarity to five prespecified SIGNAL references |
| `figures/figS1_cosmic_signal.pdf` | Catalogue-wide COSMIC and SIGNAL comparison |

`Table_S1_clone_key.csv` defines the exact analyzed and omitted sample IDs. BAP2KO3 is retained in Table S1 as an omitted library but is not included in the joint IsoMut call.

## Expected directory layout

```text
.
├── bams/
│   └── <original_clone_ID>.final.bam
├── input/
│   └── all_SNVs.isomut.gz
├── reference/
│   ├── hg38.analysisSet.fa
│   ├── hg38.analysisSet.fa.fai
│   └── hg38-blacklist.v2.bed.gz
├── Table_S1_clone_key.csv
├── isomut_config.py
├── run_samtools_qc.sh
└── wgs_signature_pipeline.R
```

The BAM and reference files are not included in the repository.

## Workflow

### 1. Run per-clone SAMtools QC

From the repository root:

```bash
chmod +x run_samtools_qc.sh
./run_samtools_qc.sh
```

The script writes one `flagstat` and one `coverage` file per BAM under `qc/`, plus `qc/mean_depth_summary.tsv`. 

### 2. Index the reference and BAM files

```bash
if [[ ! -f reference/hg38.analysisSet.fa.fai ]]; then
    samtools faidx reference/hg38.analysisSet.fa
fi

for bam in bams/*.final.bam; do
    if [[ ! -f "${bam}.bai" && ! -f "${bam%.bam}.bai" ]]; then
        samtools index "$bam"
    fi
done
```

The loop may index an excluded BAM if it remains in `bams/`; however, `isomut_config.py` passes only the 42 analyzed BAMs to IsoMut.

### 3. Run the joint IsoMut call

Edit the path constants at the top of `isomut_config.py` as required, then run it with Python 2.7 from the IsoMut package directory:

```bash
python2 isomut_config.py
```

All 42 analyzed clones must be processed in one execution. IsoMut identifies clone-specific variants by comparing each clone with the complete panel; changing the panel changes which variants are retained as clone-specific.

The downstream analysis uses the SNV output only. Compress the validated table and place it at `input/all_SNVs.isomut.gz`.

### 4. Run the complete downstream R analysis

From the repository root:

```bash
Rscript wgs_signature_pipeline.R
```

The script executes all downstream stages in one R session.

- **Step 32:** load, map, and validate the IsoMut SNV call set;
- **Step 33:** remove ENCODE blacklist overlaps and assemble per-clone counts;
- **Step 34:** validate the final 42-clone cohort and write Table S3;
- **Step 35:** construct and normalize per-clone SBS96 profiles;
- **Step 36:** calculate treatment-level mean SBS96 profiles;
- **Step 37:** load and validate SIGNAL and COSMIC reference matrices;
- **Step 38:** resolve prespecified references and calculate cosine similarities;
- **Step 39:** generate Figure 2A–C and Figure S1, and confirm the expected outputs.

## R outputs

The R script writes the following files to `results/`:

- `Table_S3_mutation_counts.csv`;
- `fig2A_burden.pdf`, `.png`, and `.svg`;
- `fig2B_spectra.pdf`, `.png`, and `.svg`;
- `fig2C_cosine.pdf`, `.png`, and `.svg`;
- `figS1_cosmic_signal.pdf`.

## Software versions used

- BWA 0.7.17
- Sambamba 0.7.1
- Picard 2.18.9
- SAMtools 1.13
- Python 2.7
- IsoMut
- R 4.6.1
- Bioconductor 3.23
- GenomicRanges 1.64.0
- rtracklayer 1.72.0
- Biostrings 2.80.1
- GenomeInfoDb 1.48.0
- MutationalPatterns 3.22.0
- BSgenome.Hsapiens.UCSC.hg38 1.4.5
- ggplot2 4.0.3
- ggdendro 0.2.0
- RColorBrewer 1.1-3

The IsoMut configuration uses the legacy Python 2 wrapper distributed with IsoMut.

## Tested operating environment

The computational workflow was executed on Ubuntu 22.04.5 LTS using the x86_64 architecture. Other operating systems were not tested.

## Data availability

Sequencing reads for the publicly deposited subset are available in the European Nucleotide Archive under accession **PRJEB102539**. Table S1 records the availability of every clone. Sequencing data for additional clones are available from the lead contact upon request.
