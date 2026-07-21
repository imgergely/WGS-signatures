# STAR Protocols — WGS mutational-signature analysis

Post-IsoMut analysis for the bioinformatics section of the STAR Protocols manuscript
(Imre & Juhász): filter the IsoMut SNV calls against the ENCODE blacklist, build SBS96
spectra, compare each clone to the SIGNAL and COSMIC catalogues, and generate Figure 2
and Figure S1.

## Repository contents

| File | Description |
|---|---|
| `run_samtools_qc.sh` | Runs per-clone `samtools flagstat` and `samtools coverage`, and calculates the reference-span-weighted mean depth used in Table S2 |
| `isomut_config.py` | Python 2.7 configuration for the joint IsoMut call across the 42 analysed clones |
| `wgs_signature_pipeline.R` | Filters IsoMut SNVs, constructs SBS96 profiles, performs cosine-similarity analyses, and generates figures and Table S3 |
| `Table_S1_clone_key.csv` | Mapping between protocol names, original clone IDs, ENA BioSample accessions, and data availability |
| `Table_S2_perclone_QC.csv` | Per-clone sequencing QC metrics |
| `Table_S3_mutation_counts.csv` | Per-clone SNV counts before and after ENCODE blacklist filtering |
| `figures/fig2A_burden.pdf` | Mutation burden by treatment |
| `figures/fig2B_spectra.pdf` | Mean treatment-level SBS96 spectra |
| `figures/fig2C_cosine.pdf` | Per-clone cosine similarity to five prespecified SIGNAL references |
| `figures/figS1_cosmic_signal.pdf` | Catalogue-wide COSMIC and SIGNAL comparison |

The exact analysed and omitted sample IDs are defined in `Table_S1_clone_key.csv`.

## Software requirements

The versions used for the protocol were:

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

The deposited IsoMut configuration uses the legacy Python 2 wrapper distributed with IsoMut.

## Inputs not included in this repository

- coordinate-sorted, duplicate-marked BAM files;
- BAM index files;
- UCSC `hg38.analysisSet.fa`;
- the corresponding FASTA index, `hg38.analysisSet.fa.fai`;
- ENCODE hg38 blacklist v2;
- the joint IsoMut SNV call set, `all_SNVs.isomut`;
- raw FASTQ sequencing files.

Sequencing reads for the publicly deposited subset are available in the European Nucleotide
Archive under accession **PRJEB102539**.

Sequencing data for additional clones and the IsoMut SNV call set are available from the
lead contact upon request.
