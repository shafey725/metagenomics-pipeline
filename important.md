# Metagenomics Pipeline (MGnify v5.0)

This repository contains a Bash script (`pipeline.sh`) that automates the processing of paired-end metagenomic sequencing data using the MGnify v5.0 reference databases.

## Features
- Automatic download & extraction of:
  - SILVA SSU and LSU databases (skips if already downloaded)
  - Rfam rRNA and other ncRNA models
- Paired-end read merging (SeqPrep)
- Read trimming (Trimmomatic)
- Quality control (FastQC)
- rRNA detection and masking (Infernal `cmsearch` + bedtools)
- Gene prediction (FragGeneScan)
- Taxonomic classification with MAPseq for SSU and LSU
- Output visualization with Krona

## Requirements
- Linux/Ubuntu environment
- Installed tools:
  - `SeqPrep`
  - `Trimmomatic`
  - `FastQC`
  - `seqkit`
  - `bedtools`
  - `cmsearch` (Infernal)
  - `FragGeneScan`
  - `MAPseq`
  - `ktImportText` (Krona)
  - `wget`, `tar`, `awk`


Missing / Future Improvements
Parallel processing for multiple samples (GNU parallel support)

Logging system for each sample's run

Automated installation of required dependencies

Error handling for missing tools

Configurable parameters (adapter file, quality thresholds, database versions)
doesn't include all pipeline v5 database only rfam/models 
## Usage
```bash
bash pipeline.sh
