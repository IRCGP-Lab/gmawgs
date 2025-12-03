GenoMycAnalyzer (GMA)

GenoMycAnalyzer (GMA) is a comprehensive WGS Bioinformatics Pipeline for Mycobacterium tuberculosis (MTB).
It automates QC, Mapping, Variant Calling, and Drug Resistance Profiling.

Features

Input: Raw FastQ files (Paired-end)

Process:

QC (Trimming)

Species Identification

Mapping (H37Rv)

Variant Calling 

Drug Resistance & Lineage 

Output: PDF Report, Integrated Excel Summary, VCF, BAM

Installation

You can install GMA via Conda:

# Add channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels shyun

# Install
conda install gma


Usage

# Run Pipeline
gmawgs -i /path/to/fastq_dir -o /path/to/output_dir -d /path/to/kraken2_db --threads 16


Citation

Kim et al. (2024). "GenoMycAnalyzer: a web-based tool..." BMC Genomics.