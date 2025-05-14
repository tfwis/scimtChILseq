### README

#### About the processing/ Directory

This directory contains scripts and auxiliary files necessary for processing sci-mtChIL-seq FASTQ files to obtain unique fragments per cell. 

- **runme.bash**: A shell script that orchestrates the entire processing workflow.
- **chilmap.bash**: A shell script that generates cell barcodes and maps the reads using FASTQ files (Read1, Read2, Index1, Index2).
- **chilrmdups.R**: An R script that removes PCR duplicates and IVT variants, and creates BED files form the BAM file.
- **chilqcs.R**: An R script that remove unassigned barcodes and double barcodes from the BED files.
- **chr.tsv**: A file required by `chilmap.bash` that specifies the chromosomes to be analyzed.

#### Required Software

The following software and tools must be installed to execute the scripts in this directory:

- **seqkit**: v0.16.1
- **umi_tools**: v1.1.1
- **Trim Galore**: v0.6.10
- **Bowtie2**: v2.3.5.1
- **SAMtools**: v1.18
- **R**: v4.4.2
  - Required R packages:
    - GenomicRanges
    - GenomicAlignments
    - furrr
    - tidyverse
    - stringdist
    - mixtools
    - rtracklayer

#### Usage

1. Place the four FASTQ files (Read1, Read2, Index1, Index2) under the `fastq/` directory. These files are uploaded to GEO.
2. Place the file containing cell barcode candidates under the `processing/` directory. This file is uploaded to Zenodo.
3. Edit the file names and paths for `cellbcs` and Bowtie2 index files in `runme.bash` as needed.
4. Execute `runme.bash`.
