README
================
Anna Quaglieri
19/12/2019

# Access repository gitbook

To read the documentation go to
<https://supporting-scripts-library-size-rna-seq.netlify.com/>.

# Aim of this document

This documentation contains the code to partially reproduce the results
in the manuscript *Finding a suitable library size to call variants in
RNA-Seq*, from the download of the raw sequenced files to variant
calling. It also contains supplementary information regarding
downsampling of the TCGA-LAML samples and the truth set used to compute
sensitivity with this cohort.

The example variant callset obtained with the Leucegene samples only
contains published somatic variants (Lavallée, 2016).

## Software used

  - `R/3.5.2`

Loaded on Unix

  - `sra-toolkit/2.8.1`
  - `fastqc/0.11.8`
  - `python/3.6.0` (for [MultiQC](https://multiqc.info/))
  - `parallel`
  - `STAR/2.5`
  - `sambamba/0.6.6`
  - `picard-tools/2.20.2`
  - `vardict/1.5.1` (loads `R3.5.2`)
  - `gatk/3.7.0`
  - `varscan/2.3.9`
  - `vcftools/0.1.13`
  - `samtools/1.6`
  - `ensembl-vep/89.0`
  - `vcflib/1.0.0-rc1`
