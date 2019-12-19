#!/bin/bash
#PBS -q submit
#PBS -l nodes=1:ppn=2,mem=40gb
#PBS -l walltime=10:00:00
#PBS -N SRX729580_test
#PBS -o ~/Supporting-scripts-library-size-RNA-Seq/funcions/SRX729580_out
#PBS -e ~/Supporting-scripts-library-size-RNA-Seq/funcions/SRX729580_err

# Run this after creating star index for the hg19 reference genome cnsidering 100bp reads

module load STAR/2.5.2
module load sambamba/0.6.6
module load R/3.5.2
module load vardict/1.5.1 # loads R 3.5.2
module load gatk/3.7.0
module load varscan/2.3.9
module load vcftools/0.1.13
module load samtools/1.6
module load ensembl-vep/89.0
module load vcflib/1.0.0-rc1
module load picard-tools/2.20.2
module load python/3.7.0 # for multiqc

set -e
 
## Setup directories and files
cd path_to_repository/

star_genome100=path_to_genome_directory/star_index_hg19_99
genome_fasta=/wehisan/home/allstaff/q/quaglieri.a/PHD_project/GEO_Leucegene_data/genomes/hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
target_regions=./data/target_regions.bed
GATK_Bundle_file1=path_to_GATK_Bundle_files/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
GATK_Bundle_file2=path_to_GATK_Bundle_files/dbsnp_138.hg19.excluding_sites_after_129.vcf \
GATK_Bundle_file3=path_to_GATK_Bundle_files/1000G_phase1.indels.hg19.sites.vcf
outdir=path_output_variants


# Specify also the path_to_vep/cache/ when necessary

#############
# STAR 1 pass
#############

# compress fastq file
pigz ./data/fastq-files/sub30*combined*.fastq

FQ1=./data/fastq-files/sub30M_SRX729580_combined_1.fastq.gz
FQ2=./data/fastq-files/sub30M_SRX729580_combined_2.fastq.gz
bamout=SRX729580
mkdir -p ./data/fastq-files/star-pass1

Rscript ./functions/run_STAR.R \
--genome_index ${star_genome100} \
--fastqfiles $FQ1,$FQ2 \
--sampleName SRX729580 \
--outdir ./data/fastq-files/star-pass1 --STARmode "1Pass"

# Check alignment stats
multiqc ./data/fastq-files/star-pass1 --interactive -n "star1pass" -o ./data/fastq-files/star-pass1
 
# concatenate splice junctions from all samples from ran in pass1
cat ./data/fastq-files/star-pass1/*SJ.out.tab > ./data/fastq-files/star-pass1/combined_sj.out.tab
# Dobin suggests to remove chrm cause they are usually False positives
awk '!/chrM/' ./data/fastq-files/star-pass1/combined_sj.out.tab > ./data/fastq-files/star-pass1/combined_sj_nochrM.out.tab
 
 
#############
# STAR 2 pass
#############

FQ1=./data/fastq-files/sub30M_SRX729580_combined_1.fastq.gz
FQ2=./data/fastq-files/sub30M_SRX729580_combined_2.fastq.gz
mkdir -p ./data/fastq-files/star-pass2

Rscript ./functions/run_STAR.R \
--genome_index $star_genome100 \
--fastqfiles $FQ1,$FQ2 \
--sampleName SRX729580 \
--RlibPath "path_to_R/3.4" \
--outdir ./data/fastq-files/star-pass2 --STARmode "2PassMulti" \
--sjfile ./data/fastq-files/star-pass1/combined_sj_nochrM.out.tab

# Check alignment stats
multiqc ./data/fastq-files/star-pass2 --interactive -n "star2pass" -o ./data/fastq-files/star-pass2
 
 
# Picard tool function to add read groups to a bamfile
AddOrReplaceReadGroups \
I=./data/fastq-files/star-pass2/SRX729580Aligned.sortedByCoord.out.bam \
O=./data/fastq-files/star-pass2/SRX729580Aligned.sortedByCoord.out.RG.bam \
RGID=SRX729580 \
RGPU=SRX729580 \
RGLB=SRX729580_L1 \
RGPL="illumina" \
RGSM=SRX729580_L1 
 
 
# Add index to final file
sambamba index ./data/fastq-files/star-pass2/SRX729580Aligned.sortedByCoord.out.RG.bam 
  
 
###################################
## Prep with GATK and call variants
###################################
 
# Vardict call variant 
 
Rscript ${functions_dir}/call_variants.R \
--reference_fasta ${genome_fasta} \
--bamfile ${fastqdir}/star-pass2/SRX729580Aligned.sortedByCoord.out.RG.bam \
--sampleName SRX729580 \
--regions ${target_regions} \
--genome_assembly 'GRCh37' \
--VarDict_dir /wehisan/home/allstaff/q/quaglieri.a/software/VarDict \
--variant_caller vardict \
--outputdir_suffix "caller_defaults" \
--VEPcall 'vep --dir_cache path_to_vep/cache/.vep --offline' \
--output_directory ${outdir}/test-variants

# GATK pre-process
Rscript ${functions_dir}/gatk_process_pipe.R \
--reference_fasta ${genome_fasta} \
--bamfile ${fastqdir}/star-pass2/SRX729580Aligned.sortedByCoord.out.RG.bam \
--sampleName SRX729580 \
--knownSites2 ${GATK_Bundle_file1} \
--knownSites1 ${GATK_Bundle_file2} \
--knownSites3 ${GATK_Bundle_file3} 

bamGATK=${fastqdir}/star-pass2/SRX729580Aligned.sortedByCoord.out.RG.split.recalibrated.bam
 
# Mutect2 call variant 
 
Rscript ${functions_dir}/call_variants.R \
--reference_fasta ${genome_fasta} \
--bamfile ${bamGATK} \
--sampleName SRX729580 \
--regions ${target_regions} \
--genome_assembly 'GRCh37' --variant_caller mutect \
--outputdir_suffix "caller_defaults" \
--VEPcall 'vep --dir_cache path_to_vep/cache/.vep --offline' \
--output_directory ${outdir}/test-variants

# Varscan2 call variant 
 
Rscript ${functions_dir}/call_variants.R \
--reference_fasta ${genome_fasta} \
--bamfile ${bamGATK} \
--sampleName SRX729580 \
--regions ${target_regions} \
--genome_assembly 'GRCh37' --variant_caller varscan \
--outputdir_suffix "caller_defaults" \
--VEPcall 'vep --dir_cache path_to_vep/cache/.vep --offline' \
--output_directory ${outdir}/test-variants
