# Download and QC


module load sra-toolkit/2.8.1
module load fastqc/0.11.8
module load python/3.6.0 # for MultiQC
module load parallel

cd path_to_repository/

# Run in a unix terminal in the repository directory
# The sra files are downloaded in your home directory under ~/ncbi/public/sra
prefetch --option-file ./data/SraAccList.txt

# Convert to Fastq
fastq-dump --split-files path_to_sra_files/*

# Move fastq files in repository directory
mkdir ./data/fastq-files
mv path_to_sra_files/*fastq ./data/fastq-files

# Multiqc to check quality of files
find ./data/fastq-files -name "SRR*.fastq" > ./data/fastq_files.txt
cat ./data/fastq_files.txt | parallel -j 10 "fastqc {}"

# Combine fastqc output and check quality in the html file
multiqc ./data/fastq-files/ --interactive -n "files-downloaded" -o ./data/fastq-files/

# Move fastqc output into a different directory to run the next part
mkdir ./data/fastqc/
mv ./data/fastq-files/*html  ./data/fastqc/
mv ./data/fastq-files/*zip  ./data/fastqc/


# Combine and QC combined fastq
cd ./data/fastq-files
python ../../functions/combined_fastqs.py ../match-SRX-SRR.txt ../sample-to-merge.txt

fastqc SRX729580_combined_1.fastq
fastqc SRX729580_combined_2.fastq
multiqc . --interactive -n "files-downloaded-merged" -o .

# Downsample
path_to_seqtk_folder/seqtk sample -s100 SRX729580_combined_1.fastq 30000000 > sub30M_SRX729580_combined_1.fastq
path_to_seqtk_folder/seqtk sample -s100 SRX729580_combined_2.fastq 30000000 > sub30M_SRX729580_combined_2.fastq

# Multiqc to check number of reads in the downsampled files
fastqc sub30M_SRX729580_combined_1.fastq
fastqc sub30M_SRX729580_combined_2.fastq
multiqc . --interactive -n "files-downsampled" -o .


