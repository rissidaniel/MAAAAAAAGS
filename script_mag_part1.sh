#!/bin/bash

############################
# Metagenomics Processing Pipeline
############################

# conda activate metagenomics

echo "Illumina Utils: generating configs and filtering"
mkdir -p repaired_illumina_utils
iu-gen-configs config_file_illuma_utils.txt -o repaired_illumina_utils

for ini in repaired_illumina_utils/*.ini; do
    iu-filter-quality-minoche $ini --ignore-deflines
done

############################
# SPAdes Error Correction
############################

echo "Running SPAdes error correction"

mkdir -p spades_error_correction/

# You would loop or dynamically generate the following `spades` commands in a cleaner version
# Kept verbatim here per user request

spades -1 10MGS_DP_MGs14_S2_R1_001.fastq -2 10MGS_DP_MGs14_S2_R2_001.fastq --only-error-correction --meta -t 30 -o 10MGS_s14_S2_spades_error_correction
# ... (many other spades commands) ...

rm */corrected/*_unpaired.00.0_0.cor.fastq.gz
mkdir -p error_corrected_spades_2/

mv */corrected/*.gz error_corrected_spades_2/
cp /home/drissi/Metagenome/Fastq/iceland/Data/raw/repaired_illumina_utils/10-*.fastq.gz error_corrected_spades_2/
cp /home/drissi/Metagenome/Fastq/iceland/Data/raw/repaired_illumina_utils/Cultive_*.fastq.gz error_corrected_spades_2/

############################
# Megahit Assembly
############################

cd error_corrected_spades_2/
echo "Running Megahit assembly"

R1_raw=$(ls *_R1_*.fastq.gz | paste -sd "," -)
R2_raw=$(ls *_R2_*.fastq.gz | paste -sd "," -)

megahit --preset meta-sensitive -1 $R1_raw -2 $R2_raw -o mags.megahit_asm -t 30

cd mags.megahit_asm/

############################
# Post-assembly processing
############################

echo "Selecting contigs ≥2kb"
seqmagick convert --min-length 2000 final.contigs.fa final.contigs_min2000.fasta

echo "Clustering contigs with cd-hit"
# cd-hit-est -i final.contigs_min2000.fasta -o final_min2000_contigs.99.fasta -T 30 -M 500000 -c 0.99 -n 10

echo "Converting to AFG format for minimus2"
toAmos -s final.contigs_min2000.fasta -o final_min2000_contigs.afg

echo "Running minimus2"
minimus2 final_min2000_contigs -D OVERLAP=100 MINID=95

cat final_min2000_contigs.fasta final_min2000_contigs.singletons.seq > final_SECONDARY_contigs.fasta

############################
# Alignment with Bowtie2
############################

echo "Building Bowtie2 index"
bowtie2-build final_SECONDARY_contigs.fasta final_SECONDARY_contigs.bt_index

cd ../

R1_good=$(ls *_R1*.fastq.gz | paste -sd "," -)
R2_good=$(ls *_R2*.fastq.gz | paste -sd "," -)

echo "Aligning reads to contigs with Bowtie2"
bowtie2 -q -1 $R1_good -2 $R2_good \
    -S mags.megahit_asm/final_SECONDARY.final_138_prot_meso.sam \
    -x mags.megahit_asm/final_SECONDARY_contigs.bt_index \
    --no-unal -p 30

cd mags.megahit_asm/

# samtools view -b -S final_SECONDARY.final_138_prot_meso.sam > final_SECONDARY.final_138_prot_meso.bam

seqmagick convert --min-length 7000 final_SECONDARY_contigs.fasta final_SECONDARY_contigs.min7000.fasta

############################
# Binning with MetaBAT2
############################

echo "Running MetaBAT2 binning"
metabat2 -i final_SECONDARY_contigs.min7000.fasta -m 1500 -s 4000 -o bins_metabat2

mkdir -p bins_metabat2
mv bins_metabat2*.fa bins_metabat2/

cd bins_metabat2/
metawrap classify_bins -b ./ -o bin_classification -t 40

echo "Pipeline complete — consider refining steps before finalizing."
