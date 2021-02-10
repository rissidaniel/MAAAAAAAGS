 ############################ rename fastq reads first
# conda activate metagenomics


        echo 'ilumina utils'
#mkdir repaired_illumina_utils
iu-gen-configs config_file_illuma_utils.txt -o repaired_illumina_utils

# Filtering step
for ini in repaired_illumina_utils/*.ini; do iu-filter-quality-minoche $ini --ignore-deflines; done


mkdir spades_error_correction/

        echo "Spades error correction"

spades -1 10MGS_DP_MGs14_S2_R1_001.fastq -2 10MGS_DP_MGs14_S2_R2_001.fastq --only-error-correction --meta -t 30 -o 10MGS_s14_S2_spades_error_correction
spades -1 10MGS_DP_MGs23_S9_R1_001.fastq -2 10MGS_DP_MGs23_S9_R2_001.fastq --only-error-correction --meta -t 30 -o 10MGS_s23_S9_spades_error_correction
spades -1 10MGS_DP_MGs25_S3_R1_001.fastq -2 10MGS_DP_MGs25_S3_R2_001.fastq --only-error-correction --meta -t 30 -o 10MGS_s25_S3_spades_error_correction
spades -1 10MGS_DP_MGs26_S10_R1_001.fastq -2 10MGS_DP_MGs26_S10_R2_001.fastq --only-error-correction --meta -t 30 -o 10MGS_s26_S10_spades_error_correction
spades -1 10MGS_DP_MGs27_S4_R1_001.fastq -2 10MGS_DP_MGs27_S4_R2_001.fastq --only-error-correction --meta -t 30 -o 10MGS_s27_S4_spades_error_correction
spades -1 10MGS_DP_MGs31_S5_R1_001.fastq -2 10MGS_DP_MGs31_S5_R2_001.fastq --only-error-correction --meta -t 30 -o 10MGS_s31_S5_spades_error_correction
spades -1 10MGS_DP_MGs3_S6_R1_001.fastq -2 10MGS_DP_MGs3_S6_R2_001.fastq --only-error-correction --meta -t 30 -o 10MGS_s3_S6_spades_error_correction
spades -1 10MGS_DP_MGs6_S7_R1_001.fastq -2 10MGS_DP_MGs6_S7_R2_001.fastq --only-error-correction --meta -t 30 -o 10MGS_s6_S7_spades_error_correction
spades -1 10MGS_DP_MGs7_S8_R1_001.fastq -2 10MGS_DP_MGs7_S8_R2_001.fastq --only-error-correction --meta -t 30 -o 10MGS_s7_S8_spades_error_correction
spades -1 10MGS_DP_MGs8_S1_R1_001.fastq -2 10MGS_DP_MGs8_S1_R2_001.fastq --only-error-correction --meta -t 30 -o 10MGS_s8_S1_spades_error_correction
spades -1 8MGs_DP_MGs_11_S7_R1_001.fastq -2 8MGs_DP_MGs_11_S7_R2_001.fastq --only-error-correction --meta -t 30 -o 8MGs_s11_S7_spades_error_correction
spades -1 8MGs_DP_MGs_12_S2_R1_001.fastq -2 8MGs_DP_MGs_12_S2_R2_001.fastq --only-error-correction --meta -t 30 -o 8MGs_s12_S2_spades_error_correction
spades -1 8MGs_DP_MGs_19_S1_R1_001.fastq -2 8MGs_DP_MGs_19_S1_R2_001.fastq --only-error-correction --meta -t 30 -o 8MGs_s19_S1_spades_error_correction
spades -1 8MGs_DP_MGs_22_S10_R1_001.fastq -2 8MGs_DP_MGs_22_S10_R2_001.fastq --only-error-correction --meta -t 30 -o 8MGs_s22_S10_spades_error_correction
spades -1 8MGs_DP_MGs_24_S5_R1_001.fastq -2 8MGs_DP_MGs_24_S5_R2_001.fastq --only-error-correction --meta -t 30 -o 8MGs_s24_S5_spades_error_correction
spades -1 8MGs_DP_MGs_28_S8_R1_001.fastq -2 8MGs_DP_MGs_28_S8_R2_001.fastq --only-error-correction --meta -t 30 -o 8MGs_s28_S8_spades_error_correction
spades -1 8MGs_DP_MGs_30_S9_R1_001.fastq -2 8MGs_DP_MGs_30_S9_R2_001.fastq --only-error-correction --meta -t 30 -o 8MGs_s30_S9_spades_error_correction
spades -1 8MGs_DP_Mgs_5_S6_R1_001.fastq -2 8MGs_DP_Mgs_5_S6_R2_001.fastq --only-error-correction --meta -t 30 -o 8MGs_s5_S6_spades_error_correction
spades -1 GreenSnow_GCCAAT_L006_R1_001.fastq -2 GreenSnow_GCCAAT_L006_R2_001.fastq --only-error-correction --meta -t 30 -o GreenSnow_001_spades_error_correction
spades -1 GreenSnow_GCCAAT_L006_R1_002.fastq -2 GreenSnow_GCCAAT_L006_R2_002.fastq --only-error-correction --meta -t 30 -o GreenSnow_002_spades_error_correction
spades -1 GreenSnow_GCCAAT_L006_R1_004.fastq -2 GreenSnow_GCCAAT_L006_R2_004.fastq --only-error-correction --meta -t 30 -o GreenSnow_004_spades_error_correction
spades -1 GreenSnow_GCCAAT_L006_R1_005.fastq -2 GreenSnow_GCCAAT_L006_R2_005.fastq --only-error-correction --meta -t 30 -o GreenSnow_005_spades_error_correction
spades -1 GreenSnow_GCCAAT_L006_R1_007.fastq -2 GreenSnow_GCCAAT_L006_R2_007.fastq --only-error-correction --meta -t 30 -o GreenSnow_007_spades_error_correction
spades -1 GreenSnow_GCCAAT_L006_R1_008.fastq -2 GreenSnow_GCCAAT_L006_R2_008.fastq --only-error-correction --meta -t 30 -o GreenSnow_008_spades_error_correction
spades -1 GreenSnow_GCCAAT_L006_R1_009.fastq -2 GreenSnow_GCCAAT_L006_R2_009.fastq --only-error-correction --meta -t 30 -o GreenSnow_009_spades_error_correction
spades -1 GreenSnow_GCCAAT_L006_R1_010.fastq -2 GreenSnow_GCCAAT_L006_R2_010.fastq --only-error-correction --meta -t 30 -o GreenSnow_010_spades_error_correction
spades -1 GreenSnow_GCCAAT_L006_R1_011.fastq -2 GreenSnow_GCCAAT_L006_R2_011.fastq --only-error-correction --meta -t 30 -o GreenSnow_011_spades_error_correction
spades -1 RedSnow_ACAGTG_L006_R1_001.fastq -2 RedSnow_ACAGTG_L006_R2_001.fastq --only-error-correction --meta -t 30 -o RedSnow_001_spades_error_correction
spades -1 RedSnow_ACAGTG_L006_R1_002.fastq -2 RedSnow_ACAGTG_L006_R2_002.fastq --only-error-correction --meta -t 30 -o RedSnow_002_spades_error_correction
spades -1 RedSnow_ACAGTG_L006_R1_003.fastq -2 RedSnow_ACAGTG_L006_R2_003.fastq --only-error-correction --meta -t 30 -o RedSnow_003_spades_error_correction
spades -1 RedSnow_ACAGTG_L006_R1_004.fastq -2 RedSnow_ACAGTG_L006_R2_004.fastq --only-error-correction --meta -t 30 -o RedSnow_004_spades_error_correction
spades -1 RedSnow_ACAGTG_L006_R1_005.fastq -2 RedSnow_ACAGTG_L006_R2_005.fastq --only-error-correction --meta -t 30 -o RedSnow_005_spades_error_correction
spades -1 RedSnow_ACAGTG_L006_R1_006.fastq -2 RedSnow_ACAGTG_L006_R2_006.fastq --only-error-correction --meta -t 30 -o RedSnow_006_spades_error_correction
spades -1 RedSnow_ACAGTG_L006_R1_007.fastq -2 RedSnow_ACAGTG_L006_R2_007.fastq --only-error-correction --meta -t 30 -o RedSnow_007_spades_error_correction
spades -1 RedSnow_ACAGTG_L006_R1_008.fastq -2 RedSnow_ACAGTG_L006_R2_008.fastq --only-error-correction --meta -t 30 -o RedSnow_008_spades_error_correction
spades -1 RedSnow_ACAGTG_L006_R1_009.fastq -2 RedSnow_ACAGTG_L006_R2_009.fastq --only-error-correction --meta -t 30 -o RedSnow_009_spades_error_correction
spades -1 RedSnow_ACAGTG_L006_R1_010.fastq -2 RedSnow_ACAGTG_L006_R2_010.fastq --only-error-correction --meta -t 30 -o RedSnow_010_spades_error_correction
spades -1 RedSnow_ACAGTG_L006_R1_011.fastq -2 RedSnow_ACAGTG_L006_R2_011.fastq --only-error-correction --meta -t 30 -o RedSnow_011_spades_error_correction
spades -1 RedSnow_ACAGTG_L006_R1_013.fastq -2 RedSnow_ACAGTG_L006_R2_013.fastq --only-error-correction --meta -t 30 -o RedSnow_013_spades_error_correction


rm */corrected/*_unpaired.00.0_0.cor.fastq.gz
mkdir error_corrected_spades_2/
mv */corrected/*.gz error_corrected_spades_2/
cp /home/drissi/Metagenome/Fastq/iceland/Data/raw/repaired_illumina_utils/10-*.fastq.gz error_corrected_spades_2/
cp /home/drissi/Metagenome/Fastq/iceland/Data/raw/repaired_illumina_utils/Cultive_*.fastq.gz error_corrected_spades_2/


cd error_corrected_spades_2/

# Assemble sequences from the same filter fraction and depth together with Megahit
	echo "Megahit Assembly"

R1_raw=`ls *_R1_*.fastq.gz | python -c 'import sys; print(",".join([x.strip() for x in sys.stdin.readlines()]))'`
R2_raw=`ls *_R2_*.fastq.gz | python -c 'import sys; print(",".join([x.strip() for x in sys.stdin.readlines()]))'`

megahit --preset meta-sensitive  -1 $R1_raw -2 $R2_raw  -o mags.megahit_asm -t 30

cd mags.megahit_asm

# Size select contigs ≥2kb in length
        echo "seqmagick --min-length 2000 "
seqmagick convert --min-length 2000 final.contigs.fa final.contigs_min2000.fasta

# cluster the contigs at 99%
        echo "clustering with cd-hit"

#################### cd-hit-est -i final.contigs_min2000.fasta -o final_min2000_contigs.99.fasta -T 30 -M 500000 -c 0.99 -n 10

# Convert file to AFG format and run minimus2 assembler

        echo "toAmos"

toAmos -s final.contigs_min2000.fasta -o final_min2000_contigs.afg

        echo "minimus2"

minimus2 final_min2000_contigs -D OVERLAP=100 MINID=95


cat final_min2000_contigs.fasta final_min2000_contigs.singletons.seq > final_SECONDARY_contigs.fasta

# Align raw sequences to SECONDARY contigs using Bowtie2

        echo "bowtie2"

bowtie2-build final_SECONDARY_contigs.fasta final_SECONDARY_contigs.bt_index

        echo "bowtie2 alignment: heavy stuff"
#Perform alignment for each sample (site, fraction size, depth) - as used in the assembly AND repeat for each of the 15 datasets
#Future iterations will then convert SAM files to BAM files

cd ../

R1_good=`ls *_R1*.fastq.gz | python -c 'import sys; print(",".join([x.strip() for x in sys.stdin.readlines()]))'`
R2_good=`ls *_R2*.fastq.gz | python -c 'import sys; print(",".join([x.strip() for x in sys.stdin.readlines()]))'`



bowtie2 -q -1 $R1s -2 $R2s -S mags.megahit_asm/final_SECONDARY.final_138_prot_meso.sam -x mags.megahit_asm/final_SECONDARY_contigs.bt_index --no-unal -p 30

cd mags.megahit_asm/

# samtools view -b -S final_SECONDARY.final_138_prot_meso.sam > final_SECONDARY.final_138_prot_meso.bam


seqmagick convert --min-length 7000 final_SECONDARY_contigs.fasta final_SECONDARY_contigs.min7000.fasta

        echo "MetaBAT2 - Binning"

metabat2 -i final_SECONDARY_contigs.min7000.fasta -m 1500 -s 4000 -o bins_metabat2
	echo "finish first part"

mkdir bins_metabat2

mv bins_metabat2*.fa bins_metabat2/

cd bins_metabat2/

metawrap classify_bins -b ./ -o bin_classification -t 40


#### add refining steps before conclude pipeline
