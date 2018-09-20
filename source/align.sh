#!/usr/bin/env bash


# =========================================================================== #
#         Perform alignment according to the GATK Best Practices              #
# =========================================================================== #


# !!!!!!!!!!!!!!! YOUR CUSTOM PATHS MIGHT BE REQUIRED HERE !!!!!!!!!!!!!!!!!!!!
readonly num_threads=2
readonly input_fastq_1="ERR174310_1.fastq"
readonly input_fastq_2="ERR174310_2.fastq"
# Output file will be $output_prefix.aln.sort.dupmark.bam
readonly output_prefix="ERR174310"
readonly gatk_bundle_path="/home/voges/tmp/GATK_bundle-2.8-b37"
# !!!!!!!!!!!!!!! YOUR CUSTOM PATHS MIGHT BE REQUIRED HERE !!!!!!!!!!!!!!!!!!!!


# -----------------------------------------------------------------------------
# Command line
# -----------------------------------------------------------------------------

readonly SCRIPT_NAME="$(basename $0)"

if [ "$#" -ne 0 ]; then
    printf "$SCRIPT_NAME: usage: $0\n"
    printf "$SCRIPT_NAME: (hint: set parameters directly in the script)\n";
    exit -1;
fi


# -----------------------------------------------------------------------------
# GATK bundle
# -----------------------------------------------------------------------------

ref_fasta="$gatk_bundle_path/human_g1k_v37.fasta"
dbsnps_vcf="$gatk_bundle_path/dbsnp_138.b37.vcf"

# Check if the reference is indexed.
if [ ! -f "${ref_fasta}.fai" ]; then
    printf "$SCRIPT_NAME: error: reference '$ref_fasta' not indexed\n"
    printf "$SCRIPT_NAME: (hint: use 'samtools faidx ref.fasta')\n"
    exit -1
fi


# -----------------------------------------------------------------------------
# Required programs
# -----------------------------------------------------------------------------

bwa="/project/omics/install/bwa-0.7.13/bwa"
gatk="/project/omics/install/gatk-4.0.8.1/gatk"
java="/usr/bin/java"
picard_jar="/project/omics/install/picard-tools-2.18.14/picard.jar"


# -----------------------------------------------------------------------------
# Do it!
# -----------------------------------------------------------------------------

# Generate BWA index.
$bwa index -a bwtsw $ref_fasta

# Map reads to reference.
$bwa mem \
    -t $num_threads \
    -M $ref_fasta \
    $input_fastq_1 \
    $input_fastq_2 \
    > $output_prefix.aln.sam

# Sort SAM file and convert to BAM.
$java -jar $picard_jar SortSam \
    I=$output_prefix.aln.sam \
    O=$output_prefix.aln.sort.bam \
    SORT_ORDER=coordinate

# Mark duplicates in the BAM file.
$java -jar $picard_jar MarkDuplicates \
    I=$output_prefix.aln.sort.bam \
    O=$output_prefix.aln.sort.dupmark.bam \
    M=$output_prefix.dedup_metrics.txt \
    ASSUME_SORTED=true
