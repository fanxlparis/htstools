#!/usr/bin/env bash


# =============================================================================
#         Perform alignment according to the GATK Best Practices
# =============================================================================


# !!!!!!!!!!!!!!! YOUR CUSTOM PATHS MIGHT BE REQUIRED HERE !!!!!!!!!!!!!!!!!!!!
readonly num_threads=2
readonly input_fastq_1="ERR174310_1.fastq"
readonly input_fastq_2="ERR174310_2.fastq"
# Output file will be: $output_prefix.aln.sort.dupmark.bam
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

readonly ref_fasta="$gatk_bundle_path/human_g1k_v37.fasta"

# Check if the reference is indexed.
if [ ! -f "${ref_fasta}.fai" ]; then
    printf "$SCRIPT_NAME: error: reference '$ref_fasta' not indexed\n"
    printf "$SCRIPT_NAME: (hint: use 'samtools faidx ref.fasta')\n"
    exit -1
fi


# -----------------------------------------------------------------------------
# Required programs
# -----------------------------------------------------------------------------

readonly bwa="/project/omics/install/bwa-0.7.13/bwa"
readonly gatk="/project/omics/install/gatk-4.0.8.1/gatk"
readonly java="/usr/bin/java"
readonly picard_jar="/project/omics/install/picard-tools-2.18.14/picard.jar"


# -----------------------------------------------------------------------------
# Do it!
# -----------------------------------------------------------------------------

# Generate BWA index.
printf "$SCRIPT_NAME: generating BWA index\n"
$bwa index -a bwtsw $ref_fasta
if [ $? -ne 0 ]; then
    printf "$SCRIPT_NAME: error: BWA index returned with non-zero status\n"
    exit -1
fi

# Map reads to reference.
printf "$SCRIPT_NAME: mapping reads to reference\n"
$bwa mem \
    -t $num_threads \
    -M $ref_fasta \
    $input_fastq_1 \
    $input_fastq_2 \
    > $output_prefix.aln.sam
if [ $? -ne 0 ]; then
    printf "$SCRIPT_NAME: error: BWA mem returned with non-zero status\n"
    exit -1
fi

# Sort SAM file and convert to BAM.
printf "$SCRIPT_NAME: sorting SAM file (while converting it to BAM)\n"
$java -jar $picard_jar SortSam \
    INPUT=$output_prefix.aln.sam \
    OUTPUT=$output_prefix.aln.sort.bam \
    SORT_ORDER=coordinate
if [ $? -ne 0 ]; then
    printf "$SCRIPT_NAME: error: Picard returned with non-zero status\n"
    exit -1
fi

# Mark duplicates in the BAM file.
printf "$SCRIPT_NAME: marking duplicates\n"
$java -jar $picard_jar MarkDuplicates \
    INPUT=$output_prefix.aln.sort.bam \
    OUTPUT=$output_prefix.aln.sort.dupmark.bam \
    METRICS_FILE=$output_prefix.dedup_metrics.txt \
    ASSUME_SORTED=true
if [ $? -ne 0 ]; then
    printf "$SCRIPT_NAME: error: Picard returned with non-zero status\n"
    exit -1
fi

# Add read group.
printf "$SCRIPT_NAME: adding read group\n"
$java -jar $picard_jar AddOrReplaceReadGroups \
    INPUT=$output_prefix.aln.sort.dupmark.bam \
    OUTPUT=$output_prefix.aln.sort.dupmark.rg.bam \
    RGID=1 \
    RGLB=library \
    RGPL=illumina \
    RGPU=platform_unit \
    RGSM=sample_name
