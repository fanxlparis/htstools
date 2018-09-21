#!/usr/bin/env bash


# =============================================================================
#              Apply BQSR according to the GATK Best Practices
# =============================================================================


# !!!!!!!!!!!!!!! YOUR CUSTOM PATHS MIGHT BE REQUIRED HERE !!!!!!!!!!!!!!!!!!!!
readonly input_bam="ERR174310.aln.sort.dupmark.bam"
# Output file will be: $output_prefix.recal.bam
readonly output_prefix="ERR174310.aln.sort.dupmark"
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
readonly dbsnp_vcf="$gatk_bundle_path/dbsnp_138.b37.vcf"

# Check if the reference is indexed.
if [ ! -f "${ref_fasta}.fai" ]; then
    printf "$SCRIPT_NAME: error: reference '$ref_fasta' not indexed\n"
    printf "$SCRIPT_NAME: (hint: use 'samtools faidx ref.fasta')\n"
    exit -1
fi


# -----------------------------------------------------------------------------
# Required programs
# -----------------------------------------------------------------------------

readonly gatk="/project/omics/install/gatk-4.0.8.1/gatk"


# -----------------------------------------------------------------------------
# Do it!
# -----------------------------------------------------------------------------

# Run the BaseRecalibrator.
printf "$SCRIPT_NAME: running GATK's BaseRecalibrator\n"
$gatk BaseRecalibrator \
    --reference $ref_fasta \
    --input $input_bam \
    --known-sites $dbsnp_vcf \
    --output $output_prefix.bam.recal_data.table
if [ $? -ne 0 ]; then
    printf "$SCRIPT_NAME: error: GATK BaseRecalibrator returned with non-zero status\n"
    exit -1
fi

# Create recalibrated BAM.
printf "$SCRIPT_NAME: creating recalibrated BAM\n"
$gatk PrintReads \
    --reference $ref_fasta \
    --input $input_bam \
    --BQSR $input_bam.recalibration_report.grp \
    --output $output_prefix.recal.bam
    --emit_original_quals # emit the OQ tag with the original base qualities
if [ $? -ne 0 ]; then
    printf "$SCRIPT_NAME: error: GATK PrintReads returned with non-zero status\n"
    exit -1
fi
