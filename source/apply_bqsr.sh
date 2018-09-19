#/bin/bash


# -----------------------------------------------------------------------------
# Constants
# -----------------------------------------------------------------------------

readonly SCRIPT_NAME="$(basename $0)"


# -----------------------------------------------------------------------------
# Functions
# -----------------------------------------------------------------------------

safe_execute ()
{
    local readonly command=("$@")
    "${command[@]}"
    # "${command[@]}" &> /dev/null
    local readonly return_code=$?
    if [ $return_code -ne 0 ]; then
        cmd="${command[@]}"
        printf "$SCRIPT_NAME: error: '$cmd' returned with non-zero status\n"
        exit -1
    fi
}


# -----------------------------------------------------------------------------
# Parse command line.
# -----------------------------------------------------------------------------

if [ "$#" -ne 0 ]; then
    printf "$SCRIPT_NAME: usage: $0\n";
    exit -1;
fi


# -----------------------------------------------------------------------------
# Input files
# -----------------------------------------------------------------------------

readonly reference_fasta=""
readonly input_bam=""
readonly dbsnp_vcf=""

if [ ! -f $reference_fasta ]; then
    printf "$SCRIPT_NAME: error: '$reference_fasta' is not a regular file\n"
    exit -1
fi

if [ ! -f $input_bam ]; then
    printf "$SCRIPT_NAME: error: '$input_bam' is not a regular file\n"
    exit -1
fi

if [ ! -f $dbsnp_vcf ]; then
    printf "$SCRIPT_NAME: error: '$dbsnp_vcf' is not a regular file\n"
    exit -1
fi


# -----------------------------------------------------------------------------
# Check if required executables are available.
# -----------------------------------------------------------------------------

# !!!!!!!!!!!!!!! YOUR CUSTOM PATHS MIGHT BE REQUIRED HERE !!!!!!!!!!!!!!!!!!!!
readonly gatk="/project/omics/bin/gatk-4.0"
# !!!!!!!!!!!!!!! YOUR CUSTOM PATHS MIGHT BE REQUIRED HERE !!!!!!!!!!!!!!!!!!!!

if [ ! -x $gatk ]; then
    printf "$SCRIPT_NAME: binary file '$gatk' is not executable\n"
    exit -1
fi


# -----------------------------------------------------------------------------
# run the BaseRecalibrator
# -----------------------------------------------------------------------------

$gatk \
    -T BaseRecalibrator \
    -R $reference_fasta \
    -I $input_bam \
    -knownSites $dbsnp_vcf \
    -o $input_bam.recal_data.table


# -----------------------------------------------------------------------------
# create recalibrated BAM
# -----------------------------------------------------------------------------

$gatk \
    -T PrintReads \
    -R $reference_fasta \
    -I $input_bam \
    -BQSR $recalibration_report_grp \
    -o $output_bam
    --emit_original_quals # emit the OQ tag with the original base qualities
