#/bin/bash

script_name="$(basename $0)"


# -----------------------------------------------------------------------------
# input files
# -----------------------------------------------------------------------------

reference_fasta=""
input_bam=""
dbsnp_vcf=""

if [ ! -f $reference_fasta ]; then
    printf "$script_name: error: '$reference_fasta' is not a regular file\n"
    exit -1
fi

if [ ! -f $input_bam ]; then
    printf "$script_name: error: '$input_bam' is not a regular file\n"
    exit -1
fi

if [ ! -f $dbsnp_vcf ]; then
    printf "$script_name: error: '$dbsnp_vcf' is not a regular file\n"
    exit -1
fi


# -----------------------------------------------------------------------------
# executables
# -----------------------------------------------------------------------------

gatk="/project/omics/bin/gatk-4.0"

if [ ! -x $gatk ]; then
    printf "$script_name: binary file '$gatk' is not executable\n"
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
