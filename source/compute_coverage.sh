#!/usr/bin/env bash


# =========================================================================== #
#                Compute coverage from a SAM/BAM/CRAM file                    #
# =========================================================================== #


# -----------------------------------------------------------------------------
# Command line
# -----------------------------------------------------------------------------

readonly SCRIPT_NAME="$(basename $0)"

if [ "$#" -ne 1 ]; then
    printf "$SCRIPT_NAME: usage: $0 alignments.[sam|bam|cram]\n"
    exit -1
fi

readonly alignments_file=$1

# Check if the input file exists
if [ ! -f $alignments_file ]; then
    printf "$SCRIPT_NAME: error: '$alignments_file' is not a regular file\n"
    exit -1
fi


# -----------------------------------------------------------------------------
# Required programs
# -----------------------------------------------------------------------------

readonly awk="/usr/bin/awk"
readonly samtools="/project/omics/install/samtools-1.3/bin/samtools"


# -----------------------------------------------------------------------------
# Do it!
# -----------------------------------------------------------------------------

printf "$SCRIPT_NAME: computing coverage\n"
$samtools depth $1 | awk '{sum+=$3} END { print "Coverage: ",sum/NR }'
if [ $? -ne 0 ]; then
    printf "$SCRIPT_NAME: error: Samtools|awk returned with non-zero status\n"
    exit -1
fi
