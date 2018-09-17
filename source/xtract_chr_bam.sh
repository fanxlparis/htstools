#!/bin/bash

script_name="$(basename $0)"


# -----------------------------------------------------------------------------
# command line
# -----------------------------------------------------------------------------

if [ "$#" -ne 2 ]; then
    printf "usage: $0 file.bam chromosome\n";
    exit -1;
fi

if [ ! -f $1 ]; then
    printf "$script_name: error: '$1' is not a regular file\n"
    exit -1
fi

root=$(echo $1 | sed 's/\.[^.]*$//') # strip .bam
chromosome=$2


# -----------------------------------------------------------------------------
# executables
# -----------------------------------------------------------------------------

samtools="/project/omics/bin/samtools-1.3"

if [ ! -x $samtools ]; then
    printf "$script_name: binary file '$samtools' is not executable\n"
    exit -1
fi


# -----------------------------------------------------------------------------
# extraction
# -----------------------------------------------------------------------------

if [ ! -f $root.bai ]; then
    printf "$script_name: creating BAM index file: $root.bai\n"
    $samtools index $root.bam $root.bai
fi

if [ ! -f $root.$chromosome.sam ]; then
    $samtools view -h $root.bam $chromosome 1> $root.$chromosome.sam
    printf "$script_name: wrote output to: $root.$chromosome.sam\n"
else
    printf "$script_name: $root.$chromosome.sam already exists\n"
fi
