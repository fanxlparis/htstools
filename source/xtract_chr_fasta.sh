#!/bin/bash

script_name="$(basename $0)"


# -----------------------------------------------------------------------------
# command line
# -----------------------------------------------------------------------------

if [ "$#" -ne 2 ]; then
    printf "usage: $0 file.[fa|fasta] chromosome\n";
    exit -1;
fi

if [ ! -f $1 ]; then
    printf "$script_name: error: '$1' is not a regular file\n"
    exit -1
fi

root=$(echo $1 | sed 's/\.[^.]*$//') # strip .fa or .fasta
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

if [ ! -f $1.fai ]; then
    printf "$script_name: creating FASTA index file: $1.fai\n"
    $samtools faidx $1
fi

if [ ! -f $root.$chromosome.fasta ]; then
    $samtools faidx $1 $chromosome 1> $root.$chromosome.fasta
    printf "$script_name: wrote output to: $root.$chromosome.fasta\n"
else
    printf "$script_name: $root.$chromosome.fasta already exists\n"
fi
