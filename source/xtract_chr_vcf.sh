#!/bin/bash

script_name="$(basename $0)"


# -----------------------------------------------------------------------------
# command line
# -----------------------------------------------------------------------------

if [ "$#" -ne 2 ]; then
    printf "usage: $0 file.vcf chromosome\n";
    exit -1;
fi

if [ ! -f $1 ]; then
    printf "$script_name: error: '$1' is not a regular file\n"
    exit -1
fi

root=$(echo $1 | sed 's/\.[^.]*$//') # strip .vcf
chromosome=$2


# -----------------------------------------------------------------------------
# executables
# -----------------------------------------------------------------------------

gatk="/project/omics/bin/gatk-4.0"
samtools="/project/omics/bin/samtools-1.3/"

if [ ! -x $gatk ]; then
    printf "$script_name: binary file '$gatk' is not executable\n"
    exit -1
fi

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

if [ ! -f $root.$chromosome.vcf ]; then
    $gatk \
        -T SelectVariants \
        -R $ref \
        -V $root.vcf \
        -o $root.$chromosome.vcf \
        -L $chromosome
    printf "$script_name: wrote output to: $root.$chromosome.vcf\n"
else
    printf "$root.$chromosome.vcf already exists (not reproducing it)\n"
fi
