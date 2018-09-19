#!/bin/bash

###############################################################################
# Perform alignment according to the GATK Best Practices             #
###############################################################################

###############################################################################
#                               Command line                                  #
###############################################################################

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 num_threads input_fastq sample"
    exit -1
fi

num_threads=$1
input_fastq=$2
sample=$3
root=$(echo $input_fastq | sed 's/\.[^.]*$//') # strip .fastq

if [ ! -f $input_fastq ]; then printf "Error: Input FASTQ file $input_fastq is not a regular file.\n"; exit -1; fi

###############################################################################
#                                GATK bundle                                  #
###############################################################################

gatk_bundle_path="/phys/intern2/muenteferi/GATK"
ref_fasta="$gatk_bundle_path/human_g1k_v37.fasta"
hapmap_vcf="$gatk_bundle_path/hapmap_3.3.b37.vcf"
omni_vcf="$gatk_bundle_path/1000G_omni2.5.b37.vcf"
KG_vcf="$gatk_bundle_path/1000G_phase1.snps.high_confidence.b37.vcf"
dbsnps_vcf="$gatk_bundle_path/dbsnp_138.b37.vcf"
mills_vcf="$gatk_bundle_path/Mills_and_1000G_gold_standard.indels.b37.vcf"
indels_vcf="$gatk_bundle_path/1000G_phase1.indels.b37.vcf"

if [ ! -f $ref_fasta ]; then printf "Error: File $ref_fasta is not a regular file.\n"; exit -1; fi
if [ ! -f $hapmap_vcf ]; then printf "Error: File $hapmap_vcf is not a regular file.\n"; exit -1; fi
if [ ! -f $omni_vcf ]; then printf "Error: File $omni_vcf is not a regular file.\n"; exit -1; fi
if [ ! -f $KG_vcf ]; then printf "Error: File $KG_vcf is not a regular file.\n"; exit -1; fi
if [ ! -f $dbsnps_vcf ]; then printf "Error: File $dbsnps_vcf is not a regular file.\n"; exit -1; fi
if [ ! -f $mills_vcf ]; then printf "Error: File $mills_vcf is not a regular file.\n"; exit -1; fi
if [ ! -f $indels_vcf ]; then printf "Error: File $indels_vcf is not a regular file.\n"; exit -1; fi

###############################################################################
#                                Executables                                  #
###############################################################################

# Binaries
bwa="/project/dna/install/bwa-0.7.13/bwa"
java="/usr/bin/java"
samtools="/project/dna/install/samtools-1.3/bin/samtools"

# JAR files
GenomeAnalysisTK_jar="/project/dna/install/gatk-3.6/GenomeAnalysisTK.jar"
picard_jar="/project/dna/install/picard-tools-2.4.1/picard.jar"

if [ ! -x $bwa ]; then printf "Error: Binary file $bwa is not executable.\n"; exit -1; fi
if [ ! -x $java ]; then printf "Error: Binary file $java is not executable.\n"; exit -1; fi
if [ ! -x $samtools ]; then printf "Error: Binary file $samtools is not executable.\n"; exit -1; fi
if [ ! -f $GenomeAnalysisTK_jar ]; then printf "Error: JAR file $GenomeAnalysisTK_jar is not a regular file.\n"; exit -1; fi
if [ ! -f $picard_jar ]; then printf "Error: JAR file $picard_jar is not a regular file.\n"; exit -1; fi

###############################################################################
#                       Map and mark duplicates                               #
###############################################################################

# Reference indexing
if [ ! -f "${ref_fasta}.fai" ]; then
    $samtools faidx $ref_fasta
fi

# Generate BWA index
$bwa index -a bwtsw $ref_fasta


#---- map to reference

# BWA MEM alignment
$bwa mem -t $num_threads -M $ref_fasta $input_fastq_1 $input_fastq_2 > $root.aln_bwa.sam

# Sort SAM file and convert to BAM
$java -jar $picard_jar SortSam I=$root.aln_bwa.sam O=$root.aln_bwa.sorted.bam SORT_ORDER=coordinate

#---- mark duplicates

# Mark duplicates in the BAM file
$java -jar $picard_jar MarkDuplicates I=$root.aln_bwa.sorted.bam O=$root.aln_bwa.sorted.dupmark.bam M=$root.dedup_metrics.txt ASSUME_SORTED=true


# Index the BAM file
$java -jar $picard_jar BuildBamIndex I=$root.aln_bwa.sorted.dupmark.rg.bam
