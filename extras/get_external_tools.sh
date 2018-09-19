#!/usr/bin/env bash

# Samtools
curl -LO https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
tar -xjf samtools-1.9.tar.bz2
cd samtools-1.9.tar
./configure --prefix=external_tools/samtools-1.9/
make
make install
cd ..
# Samtools binary: external_tools/samtools-1.9/bin/samtools
# rm samtools-1.9.tar.bz2
# rm -rf samtools-1.9

# GATK
curl -LO https://github.com/broadinstitute/gatk/releases/download/4.0.8.1/gatk-4.0.8.1.zip
unzip gatk-4.0.8.1.zip
mkdir external_tools/gatk-4.0.8.1/
cp gatk-4.0.8.1/gatk external_tools/gatk-4.0.8.1/gatk
# GATK binary: external_tools/gatk-4.0.8.1/gatk
# rm gatk-4.0.8.1.zip
# rm -rf gatk-4.0.8.1
