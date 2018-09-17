#!/usr/bin/env python

import argparse
import os
import sys


def main():
    script_name = 'replace-qual-fastq.py'

    # set up the command line parser
    parser = argparse.ArgumentParser(
        description=('replace quality scores in a FASTQ file')
    )
    parser.add_argument(
        'fastq_file_name',
        help='FASTQ file name'
    )
    parser.add_argument(
        'qual_file_name',
        help=('QUAL file name')
    )
    args = parser.parse_args()

    # retrieve the command line arguments
    fastq_file_name = args.fastq_file_name
    qual_file_name = args.qual_file_name

    # check if the user-supplied input files exists
    if not os.path.isfile(fastq_file_name):
        print('{}: error: this is not a file: {}'
              .format(script_name, fastq_file_name))
        exit(-1)
    if not os.path.isfile(qual_file_name):
        print('{}: error: this is not a file: {}'
              .format(script_name, qual_file_name))
        exit(-1)

    # check if the user-supplied FASTQ file can be suspected to be a FASTQ file
    if (not fastq_file_name.endswith(".fastq") and
            not fastq_file_name.endswith(".fq")):
        print("{}: error: FASTQ file name must end with "
              "either '.fastq' or '.fq'"
              .format(script_name))
        exit(-1)

    # open the input files
    fastq_file = open(fastq_file_name, 'r')
    qual_file = open(qual_file_name, 'r')

    # write new FASTQ contents to stdout
    while 1:
        sequence_identifier = fastq_file.readline().rstrip('\n')
        if not sequence_identifier:
            break  # reached end of file, everything ok
        sequence = fastq_file.readline().rstrip('\n')
        description = fastq_file.readline().rstrip('\n')
        quality_scores = fastq_file.readline().rstrip('\n')
        new_quality_scores = qual_file.readline().rstrip('\n')

        sys.stdout.write(sequence_identifier + "\n")
        sys.stdout.write(sequence + "\n")
        sys.stdout.write(description + "\n")
        sys.stdout.write(new_quality_scores + "\n")


if __name__ == "__main__":
    main()
