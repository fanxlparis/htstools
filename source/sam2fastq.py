#!/usr/bin/env python

import argparse
import os
import sys


def main():
    script_name = 'analyze_qual_sam.py'

    # set up the command line parser
    parser = argparse.ArgumentParser(
        description=('analyze quality scores in a SAM file')
    )
    parser.add_argument(
        'sam_file_name',
        help='SAM file name'
    )
    args = parser.parse_args()

    # retrieve the command line arguments
    sam_file_name = args.sam_file_name

    # check if the user-supplied input file exists
    if not os.path.isfile(sam_file_name):
        print('{}: error: this is not a file: {}'
              .format(script_name, sam_file_name))
        exit(-1)

    # check if the user-supplied SAM file can be suspected to be a SAM file
    if not sam_file_name.endswith(".sam"):
        print("{}: error: SAM file name must end with '.sam'"
              .format(script_name))
        exit(-1)

    # open the input file
    sam_file = open(sam_file_name, 'r')

    # write only the FASTQ components to stdout
    while 1:
        line = sam_file.readline()
        if not line:
            break
        if not line.startswith('@'):
            fields = line.split('\t')
            sys.stdout.write("@" + fields[0] + "\n")  # rname
            sys.stdout.write(fields[9] + "\n")  # seq
            sys.stdout.write("+" + "\n")
            sys.stdout.write(fields[10] + "\n")  # qual

    # close the input file
    close(sam_file)


if __name__ == "__main__":
    main()
