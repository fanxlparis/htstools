#!/usr/bin/env python

import argparse
import os
import sys


def main():
    script_name = 'xtract_part_fastq.py'

    # set up the command line parser
    parser = argparse.ArgumentParser(
        description=('extract read identifier, sequence, description or '
                     'quality scores from a FASTQ file')
    )
    parser.add_argument(
        'file_name',
        help='file name'
    )
    parser.add_argument(
        'part_id',
        help=('part ID (0 = read identifier, 1 = sequence, 2 = description, '
              '3 = quality scores)')
    )
    args = parser.parse_args()

    # retrieve the command line arguments
    file_name = args.file_name
    part_id = args.part_id

    # check if the user-supplied input file exists
    if not os.path.isfile(file_name):
        print('{}: error: this is not a file: {}'
              .format(script_name, file_name))
        exit(-1)

    # check if the user-supplied input file can be suspected to be a FASTQ file
    if not file_name.endswith(".fastq") and not file_name.endswith(".fq"):
        print("{}: error: input file name must end with "
              "either '.fastq' or '.fq'"
              .format(script_name))
        exit(-1)

    # open the input file
    file = open(file_name, 'r')

    # check if part ID is valid
    part_id = int(part_id)
    if not 0 <= part_id <= 3:
        print('{}: error: part_id must be in range [0,3]'.format(script_name))
        exit(-1)

    # write part to stdout
    while 1:
        sequence_identifier = file.readline().rstrip('\n')
        if not sequence_identifier:
            break  # reached end of file, everything ok
        sequence = file.readline().rstrip('\n')
        description = file.readline().rstrip('\n')
        quality_scores = file.readline().rstrip('\n')

        if part_id == 0:
            sys.stdout.write(sequence_identifier + "\n")
        if part_id == 1:
            sys.stdout.write(sequence + "\n")
        if part_id == 2:
            sys.stdout.write(description + "\n")
        if part_id == 3:
            sys.stdout.write(quality_scores + "\n")

    # close the input file
    close(file)


if __name__ == "__main__":
    main()
