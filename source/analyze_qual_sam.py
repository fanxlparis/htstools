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

    # initialize statistics
    qual_min = sys.maxint
    qual_max = -sys.maxint - 1
    qual_size = 0
    total_line_cnt = 0
    header_line_cnt = 0
    alignment_line_cnt = 0
    qual_dist = []
    for x in range(0, 128):
        qual_dist.append(0)

    # parse SAM file
    while 1:
        line = sam_file.readline()
        if not line:
            break
        if not line.startswith('@'):
            fields = line.split('\t')
            qual = fields[10]
            if len(qual) > 0:
                qual_size += len(qual)
                for q in qual:
                    if ord(q) > qual_max:
                        qual_max = ord(q)
                    if ord(q) < qual_min:
                        qual_min = ord(q)
                    qual_dist[ord(q)] += 1
            else:
                print("{}: error: no quality scores in line {}"
                      .format(total_line_cnt))
                exit(-1)
            alignment_line_cnt += 1
        else:
            header_line_cnt += 1
        total_line_cnt += 1

    # print statistics
    print("{}: SAM file: {}".format(script_name, sam_file_name))
    print("{}:   lines: {}".format(script_name, total_line_cnt))
    print("{}:     header lines: {}".format(script_name, header_line_cnt))
    print("{}:     alignment lines: {}"
          .format(script_name, alignment_line_cnt))
    print("{}:   quality scores: {}".format(script_name, qual_size))
    print("{}:     range (inclusive): [{},{}]"
          .format(script_name, qual_min, qual_max))
    print("{}:     distribution:".format(script_name))
    for x in range(0, 128):
        if (qual_dist[x] != 0):
            print("{}:       {}: {}".format(script_name, x, qual_dist[x]))

    # close the input file
    close(sam_file)


if __name__ == "__main__":
    main()
