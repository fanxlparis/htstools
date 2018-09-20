#!/usr/bin/env python

import argparse
import os
import sys


def main():
    script_name = 'validate_fastq.py'

    # set up the command line parser
    parser = argparse.ArgumentParser(
        description=('validate a FASTQ file')
    )
    parser.add_argument(
        'file_name',
        help='file name'
    )
    args = parser.parse_args()

    # retrieve the command line arguments
    file_name = args.file_name

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

    # set the allowed sequence alphabet
    sequence_alphabet = "ACGTNacgtn"
    print("{}: sequence alphabet: {}".format(script_name, sequence_alphabet))

    # set up counters and statistics
    line_cnt = 0
    quality_scores_min = sys.maxint
    quality_scores_max = -sys.maxint - 1
    quality_scores_cnt = 0

    while 1:
        # read a sequence identifier line
        sequence_identifier = file.readline().rstrip('\n')
        if not sequence_identifier:
            break  # reached end of file, everything ok

        # check the sequence identifier
        if len(sequence_identifier) > 0:
            if len(sequence_identifier) < 2:
                sys.exit("{}: {}:{}: "
                         "sequence identifier is too short"
                         .format(script_name, file_name, line_cnt))
            if sequence_identifier[0] != '@':
                sys.exit("{}: {}:{}: "
                         "sequence identifier does not start with '@'"
                         .format(script_name, file_name, line_cnt))
        else:
            sys.exit("{}: {}:{}: "
                     "sequence identifier is empty"
                     .format(script_name, file_name, line_cnt))

        # increment line count
        line_cnt += 1

        # read a sequence line
        sequence = file.readline().rstrip('\n')
        if not sequence:
            sys.exit("{}: {}:{}: "
                     "record is not complete"
                     .format(script_name, file_name, line_cnt))

        # check the sequence
        if len(sequence) == 0:
            sys.exit("{}: {}:{}: "
                     "sequence is empty"
                     .format(script_name, file_name, line_cnt))
        else:
            for s in sequence:
                if s not in sequence_alphabet:
                    sys.exit("{}: {}:{}: "
                             "invalid character '{}' in sequence"
                             .format(script_name, file_name, line_cnt, s))

        # increment line count
        line_cnt += 1

        # read a description line
        description = file.readline().rstrip('\n')
        if not description:
            sys.exit("{}: {}:{}: "
                     "record is not complete"
                     .format(script_name, file_name, line_cnt))

        # check the description
        if len(description) > 0:
            if description[0] != '+':
                sys.exit("{}: {}:{}: "
                         "description does not start with '+'"
                         .format(script_name, file_name, line_cnt))
            if len(description) > 1:
                if description[1:] != sequence_identifier[1:]:
                    sys.exit("{}: {}:{}: "
                             "description does not match sequence identifier"
                             .format(script_name, file_name, line_cnt))
        else:
            sys.exit("{}: {}:{}: "
                     " description is empty"
                     .format(script_name, file_name, line_cnt))

        # increment line count
        line_cnt += 1

        # read a quality scores line
        quality_scores = file.readline().rstrip('\n')
        if not quality_scores:
            sys.exit("{}: {}:{}: "
                     "record is not complete"
                     .format(script_name, file_name, line_cnt))

        # typical ranges:
        #  Sanger         Phred+33   [0,40]
        #  Solexa         Solexa+64  [-5,40]
        #  Illumina 1.3+  Phred+64   [0,40]
        #  Illumina 1.5+  Phred+64   [0,40]
        #    with
        #      0 = unused
        #      1 = unused
        #      2 = Read Segment Quality Control Indicator ('B')
        #  Illumina 1.8+  Phred+33   [0,41]

        # check the quality scores
        if len(quality_scores) > 0:
            if len(quality_scores) != len(sequence):
                sys.exit("{}: {}:{}: "
                         "quality scores length does not match sequence length"
                         .format(script_name, file_name, line_cnt))
            for q in quality_scores:
                if ord(q) < 33 or ord(q) > 126:
                    sys.exit("{}: {}:{}: "
                             "invalid character '{}'='{}' in quality string"
                             .format(script_name,
                                     file_name,
                                     line_cnt,
                                     q,
                                     ord(q)))
                if ord(q) > 104:
                    sys.exit("{}: {}:{}: "
                             "unusually high quality score '{}'='{}'"
                             .format(script_name,
                                     file_name,
                                     line_cnt,
                                     q,
                                     ord(q)))
                if ord(q) > quality_scores_max:
                    quality_scores_max = ord(q)
                if ord(q) < quality_scores_min:
                    quality_scores_min = ord(q)
            quality_scores_cnt += len(quality_scores)
        else:
            sys.exit("{}: {}:{}: "
                     "quality scores are empty"
                     .format(script_name, file_name, line_cnt))

        # increment line count
        line_cnt += 1

    print("{}: quality score range (inclusive): [{},{}]"
          .format(script_name, quality_scores_min, quality_scores_max))
    print("{}: no. of quality scores: {}"
          .format(script_name, quality_scores_cnt))
    print("{}: processed {} lines".format(script_name, line_cnt))
    print("{}: done".format(script_name))


if __name__ == "__main__":
    main()
