#!/usr/bin/env python

import argparse
import datetime
import fileinput
import re


def main():
    current_year = str(datetime.date.today().year)

    # replace 'Copyright (c) 2016-2017' with '2017-2018' (today is 2018-09-21)
    match = r"(Copyright\s\(c\)\s\d\d\d\d)-\d\d\d\d"
    replace = r"\1-" + current_year

    for line in fileinput.input(inplace=1):
        line = re.sub(match, replace, line.rstrip())
        print(line)


if __name__ == "__main__":
    main()
