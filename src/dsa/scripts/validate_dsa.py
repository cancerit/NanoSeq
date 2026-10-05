########## LICENCE ##########
# Copyright (c) 2022, 2025, 2026 Genome Research Ltd
#
# Author: Luca Barbon <lb29@sanger.ac.uk>
#
# This file is part of NanoSeq.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# 1. The usage of a range of years within a copyright statement contained within
# this distribution should be interpreted as being equivalent to a list of years
# including the first and last year specified and all consecutive years between
# them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
# 2009, 2011-2012’ should be interpreted as being identical to a statement that
# reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
# statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
# identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
# 2009, 2010, 2011, 2012’.
##########################

#!/usr/bin/env python3

from argparse import ArgumentParser
import json
import logging
import os
import subprocess
import sys


DSA_FILE_NAME = 'dsa.bed.gz'
REPORT_FILE_NAME = 'report.json'


def md5sum(fp: str) -> str:
    result = subprocess.run(
        ['md5sum', fp],
        capture_output=True,
        text=True,
        check=True,
    )

    return result.stdout.split(maxsplit=1)[0]


def main(dsa_fp: str, report_fp: str, checksum: bool = False) -> None:

    # Get DSA file size
    try:
        dsa_size = os.path.getsize(dsa_fp)
    except FileNotFoundError:
        sys.exit(f"DSA table not found at '{dsa_fp}'!")

    # Load report
    try:
        with open(report_fp) as fh:
            report = json.load(fh)
    except FileNotFoundError:
        sys.exit(f"DSA report not found at '{report_fp}'!")
    except json.JSONDecodeError:
        sys.exit(f"Invalid DSA report at '{report_fp}'!")

    # Load expected DSA file size
    compressed_meta: dict | None = report.get('compressed', None)
    if not compressed_meta or not isinstance(compressed_meta, dict):
        sys.exit("Compressed information not available!")

    dsa_size_exp: int | None = compressed_meta.get('size_bytes', None)
    if not isinstance(dsa_size_exp, int):
        sys.exit("Invalid expected size format!")

    # Compare DSA file size with the expected size
    if dsa_size != dsa_size_exp:
        sys.exit(f"Compressed size mismatch ({dsa_size} B, expected {dsa_size_exp} B)!")

    if checksum:
        logging.info("Calculating MD5...")
        dsa_md5_obs: str = md5sum(dsa_fp)
        dsa_md5_exp: str | None = compressed_meta.get('md5', None)

        if not dsa_md5_exp:
            sys.exit("Expected MD5 not found!")

        if dsa_md5_obs != dsa_md5_exp:
            sys.exit(f"Compressed MD5 does not match!")


if __name__ == '__main__':
    p = ArgumentParser(description="Validate DSA table against truncation.")
    p.add_argument('dsa', help=f"Path of {DSA_FILE_NAME} file")
    p.add_argument('report', help=f"Path of {REPORT_FILE_NAME} file")
    p.add_argument('--checksum', action='store_true')
    args = p.parse_args()

    logging.basicConfig(level=logging.INFO)
    main(args.dsa, args.report, checksum=args.checksum)
