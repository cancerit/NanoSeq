#!/usr/bin/env python3

from argparse import ArgumentParser
import json
import os
import sys


DSA_FILE_NAME = 'dsa.bed.gz'
REPORT_FILE_NAME = 'report.json'


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

    # Load expected DSA file size
    compressed_meta: dict | None = report.get('compressed', None)
    if not compressed_meta or not isinstance(compressed_meta, dict):
        sys.exit("Compressed information not available!")

    dsa_size_exp: int | None = compressed_meta.get('size_bytes', None)
    if not isinstance(dsa_size_exp, int):
        sys.exit("Invalid size!")

    # Compare DSA file size with the expected size
    if dsa_size != dsa_size_exp:
        sys.exit(f"Compressed size mismatch ({dsa_size}, expected {dsa_size_exp})!")

    if checksum:
        raise NotImplementedError


if __name__ == '__main__':
    p = ArgumentParser(description="Validate DSA table against truncation.")
    p.add_argument('dsa', help=f"Path of {DSA_FILE_NAME} file")
    p.add_argument('report', help=f"Path of {REPORT_FILE_NAME} file")
    p.add_argument('--checksum', action='store_true')
    args = p.parse_args()

    main(args.dsa, args.report, checksum=args.checksum)
