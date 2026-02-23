#!/usr/bin/env bash

# Validate DSA output directory

set -e

DSA_DIR=$1

./validate_dsa.py "${DSA_DIR}/dsa.bed.gz" "${DSA_DIR}/report.json"

set +e
