#!/usr/bin/env bash

# Validate DSA output directory

set -e

DSA_DIR=$1
SCRIPTS_DIR=$(dirname "$0")

"${SCRIPTS_DIR}/validate_dsa.py" --checksum "${DSA_DIR}/dsa.bed.gz" "${DSA_DIR}/report.json"

set +e
