import csv

from collections import defaultdict
from functools import cache

from helpers import __version__

def get_version() -> str:
    return __version__

@cache
def parse_samplesheet(samplesheet: str) -> dict[dict[str]]:
    """
    Parses a given tab-separated samplesheet into a dictionary.
    {"PD1111a": {"1": "/path/to/file.bam", "2": "/path/to/file.bam"},
    ...
    }
    """
    sample_list = defaultdict(lambda: defaultdict(str))

    with open(samplesheet) as samples:
        for sample in csv.DictReader(samples, delimiter="\t"):
            sample_list[sample["sample_id"]][sample["lane"]] = sample["path"]

    return sample_list
