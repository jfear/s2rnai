#!/usr/bin/env python
"""This script builds a coordinate set for DRSC reagent locations.

Most DRSC reagents are in the FlyBase GFF. There is a set of 12 reagents that
are not in FlyBase's files, so I had to use the DRSC UP-TORR page:

http://www.flyrnai.org/up-torr/

Returns
-------
BED formatted file ../output/drsc_coordinates.bed

"""
import gzip
import os
import urllib.request
from typing import Generator, Tuple, List

import pandas as pd

from s2rnai.io import GffRow


def main():
    coords = [parse_gff_row(row) for row in download_gff() if "DRSC_dsRNA" in row]
    coords.extend(add_missing())

    # Check for missing DRSCs
    drscs_in_gff = {x[3] for x in coords}
    drscs_in_study = get_drscs_in_study()
    missing = sorted(drscs_in_study.difference(drscs_in_gff))

    if len(missing) != 0:
        raise ValueError("Missing the following DRSCs:\n%s" % "\n".join(missing))

    # Write out table
    pd.DataFrame(coords).to_csv(snakemake.output[0], sep="\t", index=False, header=False)


def download_gff() -> Generator[str, None, None]:
    """Download and decompress GFF from FlyBase."""
    with urllib.request.urlopen(snakemake.params.url) as response:
        for line in gzip.decompress(response.read()).decode().splitlines():
            yield line


def parse_gff_row(row: str) -> Tuple[str]:
    """Parses GFF row into a tuple."""
    gff = GffRow(row)
    return gff.seqid, gff.start, gff.end, gff.parsed_attributes["Name"], gff.score, gff.strand


def get_drscs_in_study() -> set:
    """Pull out DRSCs in the study."""
    return set(pd.read_csv(snakemake.input[0], sep="\t").query('drsc != "DRSCNA"').drsc.unique())


def add_missing() -> List[Tuple[str]]:
    """Add missing drsc coordinates.

    I compared these coordinates to a list of DRSC reagents used and found some to
    be missing. I used the UP-TORR website to look up these reagent locations. I
    add them here.
    """
    return [
        ("4", "545545", "545879", "DRSC38683", ".", "+"),
        ("X", "6266464", "6266610", "DRSC40082", ".", "+"),
        ("2L", "10264095", "10264265", "DRSC40331", ".", "+"),
        ("3R", "27918291", "27918881", "DRSC40340", ".", "+"),
        ("X", "7671912", "7672385", "DRSC40438", ".", "+"),
        ("2L", "16698174", "16698323", "DRSC40694", ".", "+"),
        ("X", "2448482", "2448642", "DRSC40820", ".", "+"),
        ("2L", "15116373", "15116962", "DRSC40850", ".", "+"),
        ("2R", "21143622", "21143953", "DRSC40918", ".", "+"),
        ("3R", "25025517", "25025716", "DRSC41040", ".", "+"),
        ("2R", "12633620", "12633819", "DRSC41571", ".", "+"),
        ("2L", "9762479", "9762549", "DRSC42501", ".", "+"),
    ]


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from s2rnai.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="rnai-aln-wf",
            input="../config/sampletable.tsv",
            params=dict(
                url="ftp://ftp.flybase.net/releases/FB2019_01/dmel_r6.26/gff/dmel-all-no-analysis-r6.26.gff.gz",
            ),
        )

    main()
