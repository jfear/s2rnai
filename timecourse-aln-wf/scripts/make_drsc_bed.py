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
import tempfile
from pathlib import Path
from typing import Tuple, List
from collections import namedtuple

import pandas as pd

from s2rnai.io import GffRow

BedRecord = namedtuple("BedRecord", "chrom,start,end,name,score,strand")

def main():
    gff_file_name = download_gff()
    drsc_bed, gene_strand = split_gff(gff_file_name)
    delete_tempfile(gff_file_name)

    study_drsc_to_fbgn = (
        pd.read_table(snakemake.input[0], usecols=["drsc", "target_FBgn"])
        .drop_duplicates()
        .query("drsc != 'DRSCNA'")
        .set_index("drsc")
        .squeeze()
        .to_dict()
    )

    study_drsc = pull_study_drsc(drsc_bed, study_drsc_to_fbgn.keys())
    updated_drsc = update_drsc_strand(study_drsc, gene_strand, study_drsc_to_fbgn)

    # Write out table
    pd.DataFrame(updated_drsc).to_csv(snakemake.output[0], sep="\t", index=False, header=False)


def download_gff() -> str:
    gff_tempfile = get_tempfile_name()
    response = urllib.request.urlopen(snakemake.params.url)
    with open(gff_tempfile, "w") as file_out:
        file_out.write(gzip.decompress(response.read()).decode())
    return gff_tempfile


def get_tempfile_name() -> str:
    gff_obj = tempfile.NamedTemporaryFile(suffix=".gff", delete=False)
    gff_file_name = gff_obj.name
    gff_obj.close()
    return gff_file_name


def delete_tempfile(file_name):
    Path(file_name).unlink()


def split_gff(gff_file_name) -> Tuple[List[BedRecord], dict]:
    """Pulls out drsc and gene strand info"""
    drsc = missing_drsc()
    gene_strands = {}
    with open(gff_file_name) as fh:
        for row in fh:
            if "DRSC_dsRNA" in row:
                gff = GffRow(row)
                drsc.append(
                    BedRecord(
                        gff.seqid,
                        gff.start,
                        gff.end,
                        gff.parsed_attributes["Name"],
                        gff.score,
                        gff.strand,
                    )
                )
            elif "FlyBase\tgene" in row:
                gff = GffRow(row)
                gene_strands[gff.parsed_attributes["ID"]] = gff.strand
    return drsc, gene_strands


def missing_drsc() -> List[BedRecord]:
    """There are a handful of DRSC missing from the GFF

    Look them up missing ones here:
    https://www.flyrnai.org/up-torr/JBrowse-1.12.3/index.html to get their
    locations.
    """
    return [
        BedRecord("4", "545545", "545879", "DRSC38683", ".", "+"),
        BedRecord("X", "6266464", "6266610", "DRSC40082", ".", "+"),
        BedRecord("2L", "10264095", "10264265", "DRSC40331", ".", "+"),
        BedRecord("3R", "27918291", "27918881", "DRSC40340", ".", "+"),
        BedRecord("X", "7671912", "7672385", "DRSC40438", ".", "+"),
        BedRecord("2L", "16698174", "16698323", "DRSC40694", ".", "+"),
        BedRecord("X", "2448482", "2448642", "DRSC40820", ".", "+"),
        BedRecord("2L", "15116373", "15116962", "DRSC40850", ".", "+"),
        BedRecord("2R", "21143622", "21143953", "DRSC40918", ".", "+"),
        BedRecord("3R", "25025517", "25025716", "DRSC41040", ".", "+"),
        BedRecord("2R", "12633620", "12633819", "DRSC41571", ".", "+"),
        BedRecord("2L", "9762479", "9762549", "DRSC42501", ".", "+"),
        BedRecord("2L", "21286141", "21286386", "DRSC41069", ".", "+"),
        BedRecord("2L", "10388098", "10388500", "DRSC40784", ".", "+"),
        BedRecord("3R", "31728078", "31728479", "DRSC40439", ".", "+"),
        BedRecord("X", "15995488", "15995996", "DRSC41279", ".", "+"),
        BedRecord("2R", "13226618", "13227061", "DRSC40819", ".", "+"),
        BedRecord("3R", "16040121", "16040454", "DRSC40910", ".", "+"),
        BedRecord("3R", "23050074", "23050355", "DRSC40777", ".", "+"),  # Also found at X:11120655..11120957
        BedRecord("3L", "17561059", "17561298", "DRSC40711", ".", "+"),
        BedRecord("X", "15825153", "15825491", "DRSC40479", ".", "+"),
    ]


def pull_study_drsc(drsc_bed: List[BedRecord], drscs_in_study: List[str]) -> List[BedRecord]:
    records_in_study = [
        record
        for record in drsc_bed
        if record.name in drscs_in_study
    ]
    check_all_study_drsc_present(records_in_study, drscs_in_study)
    return records_in_study


def check_all_study_drsc_present(records: List[BedRecord], drscs_in_study: List[str]):
    present_drsc = set([record.name for record in records])
    study_set = set(drscs_in_study)
    difference = study_set.difference(present_drsc)
    if len(difference) > 0:
        raise ValueError("There are missing DRSCs:\n%s" % "\n".join(difference))


def update_drsc_strand(drsc: List[BedRecord], gene_strand: dict, drsc_to_fbgn: dict) -> List[BedRecord]:
    updated_records = []
    for record in drsc:
        drsc_id = record.name
        target_fbgn = drsc_to_fbgn[drsc_id]
        mutable_record = list(record)
        updated_strand = gene_strand[target_fbgn]
        mutable_record[5] = updated_strand
        updated_records.append(BedRecord(*mutable_record))
    return updated_records



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
