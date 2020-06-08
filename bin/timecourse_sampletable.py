#!/usr/bin/env python
"""This script creates the sample table for the S2R+ RNAi Time Course study by parsing GEO records.

The best source of information is the published GEO record (GSE89753). Using
these records I parse out all of the sample information in order to create a
sample table. This script also updates FBgn and gene symbols based on the
annotation pointed to in the config file.

Returns
-------
TSV with sample information called: ../config/sampletable_timecourse.tsv
"""
import os
import re
import sys
from tempfile import TemporaryDirectory
from xml.etree import ElementTree
from pathlib import Path

import GEOparse
import numpy as np
import pandas as pd
import yaml
from Bio import Entrez
from dotenv import load_dotenv

from s2rnai.logger import logger

load_dotenv()
Entrez.email = os.getenv("EMAIL")

PROJECT_DIR = Path(__file__).absolute().parents[1]
CONFIG = yaml.safe_load((PROJECT_DIR / "config/common.yaml").open())
GSE = "GSE89753"


def main():
    sampletable = (
        download_geo_metadata()
        .pipe(update_fbgns)
        .pipe(update_gene_symbols)
        .pipe(reorder_columns)
    )

    out_file = PROJECT_DIR / "timecourse-aln-wf/config/sampletable.tsv"
    logger.info(f"Writing out sample table: {out_file}")
    sampletable.to_csv(out_file, sep="\t", index=False)


def download_geo_metadata() -> pd.DataFrame:
    logger.info("Querying GEO for {}".format(GSE))
    tmpDir = TemporaryDirectory()
    gse = GEOparse.get_GEO(GSE, destdir=tmpDir.name, silent=True)

    ## Pull out sample attributes and build data frame
    attributes = []
    for gsm, dat in gse.gsms.items():
        try:
            # Parse sample title
            attrs = re.match(
                r"^DRSC_(?P<plate_id>plate\d)_(?P<time_point>day\d)_(?P<plate_row>[A-H])(?P<plate_column>\d+)_(?P<drsc>DRSC(\d+|NA))_(?P<fbgn>(FBgn\d+|NA))_(?P<symbol>.*?)$",
                dat.metadata["title"][0],
            ).groupdict()

            attrs["GEO"] = gsm

            # Get SRX accession from linkout
            for x in dat.metadata["relation"]:
                match = re.match(r"(\w+):.*[\/=](\w+\d+)$", x)
                if match:
                    k, v = match.groups()
                    attrs[k] = v

            # Expand out SRR accessions
            for srr in get_srrs(attrs["SRA"]):
                attrs["samplename"] = srr
                attrs["Run"] = srr
                attributes.append(attrs)

        except AttributeError:
            print(gsm, dat.metadata["title"])

    df = pd.DataFrame(attributes).rename(columns={"SRA": "SRX"}).set_index(["samplename", "SRX"])
    return df


def get_srrs(srx):
    """Generator returns all SRRs for an SRX"""
    res = Entrez.efetch(db="sra", id=srx)
    xml = res.read()
    root = ElementTree.fromstring(xml)
    for run in root.iter("RUN"):
        yield run.get("accession")
    res.close()


def update_fbgns(df: pd.DataFrame) -> pd.DataFrame:
    mapper = _fbgn_mapper()
    df["target_FBgn"] = df.fbgn.map(mapper)
    return df.drop("fbgn", axis=1)


def update_gene_symbols(df: pd.DataFrame) -> pd.DataFrame:
    mapper = _gene_symbol_mapper()
    df["target_symbol"] = df.symbol.map(mapper)
    return df.drop("symbol", axis=1)


def reorder_columns(df: pd.DataFrame) -> pd.DataFrame:
    cols = [
        "samplename",
        "Run",
        "SRX",
        "BioSample",
        "GEO",
        "drsc",
        "time_point",
        "target_FBgn",
        "target_symbol",
        "plate_id",
        "plate_row",
        "plate_column",
    ]
    return df.reset_index()[cols]


def _fbgn_mapper() -> dict:
    primary_to_secondary = (
        pd.read_table(
            PROJECT_DIR / f"lcdb-references/dmel/{CONFIG['tag']}/fb_annotation/dmel_{CONFIG['tag']}.fb_annotation"
        )
        .set_index("primary_FBgn")
        .secondary_FBgn
    )

    fbgn_mapper = {primary: primary for primary in primary_to_secondary.index}
    for primary, secondaries in primary_to_secondary.dropna().items():
        for secondary in secondaries.split(","):
            fbgn_mapper[secondary] = primary

    fbgn_mapper["NA"] = "FBgnNA"
    return fbgn_mapper


def _gene_symbol_mapper() -> dict:
    primary_to_secondary = (
        pd.read_table(
            PROJECT_DIR / f"lcdb-references/dmel/{CONFIG['tag']}/fb_synonym/dmel_{CONFIG['tag']}.fb_synonym",
            encoding="utf-8",
        )
        .query("organism_abbreviation == 'Dmel'")
        .set_index("current_symbol")
        .loc[:, "symbol_synonym(s)"]
    )

    symbol_mapper = {primary: primary for primary in primary_to_secondary.index}
    for primary, secondaries in primary_to_secondary.dropna().items():
        for secondary in secondaries.split(","):
            symbol_mapper[secondary] = primary

    symbol_mapper["NA"] = "LacZ"
    return symbol_mapper


def get_current_flybase_annotations():
    # Import config set up references
    logger.info("Loading config: {}".format(CONFIG))
    with open(CONFIG) as fh:
        config = yaml.load(fh)

    assembly = config["assembly"]
    tag = config["aligner"]["tag"]
    REF = os.path.join(os.environ["REFERENCES_DIR"], assembly, tag)

    # load flybase annotations
    FB_ANNO = os.path.join(REF, "fb_annotation/dmel_{}.fb_annotation".format(tag))
    logger.info("Loading FlyBase annotation file: {}".format(FB_ANNO))
    fb = pd.read_table(FB_ANNO)[["primary_FBgn", "gene_symbol", "secondary_FBgn"]]
    fb.rename(columns={"primary_FBgn": "FBgn", "gene_symbol": "symbol"}, inplace=True)

    ## Make map of old fbgn to current fbgn and current fbgn to current symbol
    fbgns = {}
    genes = {}
    for i, record in fb.iterrows():
        fbgn = record.FBgn
        symbol = record.symbol
        fbgn2 = record.secondary_FBgn

        fbgns[fbgn] = fbgn
        genes[fbgn] = symbol

        if isinstance(fbgn2, str):
            for f2 in fbgn2.strip().split(","):
                fbgns[f2] = fbgn


if __name__ == "__main__":
    main()
