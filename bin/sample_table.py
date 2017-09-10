#!/usr/bin/env python
"""This script creates the sample table for the S2R+ RNAi study by parsing GEO records.

The best source of information is the published GEO record (GSE81221). Using
these records I parse out all of the sample information in order to create a
sample table. This script also updates FBgn and gene symbols based on the
annotation pointed to in the config file.

Returns
-------
TSV with sample information called: ../config/sample_table.tsv
"""
oname = '../config/sampletable.tsv'
GSE = "GSE81221"
CONFIG = '../config/config.yml'
EMAIL = 'justin.fear@nih.gov'

import os
import sys
import re
from tempfile import TemporaryDirectory
import yaml
from xml.etree import ElementTree

import pandas as pd
import numpy as np
from Bio import Entrez
Entrez.email = EMAIL

import GEOparse

sys.path.insert(0, '../lib/python')
from s2rnai.logger import logger

def srr_iter(srx):
    """Create generator to return SRR given an SRX."""

    res = Entrez.efetch(db='sra', id=srx)
    xml = res.read()
    root = ElementTree.fromstring(xml)
    for run in root.iter('RUN'):
        yield run.get('accession')

    res.close()


# Import config set up references
logger.info('Loading config: {}'.format(CONFIG))
with open(CONFIG) as fh:
    config = yaml.load(fh)

assembly = config['assembly']
tag = config['aligner']['tag']
REF = os.path.join(os.environ['REFERENCES_DIR'], assembly, tag)

# load flybase annotations
FB_ANNO = os.path.join(REF, 'fb_annotation/dmel_{}.fb_annotation'.format(tag))
logger.info('Loading FlyBase annotation file: {}'.format(FB_ANNO))
fb = pd.read_table(FB_ANNO)[['primary_FBgn', 'gene_symbol', 'secondary_FBgn']]
fb.rename(
    columns={
       'primary_FBgn': 'FBgn',
       'gene_symbol': 'symbol'
    },
    inplace=True
)

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
        for f2 in fbgn2.strip().split(','):
            fbgns[f2] = fbgn


# Build sample table using GEO entry
## Query GEO
logger.info('Querying GEO for {}'.format(GSE))
tmpDir = TemporaryDirectory()
gse = GEOparse.get_GEO(GSE, destdir=tmpDir.name, silent=True)

## Pull out sample attributes and build data frame
attributes = []
for gsm, dat in gse.gsms.items():
    try:
        attrs = re.match(r'^.*_(?P<fbgn>FBgn(\d+|NA))_(?P<symbol>.*?)_.*(?P<drsc>DRSC(\d+|NA))_replicate(?P<rep>\d)(\s\[Plate(?P<plate_id>\d+)-\d_F3\]|$)', dat.metadata['title'][0]).groupdict()
    except AttributeError:
        print(gsm, dat.metadata['title'])

    attrs['GEO'] = gsm
    attrs.update(re.match(r'.*_DRSC_Plate(?P<plate_id>\d+)-\d_(?P<well_id>\w\d+)_.*',
                          dat.metadata['supplementary_file_2'][0]).groupdict())

    attrs.update(re.match(r'(?P<plate_row>\w)(?P<plate_column>\d+)', attrs['well_id']).groupdict())

    for x in dat.metadata['relation']:
        k, v = re.match(r'(\w+):.*[\/=](\w+\d+)$', x).groups()
        attrs[k] = v

    # Get SRRs
    for srr in srr_iter(attrs['SRA']):
        attrs['samplename'] = srr
        attributes.append(attrs)

df = pd.DataFrame(attributes)

df.rename(columns={'SRA': 'SRX'}, inplace=True)
df.set_index(['samplename', 'SRX'], inplace=True)

## Grab useful columns and reorder
cols = [
    'BioSample', 'GEO', 'drsc', 'fbgn', 'symbol', 'rep',
    'plate_id', 'well_id', 'plate_row', 'plate_column'
]
df = df[cols]

## Reformat gene symbols
df.symbol = df.symbol.str.replace('[', '(').str.replace(']', ')')

# Map FBgn and symbols to current version
fbgns['FBgnNA'] = 'FBgnNA'
genes['FBgnNA'] = 'LacZ'
df['curr_fbgn'] = df.apply(lambda x: fbgns[x.fbgn], axis=1)
df['curr_symbol'] = df.apply(lambda x: genes[fbgns[x.fbgn]], axis=1)

# Clean up table
cleaned = df.drop(['fbgn', 'symbol'], axis=1).rename(columns={'curr_fbgn': 'target_FBgn', 'curr_symbol': 'target_symbol'})

## Add a column for DRSC replicate (based on DRSC sort order by FBgn)
drsc = cleaned[['drsc', 'target_FBgn']].reset_index(drop=True).drop_duplicates().sort_values('drsc')

drscs = []
for g, grp in drsc.groupby('target_FBgn'):
    new = grp.copy()
    new['drsc_rep'] = range(1, grp.shape[0] + 1)
    drscs.append(new)

drsc = pd.concat(drscs, ignore_index=True)[['drsc', 'drsc_rep']]

## Merge on to cleaned dataset
cleaned = cleaned.reset_index().merge(drsc, left_on='drsc', right_on='drsc')

## Reorder columns
cols = [
    'samplename', 'SRX', 'BioSample', 'GEO',  'drsc',  'target_FBgn',  'target_symbol',
    'drsc_rep', 'rep',  'plate_id',  'well_id', 'plate_row',  'plate_column',
]

logger.info('Writing out sample table: {}'.format(oname))
cleaned[cols].to_csv(oname, sep='\t', index=False)

