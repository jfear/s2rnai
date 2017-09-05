#!/usr/bin/env python
"""This script builds a coordinate set for DRSC reagent locations.

Most DRSC reagents are in the FlyBase GFF. There is a set of 12 reagents that
are not in FlyBase's files, so I had to use the DRSC UP-TORR page:

http://www.flyrnai.org/up-torr/

Returns
-------
BED formatted file ../output/drsc_coordinates.bed

"""
oname = '../output/drsc_coordinates.bed'

import os
import sys
import urllib.request
import gzip
import yaml
import subprocess

import pandas as pd

import pybedtools
from pybedtools.featurefuncs import gff2bed

sys.path.insert(0, '../lib/python')

from s2rnai.logger import logger

if __name__ == '__main__':
    # Import config and set up references
    with open('../config/config.yml') as fh:
        config = yaml.load(fh)

    assembly = config['assembly']
    tag = config['aligner']['tag']
    REF = os.path.join(os.environ['REFERENCES_DIR'], assembly, tag)

    # Download gff file from FlyBase
    url = config['references'][config['assembly']][config['gtf']['tag']]['gtf']['url'].replace('gtf', 'gff')
    fname = os.path.join('../data/external/FlyBase', os.path.basename(url))

    if not os.path.exists(fname):
        logger.info('Createing FlyBase directory.')
        os.makedirs(os.path.dirname(fname), exist_ok=True)
        with open(fname, 'wb') as fh:
            logger.info('Downloading GFF')
            response = urllib.request.urlretrieve(url, fname)

    # Output filename to store drsc gff info
    drsc_fname = fname.replace('.gff.gz', '.drsc.gff')

    # use grep to filter our DRSC, this is the fastest way
    if not os.path.exists(drsc_fname):
        logger.info('Filtering GFF')
        cmd = 'gunzip -c {fname} | grep "DRSC_dsRNA" | grep "RNAi_reagent" > {drsc_fname}'.format(fname=fname, drsc_fname=drsc_fname)
        subprocess.run(cmd, shell=True)

    # Convert DRSC Gff to BED then to dataframe
    logger.info('Creating DRSC BED: {}'.format(oname))
    drsc_gff = pybedtools.BedTool(drsc_fname)
    drsc_bed = drsc_gff.each(gff2bed, name_field='Name').saveas()
    drsc_coords = drsc_bed.to_dataframe()

    """Add missing drsc coordinates.

    I compared these coordinates to a list of DRSC reagents used and found some to
    be missing. I used the UP-TORR website to look up these reagent locations. I
    add them here.
    """
    missing = pd.DataFrame(
        [
            ('2L', '15116372', '15116962', 'DRSC40850', '+'),
            ('2L', '16698173', '16698323', 'DRSC40694', '.', '+'),
            ('2L', '9762478', '9762549', 'DRSC42501', '.', '+'),
            ('2L', '10264094', '10264265', 'DRSC40331', '.', '+'),
            ('2R', '12633619', '12633819', 'DRSC41571', '.', '+'),
            ('2R', '21143621', '21143953', 'DRSC40918', '.', '+'),
            ('3R', '25025516', '25025716', 'DRSC41040', '.', '+'),
            ('3R', '27918290', '27918881', 'DRSC40340', '.', '+'),
            ('4', '545544', '545879', 'DRSC38683', '.', '+'),
            ('X', '7671911', '7672385', 'DRSC40438', '.', '+'),
            ('X', '6266463', '6266610', 'DRSC40082', '.', '+'),
            ('X', '2448481', '2448642', 'DRSC40820', '.', '+')
        ],
        columns=['chrom', 'start', 'end', 'name', 'score', 'strand'])

    # Combine missing to create a large dataframe.
    drsc_plus_missing = pd.concat([drsc_coords, missing])

    # Write out to a bed file.
    drsc_plus_missing.to_csv(oname, sep='\t', na_rep='.', header=False, index=False)
