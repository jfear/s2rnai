#!/usr/bin/env python
"""Adjust DRSC gene coverage.

The presence of DRSC RNAi reagents leads to a bias in counts. To figure out how
much bias there is, this script counts before and after removing the DRSC
region from the GTF.
"""
import os
import sys
from textwrap import dedent
import collections
import yaml
import argparse
from argparse import RawDescriptionHelpFormatter as Raw

import pandas as pd

import pybedtools
import HTSeq


def arguments():
    """Pulls in command line arguments."""

    DESCRIPTION = dedent("""\
    We have found that DRSC RNAi reagents were disproportionally showing up in
    the RNA-seq results. To figure out how big of a problem this is we count
    the number of reads that align to the target gene, the number of reads that
    align to the target gene after removing the reagent region, and the number
    of reads that align to the reagent region. This script also outputs a list
    of reads that align to the reagent in case we want to do some more
    downstream analysis.

    This script creates two files:
    * <ODIR>/*.bam.drsc.counts contains the counts to gene, gene w/o drsc region, and drsc region.
    * <ODIR>/*.bam.drsc.reads.fq contains the reads that map to the drsc region.
    """)

    parser = argparse.ArgumentParser(description=DESCRIPTION, formatter_class=Raw)

    parser.add_argument("--SRR", dest="srr", action='store', required=True,
                        help="Give the SRR identifier.")

    parser.add_argument("--BAM", dest="bam", type=list, action='store', required=True,
                        help="A list of BAM files to analyze")

    parser.add_argument("--out", dest="odir", action='store', required=True,
                        help="The output direcotry.")

    parser.add_argument("--config", dest="config", action='store', required=False, default='../config/config.yml',
                        help="lcdb-wf styled config. [default=../config/config.yml]")

    parser.add_argument("--sampletable", dest="stable", action='store', required=False, default='../config/sampletable.tsv',
                        help="lcdb-wf styled sampletable. [default=../config/sampletable.tsv]")

    parser.add_argument("--drsc_coords", dest="coords", action='store', required=False, default='../output/drsc_coordinates.bed',
                        help="A BED file with DRSC coordinates. [default=../output/drsc_coordinates.bed]")

    args = parser.parse_args()

    return args


def get_drsc_bed(bedName, drsc):
    """Import DRSC metadata.

    Returns
    -------
    pybedtools.BedTool: With coordinates for the current DRSC.
    """
    bed = pd.read_table(bedName, header=None)
    bed.columns = 'chrom', 'start', 'end', 'rname', 'score', 'strand'
    bed.chrom = 'chr' + bed.chrom
    mask = bed.rname == drsc
    return pybedtools.BedTool(bed[mask].values.tolist())


def get_gtf_name(config):
    """Return filename of GTF."""
    assembly = config['assembly']
    tag = config['aligner']['tag']
    return os.path.join(os.environ['REFERENCES_DIR'],
                        assembly, tag, 'gtf',
                        '{assembly}_{tag}.gtf'.format( assembly=assembly, tag=tag))


def featuretype_filter(feature, featuretype, fbgn):
    if (feature[2] == featuretype) & (feature.attrs['gene_id'] == fbgn):
        return True
    return False


def filter_gtf(config, fbgn):
    """Filter out exons for the current gene.

    Returns
    -------
    pybedtools.BedTool: With only exons for the current gene.
    """
    fname = get_gtf_name(config)
    with open(fname) as fh:
        gtf = pybedtools.BedTool([x for x in fh.readlines() if x != '\n'])

    return pybedtools.BedTool(gtf.filter(featuretype_filter,
                                         'exon', fbgn).saveas().fn)


def build_interval_from_file(fname):
    """Create HTSeq array of intervals.

    Returns
    -------
    HTSeq.GenomicArrayOfSets:
        An HTSeq array of intervals for exons.
    """
    fh = HTSeq.GFF_Reader(fname)
    interval = HTSeq.GenomicArrayOfSets("auto", stranded=True)
    for feature in fh:
        if (feature.type == 'exon'):
            interval[feature.iv] += feature.attr['gene_id']

    return interval


def build_interval_from_list(bed):
    """Create HTSeq array of intervals given a list of features.

    Returns
    -------
    HTSeq.GenomicArrayOfSets:
        An HTSeq array of intervals.
    """
    interval = HTSeq.GenomicArrayOfSets("auto", stranded=True)
    for feature in bed:
        curr_interval = HTSeq.GenomicInterval(
            chrom=feature.chrom,
            start=feature.start,
            end=feature.end,
            strand=feature.strand
        )
        interval[curr_interval] += feature.name

    return interval


def invert_strand(iv):
    iv2 = iv.copy()
    if iv2.strand == "+":
        iv2.strand = "-"
    elif iv2.strand == "-":
        iv2.strand = "+"
    else:
        raise ValueError("Illegal strand")
    return iv2


def make_fq(read):
    return dedent(""">{read_id}
                {read_seq}
                +
                {read_qual}""".format(read_id=read.name, read_seq=read.seq.decode('UTF-8'), read_qual=read.qualstr.decode('UTF-8')))


def count_algn(gene, sub, drsc, bam):
    """Count reads aligning to gene w/ drsc, gene w/o drsc, or drsc.

    Returns
    -------
    dict:
        A dicitonary with counts.
    """
    counter = {
        'gene_count': 0,
        'sub_count': 0,
        'drsc_count': 0
    }
    reads = []

    fh = HTSeq.BAM_Reader(bam)
    for algn in fh:
        # skip reads that are on different chromosomes
        if algn.iv.chrom not in list(gene.chrom_vectors.keys()):
            continue

        # Skip unaligned reads
        if not algn.aligned:
            continue

        # Count on full gene
        gene_ids = set()
        for iv, val in gene[invert_strand(algn.iv)].steps():
            gene_ids |= val

        if len(gene_ids) == 1:
            gene_id = list(gene_ids)[0]
            counter['gene_count'] += 1

        # Count on drsc subtracted gene
        gene_ids = set()
        for iv, val in sub[invert_strand(algn.iv)].steps():
            gene_ids |= val

        if len(gene_ids) == 1:
            gene_id = list(gene_ids)[0]
            counter['sub_count'] += 1

        # Count on drsc
        gene_ids = set()
        for iv, val in drsc[invert_strand(algn.iv)].steps():
            gene_ids |= val

        if len(gene_ids) == 1:
            gene_id = list(gene_ids)[0]
            counter['drsc_count'] += 1
            reads.append(make_fq(algn.read))

    return counter, reads


def get_bed_len(bed):
    """Calcuates the total length of a gene.

    Returns
    -------
    int:
        Length of a gene.
    """
    length = 0
    merged = bed.sort().merge()
    for x in merged:
        length += x.length
    return length


def main():
    # Import commandline arguments.
    args = arguments()

    # import config
    with open(args.config) as fh:
        config = yaml.load(fh)

    # import sample table information
    stable = pd.read_table(args.stable)
    srr = args.srr.strip()
    mask = stable.samplename == srr
    drsc = stable[mask].drsc.iloc[0]
    fbgn = stable[mask].target_FBgn.iloc[0]

    # Import drsc
    drsc_bed = get_drsc_bed(args.coords, drsc)

    # Import GTF
    gene_bed = filter_gtf(config, fbgn)

    # Subtract DRSC from GTF
    gene_sub = gene_bed.subtract(drsc_bed)

    # Build HTSeq Intervals
    gene_interval = build_interval_from_file(gene_bed.fn)
    gene_sub_interval = build_interval_from_file(gene_sub.fn)
    drsc_interval = build_interval_from_list(drsc_bed)

    # Count alignments
    for bam in args.bam:
        fname = os.path.basename(bam)
        counts, reads = count_algn(gene_interval, gene_sub_interval, drsc_interval, args.bam)

        # output list of reads that aligned to DRSC
        readName = os.path.join(args.odir, fname + '.drsc.reads.fq')
        with open(readName, 'w') as fh:
            fh.write('\n'.join(reads))

        # Add metadata
        counts['gene'] = fbgn
        counts['drsc'] = drsc

        # Get lengths
        counts['gene_length'] = get_bed_len(gene_bed)
        counts['sub_length'] = get_bed_len(gene_sub)
        counts['drsc_length'] = get_bed_len(drsc_bed)

        # Make DataFrame
        df = pd.DataFrame(counts, index=[srr])
        df.index.name='srr'

        # Output to standard out
        countName = os.path.join(args.odir, fname + '.drsc.counts')
        df.reset_index().to_csv(countName, sep="\t", index=False)


if __name__ == '__main__':
    main()
