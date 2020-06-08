#!/usr/bin/env python
"""Adjust DRSC gene coverage.

The presence of DRSC RNAi reagents leads to a bias in counts. To figure out how
much bias there is, this script counts before and after removing the DRSC
region from the GTF.
"""
import os
from collections import Counter
from typing import Tuple, List

import pandas as pd

import pybedtools
import HTSeq


def main():
    # Import drsc_bed and dmel GTF
    drsc_bed = get_drsc_bed(snakemake.input.drsc_bed, snakemake.params.drsc)
    dmel_bed = filter_gtf(snakemake.input.gtf, snakemake.params.fbgn)

    # Subtract DRSC from GTF
    dmel_minus_drsc = dmel_bed.subtract(drsc_bed)

    # Build HTSeq Intervals
    dmel_interval = build_interval_from_file(dmel_bed.fn)
    dmel_minus_drsc_interval = build_interval_from_file(dmel_minus_drsc.fn)
    drsc_interval = build_interval_from_list(drsc_bed)

    # Count alignments
    counts, reads = count_alignment(
        dmel_interval, dmel_minus_drsc_interval, drsc_interval, snakemake.input.bam
    )

    # output list of reads that aligned to DRSC
    with open(snakemake.output.fq, "w") as fh:
        fh.write("\n".join(reads))

    # Output counts table
    df = (
        pd.DataFrame(counts, index=pd.Index([snakemake.wildcards.sample,], name="srr"))
        .assign(FBgn=snakemake.params.fbgn)
        .assign(drsc=snakemake.params.drsc)
        .assign(gene_length=get_bed_len(dmel_bed))
        .assign(gene_minus_drsc_length=get_bed_len(dmel_minus_drsc))
        .assign(drsc_length=get_bed_len(drsc_bed))
        .set_index(["drsc", "FBgn"], append=True)
    )

    df.reset_index().to_csv(snakemake.output.counts, sep="\t", index=False)


def get_drsc_bed(bed_file_name, drsc) -> pybedtools.BedTool:
    """Import DRSC metadata.

    Returns
    -------
    pybedtools.BedTool: With coordinates for the current DRSC.
    """
    bed = pd.read_table(
        bed_file_name, header=None, names=["chrom", "start", "end", "drsc", "score", "strand"]
    )
    bed.chrom = "chr" + bed.chrom
    return pybedtools.BedTool(bed.query(f"drsc == '{drsc}'").values.tolist())


def filter_gtf(gtf_file_name, fbgn) -> pybedtools.BedTool:
    """Filter out exons for the current gene.

    Returns
    -------
    pybedtools.BedTool: With only exons for the current gene.
    """
    with open(gtf_file_name) as fh:
        gtf = pybedtools.BedTool([x for x in fh.readlines() if x != "\n"])

    return pybedtools.BedTool(gtf.filter(featuretype_filter, "exon", fbgn).saveas().fn)


def featuretype_filter(feature, featuretype, fbgn) -> bool:
    if (feature[2] == featuretype) & (feature.attrs["gene_id"] == fbgn):
        return True
    return False


def build_interval_from_file(file_name) -> HTSeq.GenomicArrayOfSets:
    """Create HTSeq array of intervals.

    Returns
    -------
    HTSeq.GenomicArrayOfSets:
        An HTSeq array of intervals for exons.
    """
    fh = HTSeq.GFF_Reader(file_name)
    interval = HTSeq.GenomicArrayOfSets("auto", stranded=True)
    for feature in fh:
        if feature.type == "exon":
            interval[feature.iv] += feature.attr["gene_id"]

    return interval


def build_interval_from_list(bed) -> HTSeq.GenomicArrayOfSets:
    """Create HTSeq array of intervals given a list of features.

    Returns
    -------
    HTSeq.GenomicArrayOfSets:
        An HTSeq array of intervals.
    """
    interval = HTSeq.GenomicArrayOfSets("auto", stranded=True)
    for feature in bed:
        curr_interval = HTSeq.GenomicInterval(
            chrom=feature.chrom, start=feature.start, end=feature.end, strand=feature.strand
        )
        interval[curr_interval] += feature.name

    return interval


def count_alignment(gene, no_drsc, drsc, bam) -> Tuple[Counter, List[str]]:
    """Count aligning reads.

    Here I count how many reads align to the different groups (gene level,
    gene minus drsc, and drsc only). I also look at reads aligning to same
    and opposite strands. These data should map to the opposite strand.

    Returns
    -------
    Counter: the various gene level counts.
    List: List of reads in FASTQ format.

    """
    counter = Counter()
    reads = []

    aln_file = HTSeq.BAM_Reader(bam)
    for read in aln_file:  # type: HTSeq._HTSeq.SAM_Alignment
        # skip reads that are on different chromosomes
        if read.iv.chrom not in list(gene.chrom_vectors.keys()):
            continue

        # Skip unaligned reads
        if not read.aligned:
            continue

        # Count on full gene same strand
        gene_ids = set()
        for iv, val in gene[read.iv].steps():
            gene_ids |= val

        if len(gene_ids) == 1:
            gene_id = list(gene_ids)[0]
            counter["gene_count_same"] += 1

        # Count on full gene opposite strand
        gene_ids = set()
        for iv, val in gene[invert_strand(read.iv)].steps():
            gene_ids |= val

        if len(gene_ids) == 1:
            gene_id = list(gene_ids)[0]
            counter["gene_count_opposite"] += 1

        # Count on drsc subtracted gene same strand
        gene_ids = set()
        for iv, val in no_drsc[read.iv].steps():
            gene_ids |= val

        if len(gene_ids) == 1:
            gene_id = list(gene_ids)[0]
            counter["gene_wo_drsc_count_same"] += 1

        # Count on drsc subtracted gene opposite strand
        gene_ids = set()
        for iv, val in no_drsc[invert_strand(read.iv)].steps():
            gene_ids |= val

        if len(gene_ids) == 1:
            gene_id = list(gene_ids)[0]
            counter["gene_wo_count_opposite"] += 1

        # Count on drsc same strand
        gene_ids = set()
        for iv, val in drsc[read.iv].steps():
            gene_ids |= val

        if len(gene_ids) == 1:
            gene_id = list(gene_ids)[0]
            counter["drsc_count_same"] += 1

        # Count on drsc opposite strand
        gene_ids = set()
        for iv, val in drsc[invert_strand(read.iv)].steps():
            gene_ids |= val

        if len(gene_ids) == 1:
            gene_id = list(gene_ids)[0]
            counter["drsc_count_opposite"] += 1
            reads.append(make_fq(read.read))

    return counter, reads


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
    return f"@{read.name}\n{read.seq.decode('UTF-8')}\n+\n{read.qualstr.decode('UTF-8')}"


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


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from s2rnai.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="rna-aln-wf",
            input=dict(
                bam="../../output/rnai-aln-wf/rnaseq_samples/SRR3486891/SRR3486891.cutadapt.bam",
                drsc_bed="../../output/drsc.bed",
                gtf="../../lcdb-references/dmel/r6-26_and_ercc/gtf/dmel_r6-26_and_ercc.gtf",
            ),
            wildcards=dict(sample="SRR3486891"),
            params=dict(drsc="DRSC07681", fbgn="FBgn0003396"),
            output=dict(counts="", fq=""),
        )

    main()
