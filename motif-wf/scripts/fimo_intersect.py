"""Intersects Fimo results with genes.

Takes Fimo mapping locations and intersects them with genes (plus 1kb
upstream). Saves results as a TSV.
"""
import os
from typing import Generator

import pandas as pd
import gffutils
from gffutils.pybedtools_integration import to_bedtool
import pybedtools

BED_HEADER = "chrom start stop name score strand".split()


def main():
    genes_plus_1k = get_genes_plus_1k(snakemake.input.db, snakemake.input.chromsizes)
    fimo_bed = get_fimo_as_bed(snakemake.input.motif)
    intersection = fimo_bed.intersect(genes_plus_1k, wo=True)
    df = wo_output_to_dataframe(intersection)

    df.to_csv(snakemake.output[0], sep="\t", index=False)


def get_genes_plus_1k(db_file, chromsizes_file) -> pybedtools.BedTool:
    db = gffutils.FeatureDB(db_file)
    genes = db.features_of_type("gene")

    def gtf_to_bed(gene: gffutils.Feature):
        return db.bed12(gene, block_featuretype="gene", name_field="gene_id")

    genes_bed = "\n".join(map(gtf_to_bed, genes))
    return pybedtools.BedTool(genes_bed, from_string=True).slop(
        l=1_000, r=0, s=True, g=chromsizes_file
    )


def get_fimo_as_bed(file_name) -> pybedtools.BedTool:
    df = (
        pd.read_table(file_name, comment="#")
        .assign(chrom=lambda x: x.sequence_name)
        .assign(name=lambda x: x.motif_id)
        .reindex(columns=BED_HEADER)
    )
    return pybedtools.BedTool(df.to_csv(sep="\t", header=None, index=False), from_string=True)


def wo_output_to_dataframe(data: pybedtools.BedTool) -> pd.DataFrame:
    return (
        data.to_dataframe(disable_auto_names=True, header=None)
        .assign(chrom=lambda x: x[0])
        .assign(start=lambda x: x[1])
        .assign(stop=lambda x: x[2])
        .assign(transcription_factor=lambda x: x[3])
        .assign(motif_score=lambda x: x[4])
        .assign(motif_strand=lambda x: x[5])
        .assign(target_FBgn=lambda x: x[9])
        .assign(target_strand=lambda x: x[11])
        .assign(database=snakemake.wildcards.db_name)
        .reindex(
            columns=[
                "chrom",
                "start",
                "stop",
                "transcription_factor",
                "motif_score",
                "motif_strand",
                "target_FBgn",
                "target_strand",
                "database",
            ]
        )
    )


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from s2rnai.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="motif-wf",
            input=dict(
                motif="../../output/motif-wf/dmmpmm2009/fimo.tsv",
                db="../../lcdb-references/dmel/r6-26/gtf/dmel_r6-26.gtf.db",
                chromsizes="../../lcdb-references/dmel/r6-26/fasta/dmel_r6-26.chromsizes",
            ),
            wildcards=dict(db_name="dmmpmm2009"),
        )

    main()
