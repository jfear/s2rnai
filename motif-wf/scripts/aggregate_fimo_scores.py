import os
from collections import ChainMap
from functools import partial

import pandas as pd


def main():
    flybase_tfs = get_tfs()
    fbgn_mapper = create_fbgn_mapper()

    motif_scores_all_db = pd.concat(
        [
            aggregate_motif_score(file_name, fbgn_mapper, flybase_tfs)
            for file_name in snakemake.input.files
        ],
        sort=True,
    )

    (
        # For each TF-Target, pull out the best score across databases
        motif_scores_all_db.groupby(["transcription_factor", "target_FBgn"])
        .motif_score.max()
        .unstack(level=0)
        .to_csv(snakemake.output[0], sep="\t")
    )


def aggregate_motif_score(file_name, fbgn_mapper, tfs):
    df = pd.read_table(file_name)
    df.transcription_factor = df.transcription_factor.map(fbgn_mapper)
    df_tfs_only = df[df.transcription_factor.isin(tfs)]
    single_gene_motifs = drop_multi_gene_motifs(df_tfs_only)
    tf_gene_score = single_gene_motifs.groupby(["transcription_factor", "target_FBgn"]).agg(
        dict(motif_score="sum", database="first")
    )
    return tf_gene_score


def drop_multi_gene_motifs(df):
    number_target_genes = df.groupby(["chrom", "start", "stop", "transcription_factor"]).size()
    single_gene_sites = (
        number_target_genes[number_target_genes == 1]
        .rename("count")
        .to_frame()
        .reset_index()
        .drop("count", axis=1)
    )
    return df.merge(
        single_gene_sites, on=["chrom", "start", "stop", "transcription_factor"], how="inner"
    )


def create_fbgn_mapper():
    fbgn2fbgn = create_fbgn2fbgn()
    symbol2fbgn = create_symbol2fbgn()
    otf2fbgn = create_otf2fbgn()
    return dict(ChainMap({}, fbgn2fbgn, symbol2fbgn, otf2fbgn))


def create_fbgn2fbgn():
    mapper = {}
    for _, row in pd.read_table(snakemake.input.annotation).iterrows():
        mapper[row.primary_FBgn] = row.primary_FBgn
        if isinstance(row.secondary_FBgn, str):
            for secondary in row.secondary_FBgn.split(","):
                mapper[secondary] = row.primary_FBgn
    return mapper


def create_symbol2fbgn():
    return (
        pd.read_table(snakemake.input.annotation, usecols=["gene_symbol", "primary_FBgn"])
        .set_index("gene_symbol")
        .squeeze()
        .to_dict()
    )


def create_otf2fbgn():
    return (
        pd.read_table(snakemake.input.otf, usecols=["name", "FBgn"])
        .set_index("name")
        .squeeze()
        .to_dict()
    )


def get_tfs():
    return pd.read_table(snakemake.input.tfs, comment="#", header=None).values.flatten()


def agg_motifs(df: pd.DataFrame):
    df.groupby(["transcription_factor", "chrom", "start", "stop"]).apply(
        lambda x: x.database.str.cat(sep="|")
    )


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from s2rnai.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="motif-wf",
            input=dict(
                files=[
                    "../../output/motif-wf/dmmpmm2009/fimo_gene_intersection.tsv",
                    "../../output/motif-wf/OnTheFly_2014_Drosophila/fimo_gene_intersection.tsv",
                ],
                annotation="../../lcdb-references/dmel/r6-26/fb_annotation/dmel_r6-26.fb_annotation",
                otf="../../output/motif-wf/onTheFlyMap.tsv",
                tfs="../../data/external/FlyBase/transcription_factors.tsv",
            ),
        )

    main()
