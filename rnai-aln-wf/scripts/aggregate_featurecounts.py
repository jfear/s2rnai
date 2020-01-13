"""Aggregates feature counts tables.

Aggregate feature counts tables into a single file. Feature counts can take a
list of BAMs and aggregates them all, but if you have too many files linux
cannot handle the line length. Instead, I am aggregating them here.

"""
import os
from pathlib import Path

import pandas as pd


def main():
    pd.concat(map(read_featurecounts, snakemake.input), axis=1, sort=True).to_csv(
        snakemake.output[0], sep="\t"
    )


def read_featurecounts(file_name):
    sample_name = Path(file_name).parent.name
    return (
        pd.read_csv(file_name, sep="\t", comment="#", index_col="Geneid")
        .iloc[:, -1]
        .rename(sample_name)
        .rename_axis("FBgn")
    )


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from s2rnai.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="rnai-aln-wf",
            input=[
                "../../output/rnai-aln-wf/rnaseq_samples/SRR3486891/SRR3486891.featurecounts.txt",
                "../../output/rnai-aln-wf/rnaseq_samples/SRR3486898/SRR3486898.featurecounts.txt",
            ],
        )

    main()
