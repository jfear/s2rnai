"""Build supplemental table describing public data used in paper.

I want to provide a supplemental table with all accession information for the
various public datasets used in the paper. These include (S2 RNA-Seq, S2
Chip-Seq, S2 RNAi RNA-Seq). I will pull metadata from my local SRA database
created with `sramongo`.

NOTE: I am using an environment with pymongo installed so I can access the
SRA database.

"""
# %%
import pandas as pd
from pymongo import MongoClient
from pymongo.collection import Collection

# %%
# Build pandas table with useful SRA metadata indexed by SRX.
def merge_srrs(sr: pd.Series):
    return sr.map(lambda x: "|".join(x))


def wrangle_attrs(sr: pd.Series):
    def pretty(attrs: dict):
        return "\n".join([f"{attr['name']}: {attr['value']}" for attr in attrs])

    return sr.map(pretty)


def wrangle_pmid(sr: pd.Series):
    def pretty(papers: dict):
        return ";".join({str(paper["accn"]) for paper in papers})

    return sr.map(pretty)


def wrangle_citation(sr: pd.Series):
    def pretty(papers: dict):
        return ";".join({paper["citation"] for paper in papers})

    return sr.map(pretty)


client = MongoClient()
db: Collection = client["sramongo"]["ncbi"]

metadata = (
    pd.DataFrame(
        db.aggregate(
            [
                {"$unwind": {"path": "$runs"}},
                {
                    "$project": {
                        "_id": False,
                        "SRX": "$srx",
                        "SRR": "$runs.srr",
                        "BioProject": "$BioProject.accn",
                        "BioSample": "$BioSample.accn",
                        "GEO": "$study.geo",
                        "sample_title": "$sample.title",
                        "sample_attributes": "$sample.attributes",
                        "PMID": "$papers",
                        "citation": "$papers",
                    }
                },
            ]
        )
    )
    .assign(sample_attributes=lambda x: wrangle_attrs(x.sample_attributes))
    .assign(PMID=lambda x: wrangle_pmid(x.PMID))
    .assign(citation=lambda x: wrangle_citation(x.citation))
)

client.close()

# %%
# S2 RNAi Experiment
rnai = (
    pd.read_table("../rnai-aln-wf/config/sampletable.tsv")
    .assign(experiment="S2 RNAi RNA-Seq")
    .rename(
        columns={
            "Run": "SRR",
            "drsc": "rnai_drsc_id",
            "target_FBgn": "rnai_target_FBgn",
            "target_symbol": "rnai_target_symbol",
        }
    )
    .drop(
        [
            "samplename",
            "BioSample",
            "GEO",
            "rep",
            "plate_id",
            "well_id",
            "plate_row",
            "plate_column",
        ],
        axis=1,
    )
)

# %%
# S2 RNA-Seq from SRA
rnaseq = (
    pd.read_csv("../data/sra_s2_rnaseq.txt", header=None, names=["SRX"])
    .assign(experiment="S2 RNA-Seq")
    .merge(metadata[["SRX", "SRR"]])
)


# %%
# S2 Chip-seq from SRA
chipseq = (
    pd.read_table("../data/sra_s2_chipseq.tsv")
    .rename(columns={"antibody": "chip_seq_antibody"})
    .drop(["replicate"], axis=1)
    .assign(experiment="S2 ChIP-Seq")
)

# %%
# merge everything together
stack = pd.concat([rnaseq, chipseq, rnai], ignore_index=True, sort=False).sort_values(
    ["SRX", "SRR"]
)
df = metadata.merge(stack).set_index(
    ["SRX", "SRR", "BioProject", "BioSample", "GEO", "PMID", "citation"]
)

# %%

df.to_csv("../output/notebook/sra_identifiers.tsv", sep="\t")


# %%
