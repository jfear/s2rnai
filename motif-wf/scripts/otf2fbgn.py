import os
from urllib.request import urlopen
import re

import pandas as pd

from s2rnai.meme import memeFile


def main():
    onTheFlyTFS = memeFile(snakemake.input[0])

    results = []
    for key, values in onTheFlyTFS.items():
        for value in values:
            # grab link to protein page from url provided
            with urlopen(value.url) as fh:
                page = fh.read().decode("UTF-8")
                name = re.findall(r"protein_entry.php\?protein_ID=(.*?)\'", page)[0]
                ID = re.findall(r"ID: (OTF\d+\.\d+)", page)[0]
            URL = (
                "https://bhapp.c2b2.columbia.edu/OnTheFly/cgi-bin/protein_entry.php?protein_ID={0}"
            )
            try:
                with urlopen(URL.format(name)) as fh:
                    fbgn = re.findall(r"FBgn\d+", fh.read().decode("utf-8"))[0]
                results.append((key, ID, value.name, fbgn))
            except:
                pass

    df = pd.DataFrame(results, columns=["motif_id", "name", "target_protein", "FBgn"])
    df.to_csv(snakemake.output[0], index=False, sep="\t")


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from s2rnai.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="motif-wf", input="../output/motif-wf/meme/OnTheFly_2014_Drosophila.meme"
        )

    main()
