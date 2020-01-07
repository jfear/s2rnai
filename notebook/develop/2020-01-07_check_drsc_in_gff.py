"""Quick check of drsc

Just checking that the drsc are the same in the GFF and the table that I used
previously.

They are!!
"""
# %%
import os
import pandas as pd

from s2rnai.io import GffRow

# %%
try:
    os.chdir(os.path.join(os.getcwd(), 'notebook'))
    print(os.getcwd())
except:
    pass


# %%
samples = pd.read_csv("../rnai-aln-wf/config/sampletable.tsv", sep='\t')

# %%
drsc = samples.drsc.unique()

# %%
gff_drsc = []
with open("/home/fearjm/Downloads/dmel-all-no-analysis-r6.26.gff") as fh:
    for row in fh.readlines():
        if "DRSC" in row:
            gff = GffRow(row)
            if gff.source == "DRSC_dsRNA":
                gff_drsc.append(gff.parsed_attributes["Name"])

# %%
for dd in drsc: 
    if dd in gff_drsc:
        continue
    print(dd)


# %%
