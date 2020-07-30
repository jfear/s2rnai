import sys
from pathlib import Path
import csv

import pandas as pd


def main():
    srx2fbgn = pd.read_csv(
        "./config/sampletable.tsv", sep="\t", index_col=0, usecols=["srx", "FBgn"], squeeze=True
    ).to_dict()

    for srx in srx2fbgn.keys():
        try:
            make_bed(srx, fbgn=srx2fbgn[srx])
        except FileNotFoundError:
            print(srx)


def make_bed(srx, fbgn):
    narrow_peak = f"../output/chipseq-wf/narrowPeak/{srx}_peaks.narrowPeak"
    bed_file = f"../output/chipseq-wf/bed/{srx}_peaks.bed"
    with open(narrow_peak, "r") as file_in, open(bed_file, "w") as file_out:
        reader = csv.reader(file_in, delimiter="\t")
        writer = csv.writer(file_out, delimiter="\t")
        for row in reader:
            peak = "_".join(row[3].split("_")[-2:])
            writer.writerow([row[0], row[1], row[2], f"{fbgn}_{srx}_{peak}", row[8]])


if __name__ == "__main__":
    main()
