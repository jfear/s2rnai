import pandas as pd

def main():
    df = pd.read_csv("../output/chipseq-wf/s2_chip.tsv", sep="\t", index_col=0)
    unique_peaks = df.groupby(["chrom", "start", "end", "tf", "target"]).score.max().reset_index()
    num_peaks_per_tf = unique_peaks.groupby("tf").size() / 1_000

    matrix = unique_peaks.groupby(["tf", "target"]).score.sum().div(num_peaks_per_tf).unstack(level=0)
    matrix.to_csv("../output/chipseq-wf/tf_weight_matrix.tsv", sep="\t")

if __name__ == "__main__":
    main()