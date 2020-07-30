import csv
from pathlib import Path


def main():
    with open("s2_chip.tsv", "w") as file_out:
        writer = csv.writer(file_out, delimiter="\t")
        writer.writerow(["peak_id", "chrom", "start", "end", "tf", "target", "score"])
        for file_name in Path("../output/chipseq-wf/intersections").iterdir():
            reader = csv.reader(file_name.open(), delimiter="\t")
            for row in reader:
                chrom = row[0]
                start = row[1]
                end = row[2]
                id_ = row[3]
                tf = id_.split("_")[0]
                srx = id_.split("_")[1]
                target = row[8]
                score = row[4]
                writer.writerow([id_, chrom, start, end, tf, target, score])


if __name__ == "__main__":
    main()
