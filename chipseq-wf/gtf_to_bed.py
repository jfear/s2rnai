import csv
import re

FBGN_PATTERN = re.compile(r"FBgn\d+")


def main():
    with open("../output/chipseq-wf/dmel-all-r6.26_genes.gtf", "r") as file_in, open(
        "../output/chipseq-wf/dmel-all-r6.26_genes.bed", "w"
    ) as file_out:
        reader = csv.reader(file_in, delimiter="\t")
        writer = csv.writer(file_out, delimiter="\t", quoting=csv.QUOTE_NONE)
        for row in reader:
            if row[0] not in ["X", "Y", "2L", "2R", "3L", "3R", "4"]:
                continue

            writer.writerow(
                [
                    f"chr{row[0]}",
                    row[3],
                    row[4],
                    re.findall(FBGN_PATTERN, row[-1])[0],
                    row[5],
                    row[6],
                ]
            )


if __name__ == "__main__":
    main()
