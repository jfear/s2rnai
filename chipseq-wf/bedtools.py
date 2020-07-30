from pathlib import Path
from pybedtools import BedTool


def main():
    out_dir = Path("../output/chipseq-wf/intersections")
    out_dir.mkdir(exist_ok=True)

    ref = "../output/chipseq-wf/dmel-all-r6.26_genes.bed"
    ref_bt = BedTool(ref).slop(b=1_000, genome="dm6")

    for file_name in Path("../output/chipseq-wf/bed").iterdir():
        file_out = out_dir / file_name.name
        bt = BedTool(file_name)
        bt.intersect(ref_bt, wb=True).saveas(file_out)


if __name__ == "__main__":
    main()
