from pathlib import Path
import shutil

def main():
    for file_name in Path("/Volumes/Promise_Pegasus/bergeric/s2cell-prior/chipseq-wf/data/chipseq_peaks/macs2/").glob("**/*.narrowPeak"):
        try:
            shutil.copy2(file_name, Path("narrowPeak", Path("../output/chipseq-wf/", file_name.name))
        except Exception:
            print(file_name.name)


if __name__ == "__main__":
    main()
