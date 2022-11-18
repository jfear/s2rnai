import os
import tarfile
from urllib.request import urlretrieve
from tempfile import TemporaryDirectory
from pathlib import Path
import shutil as sh


def main():
    with TemporaryDirectory() as tmpdir:
        # Download tar
        fname = Path(tmpdir, "db.tgz")
        urlretrieve(url=snakemake.params.url, filename=fname)

        # extract tar
        with tarfile.open(fname) as tar:
            def is_within_directory(directory, target):
                
                abs_directory = os.path.abspath(directory)
                abs_target = os.path.abspath(target)
            
                prefix = os.path.commonprefix([abs_directory, abs_target])
                
                return prefix == abs_directory
            
            def safe_extract(tar, path=".", members=None, *, numeric_owner=False):
            
                for member in tar.getmembers():
                    member_path = os.path.join(path, member.name)
                    if not is_within_directory(path, member_path):
                        raise Exception("Attempted Path Traversal in Tar File")
            
                tar.extractall(path, members, numeric_owner=numeric_owner) 
                
            
            safe_extract(tar, path=tmpdir, members=target_members(tar))

        for file_name in snakemake.output:
            name = Path(file_name).name
            sh.move(Path(tmpdir, "motif_databases/FLY", name), file_name)


def target_members(members):
    for info in members:
        if "FLY" in info.name:
            yield info


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from s2rnai.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="motif-wf",
            params=dict(
                url="http://meme-suite.org/meme-software/Databases/motifs/motif_databases.12.19.tgz"
            ),
        )

    main()
