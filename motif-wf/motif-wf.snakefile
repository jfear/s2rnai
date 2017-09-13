# ** snakemake ** #
"""The motif workflow."""
import os
import sys
from tempfile import NamedTemporaryFile, mkdtemp
from urllib.request import urlretrieve, urlopen
import tarfile

import pandas as pd

from lcdblib.snakemake import helpers
from lcdblib.utils import utils

sys.path.insert(0, '../lcdb-wf')
from lib import common

sys.path.insert(0, '../lib/python')
from s2rnai import meme

# ----------------------------------------------------------------------------
# SETUP
# ----------------------------------------------------------------------------

include: '../lcdb-wf/references.snakefile'

shell.prefix('set -euo pipefail; export TMPDIR={};'.format(common.tempdir_for_biowulf()))
shell.executable('/bin/bash')
references_dir = common.get_references_dir(config)

refdict, conversion_kwargs = common.references_dict(config)
assembly = config['assembly']

################################################################################
# Set up file naming patterns
################################################################################
patterns = {
    'meme': {
        'onTheFly': '../data/external/meme/OnTheFly_2014_Drosophila.meme',
        'flyFactor': '../data/external/meme/fly_factor_survey.meme',
        'flyReg': '../data/external/meme/flyreg.v2.meme',
        'dmmpmm2009': '../data/external/meme/dmmpmm2009.meme',
        'idmmpmm2009': '../data/external/meme/idmmpmm2009.meme',
        },
    'onTheFlyMap': 'data/onTheFlyMap.tsv',
    'fimo': {
        'xml': 'data/fimo/motif_alignments_{meme}_{ref}.xml',
        'html': 'data/fimo/motif_alignments_{meme}_{ref}.html',
        'txt': 'data/fimo/motif_alignments_{meme}_{ref}.txt',
        'gff': 'data/fimo/motif_alignments_{meme}_{ref}.gff',
        'log': 'data/fimo/motif_alignments_{meme}_{ref}.log',
        },
    'dm6': '{references_dir}/{assembly}/{tag}/fasta/{assembly}_{tag}.fasta'
}

fill = dict(meme=['flyFactor', 'onTheFly', 'flyReg', 'dmmpmm2009', 'idmmpmm2009'], ref='dm6')
targets = helpers.fill_patterns(patterns, fill)

rule targets:
    input: utils.flatten(targets)


rule download_meme:
    output: utils.flatten(patterns['meme'])
    params: url='http://meme-suite.org/meme-software/Databases/motifs/motif_databases.12.15.tgz'
    run:
        from Bio import motifs
        tmp = NamedTemporaryFile(dir=TMPDIR)

        # Download tar
        urlretrieve(url=params.url, filename=tmp.name)

        # Unpack tar to TMPDIR
        def fly(members):
            for info in members:
                if 'FLY' in info.name:
                    yield info

        with tarfile.open(tmp.name) as tar:
            tar.extractall(path=TMPDIR, members=fly(tar))

        # Iterate over outputs and move file to right spot
        for fn in output:
            name = os.path.basename(fn)
            shell('mv $TMPDIR/motif_databases/FLY/{0} {1}'.format(name, fn))

        # Remove tar file
        tmp.close()


rule map_otf_fbgn:
    input: targets['meme']['onTheFly']
    output: patterns['onTheFlyMap']
    run:
        onTheFlyTFS = meme.memeFile(input[0])

        results = []
        for key, values in onTheFlyTFS.items():
            for value in values:
                # grab link to protein page from url provided
                with urlopen(value.url) as fh:
                    page = fh.read().decode('UTF-8')
                    name = re.findall(r'protein_entry.php\?protein_ID=(.*?)\'', page)[0]
                    ID = re.findall(r'ID: (OTF\d+\.\d+)', page)[0]
                URL = 'https://bhapp.c2b2.columbia.edu/OnTheFly/cgi-bin/protein_entry.php?protein_ID={0}'
                try:
                    with urlopen(URL.format(name)) as fh:
                        fbgn = re.findall(r'FBgn\d+', fh.read().decode('utf-8'))[0]
                    results.append((key, ID, value.name, fbgn))
                except:
                    pass
        df = pd.DataFrame(results, columns=['motif_id', 'name', 'target_protein', 'FBgn'])
        df.to_csv(output[0], index=False, sep='\t')


def _fimo(wildcards):
    return [patterns['meme'][wildcards.meme], patterns[wildcards.ref]]


rule fimo:
    input: _fimo
    output:
        xml = patterns['fimo']['xml'],
        html = patterns['fimo']['html'],
        txt = patterns['fimo']['txt'],
        gff = patterns['fimo']['gff'],
        log = patterns['fimo']['log']
    run:
        fimoTMP = mkdtemp()

        shell(
            'fimo '
            '--oc {fimoTMP} '
            '{input[0]} '
            '{input[1]} '
            '2>{output.log}'
            )

        shell(
            'mv {fimoTMP}/fimo.xml {output.xml} '
            '&& mv {fimoTMP}/fimo.html {output.html} '
            '&& mv {fimoTMP}/fimo.txt {output.txt} '
            '&& mv {fimoTMP}/fimo.gff {output.gff} '
            '&& rm -rf {fimoTMP} '
            )
