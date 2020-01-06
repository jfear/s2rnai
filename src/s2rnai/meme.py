#!/usr/bin/env python
""" Quick and Dirty Meme Parser """
import re
from textwrap import dedent
import numpy as np

class memeFile(object):
    def __init__(self, filename):
        self.filename = filename
        self.tfs = {}
        self._parse()

    def _parse(self):
        if self.filename.endswith('.meme'):
            with open(self.filename) as fh:
                f = fh.read()
        else:
            raise ValueError('File must end with .meme')

        self.meme_version = int(re.findall(r'MEME version (\d+)', f)[0])
        self.meme_alphabet = re.findall(r'ALPHABET= (.*?)\n', f)[0]
        self.meme_strands = re.findall(r'strands: (.*?)\n', f)[0]
        bg = re.findall(r'A (.*?) C (.*) G (.*?) T (.*?)\n', f)[0]
        self.meme_bg = {'A': float(bg[0]), 'C': float(bg[1]), 'G': float(bg[2]), 'T': float(bg[3])}

        for m in re.findall('MOTIF .*?URL.*?\n', f, flags=re.DOTALL):
            tf = memeTF.from_string(m)
            if tf.id in self.tfs:
                self.tfs[tf.id].append(tf)
            else:
                self.tfs[tf.id] = [tf,]

    def __getitem__(self, key):
        return self.tfs[key]

    def __iter__(self):
        for key in self.tfs.keys():
            yield key

    def items(self):
        for item in self.tfs.items():
            yield item

    def keys(self):
        return self.tfs.keys()

    def count(self):
        counter = 0
        for key in self.tfs.keys():
            for entry in self.tfs[key]:
                counter += 1

        return counter

    def replace_id(self, orig, new):
        tfs = []
        for tf in self.tfs[orig]:
            tf.oldID = tf.id
            tf.id = new
            tfs.append(tf)

        self.tfs[new] = tfs
        del self.tfs[orig]



MEMETFREGEX="""MOTIF (?P<id>[A-Za-z]+\d+\.*\d*_*\d*) (?P<name>.*?)

letter-probability matrix: alength= (?P<alength>\d+) w= (?P<w>\d+) nsites= (?P<nsites>\d+) E= (?P<e>\d+)
(?P<matrix>.*?)

URL (?P<url>.*)"""

class memeTF(object):
    def __init__(self, id, name, alength, w, nsites, e, matrix, url):
        """ This is a simple parser for a meme transcription factor block. """
        self.id = id
        self.fbgn = id if id.startswith('FBgn') else None
        self.name = name
        self.alength = int(alength)
        self.w = int(w)
        self.nsites = int(nsites)
        self.e = float(e)
        self.matrix = matrix
        self.url = url

    @classmethod
    def from_string(cls, block):
        regex = re.compile(MEMETFREGEX, flags=re.DOTALL)
        match = re.match(regex, block.strip())
        try:
            return cls(**match.groupdict())
        except:
            print(block)
            raise

    @property
    def id(self):
        return self._id

    @id.setter
    def id(self, value):

        if len(value.split('_')) > 1:
            s = value.split('_')
            self._id = s[0]
            self.count = int(s[1])
        elif len(value.split('.')) > 1:
            s = value.split('.')
            self._id = s[0]
            self.count = int(s[1])
        else:
            self._id = value
            self.count = 0

    @property
    def matrix(self):
        return self._matrix

    @matrix.setter
    def matrix(self, value):
        clean = []
        for row in value.split('\n'):
            clean.append([float(x) for x in row.strip().split()])
        self._matrix = np.matrix(clean)


if __name__ == '__main__':
    pass
