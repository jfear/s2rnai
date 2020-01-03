import unittest
import sys
from s2rnai import meme

MOTIF="""MOTIF FBgn0000014 AbdA_Cell

letter-probability matrix: alength= 4 w= 7 nsites= 18 E= 0
  0.055556	  0.166667	  0.000000	  0.777778
  0.000000	  0.000000	  0.000000	  1.000000
  0.888889	  0.000000	  0.000000	  0.111111
  1.000000	  0.000000	  0.000000	  0.000000
  0.055556	  0.000000	  0.000000	  0.944444
  0.000000	  0.000000	  0.333333	  0.666667
  0.833333	  0.055556	  0.111111	  0.000000

URL http://pgfe.umassmed.edu/TFDBS/TFdetails.php?FlybaseID=FBgn0000014

"""

class testMeme(unittest.TestCase):
    def test_memeTF(self):
        tf = meme.memeTF.from_string(MOTIF)
        self.assertEqual('FBgn0000014', tf.id)
        self.assertEqual(0, tf.count)
        self.assertEqual('AbdA_Cell', tf.name)
        self.assertEqual(4, tf.alength)
        self.assertEqual(7, tf.w)
        self.assertEqual(18, tf.nsites)
        self.assertEqual(0.0, tf.e)
        self.assertEqual(0.055556, tf.matrix[0, 0])
        self.assertEqual(0.833333, tf.matrix[-1, 0])
        self.assertEqual(0.333333, tf.matrix[-2, 2])

    def test_memeParse(self):
        tfs = meme.memeFile('data/external/meme/motif_databases/FLY/fly_factor_survey.meme')
        self.assertEqual(4, tfs.meme_version)
        self.assertEqual('ACGT', tfs.meme_alphabet)
        self.assertEqual('+ -', tfs.meme_strands)
        self.assertDictEqual({'A': 0.25, 'T': 0.25, 'C': 0.25, 'G': 0.25}, tfs.meme_bg)
        self.assertEqual('AbdA_Cell', tfs['FBgn0000014'][0].name)

