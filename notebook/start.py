# Imports
import os
import sys
import yaml

from joblib import Memory
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

# Project level imports
sys.path.insert(0, '../lib')
from s2rnai.notebook import Nb

# Setup notebook
nbconfig = Nb.setup_notebook()

# Turn on cache
memory = Memory(cachedir=nbconfig.cache, verbose=0)

# Set up references
with open('../../config/config.yml') as fh:
    config = yaml.load(fh)

assembly = config['assembly']
tag = config['aligner']['tag']
REF = os.path.join(os.environ['REFERENCES_DIR'], assembly, tag)
