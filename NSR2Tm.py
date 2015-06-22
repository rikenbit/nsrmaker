"""
NSR2Tm.py

Usage:
    NSR2Tm.py (-n nsrfile)
    NSR2Tm.py -h | --help
    NSR2Tm.py -v | --version

Options:
    -n nsrfile      File of NSR sequence
    -h --help       Show this screen
    -v --version    Show version
"""

from __future__ import print_function
import pandas as pd
from Bio import SeqIO
from docopt import docopt
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('ggplot')


def calculate_Tm(seq):
    num_GC = seq.count("G") + seq.count("C")
    num_AT = seq.count("A") + seq.count("T")
    output = 4 * int(num_GC) + 2 * int(num_AT)
    return output


if __name__ == '__main__':

    NAME = "calculateTm.py"
    VERSION = "1.0.0"

    args = docopt(__doc__, version="{0} {1}".format(NAME, VERSION))
    nsrfile = args['-n']

    df_nsr = pd.read_csv(nsrfile)
    seq_nsr = df_nsr['seq'].values

    tm_nsr = [calculate_Tm(seq) for seq in seq_nsr]
    df_nsr['Tm'] = tm_nsr
    
