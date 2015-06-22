"""
NSR2Tm.py

Usage:
    NSR2Tm.py (-n nsrfile) (-o outputfile)
    NSR2Tm.py -h | --help
    NSR2Tm.py -v | --version

Options:
    -n nsrfile      File of NSR sequence
    -o outputfile   File of histogram
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
    outputfile = args['-o']

    df_nsr = pd.read_csv(nsrfile)
    seq_nsr = df_nsr['seq'].values
    Nmer = len(seq_nsr[0])

    tm_nsr = [calculate_Tm(seq) for seq in seq_nsr]
    df_nsr['Tm'] = tm_nsr
    text_str = ('Mean=' + str(round(np.mean(tm_nsr), 2)) + '\n S.D.=' +
                str(round(np.std(tm_nsr), 2)))

    # plotting
    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    ax.hist(tm_nsr, bins=28, range=(8, 36))
    ax.set_xlabel('Tm of NSR [Degree Celsius]')
    ax.set_ylabel('Number of NSR')
    ax.set_title('Histogram for Tm of NSR') 
    ax.text(0.95, 0.95, text_str, verticalalignment='top',
            horizontalalignment='right',
            transform=ax.transAxes,
            color='black', fontsize=12)
    plt.tight_layout()
    plt.savefig(outputfile)
