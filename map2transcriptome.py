"""
map2transcriptome.py

Usage:
    map2transcriptome.py (-n nsrfile) (-t transcriptomefile) (-o outputname)
    map2transcriptome.py -h | --help
    map2transcriptome.py -v | --version

Options:
    -n nsrfile              File of NSR sequence
    -t transcriptomefile    File of transcriptome sequence
    -o outputname           Name of output
    -h --help               Show this screen
    -v --version            Show version
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


if __name__ == '__main__':

    NAME = "map2transcriptome.py"
    VERSION = "1.0.0"

    args = docopt(__doc__, version="{0} {1}".format(NAME, VERSION))
    nsrfile = args['-n']
    transcriptomefile = args['-t']
    outputname = args['-o']

    # reading NSR file
    df_nsr = pd.read_csv(nsrfile)
    seq_nsr = df_nsr['seq'].values

    # downsampling of nsr
    a = len(seq_nsr) / 100.
    num_select = [100 * i for i in range(1, np.int(np.ceil(a)))]
    num_select.reverse()
    list_nsr_set = []
    list_nsr_set.append(seq_nsr)
    for i in num_select:
        this_nsr = np.random.choice(seq_nsr, i)
        list_nsr_set.append(this_nsr)

    # reading transcriptome file
    seq_transcriptome = []
    id_transcriptome = []
    for seq_record in SeqIO.parse(transcriptomefile, "fasta"):
        seq_transcriptome.append(str(
            seq_record.seq.reverse_complement()).upper())
        id_transcriptome.append(seq_record.id)

    len_transcriptome = [len(seq_transcriptome[i]) for i
                         in range(len(seq_transcriptome))]

    df_transcriptome = pd.DataFrame({'id': id_transcriptome,
                                     'seq': seq_transcriptome,
                                     'len': len_transcriptome})

    for nsr_set in list_nsr_set:
        match_transcriptome = []
        for transcript in seq_transcriptome:
            match_count = 0
            for nsr in nsr_set:
                c = transcript.count(nsr)
                match_count += c
            match_transcriptome.append(match_count)
        match_transcriptome_100 = [i * 100 for i in match_transcriptome]
        df_transcriptome['m'+str(len(nsr_set))] = (match_transcriptome_100
                                                   / df_transcriptome['len'])

    df_transcriptome.to_csv(outputname + '.csv')

    colname = ['m' + str(len(nsr_set)) for nsr_set in list_nsr_set]
    df_plot = df_transcriptome[colname]

    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    for col in colname:
        ax.hist(df_plot[col], bins=25, range=(0, 25))
    ax.set_xlabel('Number of NSR hexamers (per 100 nt)')
    ax.set_ylabel('Number of transcripts')
    ax.legend(colname)
    plt.savefig(outputname + '_hist_summary.png')

    ncols = 3
    nrows = np.int(np.ceil(len(list_nsr_set) * 1.0 / ncols))
    if nrows == 1:
        nrows = 2
    ymax_hist = len(id_transcriptome) / 2

    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(4 * ncols, 4 * nrows))
    for k in range(len(list_nsr_set)):
        i, j = divmod(k, ncols)
        axes[i, j].hist(df_plot[colname[k]], bins=25, range=(0, 25))
        axes[i, j].set_ylim([0, ymax_hist])
        axes[i, j].set_xlabel('Number of NSR hexamers (per 100 nt)')
        axes[i, j].set_ylabel('Number of transcripts')
        axes[i, j].set_title(colname[k])
    plt.tight_layout()
    plt.savefig(outputname + '_hist.png')
