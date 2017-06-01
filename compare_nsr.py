"""
nsrmaker.py

Usage:
    compare_nsr.py (--n1 nsrfile1) (--n2 nsrfile2) (-n length) (-o outprefix)
    compare_nsr.py -h | --help
    compare_nsr.py -v | --version

Options:
    --n1 nsrfile    NSR file 1
    --n2 nsrfile    NSR file 2
    -n length       Length of NSR
    -o outprefix    Prefix of output files
    -h --help       Show this screen
    -v --version    Show version
"""

from __future__ import print_function
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from docopt import docopt
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('ggplot')

# Todo ###
# comparison among species
# add many species, ex. fly, mamoset, worm, genome-read organisms
# think about meta-transcriptome using NSR
# check unbiassed distribution of designed NSR, make function



def read_nsr(nsrfile):
    df_nsr = pd.read_csv(nsrfile)
    seq_nsr = df_nsr['seq'].values
    return seq_nsr
    

def get_public_NSR(species):
    if species == 'Mmusculus':
        df_mNSR = pd.read_csv('./published_NSR/mNSR_GenomRes_suppl3.csv',
                              header=0, names=('name_published', 'seq'))
        return [df_mNSR['seq'][i] for i in range(len(df_mNSR))]
    elif species == 'Hsapiens':
        df_hNSR = pd.read_csv('./published_NSR/hNSR_nmeth_suppl2.csv',
                              header=None, names=("seq_w_tail",))
        return [df_hNSR['seq_w_tail'][i][11:] for i in range(len(df_hNSR))]
    else:
        print('Option species must be either Mmusculus or Hsapiens')
        raise OptionError()



def make_random_oligo(N):
    oligo_quad = [convert_repr_int(i, 4, N) for i in range(4 ** N)]
    dict_gatc = {'0': 'G', '1': 'A', '2': 'T', '3': 'C'}
    oligo_gatc = [replaces_str(s, dict_gatc) for s in oligo_quad]
    return oligo_gatc

def convert_repr_int(input_int, target_repr, digit):

    '''
    Return the target_repr representation of an integer
    '''

    if input_int > (target_repr ** digit - 1):
        print('Cannot change representation')
        raise OptionError()

    output = str('')
    k = input_int

    for i in reversed(range(digit)):
        j, k = divmod(k, target_repr ** i)
        output += str(j)

    return output


def replaces_str(input_str, dict):
    output = input_str
    for key in dict.keys():
        output = output.replace(key, dict[key])
    return output




class OptionError(ValueError):
    pass


if __name__ == '__main__':

    NAME = "nsrmaker.py"
    VERSION = "2.0.0"

    args = docopt(__doc__, version="{0} {1}".format(NAME, VERSION))

    nsrfile1 = args['--n1']
    nsrfile2 = args['--n2']
    N = int(args['-n'])
    outputname = args['-o']

    # reading NSR file 1, 2
    seq_nsr1 = read_nsr(nsrfile1)
    seq_nsr2 = read_nsr(nsrfile2)

    # Preparing random N-mers
    random_orig = make_random_oligo(N)

    # comparing 1 to 2

    df_c = pd.DataFrame({'seq': random_orig})
    df_c['1'] = df_c['seq'].isin(seq_nsr1)
    df_c['2'] = df_c['seq'].isin(seq_nsr2)
    TT = [df_c['1'].ix[i] == True and df_c['2'].ix[i] == True
          for i in range(len(random_orig))]
    TF = [df_c['1'].ix[i] == True and df_c['2'].ix[i] == False
          for i in range(len(random_orig))]
    FT = [df_c['1'].ix[i] == False and df_c['2'].ix[i] == True
          for i in range(len(random_orig))]
    FF = [df_c['1'].ix[i] == False and df_c['2'].ix[i] == False
          for i in range(len(random_orig))]
    print("Comparison between 1 NSR and 2 NSR: ")
    print("[both, only in 1, only in 2, neither]")
    print([pd.Series(x).sum() for x in (TT, TF, FT, FF)])

    tt = pd.Series(TT).sum()
    tf = pd.Series(TF).sum()
    ft = pd.Series(FT).sum()
    ff = pd.Series(FF).sum()
    N_all = df_c.shape[0]

    df_result = pd.DataFrame({'nsrfile1': nsrfile1,
                              'nsrfile2': nsrfile2,
                              'both': tt,
                              'only_in_1': tf,
                              'only_in_2': ft,
                              'neither': ff,
                              'all': N_all 
                          }, index=[0])
    outputfile = outputname + ".csv"
    df_result.to_csv(outputfile, index=False)
