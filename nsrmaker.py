"""
nsrmaker.py

Usage:
    nsrmaker.py (-s species) (-n length) (-e exclude) (-r rRNA...) (-t trimsd) (-o outprefix)
    nsrmaker.py fasta (-f fasta...) (-s species) (-n length) (-e exclude) (-t trimsd) (-o outprefix)
    nsrmaker.py -h | --help
    nsrmaker.py -v | --version

Options:
    -f fasta...     List of fasta files to be removed
    -s species      Species
    -n length       Length of NSR
    -e exclude      x base matching to be removed
    -r rRNA...      List of rRNA to be removed
    -t trimsd       SD value to trim NSR by its Tm
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


def get_digits_num(input_int):
    n = 0
    while True:
        d, _ = divmod(input_int, 10 ** n)
        if d == 0:
            break
        else:
            n += 1
    return n


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


def make_random_oligo(N):
    oligo_quad = [convert_repr_int(i, 4, N) for i in range(4 ** N)]
    dict_gatc = {'0': 'G', '1': 'A', '2': 'T', '3': 'C'}
    oligo_gatc = [replaces_str(s, dict_gatc) for s in oligo_quad]
    return oligo_gatc


def read_fasta(fasta_list, return_revcom=False):
    dict_rRNA = dict()
    for fasta in fasta_list:
        for seq_record in SeqIO.parse(fasta, "fasta"):
            dict_rRNA[seq_record.id] = seq_record.seq
    
    dict_rRNA_revcom = dict()
    for key in dict_rRNA.keys():
        dict_rRNA_revcom[key] = dict_rRNA[key].reverse_complement()

    if return_revcom is False:
        return dict_rRNA
    if return_revcom is True:
        return dict_rRNA_revcom
    return None


def read_rRNA(species, dir_rRNA='./rRNA', return_revcom=False,
              rRNA_subunit=('28S', '18S', '16S', '12S')):

    dict_rRNA = dict()
    for subunit in rRNA_subunit:
        for seq_record in SeqIO.parse(dir_rRNA + "/" + species +
                                      "_rRNA_" + subunit + ".fasta", "fasta"):
            dict_rRNA[subunit] = seq_record.seq

    dict_rRNA_revcom = dict()
    for key in dict_rRNA.keys():
        dict_rRNA_revcom[key] = dict_rRNA[key].reverse_complement()

    if return_revcom is False:
        return dict_rRNA
    if return_revcom is True:
        return dict_rRNA_revcom
    return None


def find_seq_in_rRNA(rRNA_seq_revcom, possible_oligo, xbasematching):
    df = pd.DataFrame()
    for key in rRNA_seq_revcom.keys():
        s = pd.Series([rRNA_seq_revcom[key].find(possible_oligo[i][
            (N-xbasematching):N]) for i in range(len(possible_oligo))])
        df[key] = s

    index_list = ([i for i in range(len(random_orig)) if all(df.ix[i] == -1)])
    return index_list, df


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


def calculate_Tm(seq):
    num_GC = seq.count("G") + seq.count("C")
    num_AT = seq.count("A") + seq.count("T")
    output = 4 * int(num_GC) + 2 * int(num_AT)
    return output


def plot_hist_tm(df_NSR, outputfile_tm):
    tm_NSR = df_NSR['tm'].values
    mean_tm = np.mean(tm_NSR)
    sd_tm = np.std(tm_NSR)
    text_str = ('Mean=' + str(round(mean_tm, 2)) + '\n S.D.=' +
                str(round(sd_tm, 2)))
    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    ax.hist(tm_NSR, bins=28, range=(8, 36))
    ax.set_xlabel('Tm of NSR [Degree Celsius]')
    ax.set_ylabel('Number of NSR')
    ax.set_title('Histogram for Tm of NSR')
    ax.text(0.95, 0.95, text_str, verticalalignment='top',
            horizontalalignment='right',
            transform=ax.transAxes,
            color='black', fontsize=12)
    plt.tight_layout()
    plt.savefig(outputfile_tm)


def trim_nsr_by_tm(df_NSR, trim_sd):
    tm_NSR = df_NSR['tm'].values
    mean_tm = np.mean(tm_NSR)
    sd_tm = np.std(tm_NSR)
    upper_tm = mean_tm + sd_tm * trim_sd
    lower_tm = mean_tm - sd_tm * trim_sd
    idx = [lower_tm < tm_NSR[i] < upper_tm for i in range(len(tm_NSR))]
    return df_NSR[idx]


def df_to_revseq_df(df_NSR):
    seq_NSR = df_NSR['seq'].values
    tm_NSR = df_NSR['tm'].values
    name_NSR = df_NSR['name'].values
    seq_NSR_rev = [str(Seq(s).reverse_complement()) for s in seq_NSR]
    name_NSR_rev = [s + "_rev" for s in name_NSR]
    tm_NSR_rev = tm_NSR
    df_NSR_rev = pd.DataFrame({'name': name_NSR_rev, 'seq': seq_NSR_rev,
                               'tm': tm_NSR_rev})
    return df_NSR_rev


class OptionError(ValueError):
    pass


if __name__ == '__main__':

    NAME = "nsrmaker.py"
    VERSION = "2.0.0"

    args = docopt(__doc__, version="{0} {1}".format(NAME, VERSION))

    species = args['-s']
    N = int(args['-n'])
    xbasematching = int(args['-e'])
    trim_sd = float(args['-t'])
    outprefix = args['-o']

    isfastamode = args['fasta']
    if isfastamode:
        list_fasta = args['-f']
    else:
        rRNA_subunit = args['-r']

#    str_rRNA_subunit = str()
#    for i in range(len(rRNA_subunit)):
#        str_rRNA_subunit = str_rRNA_subunit + str(rRNA_subunit[i])
#    outputname = ("./results/NSR_" + species + "_" + str(N) +
#                  "mer_" + str(xbasematching) + "basematchremove_"
#                  + str_rRNA_subunit)
    outputname = (outprefix + "_" + str(N) +
                  "mer_" + str(xbasematching) + "basematchremove")

    # make N6 random primer
    random_orig = make_random_oligo(N)
    digits_num = get_digits_num(4 ** N)

    # read rRNA seq
    if isfastamode:
        rRNA_seq_revcom = read_fasta(list_fasta, return_revcom=True)
    else:
        rRNA_seq_revcom = read_rRNA(species, return_revcom=True,
                                    rRNA_subunit=rRNA_subunit)

    # calculate NSR
    index_NSR, res_df = find_seq_in_rRNA(rRNA_seq_revcom, random_orig,
                                         xbasematching)
    seq_NSR = [random_orig[i] for i in index_NSR]
    name_NSR = ["NSR_" + species + "_" + str(i).zfill(digits_num)
                for i in index_NSR]

    # calculate Tm
    tm_NSR = [calculate_Tm(seq) for seq in seq_NSR]

    # export csv file
    df_NSR = pd.DataFrame({'name': name_NSR, 'seq': seq_NSR, 'tm': tm_NSR})
    outputfile_nsr = outputname + ".csv"
    df_NSR.to_csv(outputfile_nsr, index=False)

    df_NSR_rev = df_to_revseq_df(df_NSR)
    outputfile_nsr_rev = outputname + "_2nd.csv"
    df_NSR_rev.to_csv(outputfile_nsr_rev, index=False)

    # plotting histogram of NSR
    outputfile_tm = outputname + "_tm.png"
    plot_hist_tm(df_NSR, outputfile_tm)

    # trim nsr by tm_sd
    df_NSR_trimed = trim_nsr_by_tm(df_NSR, trim_sd)
    # export csv of NSR trimed
    outputfile_nsr_trimed = outputname + "_trimSD_" + str(trim_sd) + ".csv"
    df_NSR_trimed.to_csv(outputfile_nsr_trimed, index=False)

    df_NSR_trimed_rev = df_to_revseq_df(df_NSR_trimed)
    outputfile_nsr_trimed_rev = (outputname + "_trimSD_" +
                                 str(trim_sd) + "_2nd.csv")
    df_NSR_trimed_rev.to_csv(outputfile_nsr_trimed_rev, index=False)

    # plot histogram of NSR trimed
    outputfile_tm_trimed = outputname + "_trimSD_" + str(trim_sd) + "_tm.png"
    plot_hist_tm(df_NSR_trimed, outputfile_tm_trimed)

#    # comparison to published NSR seq
#    if species == 'Mmusculus' or species == 'Hsapiens':
#        seq_NSR_public = get_public_NSR(species)
#        df_c = pd.DataFrame({'seq': random_orig})
#        df_c['new'] = df_c['seq'].isin(seq_NSR)
#        df_c['public'] = df_c['seq'].isin(seq_NSR_public)
#        TT = [df_c['new'].ix[i] == True and df_c['public'].ix[i] == True
#              for i in range(len(random_orig))]
#        TF = [df_c['new'].ix[i] == True and df_c['public'].ix[i] == False
#              for i in range(len(random_orig))]
#        FT = [df_c['new'].ix[i] == False and df_c['public'].ix[i] == True
#              for i in range(len(random_orig))]
#        FF = [df_c['new'].ix[i] == False and df_c['public'].ix[i] == False
#              for i in range(len(random_orig))]
#        print("Comparison between new NSR and public NSR: " + species)
#        print("[both, only in New, only in Public, neither]")
#        print([pd.Series(x).sum() for x in (TT, TF, FT, FF)])
