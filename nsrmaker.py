"""
nsrmaker.py

Usage:
    nsrmaker.py (-s species) (-n length) (-e exclude) (-r rRNA...)
    nsrmaker.py -h | --help
    nsrmaker.py -v | --version

Options:
    -s species      Species
    -n length       Length of NSR
    -e exclude      x base matching to be removed
    -r rRNA...      List of rRNA to be removed
    -h --help       Show this screen
    -v --version    Show version
"""

from __future__ import print_function
import pandas as pd
from Bio import SeqIO
from docopt import docopt
# from Bio.Seq import Seq

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
        s = pd.Series([rRNA_seq_revcom[key].find(possible_oligo[i][(N-xbasematching):N]) for i
                      in range(len(possible_oligo))])
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


class OptionError(ValueError):
    pass


if __name__ == '__main__':

    NAME = "nsrmaker.py"
    VERSION = "1.0.0"

    args = docopt(__doc__, version="{0} {1}".format(NAME, VERSION))

    species = args['-s']
    N = int(args['-n'])
    rRNA_subunit = args['-r']
    xbasematching = int(args['-e'])

    str_rRNA_subunit = str()
    for i in range(len(rRNA_subunit)):
        str_rRNA_subunit = str_rRNA_subunit + str(rRNA_subunit[i])
    output_file = ("./results/NSR_" + species + "_" + str(N) +
                   "mer_" + str(xbasematching) + "basematchremove_" + str_rRNA_subunit + ".csv")

    # make N6 random primer
    random_orig = make_random_oligo(N)
    digits_num = get_digits_num(4 ** N)

    # read rRNA seq
    rRNA_seq_revcom = read_rRNA(species, return_revcom=True,
                                rRNA_subunit=rRNA_subunit)

    # calculate NSR
    index_NSR, res_df = find_seq_in_rRNA(rRNA_seq_revcom, random_orig, xbasematching)
    seq_NSR = [random_orig[i] for i in index_NSR]
    name_NSR = ["NSR_" + species + "_" + str(i).zfill(digits_num)
                for i in index_NSR]
    df_NSR = pd.DataFrame({'name': name_NSR, 'seq': seq_NSR})
    df_NSR.to_csv(output_file, index=False)

    # comparison to published NSR seq
    if species == 'Mmusculus' or species == 'Hsapiens':
        seq_NSR_public = get_public_NSR(species)
        df_c = pd.DataFrame({'seq': random_orig})
        df_c['new'] = df_c['seq'].isin(seq_NSR)
        df_c['public'] = df_c['seq'].isin(seq_NSR_public)
        TT = [df_c['new'].ix[i] == True and df_c['public'].ix[i] == True
              for i in range(len(random_orig))]
        TF = [df_c['new'].ix[i] == True and df_c['public'].ix[i] == False
              for i in range(len(random_orig))]
        FT = [df_c['new'].ix[i] == False and df_c['public'].ix[i] == True
              for i in range(len(random_orig))]
        FF = [df_c['new'].ix[i] == False and df_c['public'].ix[i] == False
              for i in range(len(random_orig))]
        print("Comparison between new NSR and public NSR: " + species)
        print("[both, only in New, only in Public, neither]")
        print([pd.Series(x).sum() for x in (TT, TF, FT, FF)])
