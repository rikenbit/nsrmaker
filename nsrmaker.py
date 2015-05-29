"""
nsrmaker.py

Usage:
    nsrmaker.py (-s <species>) (-n <length>) (-r <rRNA>...) (-o <output>)
    nsrmaker.py -h | --help
    nsrmaker.py -v | --version

Options:
    -s <species>    Species
    -n <length>     Length of NSR
    -r <rRNA>...    List of rRNA to be removed
    -o <output>     File name of output
    -h --help       Show this screen
    -v --version    Show version
"""

    # nsrmaker.py (-s species) [-l length] [-r rRNAs] (-o output)
from __future__ import print_function
import pandas as pd
from Bio import SeqIO
from docopt import docopt
# from Bio.Seq import Seq
# from Bio.Alphabet import IUPAC

# Todo ###
# write docopt: species, help, outputfile
# comparison among species
# add many species, ex. fly, mamoset, worm, genome-read organisms
# think about meta-transcriptome using NSR
# check unbiassed distribution of designed NSR, make function


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


def find_seq_in_rRNA(rRNA_seq_revcom, possible_oligo):
    df = pd.DataFrame()
    for key in rRNA_seq_revcom.keys():
        s = pd.Series([rRNA_seq_revcom[key].find(possible_oligo[i][0:6]) for i
                      in range(len(possible_oligo))])
        df[key] = s

    index_list = ([i for i in range(len(N6_orig)) if all(df.ix[i] == -1)])
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
    species = 'Hsapiens'

    # make N6 random primer
    N6_orig = make_random_oligo(6)

    # read rRNA seq
    rRNA_seq_revcom = read_rRNA(species, return_revcom=True,
    #                             rRNA_subunit=('5S', '23S', '16S'))
    #                            rRNA_subunit=('45S', '5S', '12S', '16S'))
    #                             rRNA_subunit=('45S', '5Ss', '5So', '12S', '16S'))
    #                             rRNA_subunit=('28S', '18S', '5p8S', '5S','12S', '16S'))
                                 rRNA_subunit=('28S', '18S', '12S', '16S'))

    # calculate NSR
    index_NSR, res_df = find_seq_in_rRNA(rRNA_seq_revcom, N6_orig)
    seq_NSR = [N6_orig[i] for i in index_NSR]
    name_NSR = ["NSR_" + species + "_" + str(i).zfill(4) for i in index_NSR]
    df_NSR = pd.DataFrame({'name': name_NSR, 'seq': seq_NSR})
    output_file = "./results/NSR_" + species + ".csv"
    df_NSR.to_csv(output_file)

    # comparison to published NSR seq
    if species == 'Mmusculus' or species == 'Hsapiens':
        seq_NSR_public = get_public_NSR(species)
        df_c = pd.DataFrame({'seq': N6_orig})
        df_c['new'] = df_c['seq'].isin(seq_NSR)
        df_c['public'] = df_c['seq'].isin(seq_NSR_public)
        TT = [df_c['new'].ix[i] == True and df_c['public'].ix[i] == True
              for i in range(len(N6_orig))]
        TF = [df_c['new'].ix[i] == True and df_c['public'].ix[i] == False
              for i in range(len(N6_orig))]
        FT = [df_c['new'].ix[i] == False and df_c['public'].ix[i] == True
              for i in range(len(N6_orig))]
        FF = [df_c['new'].ix[i] == False and df_c['public'].ix[i] == False
              for i in range(len(N6_orig))]
        print("Comparison between new NSR and public NSR: " + species)
        print("[both, only in New, only in Public, neither]")
        print([pd.Series(x).sum() for x in (TT, TF, FT, FF)])
