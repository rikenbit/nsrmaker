from __future__ import print_function
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
import pandas as pd


class OptionError(ValueError):
    pass


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
    oligo_quad = [convert_repr_int(i, 4, 6) for i in range(4 ** N)]
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

    if return_revcom is False: return dict_rRNA
    if return_revcom is True: return dict_rRNA_revcom
    return None


def find_seq_in_rRNA(rRNA_seq_revcom, possible_oligo):
    df = pd.DataFrame()
    for key in rRNA_seq_revcom.keys():
        s = pd.Series([rRNA_seq_revcom[key].find(possible_oligo[i][0:6]) for i
                      in range(len(possible_oligo))])
        df[key] = s

    return len([i for i in range(len(N6_orig)) if all(df.ix[i] == -1)])


if __name__ == '__main__':

    # make N6 random primer
    N6_orig = make_random_oligo(6)

    # read rRNA seq
    rRNA_seq_revcom = read_rRNA('Hsapiens', return_revcom=True)

    # matching
    print(find_seq_in_rRNA(rRNA_seq_revcom, N6_orig))
