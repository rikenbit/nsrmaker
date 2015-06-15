"""
map2transcriptome.py

Usage:
    map2transcriptome.py (-n nsrfile) (-t transcriptomefile) (-o outputfile)
    map2transcriptome.py -h | --help
    map2transcriptome.py -v | --version

Options:
    -n nsrfile              File of NSR sequence
    -t transcriptomefile    File of transcriptome sequence
    -o outputfile           File of output
    -h --help               Show this screen
    -v --version            Show version
"""

from __future__ import print_function
import pandas as pd
from Bio import SeqIO
from docopt import docopt
# from Bio.Seq import Seq


if __name__ == '__main__':

    NAME = "map2transcriptome.py"
    VERSION = "1.0.0"

    args = docopt(__doc__, version="{0} {1}".format(NAME, VERSION))
    nsrfile = args['-n']
    transcriptomefile = args['-t']
    outputfile = args['-o']

    # reading NSR file
    df_nsr = pd.read_csv(nsrfile)
    seq_nsr = df_nsr['seq'].values

    # reading transcriptome file
    seq_transcriptome = []
    id_transcriptome = []
    for seq_record in SeqIO.parse(transcriptomefile, "fasta"):
        seq_transcriptome.append(str(seq_record.seq))
        id_transcriptome.append(seq_record.id)

    len_transcriptome = [len(seq_transcriptome[i]) for i
                         in range(len(seq_transcriptome))]

    match_transcriptome = []
    for transcript in seq_transcriptome:
        match_count = 0
        for nsr in seq_nsr:
            c = transcript.count(nsr)
            match_count += c
        match_transcriptome.append(match_count)

    df_transcriptome = pd.DataFrame({'id': id_transcriptome,
                                     'seq': seq_transcriptome,
                                     'len': len_transcriptome,
                                     'match': match_transcriptome})

    df_transcriptome.to_csv(outputfile)
