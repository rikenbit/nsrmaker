"""
check_base_bias.py

Usage:
    check_base_bias.py -f file
    check_base_bias.py -h | --help
    check_base_bias.py -v | --version

Options:
    -f file         File of NSR
    -h --help       Show this screen
    -v --version    Show version
"""

from __future__ import print_function
import pandas as pd
from docopt import docopt
from collections import Counter
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# plt.style.use('ggplot')


def count_base(base, dict_base):
    dict_base[base] += 1


if __name__ == '__main__':

    NAME = "check_base_bias.py"
    VERSION = "0.1.0"

    args = docopt(__doc__, version="{0} {1}".format(NAME, VERSION))

    input_file = args['-f']

    df = pd.read_csv(input_file)
    seq = list(df['seq'])

    base_list = [[seq[i][j] for i in range(len(seq))]
                 for j in range(len(seq[0]))]

    base_count = [Counter(base_list[i]) for i in range(len(base_list))]

    base_df = pd.DataFrame(base_count)

    output_file = input_file.replace(".csv", "_basebias.png")

    base_df = (base_df.T * 100.0 / base_df.T.sum()).T
    base_df['base'] = range(1, 7, 1)
    base_df_melt = pd.melt(base_df, id_vars='base',
                           value_vars=['A', 'C', 'G', 'T'],
                           var_name='nucleotide', value_name='rate')

    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    sns.pointplot(x='base', y='rate', hue='nucleotide',
                  data=base_df_melt, ax=ax)
    ax.set_ylim(10, 50)
    ax.set_ylabel('rate')
    fig.savefig(output_file)
