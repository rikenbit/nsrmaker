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
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('ggplot')


def count_base(base, dict_base):
    dict_base[base] += 1


if __name__ == '__main__':

    NAME = "check_base_bias.py"
    VERSION = "0.1.0"

    args = docopt(__doc__, version="{0} {1}".format(NAME, VERSION))

    input_file = args['-f']

    df = pd.read_csv(input_file)
    seq = list(df['seq'])

    base_list = [[seq[i][j] for i in range(len(seq))] for j in range(len(seq[0]))]

    base_count = [Counter(base_list[i]) for i in range(len(base_list))]

    base_df = pd.DataFrame(base_count)

    base_df.plot()

    output_file = input_file.replace(".csv", "_basebias.png")
    plt.savefig(output_file)
