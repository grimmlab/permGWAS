import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pathlib
import seaborn as sns
import scipy.stats as stats
from utils import helper_functions

def manhattan(df: pd.DataFrame, output: pathlib.Path, sig_level: int=5, min_p_val: np.array=None):
    """
    Create manhattan plot with Bonferroni threshold, if df_min is given, additionally with permutation based threshold.
    :param df: pandas DataFrame containing p-values
    :param output: filename to save plot
    :param min_p_val: min p-values for permutation based threshold
    """
    data = df[df['p_value'] <= 0.1].copy()
    data['-log10'] = -np.log10(data['p_value'])
    running_pos = 0
    cumul_pos = []
    for chrom, group_df in data.groupby('CHR'):
        cumul_pos.append(group_df['POS'] + running_pos)
        running_pos += group_df['POS'].max()
    data['Cumul_Pos'] = pd.concat(cumul_pos)

    g = sns.relplot(
        data=data,
        x='Cumul_Pos',
        y='-log10',
        aspect=4,
        hue='CHR',
        palette='Set1',
        linewidth=0,
        s=6,
        legend=None)

    g.ax.set_xlabel('Chromosome')
    g.ax.set_ylabel('-log10(p-value)')
    g.ax.set_xticks(data.groupby('CHR')['Cumul_Pos'].median())
    labels = np.unique(data['CHR'])
    g.ax.set_xticklabels(labels)

    g.axes[0, 0].axhline(-np.log10(helper_functions.compute_bonf_threshold(len(df), sig_level)), color='green',
                         label='Bonferroni')

    if min_p_val is not None:
        g.axes[0, 0].axhline(-np.log10(helper_functions.compute_perm_threshold(min_p_val, sig_level)), color='blue',
                             label='perm threshold')
    g.ax.legend(loc="upper right")
    plt.savefig(output)
    plt.clf()


def qq_plot(df: pd.DataFrame, output: pathlib.Path):
    """
    save QQ-plot
    :param df: dataframe with p-values
    :param output: filename to save plot
    """
    observed_p = -np.log10(np.sort(df['p_value']))
    expected_p = -np.log10(np.arange(1.0 / float(len(df)), 1, 1.0 / float(len(df) + 1)))
    inflation_factor = np.median(stats.chi2.isf(df['p_value'], 1)) / 0.456

    plt.figure(figsize=(4, 4))
    plt.plot(expected_p, observed_p, '.', markersize=4, markeredgewidth=0, alpha=0.8)
    plt.plot(expected_p, expected_p, 'k--', linewidth=0.75)
    plt.text(4, 0.5, "$\lambda=%.2f$" % inflation_factor)
    plt.xlabel('Expected $-log10(p-value)$')
    plt.ylabel('Observed $-log10(p-value)$')
    plt.savefig(output)
    plt.clf()