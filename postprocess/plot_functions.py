import pathlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats

from utils import helper_functions


def manhattan_plot(df: pd.DataFrame, data_dir: pathlib.Path, filename: str, min_p_values: np.array = None,
                   sig_level: int = 5):
    """
    Save Manhattan plot as manhattan_FILENAME.png to data_dir

    :param df: DataFrame containing chromosome (CHR) and position (POS) identifiers, and corresponding p_values
    :param data_dir: full path to save directory
    :param filename: name of file
    :param min_p_values: array containing minimal p_values to compute permutation-based threshold
    :param sig_level: significance level for Bonferroni and perm thresholds, default is 5
    """
    if not {'CHR', 'POS', 'p_value'}.issubset(df.columns):
        raise Exception('Cannot create Manhattan plot; need CHR, POS and p_value in DataFrame.')
    n_snps = len(df)
    df = df[df['p_value'] <= 0.1].copy()
    df['-log10'] = -np.log10(df['p_value'])
    running_pos = 0
    cumul_pos = []
    for chrom, group_df in df.groupby('CHR'):
        cumul_pos.append(group_df['POS'] + running_pos)
        running_pos += group_df['POS'].max()
    df['cumul_pos'] = pd.concat(cumul_pos)

    g = sns.relplot(
        data=df,
        x='cumul_pos',
        y='-log10',
        aspect=4,
        hue='CHR',
        palette='Set1',
        linewidth=0,
        s=6,
        legend=None)

    g.ax.set_xlabel('Chromosome')
    g.ax.set_ylabel('-log10(p-value)')
    g.ax.set_xticks(df.groupby('CHR')['cumul_pos'].median())
    labels = np.unique(df['CHR'])
    g.ax.set_xticklabels(labels)

    g.axes[0, 0].axhline(-np.log10(helper_functions.compute_bonf_threshold(n_snps, sig_level)),
                         color='green', label='Bonferroni')

    if min_p_values is not None:
        g.axes[0, 0].axhline(-np.log10(helper_functions.compute_perm_threshold(min_p_values, sig_level)),
                             color='blue', label='perm threshold')
    g.ax.legend(loc="upper right")
    plt.savefig(data_dir.joinpath('manhattan_' + pathlib.Path(filename).with_suffix('.png').as_posix()))
    plt.clf()


def qq_plot(p_values: np.array, data_dir: pathlib.Path, filename: str):
    """
    Save QQ-plot as qq_plot_FILENAME.png to data_dir

    :param p_values: array containing p_values
    :param data_dir: full path to save directory
    :param filename: name of file
    """
    n_snps = len(p_values)
    observed_p = -np.log10(np.sort(p_values))
    expected_p = -np.log10(np.arange(1.0 / float(n_snps), 1, 1.0 / float(n_snps + 1)))
    inflation_factor = np.median(stats.chi2.isf(p_values, 1)) / 0.456

    plt.figure(figsize=(4, 4))
    plt.plot(expected_p, observed_p, '.', markersize=4, markeredgewidth=0, alpha=0.8)
    plt.plot(expected_p, expected_p, 'k--', linewidth=0.75)
    plt.text(3.5, 0.5, "$\lambda=%.2f$" % inflation_factor)
    plt.xlabel('Expected $-log10(p-value)$')
    plt.ylabel('Observed $-log10(p-value)$')
    plt.savefig(data_dir.joinpath('qq_plot_' + pathlib.Path(filename).with_suffix('.png').as_posix()))
    plt.clf()
