import pathlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
plt.rc('axes', axisbelow=True)
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['axes.titlesize'] = 20

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
    df = df[df['p_value'] <= 0.01].copy()
    if isinstance(df['CHR'].values[0], str):
        try:
            df['CHR'] = [int(x.replace('Chr', '')) for x in df['CHR']]
        except Exception as exc:
            print("Chromosome identifier might be wrong. Use the chromosome number.")
            print(exc)
    running_pos = 0
    cumul_pos = []
    for chrom, group_df in df.groupby('CHR'):
        cumul_pos.append(group_df['POS'] + running_pos)
        running_pos += group_df['POS'].max()
    df['cumul_pos'] = pd.concat(cumul_pos)

    fig, ax = plt.subplots(1, 1, figsize=(20, 5), constrained_layout=True)
    sns.scatterplot(ax=ax, data=df, x='cumul_pos', y='p_value', hue='CHR', palette='colorblind', linewidth=0, s=20,
                    legend=None)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_yscale("log")
    ax.invert_yaxis()
    ax.minorticks_off()
    ax.set_xlabel('Chromosome')
    ax.set_ylabel(r'$-log_{10}$(p-value)')
    ax.set_xticks(df.groupby('CHR')['cumul_pos'].median())
    ax.set_xticklabels(np.unique(df['CHR']))

    if min_p_values is not None:
        ax.axhline(helper_functions.compute_perm_threshold(min_p_values, sig_level), linewidth=1.5, color='blue',
                   label='permGWAS2')
    ax.axhline(helper_functions.compute_bonf_threshold(n_snps, sig_level), linewidth=1.5, color='red',
               label='Bonferroni')
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.13), fancybox=True, ncol=2, frameon=True)
    fig.savefig(data_dir.joinpath('manhattan_' + pathlib.Path(filename).with_suffix('.png').as_posix()))
    fig.clf()


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
