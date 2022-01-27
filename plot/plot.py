import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def manhattan(df, name, output, df_min=None):
    df['-log10'] = -np.log10(df[name])
    running_pos = 0
    cumul_pos = []
    for chrom, group_df in df.groupby('CHR'):
        cumul_pos.append(group_df['POS'] + running_pos)
        running_pos += group_df['POS'].max()
    df['Cumul_Pos'] = pd.concat(cumul_pos)

    g = sns.relplot(
        data=df,
        x='Cumul_Pos',
        y='-log10',
        aspect=4,
        hue='CHR',
        palette='Set1',
        linewidth=0,
        s=6,
        legend=None)

    g.ax.set_xlabel('Position')
    g.ax.set_ylabel('-log10(p-value)')
    g.ax.set_xticks(df.groupby('CHR')['Cumul_Pos'].median())
    g.ax.set_xticklabels([])

    g.axes[0, 0].axhline(-np.log10(0.05 / len(df)), color='green', label='Bonferroni')

    if df_min is not None:
        threshold = np.percentile(df_min['min_p_val'], 5)
        g.axes[0, 0].axhline(-np.log10(threshold), color='blue', label='perm threshold')
    g.ax.legend(loc="upper right")
    plt.savefig(output)
    plt.clf()


def create_manhattan(df, arg, df_min=None):
    """
    Create manhattan plot with Bonferroni threshold, if df_min is given, additionally with permutation based threshold.
    if perm>0, additionally create manhattan plot for adjusted p-values with Bonferroni threshold
    :param df: pandas DataFrame containing p-values
    :param arg: argparse.Namespace
    :param df_min: pandas DataFrame containing min p-values for permutation based threshold
    :return: save manhattan plot
    """
    output = arg.out_dir.joinpath(arg.out_file + '_manhattan.png')
    manhattan(df, 'p_value', output, df_min)
    if arg.perm is not None:
        output = arg.out_dir.joinpath(arg.out_file + '_manhattan_adjusted.png')
        manhattan(df, 'adjusted_p_val', output)
