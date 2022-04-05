import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def manhattan(df: pd.DataFrame, name: str, output: str, df_min=None):
    """
    Create manhattan plot with Bonferroni threshold, if df_min is given, additionally with permutation based threshold.
    :param df: pandas DataFrame containing p-values
    :param name: name of column of DataFrame which contains the p-values
    :param output: name to save plot
    :param df_min: pandas DataFrame containing min p-values for permutation based threshold
    :return: save manhattan plot
    """
    data = df[df['p_value'] <= 0.1].copy()
    data['-log10'] = -np.log10(data[name])
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

    g.axes[0, 0].axhline(-np.log10(0.05 / len(data)), color='green', label='Bonferroni')

    if df_min is not None:
        threshold = np.percentile(df_min['min_p_val'], 5)
        g.axes[0, 0].axhline(-np.log10(threshold), color='blue', label='perm threshold')
    g.ax.legend(loc="upper right")
    plt.savefig(output)
    plt.clf()
