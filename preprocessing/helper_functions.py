import numpy as np
import pathlib
import argparse
import time

def estimate_heritability(v_g: float, v_e: float):
    """
    compute narrow sense heritability
    :param v_g: genetic variance component
    :param v_e: residual variance component
    :return: narrow sense heritability
    """
    return v_g / (v_g + v_e)


def compute_perm_threshold(min_p_val: np.array, sig_level: int):
    """
    Compute permutation-based threshold
    :param min_p_val: array with minimal p-values
    :param sig_level: significance level as percentage value
    :return: threshold
    """
    return np.percentile(min_p_val, sig_level)


def compute_bonf_threshold(number_snps: int, sig_level: int):
    """
    Compute Bonferroni threshold
    :param number_snps: number of SNPs
    :param sig_level: significance level as percentage value
    :return: threshold
    """
    return (sig_level / 100) / number_snps


def print_summary_stats(arguments: argparse.Namespace, samples: int, snps: int, v_g: float, v_e: float,
                       h2: float, bonf1: float, bonf5: float, perm1: float, perm5: float):
    """
    Print summary statistics
    :param arguments: user input arguments
    :param samples: number of samples used
    :param snps: number of SNPs used
    :param v_g: genetic variance component
    :param v_e: residual variiance component
    :param h2: narrow-sense heritability
    :param bonf1: Bonferroni threshold significance level 1%
    :param bonf5: Bonferroni threshold significance level 5%
    :param perm1: permutation-based threshold significance level 1%
    :param perm5: permutation-based threshold significance level 5%
    """
    print('\n')
    print('+++++++++ Summary Statistics +++++++++')
    print('## Genotype file: ' + str(arguments.x))
    print('## Phenotype file: ' + str(arguments.y))
    print('## Phenotype: ' + arguments.y_name)
    if arguments.cov_file is not None:
        print('## Covariate file: ' + str(arguments.cov_file))
    if arguments.k is not None:
        print('## Kinship file: ' + str(arguments.k))
    print('## Number of individuals: ' + str(samples))
    print('## Number of SNPs: ' + str(snps))
    print('## Number of permutations: ' + str(arguments.perm))
    print('## v_g estimate in null model: ' + str(v_g))
    print('## v_e estimate in null model: ' + str(v_e))
    print('## Narrow-sense heritability estimate: ' + str(h2))
    print('## Bonferroni threshold (1% significance level): ' + str(bonf1))
    print('## Bonferroni threshold (5% significance level): ' + str(bonf5))
    if perm1 is not None:
        print('## Permutation-based threshold (1% significance level): ' + str(perm1))
        print('## Permutation-based threshold (5% significance level): ' + str(perm5))
    print('+++++++++++++++++++++++++++')


def write_summary_stats(arguments: argparse.Namespace, samples: int, snps: int, v_g: float, v_e: float,
                       h2: float, bonf1: float, bonf5: float, perm1: float, perm5: float):
    """
    Save summary statistics to txt file
    :param arguments: user input arguments
    :param samples: number of samples used
    :param snps: number of SNPs used
    :param v_g: genetic variance component
    :param v_e: residual variiance component
    :param h2: narrow-sense heritability
    :param bonf1: Bonferroni threshold significance level 1%
    :param bonf5: Bonferroni threshold significance level 5%
    :param perm1: permutation-based threshold significance level 1%
    :param perm5: permutation-based threshold significance level 5%
    """
    filename = arguments.out_dir.joinpath('summary_statistics_' +
                                          pathlib.Path(arguments.out_file).with_suffix('.txt').as_posix())
    with open(filename, 'w') as f:
        f.write('Summary Statistics:\n')
        f.write('## Genotype file:\t' + str(arguments.x) + '\n')
        f.write('## Phenotype file:\t' + str(arguments.y) + '\n')
        f.write('## Phenotype:\t' + arguments.y_name + '\n')
        if arguments.cov_file is not None:
            f.write('## Covariate file: ' + str(arguments.cov_file) + '\n')
        if arguments.k is not None:
            f.write('## Kinship file: ' + str(arguments.k) + '\n')
        f.write('## Number of individuals:\t' + str(samples) + '\n')
        f.write('## Number of SNPs:\t' + str(snps) + '\n')
        f.write('## Number of permutations:\t' + str(arguments.perm) + '\n')
        f.write('## v_g estimate in null model:\t' + str(v_g) + '\n')
        f.write('## v_e estimate in null model:\t' + str(v_e) + '\n')
        f.write('## Narrow-sense heritability estimate:\t' + str(h2) + '\n')
        f.write('## Bonferroni threshold (1% significance level):\t' + str(bonf1) + '\n')
        f.write('## Bonferroni threshold (5% significance level):\t' + str(bonf5) + '\n')
        if perm1 is not None:
            f.write('## Permutation-based threshold (1% significance level):\t' + str(perm1) + '\n')
            f.write('## Permutation-based threshold (5% significance level):\t' + str(perm5))


def get_summary_stats(arguments: argparse.Namespace, samples: int, snps: int, v_g: float, v_e: float,
                      min_p_val: np.array):
    """
    Compute summary statistics, print and save them to file
    :param arguments: user input arguments
    :param samples: number of samples used
    :param snps: number of SNPs used
    :param v_g: genetic variance component
    :param v_e: residual variiance component
    :param min_p_val: minimal p-values
    """
    h2 = estimate_heritability(v_g=v_g, v_e=v_e)
    bonf1 = compute_bonf_threshold(number_snps=snps, sig_level=1)
    bonf5 = compute_bonf_threshold(number_snps=snps, sig_level=5)
    if min_p_val is not None:
        perm1 = compute_perm_threshold(min_p_val=min_p_val, sig_level=1)
        perm5 = compute_perm_threshold(min_p_val=min_p_val, sig_level=5)
    else:
        perm1 = None
        perm5 = None

    write_summary_stats(arguments=arguments, samples=samples, snps=snps, v_g=v_g, v_e=v_e, h2=h2,
                        bonf1=bonf1, bonf5=bonf5, perm1=perm1, perm5=perm5)
    print_summary_stats(arguments=arguments, samples=samples, snps=snps, v_g=v_g, v_e=v_e, h2=h2,
                        bonf1=bonf1, bonf5=bonf5, perm1=perm1, perm5=perm5)