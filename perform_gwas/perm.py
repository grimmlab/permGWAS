# functions for permutations
import numpy as np
import torch
from perform_gwas.model import get_p_value


def permute_phenotype(y, perm):
    """
    compute permutations of phenotype vector
    :param y: phenotype vector
    :param perm: number of permutations
    :return: matrix with different permutations as columns and vector with permutation seeds
    """
    y = y.to(torch.device("cpu"))
    y_perm = []
    perm_seeds = []
    for j in range(perm):
        my_seed = np.ndarray.item(np.random.randint(9999, size=1))
        perm_seeds.append(my_seed)
        rng = np.random.default_rng(my_seed)
        y_perm.append(rng.permutation(y))
    y_perm = torch.tensor(np.array(y_perm))
    return torch.t(y_perm), perm_seeds


def get_k_batch(v, K):
    """
    :param v: vector of variance components of length p (for each permutation)
    :param K: matrix of dim (n,n)
    :return: batch of matrices v*K of dim (p,n,n)
    """
    return v.unsqueeze(1).unsqueeze(2)*K.repeat(v.shape[0], 1, 1)


def get_perm_kinships(K, var_comps):
    """
    :param K: kinship matrix of dim (n,n)
    :param var_comps: matrix with variance components v_g in first and v_e in second col, dim (p,2)
    :return: batch of matrices v_g*K+v_e*I of dim (p,n,n)
    """
    n = K.shape[0]
    return get_k_batch(var_comps[:, 0], K)+get_k_batch(var_comps[:, 1], torch.eye(n, dtype=K.dtype, device=K.device))


def get_perm_p_value(perm_test_stats, true_test_stats):
    """
    :param perm_test_stats: matrix containing test-statistics for all permutations and SNPs, dim (p,m)
    :param true_test_stats: vector containing test-statistics of true observations for all SNPs, length m
    :return: adjusted p-values
    """
    sorted_test_stats, ind = torch.sort(perm_test_stats.flatten())
    n = sorted_test_stats.shape[0]
    test_stats_ind = torch.searchsorted(sorted_test_stats.contiguous(), true_test_stats.contiguous(), right=True)
    if test_stats_ind == n:
        adj_p_value = 1 / n
    else:
        adj_p_value = (n - test_stats_ind) / n
    return adj_p_value


def get_min_p_value(test_stats, n, freedom_deg):
    """
    :param test_stats: matrix containing test-statistics for all permutations and SNPs, dim (p,m)
    :param n: number of samples
    :return: vector containing the minimal p-value for each permutation
    """
    max_test_stats, _ = torch.max(test_stats, dim=1)
    min_p_val = []
    for test in max_test_stats:
        min_p_val.append(get_p_value(test, n, freedom_deg))
    return torch.tensor(min_p_val)
