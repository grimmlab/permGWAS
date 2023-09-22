# perform gwas, use sub functions
import argparse
import torch
import numpy as np
from itertools import repeat
from perform_gwas import model, kin, perm
from preprocessing import load_files, prepare_data
import time


def gwas(arguments: argparse.Namespace, X: torch.tensor, y: torch.tensor, K: torch.tensor, covs: torch.tensor,
         X_index: np.array, m: int):
    """
    Function to perform batch-wise computation of univariate test:
        (1) estimate variance components
        (2) calculate Cholesky decomposition
        (3) linearly transform data
        (4) calculate residual sum of squares of null model
        (5) batch-wise:
            (a) linearly transform marker
            (b) calculate effect size, residual sum of squares and standard error
            (c) calculate test statistic
        (6) calculate p-values
    :param arguments: user input
    :param X: genotype matrix of shape (n,m)
    :param y: phenotype vector of shape (n)
    :param K: kinship matrix of shape (n,n)
    :param covs: vector/matrix of covariates of shape (n,c), optional
    :param X_index: indices of genotype matrix samples to load in batches
    :param m: total number of SNPs to work on
    :return: p-value, test statistic, standard error and effect size of each SNP, variance components and minor allele
    frequencies if X was loaded batch-wise
    """
    start = time.time()
    n = y.shape[0]
    fixed = model.get_fixed_effects(arguments, covs, n)  # shape: (n,1) or (n,c+1)
    # estimate variance components
    v_g, v_e = kin.estimate_variance_components(y, K, fixed)  # shape: (2)
    # Cholesky decomposition
    C = kin.get_cholesky(v_g * K + v_e * torch.eye(n, dtype=K.dtype, device=arguments.device))  # shape: (n,n)
    # linearly transform data, i.e. solve y = Cb
    y_trans = kin.transform_input(y, C)  # shape: (n)
    fixed = kin.transform_input(fixed, C)  # shape: (n) or (n,c+1)
    # calculate residual sum of squares of null model
    RSS_0 = model.get_rss_h0(y_trans, fixed)  # shape: (1)
    freedom_deg = fixed.shape[1] + 1
    tmp = []
    freq = []
    # in batches:
    for i in range(int(np.ceil(m / arguments.batch))):
        lower_bound = i * arguments.batch
        upper_bound = (i + 1) * arguments.batch
        if upper_bound > m:
            upper_bound = m
        if X is None:
            # load X batch-wise
            X_batch = load_files.load_genotype_matrix(arguments.x, sample_index=X_index, snp_lower_index=lower_bound,
                                                      snp_upper_index=upper_bound, not_add=arguments.not_add)
            # calculate minor allele frequencies
            freq.append(prepare_data.get_maf(X_batch))
            # linearly transform X
            X_batch = kin.transform_input(X_batch.to(arguments.device), C)  # shape: (n,b)
        else:
            # linearly transform X_batch if X was completely loaded before
            X_batch = kin.transform_input(X[:, lower_bound:upper_bound].to(arguments.device), C)  # shape: (n,b)
        X_batch = model.get_x_batch(X_batch, fixed, lower_bound, upper_bound)  # shape: (b,n,2) or (b,n,c+2)
        # calculate effect size, residual sum of squares and standard error
        RSS_1, SE, effSize = model.get_rss_and_se(X_batch, y_trans)
        # calculate test statistic
        F_score = model.get_f_score(RSS_0, RSS_1, n, freedom_deg)
        tmp.append(torch.stack((F_score, SE, effSize), dim=1).to(torch.device("cpu")))
        if arguments.device.type != "cpu":
            with torch.cuda.device(arguments.device):
                del RSS_1
                del F_score
                del X_batch
                del SE
                del effSize
                torch.cuda.empty_cache()
    output = torch.cat(tmp, dim=0)  # shape(m,3)
    output = output.to(torch.device("cpu"))
    if X is None:
        freq = torch.cat(freq, dim=0)
    time_test_stats = time.time()
    print("Have test statistics of %d SNPs. Elapsed time: %f" % (m, time_test_stats-start))
    print("Calculate P-values now")
    # compute p-values
    p_val = list(map(model.get_p_value, output[:, 0], repeat(n), repeat(freedom_deg)))
    print("Have P-values. Elapsed time: ", time.time()-time_test_stats)
    return torch.cat((torch.tensor(p_val).unsqueeze(1), output), dim=1), v_g, v_e, freq


def perm_gwas(arguments: argparse.Namespace, X: torch.tensor, y: torch.tensor, K: torch.tensor,
              true_test_stats: torch.tensor, covs: torch.tensor, X_index: np.array, m: int):
    """
    Function to perform batch-wise computation of permutation-based test:
        (1) compute permutations of phenotype
        (2) estimate variance components for each permutation
        (3) calculate Cholesky decomposition for each permutation
        (4) linearly transform data
        (5) calculate residual sum of squares of null model for each permutation
        (6) batch-wise:
            (a) linearly transform marker
            (b) calculate residual sum of squares
            (c) calculate test statistic
        (7) calculate permutation-based p-values
        (8) calculate Westfall-Young permutation-based threshold
    :param arguments: user input
    :param X: genotype matrix of shape (n,m)
    :param y: phenotype vector of shape (n)
    :param K: kinship matrix of shape (n,n)
    :param true_test_stats: test statistics of true observations for each SNP
    :param covs: vector/matrix of covariates of shape (n,c), optional
    :param X_index: indices of genotype matrix samples to load in batches
    :param m: total number of SNPs to work on
    :return: adjusted p-value for each SNP, minimal p-value for each permutation, permutation seeds
    """
    start = time.time()
    n = y.shape[0]
    fixed = model.get_fixed_effects_perm(arguments, covs, n)  # shape: (p,n,1) or (p,n,c+1)
    # compute permutations of phenotype
    y_perm, my_seeds = perm.permute_phenotype(y, arguments.perm)
    y_perm = y_perm.to(arguments.device)  # shape: (n,p)
    var_comps = []
    # estimate variance components for each permutation
    for i in range(arguments.perm):
        var_comps.append(kin.estimate_variance_components(y_perm[:, i], K, fixed[0, :, :]))  # shape: (p,2)
    var_comps = torch.tensor(np.array(var_comps), device=arguments.device)
    # calculate Cholesky decomposition for each permutation
    C_perm = perm.get_perm_kinships(K, var_comps)  # shape: (p,n,n)
    C_perm = kin.get_cholesky(C_perm)  # shape: (p,n,n)
    # linearly transform data
    y_perm = kin.transform_input(torch.unsqueeze(torch.t(y_perm), 2), C_perm)  # shape: (p,n,1)
    fixed = kin.transform_input(fixed, C_perm)  # shape: (p,n,1) or (p,n,c+1)
    # calculate residual sum of squares of null model for each permutation
    RSS0 = model.get_rss_h0_perm(y_perm, fixed)  # shape: (p)
    if arguments.device.type != "cpu":
        with torch.cuda.device(arguments.device):
            del var_comps
            del K
            del y
            del covs
            torch.cuda.empty_cache()
    freedom_deg = fixed.shape[2]+1
    test_stats = []
    # in batches:
    for i in range(int(np.ceil(m / arguments.batch_perm))):
        lower_bound = i * arguments.batch_perm
        upper_bound = (i + 1) * arguments.batch_perm
        if upper_bound > m:
            upper_bound = m
        print("Calculate perm test statistics for SNPs %d to %d" % (lower_bound, upper_bound))
        if X is None:
            # load X batch-wise and linearly transform it
            X_batch = load_files.load_genotype_matrix(arguments.x, sample_index=X_index, snp_lower_index=lower_bound,
                                                      snp_upper_index=upper_bound, not_add=arguments.not_add)
            X_batch = kin.transform_input(model.get_v_batch(X_batch.to(arguments.device),
                                                            arguments.perm), C_perm)  # shape: (p,n,b)
        else:
            # if X was already loaded completely, linearly transform it
            X_batch = kin.transform_input(model.get_v_batch(X[:, lower_bound:upper_bound].to(arguments.device),
                                                            arguments.perm), C_perm)  # shape: (p,n,b)
        X_batch = model.get_x_batch_perm(X_batch, fixed, lower_bound, upper_bound)  # shape: (p,b,n,2) or (p,b,n,c+2)
        # calculate residual sum of squares
        RSS = model.get_rss_perm(X_batch, y_perm)  # shape: (p,b)
        # calculate test statistics
        F_score = model.get_f_score(torch.t(RSS0.repeat(upper_bound-lower_bound, 1)), RSS, n, freedom_deg)  # shape: (p,b)
        test_stats.append(F_score.to(torch.device("cpu")))
        if arguments.device.type != "cpu":
            with torch.cuda.device(arguments.device):
                del RSS
                del F_score
                del X_batch
                torch.cuda.empty_cache()
    test_stats = torch.cat(test_stats, dim=1)  # shape: (p,m)
    time_test_stats = time.time()
    print("Have perm test statistics. Elapsed time: ", time_test_stats-start)
    test_stats = test_stats.to(torch.device("cpu"))
    # calculate permutation-based p-values
    perm_p_val = perm.get_perm_p_value(test_stats, true_test_stats)  # shape: (m)
    # calculate Westfall-Young permutation-based threshold
    min_p_val = perm.get_min_p_value(test_stats, n, freedom_deg)  # shape: (p)
    print("Have adjusted p-values and minimal p-values. Elapsed time: ", time.time()-time_test_stats)
    return perm_p_val, min_p_val, my_seeds
