# perform gwas, use sub functions
import torch
import numpy as np
from itertools import repeat
from perform_gwas import model, kin, perm
import time


def gwas(X, y, K, batchsize, covs):
    """
    perform gwas use subfunctions
    :param X: genotype matrix of shape (n,m)
    :param y: phenotype vector of shape (n)
    :param K: kinship matrix of shape (n,n)
    :param batchsize: number of SNPs to perform GWAS on simultaneously
    :param covs: vector/matrix of covariates of shape (n,c), optional
    :return: p-value, F-score, standard error and effect size of each SNP
    """
    start = time.time()
    n = X.shape[0]
    m = X.shape[1]
    fixed = model.get_fixed_effects(covs, n, y.device)  # shape: (n,1) or (n,c+1)
    K_N = kin.normalize_kinship(K)  # shape: (n,n)
    v_g, v_e = kin.estimate_variance_components(y, K_N, fixed)  # shape: (2)
    C = kin.get_cholesky(v_g * K_N + v_e * torch.eye(n, dtype=K.dtype, device=K.device))  # shape: (n,n)
    y_trans = kin.transform_input(y, C)  # shape: (n)
    fixed_trans = kin.transform_input(fixed, C)  # shape: (n) or (n,c+1)
    RSS_0 = model.get_rss_h0(y_trans, fixed_trans)  # shape: (1)
    freedom_deg = fixed.shape[1] + 1
    tmp = []
    for i in range(int(np.ceil(m / batchsize))):
        lower_bound = i * batchsize
        upper_bound = (i + 1) * batchsize
        if upper_bound > m:
            upper_bound = m
        X_batch = kin.transform_input(X[:, lower_bound:upper_bound].to(C.device), C)  # shape: (n,b)
        X_batch = model.get_x_batch(X_batch, fixed_trans, lower_bound, upper_bound)  # shape: (b,n,2) or (b,n,c+2)
        RSS_1, SE, effSize = model.get_rss_and_se(X_batch, y_trans)
        F_score = model.get_f_score(RSS_0, RSS_1, n, freedom_deg)
        tmp.append(torch.stack((F_score, SE, effSize), dim=1))
        del RSS_1
        del F_score
        del X_batch
        del SE
        del effSize
        torch.cuda.empty_cache()
    output = torch.cat(tmp, dim=0)  # shape(m,3)
    output = output.to(torch.device("cpu"))
    have_test_stats_time = time.time()
    print("Have test statistics of %d SNPs. Elapsed time: %f" % (m, have_test_stats_time-start))
    print("Calculate P-values now")
    p_val = list(map(model.get_p_value, output[:, 0], repeat(n), repeat(freedom_deg)))
    print("Have P-values. Elapsed time: ", time.time()-have_test_stats_time)
    return torch.cat((torch.tensor(p_val).unsqueeze(1), output), dim=1)
# TODO threshold for p-values


def perm_gwas(X, y, K, true_test_stats, perms, batchsize, covs):
    """
    perform gwas with permutations, use subfunctions
    :param X: genotype matrix of shape (n,m)
    :param y: phenotype vector of shape (n)
    :param K: kinship matrix of shape (n,n)
    :param true_test_stats: test statistics of true observations for each SNP
    :param perms: number of permutations
    :param batchsize: number of SNPs to perform GWAS on simultaneously
    :param covs: vector/matrix of covariates of shape (n,c), optional
    :return: adjusted p-value for each SNP, minimal p-value for each permutation, permutation seeds
    """
    start = time.time()
    n = X.shape[0]
    m = X.shape[1]
    fixed = model.get_fixed_effects_perm(covs, n, perms, y.device)  # shape: (p,n,1) or (p,n,c+1)
    K_N = kin.normalize_kinship(K)
    y_perm, my_seeds = perm.permute_phenotype(y, perms)
    y_perm = y_perm.to(y.device)
    var_comps = []
    for i in range(perms):
        var_comps.append(kin.estimate_variance_components(y_perm[:, i], K_N, fixed[0, :, :]))  # shape: (p,2)
    var_comps = torch.tensor(np.array(var_comps), device=K.device)
    C_perm = perm.get_perm_kinships(K_N, var_comps)  # shape: (p,n,n)
    C_perm = kin.get_cholesky(C_perm)  # shape: (p,n,n)
    y_perm = kin.transform_input(torch.unsqueeze(torch.t(y_perm), 2), C_perm)  # shape: (p,n,1)
    fixed = kin.transform_input(fixed, C_perm)  # shape: (p,n,1) or (p,n,c+1)
    RSS0 = model.get_rss_h0_perm(y_perm, fixed)  # shape: (p)
    del K_N
    del var_comps
    del K
    del y
    torch.cuda.empty_cache()
    freedom_deg = fixed.shape[2]+1
    test_stats = []
    for i in range(int(np.ceil(m / batchsize))):
        lower_bound = i * batchsize
        upper_bound = (i + 1) * batchsize
        if upper_bound > m:
            upper_bound = m
        print("Calculate perm test statistics for SNPs %d to %d" % (lower_bound, upper_bound))
        X_batch = kin.transform_input(model.get_v_batch(X[:, lower_bound:upper_bound].to(C_perm.device), perms), C_perm)  # shape: (p,n,b)
        X_batch = model.get_x_batch_perm(X_batch, fixed, lower_bound, upper_bound)  # shape: (p,b,n,2) or (p,b,n,c+2)
        RSS = model.get_rss_perm(X_batch, y_perm)  # shape: (p,b)
        F_score = model.get_f_score(torch.t(RSS0.repeat(upper_bound-lower_bound, 1)), RSS, n, freedom_deg)  # shape: (p,b)
        test_stats.append(F_score)
        del RSS
        del F_score
        del X_batch
        torch.cuda.empty_cache()
    test_stats = torch.cat(test_stats, dim=1)  # shape: (p,m)
    have_test_stats_time = time.time()
    print("Have perm test statistics. Elapsed time: ", have_test_stats_time-start)
    test_stats = test_stats.to(torch.device("cpu"))
    perm_p_val = perm.get_perm_p_value(test_stats, true_test_stats)  # shape: (m)
    min_p_val = perm.get_min_p_value(test_stats, n, freedom_deg)  # shape: (p)
    print("Have adjusted p-values and minimal p-values. Elapsed time: ", time.time()-have_test_stats_time)
    return perm_p_val, min_p_val, my_seeds
