# functions for GWAS with EMMAX and EMMAXperm
import argparse

import torch
import scipy.stats as stats
from perform_gwas.kin import transform_input


# functions for EMMAX
def get_fixed_effects(arguments: argparse.Namespace, covs: torch.tensor, n: int):
    """
    check for covariates and create fixed effects vector
    :param arguments: user input
    :param covs: vector or matrix of covariates or None
    :param n: number of samples
    :return: vector of ones if no covs are used, matrix containing column of ones and covs, dim (n,c+1)
    """
    if covs is None:
        return torch.ones((n, 1), dtype=torch.float64, device=arguments.device)
    elif covs.ndim == 1:
        return torch.stack((torch.ones(n, dtype=torch.float64, device=arguments.device), covs), dim=1)
    else:
        return torch.cat((torch.ones((n, 1), dtype=torch.float64, device=arguments.device), covs), dim=1)


def get_v_batch(v: torch.tensor, batchsize: int):
    """
    create copies of input tensor
    :param v: vector of length n or matrix of dim (n,c)
    :param batchsize:
    :return: tensor of copies of v with dim (batchsize,n,1) or (batchsize,n,c)
    """
    if v.ndim == 1:
        return torch.unsqueeze(v.repeat(batchsize, 1), 2)
    if v.ndim == 2:
        return v.repeat(batchsize, 1, 1)


def get_x_batch(X: torch.tensor, fixedEff: torch.tensor, lower_bound: int, upper_bound: int):
    """
    create 3D tensor where each matrix in the tensor contains the same fixed effects and a different SNP
    :param X: nxm genotype matrix
    :param fixedEff: nxc matrix of fixed effects
    :param lower_bound: lower bound of batch
    :param upper_bound: upper bound of batch
    :return: tensor of dim (upper_bound-lower_bound, n, c+1)
    """
    b = get_v_batch(fixedEff, upper_bound-lower_bound)
    Xb = torch.transpose(torch.unsqueeze(X, 0), 0, 2)
    return torch.cat((b, Xb), dim=2)


def get_rss_h0(y: torch.tensor, fixedEff: torch.tensor):
    """
    compute residual sum of squares of H0: marker has no effect on phenotype
    :param y: phenotype vector
    :param fixedEff: vector or matrix of fixed effects
    :return: residual sum of squares
    """
    fixedT = torch.t(fixedEff)
    beta = transform_input(torch.matmul(fixedT, y), torch.matmul(fixedT, fixedEff))
    dif = y - torch.matmul(fixedEff, beta)
    return torch.matmul(torch.t(dif), dif)


def get_rss_and_se(X: torch.tensor, y: torch.tensor):
    """
    in batches:
    First estimate coefficients beta
    Then compute residual sum of squares of alternative hypothesis: marker has effect on phenotype
    Also compute standard error SE
    :param X: fixed effects matrix containing intercept, covariates and several markers in batches
    :param y: phenotype values
    :return: residual sum of squares, standard error and effect size in batches
    """
    y_batch = get_v_batch(y, X.shape[0])
    XT = torch.transpose(X, 1, 2)
    XTX_inv = torch.linalg.pinv(torch.matmul(XT, X), hermitian=True)
    beta = torch.matmul(torch.matmul(XTX_inv, XT), y_batch)
    dif = y_batch-torch.matmul(X, beta)
    diag = torch.diagonal(XTX_inv, dim1=1, dim2=2)[:, 1]
    rss = torch.squeeze(torch.matmul(torch.transpose(dif, 1, 2), dif))
    sigma2 = rss / (y.shape[0] - X.shape[2])
    se = torch.sqrt(sigma2 * diag)
    return rss, se, torch.squeeze(beta[:, -1])


def get_f_score(rss0: torch.tensor, rss1: torch.tensor, n: int, freedom_deg: int):
    """
    compute test statistics in batches
    :param rss0: residual sum of squares of H0: marker has no effect on phenotype
    :param rss1: residual sum of squares of H1: marker has effect on phenotype
    :param n: number of samples
    :param freedom_deg: degrees of freedom
    :return: F1 score
    """
    return (n-freedom_deg)*(rss0-rss1)/rss1


def get_p_value(f1: float, n: int, freedom_deg: int):
    """
    compute p-value using survival function of f distribution
    :param f1: F1 score
    :param n: number of samples
    :param freedom_deg: degrees of freedom
    :return: p-value
    """
    return stats.f.sf(f1, 1, n-freedom_deg)


# functions for EMMAXperm
def get_fixed_effects_perm(arguments: argparse.Namespace, covs: torch.tensor, n: int):
    """
        check for covariates and create fixed effects tensor
        :param arguments: user input
        :param covs: vector or matrix of covariates or None
        :param n: number of samples
        :return: tensor of ones if no covs are used dim (p,n,1),tensor containing ones and covs, dim (p,n,c+1)
        """
    if covs is None:
        return torch.ones((arguments.perm, n, 1), dtype=torch.float64, device=arguments.device)
    else:
        covs_batch = get_v_batch(covs, arguments.perm)
        return torch.cat((torch.ones((arguments.perm, n, 1), dtype=torch.float64, device=arguments.device), covs_batch),
                         dim=2)


def get_v_batch_perm(v: torch.tensor, batchsize: int):
    """
    create copies of input tensor
    :param v: tensor of dim (p,n,1) or (p,n,c)
    :param batchsize:
    :return: tensor of copies of v with dim (p,b,n,1) or (p,b,n,c)
    """
    return torch.transpose(v.repeat(batchsize, 1, 1, 1), dim0=0, dim1=1)


def get_x_batch_perm(X: torch.tensor, fixedEff: torch.tensor, lower: int, upper: int):
    """
    create 4D tensor that contains copies of a 3D tensors that contains matrices with the same fixed effects and a
    different SNP each
    :param X: tensor of dim (p,n,m) containing a matrix X for each permutation
    :param fixedEff: tensor of dim (p,n,1) or (p,n,c)
    :param lower: boundary
    :param upper: boundary
    :return: tensor of dim (p,upper-lower,n,2) or (p,upper-lower,n,c+1)
    """
    b = get_v_batch_perm(fixedEff, upper-lower)
    Xb = torch.unsqueeze(torch.transpose(X, dim0=1, dim1=2), 3)
    return torch.cat((b, Xb), dim=3)


def get_rss_h0_perm(y: torch.tensor, fixedEff: torch.tensor):
    """
    compute residual sum of squares of H0 with permutations
    :param y: phenotype matrix, containing permutations in columns
    :param fixedEff: vector or matrix of fixed effects
    :return: residual sum of squares
    """
    fixedT = torch.transpose(fixedEff, dim0=1, dim1=2)
    beta = transform_input(torch.matmul(fixedT, y), torch.matmul(fixedT, fixedEff))
    dif = y - torch.matmul(fixedEff, beta)
    return torch.squeeze(torch.sum(torch.mul(dif, dif), dim=1))


def get_rss_perm(X: torch.tensor, y: torch.tensor):
    """
    compute residual sum of squares of alternative hypothesis with permutations
    :param X: genotype matrix
    :param y: phenotype vector
    :return: RSS
    """
    y_batch = get_v_batch_perm(y, X.shape[1])
    XT = torch.transpose(X, dim0=2, dim1=3)
    beta = transform_input(torch.matmul(XT, y_batch), torch.matmul(XT, X))
    dif = y_batch-torch.matmul(X, beta)
    return torch.squeeze(torch.sum(torch.mul(dif, dif), dim=2))
