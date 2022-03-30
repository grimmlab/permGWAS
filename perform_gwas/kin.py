# functions for kinship matrix
import torch
from numpy_sugar.linalg import economic_qs
from numpy_sugar import is_all_finite
from glimix_core.lmm import LMM


def estimate_variance_components(y: torch.tensor, K: torch.tensor, fixed_eff: torch.tensor, verbose=False):
    """
    Function to estimate variance components, using the packages numpy_sugar and glimix_core
    :param y: phenotype vector
    :param K: Kinship matrix
    :param fixed_eff: matrix of fixed effects containing intercept and covariates
    :param verbose:
    :return: variance components v_g and v_e
    """
    y_c = y.to(torch.device("cpu"))
    K_c = K.to(torch.device("cpu"))
    fixed_c = fixed_eff.to(torch.device("cpu"))
    if not is_all_finite(y_c):
        raise ValueError("Outcome must have finite values only.")
    if not is_all_finite(K_c):
        raise ValueError("Outcome must have finite values only.")
    if K is not None:
        q_s = economic_qs(K_c)
    else:
        q_s = None

    method = LMM(y_c, fixed_c, q_s, restricted=True)
    method.fit(verbose=verbose)

    v_g = method.scale * (1 - method.delta)
    v_e = method.scale * method.delta

    return v_g, v_e


def get_cholesky(K: torch.tensor):
    """
    compute the cholesky decomposition K = CC^T
    :param K:  matrix
    :return: lower triangle matrix C
    """
    return torch.linalg.cholesky(K)


def transform_input(data: torch.tensor, C: torch.tensor):
    """
    solve: data = Cx
    :param data: input vector
    :param C: matrix
    :return: x
    """
    return torch.linalg.solve(C, data)
