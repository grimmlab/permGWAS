# functions for kinship matrix
import torch
from numpy_sugar.linalg import economic_qs
from numpy_sugar import is_all_finite
from glimix_core.lmm import LMM


def normalize_kinship(K):
    """
    normalize K using a Gower's centered matrix
    :param K: kinship matrix
    :return: normalized kinship matrix
    """
    n = K.shape[0]
    P = (torch.eye(n, dtype=K.dtype, device=K.device) - torch.ones(n, n, dtype=K.dtype, device=K.device)/n)
    return (n-1)/torch.sum(torch.mul(P, K))*K


# TODO without limix
def estimate_variance_components(y, K, fixed_eff, verbose=False):
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


def get_cholesky(K):
    """
    compute the cholesky decomposition K = CC^T
    :param K:  matrix
    :return: lower triangle matrix C
    """
    return torch.linalg.cholesky(K)


def transform_input(data, C):
    """
    solve: data = Cx
    :param data: input vector
    :param C: matrix
    :return: x
    """
    return torch.linalg.solve(C, data)
