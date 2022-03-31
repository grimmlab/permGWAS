# functions to load and match genotype and phenotype
import argparse
import sys
import numpy as np
import torch
from preprocessing import load_files


def filter_non_informative_snps(X: torch.tensor, pos: np.array, chrom: np.array):
    """
    Remove constant SNPs from genotype matrix
    :param X: genotype matrix
    :param pos: vector containing positions
    :param chrom: vector containing chromosomes
    :return: X, positions, chromosomes
    """
    tmp = np.where(X.std(axis=0) == 0)[0]
    X = np.delete(X, tmp, axis=1)
    pos = np.delete(pos, tmp, axis=0)
    chrom = np.delete(chrom, tmp, axis=0)
    return X, pos, chrom


def standardize_matrix(matrix: torch.tensor):
    """
    :param matrix: genotype matrix
    :return: matrix with zero mean and unit variance
    """
    return (matrix-matrix.mean(axis=0))/matrix.std(axis=0)


def get_kinship(X: torch.tensor):
    """
    compute realized relationship matrix
    :param X: genotype matrix in additive encoding
    :return: kinship matrix
    """
    X_stand = standardize_matrix(X)
    K = torch.matmul(X_stand, torch.t(X_stand))/X.shape[1]
    # set negative values in K to zero
    return torch.where(K > 0, K, 0.)


def normalize_kinship(K: torch.tensor):
    """
    normalize K using a Gower's centered matrix
    :param K: kinship matrix
    :return: normalized kinship matrix
    """
    n = K.shape[0]
    P = (torch.eye(n, dtype=K.dtype, device=K.device) - torch.ones(n, n, dtype=K.dtype, device=K.device)/n)
    return (n-1)/torch.sum(torch.mul(P, K))*K


def load_and_prepare_data(arguments: argparse.Namespace):
    """
    load and match genotype and phenotype and covariates, load/create kinship matrix
    If kinship is provided and genotype file is in (.h5, .h5py, .hdf5) it is possible to load genotype matrix batch wise
    during computations to save memory.
    :param arguments: user input
    :return: genotype matrix, phenotype vector, kinship matrix, vector with SNP positions and corresponding chromosomes
    """
    # load and match genotype ids and phenotype
    sample_ids, pos, chrom = load_files.load_genotype_ids(arguments.x)
    y = load_files.load_phenotype(arguments)
    y_ids = np.asarray(y.index, dtype=sample_ids.dtype).flatten()
    sample_index = (np.reshape(y_ids, (y_ids.shape[0], 1)) == sample_ids).nonzero()
    if len(sample_index[0]) == 0:
        print("Samples of genotype and phenotype do not match.")
        sys.exit()
    y = torch.tensor(y.values, dtype=torch.float64).flatten()[sample_index[0]]

    # check if kinship is provided if not, load genotype and create kinship
    if arguments.k is None:
        X = load_files.load_genotype_matrix(arguments.x, sample_index=sample_index[1])
        X, pos, chrom = filter_non_informative_snps(X, pos, chrom)
        K = get_kinship(X)
    # load genotype if needed
    else:
        if arguments.load_genotype or (arguments.x.suffix not in ('.h5', '.hdf5', '.h5py')):
            X = load_files.load_genotype_matrix(arguments.x, sample_index=sample_index[1])
            X, pos, chrom = filter_non_informative_snps(X, pos, chrom)
        else:
            X = None
        # load kinship from file if needed
        K, K_ids = load_files.load_kinship(arguments)
        K_index = (np.reshape(sample_ids[sample_index[1]],
                              (sample_ids[sample_index[1]].shape[0], 1)) == K_ids).nonzero()
        if len(K_index[1]) == len(sample_index[1]):
            K = K[K_index[1], :][:, K_index[1]]
        else:
            print("Sample ids of genotype and kinship matrix do not match.")
            sys.exit()
    # normalize kinship matrix
    K = normalize_kinship(K)
    # load and match covariates if provided
    if arguments.cov_file is None:
        covs = None
    else:
        covs = load_files.load_covariates(arguments)
        covs_ids = np.asarray(covs.index, dtype=y_ids.dtype).flatten()
        covs_index = (np.reshape(y_ids[sample_index[0]], (y_ids[sample_index[0]].shape[0], 1)) == covs_ids).nonzero()
        if len(covs_index[1]) == len(sample_index[0]):
            covs = torch.tensor(covs.values, dtype=torch.float64).flatten()[covs_index[1]]
        else:
            print('Sample ids of covariates and phenotype do not match.')
            sys.exit()
    return X, y, K, covs, pos, chrom, sample_index[1]


def get_maf(X: torch.tensor):
    """
    Function to calculate minor allele frequencies of each SNP
    :param X: genotype matrix
    :return: vector containing frequencies
    """
    freq = (torch.sum(X, 0)) / (2 * X.shape[0])
    return freq


def use_maf_filter(X: torch.tensor, positions: np.array, chrom: np.array, maf: int):
    """
    filter genotype by minor allele frequency
    :param X: genotype matrix
    :param positions: SNP positions
    :param chrom: SNP chromosomes
    :param maf: maf threshold
    :return: filtered genotype matrix and vectors with SNP positions and chromosomes
    """
    freq = get_maf(X)
    tmp = np.where(freq <= maf/100)[0]
    X = np.delete(X, tmp, axis=1)
    positions = np.delete(positions, tmp, axis=0)
    chrom = np.delete(chrom, tmp, axis=0)
    freq = np.delete(freq, tmp, axis=0)
    return X, positions, chrom, freq
