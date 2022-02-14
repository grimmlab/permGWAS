import argparse

import numpy as np
import torch
import pandas as pd
import h5py
from pandas_plink import read_plink1_bin


def load_genotype(arguments: argparse.Namespace):
    """
    load genotype takes PLink files, binary PLINK files, .csv and .h5, .hdf5, .h5py files
    for binary PLINK files: need either "NAME.bed", "NAME.bim" or "NAME.fam", all need to be in same folder
    :param arguments: user input
    :return: genotype in additive encoding, with sample ids and SNP positions
    """
    suffix = arguments.x.suffix
    if suffix in ('.bed', '.bim', '.fam'):
        x_file = arguments.x.with_suffix('').as_posix()
        gt = read_plink1_bin(x_file + '.bed', x_file + '.bim',
                             x_file + '.fam', ref="a0", verbose=False)
        sample_ids = np.array(gt['fid'], dtype=np.int).flatten()
        positions = np.array(gt['pos']).flatten()
        chromosomes = np.array(gt['chrom']).flatten()
        X = torch.tensor(gt.values, dtype=torch.float64)
    elif suffix in ('.h5', '.hdf5', '.h5py'):
        with h5py.File(arguments.x, "r") as gt:
            chromosomes = gt['chr_index'][:].astype(str)
            positions = gt['position_index'][:].astype(int)
            sample_ids = gt['sample_ids'][:].astype(int)
            X = torch.tensor(gt['snps'][:], dtype=torch.float64)
    elif suffix == '.csv':
        gt = pd.read_csv(arguments.x, index_col=0)
        identifiers = np.array(list(map(lambda a: a.split("_"), gt.columns.values)))
        chromosomes = identifiers[:, 0]
        positions = identifiers[:, 1]
        sample_ids = np.asarray(gt.index)
        X = torch.tensor(gt.values, dtype=torch.float64)
    elif suffix in ('map', 'ped'):
        x_file = arguments.x.with_suffix('').as_posix()
        with open(x_file + '.map', 'r') as f:
            chromosomes = []
            positions = []
            for line in f:
                tmp = line.strip().split(" ")
                chromosomes.append(tmp[0].strip())
                positions.append(tmp[-1].strip())
        chromosomes = np.array(chromosomes)
        positions = np.array(positions)
        iupac_map = {"AA": "A", "GG": "G", "TT": "T", "CC": "C", "AG": "R", "GA": "R", "CT": "Y", "TC": "Y", "GC": "S",
                     "CG": "S", "AT": "W", "TA": "W", "GT": "K", "TG": "K", "AC": "M", "CA": "M"}
        with open(x_file + '.ped', 'r') as f:
            sample_ids = []
            raw = []
            for line in f:
                tmp = line.strip().split(" ")
                sample_ids.append(int(tmp[1].strip()))
                snps = []
                j = 6
                while j < len(tmp) - 1:
                    snps.append(iupac_map[tmp[j] + tmp[j + 1]])
                    j += 2
                raw.append(snps)
        sample_ids = np.array(sample_ids)
        raw = np.array(raw)
        X = encode_homozygous(raw)
        # TODO encode heterozygous
    else:
        raise NotImplementedError('Only accept .h5, .hdf5, .h5py, .csv and binary PLINK genotype files')
    return X, sample_ids, positions, chromosomes


def encode_homozygous(matrix: np.array):
    """
    :param matrix:
    :return: get additive encoding of genotype matrix
    """
    maj_min = []
    index_arr = []
    for col in np.transpose(matrix):
        _, inv, counts = np.unique(col, return_counts=True, return_inverse=True)
        tmp = np.where(counts == np.max(counts), 0., 2.)
        maj_min.append(tmp)
        index_arr.append(inv)
    maj_min = np.transpose(np.array(maj_min))
    ind_arr = np.transpose(np.array(index_arr))
    cols = np.arange(maj_min.shape[1])
    return torch.tensor(maj_min[ind_arr, cols])


def load_phenotype(arguments: argparse.Namespace):
    """
    load phenotype
    Accept .csv, .pheno, .txt files. For .txt files assume that separator is a single space. First column should contain
     sample_ids. Name of phenotype should be column name
    :param arguments: user input
    :return: pandas DataFrame with sample ids as index
    """
    suffix = arguments.y.suffix
    if suffix == ".csv":
        y = pd.read_csv(arguments.y)
    elif suffix in (".pheno", ".txt"):
        y = pd.read_csv(arguments.y, sep=" ")
    else:
        raise NotImplementedError('Only accept .csv, .pheno and .txt phenotype files')
    y = y.sort_values(y.columns[0]).groupby(y.columns[0]).mean()
    if arguments.y_name not in y.columns:
        raise Exception('Phenotype ' + arguments.y_name + ' is not in phenotype file ' + str(arguments.y))
    else:
        y = y[[arguments.y_name]].dropna()
    return y


def load_covariates(arguments: argparse.Namespace):
    """
    Only take csv-files: sample ids have to be in first column
    :param arguments: user input
    :return: pandas DataFrame with sample ids as index
    """
    if arguments.cov_file.suffix == ".csv":
        covs = pd.read_csv(arguments.cov_file)
        covs = covs.sort_values(covs.columns[0]).groupby(covs.columns[0]).mean().dropna()
    else:
        raise NotImplementedError('Only accept .csv covariates files')
    return covs


def load_kinship(arguments: argparse.Namespace):
    """
    load kinship matrix fom file. Only take .csv files. Sample ids have to be in first column
    :param arguments: user input
    :return: kinship matrix and sample ids
    """
    if arguments.k.suffix == ".csv":
        kin = pd.read_csv(arguments.k, index_col=0)
        K = torch.tensor(kin.values)
        sample_ids = np.array(kin.index)
    else:
        raise NotImplementedError('Only accept .csv kinship files')
    return K, sample_ids
