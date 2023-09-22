import argparse

import numpy as np
import torch
import pandas as pd
import h5py
from pandas_plink import read_plink1_bin


def load_genotype_ids(genotype_file):
    """
    Function to load sample ids and chromosome and position identifiers of genotype.
    Takes PLINK files, binary PLINK files, CSV and H5, HDF5, H5PY files as input
    for binary PLINK files: need either "NAME.bed", "NAME.bim" or "NAME.fam", all need to be in same folder
    Raises Exception if not the right file type
    H5, HDF5, H5PY files need to have the following structure:
        snps:           genotype values in additive encoding (012) with samples as rows and markers as columns
        sample_ids:     sample identifier in the same order as the rows of the genotype matrix
        chr_index:      chromosome identifier in the same order as the columns of the genotype matrix
        position_index: position number (integer) in the same order as the columns of the genotype matrix
    CSV files need to have the following structure:
        sample identifier in the first column
        SNP identifier in the first row need to be formatted as CHROMOSOME_POSITION
        genotype values in additive encoding in the remaining rows and columns
    For PLINK binary files GENOTYPE.bed, GENOTYPE.bim, GENOTYPE.fam need to be in same folder
    For PLINK files GENOTYPE.ped and GENOTYPE.map files need to be in same folder
    :param genotype_file: path to genotype file
    :return: sample ids, SNP positions and chromosomes
    """
    suffix = genotype_file.suffix
    # load binary PLINK
    if suffix in ('.bed', '.bim', '.fam'):
        prefix = genotype_file.with_suffix('').as_posix()
        gt = read_plink1_bin(prefix + '.bed', prefix + '.bim', prefix + '.fam', ref="a0", verbose=False)
        sample_ids = np.array(gt['fid'], dtype=str).flatten()
        positions = np.array(gt['pos']).flatten()
        chromosomes = np.array(gt['chrom']).flatten()
    # load H5, HDF5, H5PY
    elif suffix in ('.h5', '.hdf5', '.h5py'):
        with h5py.File(genotype_file, "r") as gt:
            chromosomes = gt['chr_index'][:].astype(str)
            positions = gt['position_index'][:].astype(int)
            sample_ids = gt['sample_ids'][:].astype(str)
    # load CSV
    elif suffix == '.csv':
        gt = pd.read_csv(genotype_file, index_col=0)
        identifiers = np.array(list(map(lambda a: a.split("_"), gt.columns.values)))
        chromosomes = identifiers[:, 0]
        positions = identifiers[:, 1].astype(int)
        sample_ids = np.asarray(gt.index, dtype=str)
    # load PLINK
    elif suffix in ('.map', '.ped'):
        prefix = genotype_file.with_suffix('').as_posix()
        with open(prefix + '.map', 'r') as f:
            chromosomes = []
            positions = []
            for line in f:
                tmp = line.strip().split(" ")
                chromosomes.append(tmp[0].strip())
                positions.append(tmp[-1].strip())
        chromosomes = np.array(chromosomes)
        positions = np.array(positions)
        with open(prefix + '.ped', 'r') as f:
            sample_ids = []
            for line in f:
                tmp = line.strip().split(" ")
                sample_ids.append(int(tmp[1].strip()))
        sample_ids = np.array(sample_ids, dtype=str)
    else:
        raise NotImplementedError('Only accept H5, HDF5, H5PY, CSV, PLINK and binary PLINK genotype files')
    return sample_ids, positions, chromosomes


def load_genotype_matrix(genotype_file, sample_index=None, snp_lower_index=None, snp_upper_index=None, not_add=False):
    """
    Load genotype matrix. Accepts PLINK files, binary PLINK files, CSV and H5, HDF5, H5PY files.
    For H5, HDF5, H5PY files it is possible to only load certain samples and SNPs batch wise.
    For all other files genotype will be loaded completely but only needed samples will be returned
    H5, HDF5, H5PY files need to have the following structure:
        snps:           genotype values in additive encoding (012) with samples as rows and markers as columns
        sample_ids:     sample identifier in the same order as the rows of the genotype matrix
        chr_index:      chromosome identifier in the same order as the columns of the genotype matrix
        position_index: position number (integer) in the same order as the columns of the genotype matrix
    CSV files need to have the following structure:
        sample identifier in the first column
        SNP identifier in the first row need to be formatted as CHROMOSOME_POSITION
        genotype values in additive encoding in the remaining rows and columns
    For PLINK binary files GENOTYPE.bed, GENOTYPE.bim, GENOTYPE.fam need to be in same folder
    For PLINK files GENOTYPE.ped and GENOTYPE.map files need to be in same folder
    :param genotype_file: path to genotype file
    :param sample_index: either a list/np.array containing the indices of the samples to load from the genotype matrix,
     or a single integer to load the genotype matrix from row 0 to row sample_index
    :param snp_lower_index: lower bound of batch
    :param snp_upper_index: upper bound of batch
    :return: X
    """
    suffix = genotype_file.suffix
    # only load needed samples and markers if H5, HDF5, H5PY files
    if suffix in ('.h5', '.hdf5', '.h5py'):
        with h5py.File(genotype_file, "r") as gt:
            if isinstance(sample_index, (np.ndarray, list)):
                indices, inverse = np.unique(sample_index, return_inverse=True)
                X = gt['snps'][indices, snp_lower_index:snp_upper_index]
                X = torch.tensor(X[inverse, :], dtype=torch.float64)
            else:
                X = torch.tensor(gt['snps'][:sample_index, snp_lower_index:snp_upper_index], dtype=torch.float64)
    else:
        # load binary PLINK
        if suffix in ('.bed', '.bim', '.fam'):
            prefix = genotype_file.with_suffix('').as_posix()
            gt = read_plink1_bin(prefix + '.bed', prefix + '.bim', prefix + '.fam', ref="a0", verbose=False)
            X = torch.tensor(gt.values, dtype=torch.float64)
        # load CSV
        elif suffix == '.csv':
            gt = pd.read_csv(genotype_file, index_col=0)
            X = torch.tensor(gt.values, dtype=torch.float64)
        # load PLINK
        elif suffix in ('.map', '.ped'):
            prefix = genotype_file.with_suffix('').as_posix()
            iupac_map = {"AA": "A", "GG": "G", "TT": "T", "CC": "C", "AG": "R", "GA": "R", "CT": "Y", "TC": "Y",
                         "GC": "S", "CG": "S", "AT": "W", "TA": "W", "GT": "K", "TG": "K", "AC": "M", "CA": "M"}
            with open(prefix + '.ped', 'r') as f:
                raw = []
                for line in f:
                    tmp = line.strip().split(" ")
                    snps = []
                    j = 6
                    while j < len(tmp) - 1:
                        snps.append(iupac_map[tmp[j] + tmp[j + 1]])
                        j += 2
                    raw.append(snps)
            raw = np.array(raw)
            X = get_additive_encoding(raw)
        else:
            raise NotImplementedError('Only accept H5, HDF5, H5PY, CSV, PLINK and binary PLINK genotype files')
        if sample_index is not None:
            X = X[sample_index, :]
    unique = np.unique(X)
    if not_add:
        print('Genotype might not be additive. Will not check encoding of genotype.')
    else:
        if not all(z in [0, 1, 2] for z in unique):
            raise Exception('Genotype not in additive encoding.')
    return X


def get_additive_encoding(X: np.array):
    """
    Function to compute additive encoding of genotype matrix with
        0: homozygous major allele
        1: heterozygous
        2: homozygous minor allele
    :param X: genotype matrix containing raw nucleotides in iupac single nucleotide notation
    :return: torch tensor with genotype matrix in additive encoding
    """
    alleles = []
    index_arr = []
    pairs = [['A', 'C'], ['A', 'G'], ['A', 'T'], ['C', 'G'], ['C', 'T'], ['G', 'T']]
    heterozygous_nuc = ['M', 'R', 'W', 'S', 'Y', 'K']
    for i, col in enumerate(np.transpose(X)):
        unique, inv, counts = np.unique(col, return_counts=True, return_inverse=True)
        unique = unique.astype(str)
        boolean = (unique == 'A') | (unique == 'T') | (unique == 'C') | (unique == 'G')
        tmp = np.zeros(3)
        if len(unique) > 3:
            raise Exception('More than two alleles encountered at snp ' + str(i))
        elif len(unique) == 3:
            hetero = unique[~boolean][0]
            homozygous = unique[boolean]
            for j, pair in enumerate(pairs):
                if all(h in pair for h in homozygous) and hetero != heterozygous_nuc[j]:
                    raise Exception('More than two alleles encountered at snp ' + str(i))
            tmp[~boolean] = 1.0
            tmp[np.argmin(counts[boolean])] = 2.0
        elif len(unique) == 2:
            if list(unique) in pairs:
                tmp[np.argmin(counts)] = 2.0
            else:
                tmp[(~boolean).nonzero()] = 1.0
        else:
            if unique[0] in heterozygous_nuc:
                tmp[0] = 1.0
        alleles.append(tmp)
        index_arr.append(inv)
    alleles = np.transpose(np.array(alleles))
    index_arr = np.transpose(np.array(index_arr))
    cols = np.arange(alleles.shape[1])
    return torch.tensor(alleles[index_arr, cols], dtype=torch.float64)


def load_phenotype(arguments: argparse.Namespace):
    """
    Function to load phenotype matrix
    Accept CSV, PHENO, TXT files. For .txt files assume that separator is a single space. First column should contain
    sample_ids. Name of phenotype needs to be a column name
    Raises Exception if not the right file type
    :param arguments: user input
    :return: pandas DataFrame with phenotype values and sample ids as index
    """
    suffix = arguments.y.suffix
    # load CSV
    if suffix == ".csv":
        y = pd.read_csv(arguments.y)
    # load PHENO or TXT
    elif suffix == ".txt":
        y = pd.read_csv(arguments.y, sep=" ")
    elif suffix == ".pheno":
        y = pd.read_csv(arguments.y, sep=" ")
        if {'FID','IID'}.issubset(set(y.columns)):
            y.drop(columns='FID', inplace=True)
    else:
        raise NotImplementedError('Only accept CSV, PHENO and TXT phenotype files')
    y = y.sort_values(y.columns[0]).groupby(y.columns[0]).mean()
    if arguments.y_name not in y.columns:
        raise Exception('Phenotype ' + arguments.y_name + ' is not in phenotype file ' + str(arguments.y))
    else:
        y = y[[arguments.y_name]].dropna()
    return y


def load_covariates(arguments: argparse.Namespace):
    """
    Only take CSV files: sample ids have to be in first column, will load all available columns
    Raises Exception if not the right file type
    :param arguments: user input
    :return: pandas DataFrame with sample ids as index
    """
    if arguments.cov_file.suffix == ".csv":
        covs = pd.read_csv(arguments.cov_file)
        covs = covs.sort_values(covs.columns[0]).groupby(covs.columns[0]).mean().dropna()
    else:
        raise NotImplementedError('Only accept CSV covariates files')
    return covs


def load_kinship(arguments: argparse.Namespace):
    """
    load kinship matrix from file. Only take CSV, H5, HDF5, H5PY files.
    For CSV files sample ids have to be in first column, H5, HDF5, H5PY files need to contain the kinship matrix
    with key 'kinship' and the corresponding sample ids with key 'sample_ids'.
    Raises Exception if not the right file type
    :param arguments: user input
    :return: kinship matrix and sample ids
    """
    # load CSV
    if arguments.k.suffix == ".csv":
        kin = pd.read_csv(arguments.k, index_col=0)
        K = torch.tensor(kin.values)
        sample_ids = np.array(kin.index, dtype=str)
    # load H5, HDF5, H5PY
    elif arguments.k.suffix in (".h5", ".hdf5", ".h5py"):
        with h5py.File(arguments.k, "r") as f:
            K = torch.tensor(f['kinship'][:], dtype=torch.float64)
            sample_ids = f['sample_ids'][:].astype(str)
    else:
        raise NotImplementedError('Only accept .csv, .h5, .hdf5, .h5py kinship files')
    return K, sample_ids
