import numpy as np
import torch
import pandas as pd
import h5py
from pandas_plink import read_plink1_bin


def load_genotype(gt_file):
    """
    load genotype
    for binary PLINK files: need either "NAME.bed", "NAME.bim" or "NAME.fam", all need to be in same folder
    :param gt_file: takes PLink files, binary PLINK files, .csv and .h5, .hdf5, .h5py files
    :return: genotype in additive encoding, with sample ids and SNP positions
    """
    suffix = gt_file.suffix
    if suffix in ('.bed', '.bim', '.fam'):
        x_file = gt_file.with_suffix('').as_posix()
        gt = read_plink1_bin(x_file + '.bed', x_file + '.bim',
                             x_file + '.fam', ref="a0", verbose=False)
        sample_ids = np.array(gt['fid'], dtype=np.int).flatten()
        positions = np.array(gt['pos']).flatten()
        chromosomes = np.array(gt['chrom']).flatten()
        X = torch.tensor(gt.values, dtype=torch.float64)
    elif suffix in ('.h5', '.hdf5', '.h5py'):
        with h5py.File(gt_file, "r") as gt:
            chromosomes = gt['chr_index'][:].astype(str)
            positions = gt['position_index'][:].astype(int)
            sample_ids = gt['sample_ids'][:].astype(int)
            X = torch.tensor(gt['snps'][:], dtype=torch.float64)
    elif suffix == '.csv':
        gt = pd.read_csv(gt_file, index_col=0)
        identifiers = np.array(list(map(lambda a: a.split("_"), gt.columns.values)))
        chromosomes = identifiers[:, 0]
        positions = identifiers[:, 1]
        sample_ids = np.asarray(gt.index)
        X = torch.tensor(gt.values, dtype=torch.float64)
    elif suffix in ('map', 'ped'):
        x_file = gt_file.with_suffix('').as_posix()
        with open(x_file + '.map', 'r') as f:
            chromosomes = []
            positions = []
            for line in f:
                tmp = line.strip().split("\t")
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


def encode_homozygous(matrix):
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


def load_phenotype(pt_file, pt_name):
    """
    load phenotype
    for csv-files: column name of sample ids has to be "accession_id"
    :param pt_file: takes csv-files
    :param pt_name: name of phenotype (column) to be used
    :return: pandas DataFrame with sample ids as index
    """
    if pt_file.suffix == ".csv":
        y = pd.read_csv(pt_file).sort_values(['accession_id']).groupby('accession_id').mean()
        if pt_name not in y.columns:
            raise Exception('Phenotype ' + pt_name + ' is not in phenotype file ' + str(pt_file))
        else:
            y = y[[pt_name]].dropna()
    else:
        raise NotImplementedError('Only accept .csv phenotype files')
    return y


def load_covariates(cov_file):
    """
    for csv-files: column name of sample ids has to be "sample_id"
    :param cov_file: file containing covariates, take .csv
    :return: pandas DataFrame with sample ids as index
    """
    if cov_file.suffix == ".csv":
        covs = pd.read_csv(cov_file).sort_values(['accession_id']).groupby('accession_id').mean().dropna()
    else:
        raise NotImplementedError('Only accept .csv covariates files')
    return covs


def load_kinship(k_file):
    """
    :param k_file: file containing kinship matrix, take .csv
    :return: kinship matrix and sample ids
    """
    if k_file.suffix == ".csv":
        kin = pd.read_csv(k_file, index_col=0)
        K = torch.tensor(kin.values)
        sample_ids = np.array(kin.index)
    else:
        raise NotImplementedError('Only accept .csv kinship files')
    return K, sample_ids
