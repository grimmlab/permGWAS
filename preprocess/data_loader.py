import numpy as np
import torch
import pandas as pd
import h5py
import pathlib
from pandas_plink import read_plink1_bin


class Genotype:
    """
    Class for loading of genotype data.

    **Attributes**

        - genotype_file (*pathlib.Path*): full path to genotype file for data loading
        - X (*torch.tensor*): matrix containing genotype values
        - sample_ids (*numpy.array*): ids of genotype samples
        - chromosomes (*numpy.array*): chromosome identifier of SNPs
        - positions (*numpy.array*): position identifier of SNPs
        - maf (*torch.tensor*): vector containing minor allele frequencies
        - sample_index (*numpy.array*): indices of the samples to load from the genotype matrix
        - n_samples (*int*): number of samples
        - n_snps (*int*): number of SNPs
        - maf_threshold (*int*): threshold for minor allele frequency filtering

    **Functions**

        -   load_genotype_ids(load_genotype): load sample_ids from .h5/.hdf5/.h5py file
        -   load_genotype_data(): load and encode genotype data from file, calls the following functions:
            -   load_genotype_hdf5_file(sample_index, snp_lower_index, snp_upper_index): load genotype data from
                .h5/.hdf5/.h5py files
            -   load_genotype_csv_file(): load genotype data from .csv files
            -   load_genotype_binary_plink_file(): load genotype data from binary PLINK files
            -   load_genotype_plink_file(): load genotype data from PLINK files
            -   encode_genotype(): check encoding of genotype, change to additive if necessary, create torch.tensor,
                calls the following functions:
                -   check_encoding()
                -   get_additive_encoding()
        -   load_genotype_batch_wise(maf_threshold, snp_lower_index, snp_upper_index): batch-wise loading and filtering
            of genotype data
        -   filter_monomorphic_snps(): remove monomorphic SNPs
        -   get_minor_allele_freq(): compute minor allele frequencies
        -   use_maf_filter(maf_threshold): filter for minor allele frequency
        -   save_genotype_hdf5(filename): save genotype data as .h5 file
        -   reset_genotype(): delete X for batch-wise loading
        -   get_matched_data(data, row_index): filter samples of data

    :param genotype_file: full path to genotype file
    :param maf_threshold: threshold for minor allele frequency filtering
    :param not_add: use if genotype has different / not additive encoding
    """

    def __init__(self, genotype_file: pathlib.Path, maf_threshold: int = 0, not_add: bool = False):
        self.genotype_file = genotype_file
        self.maf_threshold = maf_threshold
        self.not_add = not_add
        self.sample_ids = None
        self.chromosomes = None
        self.positions = None
        self.X = None
        self.maf = None
        self.sample_index = None
        self.n_samples = None
        self.n_snps = None

    def load_genotype_ids(self, load_genotype: bool = False) -> np.array:
        """
        Load sample_ids from .h5/.hdf5/.h5py genotype file.
        """
        if self.genotype_file.suffix not in ('.h5', '.hdf5', '.h5py'):
            raise Exception('Can only load genotype IDs from .h5/.hdf5/.h5py files.')
        with h5py.File(self.genotype_file, "r") as gt:
            self.sample_ids = gt['sample_ids'][:].astype(str)
            if not load_genotype:
                self.n_snps = len(gt['position_index'][:])

    def load_genotype_data(self):
        """
        Load and encode genotype data. Accepts PLINK files, binary PLINK files, .csv and .h5, .hdf5, .h5py files.
        For .h5/.hdf5/.h5py files only load needed samples defined in self.sample_index.
        After loading check encoding of genotype and change to additive if necessary.
        Return genotype matrix as torch.tensor, chromosomes, positions and sample_ids as np.arrays.
        """
        suffix = self.genotype_file.suffix
        if suffix in ('.h5', '.hdf5', '.h5py'):
            self.X, self.chromosomes, self.positions = self.load_genotype_hdf5_file()
        elif suffix == '.csv':
            self.X, self.sample_ids, self.chromosomes, self.positions = self.load_genotype_csv_file()
        elif suffix in ('.bed', '.bim', '.fam'):
            self.X, self.sample_ids, self.chromosomes, self.positions = self.load_genotype_binary_plink_file()
        elif suffix in ('.map', '.ped'):
            self.X, self.sample_ids, self.chromosomes, self.positions = self.load_genotype_plink_file()
        # check if genotype is in additive encoding, change encoding if not
        # change X from np.array to torch.tensor
        self.encode_genotype()
        self.n_samples = len(self.sample_ids)
        self.n_snps = len(self.positions)

    def load_genotype_batch_wise(self, device: torch.device = torch.device("cpu"), save_meta: bool = True,
                                 snp_lower_index: int = None, snp_upper_index: int = None):
        """
        Load and encode genotype data batch-wise. After loading filter for monomorphic snps and minor allele frequency.
        Only accept .h5(.hdf5/.h5py files.

        :param device: device (cpu/gpu) for computations
        :param save_meta: save chromosome and position identifiers if True
        :param snp_lower_index: lower bound of batch
        :param snp_upper_index: upper bound of batch
        """
        self.X, chromosomes, positions = self.load_genotype_hdf5_file(snp_lower_index=snp_lower_index,
                                                                      snp_upper_index=snp_upper_index)
        self.encode_genotype()
        chromosomes, positions = self.filter_monomorphic_snps(chromosomes=chromosomes, positions=positions)
        maf = self.get_minor_allele_freq()
        if self.maf_threshold != 0:
            maf, chromosomes, positions = self.use_maf_filter(maf=maf, chromosomes=chromosomes, positions=positions)
        self.X = self.X.to(device)

        if save_meta:
            if self.chromosomes is None:
                self.chromosomes = chromosomes
                self.positions = positions
                self.maf = maf
            else:
                self.chromosomes.append(chromosomes)
                self.positions.append(positions)
                self.maf.append(maf)

    def load_genotype_hdf5_file(self, snp_lower_index: int = None, snp_upper_index: int = None) -> tuple:
        """
        Load genotype matrix from .h5/.hdf5/.h5py file.
        Only load needed samples and SNPs batch wise:
            will only load specified samples given in sample_index
            if snp_upper_bound/snp_lower_bound is given, will load SNPs batch-wise, else will load all SNPs
        H5, HDF5, H5PY files need to have the following structure:
            snps:           genotype matrix either in additive encoding or in raw nucleotide encoding (biallelic
                            notation (i.e. 'AA', 'AT', ...) or iupac notation (i.e. 'A', 'W', ...)) with samples as
                            rows and markers as columns
            sample_ids:     sample identifier in the same order as the rows of the genotype matrix
            chr_index:      chromosome identifier in the same order as the columns of the genotype matrix
            position_index: position number (integer) in the same order as the columns of the genotype matrix

        :param snp_lower_index: lower bound of batch
        :param snp_upper_index: upper bound of batch

        :return: Genotype values, chromosomes and positions and sample_ids if no sample_index is specified
        """
        with h5py.File(self.genotype_file, "r") as gt:
            chromosomes = gt['chr_index'][snp_lower_index:snp_upper_index].astype(str)
            positions = gt['position_index'][snp_lower_index:snp_upper_index].astype(int)
            if isinstance(self.sample_index, (np.ndarray, list)):
                # using sample indices directly does not work for h5py --> use workaround
                indices, inverse = np.unique(self.sample_index, return_inverse=True)
                X = gt['snps'][indices, snp_lower_index:snp_upper_index]
                X = X[inverse, :]
                return X, chromosomes, positions
            else:
                raise Exception('sample_index needs to be a list in order to load certain genotype samples only.')

    def load_genotype_csv_file(self) -> (np.array, np.array, np.array, np.array):
        """
        Load .csv genotype file. File must have the following structure:
        First column must contain the sample ids, the column names should be the SNP ids as CHROMOSOME_POSITION.
        The values should be the genotype matrix either in additive encoding or in raw nucleotide encoding (biallelic
        notation (i.e. 'AA', 'AT', ...) or iupac notation (i.e. 'A', 'W', ...)).

        :return: Genotype values, sample_ids, chromosomes and positions
        """
        gt = pd.read_csv(self.genotype_file, index_col=0)
        snp_ids = np.array(list(map(lambda a: a.split("_"), gt.columns.values)))
        chromosomes = snp_ids[:, 0]
        positions = snp_ids[:, 1].astype(int)
        sample_ids = np.asarray(gt.index, dtype=str)
        X = np.asarray(gt.values)
        return X, sample_ids, chromosomes, positions

    def load_genotype_binary_plink_file(self) -> (np.array, np.array, np.array, np.array):
        """
        Load binary PLINK file, .bim, .fam, .bed files with same prefix need to be in same folder.

        :return: Genotype values, sample_ids, chromosomes and positions
        """
        prefix = self.genotype_file.with_suffix('').as_posix()
        gt = read_plink1_bin(prefix + '.bed', prefix + '.bim', prefix + '.fam', ref="a0", verbose=False)
        sample_ids = np.array(gt['fid'], dtype=str).flatten()
        positions = np.array(gt['pos']).flatten()
        chromosomes = np.array(gt['chrom']).flatten()
        X = np.asarray(gt.values)
        return X, sample_ids, chromosomes, positions

    def load_genotype_plink_file(self) -> (np.array, np.array, np.array, np.array):
        """
        Load PLINK files, .map and .ped file with same prefix need to be in same folder.
        Accepts GENOTYPENAME.ped and GENOTYPENAME.map as input

        :return: Genotype values, sample_ids, chromosomes and positions
        """
        prefix = self.genotype_file.with_suffix('').as_posix()
        with open(prefix + '.map', 'r') as f:
            chromosomes = []
            positions = []
            for line in f:
                tmp = line.strip().split(" ")
                chromosomes.append(int(tmp[0].strip()))
                positions.append(int(tmp[-1].strip()))
        chromosomes = np.array(chromosomes)
        positions = np.array(positions)
        iupac_map = {"AA": "A", "GG": "G", "TT": "T", "CC": "C", "AG": "R", "GA": "R", "CT": "Y", "TC": "Y", "GC": "S",
                     "CG": "S", "AT": "W", "TA": "W", "GT": "K", "TG": "K", "AC": "M", "CA": "M"}
        with open(prefix + '.ped', 'r') as f:
            sample_ids = []
            X = []
            for line in f:
                tmp = line.strip().split(" ")
                sample_ids.append(int(tmp[1].strip()))
                snps = []
                j = 6
                while j < len(tmp) - 1:
                    snps.append(iupac_map[tmp[j] + tmp[j + 1]])
                    j += 2
                X.append(snps)
        sample_ids = np.array(sample_ids, dtype=str)
        X = np.array(X)
        return X, sample_ids, chromosomes, positions

    def encode_genotype(self):
        """
        first check encoding of genotype, then change to additive if necessary, finally change X from np.array
        to torch.tensor
        """
        if self.not_add:
            print('Genotype might not be in additive encoding. Will not check encoding of genotype.')
            self.X = torch.tensor(self.X, dtype=torch.float64)
        else:
            enc_of_X = self.check_encoding()
            # if genotype in biallelic notation, will change to iupac notation and then encode additively
            if enc_of_X == 'biallelic':
                iupac_map = {"AA": "A", "GG": "G", "TT": "T", "CC": "C", "AG": "R", "GA": "R", "CT": "Y", "TC": "Y",
                             "GC": "S", "CG": "S", "AT": "W", "TA": "W", "GT": "K", "TG": "K", "AC": "M", "CA": "M"}
                self.X = np.vectorize(iupac_map.__getitem__)(self.X.astype(str))
                enc_of_X = 'iupac'
            if enc_of_X == 'iupac':
                self.X = torch.tensor(self.get_additive_encoding(), dtype=torch.float64)
            elif enc_of_X == 'additive':
                self.X = torch.tensor(self.X, dtype=torch.float64)
            else:
                raise Exception('Genotype in wrong encoding. Can only deal with additive, iupac and biallelic '
                                'encoding. If you want to use different encoding use flag -not_add.')

    def check_encoding(self):
        """
        Check the encoding of the genotype matrix

        :return: encoding of the genotype matrix
        """
        if self.X[0, 0].astype(str) in ['A', 'C', 'G', 'T', 'M', 'R', 'W', 'S', 'Y', 'K']:
            return 'iupac'
        elif self.X[0, 0] in [0, 1, 2]:
            return 'additive'
        elif self.X[0, 0] in ["AA", "GG", "TT", "CC", "AG", "GA", "CT", "TC", "GC", "CG", "AT", "TA", "GT", "TG",
                              "AC", "CA"]:
            return 'biallelic'
        else:
            raise Exception('Genotype in wrong encoding. Can only deal with additive, iupac and biallelic encoding. '
                            'Please check again.')

    def get_additive_encoding(self):
        """
        Function to compute additive encoding of genotype matrix with
            0: homozygous major allele
            1: heterozygous
            2: homozygous minor allele

        :return: gnotype in additive encoding
        """
        alleles = []
        index_arr = []
        pairs = [['A', 'C'], ['A', 'G'], ['A', 'T'], ['C', 'G'], ['C', 'T'], ['G', 'T']]
        heterozygous_nuc = ['M', 'R', 'W', 'S', 'Y', 'K']
        for i, col in enumerate(np.transpose(self.X)):
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
        return alleles[index_arr, cols]

    def filter_monomorphic_snps(self, chromosomes: np.array = None, positions: np.array = None) -> (np.array, np.array):
        """
        Remove monomorphic SNPs, i.e., SNPs that are constant

        :param chromosomes: vector with chromosome identifiers
        :param positions: vector with position identifiers

        :return filtered chromosomes and positions
        """
        tmp = self.X == self.X[0, :]
        self.X = self.X[:, ~tmp.all(0)]
        if chromosomes is None:
            self.chromosomes = self.chromosomes[~tmp.all(0)]
            self.positions = self.positions[~tmp.all(0)]
        else:
            return chromosomes[~tmp.all(0)], positions[~tmp.all(0)]

    def get_minor_allele_freq(self):
        """
        Function to calculate minor allele frequencies of each SNP

        :return: vector containing frequencies
        """

        return (torch.sum(self.X, 0)) / (2 * self.X.shape[0])

    def use_maf_filter(self, maf: torch.tensor = None, chromosomes: np.array = None, positions: np.array = None) \
            -> (torch.tensor, np.array, np.array):
        """
        filter genotype by minor allele frequency

        :param maf: vector containing minor allele frequencies
        :param chromosomes: vector with chromosome identifiers
        :param positions: vector with position identifiers

        :return: tensor with filtered maf frequencies, chromosomes and positions
        """
        if maf is None:
            tmp = self.maf > (self.maf_threshold / 100)
            self.X = self.X[:, tmp]
            self.chromosomes = self.chromosomes[tmp]
            self.positions = self.positions[tmp]
            self.maf = self.maf[tmp]
        else:
            # for batch-wise loading
            tmp = maf > (self.maf_threshold / 100)
            self.X = self.X[:, tmp]
            return maf[tmp], chromosomes[tmp], positions[tmp]

    def save_genotype_hdf5(self, filename: pathlib.Path):
        """
        Save genotype data to .h5 file

        :param filename: Full path to new genotype file
        """
        if any(elem is None for elem in [self.X, self.sample_ids, self.chromosomes, self.positions]):
            raise Exception('Cannot save genotype file. Some values are None, please check again.')
        print('Save genotype data as ' + filename.as_posix() + '.\nThis might take some time.')
        with h5py.File(filename.with_suffix('.h5'), 'w') as f:
            f.create_dataset('sample_ids', data=self.sample_ids.astype(bytes), chunks=True, compression="gzip")
            f.create_dataset('chr_index', data=self.chromosomes.astype(bytes), chunks=True, compression="gzip")
            f.create_dataset('position_index', data=self.positions.astype(int), chunks=True, compression="gzip")
            f.create_dataset('snps', data=self.X, chunks=True, compression="gzip", compression_opts=7)
        print('Done saving H5 file.')

    def reset_genotype(self):
        """
        Delete X for batchwise loading
        """
        self.X = None

    @staticmethod
    def get_matched_data(data, row_index: np.array):
        """
        Get rows of data specified in index array

        :param data: data to match, either np.array or torch.tensor
        :param row_index: row-index array for filtering / matching
        """
        if data.ndim == 2:
            return data[row_index, :]
        if data.ndim == 1:
            return data[row_index]
        else:
            raise Exception('Cannot match data, dimensions are wrong. Expected dimension 1 or 2 but got '
                            + str(data.ndim) + ' instead. Please check again.')


class Dataset(Genotype):
    """
    Class for loading and preparation of genotype, phenotype, kinship and covariates.

    **Attributes**

        - genotype_file (*pathlib.Path*): full path to genotype file for data loading
        - X (*torch.tensor*): matrix containing genotype values
        - sample_ids (*numpy.array*): ids of genotype samples
        - chromosomes (*numpy.array*): chromosome identifier of SNPs
        - positions (*numpy.array*): position identifier of SNPs
        - y (*torch.tensor*): tensor containing phenotypic values
        - K (*torch.tensor*): kinship matrix
        - fixed (*torch.tensor*): matrix containing fixed effects, i.e. vector of ones and covariates if available
        - maf (*torch.tensor*): vector containing minor allele frequencies
        - sample_index (*np.array*): vector containing sample indices for batch-wise loading of X
        - n_samples (*int*): number of samples
        - n_snps (*int*): number of SNPs
        - maf_threshold (*int*): threshold for minor allele frequency filtering

    **Functions**

        -   load_and_prepare_data(): load load and match data, calls the following functions:
            -   see class Genotype for all genotype specific functions
            -   load_phenotype(phenotype_file, trait): load phenotype fom file
            -   load_kinship(kinship_file): load kinship matrix from file
            -   compute_rrk_kinship(): compute realized relationship kernel
            -   normalize_kinship(): normalize kinship matrix using a Gower's centered matrix
            -   load_covariates(covariates_file, column_list): load covariates from file
            -   get_fixed_effects(): create fixed effects vector/matrix
            -   match_data(data_ids1, data_ids2): match ids of two datasets
        -   to_device(device): move tensors to device

    :param genotype_file: full path to genotype file
    :param phenotype_file: full path to phenotype file
    :param trait: name of phenotypic trait to use
    :param maf_threshold: minor allele frequency threshold to use for SNP filtering, default is 0 (no filtering)
    :param load_genotype: bool, if False load genotype batch-wise during computations, default is False
    :param kinship_file: full path to kinship file, optional, if missing, compute rrk kinship
    :param covariate_file: full path to covariate file, optional
    :param covariate_list: list of covariates to use, optional
    :param not_add: use if genotype has different / not additive encoding
    """

    def __init__(self, genotype_file: pathlib.Path, phenotype_file: pathlib.Path, trait: str, maf_threshold: int = 0,
                 load_genotype: bool = False, kinship_file: pathlib.Path = None, covariate_file: pathlib.Path = None,
                 covariate_list: list = None, not_add: bool = False):
        super().__init__(genotype_file=genotype_file, maf_threshold=maf_threshold, not_add=not_add)

        self.y = None
        self.K = None
        self.fixed = None
        self.load_and_prepare_data(phenotype_file=phenotype_file, trait=trait, load_genotype=load_genotype,
                                   kinship_file=kinship_file, covariate_file=covariate_file,
                                   covariate_list=covariate_list)

    def load_and_prepare_data(self, phenotype_file: pathlib.Path, trait: str, load_genotype: bool = False,
                              kinship_file: pathlib.Path = None, covariate_file: pathlib.Path = None,
                              covariate_list: list = None):
        """
        Load and match genotype, phenotype, kinship and covariates.
        1. Load phenotype from file.
        2. Load genotype and match with pheno:
            If load_genotype is False, only load geno sample_ids from file and match data
            Load genotype sample_ids, match with pheno and load geno data only for needed samples
        3. Filter genotype for monomorphic SNPs and minor allele frequency
        4. Load kinship from file and match with geno, or compute kinship from geno data
        5. if available load covariates from file

        :param phenotype_file: full path to phenotype file
        :param trait: name of phenotypic trait to use
        :param load_genotype: bool, if False load genotype batch-wise during computations, default is False
        :param kinship_file: full path to kinship file, optional, if missing, compute rrk kinship
        :param covariate_file: full path to covariate file, optional
        :param covariate_list: list of covariates to use, optional
        """
        # load phenotype
        y, y_ids = self.load_phenotype(phenotype_file=phenotype_file, trait=trait)
        # load genotype
        if not load_genotype:
            # only load and match sample ids of genotype, values will be loaded batch-wise during computations
            self.load_genotype_ids(load_genotype=False)
            pheno_index, self.sample_index = self.match_data(data_ids1=y_ids, data_ids2=self.sample_ids)
            if len(pheno_index) == 0:
                raise Exception("Samples of genotype and phenotype do not match.")
        else:
            if self.genotype_file.suffix in ('.h5', '.hdf5', '.h5py'):
                # load genotype sample ids, match data and only load genotype values for needed samples
                self.load_genotype_ids()
                pheno_index, self.sample_index = self.match_data(data_ids1=y_ids, data_ids2=self.sample_ids)
                if len(pheno_index) == 0:
                    raise Exception("Samples of genotype and phenotype do not match.")
                self.load_genotype_data()
            else:
                self.load_genotype_data()
                pheno_index, self.sample_index = self.match_data(data_ids1=y_ids, data_ids2=self.sample_ids)
                if len(pheno_index) == 0:
                    raise Exception("Samples of genotype and phenotype do not match.")
                self.X = self.get_matched_data(data=self.X, row_index=self.sample_index)
            self.filter_monomorphic_snps()
            self.maf = self.get_minor_allele_freq()
            if self.maf_threshold != 0:
                self.use_maf_filter()
            self.n_snps = len(self.positions)
        self.y = self.get_matched_data(data=y, row_index=pheno_index)
        self.sample_ids = self.get_matched_data(data=self.sample_ids, row_index=self.sample_index)
        self.n_samples = len(self.y)
        # kinship
        if kinship_file is None:
            # compute kinship matrix
            self.K = self.compute_rrk_kinship()
        else:
            # load kinship from file
            self.K, K_ids = self.load_kinship(kinship_file=kinship_file)
            _, K_index = self.match_data(data_ids1=self.sample_ids, data_ids2=K_ids)
            if len(K_index) == len(self.sample_ids):
                self.K = self.K[K_index, :][:, K_index]
            else:
                raise Exception("Sample ids of genotype and kinship matrix do not match. Please check again")
        self.normalize_kinship()
        # fixed effects
        if covariate_file is not None:
            # load covariates from file
            cov = self.load_covariates(covariate_file=covariate_file, covariate_list=covariate_list)
            cov_ids = np.asarray(cov.index, dtype=y_ids.dtype).flatten()
            _, cov_index = self.match_data(data_ids1=self.sample_ids, data_ids2=cov_ids)
            if len(cov_index) == len(self.sample_ids):
                self.fixed = torch.tensor(cov.values, dtype=torch.float64).flatten()[cov_index]
            else:
                raise Exception('Sample ids of covariates and phenotype do not match.')
        self.get_fixed_effects()

    def load_phenotype(self, phenotype_file: pathlib.Path, trait: str) -> (torch.Tensor, np.array):
        """
        Load phenotype from file. Accept .csv and single white space separated .txt and .pheno files.
        Phenotype data needs to contain sample identifiers as first column and phenotypic traits as remaining columns.
        The trait name should be the respective column name. Can contain more than one phenotype columns.
        Will drop NAN values during preparation and compute mean over replicates.

        :param phenotype_file: full path to phenotype file
        :param trait: name of phenotypic trait / column to use

        :return: tensor containing phenotypic traits and array containing respective sample_ids
        """

        suffix = phenotype_file.suffix
        # load CSV
        if suffix == ".csv":
            y = pd.read_csv(phenotype_file)
        # load PHENO or TXT
        elif suffix == ".txt":
            y = pd.read_csv(phenotype_file, sep=" ")
        elif suffix == ".pheno":
            y = pd.read_csv(phenotype_file, sep=" ")
            if {'FID', 'IID'}.issubset(set(y.columns)):
                y.drop(columns='FID', inplace=True)
        else:
            raise NotImplementedError('Only accept CSV, PHENO and TXT phenotype files')
        # account for replicates
        y = y.sort_values(y.columns[0]).groupby(y.columns[0]).mean()
        if trait not in y.columns:
            raise Exception('Phenotype ' + trait + ' is not in phenotype file ' + phenotype_file.as_posix())
        else:
            y = y[[trait]].dropna()
        return torch.tensor(y.values, dtype=torch.float64).flatten(), np.asarray(y.index, dtype=str).flatten()

    def load_kinship(self, kinship_file: pathlib.Path) -> (torch.tensor, np.array):
        """
        load kinship matrix from file. Only take .csv or .h5/.hdf5/.h5py files.
        For .csv files sample ids have to be in first column, .h5/.hdf5/.h5py files need to contain the kinship matrix
        with key 'kinship' and the corresponding sample ids with key 'sample_ids'.

        :param kinship_file: full path to kinship file

        :return: torch.tensor containing kinship matrix and array with sample ids
        """
        # load .csv
        suffix = kinship_file.suffix
        if suffix == ".csv":
            kin = pd.read_csv(kinship_file, index_col=0)
            K = torch.tensor(kin.values)
            sample_ids = np.array(kin.index, dtype=str)
        # load .h5/.hdf5/.h5py
        elif suffix in (".h5", ".hdf5", ".h5py"):
            with h5py.File(kinship_file, "r") as f:
                K = torch.tensor(f['kinship'][:], dtype=torch.float64)
                sample_ids = f['sample_ids'][:].astype(str)
        else:
            raise NotImplementedError('Only accept .csv, .h5, .hdf5, .h5py kinship files')
        return K, sample_ids

    def compute_rrk_kinship(self) -> torch.tensor:
        """
        compute realized relationship kernel as kinship matrix

        :return: kinship matrix
        """
        if self.X is None:
            raise Exception('Cannot compute kinship matrix, no genotype matrix available.')
        X_stand = (self.X - self.X.mean(axis=0)) / self.X.std(axis=0)
        K = torch.matmul(X_stand, torch.t(X_stand)) / self.X.shape[1]
        # set negative values in K to zero
        return torch.where(K > 0, K, 0.)

    def normalize_kinship(self):
        """
        normalize kinship matrix using a Gower's centered matrix
        """
        n = self.K.shape[0]
        P = (torch.eye(n, dtype=self.K.dtype, device=self.K.device) -
             torch.ones(n, n, dtype=self.K.dtype, device=self.K.device) / n)
        self.K = (n - 1) / torch.sum(torch.mul(P, self.K)) * self.K

    def load_covariates(self, covariate_file: pathlib.Path, covariate_list: list = None) -> torch.tensor:
        """
        Only take .csv files: sample ids have to be in first column, if column_list is available, will load all columns
        specified, else will load all available columns

        :param covariate_file: full path to covariates file
        :param covariate_list: list containing column names/headers of covariates to load

        :return: pandas DataFrame containing covariates with sample ids as index
        """
        if covariate_file.suffix == ".csv":
            covs = pd.read_csv(covariate_file)
            covs = covs.sort_values(covs.columns[0]).groupby(covs.columns[0]).mean().dropna()
            if covariate_list is not None:
                if set(covariate_list).issubset(set(covs.columns)):
                    covs = covs[covariate_list]
                else:
                    raise Exception('Specified covariates are not available in covariate file. Please check again.')
        else:
            raise NotImplementedError('Only accept .csv covariates files')
        return covs

    def get_fixed_effects(self):
        """
        Check for covariates and create fixed effects matrix with ones as first column and covariates as remaining
        columns if available --> dim: (n, c+1)
        """
        if self.fixed is None:
            self.fixed = torch.ones((len(self.y), 1), dtype=torch.float64)
        elif self.fixed.ndim == 1:
            self.fixed = torch.stack((torch.ones(len(self.y), dtype=torch.float64), self.fixed), dim=1)
        else:
            self.fixed = torch.cat((torch.ones((len(self.y), 1), dtype=torch.float64), self.fixed), dim=1)

    def to_device(self, device: torch.device):
        """
        move data to device

        :param device: cpu or cuda
        """
        self.y = self.y.to(device)
        self.K = self.K.to(device)
        self.fixed = self.fixed.to(device)

    @staticmethod
    def match_data(data_ids1: np.array, data_ids2: np.array) -> (np.array, np.array):
        """
        match two datasets

        :param data_ids1: ids of first dataset
        :param data_ids2: ids of second dataset

        :return: two arrays with indices of matched data
        """
        return (np.reshape(data_ids1, (data_ids1.shape[0], 1)) == data_ids2.astype(data_ids1.dtype)).nonzero()
