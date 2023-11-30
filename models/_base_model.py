import abc
import torch
import pathlib
import pandas as pd
import numpy as np

from preprocess import data_loader
from postprocess import plot_functions


class BaseModel(abc.ABC):

    def __init__(self, dataset: data_loader.Dataset, batch_size: int, device: torch.device, perm: int = None,
                 perm_batch_size: int = None):
        self.dataset = dataset
        self.batch_size = batch_size
        self.device = device
        self.perm_batch_size = perm_batch_size
        self.perm = perm
        self.v_g = None  # genetic variance component for LMM
        self.v_e = None  # residual variance component for LMM
        self.delta = None  # v_e/v_g
        self.effect_size = None  # effect sizes for all SNPs
        self.SE = None  # standard errors for all SNPs
        self.test_stat = None  # tests statistics for all SNPs
        self.p_value = None  # p_values for all SNPs
        self.seeds = None  # seeds for permutation with numpy generator
        self.perm_p_val = None  # permutation-based p-values
        self.min_p_value = None  # minimal p-values for all permutations

    @abc.abstractmethod
    def gwas(self):
        """
        Function to perform batch-wise computation of univariate test

        """

    @abc.abstractmethod
    def perm_gwas(self, **kwargs):
        """
        Function to perform batch-wise computation of permutation-based test

        """

    # general methods
    def perm_seeds(self) -> np.array:
        """
        get seeds for permutations

        :return: array with seeds
        """
        rng = np.random.default_rng()
        return rng.choice(1000000, self.perm, replace=False)

    def permute(self, data: torch.tensor) -> torch.tensor:
        """
        Create tensor with permutations of input data

        :param data: input data to permute of shape (n,c) or (n)

        :return: tensor with permuted data of shape (p,n,c) or (n,p)
        """
        data = data.to(torch.device("cpu"))
        x_perm = []
        for seed in self.seeds:
            tmp = np.random.default_rng(seed=seed)
            x_perm.append(tmp.permutation(data, axis=0))
        if data.ndim == 1:
            return torch.t(torch.tensor(np.array(x_perm), dtype=torch.float64, device=self.device))
        else:
            return torch.tensor(np.array(x_perm), dtype=torch.float64, device=self.device)

    def save_results(self, data_dir: pathlib.Path, filename: str):
        """
        Save p-values results to csv file as p_values_filename. If permutations were computed, also save
        minimal p-values as min_p_values_filename.

        :param data_dir: full path to results directory
        :param filename: name of results file
        """
        df = pd.DataFrame({'CHR': self.dataset.chromosomes,
                           'POS': self.dataset.positions,
                           'p_value': self.p_value,
                           'test_stat': self.test_stat,
                           'maf': self.dataset.maf,
                           'SE': self.SE,
                           'effect_size': self.effect_size})
        if self.perm_p_val is not None:
            df['adjusted_p_val'] = self.perm_p_val
        df.to_csv(data_dir.joinpath('p_values_' + filename), index=False)
        if self.min_p_value is not None:
            df_min = pd.DataFrame({'seed': self.seeds,
                                   'min_p_val': self.min_p_value})
            df_min.to_csv(data_dir.joinpath('min_p_values_' + filename), index=False)

    def manhattan_plot(self, data_dir: pathlib.Path, filename: str, sig_level: int = 5):
        """
        Save Manhattan plot as manhattan_FILENAME.png to data_dir

        :param data_dir: full path to save directory
        :param filename: name of file
        :param sig_level: significance level for Bonferroni and perm thresholds, default is 5
        """
        df = pd.DataFrame({'CHR': self.dataset.chromosomes,
                           'POS': self.dataset.positions,
                           'p_value': self.p_value})

        plot_functions.manhattan_plot(df=df, data_dir=data_dir, filename=filename, min_p_values=self.min_p_value,
                                      sig_level=sig_level)

    def qq_plot(self, data_dir: pathlib.Path, filename: str):
        """
        Save QQ-plot as qq_plot_FILENAME.png to data_dir

        :param data_dir: full path to save directory
        :param filename: name of file
        """
        plot_functions.qq_plot(p_values=self.p_value, data_dir=data_dir, filename=filename)
