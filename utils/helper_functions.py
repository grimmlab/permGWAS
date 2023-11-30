import yaml
import pathlib
import importlib
import inspect
import numpy as np

import models


def parse_config_file(args: dict) -> dict:
    """
    Read yaml config file to update all user specified arguments

    :param args: dict with user specified arguments

    :return: updated dict with arguments
    """
    config_path = pathlib.Path(args["config_file"])
    if not config_path.is_file():
        raise FileNotFoundError('Specified config file does not exist. Please check again.')
    if config_path.suffix not in ['.yaml', '.yml']:
        raise Exception('Only accept yaml config files. Please check again.')
    config = yaml.safe_load(open(config_path))
    args.update(config)
    return args


def get_model_class_name(model_name: str = 'lmm'):
    """
    Get class name of model for user input

    :param model_name: user input of model name
    :return: model class name
    """
    if model_name in models.__all__:
        model_name = 'models.' + model_name
        for name, cls in inspect.getmembers(importlib.import_module(model_name), inspect.isclass):
            if cls.__module__ == model_name:
                return cls
            else:
                raise NotImplementedError('No class named ', model_name)
    else:
        raise NotImplementedError('No class named ', model_name)


def estimate_heritability(v_g: float, v_e: float) -> float:
    """
    compute narrow sense heritability
    :param v_g: genetic variance component
    :param v_e: residual variance component
    :return: narrow sense heritability
    """
    return v_g / (v_g + v_e)


def compute_perm_threshold(min_p_val: np.array, sig_level: int) -> float:
    """
    Compute permutation-based threshold
    :param min_p_val: array with minimal p-values
    :param sig_level: significance level as percentage value
    :return: threshold
    """
    return np.percentile(min_p_val, sig_level)


def compute_bonf_threshold(number_snps: int, sig_level: int) -> float:
    """
    Compute Bonferroni threshold
    :param number_snps: number of SNPs
    :param sig_level: significance level as percentage value
    :return: threshold
    """
    return (sig_level / 100) / number_snps


def print_summary_stats(genotype_file: pathlib.Path, phenotype_file: pathlib.Path, trait: str, samples: int, snps: int,
                        model: str, maf_threshold: int, perm: int, v_g: float, v_e: float, h2: float, bonf1: float,
                        bonf5: float, perm1: float, perm5: float, time: float, kinship_file: pathlib.Path = None,
                        covariate_file: pathlib.Path = None, covariate_list: list = None, perm_method: str = None):
    """
    Print summary statistics

    :param genotype_file:
    :param phenotype_file:
    :param trait: name of phenotypic trait
    :param samples: number of samples used
    :param snps: number of SNPs used
    :param model: model used for GWAS
    :param maf_threshold: threshold used for maf filtering
    :param perm: number of permutations
    :param v_g: genetic variance component
    :param v_e: residual variiance component
    :param h2: narrow-sense heritability
    :param bonf1: Bonferroni threshold significance level 1%
    :param bonf5: Bonferroni threshold significance level 5%
    :param perm1: permutation-based threshold significance level 1%
    :param perm5: permutation-based threshold significance level 5%
    :param kinship_file:
    :param covariate_file:
    :param covariate_list: list containing covariates
    :param perm_method: method used for permutations
    """
    print('\n')
    print('+++++++++ Summary Statistics +++++++++')
    print('## Genotype file: ' + genotype_file.as_posix())
    print('## Phenotype file: ' + phenotype_file.as_posix())
    print('## Phenotype: ' + trait)
    if covariate_file is not None:
        print('## Covariate file: ' + covariate_file.as_posix())
        if covariate_list is not None:
            print('## Used covariates: ' + ",".join(covariate_list))
        else:
            print('## Used all available covariates')
    if kinship_file is not None:
        print('## Kinship file: ' + kinship_file.as_posix())
    print('## Number of individuals: ' + str(samples))
    print('## Number of SNPs: ' + str(snps))
    print('## Model: ' + model)
    print('## MAF threshold: ' + str(maf_threshold))
    print('## Number of permutations: ' + str(perm))
    if perm_method is not None:
        print('## permutation method: ' + perm_method)
    if model == 'lmm':
        print('## v_g estimate in null model: ' + str(v_g))
        print('## v_e estimate in null model: ' + str(v_e))
        print('## Narrow-sense heritability estimate: ' + str(h2))
    print('## Bonferroni threshold (1% significance level): ' + str(bonf1))
    print('## Bonferroni threshold (5% significance level): ' + str(bonf5))
    if perm1 is not None:
        print('## Permutation-based threshold (1% significance level): ' + str(perm1))
        print('## Permutation-based threshold (5% significance level): ' + str(perm5))
    print('## Total time: %.2f s' %time)
    print('+++++++++++++++++++++++++++')
    print('\n')


def write_summary_stats(out_dir: pathlib.Path, out_file: str, genotype_file: pathlib.Path, phenotype_file: pathlib.Path,
                        trait: str, samples: int, snps: int, model: str, maf_threshold: int, perm: int, v_g: float,
                        v_e: float, h2: float, bonf1: float, bonf5: float, perm1: float, perm5: float, time: float,
                        kinship_file: pathlib.Path = None, covariate_file: pathlib.Path = None,
                        covariate_list: list = None, perm_method: str = None):
    """
    Save summary statistics to txt file

    :param out_dir:
    :param out_file:
    :param genotype_file:
    :param phenotype_file:
    :param trait: name of phenotypic trait
    :param samples: number of samples used
    :param snps: number of SNPs used
    :param model: model used for GWAS
    :param maf_threshold: threshold used for maf filtering
    :param perm: number of permutations
    :param v_g: genetic variance component
    :param v_e: residual variance component
    :param h2: narrow-sense heritability
    :param bonf1: Bonferroni threshold significance level 1%
    :param bonf5: Bonferroni threshold significance level 5%
    :param perm1: permutation-based threshold significance level 1%
    :param perm5: permutation-based threshold significance level 5%
    :param kinship_file:
    :param covariate_file:
    :param covariate_list: list containing covariates
    :param perm_method: method used for permutations
    """
    filename = out_dir.joinpath('summary_statistics_' + pathlib.Path(out_file).with_suffix('.txt').as_posix())
    with open(filename, 'w') as f:
        f.write('Summary Statistics:\n')
        f.write('## Genotype file:\t' + genotype_file.as_posix() + '\n')
        f.write('## Phenotype file:\t' + phenotype_file.as_posix() + '\n')
        f.write('## Phenotype:\t' + trait + '\n')
        if covariate_file is not None:
            f.write('## Covariate file:\t' + covariate_file.as_posix() + '\n')
            if covariate_list is not None:
                f.write('## Used covariates:\t' + ",".join(covariate_list) + '\n')
            else:
                f.write('## Used all available covariates' + '\n')
        if kinship_file is not None:
            f.write('## Kinship file:\t' + kinship_file.as_posix() + '\n')
        f.write('## Number of individuals:\t' + str(samples) + '\n')
        f.write('## Number of SNPs:\t' + str(snps) + '\n')
        f.write('## Model:\t' + model + '\n')
        f.write('## MAF threshold:\t' + str(maf_threshold) + '\n')
        f.write('## Number of permutations:\t' + str(perm) + '\n')
        if perm_method is not None:
            f.write('## permutation method:\t' + perm_method + '\n')
        if model == 'lmm':
            f.write('## v_g estimate in null model:\t' + str(v_g) + '\n')
            f.write('## v_e estimate in null model:\t' + str(v_e) + '\n')
            f.write('## Narrow-sense heritability estimate:\t' + str(h2) + '\n')
        f.write('## Bonferroni threshold (1% significance level):\t' + str(bonf1) + '\n')
        f.write('## Bonferroni threshold (5% significance level):\t' + str(bonf5) + '\n')
        if perm1 is not None:
            f.write('## Permutation-based threshold (1% significance level):\t' + str(perm1) + '\n')
            f.write('## Permutation-based threshold (5% significance level):\t' + str(perm5) + '\n')
        f.write('## Total time:\t' + str(time) + ' s\n')


def get_summary_stats(out_dir: pathlib.Path, out_file: str, genotype_file: pathlib.Path, phenotype_file: pathlib.Path,
                      trait: str, samples: int, snps: int, model: str, maf_threshold: int, perm: int, v_g: float,
                      v_e: float, min_p_val: np.array, time: float, kinship_file: pathlib.Path = None,
                      covariate_file: pathlib.Path = None, covariate_list: list = None, perm_method: str = None):
    """
    Compute summary statistics, print and save them to file

    :param out_dir:
    :param out_file:
    :param genotype_file:
    :param phenotype_file:
    :param trait: name of phenotypic trait
    :param samples: number of samples used
    :param snps: number of SNPs used
    :param model: model used for GWAS
    :param maf_threshold: threshold used for maf filtering
    :param perm: number of permutations
    :param v_g: genetic variance component
    :param v_e: residual variance component
    :param min_p_val: minimal p-values
    :param kinship_file:
    :param covariate_file:
    :param covariate_list: list containing covariates
    :param perm_method: method used for permutations
    """
    if model == 'lmm':
        h2 = estimate_heritability(v_g=v_g, v_e=v_e)
    else:
        h2 = None
    bonf1 = compute_bonf_threshold(number_snps=snps, sig_level=1)
    bonf5 = compute_bonf_threshold(number_snps=snps, sig_level=5)
    if min_p_val is not None:
        perm1 = compute_perm_threshold(min_p_val=min_p_val, sig_level=1)
        perm5 = compute_perm_threshold(min_p_val=min_p_val, sig_level=5)
    else:
        perm1 = None
        perm5 = None

    write_summary_stats(out_dir=out_dir, out_file=out_file, genotype_file=genotype_file, phenotype_file=phenotype_file,
                        trait=trait, samples=samples, snps=snps, model=model, maf_threshold=maf_threshold, perm=perm,
                        v_g=v_g, v_e=v_e, h2=h2, bonf1=bonf1, bonf5=bonf5, perm1=perm1, perm5=perm5, time=time,
                        kinship_file=kinship_file, covariate_file=covariate_file, covariate_list=covariate_list,
                        perm_method=perm_method)
    print_summary_stats(genotype_file=genotype_file, phenotype_file=phenotype_file,
                        trait=trait, samples=samples, snps=snps, model=model, maf_threshold=maf_threshold, perm=perm,
                        v_g=v_g, v_e=v_e, h2=h2, bonf1=bonf1, bonf5=bonf5, perm1=perm1, perm5=perm5, time=time,
                        kinship_file=kinship_file, covariate_file=covariate_file, covariate_list=covariate_list,
                        perm_method=perm_method)
