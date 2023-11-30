import pathlib
import time
import torch

from preprocess import data_loader
from utils import helper_functions


def run(genotype_file: pathlib.Path, phenotype_file: pathlib.Path, model: str, trait: str = 'phenotype_value',
        kinship_file: pathlib.Path = None, covariate_file: pathlib.Path = None, covariate_list: list = None,
        maf_threshold: int = 0, load_genotype: bool = False,
        out_dir: pathlib.Path = pathlib.Path.cwd().joinpath('results'), out_file: str = None,
        device: torch.device = torch.device('cpu'), perm: int = 0, perm_method: str = 'x',
        adj_p_value: bool = False, batch_size: int = 50000, perm_batch_size: int = 1000, manhattan: bool = False,
        qqplot: bool = False, not_add: bool = False):
    # check user specified arguments
    start = time.time()
    print('Start loading data now')

    # load data
    dataset = data_loader.Dataset(genotype_file=genotype_file, phenotype_file=phenotype_file, trait=trait,
                                  maf_threshold=maf_threshold, load_genotype=load_genotype, kinship_file=kinship_file,
                                  covariate_file=covariate_file, covariate_list=covariate_list, not_add=not_add)
    dataset.to_device(device=device)
    have_data = time.time()
    print('Loaded data, elapsed time: %f s.' % (have_data - start))
    print('Start performing GWAS on phenotype %s for %d samples and %d SNPs.'
          % (trait, dataset.n_samples, dataset.n_snps))

    # perform GWAS
    gwas_model = helper_functions.get_model_class_name(model_name=model)(dataset=dataset, batch_size=batch_size,
                                                                         device=device, perm=perm,
                                                                         perm_batch_size=perm_batch_size)
    gwas_model.gwas()
    done_gwas = time.time()
    print('Done performing GWAS on phenotype %s for %d samples and %d SNPs.\n'
          'Elapsed time: %f s' % (trait, dataset.n_samples, len(dataset.positions), done_gwas - have_data))

    # perform GWAS with permutations
    if perm > 0:
        print('Start performing GWAS with %d permutations.' % perm)
        gwas_model.perm_gwas(perm_method=perm_method, adj_p_value=adj_p_value)
        done_perm = time.time()
        print('Done performing GWAS with %d permutations.\n'
              'Elapsed time: %f s' % (perm, done_perm - done_gwas))

    # save results
    print('Save results.')
    gwas_model.save_results(data_dir=out_dir, filename=out_file)
    total_time = time.time() - start
    print('Total time: ', total_time)

    # plots
    if manhattan:
        print('Save Manhattan plot with significance level of 5%.')
        gwas_model.manhattan_plot(data_dir=out_dir, filename=out_file, sig_level=5)
        total_time = time.time() - start
    if qqplot:
        print('Save QQ-plot.')
        gwas_model.qq_plot(data_dir=out_dir, filename=out_file)
        total_time = time.time() - start

    # summary statistics
    if not load_genotype:
        # reset number of SNPs in case of batch-wise loading
        dataset.n_snps = len(dataset.positions)
    helper_functions.get_summary_stats(out_dir=out_dir, out_file=out_file, genotype_file=genotype_file,
                                       phenotype_file=phenotype_file, trait=trait, samples=dataset.n_samples,
                                       snps=dataset.n_snps, model=model, maf_threshold=maf_threshold, perm=perm,
                                       v_g=gwas_model.v_g.item(), v_e=gwas_model.v_e.item(),
                                       min_p_val=gwas_model.min_p_value, time=total_time, kinship_file=kinship_file,
                                       covariate_file=covariate_file, covariate_list=covariate_list,
                                       perm_method=perm_method)
