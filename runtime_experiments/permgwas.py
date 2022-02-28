import torch
from preprocessing import prepare_data as prep
from preprocessing import load_files
from perform_gwas import gwas
import pandas as pd
import time


def permgwas(args, number_of_samples=None, number_of_snps=None):
    """
    Function for runtime experiments. Perform complete gwas including data load, normal gwas and permutation based gwas.
    :param args:
    :param number_of_samples:
    :param number_of_snps:
    :return: runtime
    """
    # load data
    print('Start loading data.')
    start = time.time()
    X, y, K, covs, positions, chrom, X_index = prep.load_and_prepare_data(args, number_of_samples, number_of_snps)
    if args.maf > 0:
        if X is None:
            X = load_files.load_data(args, sample_index=X_index, snp_upper_index=number_of_snps)
        X, positions, chrom, freq = prep.maf_filter(X, positions, chrom, args.maf)
    y = y.to(args.device)
    K = K.to(args.device)
    if covs is not None:
        covs = covs.to(args.device)
    time_loaddata = time.time()
    runtime = [['load-data', time_loaddata - start]]
    print('Loaded data, elapsed time: %f s.' % (time_loaddata - start))

    # perform GWAS
    m = len(positions)
    print('Start performing GWAS on phenotype %s for %d SNPs and %d samples.' % (args.y_name, m, len(y)))
    start_gwas = time.time()
    output, runtime_gwas = gwas.gwas(args, X, y, K, covs, X_index, m)
    df = pd.DataFrame({'CHR': chrom,
                       'POS': positions,
                       'p_value': output[:, 0],
                       'test_stat': output[:, 1],
                       'SE': output[:, 2],
                       'effect_size': output[:, 3]})
    runtime = runtime + runtime_gwas
    time_gwas = time.time()
    runtime.append(['gwas_full', time_gwas - start_gwas])
    print('Done performing GWAS on phenotype %s for %d SNPs.\n '
          'Elapsed time: %f s' % (args.y_name, len(positions), time_gwas - start_gwas))
    if args.device.type != "cpu":
        with torch.cuda.device(args.device):
            torch.cuda.empty_cache()

    # perform GWAS with permutations
    if args.perm > 0:
        print('Start performing GWAS with %d permutations.' % args.perm)
        start_perm = time.time()
        adjusted_p_val, min_p_val, my_seeds, runtime_gwas = gwas.perm_gwas(args, X, y, K, output[:, 1], covs, X_index, m)
        df['adjusted_p_val'] = adjusted_p_val
        df_min = pd.DataFrame({'seed': my_seeds,
                               'min_p_val': min_p_val})
        df_min.to_csv(args.out_dir.joinpath(args.out_file + '_min_p_values.csv'), index=False)
        runtime = runtime + runtime_gwas
        time_perm = time.time()
        runtime.append(['perm_gwas_full', time_perm - start_perm])
        print('Done performing GWAS with %d permutations.'
              'Elapsed time: %f s' % (args.perm, time_perm - start_perm))

    # save p values
    df.to_csv(args.out_dir.joinpath(args.out_file + '_p_values.csv'), index=False)
    runtime.append(['total', time.time() - start])
    print('Total time: ', time.time() - start)
    return runtime
