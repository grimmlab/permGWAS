import argparse
import torch
from pathlib import Path
from preprocessing import prepare_data as prep
from preprocessing import check_functions as check
from preprocessing import load_files
from preprocessing import helper_functions
from perform_gwas import gwas
from plot import plot
import pandas as pd
import time


if __name__ == "__main__":

    # Input parameters
    parser = argparse.ArgumentParser()
    parser.add_argument('-x', '-x_file', '-genotype', type=str,
                        help='specify the name of the genotype file, absolute and relative paths are accepted, '
                             'only accept .h5, .hdf5, .h5py, .csv, PLINK and binary PLINK files, '
                             'PLINK and binary PLINK: all required files must be in the same folder with same prefix,'
                             'for format of .h5, .hdf5, .h5py and .csv files check documentation')
    parser.add_argument('-y', '-y_file', '-phenotype', type=str,
                        help='specify the name of the phenotype file, absolute and relative paths are accepted, '
                             'currently only accept .csv files')
    parser.add_argument('--y_name', '--pt_name', type=str, default='phenotype_value',
                        help='specify name of phenotype (column) to be used in phenotype file,'
                             'default is "phenotype_value"')
    parser.add_argument('--k', '--k_file', '--kinship', type=str,
                        help='optional, if not provided realized relationship kernel will be calculated'
                             'specify the name of the kinship file, absolute and relative paths are accepted,'
                             'currently only accept .csv files')
    parser.add_argument('--cov_file', '--cov', type=str,
                        help='optional, if not provided only intercept will be used as fixed effect, '
                             'specify the name of the covariates file, absolute and relative paths are accepted,'
                             'currently only accept .csv files')
    parser.add_argument('--maf', type=int, choices=range(0, 30), default=0,
                        help='specify minor allele frequency threshold as percentage value,'
                             'optional, if not provided no maf filtering will be performed')
    parser.add_argument('--load_genotype', '--load_x', action='store_true',
                        help='If True, genotype matrix will be loaded completely during data load. If False, genotype '
                             'matrix will be loaded batch-wise during computations of test statistics. '
                             'Default is False')
    parser.add_argument('--perm', type=int,
                        help='specify the number of permutations (integer value) to be performed,'
                             'optional, if not provided no permutations will be performed')
    parser.add_argument('--out_dir', type=str, default=Path.cwd().joinpath('results'),
                        help='optional, if not provided, files will be stored in folder "results" in current directory,'
                             'specify the name of the directory result-files should be stored in,'
                             'absolute and relative paths are accepted')
    parser.add_argument('--out_file', type=str,
                        help='specify NAME of result files, will be stored as p_values_NAME and min_p_values_NAME,'
                             'optional, if not provided name of phenotype will be used')
    parser.add_argument('--device', type=int, default=0,
                        help='specify GPU device to be used, default is 0')
    parser.add_argument('--batch', type=int, default=50000,
                        help='specify number of SNPs to work on simultaneously, default is 50000')
    parser.add_argument('--batch_perm', type=int, default=1000,
                        help='specify number of SNPs to work on simultaneously while using permutations, '
                             'default is 1000')
    parser.add_argument('--plot', action='store_true',
                        help='optional, creates manhattan plot')

    args = parser.parse_args()

    '''check if all files and directories exist'''
    if args.out_file is None:
        args.out_file = args.y_name + '.csv'
    args.x = check.check_file_paths(args.x)
    args.y = check.check_file_paths(args.y)
    args.k = check.check_file_paths(args.k)
    args.cov_file = check.check_file_paths(args.cov_file)
    args.out_dir, args.out_file = check.check_dir_paths(args.out_dir, args.out_file, check_again=True)

    '''check if GPU is available'''
    if torch.cuda.is_available():
        dev = "cuda:"+str(args.device)
        print('GPU is available. Perform computations on device ', dev)
    else:
        dev = "cpu"
        print('GPU is not available. Perform computations on device ', dev)
    args.device = torch.device(dev)

    '''load data'''
    print('Checked if all specified files exist. Start loading data.')
    start = time.time()
    X, y, K, covs, positions, chrom, X_index = prep.load_and_prepare_data(args)
    if args.maf > 0:
        if X is None:
            X = load_files.load_genotype_matrix(args.x, sample_index=X_index)
        X, positions, chrom, freq = prep.use_maf_filter(X, positions, chrom, args.maf)
    elif X is not None:
        freq = prep.get_maf(X)
    print('Loaded data, elapsed time: %f s.' % (time.time()-start))
    y = y.to(args.device)
    K = K.to(args.device)
    if covs is not None:
        covs = covs.to(args.device)

    n_samples = len(y)
    n_snps = len(positions)

    '''perform GWAS'''
    print('Start performing GWAS on phenotype %s for %d SNPs and %d samples.' % (args.y_name, n_snps, n_samples))
    start_gwas = time.time()
    if X is not None:
        output, v_g, v_e, _ = gwas.gwas(args, X, y, K, covs, X_index, n_snps)
    else:
        output, v_g, v_e, freq = gwas.gwas(args, X, y, K, covs, X_index, n_snps)
    df = pd.DataFrame({'CHR': chrom,
                       'POS': positions,
                       'p_value': output[:, 0],
                       'test_stat': output[:, 1],
                       'maf': freq,
                       'SE': output[:, 2],
                       'effect_size': output[:, 3]})
    print('Done performing GWAS on phenotype %s for %d SNPs.\n'
          'Elapsed time: %f s' % (args.y_name, n_snps, time.time()-start_gwas))
    if args.device.type != "cpu":
        with torch.cuda.device(args.device):
            torch.cuda.empty_cache()

    '''perform GWAS with permutations'''
    if args.perm is not None:
        print('Start performing GWAS with %d permutations.' % args.perm)
        start_perm = time.time()
        adjusted_p_val, min_p_val, my_seeds = gwas.perm_gwas(args, X, y, K, output[:, 1], covs, X_index, n_snps)
        df['adjusted_p_val'] = adjusted_p_val
        df_min = pd.DataFrame({'seed': my_seeds,
                               'min_p_val': min_p_val})
        df_min.to_csv(args.out_dir.joinpath('min_p_values_' + args.out_file), index=False)
        print('Done performing GWAS with %d permutations.'
              'Elapsed time: %f s' % (args.perm, time.time()-start_perm))
    else:
        min_p_val = None
    '''save p values'''
    df.to_csv(args.out_dir.joinpath('p_values_' + args.out_file), index=False)
    print('Total time: ', time.time()-start)

    '''create manhattan plot'''
    if args.plot is True:
        if args.perm is None:
            plot.manhattan(df, 'p_value',
                           args.out_dir.joinpath('manhattan_' + Path(args.out_file).with_suffix('.png').as_posix()))
        else:
            plot.manhattan(df, 'p_value',
                           args.out_dir.joinpath('manhattan_' + Path(args.out_file).with_suffix('.png').as_posix()),
                           df_min)

    '''get summary statistics'''
    helper_functions.get_summary_stats(arguments=args, samples=n_samples, snps=n_snps, v_g=v_g, v_e=v_e,
                                       min_p_val=min_p_val)