import argparse
import torch
import numpy as np
from pathlib import Path
from preprocessing import check_functions as check
from runtime_experiments import permgwas
import pandas as pd


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
    parser.add_argument('--load_genotype', '--load_x', type=bool, default=False,
                        help='If True, genotype matrix will be loaded completely during data load. If False, genotype '
                             'matrix will be loaded batch-wise during computations of test statistics. '
                             'Default is False')
    parser.add_argument('--maf', type=int, choices=range(0, 30), default=0,
                        help='specify minor allele frequency threshold as percentage value,'
                             'optional, if not provided no maf filtering will be performed')
    parser.add_argument('--perm', type=int,
                        help='specify the number of permutations (integer value) to be performed,'
                             'optional, if not provided no permutations will be performed')
    parser.add_argument('--out_dir', type=str, default=Path.cwd().joinpath('results/runtime'),
                        help='optional, if not provided, files will be stored in folder "results/runtime" in current '
                             'directory, specify the name of the directory result-files should be stored in,'
                             'absolute and relative paths are accepted')
    parser.add_argument('--out_file', type=str,
                        help='specify NAME of result files, will be stored as NAME_p_values and NAME_min_p_values,'
                             'optional, if not provided name of phenotype will be used')
    parser.add_argument('--device', type=int, default=0,
                        help='specify GPU device to be used, default is 0')
    parser.add_argument('--batch', type=int, default=450000,
                        help='specify number of SNPs to work on simultaneously, default is 50000')
    parser.add_argument('--batch_perm', type=int, default=5500,
                        help='specify number of SNPs to work on simultaneously while using permutations, '
                             'default is 1000')
    parser.add_argument('--experiment', type=str, default='snps', help='specify type of runtime experiment, '
                        'choose "snps" to run experiments with different numbers of snps and 1000 samples,'
                        'choose "samples" to run experiments with different numbers of samples and 1mio snps, '
                        'choose "perm" to run exp with 1000 samples, 1mio snps and different number of permutations.'
                        'Default is "snps".')

    args = parser.parse_args()

    '''check if all files and directories exist'''
    if args.out_file is None:
        args.out_file = 'runtime_experiment'
    args.x = check.check_file_paths(args.x)
    args.y = check.check_file_paths(args.y)
    args.k = check.check_file_paths(args.k)
    args.cov_file = check.check_file_paths(args.cov_file)
    args.out_dir = check.make_dir_path(arg=args)
    print('Checked if all specified files exist.')

    '''check if GPU is available'''
    if torch.cuda.is_available():
        dev = "cuda:"+str(args.device)
        print('GPU is available. Perform computations on device ', dev)
    else:
        dev = "cpu"
        print('GPU is not available. Perform computations on device ', dev)
    args.device = torch.device(dev)

    if dev != "cpu":
        dev = "gpu"  # variable needed to save experiments

    if args.experiment == 'snps':
        print('Perform runtime experiment on number of snps.')
        # set parameters for runtime experiment
        list_of_snps = np.append(np.append(np.arange(10000, 100000, 5000), np.arange(100000, 1000000, 50000)),
                                 np.arange(500000, 5500000, 500000))
        number_of_samples = 1000
        for snps in list_of_snps:
            print('Start experiment with %s samples and %s snps' % (number_of_samples, snps))
            # TODO adjust batch and batch_perm to size of GPU
            runtime = []
            for run in range(3):
                print('start run ' + str(run + 1) + ' of 3')
                runtime = runtime + permgwas.permgwas(args, number_of_samples, snps)
            runtime = np.array(runtime)
            df_run = pd.DataFrame({'label': runtime[:, 0],
                                   'time': runtime[:, 1]})
            df_run.to_csv(
                args.out_dir.joinpath('runtime_snps_' + str(number_of_samples) + '_samples_' + str(snps) + '_snps_' +
                                      str(args.perm) + '_permutations_' + dev + '_.csv'), index=False)
    elif args.experiment == 'samples':
        print('Perform runtime experiment on number of samples.')
        # set parameters for runtime experiment
        list_of_samples = np.append(np.arange(100, 1000, 50), np.arange(1000, 10500, 500))
        number_of_snps = 1000000
        for samples in list_of_samples:
            print('Start experiment with %s samples and %s snps' % (samples, number_of_snps))
            args.batch = int(450000000/samples)  # TODO adjust to size of GPU
            runtime = []
            for run in range(3):
                print('start run ' + str(run + 1) + ' of 3')
                runtime = runtime + permgwas.permgwas(args, samples, number_of_snps)
            runtime = np.array(runtime)
            df_run = pd.DataFrame({'label': runtime[:, 0],
                                   'time': runtime[:, 1]})
            df_run.to_csv(
                args.out_dir.joinpath('runtime_samples_' + str(samples) + '_samples_' + str(number_of_snps) + '_snps_' +
                                      str(args.perm) + '_permutations_' + dev + '_.csv'), index=False)

    elif args.experiment == 'perm':
        print('Perform runtime experiment on number of permutations.')
        # set parameters for runtime experiment
        list_of_perm = np.append(np.arange(10, 100, 10), np.arange(100, 1000, 50))
        number_of_samples = 1000
        number_of_snps = 1000000
        for perm in list_of_perm:
            print('Start experiment with %s samples, %s snps and %s permutations'
                  % (number_of_samples, number_of_snps, perm))
            args.perm = perm
            args.batch_perm = int(550000000/(perm*number_of_samples))  # TODO adjust to size of GPU
            runtime = []
            for run in range(3):
                print('start run ' + str(run + 1) + ' of 3')
                runtime = runtime + permgwas.permgwas(args, number_of_samples, number_of_snps)
            runtime = np.array(runtime)
            df_run = pd.DataFrame({'label': runtime[:, 0],
                                   'time': runtime[:, 1]})
            df_run.to_csv(
                args.out_dir.joinpath('runtime_perm_' + str(number_of_samples) + '_samples_' + str(number_of_snps) +
                                      '_snps_' + str(args.perm) + '_permutations_' + dev + '_.csv'), index=False)
