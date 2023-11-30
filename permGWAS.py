# run the script here
import argparse
import pathlib

from utils import check_functions
import perform_gwas


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-x', '--genotype_file', type=str, default=None,
                        help='Specify the full path to the genotype file, absolute and relative paths are accepted, '
                             'only accept .h5, .hdf5, .h5py, .csv, PLINK and binary PLINK files, '
                             'PLINK and binary PLINK: all required files must be in the same folder with same prefix. '
                             'See documentation for correct format.')
    parser.add_argument('-y', '--phenotype_file', type=str, default=None,
                        help='Specify the full path to the phenotype file, absolute and relative paths are '
                             'accepted, only accept .csv, .txt and .pheno files. See documentation for correct format.')
    parser.add_argument('-trait', '--trait', '--y_name', nargs='+', type=str, default=['phenotype_value'],
                        help='Specify the name of phenotype (column) to be used in phenotype file,'
                             'default is "phenotype_value". You can run permGWAS on several phenotypes one after '
                             'another if they are in the same phenotype_file. Juste name the phenotypes, '
                             'e.g. --trait pheno1 pheno2 if you want to use all available traits use --trait all')
    parser.add_argument('-k', '--kinship_file', '--k', '--kinship', type=str, default=None,
                        help='Specify the the full path to the kinship file, absolute and relative paths are accepted,'
                             'only accept .csv and .h5/.h5py/.hdf5 files. See documentation for correct format. '
                             'Optional, if not provided realized relationship kernel will be calculated')
    parser.add_argument('-cov', '--covariate_file', '--cov', '--cov_file', type=str, default=None,
                        help='Specify the full path to the covariates file, absolute and relative paths are accepted,'
                             'currently only accept .csv files. Optional, if not provided only intercept will be used '
                             'as fixed effect.')
    parser.add_argument('-cov_list', '--covariate_list', nargs='+', type=str, default=None,
                        help='Specify the covariates (column headers) to use from the covariates file. Optional, if '
                             'not provided, will use all available columns as covariates.')
    parser.add_argument('-maf', '--maf_threshold', '--maf', type=int, choices=range(0, 31), default=0,
                        help='Specify minor allele frequency threshold as percentage value. '
                             'Optional, if not provided no maf filtering will be performed.')
    parser.add_argument('-load_genotype', action='store_true',
                        help='If used, genotype matrix will be completely loaded from file during preprocessing. '
                             'Otherwise load genotype batch-wise during computations of test statistics. '
                             'Batch-wise loading is only possible, if kinship file is provided. Default is False')
    parser.add_argument('-config', '--config_file', type=str, default=None,
                        help='Specify the full path to the yaml config file. Specify all required arguments to use in '
                             'this config file and just give the config file instead of all required parameters. '
                             'For more info regarding the required format see the documentation.')
    parser.add_argument('-model', type=str, default='lmm',
                        help='Specify the model to use for GWAS. Currently only lmm (linear mixed model) is '
                             'implemented.')
    parser.add_argument('-out_dir', '--out_dir', type=str, default=pathlib.Path.cwd().joinpath('results'),
                        help='Specify the name of the directory result-files should be stored in,'
                             'absolute and relative paths are accepted. Optional, if not provided, files will be '
                             'stored in folder "results" in current directory,')
    parser.add_argument('-out_file', '--out_file', type=str, default=None,
                        help='Specify NAME of result files, will be stored as p_values_NAME and min_p_values_NAME,'
                             'optional, if not provided name of phenotype will be used. If you run permGWAS with '
                             'several phenotypes, will always use name of phenotype.')
    parser.add_argument('-disable_gpu', action='store_true',
                        help='If used, GPUs will be disabled and only CPUs will be used for computations.')
    parser.add_argument('-device', '--device', type=int, default=0,
                        help='Specify GPU device to be used, default is 0.')
    parser.add_argument('-perm', '--perm', type=int, default=0,
                        help='Specify the number of permutations (integer value) to be performed, optional, if not '
                             'provided no permutations will be performed')
    parser.add_argument('-perm_method', type=str, default='x',
                        help='Specify the method to use for permutations: x or y,'
                             'for x permute fixed effects matrix including SNP of interest, which is equivalent to '
                             'permuting the phenotype and the covariance matrix; for y permute only the phenotype '
                             'vector as in permGWAS Version1. Default is x.')
    parser.add_argument('-adj_p_value', action='store_true',
                        help='If used, will additionally compute adjusted permutation-based p-values for each SNP.')
    parser.add_argument('-batch', '--batch_size', '--batch', type=int, default=50000,
                        help='Specify number of SNPs to work on simultaneously, default is 50000')
    parser.add_argument('-batch_perm', '--perm_batch_size', '--batch_perm', type=int, default=1000,
                        help='Specify number of SNPs to work on simultaneously, default is 1000')
    parser.add_argument('-mplot', '--manhattan', '--plot', action='store_true',
                        help='optional, creates manhattan plot')
    parser.add_argument('-qqplot', '--qqplot', action='store_true',
                        help='optional, creates QQ-plot')
    parser.add_argument('-not_add', '--not_add', action='store_true',
                        help='optional, use if genotype has different encoding.')
    args = vars(parser.parse_args())
    # check config file
    args = check_functions.check_all_arguments(args=args)
    phenotypes = args["trait"]

    # run pipeline
    for trait in phenotypes:
        print('Working on phenotype ', trait)
        args["trait"] = trait
        args = check_functions.check_output_files(args=args)
        print('Checked if all specified files exist.')
        try:
            perform_gwas.run(**args)
            args["out_file"] = None
        except Exception as exc:
            print("Failure when running permGWAS2.0")
            print(exc)
            continue
