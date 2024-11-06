import pathlib
import torch
import pandas as pd
from utils import helper_functions
import models


def check_all_arguments(args: dict) -> dict:
    """
    Check user specified arguments for plausibility and turn all file paths to pathlib.Path objects
    :param args:
    :return:
    """
    if args["config_file"] is not None:
        args = helper_functions.parse_config_file(args=args)
    del args["config_file"]
    # check if specified files exist
    args["genotype_file"] = check_file(filepath=args["genotype_file"])
    args["phenotype_file"] = check_file(filepath=args["phenotype_file"])
    args["kinship_file"] = check_file(filepath=args["kinship_file"])
    args["covariate_file"] = check_file(filepath=args["covariate_file"])
    if args["trait"] is None:
        args["trait"] = 'phenotype_value'
    elif (args["trait"] == 'all') or (args["trait"] == ['all']):
        print('Will perform computations on all available phenotypes.')
        args["out_file"] = None
        suffix = args["phenotype_file"].suffix
        if suffix == ".csv":
            df = pd.read_csv(args["phenotype_file"], index_col=0)
        # load PHENO or TXT
        elif suffix == ".txt":
            df = pd.read_csv(args["phenotype_file"], index_col=0, sep=" ")
        elif suffix == ".pheno":
            df = pd.read_csv(args["phenotype_file"], index_col=0, sep=" ")
            if 'FID' in df.columns:
                df.drop(columns='FID', inplace=True)
            if 'IID' in df.columns:
                df.drop(columns='IID', inplace=True)
        else:
            raise Exception('Only accept .txt, .pheno or .csv phenotype files.')
        args["trait"] = df.columns.tolist()
    elif isinstance(args["trait"], str):
        args["trait"] = [args["trait"]]
    elif isinstance(args["trait"], list):
        args["out_file"] = None
    else:
        raise Exception('Something is wrong with the trait name. Please check again.')
    # sanity checks for fast loading and batch-wise loading
    if args["kinship_file"] is None:
        args["load_genotype"] = True
    if args["genotype_file"].suffix not in ('.h5', '.hdf5', '.h5py'):
        args["load_genotype"] = True
    # check gpu
    if torch.cuda.is_available() and not args["disable_gpu"]:
        dev = "cuda:" + str(args["device"])
        print('GPU is available. Perform computations on device ', dev)
    else:
        dev = "cpu"
        print('GPU is not available. Perform computations on device ', dev)
    del args["disable_gpu"]
    args["device"] = torch.device(dev)
    # check model
    if args["model"] is None:
        args["model"] = 'lmm'
    if args["model"] not in models.__all__:
        raise NotImplementedError('Specified model not implemented')

    # sanity checks
    if args["maf_threshold"] is None:
        args["maf_threshold"] = 0
    if isinstance(args["covariate_list"], str):
        args["covariate_list"] = [args["covariate_list"]]
    # check permutation method
    if args["perm"] is None:
        args["perm"] = 0
    if args["perm"] > 0:
        if args["perm_method"] not in ('x', 'y'):
            raise NotImplementedError(' Can only perform permutation methods x and y. Please check again.')
    if args["adj_p_value"] and args["perm"] == 0:
        raise Exception('Can not compute adjusted p-values with 0 permutations. Please check again.')
    return args


def check_output_files(args: dict) -> dict:
    # check output directory and file
    if args["out_file"] is None:
        args["out_file"] = args["trait"] + '.csv'
    if args["out_dir"] is None:
        args["out_dir"] = pathlib.Path.cwd().joinpath('results')
    args["out_dir"], args["out_file"] = check_dir_paths(args["out_dir"], args["out_file"])
    return args


def check_file(filepath: str):
    """
    Check if specified file exists

    :param filepath: full path to file

    :return: path to file as Path object
    """
    if filepath is None:
        return None
    else:
        filepath = pathlib.Path(filepath)
        if filepath.is_file():
            return filepath
        else:
            raise FileNotFoundError('There is no file ', filepath.as_posix())


def check_dir_paths(out_dir: str, out_file: str, prefix: str = 'p_values_') -> (pathlib.Path, pathlib.Path):
    """
    Check if directory for result files exists, if not, create directory.
    Then check if result files already exist, if they already exist, rename result file by adding (i) to the
    end of the file

    :param out_dir: directory to save result files
    :param out_file: result file
    :param prefix: prefix to use when checking for existing files, default is p_values_

    :return: path object
    """
    my_path = pathlib.Path(out_dir)
    if prefix in ('manhattan_', 'qq_plot_'):
        suffix = '.png'
    elif prefix == '':
        suffix = '.h5'
    else:
        suffix = '.csv'
    out_file = pathlib.Path(out_file).with_suffix(suffix).as_posix()
    if my_path.is_dir():
        if my_path.joinpath(prefix + out_file).exists():
            if suffix == '.h5':
                raise Exception('File %s already exists in chosen directory %s.' % (out_file, out_dir))
            i = 1
            new_file = pathlib.Path(out_file).with_suffix('').as_posix() + '(' + str(i) + ')' + suffix
            new_path = my_path.joinpath(prefix + new_file)
            while new_path.exists():
                i += 1
                new_file = pathlib.Path(out_file).with_suffix('').as_posix() + '(' + str(i) + ')' + suffix
                new_path = my_path.joinpath(prefix + new_file)
            print('The file %s already exists in chosen directory %s. Changed filename to %s.'
                  % (prefix + out_file, out_dir, prefix + new_file))
        else:
            new_file = out_file
    else:
        new_file = out_file
        my_path.mkdir(parents=True, exist_ok=True)
    return my_path, new_file
