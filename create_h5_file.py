import argparse
import pathlib
from preprocess import data_loader
from utils import check_functions

if __name__ == "__main__":
    # Input parameters
    parser = argparse.ArgumentParser()
    parser.add_argument('-x', '--genotype_file', type=str,
                        help='specify the name of the genotype file, absolute and relative paths are accepted, '
                             'only accept CSV, PLINK and binary PLINK files, '
                             'PLINK and binary PLINK: all required files must be in the same folder with same prefix,'
                             'for format CSV files check documentation')
    parser.add_argument('-sd', '--save_dir', type=str, default=None,
                        help='specify a directory to save newly generated H5 file. Optional, if None is specified, '
                             'H5 file will be saved in same directory as original genotype file.')

    args = vars(parser.parse_args())
    args["genotype_file"] = check_functions.check_file(args["genotype_file"])
    if pathlib.Path(args["genotype_file"]).suffix in ('.h5', '.hdf5', '.h5py'):
        raise Exception('Genotype file is already in HDF5, H5, H5PY')
    if args["save_dir"] is None:
        args["save_dir"] = pathlib.Path(args["genotype_file"]).parent
    out_file = pathlib.Path(args["genotype_file"]).with_suffix('.h5').stem
    args["save_dir"], out_file = check_functions.check_dir_paths(out_dir=args["save_dir"], out_file=out_file, prefix='')

    # load data from file
    print('Load data from file ' + str(args["genotype_file"]))
    dataset = data_loader.Genotype(genotype_file=args["genotype_file"])
    dataset.load_genotype_data()

    # save data as H5
    dataset.save_genotype_hdf5(filename=args["save_dir"].joinpath(out_file))
