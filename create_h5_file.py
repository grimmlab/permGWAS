import argparse
import h5py
from preprocessing import prepare_data as prep
from preprocessing import check_functions as check


if __name__ == "__main__":
    # Input parameters
    parser = argparse.ArgumentParser()
    parser.add_argument('-x', '-x_file', '-genotype', type=str,
                        help='specify the name of the genotype file, absolute and relative paths are accepted, '
                             'only accept CSV, PLINK and binary PLINK files, '
                             'PLINK and binary PLINK: all required files must be in the same folder with same prefix,'
                             'for format CSV files check documentation')
    parser.add_argument('-sd', '-save_dir', type=str, default=None,
                        help='specify a directory to save newly generated H5 file. Optional, if None is specified, '
                             'H5 file will be saved in same directory as original genotype file.')

    args = parser.parse_args()
    args.x = check.check_file_paths(args.x)
    if args.x.suffix in ('.h5', '.hdf5', '.h5py'):
        raise Exception('Genotype file is already in HDF5, H5, H5PY')
    if args.sd is None:
        x_file = args.x.with_suffix('').as_posix() + '.h5'
    else:
        args.sd, _ = check.check_dir_paths(args.sd, args.x.stem + '.h5')
        x_file = args.sd.joinpath(args.x.stem + '.h5')

    '''load data from file'''
    print('Load data from file ' + str(args.x))
    sample_ids, positions, chromosomes = prep.load_files.load_genotype_ids(args.x)
    X = prep.load_files.load_genotype_matrix(args.x)

    '''save data as H5'''
    print('Have data.\nSave data as ' + str(x_file) + '.\nThis might take some time.')
    with h5py.File(x_file, 'w') as f:
        f.create_dataset('sample_ids', data=sample_ids.astype(bytes), chunks=True, compression="gzip")
        f.create_dataset('chr_index', data=chromosomes.astype(bytes), chunks=True, compression="gzip")
        f.create_dataset('position_index', data=positions.astype(int), chunks=True, compression="gzip")
        f.create_dataset('snps', data=X, chunks=True, compression="gzip", compression_opts=7)
    print('Done saving H5 file.')
