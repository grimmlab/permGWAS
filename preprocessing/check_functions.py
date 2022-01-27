from pathlib import Path


def check_file_paths(file_path):
    """
    check if specified files exist
    :param file_path: absolute or relative path
    :return: path object
    """
    if file_path is None:
        return file_path
    else:
        my_path = Path(file_path)
        if my_path.is_file():
            return my_path
        else:
            raise FileNotFoundError('There is no file ', file_path)


def check_dir_paths(arg):
    """
    check if directory for result files exists, if not, create directory, check if result files already exist
    :param arg: argparse.Namespace
    :return: path object
    """
    my_path = Path(arg.out_dir)
    if my_path.is_dir():
        if my_path.joinpath(arg.out_file+'_p_values.csv').exists():
            raise FileExistsError('The file %s already exists in chosen directory %s'
                                  % (arg.out_file+'_p_values.csv', arg.out_dir))
        else:
            if my_path.joinpath(arg.out_file+'_min_p_values.csv').exists() and arg.perm is not None:
                raise FileExistsError('The file %s already exists in chosen directory %s'
                                      % (arg.out_file+'_min_p_values.csv', arg.out_dir))
    else:
        my_path.mkdir(parents=True, exist_ok=True)
    return my_path
