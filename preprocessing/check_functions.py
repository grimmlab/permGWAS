from pathlib import Path


def check_file_paths(file_path):
    """
    check if specified input file exists, raises Exception if not
    :param file_path: absolute or relative path
    :return: path object if file exists
    """
    if file_path is None:
        return file_path
    else:
        my_path = Path(file_path)
        if my_path.is_file():
            return my_path
        else:
            raise FileNotFoundError('There is no file ', file_path)


def check_dir_paths(out_dir: str, out_file: str, check_again=False):
    """
    check if directory for result files exists, if not, create directory.
    Then check if result files already exist, if they already exist, raises Exception
    :param out_dir: directory to save result files
    :param out_file: result file
    :param check_again: if True, rename result file via adding (i) to the end of the file
    :return: path object
    """
    my_path = Path(out_dir)
    if my_path.is_dir():
        if my_path.joinpath(out_file).exists():
            if not check_again:
                raise FileExistsError('The file %s already exists in chosen directory %s' % (out_file, out_dir))
            else:
                i = 1
                new_file = out_file.split('.')[0] + '(' + str(i) + ').csv'
                new_path = my_path.joinpath(new_file)
                while new_path.exists():
                    i += 1
                    new_file = out_file.split('.')[0] + '(' + str(i) + ').csv'
                    new_path = my_path.joinpath(new_file)
                print('The file %s already exists in chosen directory %s. Changed filename to %s.'
                      % (out_file, out_dir, new_file))
        else:
            new_file = out_file
    else:
        new_file = out_file
        my_path.mkdir(parents=True, exist_ok=True)
    return my_path, new_file
