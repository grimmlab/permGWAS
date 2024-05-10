# create Manhattan and QQ-plots
import pandas as pd
import pathlib
import argparse

from utils import check_functions
from postprocess import plot_functions

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-p_val', '--p_value_file', type=str, default=None,
                        help='Specify the full path to the p_value file, absolute and relative paths are accepted, '
                             'only accept .csv files. p_value files must at least contain chromosome ids (CHR), '
                             'position ids (POS) and corresponding p_values (p_value).')
    parser.add_argument('-min_p_val', '--min_p_value_file', type=str, default=None,
                        help='Optional, specify the full path to the file containing minimal p-values in order to '
                             'compute permutation-based thresholds, absolute and relative paths are accepted, '
                             'only accept .csv files.')
    parser.add_argument('-mplot', '--manhattan', action='store_true',
                        help='optional, creates manhattan plot')
    parser.add_argument('-qqplot', action='store_true',
                        help='optional, creates QQ-plot')
    parser.add_argument('-out_dir', type=str, default=None,
                        help='Specify the name of the directory plots should be stored in,'
                             'absolute and relative paths are accepted. Optional, if not provided, files will be '
                             'stored in same folder as p_value file.')
    parser.add_argument('-out_file', type=str, default=None,
                        help='Specify NAME of plots, will be stored as manhattan_NAME.png or qq_plot_NAME.png,'
                             'optional, if not provided name of p_value file will be used.')
    parser.add_argument('-sig_level', type=int, default=5,
                        help='Significance level (percentage values) to compute threshold for Manhattan plot. '
                             'Optional, default is 5.')
    args = vars(parser.parse_args())

    args["p_value_file"] = check_functions.check_file(args["p_value_file"])
    if args["min_p_value_file"] is not None:
        args["min_p_value_file"] = check_functions.check_file(args["min_p_value_file"])
    if args["out_dir"] is None:
        args["out_dir"] = pathlib.Path(args["p_value_file"]).parent
    if args["out_file"] is None:
        args["out_file"] = pathlib.Path(args["p_value_file"]).stem

    df = pd.read_csv(args["p_value_file"])
    if not {'CHR', 'POS', 'p_value'}.issubset(df.columns):
        raise Exception('Cannot create Manhattan plot; need CHR, POS and p_value in DataFrame.')

    if args["manhattan"]:
        out_dir, out_file = check_functions.check_dir_paths(out_dir=args["out_dir"], out_file=args["out_file"],
                                                            prefix='manhattan_')
        print('Save Manhattan plot with significance level of %d.' % args["sig_level"])
        if args["min_p_value_file"] is not None:
            df_min = pd.read_csv(args["min_p_value_file"])
            if not 'min_p_val' in df_min.columns:
                raise Exception('Cannot compute permutation-based threshold, need min_p_val in DataFrame.')
            min_p_val = df_min['min_p_val'].values
        else:
            min_p_val = None
        plot_functions.manhattan_plot(df=df, data_dir=out_dir, filename=out_file,
                                      min_p_values=min_p_val, sig_level=args["sig_level"])

    if args["qqplot"]:
        out_dir, out_file = check_functions.check_dir_paths(out_dir=args["out_dir"], out_file=args["out_file"],
                                                            prefix='qq_plot_')
        print('Save QQ-plot.')
        plot_functions.qq_plot(p_values=df['p_value'].values, data_dir=out_dir, filename=out_file)
