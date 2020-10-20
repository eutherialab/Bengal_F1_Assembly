"""
Author: Andrew Harris
Python 3.8.3

This script combines csv files into a single lst file.
The .lst file is used with Seqtk to subset out reads.
"""
import argparse
from pathlib import Path

import pandas as pd


def main():
    # ---- Argparse Variables ----
    parser = argparse.ArgumentParser(
        description="Concatenate subtype read name files.")
    # Required input arguments
    parser.add_argument('--cross', type=str, required=True,
                        action='store', help='Input file cross name')
    parser.add_argument('--inputdir', type=str, required=True, action='store',
                        help="Pathway to directory containing input files (can be one or more)")
    args = parser.parse_args()

    CROSS = args.cross
    INPUT_DIR = Path(args.inputdir)

    # ---- Main body of script ----

    # Collect files in INPUT_DIR
    read_name_files = [f for f in INPUT_DIR.iterdir() if f.is_file() and (
        f.stem[0] != ".") and (f.suffix != ".lst")]

    # DataFrame collection list
    file_dfs = []

    # Load files into pandas DataFrames and add to file_dfs = []
    for rn_file in read_name_files:
        df = pd.read_csv(rn_file, sep=",", header=None)
        file_dfs.append(df)
        continue

    # Concatenate all dataframes
    concat_dfs = pd.concat(file_dfs)
    concat_dfs.drop_duplicates(keep='first', inplace=True)

    # Output to a single file
    output_file = INPUT_DIR / f"{CROSS}.lst"
    concat_dfs.to_csv(output_file, sep=",", index=False, header=False)
    print(f"Number of reads in output file: {len(concat_dfs)}")
    return


if __name__ == "__main__":
    main()
