"""
Author: Andrew Harris
Python Version: 3.6, 3.7, 3.8

This script will take in pre-windowed multi-alignment files of any format specified at (https://biopython.org/wiki/AlignIO) 
and will calculate the p-distance of each sample from the indicated reference sample provided.

*** Make sure to install dependencies Pandas and Biopython in your python enviornment.
"""
import glob
import sys
import os
from pathlib import Path
import argparse
import time
import datetime
import logging
import pathlib

import numpy as np
import pandas as pd
from Bio import AlignIO

# ------------------- Helper functions -------------------

def _get_chromosome_dirs(input_directory):
    """Collect chromosome directories"""
    dirs = []
    for d in input_directory.iterdir():
        if not d.is_dir():
            continue
        # Just in case user re-runs and
        # does not delete output files
        elif d.name == 'logs':
            continue
        elif d.name == 'p_distance_output':
            continue
        else:
            dirs.append(d)
    return dirs


def base_check(base):
        """Check if current base is A, T, C, or G. 
        All other symbols will be concidered as N """
        bases = ['A', 'T', 'C', 'G']
        if base not in bases:
            return False
        else:
            return True


def p_distance_calc(current_sample_seq, reference_sample_seq, missing_data_threshold):
    """This function takes in the current sample sequence and the reference 
    sample sequence and calculates the p-distance per window"""
    snp_count = 0
    same_base = 0
    nan_count = 0

    logging.debug(f"Current Sample: {current_sample_seq[:10]}--{current_sample_seq[-10:]}")
    logging.debug(f"Reference Sample: {reference_sample_seq[:10]}--{reference_sample_seq[-10:]}\n")

    for curr_samp_base, ref_base in zip(current_sample_seq, reference_sample_seq):
        curr_samp_check = base_check(curr_samp_base)
        ref_base_check = base_check(ref_base)
        if not curr_samp_check or not ref_base_check:
            nan_count += 1
        elif curr_samp_base != ref_base:
            snp_count += 1
        else:
            same_base += 1
            continue
    assert snp_count + same_base + nan_count == len(current_sample_seq)
    if nan_count / (snp_count + same_base + nan_count) >= missing_data_threshold:
        p_distance = np.nan
    else:
        p_distance = float(snp_count / len(current_sample_seq))

    return p_distance


def per_sample_calc(
    reference_sample, 
    alignment_df, 
    file_name, 
    file_df,
    missing_data_threshold,
    outfile_index
):
    # ["Chromosome", "Start", "Stop"]
    split_file_name = file_name.split("_")
    chromosome = str(split_file_name[0])
    start = int(split_file_name[1])
    stop = int(split_file_name[2])

    # Set Chromosome, Start, and Stop in file_df
    file_df.at[outfile_index, "Chromosome"] = chromosome
    file_df.at[outfile_index, "Start"] = start
    file_df.at[outfile_index, "Stop"] = stop

    reference_sample_sequence = list(alignment_df[reference_sample])

    logging.debug(f"Samples: {alignment_df.columns.to_list()}")

    for current_sample in alignment_df.columns:
        if current_sample == reference_sample:
            file_df.at[outfile_index, str(reference_sample)] = 0
            continue
        else:
            current_sample_seq = list(alignment_df[current_sample])
            p_distance = p_distance_calc(
                current_sample_seq, 
                reference_sample_sequence,
                missing_data_threshold,
            )
            file_df.at[outfile_index, str(current_sample)] = p_distance
            continue

    return file_df


# ------------------------- BODY -------------------------

def p_distance_calculator(
        input_directory=None,
        output_directory=None,
        reference_sample=None,
        project_id=None,
        missing_data_threshold=None,
        window_size=None,
):

    # set Paths()
    input_directory = Path(input_directory)
    if not output_directory:
        output_directory = Path(input_directory) / "p_distance_output"
    else:
        output_directory = Path(output_directory)

    # Make output directory
    output_directory.mkdir(parents=True, exist_ok=True)

    # Set Date
    DATE = str(datetime.date.today()).replace('-', '_')

    # Initialize log file
    log_path = output_directory / "logs/"
    log_path.mkdir(parents=True, exist_ok=True)
    log_file = log_path / f"{project_id}_{DATE}.log"
    logging.basicConfig(
        filename=log_file,
        level=logging.DEBUG,
        filemode='w',
        format='%(asctime)s %(levelname)-8s %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    """Collect chromosome directories into list"""
    chromosome_dirs = _get_chromosome_dirs(input_directory)

    """Iterate through each chromosome directory in reverse"""
    """and calculate p-distance per-file"""
    for chrom_dir in reversed(chromosome_dirs):
        print(f"---- Running {chrom_dir.name} ----")
        logging.debug(f"---- Running {chrom_dir.name} ----")

        windowed_files = [f for f in chrom_dir.iterdir() if f.is_file()]
        file_df = None
        outfile_index = 0
        columns_flag = True
        file_tracker = 1

        """ Iterate through each file in the chrom directory"""
        for f in windowed_files:
            if file_tracker % 10 == 0:
                print(f"-- Completed {file_tracker} of {len(windowed_files)} files --")
            try:
                read_alignment = AlignIO.read(f, 'fasta')
            except ValueError:
                raise ValueError("Faulty input file -- exit")

            alignment_dict = {str(read.id): list(read.seq)
                              for read in read_alignment}
            sample_names = [sample_name.name for sample_name in read_alignment]
            alignment_df = pd.DataFrame(
                data=alignment_dict, columns=sample_names)

            # Set up output_df
            if columns_flag:
                column_names = list(alignment_df)
                column_names = ["Chromosome", "Start", "Stop"] + column_names
                file_df = pd.DataFrame(columns=column_names)
                columns_flag = False

            logging.debug(f"-- Calculate p-distance for file {f.name} --")

            file_df = per_sample_calc(
                reference_sample,
                alignment_df,
                f.stem,
                file_df,
                missing_data_threshold,
                outfile_index
            )
            outfile_index += 1
            file_tracker += 1
            continue

        print(f"---- Completed {file_tracker} of {len(windowed_files)} files ----")

        # sort output file_df
        file_df.sort_values(by=['Start'], inplace=True)

        # Output file_df
        output_file = output_directory / f"{chrom_dir.name}_{project_id}_{DATE}.csv"
        file_df.to_csv(output_file, sep=',', index=False)

    return


def concat_output_files(
    input_directory=None,
    output_directory=None,
    reference_sample=None,
    project_id=None,
):
    """This function will concatenate all chromosome output files 
    into a single file for viewing on P-distance Graphing Tool.
    It will also remove the reference sample data as it is only 
    0's and uniformative."""
    logging.debug(f"Starting concatenation process")

    DATE = str(datetime.date.today()).replace('-', '_')

    # set Paths()
    input_directory = Path(input_directory)
    if not output_directory:
        output_directory = Path(input_directory) / "p_distance_output"
    else:
        output_directory = Path(output_directory)

    # Collect output files
    output_files = [f for f in output_directory.iterdir() if f.is_file()]
    pprint_list = "\n".join([str(f) for f in output_files])
    logging.debug(f"Files to concatenate:\n{pprint_list}\n")

    # Read files into dataframe
    dfs_to_concat = [pd.read_csv(f) for f in output_files]

    # Concatenate the files into a single directory
    concat_dfs = pd.concat(dfs_to_concat)

    # Drop reference column from final dataset
    concat_dfs.drop(columns=reference_sample, inplace=True)

    # Create output filename and output to file
    output_filename = output_directory / f"{reference_sample}_{project_id}_{DATE}.csv"
    concat_dfs.to_csv(output_filename, index=False)
    return


def main():
    # ---- Argeparse Parsing ----
    parser = argparse.ArgumentParser(
        description='P-Distance calculation and graphing.')
    parser.add_argument(
        '--pdistance',
        action='store_true', 
        default=False, 
        help="Run just p-distance calculation"
    )
    parser.add_argument(
        '--concat',
        action='store_true', 
        default=False, 
        help="Run just file concatenation"
    )
    parser.add_argument(
        '-r',
        '--reference',
        type=str,
        action='store', 
        required=True, 
        help="Reference sample name."
    )
    parser.add_argument(
        '-i',
        '--input',
        type=str,
        action='store', 
        required=True, 
        help="Directory containing subdirectories of per-chromosome windowed fasta files (output of fasta_windower.py)"
    )
    parser.add_argument(
        '-p', 
        '--project_id', 
        type=str,
        action='store', 
        required=True, 
        help='Provide a project ID to track this specific run through the Dash Graphing application'
    )
    parser.add_argument(
        '-c', 
        '--cutoff', 
        type=float,
        action='store', 
        required=True, 
        help='Provide the maximum frequency of missing data points allowed in a given window. (0.0 - 1.0)'
    )
    parser.add_argument(
        '-w',
        '--window_size',
        type=int,
        action='store', 
        required=True,
        help='Provide the window size used.'
    )
    # Optional Arguments
    parser.add_argument(
        '-o',
        '--output',
        type=str,
        action='store',
        help='Indicate alternate output location if you do not want the output to be placed in the current working directory.'
    )
    args = parser.parse_args()

    # ---- Argparse Variables ----
    run_pdist = args.pdistance
    run_concat = args.concat
    input_directory = args.input
    output_directory = args.output
    reference_sample = args.reference
    project_id = args.project_id
    missing_data_threshold = args.cutoff
    window_size = args.window_size
    
    # If --pdistance argument provided, just run p-dist-calc
    if run_pdist:
        p_distance_calculator(
            input_directory,
            output_directory,
            reference_sample,
            project_id,
            missing_data_threshold,
            window_size,
        )

    # If --concat argument provided, just run p-dist-calc
    if run_concat:
        concat_output_files(
            input_directory,
            output_directory,
            reference_sample,
            project_id,
        )
    if not run_pdist and not run_concat:
        # Step 1: Run p-distance calculations on windowed files
        p_distance_calculator(
            input_directory,
            output_directory,
            reference_sample,
            project_id,
            missing_data_threshold,
            window_size,
        )
        # Step 2: Concatenate output files
        concat_output_files(
            input_directory,
            output_directory,
            reference_sample,
            project_id,
        )
    return


if __name__ == "__main__":
    main()
