"""
Author: Andrew Harris
Version: Python 3.6.4
"""
import argparse
from pathlib import Path
import os
import time
import textwrap
import logging

from pyfaidx import Fasta


def parse_fasta(f, output, window_size, chromosome, line_len):
    # ----  Load in input information ----
    fasta_file = Fasta(str(f))
    headers = fasta_file.keys()

    flag = True
    for header in list(headers):
        print( f"Parsing out {header} into {window_size} base pair windows")

        sample_seq_df = {str(header): list(fasta_file[str(header)][:].seq)}

        number_of_out_seqs = (
            len(sample_seq_df[str(header)]) // window_size + 1)

        start_pos = 0
        end_pos = window_size

        for _ in range(number_of_out_seqs):
            current_file_path = output / f"{chromosome}/{chromosome}_{start_pos}_{end_pos}.fasta"
            if flag:
                if current_file_path.is_file():
                    os.remove(current_file_path)
                else:
                    pass

            with open(current_file_path, 'a') as current_file:
                seq = textwrap.wrap("".join(list(sample_seq_df[str(header)][start_pos:end_pos])), line_len)
                if len(seq) == 0:
                    continue
                current_file.write(">{}\n".format(header))
                current_file.write("{}\n".format("\n".join(seq)))
                start_pos += window_size
                end_pos += window_size
        flag = False
    return


def main():
    # ---- Argparse Parsing ----
    parser = argparse.ArgumentParser(
        description='Break fasta file into windows of a given size.')
    parser.add_argument(
        '-i', 
        '--indir', 
        type=str, 
        action='store',
        required=True,
        help="Directory containting fasta files to run.",
    )
    parser.add_argument(
        '-o', 
        '--outdir',
        type=str,
        action='store', 
        required=True, 
        help="Directory to write output.",
    )
    parser.add_argument(
        '-w', 
        '--window_size', 
        type=int,
        action='store', 
        required=True, 
        help='Size of window to parse fasta file into.',
    )
    parser.add_argument(
        '-l', 
        '--line_length', 
        type=int,
        action='store', required=False,
        help='Fasta file line length (optional; default=80)',
        default=80,
    )
    args = parser.parse_args()

    # ---- Argument Variables ----
    INPUT_DIR = Path(args.indir)
    OUTPUT_DIR = Path(args.outdir)
    WINDOW_SIZE = args.window_size
    FASTA_LINE_LENGTH = args.line_length

    # Ensure output dir is made
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Collect files
    files_to_run = [f for f in INPUT_DIR.iterdir() if (
        f.is_file()) and (f.suffix != ".fai")]

    # Iterate through each file and parse
    for f in files_to_run:
        # Identify chromosome to be run
        CHR = f.stem.split("_")[0]

        # Make output directory before calling main function
        CHROM_OUTPUT = OUTPUT_DIR / f"{CHR}/"
        CHROM_OUTPUT.mkdir(parents=True, exist_ok=True)

        print(f"Initiating windowing for Chromosome {CHR}")

        parse_fasta(f, OUTPUT_DIR, WINDOW_SIZE, CHR, FASTA_LINE_LENGTH)
    return


if __name__ == "__main__":
    main()
