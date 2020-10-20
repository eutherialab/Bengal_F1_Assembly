"""
Author: Andrew Harris
Python 3.8.3

This script takes in the RepeatMasker.masked output file 
and provides summary statistics on repetitive content,
and average read lengths. 

Input file: RepeatMasker output .masked file
"""
import argparse
import logging
from pathlib import Path

from pyfaidx import Fasta


def check_input_filetype(input_file):
    """Make sure that input file is a .masked file"""
    try:
        assert(input_file.suffix == ".masked")
    except AssertionError:
        print("Input file is not a .masked file, please correct and run again.")
        exit()


def output_nonrepetitive_read_names(outputdir, cross, nonrepetitive_read_name_list):
    nr_output_file = outputdir / f"{cross}_non-repetitive_read_names.lst"
    with open(nr_output_file, "w") as oh:
        for name in nonrepetitive_read_name_list:
            oh.writelines(f"{name}\n")
            continue
    return

def main():
    parser = argparse.ArgumentParser(description='Read statistic script for PHA')
    parser.add_argument('--input', type=str, action='store', required=True, help="RepeatMasker.masked file")
    parser.add_argument('--outdir', type=str, action='store', required=True, help='Output directory')
    parser.add_argument('--cross', type=str, action='store', required=True, help='Input data Cross (i.e. LilBubxPbe14)')
    args = parser.parse_args()

    # Set args to variables
    input_file = Path(args.input)
    outputdir = Path(args.outdir)
    cross = args.cross

    check_input_filetype(input_file)

    # Make output dir if not already created
    outputdir.mkdir(parents=True, exist_ok=True)

    # Initital log file
    logfile_name = outputdir / f"{cross}_count_N_in_masked_output.log"
    logging.basicConfig(filename=logfile_name, level=0, filemode='w', format='')

    # Create collection list dict={"read_name": int(read_length),}
    below_50_percent = {}
    above_50_percent = {}
    non_repetitive = {}

    # Collect # reads below 10kb and 5kb
    below_five_kb = 0
    below_ten_kb = 0

    # Stats collection
    total_reads = 0
    total_nonrepetitive = 0
    repetitive_read_length_list = []
    nonrepetitive_read_length_list = []
    nonrepetitive_read_name_list = []
    all_read_lengths = []

    # Load into pyfaidx Fasta object
    input_fasta = Fasta(str(input_file), "fasta")

    # Iterate through each read and record stats
    for read_name in input_fasta.keys():
        read_length = len(input_fasta[read_name])
        all_read_lengths.append(read_length)
        # below 5kb or 10kb
        if read_length < 5000:
            below_five_kb += 1
            below_ten_kb += 1
            pass
        elif read_length < 10000:
            below_ten_kb += 1
        n_count = str(input_fasta[read_name]).count("N")
        percent_N = (n_count / read_length) * 100
        
        if below_50_percent.get(read_name) != None:
            continue
        elif above_50_percent.get(read_name) != None:
            continue
        elif non_repetitive.get(read_name) != None:
            continue
        elif percent_N == 0:
            total_nonrepetitive += 1
            total_reads += 1
            below_50_percent[read_name] = (percent_N)
            non_repetitive[read_name] = (percent_N)
            nonrepetitive_read_length_list.append(read_length)
            nonrepetitive_read_name_list.append(read_name)
            continue
        elif percent_N < 50:
            total_reads += 1
            below_50_percent[read_name] = (percent_N)
            repetitive_read_length_list.append(read_length)
            continue
        elif percent_N >= 50:
            total_reads += 1
            above_50_percent[read_name] = (percent_N)
            repetitive_read_length_list.append(read_length)
        else:
            raise AssertionError

    # Make sure all reads were accounted for
    assert len(input_fasta.keys()) == total_reads

    # Percent total read breakdown
    percent_reads_below_50 = (len(below_50_percent.keys())) / (len(input_fasta.keys())) * 100
    percent_reads_above_50 = (len(above_50_percent.keys())) / (len(input_fasta.keys())) * 100

    # Mean % N per-type
    mean_n_below_50 = (sum(below_50_percent.values()) / len(below_50_percent.values()))
    mean_n_above_50 = (sum(above_50_percent.values()) / len(above_50_percent.values()))

    # Avg. read length
    avg_repetitive_read_len = round((sum(repetitive_read_length_list) / len(repetitive_read_length_list)), 2)
    avg_nonrepetitive_read_len = round((sum(nonrepetitive_read_length_list) / len(nonrepetitive_read_length_list)), 2)
    avg_all_read_length = round((sum(all_read_lengths) / len(all_read_lengths)), 2)

    # Output nonrepetitive read name list
    output_nonrepetitive_read_names(outputdir, cross, nonrepetitive_read_name_list)
    
    # Print and log stats
    print("===============================")
    print("---- .masked file stats ----")
    print(f"Total number of reads in file: {total_reads:,}")
    print(f"Total repetitive reads in file: {(total_reads - total_nonrepetitive):,}")
    print(f"Total non-repetitive reads in file: {total_nonrepetitive:,}")
    print(f"Total avg. read length: {avg_all_read_length:,}")
    print(f"Avg. repetitive read length: {avg_repetitive_read_len:,}")
    print(f"Avg. non-repetitive read length: {avg_nonrepetitive_read_len:,}")
    print("")
    print("===============================")
    print("---- Reads Below 10kb in Length ----")
    print(f"% below 10kb: {round((below_ten_kb/total_reads)*100, 2)}%")
    print("")
    print("===============================")
    print("---- Reads Below 5kb in Length ----")
    print(f"% below 5kb: {round((below_five_kb/total_reads)*100, 2)}%")
    print("")
    print("===============================")
    print("---- Below 50% Repetitive ----")
    print(f"Number of Reads: {len(below_50_percent.keys()):,}")
    print(f"Percent of total reads: {round(percent_reads_below_50, 2)}%")
    print(f"Mean % Repetitive: {round(mean_n_below_50, 2)}")
    print("")
    print("===============================")
    print("---- Above 50% Repetitive ----")
    print(f"Number of Reads: {len(above_50_percent.keys()):,}")
    print(f"Percent of total reads: {round(percent_reads_above_50, 2)}%")
    print(f"Mean % Repetitive: {round(mean_n_above_50, 2)}")
    print("")
 

    logging.info("===============================")
    logging.info("---- .masked file stats ----")
    logging.info(f"Total number of reads in file: {total_reads:,}")
    logging.info(f"Total repetitive reads in file: {(total_reads - total_nonrepetitive):,}")
    logging.info(f"Total non-repetitive reads in file: {total_nonrepetitive:,}")
    logging.info(f"Total avg. read length: {avg_all_read_length:,}")
    logging.info(f"Avg. repetitive read length: {avg_repetitive_read_len:,}")
    logging.info(f"Avg. non-repetitive read length: {avg_nonrepetitive_read_len:,}")
    logging.info("")
    logging.info("===============================")
    logging.info("---- Reads Below 10kb in Length ----")
    logging.info(f"% below 10kb: {round((below_ten_kb/total_reads)*100, 2)}%")
    logging.info("")
    logging.info("===============================")
    logging.info("---- Reads Below 5kb in Length ----")
    logging.info(f"% below 5kb: {round((below_five_kb/total_reads)*100, 2)}%")
    logging.info("")
    logging.info("===============================")
    logging.info("---- 0-49% Repetitive ----")
    logging.info(f"Number of Reads: {len(below_50_percent.keys()):,}")
    logging.info(f"Percent of total reads: {round(percent_reads_below_50, 2)}%")
    logging.info(f"Mean % Repetitive: {round(mean_n_below_50, 2)}")
    logging.info("")
    logging.info("===============================")
    logging.info("---- 50-69% Repetitive ----")
    logging.info(f"Number of Reads: {len(above_50_percent.keys()):,}")
    logging.info(f"Percent of total reads: {round(percent_reads_above_50, 2)}%")
    logging.info(f"Mean % Repetitive: {round(mean_n_above_50, 2)}")
    logging.info("")
    return


if __name__ == "__main__":
    main()
