"""
Author: Andrew Harris
Python 3.8.3
"""
import argparse
import os
from pathlib import Path
import sys

from pyfaidx import Fasta
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly


# Set working directory to current directory
sys.path.insert(0, os.path.realpath(os.path.dirname(__file__)))
os.chdir(os.path.realpath(os.path.dirname(__file__)))


def phased_haplotype_analysis():
    """
    This script compares the TrioCanu output's from a reference cross 
    (bio-parent x bio-parent) to a replacement cross of one or both bio-parents.
    It provides a breakdown of the number of reads sorted to the same bins in both
    crosses (correctly sorted) and the number of reads that were sorted differently
    (incorrectly sorted) from the reference cross. It also further breaks down
    the incorrectly sorted reads into bins describing where the read was
    incorrectly sorted to in the replacement cross.
    """
    # ---------------- Argparse inputs ----------------
    parser = argparse.ArgumentParser(description="Phased Haplotype Analysis (PHA)")
    # Required input arguments
    parser.add_argument('--bio_mom', type=str, required=True,
                        action='store', help='Pathway to biological mothers fasta file')
    parser.add_argument('--bio_dad', type=str, required=True,
                        action='store', help='Pathway to biological fathers fasta file')
    parser.add_argument('--bio_unknown', type=str, required=True,
                        action='store', help='Pathway to biological unknown fasta file')
    parser.add_argument('--nonbio_mom', type=str, required=True,
                        action='store', help='Pathway to non-biological mothers fasta file')
    parser.add_argument('--nonbio_dad', type=str, required=True,
                        action='store', help='Pathway to non-biological mothers fasta file')
    parser.add_argument('--nonbio_unknown', type=str, required=True,
                        action='store', help='Pathway to non-biological unknown fasta file')
    parser.add_argument('--cross', type=str, required=True, action='store',
                        help="Name of replacement cross (i.e. LilBubxPbe14).")
    parser.add_argument('--biocross', type=str, required=True, action='store',
                        help="Name of replacement cross (i.e. LilBubxPbe14).")
    parser.add_argument('--excel', type=str, required=True, action='store',
                        help="Excel file with TrioCanu output read counts and coverages.")
    # Optional - defaults to False
    parser.add_argument('--show_graph', action='store_True',
                        default=False, help='Automatically load graph into browser')
    args = parser.parse_args()

    # input variables
    bio_mom = args.bio_mom
    bio_dad = args.bio_dad
    bio_unknown = args.bio_unknown
    non_bio_mom = args.nonbio_mom
    non_bio_dad = args.nonbio_dad
    non_bio_unknown = args.nonbio_unknown
    cross_used = args.cross
    bio_parents = args.biocross
    auto_open = args.show_graph
    excel_file_pathway = args.excel

    # ---------------- Helper Functions ----------------
    def output_directory_setup():

        # Create cross used output directory
        Path(f"./read_names/{cross_used}").mkdir(parents=True, exist_ok=True)

        # Create data output directory
        isr_pathway = f"./data_output/{cross_used}/incorrectly_sorted_reads/"
        Path(isr_pathway).mkdir(parents=True, exist_ok=True)

        return

    def self_check(df, individual):
        """
        This function will perform a duplicate drop test on each file to see if a particular individual contains duplicate reads.
        """
        df_copy = df.copy()
        df_copy.sort_values(by='names', inplace=True)
        before_drop_length = len(df_copy)

        df_copy.drop_duplicates(subset='names', keep=False, inplace=True)

        if before_drop_length - len(df_copy) == 0:
            print("--------------------------------------------")
            print("\t-- Duplicate Self Check --")
            print(
                "|-- {} input file **DOES NOT** contain any duplicate reads".format(individual))
            return True
        else:
            print(
                "|-- {} input file **DOES** contain any duplicate reads".format(individual))
            return False

    def read_names_to_csv(file_pathway, sample_name):
        """
        This function takes in the fasta files and pulls the read names and then outputs them into csv files.
        :param file_pathway: Pathway to Fasta files
        :param sample_name: Name of current sample being loaded (i.e. bio_mom)
        :return: The output filename
        """
        output_filename = f"./read_names/{cross_used}/{sample_name}.csv"
        output_df = pd.DataFrame()
        print(f"-- Loading {sample_name} Fasta file --")
        file_read = Fasta(file_pathway)
        print(f"-- Loaded {sample_name} Fasta file --")
        file_read_names_list = [name for name in file_read.keys()]
        file_read_names_list.sort()
        output_df['names'] = file_read_names_list
        print(f"-- Output {sample_name} read names to {sample_name}.csv --")
        output_df.to_csv(output_filename, sep=',', index=False)
        print(f'-- {sample_name} finished --')
        return output_filename

    def mom_check(bio_mom, nonbio_mom):
        """
        This function takes in the bio and non-bio mom reads and...
        1. Merges the two list together
        2. Drops any duplicate (correctly sorted reads)
        3. Separates the remaining reads back into bio and non-bio lists.

        The remaining bio-mom reads are the reads lost in the sorting.
        The remaining nonbio-mom reads are the reads that were gained in the sorting.
        """
        # Run self check to see if any individuals contain duplicate reads
        self_check(bio_mom, 'bio_mom')
        self_check(nonbio_mom, 'nonbio_mom')

        concat_df = pd.concat([bio_mom, nonbio_mom], sort=True)
        concat_df = concat_df.reset_index(drop=True)
        concat_df.sort_values(by='names', inplace=True)
        before_length = len(concat_df)

        print("--------------------------------------------")
        print('|--')
        print("|-- Maternal read count before drop:", len(concat_df))

        concat_df.drop_duplicates(subset='names', keep=False, inplace=True)

        print("|-- Maternal read count after drop:", len(concat_df))
        print("|-- Maternal read count difference:",
              (before_length - len(concat_df)))
        print('|--')

        num_correctly_sorted = (before_length - len(concat_df)) / 2

        nonbio_mom_tally = []
        [nonbio_mom_tally.append(
            1) for relationship in concat_df["parental_relationships"] if relationship == 'nonbio_mom']

        bio_mom_tally = []
        [bio_mom_tally.append(
            1) for relationship in concat_df["parental_relationships"] if relationship == 'bio_mom']

        lost_mom_reads = concat_df[concat_df.isin(
            {'bio_mom'}).any(1)]  # removed ['names']
        gained_mom_reads = concat_df[concat_df.isin(
            {'nonbio_mom'}).any(1)]  # removed ['names']

        lost_mom_reads.reset_index(inplace=True, drop=True)
        gained_mom_reads.reset_index(inplace=True, drop=True)

        print("--------------------------------------------")
        print('|--')
        print('|-- Total maternal reads correctly sorted:', num_correctly_sorted)
        print('|-- Total maternal reads gained:', sum(nonbio_mom_tally))
        print('|-- Total maternal reads lost:', sum(bio_mom_tally))
        print('|-- Total gained+lost:',
              sum(bio_mom_tally) + sum(nonbio_mom_tally))
        print('|--')
        print("--------------------------------------------")
        print("------ Maternal Read Sorting Complete ------")
        print("--------------------------------------------")
        print()
        print("\t-- Moving to Paternal Reads --\t")

        return lost_mom_reads, gained_mom_reads, num_correctly_sorted

    def dad_check(bio_dad, nonbio_dad):
        """
        This function takes in the bio and non-bio dad reads and...
        1. Merges the two list together
        2. Drops any duplicate (correctly sorted reads)
        3. Separates the remaining reads back into bio and non-bio lists.

        The remaining bio-dad reads are the reads lost in the sorting.
        The remaining nonbio-dad reads are the reads that were gained in the sorting.
        """

        # Run self check to see if any individuals contain duplicate reads
        self_check(bio_dad, 'bio_dad')
        self_check(nonbio_dad, 'nonbio_dad')

        concat_df = pd.concat([bio_dad, nonbio_dad], sort=True)
        concat_df = concat_df.reset_index(drop=True)
        concat_df.sort_values(by='names', inplace=True)
        before_length = len(concat_df)

        print("--------------------------------------------")
        print('|--')
        print("|-- Paternal before drop count:", len(concat_df))

        concat_df.drop_duplicates(subset='names', keep=False, inplace=True)

        print("|-- Paternal after drop count:", len(concat_df))
        print("|-- Paternal read count difference:",
              (before_length - len(concat_df)))
        print('|--')

        # Divide by two because we dropped two copies of each read
        num_correctly_sorted = (before_length - len(concat_df)) / 2

        nonbio_dad_tally = []
        [nonbio_dad_tally.append(
            1) for relationship in concat_df["parental_relationships"] if relationship == 'nonbio_dad']

        bio_dad_tally = []
        [bio_dad_tally.append(
            1) for relationship in concat_df["parental_relationships"] if relationship == 'bio_dad']

        lost_dad_reads = concat_df[concat_df.isin({'bio_dad'}).any(1)]
        gained_dad_reads = concat_df[concat_df.isin({'nonbio_dad'}).any(1)]

        lost_dad_reads.reset_index(inplace=True, drop=True)
        gained_dad_reads.reset_index(inplace=True, drop=True)

        print("--------------------------------------------")
        print('|--')
        print('|-- Total paternal reads correctly sorted:', num_correctly_sorted)
        print('|-- Total paternal reads gained:', sum(nonbio_dad_tally))
        print('|-- Total paternal reads lost:', sum(bio_dad_tally))
        print('|-- Total gained+lost:',
              sum(bio_dad_tally) + sum(nonbio_dad_tally))
        print('|--')
        print("--------------------------------------------")
        print("------ Paternal Read Sorting Complete ------")
        print("--------------------------------------------")
        print()
        print("\t-- Moving to Unknown Reads --\t")

        return lost_dad_reads, gained_dad_reads, num_correctly_sorted

    def unknown_check(bio_unknown, nonbio_unknown):
        """
        This function takes in the bio and non-bio unknown reads and...
        1. Merges the two list together
        2. Drops any duplicate (correctly sorted reads)
        3. Separates the remaining reads back into bio and non-bio lists.

        The remaining bio-unknown reads are the reads lost in the sorting.
        The remaining nonbio-unknown reads are the reads that were gained in the sorting.
        """

        # Run self check to see if any individuals contain duplicate reads
        self_check(bio_unknown, 'bio_unknown')
        self_check(nonbio_unknown, 'nonbio_unknown')

        concat_df = pd.concat([bio_unknown, nonbio_unknown], sort=True)
        concat_df = concat_df.reset_index(drop=True)
        concat_df.sort_values(by='names', inplace=True)
        before_length = len(concat_df)

        print("--------------------------------------------")
        print('|--')
        print("|-- Unknown before drop count:", len(concat_df))

        concat_df.drop_duplicates(subset='names', keep=False, inplace=True)
        print("|-- Unknown after drop count:", len(concat_df))
        print("|-- Unknown read count difference:",
              (before_length - len(concat_df)))
        print('|--')

        num_correctly_sorted = (before_length - len(concat_df)) / 2

        nonbio_unknown_tally = []
        [nonbio_unknown_tally.append(
            1) for relationship in concat_df["parental_relationships"] if relationship == 'nonbio_unknown']

        bio_unknown_tally = []
        [bio_unknown_tally.append(
            1) for relationship in concat_df["parental_relationships"] if relationship == 'bio_unknown']

        lost_unknown_reads = concat_df[concat_df.isin({'bio_unknown'}).any(1)]

        gained_unknown_reads = concat_df[concat_df.isin(
            {'nonbio_unknown'}).any(1)]

        lost_unknown_reads.reset_index(inplace=True, drop=True)
        gained_unknown_reads.reset_index(inplace=True, drop=True)

        print("--------------------------------------------")
        print('|--')
        print('|-- Total unknown reads correctly sorted:', num_correctly_sorted)
        print('|-- Total unknown reads gained:', sum(nonbio_unknown_tally))
        print('|-- Total unknown reads lost:', sum(bio_unknown_tally))
        print('|-- Total gained+lost:',
              sum(bio_unknown_tally) + sum(nonbio_unknown_tally))
        print('|--')
        print("--------------------------------------------")
        print("------ Unknown Read Sorting Complete ------")
        print("--------------------------------------------")
        print()
        print("\t*** -- Initiating Phase 2 -- ***\t\n")

        return lost_unknown_reads, gained_unknown_reads, num_correctly_sorted

    def maternal_load_n_run(bio_mom, non_bio_mom):
        bio_mom_read_to_csv = read_names_to_csv(bio_mom, 'bio_mom')
        bio_mom_read_names = pd.read_csv(
            bio_mom_read_to_csv, dtype={'names': str})
        bm_df = pd.DataFrame(
            [value for value in bio_mom_read_names['names']], columns=['names'])
        bm_df['parental_relationships'] = [
            'bio_mom' for _ in range(len(bm_df.index))]
        bio_mom_read_count = len(bm_df)
        stats_collecting_df.at["biomom_initial_read_count", str(
            cross_used)] = bio_mom_read_count

        non_bio_mom_read_to_csv = read_names_to_csv(non_bio_mom, 'nonbio_mom')
        non_bio_mom_read_names = pd.read_csv(
            non_bio_mom_read_to_csv, dtype={'names': str})
        nbm_df = pd.DataFrame(
            [value for value in non_bio_mom_read_names['names']], columns=['names'])
        nbm_df['parental_relationships'] = [
            'nonbio_mom' for _ in range(len(nbm_df.index))]

        print("--------------------------------------------")
        print("-- Maternal gained/lost read information --")
        print("--------------------------------------------")
        print('|--')
        print('|-- Number of bio mom reads:', len(bm_df))
        print('|-- Number of nonbio mom reads:', len(nbm_df))
        print('|--')

        mom_check_output = mom_check(bm_df, nbm_df)

        return mom_check_output, bio_mom_read_count

    def paternal_load_n_run(bio_dad, non_bio_dad):
        bio_dad_read_to_csv = read_names_to_csv(bio_dad, 'bio_dad')
        bio_dad_read_names = pd.read_csv(
            bio_dad_read_to_csv, dtype={'names': str})
        bd_df = pd.DataFrame(
            [value for value in bio_dad_read_names['names']], columns=['names'])
        bd_df['parental_relationships'] = [
            'bio_dad' for _ in range(len(bd_df.index))]
        bio_dad_read_count = len(bd_df)
        stats_collecting_df.at["biodad_initial_read_count", str(
            cross_used)] = bio_dad_read_count

        non_bio_dad_read_to_csv = read_names_to_csv(non_bio_dad, 'nonbio_dad')
        non_bio_dad_read_names = pd.read_csv(
            non_bio_dad_read_to_csv, dtype={'names': str})
        nbd_df = pd.DataFrame(
            [value for value in non_bio_dad_read_names['names']], columns=['names'])
        nbd_df['parental_relationships'] = [
            'nonbio_dad' for _ in range(len(nbd_df.index))]
        print()
        print("--------------------------------------------")
        print("-- Paternal gained/lost read information --")
        print("--------------------------------------------")
        print('|--')
        print('|-- Number of bio dad reads:', len(bd_df))
        print('|-- Number of nonbio dad reads:', len(nbd_df))
        print('|--')

        dad_check_output = dad_check(bd_df, nbd_df)

        return dad_check_output, bio_dad_read_count

    def unknown_load_n_run(bio_unknown, non_bio_unknown):
        # Load Bio Unknown and py
        bio_unknown_read_to_csv = read_names_to_csv(bio_unknown, 'bio_unknown')
        bio_unknown_read_names = pd.read_csv(
            bio_unknown_read_to_csv, dtype={'names': str})
        bu_df = pd.DataFrame(
            [value for value in bio_unknown_read_names['names']], columns=['names'])
        bu_df['parental_relationships'] = [
            'bio_unknown' for _ in range(len(bu_df.index))]
        bio_unknown_read_count = len(bu_df)
        stats_collecting_df.at["biounknown_initial_read_count", str(
            cross_used)] = bio_unknown_read_count

        non_bio_unknown_read_to_csv = read_names_to_csv(
            non_bio_unknown, 'nonbio_unknown')
        non_bio_unknown_read_names = pd.read_csv(
            non_bio_unknown_read_to_csv, dtype={'names': str})
        nbu_df = pd.DataFrame(
            [value for value in non_bio_unknown_read_names['names']], columns=['names'])
        nbu_df['parental_relationships'] = [
            'nonbio_unknown' for _ in range(len(nbu_df.index))]
        print()
        print("--------------------------------------------")
        print("----Unknown gained/lost read information----")
        print("--------------------------------------------")
        print('|--')
        print('|-- Number of bio unknown reads:', len(bu_df))
        print('|-- Number of nonbio unknown reads:', len(nbu_df))

        unknown_check_output = unknown_check(bu_df, nbu_df)

        return unknown_check_output, bio_unknown_read_count

    def phase_2(maternal_lists, paternal_lists, unknown_lists):
        """Indentify number and direction of incorrect sorting"""
        # Lists are ordered [lost reads (bio), gained reads (nonbio)]
        maternal_lost_reads = maternal_lists[0]
        maternal_gained_reads = maternal_lists[1]

        paternal_lost_reads = paternal_lists[0]
        paternal_gained_reads = paternal_lists[1]

        unknown_lost_reads = unknown_lists[0]
        unknown_gained_reads = unknown_lists[1]

        def maternal_to_other_stats(maternal_lost, paternal_gained, unknown_gained):
            def paternal_gain_check(maternal_lost, paternal_gain_list):
                concat = pd.concat([maternal_lost, paternal_gain_list])

                before_length = len(concat)

                before_drop_df = concat.copy()
                before_drop_df.drop_duplicates(
                    subset='names', keep='first', inplace=True)

                concat.drop_duplicates(
                    subset='names', keep=False, inplace=True)

                difference = before_length - len(concat)
                number_of_reads = int(difference / 2)

                duplicate_df = pd.concat([concat, before_drop_df])
                duplicate_df.sort_values(by='names', inplace=True)

                duplicate_df.drop_duplicates(
                    subset='names', keep=False, inplace=True)
                duplicate_df['names'].to_csv(
                    "./data_output/{}/incorrectly_sorted_reads/mom_to_dad_reads.csv".format(cross_used), sep=',', index=False)

                return number_of_reads

            def unknown_gain_check(maternal_lost, unknown_gain_list):
                concat = pd.concat([maternal_lost, unknown_gain_list])

                before_length = len(concat)

                before_drop_df = concat.copy()
                before_drop_df.drop_duplicates(
                    subset='names', keep='first', inplace=True)

                concat.drop_duplicates(
                    subset='names', keep=False, inplace=True)

                difference = before_length - len(concat)
                number_of_reads = int(difference / 2)

                duplicate_df = pd.concat([concat, before_drop_df])
                duplicate_df.sort_values(by='names', inplace=True)

                duplicate_df.drop_duplicates(
                    subset='names', keep=False, inplace=True)
                duplicate_df['names'].to_csv(
                    "./data_output/{}/incorrectly_sorted_reads/mom_to_unknown_reads.csv".format(cross_used), sep=',', index=False)

                return number_of_reads

            count_to_dad = paternal_gain_check(maternal_lost, paternal_gained)

            count_to_unknown = unknown_gain_check(
                maternal_lost, unknown_gained)

            count = count_to_dad + count_to_unknown

            print("--------------------------------------------")
            print('\t-- Maternal Check Complete --')
            print('|-- Num. Reads 2 Dad = {}'. format(count_to_dad))
            print('|-- Num. Reads 2 Unknown = {}'. format(count_to_unknown))
            print('|-- Total Moved Reads = {}'. format(count))
            print("--------------------------------------------")

            return count_to_dad, count_to_unknown

        def paternal_to_other_stats(paternal_lost, maternal_gained, unknown_gained):
            def maternal_gain_check(paternal_lost, maternal_gained):
                concat = pd.concat([paternal_lost, maternal_gained])
                before_length = len(concat)

                before_drop_df = concat.copy()
                before_drop_df.drop_duplicates(
                    subset='names', keep='first', inplace=True)

                concat.drop_duplicates(
                    subset='names', keep=False, inplace=True)

                difference = before_length - len(concat)
                number_of_reads = int(difference / 2)

                duplicate_df = pd.concat([concat, before_drop_df])
                duplicate_df.sort_values(by='names', inplace=True)

                duplicate_df.drop_duplicates(
                    subset='names', keep=False, inplace=True)
                duplicate_df['names'].to_csv(
                    "./data_output/{}/incorrectly_sorted_reads/dad_to_mom_reads.csv".format(cross_used), sep=',', index=False)

                return number_of_reads

            def unknown_gain_check(paternal_lost, unknown_gained):
                concat = pd.concat([paternal_lost, unknown_gained])
                before_length = len(concat)

                before_drop_df = concat.copy()
                before_drop_df.drop_duplicates(
                    subset='names', keep='first', inplace=True)

                concat.drop_duplicates(
                    subset='names', keep=False, inplace=True)

                difference = before_length - len(concat)
                number_of_reads = int(difference / 2)

                duplicate_df = pd.concat([concat, before_drop_df])
                duplicate_df.sort_values(by='names', inplace=True)

                duplicate_df.drop_duplicates(
                    subset='names', keep=False, inplace=True)
                duplicate_df['names'].to_csv(
                    "./data_output/{}/incorrectly_sorted_reads/dad_to_unknown_reads.csv".format(cross_used), sep=',', index=False)

                return number_of_reads

            count_to_mom = maternal_gain_check(paternal_lost, maternal_gained)

            count_to_unknown = unknown_gain_check(
                paternal_lost, unknown_gained)

            count = count_to_mom + count_to_unknown

            print("--------------------------------------------")
            print('\t-- Paternal Check Complete --')
            print('|-- Num. Reads 2 Mom = {}'. format(count_to_mom))
            print('|-- Num. Reads 2 Unknown = {}'. format(count_to_unknown))
            print('|-- Total Moved Reads = {}'. format(count))
            print("--------------------------------------------")

            return count_to_mom, count_to_unknown

        def unknown_to_other_stats(unknown_lost, maternal_gained, paternal_gained):
            def mom_gain_check(unknown_lost, maternal_gained):
                concat = pd.concat([unknown_lost, maternal_gained])
                before_length = len(concat)

                before_drop_df = concat.copy()
                before_drop_df.drop_duplicates(
                    subset='names', keep='first', inplace=True)

                concat.drop_duplicates(
                    subset='names', keep=False, inplace=True)

                difference = before_length - len(concat)
                number_of_reads = int(difference / 2)

                duplicate_df = pd.concat([concat, before_drop_df])
                duplicate_df.sort_values(by='names', inplace=True)

                duplicate_df.drop_duplicates(
                    subset='names', keep=False, inplace=True)
                duplicate_df['names'].to_csv(
                    "./data_output/{}/incorrectly_sorted_reads/unknown_to_mom_reads.csv".format(cross_used), sep=',', index=False)

                return number_of_reads

            def dad_gain_check(unknown_lost, paternal_gained):
                concat = pd.concat([unknown_lost, paternal_gained])
                before_length = len(concat)

                before_drop_df = concat.copy()
                before_drop_df.drop_duplicates(
                    subset='names', keep='first', inplace=True)

                concat.drop_duplicates(
                    subset='names', keep=False, inplace=True)

                difference = before_length - len(concat)
                number_of_reads = int(difference / 2)

                duplicate_df = pd.concat([concat, before_drop_df])
                duplicate_df.sort_values(by='names', inplace=True)

                duplicate_df.drop_duplicates(
                    subset='names', keep=False, inplace=True)
                duplicate_df['names'].to_csv(
                    "./data_output/{}/incorrectly_sorted_reads/unknown_to_dad_reads.csv".format(cross_used), sep=',', index=False)

                return number_of_reads

            count_to_mom = mom_gain_check(unknown_lost, maternal_gained)

            count_to_dad = dad_gain_check(unknown_lost, paternal_gained)

            count = count_to_dad + count_to_mom

            print("--------------------------------------------")
            print('\t-- Unknown Check Complete --')
            print('|-- Num. Reads 2 Mom = {}'. format(count_to_mom))
            print('|-- Num. Reads 2 Dad = {}'. format(count_to_dad))
            print('|-- Total Moved Reads = {}'. format(count))
            print("--------------------------------------------")

            return count_to_mom, count_to_dad

        # Stats function calls
        maternal_stats = maternal_to_other_stats(
            maternal_lost_reads,
            paternal_gained_reads,
            unknown_gained_reads,
        )

        paternal_stats = paternal_to_other_stats(
            paternal_lost_reads,
            maternal_gained_reads,
            unknown_gained_reads,
        )

        unknown_stats = unknown_to_other_stats(
            unknown_lost_reads,
            maternal_gained_reads,
            paternal_gained_reads,
        )

        """
        Return will look like: ((Maternal Stats), (Paternal Stats), (Unknown Stats))
        Maternal Stats = ((Num. reads to Dad), (Num. reads to Unknown))
        Paternal Stats = ((Num. reads to Mom), (Num. reads to Unknown))
        Unknown Stats = ((Num. reads to Mom), (Num. reads to Dad))
        """
        return maternal_stats, paternal_stats, unknown_stats

    def graph_results(
        parental_stats, 
        mom_correct_read_count, 
        dad_correct_read_count, 
        unknown_correct_read_count, 
        ef, 
        cu,
    ):
        """
        Input data reference:
        ((250, 0), (0, 100), (0, 0))
        ((mom2dad=[0][0], mom2unknown=[0][1]), (dad2mom=[1][0], dad2unknown=[1][1]), (unknown2mom=[2][0], unknown2dad=[2][1]))
        """
        # Read count graph data setup
        mom_to_dad = parental_stats[0][0]
        mom_to_unknown = parental_stats[0][1]

        dad_to_mom = parental_stats[1][0]
        dad_to_unknown = parental_stats[1][1]

        unknown_to_mom = parental_stats[2][0]
        unknown_to_dad = parental_stats[2][1]

        stats_collecting_df.at["mom_to_dad_count",
                               str(cross_used)] = mom_to_dad
        stats_collecting_df.at["mom_to_unknown_count",
                               str(cross_used)] = mom_to_unknown

        stats_collecting_df.at["dad_to_mom_count",
                               str(cross_used)] = dad_to_mom
        stats_collecting_df.at["dad_to_unknown_count",
                               str(cross_used)] = dad_to_unknown

        stats_collecting_df.at["unknown_to_mom_count",
                               str(cross_used)] = unknown_to_mom
        stats_collecting_df.at["unknown_to_dad_count",
                               str(cross_used)] = unknown_to_dad

        read_count_coverage_y_max = max([mom_to_dad, mom_to_unknown, dad_to_mom, dad_to_unknown, unknown_to_dad, unknown_to_mom, mom_correct_read_count,
                                         dad_correct_read_count, unknown_correct_read_count, bio_mom_read_count, bio_dad_read_count, bio_unknown_read_count])

        # Coverage graph data setup

        # Load in excel file
        read_excel = pd.read_excel(ef, index_col='Cross')

        # clean dataset
        read_excel.drop(columns=['avg. len of genome', 'Maternal Base Count',
                                 'Paternal Base Count', 'Unknown Base Count'], inplace=True)

        # Round values to 2 decimals
        read_excel = read_excel.round(
            {'Maternal Coverage': 0, 'Paternal Coverage': 0, 'Unknown Coverage': 3})

        # Set y-max for coverage graphs
        coverage_y_max = []
        [coverage_y_max.append(value)
         for value in read_excel['Maternal Coverage']]
        [coverage_y_max.append(value)
         for value in read_excel['Paternal Coverage']]
        [coverage_y_max.append(value)
         for value in read_excel['Unknown Coverage']]

        fig = make_subplots(
            rows=2,
            cols=3,
            # subplot_titles=('', 'Read Count Analysis: {}'.format(cross_used), '', '', 'Coverage Analysis: {}'.format(cross_used), '')
        )
        # Read Count traces
        # Maternal Read Counts
        fig.add_trace(go.Bar(x=['Maternal Read Count'], y=[bio_mom_read_count], name='Bio Maternal Count', text=str(
            bio_mom_read_count), legendgroup="bio", marker=dict(color=['blue'])), row=1, col=1)
        fig.add_trace(go.Bar(x=['Maternal Read Count'], y=[mom_correct_read_count], name='Correctly Sorted Mom', text=str(
            mom_correct_read_count), legendgroup="correctly_sorted", marker=dict(color=['purple'])), row=1, col=1)
        fig.add_trace(go.Bar(x=['Maternal Read Count'], y=[mom_to_dad], name='mom -> dad', text=str(
            mom_to_dad), legendgroup="maternal", marker=dict(color=['pink'])), row=1, col=1)
        fig.add_trace(go.Bar(x=['Maternal Read Count'], y=[mom_to_unknown], name='mom -> unknown', text=str(
            mom_to_unknown), legendgroup="maternal", marker=dict(color=['pink'])), row=1, col=1)

        # Paternal Read Counts
        fig.add_trace(go.Bar(x=['Paternal Read Count'], y=[bio_dad_read_count], name='Bio Paternal Count', text=str(
            bio_dad_read_count), legendgroup="bio", marker=dict(color=['blue'])), row=1, col=2)
        fig.add_trace(go.Bar(x=['Paternal Read Count'], y=[dad_correct_read_count], name='Correctly Sorted Dad', text=str(
            dad_correct_read_count), legendgroup="correctly_sorted", marker=dict(color=['purple'])), row=1, col=2)
        fig.add_trace(go.Bar(x=['Paternal Read Count'], y=[dad_to_mom], name='dad -> mom', text=str(
            dad_to_mom), legendgroup="paternal", marker=dict(color=['green'])), row=1, col=2)
        fig.add_trace(go.Bar(x=['Paternal Read Count'], y=[dad_to_unknown], name='dad -> unknown', text=str(
            dad_to_unknown), legendgroup="paternal", marker=dict(color=['green'])), row=1, col=2)

        # Unknown Read Counts
        fig.add_trace(go.Bar(x=['Unknown Read Count'], y=[bio_unknown_read_count], name='Bio Unknown Count', text=str(
            bio_unknown_read_count), legendgroup="bio", marker=dict(color=['blue'])), row=1, col=3)
        fig.add_trace(go.Bar(x=['Unknown Read Count'], y=[unknown_correct_read_count], name='Correctly Sorted Unknown', text=str(
            unknown_correct_read_count), legendgroup="correctly_sorted", marker=dict(color=['purple'])), row=1, col=3)
        fig.add_trace(go.Bar(x=['Unknown Read Count'], y=[unknown_to_mom], name='unknown -> mom', text=str(
            unknown_to_mom), legendgroup="unknown", marker=dict(color=['grey'])), row=1, col=3)
        fig.add_trace(go.Bar(x=['Unknown Read Count'], y=[unknown_to_dad], name='unknown -> dad', text=str(
            unknown_to_dad), legendgroup="unknown", marker=dict(color=['grey'])), row=1, col=3)

        # Coverage Traces

        # Bio-Parents
        fig.add_trace(go.Bar(x=['Maternal Coverage'],
                             y=[read_excel.at[bio_parents, 'Maternal Coverage']],
                             name='Bio Maternal Coverage',
                             text='{}'.format(
                                 read_excel.at[bio_parents, 'Maternal Coverage']),
                             legendgroup="Bio-Parents", marker=dict(color=['red'])), row=2, col=1
                      )

        fig.add_trace(go.Bar(x=['Paternal Coverage'],
                             y=[read_excel.at[bio_parents, 'Paternal Coverage']],
                             name='Bio Paternal Coverage',
                             text='{}'.format(
                                 read_excel.at[bio_parents, 'Paternal Coverage']),
                             legendgroup="Bio-Parents", marker=dict(color=['red'])), row=2, col=2
                      )

        fig.add_trace(go.Bar(x=['Unknown Coverage'],
                             y=[read_excel.at[bio_parents, 'Unknown Coverage']],
                             name='Bio Unknown Coverage',
                             text='{}'.format(
                                 read_excel.at[bio_parents, 'Unknown Coverage']),
                             legendgroup="Bio-Parents", marker=dict(color=['red'])), row=2, col=3
                      )

        # Non-bio Crosses
        fig.add_trace(go.Bar(x=['Maternal Coverage'], y=[read_excel.at["{}".format(cu), 'Maternal Coverage']], name='{} Maternal Coverage'.format(
            cu), text='{}'.format(read_excel.at["{}".format(cu), 'Maternal Coverage']), legendgroup="{}".format(cu), marker=dict(color=['orange'])), row=2, col=1)
        fig.add_trace(go.Bar(x=['Paternal Coverage'], y=[read_excel.at["{}".format(cu), 'Paternal Coverage']], name='{} Paternal Coverage'.format(
            cu), text='{}'.format(read_excel.at["{}".format(cu), 'Paternal Coverage']), legendgroup="{}".format(cu), marker=dict(color=['orange'])), row=2, col=2)
        fig.add_trace(go.Bar(x=['Unknown Coverage'], y=[read_excel.at["{}".format(cu), 'Unknown Coverage']], name='{} Unknown Coverage'.format(
            cu), text='{}'.format(read_excel.at["{}".format(cu), 'Unknown Coverage']), legendgroup="{}".format(cu), marker=dict(color=['orange'])), row=2, col=3)

        # Update Count Graph Properties
        fig.update_layout(
            title_text='TrioCanu Phased Haplotype Analysis: {}'.format(
                cross_used),
            template="plotly_dark"
        )

        fig.update_traces(
            texttemplate='%{text:.4s}', textposition='outside', row=1)

        fig.update_yaxes(title_text='Number of Reads', title_font=dict(size=18), range=[
                         0, (read_count_coverage_y_max + (read_count_coverage_y_max * .1))], row=1, col=1)
        fig.update_yaxes(title_text='Number of Reads', title_font=dict(size=18), range=[
                         0, (read_count_coverage_y_max + (read_count_coverage_y_max * .1))], row=1, col=2)
        fig.update_yaxes(title_text='Number of Reads', title_font=dict(
            size=18), range=[0, 20000], row=1, col=3)

        fig.update_traces(
            texttemplate='%{text}X', textposition='outside', row=2)

        # fig.update_yaxes(title_text='Coverage (X)', title_font=dict(size=18), range=[0, (coverage_yaxis_ma + (coverage_yaxis_ma * .1))], row=2, col=1)
        # fig.update_yaxes(title_text='Coverage (X)', title_font=dict(size=18), range=[0, (coverage_yaxis_ma + (coverage_yaxis_ma * .1))], row=2, col=2)
        fig.update_yaxes(title_text='Coverage (X)', title_font=dict(
            size=18), range=[0, 50], row=2, col=1)
        fig.update_yaxes(title_text='Coverage (X)', title_font=dict(
            size=18), range=[0, 50], row=2, col=2)
        fig.update_yaxes(title_text='Coverage (X)', title_font=dict(
            size=18), range=[0, 0.05], row=2, col=3)

        fig.update_xaxes(title_font=dict(size=20), row=2, col=1)
        fig.update_xaxes(title_font=dict(size=20), row=2, col=2)
        fig.update_xaxes(title_font=dict(size=20), row=2, col=3)

        if auto_open:
            plotly.offline.plot(fig, filename="./data_output/{}/nonbio_parents_{}.html".format(
                cross_used, cross_used), auto_open=True)
        else:
            plotly.offline.plot(fig, filename="./data_output/{}/nonbio_parents_{}.html".format(
                cross_used, cross_used), auto_open=False)
        return

    # ---------------- Body ----------------

    # Set up out directories
    output_directory_setup()

    # Set up stats collection DataFrame
    stats_collecting_df = pd.DataFrame(
        index=[
            "biomom_initial_read_count",
            "biodad_initial_read_count",
            "biounknown_initial_read_count",
            "correctly_sorted_to_mom",
            "mom_to_dad_count",
            "mom_to_unknown_count",
            "correctly_sorted_to_dad",
            "dad_to_mom_count",
            "dad_to_unknown_count",
            "correctly_sorted_to_unknown",
            "unknown_to_mom_count",
            "unknown_to_dad_count",
        ],
        columns=[str(cross_used)]
    )

    # Print Phase 1 statement
    print()
    print(f"\t*** -- Cross: {cross_used} -- ***\t")
    print("\t*** -- Starting Phase 1 -- ***\t")
    print()

    # Run read name comparisons
    maternal_gain_loss_list = maternal_load_n_run(bio_mom, non_bio_mom)
    paternal_gain_loss_list = paternal_load_n_run(bio_dad, non_bio_dad)
    unknown_gain_loss_list = unknown_load_n_run(bio_unknown, non_bio_unknown)

    bio_mom_read_count = maternal_gain_loss_list[1]
    bio_dad_read_count = paternal_gain_loss_list[1]
    bio_unknown_read_count = unknown_gain_loss_list[1]

    stats_collecting_df.at["correctly_sorted_to_mom",
                           str(cross_used)] = maternal_gain_loss_list[0][2]
    stats_collecting_df.at["correctly_sorted_to_dad",
                           str(cross_used)] = paternal_gain_loss_list[0][2]
    stats_collecting_df.at["correctly_sorted_to_unknown",
                           str(cross_used)] = unknown_gain_loss_list[0][2]

    # Run Phase 2
    stats = phase_2(
        maternal_gain_loss_list[0][:2], paternal_gain_loss_list[0][:2], unknown_gain_loss_list[0][:2])

    graph_results(stats, maternal_gain_loss_list[0][2], paternal_gain_loss_list[0]
                  [2], unknown_gain_loss_list[0][2], excel_file_pathway, cross_used)

    stats_collecting_df.to_csv(
        "./data_output/{}/{}_stats.csv".format(cross_used, cross_used), sep=',')
    return


if __name__ == "__main__":
    phased_haplotype_analysis()
