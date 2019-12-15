# -*- coding: utf-8 -*-
"""Main module."""

import sys
import pickle
import click
import logging
from send2trash import send2trash
from binomial_p import compute_binom_p
from crispector_algorithm import CrispectorAlgorithm
from exceptions import FastpRunTimeError, NoneValuesInAmpliconsCSV, SgRNANotInReferenceSequence, \
    ConfiguratorIsCalledBeforeInitConfigPath, PriorPositionHasWrongLength, \
    UnknownAlignmentChar, Bowtie2RunTimeError, Bowtie2BuildRunTimeError, CantOpenDemultiplexedSamFile, \
    AlignerSubstitutionDoesntExist, ClassificationFailed, BadSgRNAChar, BadReferenceAmpliconChar
from constants_and_types import ExpType, Path, FASTP_DIR, welcome_msg, FREQ, TX_READ_NUM, MOCK_READ_NUM, EDIT_PERCENT, \
    SUMMARY_RESULTS_TITLES, SITE_NAME, REFERENCE, SGRNA, ON_TARGET, CUT_SITE, AlgResult, R_PRIMER, F_PRIMER
from input_processing import InputProcessing
import traceback
from utils import Logger, Configurator, plot_editing_activity, create_reads_statistics_report, summary_result_to_excel, \
    discarded_sites_text
import os
import pandas as pd
from modification_tables import ModificationTables
from typing import Dict
from modification_types import ModificationTypes

# TODO - coin. when not tx == mock number of reads
def run(tx_in1: Path, tx_in2: Path, mock_in1: Path, mock_in2: Path, output: Path, amplicons_csv: Path,
        fastp_options_string: str, override_fastp: bool, keep_fastp_output: bool, verbose: bool, min_num_of_reads: int,
        cut_site_position: int, amplicon_min_score: float, min_read_length: int, config: Path,
        override_binomial_p: bool, confidence_interval: float, editing_threshold: float, suppress_site_output: bool,
        experiment_name: str, fastp_threads: int, allow_translocations: bool, max_error_on_primer: int,
        enable_substitutions: bool, ambiguous_cut_site_detection: bool, debug: bool, alignment_input, table_input):

    # Init the logger
    Logger.set_output_dir(output)
    if verbose:
        Logger.set_logger_level(logging.DEBUG)
    else:
        Logger.set_logger_level(logging.INFO)
    logger = Logger.get_logger()

    try:
        # Display welcome msg
        click.echo(welcome_msg)

        # Set config path and get configuration
        Configurator.set_cfg_path(config)
        cfg = Configurator.get_cfg()

        # Convert the amplicon_csv to a pandas.DataFrame
        ref_df = pd.read_csv(amplicons_csv, names=[SITE_NAME, REFERENCE, SGRNA, ON_TARGET, F_PRIMER, R_PRIMER])
        if ref_df[[SITE_NAME, REFERENCE, SGRNA, ON_TARGET]].isnull().values.any():
            raise NoneValuesInAmpliconsCSV()
        ref_df.index = ref_df[SITE_NAME]

        # Convert all bases to upper case
        ref_df[REFERENCE] = ref_df[REFERENCE].apply(lambda x : x.upper())
        if ref_df[REFERENCE].str.contains('[^ATGC]',regex=True).any():
            raise BadReferenceAmpliconChar()

        ref_df[SGRNA] = ref_df[SGRNA].apply(lambda x: x.upper())
        if ref_df[SGRNA].str.contains('[^ATGC]', regex=True).any():
            raise BadSgRNAChar()
        # TODO - HDR - Change to default coin, delete on-target results and add warning. read HDR and remove.

        # Create InputProcessing instance
        input_processing = InputProcessing(ref_df, output, amplicon_min_score, min_read_length, max_error_on_primer)

        # Find expected cut-site position
        input_processing.convert_sgRNA_to_cut_site_position(cut_site_position)

        # Detect ambiguous cut-sites
        input_processing.detect_ambiguous_cut_site(cut_site_position, ambiguous_cut_site_detection)

        # Filter low quality reads and merge pair-end reads with fastp
        if not override_fastp:
            tx_merged = input_processing.fastp(tx_in1, tx_in2, fastp_options_string, fastp_threads, ExpType.TX)
            mock_merged = input_processing.fastp(mock_in1, mock_in2, fastp_options_string, fastp_threads, ExpType.MOCK)
        else:
            logger.info("Skip merging with fastp and read merged files from input.")
            tx_merged, mock_merged = tx_in1, mock_in1

        # TODO - delete dump to pickle files.
        if alignment_input:
            with open(os.path.join(alignment_input, "tx_reads_d.pkl"), "rb") as file:
                tx_reads_d = pickle.load(file)
            with open(os.path.join(alignment_input, "mock_reads_d.pkl"), "rb") as file:
                mock_reads_d = pickle.load(file)
            tx_trans_df = pd.read_csv(os.path.join(alignment_input, "tx_translocation_reads.csv"))
            mock_trans_df = pd.read_csv(os.path.join(alignment_input, "mock_translocation_reads.csv"))
        else:
            # Demultiplexing reads
            tx_reads = input_processing.demultiplex_reads(tx_merged, ExpType.TX)
            mock_reads = input_processing.demultiplex_reads(mock_merged, ExpType.MOCK)

            # Align reads to the amplicon reference sequences
            tx_reads_d, tx_trans_df = input_processing.align_reads(tx_reads, ExpType.TX, allow_translocations, debug)
            mock_reads_d, mock_trans_df = input_processing.align_reads(mock_reads, ExpType.MOCK, allow_translocations, debug)

            # TODO - remove
            tx_trans_df.to_csv(os.path.join(output, "tx_translocation_reads.csv"), index=False)
            mock_trans_df.to_csv(os.path.join(output, "mock_translocation_reads.csv"), index=False)

            with open(os.path.join(output, "tx_reads_d.pkl"), "wb") as file:
                pickle.dump(tx_reads_d, file)
            with open(os.path.join(output, "mock_reads_d.pkl"), "wb") as file:
                pickle.dump(mock_reads_d, file)

       # TODO - uneven reads in Tx & Mock - put warning (and in final report)

        # Get modification types and positions priors
        modifications = ModificationTypes.init_from_cfg(enable_substitutions)

        # Convert alignment to modification tables
        logger.info("Convert alignment tables to algorithm input")
        tables_d: Dict[str, ModificationTables] = dict()

        # TODO - delete dump to pickle files.
        if table_input:
            with open(os.path.join(table_input, "tables_d.pkl"), "rb") as file:
                tables_d = pickle.load(file)
        else:
            for site, row in ref_df.iterrows():
                tx_reads_num = tx_reads_d[site][FREQ].sum()
                mock_reads_num = mock_reads_d[site][FREQ].sum()
                if min(tx_reads_num, mock_reads_num) < min_num_of_reads:
                    logger.info("Site {} - Discarded from evaluation due to low number of reads (treatment={}, "
                                "mock={}).".format(site, tx_reads_num, mock_reads_num))
                else:
                    tables_d[site] = ModificationTables(tx_reads_d[site], mock_reads_d[site], modifications, row)
                    logger.debug(
                        "Site {} - Converted. Number of reads (treatment={}, mock={}).".format(site, tx_reads_num,
                                                                                               mock_reads_num))
            with open(os.path.join(output, "tables_d.pkl"), "wb") as file:
                pickle.dump(tables_d, file)

        # Compute binomial coin for all modification types
        binom_p_d = compute_binom_p(tables_d, modifications, override_binomial_p, ref_df)

        result_summary_d: AlgResult = dict()  # Algorithm result dictionary
        # Run crispector algorithm on all sites
        for site, row in ref_df.iterrows():
            cut_site = row[CUT_SITE]

            # Continue if site was discarded
            if site not in tables_d:
                # Log the following in the result dict
                tx_reads_num = tx_reads_d[site][FREQ].sum()
                mock_reads_num = mock_reads_d[site][FREQ].sum()
                result_summary_d[site] = {TX_READ_NUM: tx_reads_num, MOCK_READ_NUM: mock_reads_num, ON_TARGET: row[ON_TARGET]}
                continue

            # Create output folder
            if not suppress_site_output:
                site_output = os.path.join(output, site)
                if not os.path.exists(site_output):
                    os.makedirs(site_output)
            else:
                site_output = None

            logger.info("Site {} - Evaluate editing activity".format(site))
            algorithm = CrispectorAlgorithm(site, experiment_name, cut_site, modifications, binom_p_d[site],
                                            confidence_interval, row[ON_TARGET], site_output)
            result_summary_d[site] = algorithm.evaluate(tables_d[site])
            result_summary_d[site][ON_TARGET] = row[ON_TARGET]
            logger.debug("Site {} - Done! Editing activity is {:.2f}".format(site,
                                                                             result_summary_d[site][EDIT_PERCENT]))


        # Dump summary results
        summary_result_df = pd.DataFrame.from_dict(result_summary_d, orient='index')
        summary_result_df[SITE_NAME] = summary_result_df.index
        summary_result_df = summary_result_df.reindex(ref_df.index, columns=SUMMARY_RESULTS_TITLES)
        summary_result_to_excel(summary_result_df, confidence_interval, output)

        # Create a text file with all discarded sites
        discarded_sites_text(summary_result_df, min_num_of_reads, output)

        # Create bar plot for editing activity
        plot_editing_activity(summary_result_df, confidence_interval, editing_threshold, output, experiment_name)

        # Dump reads statistics
        tx_input_n, tx_merged_n, tx_aligned_n = input_processing.read_numbers(ExpType.TX)
        mock_input_n, mock_merged_n, mock_aligned_n = input_processing.read_numbers(ExpType.MOCK)
        create_reads_statistics_report(summary_result_df, min_num_of_reads, tx_input_n, tx_merged_n, tx_aligned_n,
                                       mock_input_n, mock_merged_n, mock_aligned_n, output, experiment_name)

        # Remove fastp files
        if not keep_fastp_output:
            for exp_type in ExpType:
                fastp_output = os.path.join(output, FASTP_DIR[exp_type])
                if os.path.exists(fastp_output):
                    send2trash(fastp_output)

    except NoneValuesInAmpliconsCSV:
        logger.info("{}".format(traceback.format_exc()))
        logger.error("amplicons_csv file contains None values! Check the input file.")
        return 1
    except SgRNANotInReferenceSequence as e:
        logger.info("{}".format(traceback.format_exc()))
        logger.error("Site={} - sgRNA not in reference sequence!".format(e.site_name))
        return 2
    except FastpRunTimeError:
        logger.info("{}".format(traceback.format_exc()))
        logger.error("fastp failed ! check log for more details")
        return 3
    except CantOpenDemultiplexedSamFile as e:
        logger.info("{}".format(traceback.format_exc()))
        logger.error("Can't read demultiplexed sam file at {}. Please check path and permissions".format(e.path))
        return 4
    except ConfiguratorIsCalledBeforeInitConfigPath:
        logger.info("{}".format(traceback.format_exc()))
        logger.error("Configurator is called before set config path")
        return 5
    except PriorPositionHasWrongLength as e:
        logger.info("{}".format(traceback.format_exc()))
        logger.error("Configuration for modification type=({} size {}-{}) has a bad pos_prior length {} (expected={})."
                     .format(e.indel_type.name, e.indel_min, e.indel_max, e.actual_len, e.expected_len))
        return 6
    except UnknownAlignmentChar:
        logger.info("{}".format(traceback.format_exc()))
        logger.error("Unknown alignment character from Bio-python-Alignment. Please open an issue on GitHub.")
    except Bowtie2BuildRunTimeError:
        logger.info("{}".format(traceback.format_exc()))
        logger.error("Bowtie2-build failed ! check log for more details")
        return 7
    except Bowtie2RunTimeError:
        logger.info("{}".format(traceback.format_exc()))
        logger.error("Bowtie2 failed ! check log for more details")
        return 8
    except AlignerSubstitutionDoesntExist as e:
        logger.info("{}".format(traceback.format_exc()))
        logger.error("Configuration aligner substitution matrix {} Doesn't exist!".format(e.name))
        return 9
    except ClassificationFailed:
        logger.info("{}".format(traceback.format_exc()))
        logger.error("Classification failed!!! Please open an issue on GitHub.")
        return 10
    except BadSgRNAChar:
        logger.info("{}".format(traceback.format_exc()))
        logger.error("Bad character in sgRNA column! Check input!")
        return 11
    except BadReferenceAmpliconChar:
        logger.info("{}".format(traceback.format_exc()))
        logger.error("Bad character in REFERENCE amplicon column! Check input!")
        return 12
    except Exception:
        logger.info("{}".format(traceback.format_exc()))
        logger.error("Unknown Error. Open an issue on GitHub (#######)")  # Add GitHub link

        return -1
    # exit with success
    return 0
