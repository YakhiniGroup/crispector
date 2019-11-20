# -*- coding: utf-8 -*-
"""Main module."""

import sys
import pickle
import click
import logging
from send2trash import send2trash #TODO - add to poroject requiremnts
from algorithm_utils import compute_binom_p
from crispector_algorithm import CrispectorAlgorithm
from exceptions import FastpRunTimeError, NoneValuesInAmpliconsCSV, SgRNANotInReferenceSequence, \
    ConfiguratorIsCalledBeforeInitConfigPath, PriorPositionHasWrongLength, \
    UnknownAlignmentChar, Bowtie2RunTimeError, Bowtie2BuildRunTimeError, CantOpenDemultiplexedSamFile
from constants_and_types import ExpType, Path, FASTP_DIR, welcome_msg, FREQ, TX_READ_NUM, MOCK_READ_NUM, EDIT_PERCENT, \
    SUMMARY_RESULTS_TITLES, SITE_NAME, REFERENCE, SGRNA, ON_TARGET, CUT_SITE, AlgResult
from input_processing import InputProcessing
import traceback
from utils import Logger, Configurator, plot_editing_activity, create_reads_statistics_report, summary_result_to_excel, \
    discarded_sites_text
import os
import pandas as pd #TODO - add to project requiremnts
from modification_tables import ModificationTables
from typing import Dict
# TODO - add bowtie2 to project requirements
from modification_types import ModificationTypes
# TODO simplesam - add t project - notice, not conda or bioconda for this. I think we can copy the code


def run(tx_in1: Path, tx_in2: Path, mock_in1: Path, mock_in2: Path, output: Path, amplicons_csv: Path,
        fastp_options_string: str, override_fastp: bool, keep_fastp_output: bool, verbose: bool, min_num_of_reads: int,
        cut_site_position: int, amplicon_min_alignment_score: float, override_alignment: bool, config: Path,
        override_binomial_p: bool, confidence_interval: float, editing_threshold: float, suppress_site_output: bool,
        experiment_name: str, fastp_threads: int, bowtie2_threads: int, bowtie2_options_string: str,
        override_bowtie2: bool):

    # TODO - of logger not on debug don't print fastp and bowtie2 messages.

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

        # TODO - Add code that check amplicon have only correct bases (and make everything upper).
        # TODO - Same for guides.

        # Convert the amplicon_csv to a pandas.DataFrame
        ref_df = pd.read_csv(amplicons_csv, names=[SITE_NAME, REFERENCE, SGRNA, ON_TARGET])
        if ref_df.isnull().values.any():
            raise NoneValuesInAmpliconsCSV()
        ref_df = InputProcessing.convert_sgRNA_to_cut_site_position(ref_df, cut_site_position)
        ref_df.index = ref_df[SITE_NAME]
        # canonial (index based) site names.
        canonial_names = {name: "site_{}".format(idx) for idx, name in enumerate(ref_df.index)}

        # Filter low quality reads and merge pair-end reads with fastp
        if not override_fastp:
            tx_merged, tx_input_n, tx_merged_n = InputProcessing.fastp(tx_in1, tx_in2, fastp_options_string,
                                                                       fastp_threads, output, ExpType.TX)

            mock_merged, mock_input_n, mock_merged_n = InputProcessing.fastp(mock_in1, mock_in2, fastp_options_string,
                                                                             fastp_threads, output, ExpType.MOCK)
        else:
            logger.info("Skip merging with fastp and read merged files from input.")
            tx_merged, mock_merged = tx_in1, mock_in1
            tx_input_n, mock_input_n = -1, -1
            tx_merged_n, mock_merged_n = -1, -1

        # Demultiplexing reads with Bowtie2
        if not override_bowtie2:
            tx_reads = InputProcessing.demultiplex_reads(tx_merged, ref_df, ExpType.TX)
            mock_reads = InputProcessing.demultiplex_reads(mock_merged, ref_df, ExpType.MOCK)
        else:
            tx_reads = tx_in1
            mock_reads = mock_in1

        allow_translocations = True # TODO - Change to user input

        # Align reads to the amplicon reference sequences
        tx_reads_d, tx_aligned_n = InputProcessing.align_reads(tx_reads, ref_df,
                                                               amplicon_min_alignment_score, output, ExpType.TX,
                                                               override_alignment)

        mock_reads_d, mock_aligned_n = InputProcessing.align_reads(mock_reads, ref_df,
                                                                   amplicon_min_alignment_score, output, ExpType.MOCK,
                                                                   override_alignment)

        # Get modification types and positions priors
        modifications = ModificationTypes.init_from_cfg(cfg)

        # Convert alignment to modification tables
        logger.info("Convert alignment tables to algorithm input")
        tables_d: Dict[str, ModificationTables] = dict()

        # TODO - delete dump to pickle files.
        if os.path.exists(os.path.join(output, "tables.pkl")) and override_alignment:
            with open(os.path.join(output, "tables.pkl"), "rb") as file:
                tables_d = pickle.load(file)
        else:
            for site, row in ref_df.iterrows():
                amplicon = row[REFERENCE]
                tx_reads_num = tx_reads_d[site][FREQ].sum()
                mock_reads_num = mock_reads_d[site][FREQ].sum()
                if min(tx_reads_num, mock_reads_num) < min_num_of_reads:
                    logger.info("Site {} - Discarded from evaluation due to low number of reads (treatment={}, "
                                "mock={}).".format(site, tx_reads_num, mock_reads_num))
                else:
                    tables_d[site] = ModificationTables(tx_reads_d[site], mock_reads_d[site], modifications, amplicon)
                    logger.debug(
                        "Site {} - Converted. Number of reads (treatment={}, mock={}).".format(site, tx_reads_num,
                                                                                               mock_reads_num))

        # Compute binomial coin for all modification types
        # TODO -delete readDF from inputs, output as well
        binom_p_l = compute_binom_p(tables_d, modifications, override_binomial_p, ref_df, output)

        result_summary_d: AlgResult = dict()  # Algorithm result dictionary
        # Run crispector algorithm on all sites
        for site, row in ref_df.iterrows():
            cut_site = row[CUT_SITE]

            # Continue if site was discarded
            if site not in tables_d:
                # Log the following in the result dict
                tx_reads_num = tx_reads_d[site][FREQ].sum()
                mock_reads_num = mock_reads_d[site][FREQ].sum()
                result_summary_d[site] = {TX_READ_NUM: tx_reads_num, MOCK_READ_NUM: mock_reads_num,
                                          ON_TARGET: row[ON_TARGET]}
                continue

            # Create output folder
            if not suppress_site_output:
                site_output = os.path.join(output, site)
                if not os.path.exists(site_output):
                    os.makedirs(site_output)
            else:
                site_output = None

            logger.info("Site {} - Evaluate editing activity".format(site))
            algorithm = CrispectorAlgorithm(site, experiment_name, cut_site, modifications, binom_p_l,
                                            confidence_interval, row[ON_TARGET], site_output)
            result_summary_d[site] = algorithm.evaluate(tables_d[site])
            result_summary_d[site][ON_TARGET] = row[ON_TARGET]
            logger.debug("Site {} - Done! Editing activity is {:.2f}".format(site,
                                                                             result_summary_d[site][EDIT_PERCENT]))

        with open(os.path.join(output, "tables.pkl"), "wb") as file:
            pickle.dump(tables_d, file)

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
        create_reads_statistics_report(summary_result_df, min_num_of_reads, tx_input_n, tx_merged_n, tx_aligned_n,
                                       mock_input_n, mock_merged_n, mock_aligned_n, output, experiment_name)

        # Remove fastp files
        if not keep_fastp_output:
            for exp_type in ExpType:
                fastp_output = os.path.join(output, FASTP_DIR[exp_type])
                if os.path.exists(fastp_output):
                    send2trash(fastp_output)

    except NoneValuesInAmpliconsCSV:
        if verbose:
            traceback.print_exc(file=sys.stdout)
        logger.error("amplicons_csv file contains None values! Check the input file.")
        return 1
    except SgRNANotInReferenceSequence as e:
        if verbose:
            traceback.print_exc(file=sys.stdout)
        logger.error("Site={} - sgRNA not in reference sequence!".format(e.site_name))
        return 2
    except FastpRunTimeError:
        if verbose:
            traceback.print_exc(file=sys.stdout)
        logger.error("fastp failed ! check log for more details")
        return 3
    except CantOpenDemultiplexedSamFile as e:
        if verbose:
            traceback.print_exc(file=sys.stdout)
        logger.error("Can't read demultiplexed sam file at {}. Please check path and permissions".format(e.path))
        return 4
    except ConfiguratorIsCalledBeforeInitConfigPath:
        if verbose:
            traceback.print_exc(file=sys.stdout)
        logger.error("Configurator is called before set config path")
        return 5
    except PriorPositionHasWrongLength as e:
        if verbose:
            traceback.print_exc(file=sys.stdout)
        logger.error("Configuration for modification type=({} size {}-{}) has a bad pos_prior length {} (expected={})."
                     .format(e.indel_type.name, e.indel_min, e.indel_max, e.actual_len, e.expected_len))
        return 6
    except UnknownAlignmentChar:
        if verbose:
            traceback.print_exc(file=sys.stdout)
        logger.error("Unknown alignment character from Bio-python-Alignment. Please open a defect on GitHub.")
    except Bowtie2BuildRunTimeError:
        if verbose:
            traceback.print_exc(file=sys.stdout)
        logger.error("Bowtie2-build failed ! check log for more details")
        return 7
    except Bowtie2RunTimeError:
        if verbose:
            traceback.print_exc(file=sys.stdout)
        logger.error("Bowtie2 failed ! check log for more details")
        return 8
    except Exception:
        if verbose:
            traceback.print_exc(file=sys.stdout)
        logger.error("Unknown Error. Open an issue on GitHub (#######)")  # Add GitHub link
        return -1
    # exit with success
    return 0
