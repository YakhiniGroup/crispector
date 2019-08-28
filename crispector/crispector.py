# -*- coding: utf-8 -*-
"""Main module."""

import sys
import pickle
import click
import logging
from send2trash import send2trash #TODO - add to poroject requiremnts
from algorithm_utils import compute_binom_p
from crispector_algorithm import CrispectorAlgorithm
from crispector_constants import welcome_msg, FASTP_DIR, FREQ, CUT_SITE, EDIT_PERCENT, TX_READ_NUM, MOCK_READ_NUM
from crispector_exceptions import FastpRunTimeError, NoneValuesInAmpliconsCSV, SgRNANotInReferenceSequence, \
    CantOpenMergedFastqFile, ConfiguratorIsCalledBeforeInitConfigPath, PriorPositionHasWrongLength
from crispector_types import ExpType, Path
from input_processing import InputProcessing
import traceback
from crispector_utils import Logger, Configurator
import os
import pandas as pd #TODO - add to poroject requiremnts
from crispector_constants import SITE_NAME, REFERENCE, SGRNA, ON_TARGET
from modification_tables import ModificationTables
from typing import Dict

from modification_types import ModificationTypes


def run(tx_in1: Path, tx_in2: Path, mock_in1: Path, mock_in2: Path, output: Path, amplicons_csv: Path,
        fastp_options_string: str, override_fastp: bool, keep_fastp_output: bool, verbose: bool, min_num_of_reads: int,
        cut_site_position: int, amplicon_min_alignment_score: float, override_alignment: bool, config: Path,
        override_binomial_p: bool, confidence_interval: float):

    # Init the logger
    Logger.set_log_path(output)
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
        ref_df = pd.read_csv(amplicons_csv, names=[SITE_NAME, REFERENCE, SGRNA, ON_TARGET])
        if ref_df.isnull().values.any():
            raise NoneValuesInAmpliconsCSV()
        ref_df = InputProcessing.convert_sgRNA_to_cut_site_position(ref_df, cut_site_position)

        # Filter low quality reads and merge pair-end reads with fastp
        if not override_fastp:
            tx_merged = InputProcessing.fastp(tx_in1, tx_in2, fastp_options_string, output, ExpType.TX)
            mock_merged = InputProcessing.fastp(mock_in1, mock_in2, fastp_options_string, output, ExpType.MOCK)
        else:
            logger.info("Skip merging with fastp and read merged files from input.")
            tx_merged, mock_merged = tx_in1, mock_in1

        # Align reads to the amplicon reference sequences
        tx_reads_d, tx_discarded_reads_num = InputProcessing.split_read_and_align(tx_merged, ref_df,
                                                                                  amplicon_min_alignment_score, output,
                                                                                  ExpType.TX, override_alignment)
        mock_reads_d, mock_discarded_reads_num = InputProcessing.split_read_and_align(mock_merged, ref_df,
                                                                                      amplicon_min_alignment_score,
                                                                                      output, ExpType.MOCK,
                                                                                      override_alignment)
        # Get modification types and positions priors
        modifications = ModificationTypes.init_from_cfg(cfg)

        # Convert alignment to modification tables
        logger.info("Convert alignment tables to algorithm input")
        tables_d: Dict[str, ModificationTables] = dict()
        summary_result_d = dict()  # Algorithm result dictionary

        # TODO - delete dump to pickle files.
        if os.path.exists(os.path.join(output, "tables.pkl")):
            with open(os.path.join(output, "tables.pkl"), "rb") as file:
                tables_d = pickle.load(file)
        else:
            for _, row in ref_df.iterrows():
                site = row[SITE_NAME]
                amplicon = row[REFERENCE]
                tx_reads_num = tx_reads_d[site][FREQ].sum()
                mock_reads_num = mock_reads_d[site][FREQ].sum()
                if min(tx_reads_num, mock_reads_num) < min_num_of_reads:
                    summary_result_d[site] = {TX_READ_NUM: tx_reads_num, MOCK_READ_NUM: mock_reads_num}
                    logger.info("Site {} - Discarded from evaluation due to low number of reads (treatment={}, "
                                "mock={}).".format(site, tx_reads_num, mock_reads_num))
                else:
                    tables_d[site] = ModificationTables(tx_reads_d[site], mock_reads_d[site], modifications, amplicon)
                    logger.debug(
                        "Site {} - Converted. Number of reads (treatment={}, mock={}).".format(site, tx_reads_num,
                                                                                               mock_reads_num))

            with open(os.path.join(output, "tables.pkl"), "wb") as file:
                pickle.dump(tables_d, file)

        # Compute binomial coin for all modification types
        # TODO -delete readDF from inputs
        binom_p_l = compute_binom_p(tables_d, modifications, override_binomial_p, ref_df)

        # Run crispector algorithm on all sites
        for _, row in ref_df.iterrows():
            site = row[SITE_NAME]
            cut_site = row[CUT_SITE]
            # Continue if site was discarded
            if site not in tables_d:
                continue
            # Create output folder
            site_output = os.path.join(output, site)
            if not os.path.exists(site_output):
                os.makedirs(site_output)

            logger.info("Site {} - Evaluate editing activity".format(site))
            algorithm = CrispectorAlgorithm(site, cut_site, modifications, binom_p_l, confidence_interval, site_output)
            summary_result_d[site] = algorithm.evaluate(tables_d[site])
            logger.debug("Site {} - Done! Editing activity is {:.2f}".format(site,
                                                                             summary_result_d[site][EDIT_PERCENT]))

        # Dump summary results TODO - create function and reorder columns
        summary_result_df = pd.DataFrame.from_dict(summary_result_d, orient='index')
        summary_result_df[SITE_NAME] = summary_result_df.index
        summary_result_df.sort_values(by=EDIT_PERCENT, ascending=False)
        summary_result_df.to_csv(os.path.join(output,"summary_results.csv"),
                                            index=False, float_format='%.4f')

        # TODO - Create final bar plot.
        # TODO - Add final report with numbers.

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
        logger.error("fastp failed ! check fastp log.")  # TODO - fastp has a log?
        return 3
    except CantOpenMergedFastqFile as e:
        if verbose:
            traceback.print_exc(file=sys.stdout)
        logger.error("Can't read merged fastq file at {}. Please check path and permissions".format(e.path))
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
    except Exception:
        if verbose:
            traceback.print_exc(file=sys.stdout)
        logger.error("Unknown Error. Open an issue on GitHub (#######)")  # Add GitHub link
        return -1
    # exit with success
    return 0
