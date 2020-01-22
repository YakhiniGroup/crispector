# -*- coding: utf-8 -*-
"""Main module."""

import pickle
import click
import logging
from algorithm.binomial_probability import compute_binom_p
from algorithm.core_algorithm import CoreAlgorithm
from input_processing.utils import read_exp_config_and_check_input
from utils.exceptions import FastpRunTimeError, SgRNANotInReferenceSequence, ConfiguratorIsCalledBeforeInitConfigPath, \
    PriorPositionHasWrongLength, UnknownAlignmentChar, \
    AlignerSubstitutionDoesntExist, ClassificationFailed, BadSgRNAChar, BadReferenceAmpliconChar, BadInputError
from utils.constants_and_types import Path, welcome_msg, FREQ, TX_READ_NUM, MOCK_READ_NUM, EDIT_PERCENT, \
    SITE_NAME, ON_TARGET, CUT_SITE, AlgResult, OUTPUT_DIR, SUMMARY_RESULTS_TITLES, \
    AlgResultDf, DONOR
from report.html_report import create_final_html_report
from input_processing.input_processing import InputProcessing
import traceback
from algorithm.translocations import translocations_test
from utils.logger import Logger
from utils.configurator import Configurator
from report.visualization_and_output import create_site_output, create_experiment_output
import os
import pandas as pd
from modifications.modification_tables import ModificationTables
from typing import Dict
from modifications.modification_types import ModificationTypes

def run(tx_in1: Path, tx_in2: Path, mock_in1: Path, mock_in2: Path, report_output: Path, experiment_config: Path,
        fastp_options_string: str, verbose: bool, min_num_of_reads: int,
        cut_site_position: int, amplicon_min_score: float, translocation_amplicon_min_score: float,
        min_read_length: int, crispector_config: Path, override_binomial_p: bool, confidence_interval: float,
        editing_threshold: float, translocation_p_value: float, suppress_site_output: bool,
        disable_translocations: bool, enable_substitutions: bool,
        ambiguous_cut_site_detection: bool, debug: bool, table_input: Path, keep_fastp_output: bool):

    output = os.path.join(report_output, OUTPUT_DIR)
    # Create output folder
    if not os.path.exists(output):
        os.makedirs(output)

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
        Configurator.set_cfg_path(crispector_config)

        ref_df = read_exp_config_and_check_input(experiment_config, tx_in1, tx_in2, mock_in1, mock_in2)

        donor = ref_df[DONOR].notnull().any() # Donor experiment flag
        if donor:
            logger.warning("Donor experiment (HDR). NHEJ activity or translocations won't be evaluated for on-target site.")

        # Create site output folders
        for site, row in ref_df.iterrows():
            site_output = os.path.join(output, site)
            if not os.path.exists(site_output):
                os.makedirs(site_output)

        # Create InputProcessing instance
        input_processing = InputProcessing(ref_df, output, amplicon_min_score, translocation_amplicon_min_score,
                                           min_read_length, cut_site_position, ambiguous_cut_site_detection,
                                           disable_translocations, fastp_options_string, debug, keep_fastp_output)
        if table_input is None:
            tx_reads_d, mock_reads_d, tx_trans_df, mock_trans_df = input_processing.run(tx_in1, tx_in2, mock_in1, mock_in2)

        # Get modification types and positions priors
        modifications = ModificationTypes.init_from_cfg(enable_substitutions)

        # Convert alignment to modification tables
        logger.info("Convert alignment tables to algorithm input")
        tables_d: Dict[str, ModificationTables] = dict()


        if table_input is not None: # TODO - delete dump to pickle files.
            with open(os.path.join(table_input, "tables_d.pkl"), "rb") as file:
                tables_d = pickle.load(file)
            if os.path.exists(os.path.join(table_input, "tx_translocations_reads.csv")):
                tx_trans_df = pd.read_csv(os.path.join(table_input, "tx_translocations_reads.csv"))
            else:
                tx_trans_df = pd.DataFrame()
            if os.path.exists(os.path.join(table_input, "mock_translocations_reads.csv")):
                mock_trans_df = pd.read_csv(os.path.join(table_input, "mock_translocations_reads.csv"))
            else:
                mock_trans_df  = pd.DataFrame()
        else:
            for site, row in ref_df.iterrows():
                tx_reads_num = tx_reads_d[site][FREQ].sum().astype(int)
                mock_reads_num = mock_reads_d[site][FREQ].sum().astype(int)
                if donor and row[ON_TARGET]:
                    logger.info("Site {} - Discarded from evaluation. On-target site with donor (HDR) evaluation is not"
                                "supported by crispector.".format(site))
                elif min(tx_reads_num, mock_reads_num) < min_num_of_reads:
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

        # Run crispector core algorithm on all sites
        logger.info("Start Evaluating editing activity for all sites")
        result_summary_d: AlgResult = dict()  # Algorithm result dictionary
        algorithm_d: Dict[str, CoreAlgorithm] = dict()
        for site, row in ref_df.iterrows():
            cut_site = row[CUT_SITE]
            # Continue if site was discarded
            if site not in tables_d:
                # Log the following in the result dict
                tx_reads_num = 0 #tx_reads_d[site][FREQ].sum().astype(int) TODO - change back
                mock_reads_num = 0 #mock_reads_d[site][FREQ].sum().astype(int) TODO - change back
                result_summary_d[site] = {TX_READ_NUM: tx_reads_num, MOCK_READ_NUM: mock_reads_num, ON_TARGET: row[ON_TARGET]}
                continue

            algorithm_d[site] = CoreAlgorithm(cut_site, modifications, binom_p_d[site], confidence_interval, row[ON_TARGET])
            result_summary_d[site] = algorithm_d[site].evaluate(tables_d[site])
            result_summary_d[site][ON_TARGET] = row[ON_TARGET]
            logger.debug("Site {} - Editing activity is {:.2f}".format(site, result_summary_d[site][EDIT_PERCENT]))

        # Convert result_summary dict to DataFrame
        summary_df: AlgResultDf = pd.DataFrame.from_dict(result_summary_d, orient='index')
        summary_df[SITE_NAME] = summary_df.index
        summary_df = summary_df.reindex(ref_df.index, columns=SUMMARY_RESULTS_TITLES)
        summary_df = summary_df.reset_index(drop=True)

        logger.info("Evaluating editing activity for all sites - Done!")

        # Run translocations test and call translocations reads
        logger.info("Translocations - Run HG tests")
        trans_result_df = translocations_test(summary_df, tx_trans_df , mock_trans_df, translocation_p_value,
                                              editing_threshold)
        logger.debug("Translocations - HG test - Done!")

        # Create plots & tables for all sites
        logger.info("Start creating experiment plots and tables")
        exp_param_d = create_experiment_output(summary_df,  tx_trans_df, mock_trans_df, trans_result_df,
                                               input_processing, min_num_of_reads, confidence_interval,
                                               editing_threshold, translocation_p_value, output)
        site_param_d = dict()
        if not suppress_site_output:
            for site, algorithm in algorithm_d.items():
                logger.debug("Site {} - Start creating plots and tables".format(site))

                # Create plots and tables
                site_output = os.path.join(output, site)
                site_param_d[site] = create_site_output(algorithm, modifications, tables_d[site], result_summary_d[site],
                                                        site, site_output)
        logger.info("Creating experiment plots and tables - Done!")

        # Create final HTML report
        with open(os.path.join(output, "exp_param_d.pkl"), "wb") as file:
            pickle.dump(exp_param_d, file)
        with open(os.path.join(output, "site_param_d.pkl"), "wb") as file:
            pickle.dump(site_param_d, file)
        create_final_html_report(exp_param_d, site_param_d, report_output)

    # Catch exceptions
    except BadInputError as e:
        logger.info("{}".format(traceback.format_exc()))
        logger.error("{}".format(e.msg))
        return 1
    except BadReferenceAmpliconChar:
        logger.info("{}".format(traceback.format_exc()))
        logger.error("Bad character in AmpliconReference or DonorReference columns! Check experiment_config.csv!")
        return 2
    except BadSgRNAChar:
        logger.info("{}".format(traceback.format_exc()))
        logger.error("Bad character in sgRNA column! Check experiment_config.csv!")
        return 3
    except SgRNANotInReferenceSequence as e:
        logger.info("{}".format(traceback.format_exc()))
        logger.error("Site={} - sgRNA is not in reference sequence!".format(e.site_name))
        return 4
    except FastpRunTimeError:
        logger.info("{}".format(traceback.format_exc()))
        logger.error("fastp failed! check log for more details")
        return 5
    except ConfiguratorIsCalledBeforeInitConfigPath:
        logger.info("{}".format(traceback.format_exc()))
        logger.error("Configurator is called before set config path")
        return 6
    except PriorPositionHasWrongLength as e:
        logger.info("{}".format(traceback.format_exc()))
        logger.error("Configuration for modification type=({} size {}-{}) has a bad pos_prior length {} (expected={})."
                     .format(e.indel_type.name, e.indel_min, e.indel_max, e.actual_len, e.expected_len))
        return 7
    except UnknownAlignmentChar:
        logger.info("{}".format(traceback.format_exc()))
        logger.error("Unknown alignment character from Bio-python-Alignment. Please contact crispector package author.")
    except AlignerSubstitutionDoesntExist as e:
        logger.info("{}".format(traceback.format_exc()))
        logger.error("Configuration aligner substitution matrix {} Doesn't exist!".format(e.name))
        return 8
    except ClassificationFailed:
        logger.info("{}".format(traceback.format_exc()))
        logger.error("Classification failed!!! Please contact crispector package author.")
        return 9
    except Exception:
        logger.info("{}".format(traceback.format_exc()))
        logger.error("Unknown Error. Please contact crispector package author")
        return -1
    # exit with success
    return 0
