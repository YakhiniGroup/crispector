# -*- coding: utf-8 -*-
"""Main module."""

import sys
import click
import logging
from send2trash import send2trash #TODO - add to poroject requiremnts
from crispector_constants import welcome_msg, FASTP_DIR, FREQ
from crispector_exceptions import FastpRunTimeError, NoneValuesInAmpliconsCSV, SgRNANotInReferenceSequence, \
    CantOpenMergedFastqFile, ConfiguratorIsCalledBeforeInitConfigPath
from crispector_types import ExpType, Path
from input_process_utils import InputProcess
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
        cut_site_position: int, amplicon_min_alignment_score: float, override_alignment: bool):

    # TODO - change all str to DNASeq a
    # TODO - change to option
    cfg_path = ""

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
        Configurator.set_cfg_path(cfg_path)
        cfg = Configurator.get_cfg()
        # TODO - Handle priors
        modifications = ModificationTypes.init_from_cfg(cfg)

        # Convert the amplicon_csv to a pandas.DataFrame
        ref_df = pd.read_csv(amplicons_csv, names=[SITE_NAME, REFERENCE, SGRNA, ON_TARGET])
        if ref_df.isnull().values.any():
            raise NoneValuesInAmpliconsCSV()
        ref_df = InputProcess.convert_sgRNA_to_cut_site_position(ref_df, cut_site_position)

        # Filter low quality reads and merge pair-end reads with fastp
        if not override_fastp:
            tx_merged = InputProcess.fastp(tx_in1, tx_in2, fastp_options_string, output, ExpType.TX)
            mock_merged = InputProcess.fastp(mock_in1, mock_in2, fastp_options_string, output, ExpType.MOCK)
        else:
            logger.info("Skip merging with fastp and read merged files from input.")
            tx_merged, mock_merged = tx_in1, mock_in1

        # Align reads to the amplicon reference sequences
        tx_reads_d, tx_discarded_reads_num = InputProcess.split_read_and_align(tx_merged, ref_df,
                                                                               amplicon_min_alignment_score, output,
                                                                               ExpType.TX,
                                                                               override_alignment)
        mock_reads_d, mock_discarded_reads_num = InputProcess.split_read_and_align(mock_merged, ref_df,
                                                                                   amplicon_min_alignment_score,
                                                                                   output, ExpType.MOCK,
                                                                                   override_alignment)
        # Convert alignment to modification tables
        site_d = {k: v for (k, v) in zip(ref_df[SITE_NAME].values, ref_df[REFERENCE].values)}
        tables_d: Dict[ModificationTables] = dict()

        for site, amplicon in site_d.items():
            tx_reads_num = tx_reads_d[site][FREQ].sum()
            mock_reads_num = mock_reads_d[site][FREQ].sum()
            if min(tx_reads_num, mock_reads_num) < min_num_of_reads:
                # TODO - add information in final report?
                logger.info("Site {} - Discarded from evaluation due to low number of reads (treatment={}, "
                            "mock={}).".format(site, tx_reads_num, mock_reads_num))
            else:
                logger.debug("Site {} - Number of reads (treatment={}, mock={}).".format(site, tx_reads_num,
                                                                                         mock_reads_num))
                tables_d[site] = ModificationTables(tx_reads_d[site], mock_reads_d[site], modifications, amplicon)

        # Remove fastp files
        if not keep_fastp_output:
            for exp_type in ExpType:
                fastp_output = os.path.join(output, FASTP_DIR[exp_type])
                if os.path.exists(fastp_output):
                    send2trash(fastp_output)

        return 0
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
    except Exception:
        if verbose:
            traceback.print_exc(file=sys.stdout)
        logger.error("Unknown Error. Open an issue on GitHub (#######)")  # Add GitHub link
        return -1

    # exit with success
    return 0
