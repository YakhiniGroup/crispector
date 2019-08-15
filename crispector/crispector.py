# -*- coding: utf-8 -*-

"""Main module."""
import sys
import click
import logging

from send2trash import send2trash
from constants import welcome_msg, TX_FASTP_DIR, MOCK_FASTP_DIR
from exceptions import FastpRunTimeError, NoneValuesInAmpliconsCSV, GeneralReadAmpliconsCSV, SgRNANotInReferenceSequence
from input_process_utils import fastp, convert_sgRNA_to_cut_site_position
import traceback
from utils import Logger
import os
import pandas as pd
from constants import SITE_NAME, REFERENCE, SGRNA, ON_TARGET


def run(tx_in1: str, tx_in2: str, mock_in1: str, mock_in2: str, output: str, amplicons_csv: str,
        fastp_options_string: str, override_fastp: bool, keep_fastp_output: bool, verbose: bool, min_num_of_reads: int,
        cut_site_position: int):

    # init the logger
    Logger.set_log_path(output)
    if verbose:
        Logger.set_logger_level(logging.DEBUG)
    else:
        Logger.set_logger_level(logging.INFO)
    logger = Logger.get_logger()

    try:
        # Display welcome msg
        click.echo(welcome_msg)

        # Convert amplicon_csv to pandas.DataFrame
        try:
            ref_df = pd.read_csv(amplicons_csv, names=[SITE_NAME, REFERENCE, SGRNA, ON_TARGET])
            if ref_df.isnull().values.any():
                raise NoneValuesInAmpliconsCSV()
            ref_df = convert_sgRNA_to_cut_site_position(ref_df, cut_site_position)
        except Exception:
            raise GeneralReadAmpliconsCSV()

        # Filter low quality reads and merge pair-end reads with fastp
        if not override_fastp:
            tx_merged = fastp(tx_in1, tx_in2, fastp_options_string, output, fastq_type="Treatment")
            mock_merged = fastp(mock_in1, mock_in2, fastp_options_string, output, fastq_type="Mock")
        else:
            tx_merged, mock_merged = tx_in1, mock_in1

        # Align reads to references sequences



        # Remove fastp files
        if not keep_fastp_output:
            for dir in [TX_FASTP_DIR, MOCK_FASTP_DIR]:
                fastp_output = os.path.join(output, dir)
                if os.path.exists(fastp_output):
                    send2trash(fastp_output)

        return 0
    except NoneValuesInAmpliconsCSV:
        if verbose:
            traceback.print_exc(file=sys.stdout)
        logger.error("amplicons_csv file contains None values! Check the input file.")
        return 1
    except GeneralReadAmpliconsCSV:
        if verbose:
            traceback.print_exc(file=sys.stdout)
        logger.error("amplicons_csv file read error! Check the input file.")
        return 2
    except SgRNANotInReferenceSequence as e:
        if verbose:
            traceback.print_exc(file=sys.stdout)
        logger.error("Site={} - sgRNA not in reference sequence!".format(e.site_name))
        return 3
    except FastpRunTimeError:
        if verbose:
            traceback.print_exc(file=sys.stdout)
        logger.error("fastp failed ! check fastp log.")  # TODO - fastp has a log?
        return 4
    except Exception:
        if verbose:
            traceback.print_exc(file=sys.stdout)
        logger.error("Unknown Error. Open an issue on GitHub (#######)")  # Add GitHub link
        return -1

    # exit with success
    return 0
