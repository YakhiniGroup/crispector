from crispector.utils.constants_and_types import COMPLEMENT, IndelType, Path, DNASeq, CIGAR_D, CIGAR_I, CIGAR_S, CIGAR_M, \
    AmpliconDf, SITE_NAME, REFERENCE, SGRNA, ON_TARGET, F_PRIMER, R_PRIMER, TX_IN1, TX_IN2, MOCK_IN1, MOCK_IN2, DONOR
from typing import List, Tuple
import re
from crispector.utils.exceptions import CantOpenMergedFastqFile, BadInputError, BadReferenceAmpliconChar, BadSgRNAChar
import pandas as pd
import os
import binascii
import gzip


def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return binascii.hexlify(test_f.read(2)) == b'1f8b'


def read_DNA_lines(fastq_file) -> List[DNASeq]:
    sequences = []
    line = fastq_file.readline()

    while line:
        line = fastq_file.readline()
        if re.match("[ACGT]+\Z", line[:-1]):
            sequences.append(line[:-1])
    return sequences


def reverse_complement(seq: DNASeq) -> DNASeq:
    return "".join(COMPLEMENT.get(base, base) for base in reversed(seq))


def parse_fastq_file(file_name: Path) -> List[DNASeq]:
    """
    process fastq_file to a list of DNA strings
    :param file_name: fastq file name
    :return: list of DNA strings
    """

    try:
        if is_gz_file(file_name):
            fastq_file = gzip.open(file_name, 'rt')
            return read_DNA_lines(fastq_file)
        else:
            with open(file_name) as fastq_file:
                return read_DNA_lines(fastq_file)
    except IOError:
        raise CantOpenMergedFastqFile(file_name)


def parse_cigar(cigar: str) -> Tuple[int, IndelType]:
    """
    Generator function for cigar path.
    Yield each iteration the length of the of the indel and indel type
    :param cigar:
    :return: Yield each iteration the length of the of the indel and indel type
    """
    for length, indel in re.findall(r'(\d+)([{}{}{}{}])'.format(CIGAR_D, CIGAR_I, CIGAR_S, CIGAR_M), cigar):
        length = int(length)
        indel = IndelType.from_cigar(indel)
        yield length, indel


def parse_cigar_with_mixed_indels(cigar: str) -> List[Tuple[int, int, IndelType, List[Tuple[int, IndelType]]]]:
    """
    Function returns a list of Tuple[indel length, indel length without insertion, IndelType, Mixed_list] where
    adjacent indels are aggregated into indel_type = Mixed.
    :param cigar:
    :return: list of Tuple[indel length, indel length without insertion, IndelType, Mixed list]
    """
    indel_list = []
    prev_indel = IndelType.MATCH
    prev_length = 0
    prev_length_wo_ins = 0  # used to understand with positions are relevant for this modification
    mixed_list = []  # List of all indels comprising current mixed
    mixed_count = 0  # Current mixed indel count of comprising modifications
    for length, indel in re.findall(r'(\d+)([{}{}{}{}])'.format(CIGAR_D, CIGAR_I, CIGAR_S, CIGAR_M), cigar):
        indel = IndelType.from_cigar(indel)
        length = int(length)
        length_wo_ins = int(length) if indel != IndelType.INS else 0
        mixed_list.append((length, indel))
        if (indel != IndelType.MATCH) and (prev_indel != IndelType.MATCH):
            mixed_count += 1
            length += prev_length
            length_wo_ins += prev_length_wo_ins
            indel_list[-1] = (length, length_wo_ins, IndelType.MIXED, mixed_list[-(mixed_count + 1):])
        else:
            indel_list.append((length, length_wo_ins, indel, []))
            mixed_count = 0

        # update prev
        prev_indel = indel
        prev_length = length
        prev_length_wo_ins = length_wo_ins

    return indel_list

def read_exp_config_and_check_input(experiment_config: Path, tx_in1: Path, tx_in2: Path, mock_in1: Path,
                                    mock_in2: Path) -> AmpliconDf:
    # Convert the experiment_config.csv to a pandas.DataFrame
    ref_df: AmpliconDf = pd.read_csv(experiment_config, names=[SITE_NAME, REFERENCE, SGRNA, ON_TARGET, F_PRIMER,
                                                               R_PRIMER, TX_IN1, TX_IN2, MOCK_IN1, MOCK_IN2, DONOR],
                                     skiprows=1, engine='python')

    if ((tx_in1 is None) and (mock_in1 is not None)) or ((tx_in1 is not None) and (mock_in1 is None)):
        raise BadInputError("--tx_in1 & --mock_in1 must be specified together")

    if ((tx_in2 is None) and (mock_in2 is not None)) or ((tx_in2 is not None) and (mock_in2 is None)):
        raise BadInputError("--tx_in2 & --mock_in2 must be specified together")

    if (tx_in2 is not None) and (tx_in1 is None):
        raise BadInputError("--tx_in1 & --mock_in1 must be specified when --tx_in2 & --mock_in2 are")

    for col in [SITE_NAME, REFERENCE, SGRNA, ON_TARGET]:
        if ref_df[col].isna().any():
            raise BadInputError("experiment_config.csv file has missing values in column {}".format(col))

    if tx_in1 is None:
        for col in [TX_IN1, MOCK_IN1]:
            if ref_df[col].isna().any():
                raise BadInputError("experiment_config.csv file has missing values in column {}. But --tx_in1 is not "
                                    "specified!!!".format(col))
            for idx, path in enumerate(ref_df[col]):
                if not os.path.exists(path):
                    raise BadInputError("experiment_config.csv, at row={}, column={}  - Path doesn't exist!"
                                        "\nPath={}".format(idx, col, path))
        if not (ref_df[TX_IN2].isna().all() and ref_df[MOCK_IN2].isna().all()):
            for col in [TX_IN2, MOCK_IN2]:
                if ref_df[col].isna().any():
                    raise BadInputError("experiment_config.csv file has missing values in column {}. But --tx_in1 is "
                                        "not specified!!!".format(col))
                for idx, path in enumerate(ref_df[col]):
                    if not os.path.exists(path):
                        raise BadInputError("experiment_config.csv, at row={}, column={}  - Path doesn't exist!"
                                            "\nPath={}".format(idx, col, path))
    else:
        for col in [TX_IN1, TX_IN2, MOCK_IN1, MOCK_IN2]:
            if ref_df[col].notnull().any():
                raise BadInputError("experiment_config.csv file contains values in column {}. But --tx_in1 is "
                                    "specified!!! leave this column empty".format(col))

    # Convert all bases to upper case
    ref_df[REFERENCE] = ref_df[REFERENCE].apply(lambda x: x.upper())
    if ref_df[REFERENCE].str.contains('[^ATGC]', regex=True).any():
        raise BadReferenceAmpliconChar()

    ref_df[SGRNA] = ref_df[SGRNA].apply(lambda x: x.upper())
    if ref_df[SGRNA].str.contains('[^ATGC]', regex=True).any():
        raise BadSgRNAChar()

    # Donor
    if ref_df[DONOR].notnull().any():
        ref_df.loc[ref_df[DONOR].notnull(), DONOR] = ref_df.loc[ref_df[DONOR].notnull(), DONOR].apply(lambda x: x.upper())
        if ref_df.loc[ref_df[DONOR].notnull(), DONOR].str.contains('[^ATGC]', regex=True).any():
            raise BadReferenceAmpliconChar()

    ref_df.index = ref_df[SITE_NAME]
    return ref_df

