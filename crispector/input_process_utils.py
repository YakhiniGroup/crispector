import os
import subprocess

from constants import MOCK_FASTP_DIR, TX_FASTP_DIR, REFERENCE, SGRNA, COMPLEMENT, SITE_NAME, CUT_SITE
from exceptions import FastpRunTimeError, SgRNANotInReferenceSequence
from types import AmpliconDf
from utils import Logger
import edlib      # TODO - add to install requirments
from typing import List
import re



def fastp(in1: str, in2: str, fastp_options_string: str, output_dir: str, fastq_type: str) -> str:
    """
    Wrapper for fastp SW.
    :param in1: read1 input
    :param in2: read2 input
    :param fastp_options_string: any additional fastp options
    :param output_dir: output path
    :param fastq_type: Can be either Treatment or Mock
    :return: path to merged fastq file
    """
    logger = Logger.get_logger()

    if fastq_type == "Treatment":
        fastp_output = os.path.join(output_dir, TX_FASTP_DIR)
    else:
        fastp_output = os.path.join(output_dir, MOCK_FASTP_DIR)

    # Create output folder
    if not os.path.exists(fastp_output):
        os.makedirs(fastp_output)

    merged_path = os.path.join(fastp_output, "merged_reads.fastq")

    command = ["fastp", "-i", in1, "-I", in2, "-o", os.path.join(fastp_output, "r1_filtered_reads.fastq"),
               "-O", os.path.join(fastp_output, "r2_filtered_reads.fastq"), "-m", "--merged_out", merged_path,
               "-j", os.path.join(fastp_output, "fastp.json"), "-h", os.path.join(fastp_output, "fastp.html")]
    command = " ".join(command) + fastp_options_string

    logger.debug("fastp on {} - Command {}".format(fastq_type, command))
    logger.info("fastp on {} - Run.".format(fastq_type))
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError:
        raise FastpRunTimeError()

    logger.info("fastp on {} - Done.".format(fastq_type))

    return merged_path


# parse fastq file
def parse_fastq_file(file_name: str) -> List[str]:
    """
    process fastq_file to a list of strings
    :param file_name: fastq file name
    :return:
    """
    sequences = []
    with open(file_name) as fastq_file:
        line = fastq_file.readline()

        while line:
            line = fastq_file.readline()
            if re.match("[ACGT]+\Z", line[:-1]):
                sequences.append(line[:-1])

    return sequences


def reverseAndComplement(seq: str):
    return "".join(COMPLEMENT.get(base, base) for base in reversed(seq))


def get_expected_cut_site(reference: str , sgRNA: str, cut_site_position: int, site_name:str =""):
    """
    Find sgRNA (or the reverse complement) inside the reference and return the expected cut-site.
    The cut-site is LEFT to returned index
    :param reference: reference sequence
    :param sgRNA: sgRNA sequence
    :param cut_site_position: position relative to the PAM
    :param site_name: site name
    :return: expected cut-site
    """
    sgRNA_start_idx = reference.find(sgRNA)
    if sgRNA_start_idx == -1:
        sgRNA_start_idx = reference.find(reverseAndComplement(sgRNA))
        if sgRNA_start_idx != -1:
            raise SgRNANotInReferenceSequence(site_name)
        else:
            cut_site = sgRNA_start_idx - cut_site_position
    else:
        cut_site = sgRNA_start_idx + len(sgRNA) + cut_site_position

    return cut_site


def convert_sgRNA_to_cut_site_position(ref_df: AmpliconDf, cut_site_position: int) -> AmpliconDf:
    """
    :param ref_df: amplicons reference sequence df
    :param cut_site_position: position relative to the PAM
    :return: ref_df with new column - CUT_SITE
    """
    ref_df[CUT_SITE] = ref_df.apply(lambda row: get_expected_cut_site(row[REFERENCE], row[SGRNA],
                                                                      cut_site_position, row[SITE_NAME]), axis=1)
    return ref_df

"""

## TODO - continue writing alignment!!!!!
reads = prase_fastq_file(reads_file)

# Create the reads dataframe.
reads_df = df.DataFrame(data=reads, columns=['sequence'])

# Group similar reads togehter and count them.
self.seq_df = self.seq_df.groupby('z').size().to_frame('count').reset_index()

# Read design
design_df = pd.read_csv(desgin_file)
def align():

"""
