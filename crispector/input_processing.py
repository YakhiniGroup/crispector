import os
import subprocess
from crispector_constants import REFERENCE, SGRNA, COMPLEMENT, SITE_NAME, CUT_SITE, READ, \
    FREQ, FASTP_DIR, NONE, CIGAR, ALIGNMENT_W_INS, ALIGNMENT_W_DEL
from crispector_exceptions import FastpRunTimeError, SgRNANotInReferenceSequence, CantOpenMergedFastqFile
from crispector_types import AmpliconDf, ReadsDict, ExpType, ReadsDf, IndelType, Path, DNASeq
from crispector_utils import Logger
import edlib      # TODO - add to install requirements
from typing import List, Tuple
import re
import pandas as pd  # TODO - add to install requirements


class InputProcessing:
    """
    A container static class with all relevant function to process crispector input.
    """

    @staticmethod
    def fastp(in1: Path, in2: Path, fastp_options_string: Path, output_dir: Path, fastq_type: ExpType) -> Path:
        """
        Wrapper for fastp SW.
        :param in1: read1 input
        :param in2: read2 input
        :param fastp_options_string: any additional fastp options
        :param output_dir: output path
        :param fastq_type: ExpType
        :return: path to merged fastq file
        """
        logger = Logger.get_logger()

        # Create output folder
        fastp_output = os.path.join(output_dir, FASTP_DIR[fastq_type])
        if not os.path.exists(fastp_output):
            os.makedirs(fastp_output)

        merged_path = os.path.join(fastp_output, "merged_reads.fastq")

        command = ["fastp", "-i", in1, "-I", in2, "-o", os.path.join(fastp_output, "r1_filtered_reads.fastq"),
                   "-O", os.path.join(fastp_output, "r2_filtered_reads.fastq"), "-m", "--merged_out", merged_path,
                   "-j", os.path.join(fastp_output, "fastp.json"), "-h", os.path.join(fastp_output, "fastp.html")]
        command = " ".join(command) + fastp_options_string

        logger.debug("fastp for {} - Command {}".format(fastq_type.name(), command))
        logger.info("fastp for {} - Run (may take a few minutes).".format(fastq_type.name()))
        try:
            subprocess.run(command, shell=True, check=True)
        except subprocess.CalledProcessError:
            raise FastpRunTimeError()

        logger.info("fastp for {} - Done.".format(fastq_type.name()))

        return merged_path

    @classmethod
    def split_read_and_align(cls, merged_fastq: Path, ref_df: AmpliconDf, min_score: float, output: Path,
                             exp_type: ExpType, override_alignment: bool) -> Tuple[ReadsDict, int]:
        # TODO - delete override_alignment
        """
        Split all reads and align them to their amplicon sequence:
        1. Find for every read in merged_fastq the alignment with min_edit_distance.
        2. If the min_edit_distance is below the min_distance threshold - Mark the read with this amplicon.
        3. For all matched reads - Compute cigar path.
        4. Dump unmatched to fastq file.
        :param merged_fastq: merged fastq path
        :param ref_df: AmpliconDf type
        :param min_score: user min alignment score (0-100)
        :param output: output path
        :param exp_type: ExpType
        :return: A tuple of (dictionary with key=SITE_NAME and value=ReadsDf, number of unaligned reads)
        """

        logger = Logger.get_logger()

        amplicons = ref_df[REFERENCE].values
        amplicon_names = ref_df[SITE_NAME].values
        amplicon_d = {k: v for (k, v) in zip(amplicon_names, amplicons)}
        unaligned_reads_num = 0

        # TODO - remove when delete override_alignment
        if override_alignment:
            logger.debug("Alignment for {} - Skip alignment and read from input.".format(exp_type.name()))
            reads_df = pd.read_csv(merged_fastq)
        else:
            # Read file and convert to pandas
            logger.debug("Alignment for {} - Read merged fastq from {}".format(exp_type.name(), merged_fastq))
            reads = cls._parse_fastq_file(merged_fastq)
            reads_df = pd.DataFrame(data=reads, columns=[READ])

            # Group identical reads together
            reads_df = reads_df.groupby(READ).size().to_frame(FREQ).reset_index()

            # Find amplicons and reads match
            logger.info("Alignment for {} - Start alignment for {} sites.".format(exp_type.name(), ref_df.shape[0]))
            reads_df[SITE_NAME] = reads_df.apply((lambda row: cls._match_by_edit_distance(row[READ], amplicons,
                                                  amplicon_names, min_score)), axis=1)

            # Dump unaligned reads
            unaligned_reads_num = cls._dump_unaligned_reads(reads_df, output, exp_type)
            reads_df.drop(reads_df.loc[reads_df[SITE_NAME] == NONE].index, inplace=True)

            # Compute cigar path
            reads_df[CIGAR] = reads_df.apply(
                (lambda row: cls._get_cigar_path(row[READ], amplicon_d[row[SITE_NAME]])), axis=1)

            # Convert cigar path to full alignment
            reads_df["alignment"] = reads_df.apply(
                (lambda row: cls._cigar_to_full_alignment(amplicon_d[row[SITE_NAME]], row[READ], row[CIGAR])), axis=1)
            reads_df[[ALIGNMENT_W_INS, ALIGNMENT_W_DEL]] = pd.DataFrame(reads_df['alignment'].tolist(), index=reads_df.index)
            reads_df.drop(columns=["alignment"], inplace=True)

            # TODO - move from here to the end of CRISPECTOR run
            aligned_read_path = os.path.join(output, "{}_aligned.csv".format(exp_type.name()))
            reads_df.to_csv(aligned_read_path, index=False)

        # Split read_df to all the different sites
        read_d: ReadsDict = dict()
        for site in amplicon_names:
            read_d[site] = reads_df.loc[reads_df[SITE_NAME] == site].sort_values(by=[FREQ], ascending=False).reset_index()

        logger.info("Alignment for {} - Done.".format(exp_type.name(), ref_df.shape[0]))
        return read_d, unaligned_reads_num

    @staticmethod
    def _dump_unaligned_reads(reads_df: ReadsDf, output: Path, exp_type: ExpType) -> int:
        """
        Store all unaligned reads in a fasta format.
        :param reads_df: aligned reads df.
        :param output: output_path
        :param exp_type:
        :return: number of unaligned reads.
        """
        logger = Logger.get_logger()

        total_reads_num = reads_df[FREQ].sum()
        unaligned_df = reads_df.loc[reads_df[SITE_NAME] == NONE]
        unaligned_reads_num = unaligned_df[FREQ].sum()

        logger.info("Alignment for {} - {} reads weren't aligned ({:.2f}% of all reads)".format(exp_type.name(),
                    unaligned_reads_num, unaligned_reads_num/total_reads_num))

        with open(os.path.join(output,"{}_unaligned_reads.fasta".format(exp_type.name())), 'w') as file:
            for index, row in unaligned_df.iterrows():
                file.write("> unaligned read index {} - read has {} copies in the original fastq file\n".format(index,
                           row[FREQ]))
                file.write("{}\n".format(row[READ]))

        return unaligned_reads_num

    @staticmethod
    def _match_by_edit_distance(read: DNASeq, amplicons: List[DNASeq], amplicon_names: List[str],
                                min_score: float) -> str:
        """
        Find for every read the alignment with min_edit_distance.
        If the min_edit_distance is below the min_distance threshold - Mark the read with this amplicon.
        :param read:  reads
        :param amplicons: List of reference sequences
        :param amplicon_names: List of site names
        :param min_score: user min alignment score (0-100)
        :return: site_name (or NONE string)
        """
        min_name = amplicon_names[0]
        min_amplicon = amplicons[0]
        min_dist = edlib.align(read, amplicons[0])['editDistance']

        for amplicon, name in zip(amplicons[1:], amplicon_names[1:]):
            d = edlib.align(read, amplicon, k=min_dist)['editDistance']
            if d < min_dist and d != -1:
                min_dist = d
                min_name = name
                min_amplicon = amplicon

        # check if min distance above min_distance threshold
        score = 100*((len(min_amplicon) - min_dist) / len(min_amplicon))

        if score < min_score:
            return NONE
        else:
            return min_name

    @staticmethod
    def _get_cigar_path(read: DNASeq, reference_sequence: DNASeq) -> str:
        """
        Compute cigar path for a read.
        :param read:  reads
        :param reference_sequence:  amplicon reference sequences
        :return: cigar path (string)
        """
        return edlib.align(read, reference_sequence, task='path')['cigar']

    @staticmethod
    def _parse_fastq_file(file_name: Path) -> List[DNASeq]:
        """
        process fastq_file to a list of DNA strings
        :param file_name: fastq file name
        :return: list of DNA strings
        """
        sequences = []
        try:
            with open(file_name) as fastq_file:
                line = fastq_file.readline()

                while line:
                    line = fastq_file.readline()
                    if re.match("[ACGT]+\Z", line[:-1]):
                        sequences.append(line[:-1])
        except IOError:
            raise CantOpenMergedFastqFile(file_name)

        return sequences

    @staticmethod
    def reverse_and_complement(seq: DNASeq):
        return "".join(COMPLEMENT.get(base, base) for base in reversed(seq))

    @classmethod
    def _get_expected_cut_site(cls, reference: DNASeq, sgRNA: DNASeq, cut_site_position: int, site_name: str = '') \
        -> int:
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
            sgRNA_start_idx = reference.find(cls.reverse_and_complement(sgRNA))
            if sgRNA_start_idx == -1:
                raise SgRNANotInReferenceSequence(site_name)
            else:
                cut_site = sgRNA_start_idx - cut_site_position
        else:
            cut_site = sgRNA_start_idx + len(sgRNA) + cut_site_position

        return cut_site

    @classmethod
    def convert_sgRNA_to_cut_site_position(cls, ref_df: AmpliconDf, cut_site_position: int) -> AmpliconDf:
        """
        :param ref_df: amplicons reference sequence df
        :param cut_site_position: position relative to the PAM
        :return: ref_df with new column - CUT_SITE
        """
        ref_df[CUT_SITE] = ref_df.apply(lambda row: cls._get_expected_cut_site(row[REFERENCE], row[SGRNA],
                                        cut_site_position, row[SITE_NAME]), axis=1)
        return ref_df

    @classmethod
    def _cigar_to_full_alignment(cls, ref: DNASeq, read: DNASeq, cigar: str):
        """
        convert cigar path to full alignment
        :param ref: reference
        :param read: read
        :param cigar: cigar path
        :return:
        """
        read_idx = 0
        ref_idx = 0
        align_read = read
        align_ref = ref
        for length, indel in cls.parse_cigar(cigar):
            if indel in [IndelType.MATCH, IndelType.SUB]:
                read_idx += length
                ref_idx += length
            elif indel == IndelType.INS:
                read_idx += length
                align_ref = cls.insert_dash(align_ref, ref_idx, length)
            elif indel == IndelType.DEL:
                align_read = cls.insert_dash(align_read, read_idx, length)
                ref_idx += length

        return align_ref, align_read

    @staticmethod
    def parse_cigar(cigar: str) -> Tuple[int, IndelType]:
        """
        Generator function for cigar path.
        Yield each iteration the length of the of the indel and indel type
        :param cigar:
        :return: Yield each iteration the length of the of the indel and indel type
        """
        for length, indel in re.findall(r'(\d+)([DIX=])', cigar):
            length = int(length)
            indel = IndelType.from_cigar(indel)
            yield length, indel

    @staticmethod
    def insert_dash(string: str, index: int, length: int) -> str:
        """
        insert dashes at the index position
        :param string:
        :param index:
        :param length:
        :return:
        """
        return string[:index] + length*'-' + string[index:]

