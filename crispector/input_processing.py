import os
import subprocess
from constants import REFERENCE, SGRNA, COMPLEMENT, SITE_NAME, CUT_SITE, READ, \
    FREQ, CIGAR, ALIGNMENT_W_INS, ALIGNMENT_W_DEL, INDEL_COLS, DEL_LEN, DEL_START, DEL_END, SUB_POS, \
    SUB_CNT, INS_LEN, INS_POS, CIGAR_D, CIGAR_I, CIGAR_S, CIGAR_M, SCORE, ALIGN_SCORE
from exceptions import FastpRunTimeError, SgRNANotInReferenceSequence, CantOpenMergedFastqFile
from enum_types import AmpliconDf, ReadsDict, ExpType, ReadsDf, IndelType, Path, DNASeq, CigarPath, FASTP_DIR
from utils import Logger, Configurator
import edlib      # TODO - add to install requirements
from typing import List, Tuple, Dict
import re
import pandas as pd  # TODO - add to install requirements
from Bio import Align  # TODO - add to install requirements
from Bio.SubsMat import MatrixInfo
from collections import defaultdict

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
        amplicon_d = {k: v for (k, v) in zip(ref_df[SITE_NAME].values, ref_df[REFERENCE].values)}
        cut_site_d = {k: v for (k, v) in zip(ref_df[SITE_NAME].values, ref_df[CUT_SITE].values)}
        unaligned_reads_num = 0
        cfg = Configurator.get_cfg()
        window_size = cfg["window_size"]

        # TODO - remove when delete override_alignment
        if override_alignment:
            logger.debug("Alignment for {} - Skip alignment and read from input.".format(exp_type.name()))
            reads_df = pd.read_csv(merged_fastq)
        else:
            # Read file and convert to pandas
            logger.debug("Alignment for {} - Read merged fastq from {}".format(exp_type.name(), merged_fastq))
            reads = cls._parse_fastq_file(merged_fastq)
            reads_df: ReadsDf = pd.DataFrame(data=reads, columns=[READ])

            # Group identical reads together
            reads_df = reads_df.groupby(READ).size().to_frame(FREQ).reset_index()

            # Match reads and amplicons
            logger.info("Alignment for {} - Start alignment for {} sites.".format(exp_type.name(), ref_df.shape[0]))
            logger.debug("Alignment for {} - Demultiplexing reads for {} sites.".format(exp_type.name(),
                                                                                        ref_df.shape[0]))
            # prepare all amplicons and reverse amplicons
            reference_reverse = ref_df.apply(lambda row: cls.reverse_and_complement(row[REFERENCE]), axis=1)
            all_amplicons = list(ref_df[REFERENCE].values) + list(reference_reverse)
            all_amplicon_names = list(ref_df[SITE_NAME].values) + ["rev_" + x for x in ref_df[SITE_NAME].values]
            # Run matching
            reads_df[SITE_NAME] = reads_df.apply((lambda row: cls._match_by_edit_distance(row[READ], all_amplicons,
                                                  all_amplicon_names)), axis=1)

            # Reverse all the reversed reads
            rev_reads = reads_df.loc[reads_df[SITE_NAME].str.startswith("rev_"), READ]
            reads_df.loc[rev_reads.index, READ] = rev_reads.apply(cls.reverse_and_complement)
            rev_names = reads_df.loc[rev_reads.index, SITE_NAME].values
            reads_df.loc[rev_reads.index, SITE_NAME] = [name[4:] for name in rev_names]

            logger.debug("Alignment for {} - Demultiplexing Done!".format(exp_type.name()))

            # Align reads to their amplicon
            logger.debug("Alignment for {} - Start Needleman-Wunsch alignment for all reads.".format(exp_type.name()))

            cls._align_reads_to_amplicon(reads_df, amplicon_d)

            # Filter reads with low alignment score
            unaligned_reads_num = cls._filter_low_score_reads(reads_df, output, exp_type, min_score, amplicon_d)
            logger.debug("Alignment for {} - Needleman-Wunsch alignment done.".format(exp_type.name()))

            # TODO - shift reads into cut-site, same score.

            # TODO - Move this code Convert to the function of shift reads into cut-site
            # reads_df["alignment"] = reads_df.apply(
            #     (lambda row: cls._cigar_to_full_alignment(amplicon_d[row[SITE_NAME]], row[READ], row[CIGAR])), axis=1)
            # reads_df[[ALIGNMENT_W_INS, ALIGNMENT_W_DEL]] = pd.DataFrame(reads_df['alignment'].tolist(), index=reads_df.index)
            # reads_df.drop(columns=["alignment"], inplace=True)

            # TODO - move from here to the end of CRISPECTOR run
            aligned_read_path = os.path.join(output, "{}_aligned.csv".format(exp_type.name()))
            reads_df.to_csv(aligned_read_path, index=False)

        # Split read_df to all the different sites
        read_d: ReadsDict = dict()
        for site in amplicon_d.keys():
            read_d[site] = reads_df.loc[reads_df[SITE_NAME] == site].sort_values(by=[FREQ],
                                                                                 ascending=False).reset_index(drop=True)
            # Add indels columns to reads df
            cls._add_indels_columns(read_d[site], cut_site_d[site], window_size)

        logger.info("Alignment for {} - Done.".format(exp_type.name(), ref_df.shape[0]))
        return read_d, unaligned_reads_num

    @staticmethod
    def _add_indels_columns(reads: ReadsDf, cut_site: int, window_size: int):
        """
        Add columns with info about indels position and lengths in the qualification window.
        :param reads: The site aggregated reads
        :param cut_site: The site cut-site
        :param window_size: the window size
        :return:
        """

        start_idx = cut_site - window_size  # Start index to include indel
        end_idx = cut_site + window_size  # end index to include indel
        new_col_d = dict()
        for col in INDEL_COLS:
            new_col_d[col] = reads.shape[0] * [None]

        for row_idx, row in reads.iterrows():
            pos_idx = 0  # position index
            for length, indel_type in InputProcessing.parse_cigar(row[CIGAR]):
                # If outside the qualification window then move to the next read
                if pos_idx > end_idx:
                    break

                # For a match - continue
                if indel_type == IndelType.MATCH:
                    pos_idx += length

                # Deletions
                elif indel_type == IndelType.DEL:
                    if (pos_idx + length > start_idx) and (pos_idx < end_idx):
                        # First deletion
                        if new_col_d[DEL_LEN][row_idx] is None:
                            new_col_d[DEL_LEN][row_idx] = str(length)
                            new_col_d[DEL_START][row_idx] = str(pos_idx)
                            new_col_d[DEL_END][row_idx] = str(pos_idx + length - 1)
                        else:
                            new_col_d[DEL_LEN][row_idx] += ", {}".format(length)
                            new_col_d[DEL_START][row_idx] += ", {}".format(pos_idx)
                            new_col_d[DEL_END][row_idx] += ", {}".format(pos_idx + length - 1)

                    pos_idx += length

                # Substations
                elif indel_type == IndelType.SUB:
                    if (pos_idx + length > start_idx) and (pos_idx < end_idx):
                        # First snp
                        if new_col_d[SUB_CNT][row_idx] is None:
                            new_col_d[SUB_CNT][row_idx] = int(length)
                            new_col_d[SUB_POS][row_idx] = str(pos_idx)
                            new_col_d[SUB_POS][row_idx] += "".join([", {}".format(pos_idx+i) for i in range(1, length)])
                        else:
                            new_col_d[SUB_CNT][row_idx] += int(length)
                            new_col_d[SUB_POS][row_idx] += "".join([", {}".format(pos_idx+i) for i in range(length)])

                    pos_idx += length

                # Insertions
                elif indel_type == IndelType.INS:
                    if pos_idx >= start_idx:
                        # First Insertion
                        if new_col_d[INS_LEN][row_idx] is None:
                            new_col_d[INS_LEN][row_idx] = str(length)
                            new_col_d[INS_POS][row_idx] = str(pos_idx)
                        else:
                            new_col_d[INS_LEN][row_idx] += ", {}".format(length)
                            new_col_d[INS_POS][row_idx] += ", {}".format(pos_idx)

        # Add the cut-site
        reads[CUT_SITE] = reads.shape[0] * [cut_site]

        # Add the new
        for col in INDEL_COLS:
            reads[col] = new_col_d[col]

    @staticmethod
    def _filter_low_score_reads(reads_df: ReadsDf, output: Path, exp_type: ExpType, min_score: float, amplicon_d: Dict)\
        -> int:
        """
        - Filter low alignment score reads
        - Store all unaligned reads in a fasta format.
        :param reads_df: aligned reads df.
        :param output: output_path
        :param exp_type:
        :param min_score: user min alignment score (0-100)
        :param amplicon_d:
        :return: number of unaligned reads.
        """
        logger = Logger.get_logger()

        # Find all indexes with lower score than the threshold
        unaligned_indexes = []
        for site in amplicon_d.keys():
            site_df = reads_df.loc[reads_df[SITE_NAME] == site]
            max_score = site_df[ALIGN_SCORE].max()
            score_threshold = (min_score/100)*max_score
            unaligned_indexes += list(site_df.loc[site_df[ALIGN_SCORE] < score_threshold].index)

        total_reads_num = reads_df[FREQ].sum()
        unaligned_df = reads_df.loc[unaligned_indexes]
        unaligned_reads_num = unaligned_df[FREQ].sum()

        logger.info("Alignment for {} - {} reads weren't aligned ({:.2f}% of all reads)".format(exp_type.name(),
                    unaligned_reads_num, unaligned_reads_num/total_reads_num))

        with open(os.path.join(output, "{}_unaligned_reads.fasta".format(exp_type.name())), 'w') as file:
            for index, row in unaligned_df.iterrows():
                file.write("> unaligned read index {} - read has {} copies in the original fastq file and it's most "
                           "similar to site {}\n".format(index, row[FREQ], row[SITE_NAME]))
                file.write("{}\n".format(row[READ]))

        reads_df.drop(unaligned_df.index, inplace=True)

        return unaligned_reads_num

    @staticmethod
    def _match_by_edit_distance(read: DNASeq, amplicons: List[DNASeq], amplicon_names: List[str]) -> str:
        """
        Find for every read the alignment with min_edit_distance.
        If the min_edit_distance is below the min_distance threshold - Mark the read with this amplicon.
        :param read:  reads
        :param amplicons: List of reference sequences
        :param amplicon_names: List of site names
        :return: site_name (or NONE string)
        """
        min_name = amplicon_names[0]
        min_dist = edlib.align(read, amplicons[0])['editDistance']

        for amplicon, name in zip(amplicons[1:], amplicon_names[1:]):
            d = edlib.align(read, amplicon, k=min_dist)['editDistance']
            if d < min_dist and d != -1:
                min_dist = d
                min_name = name

        return min_name

    @classmethod
    def _align_reads_to_amplicon(cls, reads: ReadsDf, amplicons: Dict):
        """
        - Align all reads to their amplicon.
        - Compute cigar path.
        - Compute score
        :param reads: all reads
        :param amplicons: amplicons dict
        :return:
        """
        # TODO - Consider to delete one of the alignments score types.
        # Get alignment config
        cfg = Configurator.get_cfg()
        align_cfg = cfg["alignment"]

        # Init biopython aligner
        aligner = Align.PairwiseAligner()
        aligner.match_score = align_cfg["match_score"]
        aligner.mismatch_score = align_cfg["mismatch_score"]
        aligner.open_gap_score = align_cfg["open_gap_score"]
        aligner.extend_gap_score = align_cfg["extend_gap_score"]
        if align_cfg["substitution_matrix"] != "":
            if align_cfg["substitution_matrix"] in MatrixInfo.__dict__:
                aligner.substitution_matrix = MatrixInfo.__dict__[align_cfg["substitution_matrix"]]
            else:
                assert 0, "add an exception"

        new_cols_d = defaultdict(list)
        for row_idx, row in reads.iterrows():
            ref_w_ins, read_w_del, cigar, score, align_score = cls._compute_needle_wunsch_alignment(
                reference=amplicons[row[SITE_NAME]], read=row[READ], aligner=aligner)

            new_cols_d[ALIGNMENT_W_INS].append(ref_w_ins)
            new_cols_d[ALIGNMENT_W_DEL].append(read_w_del)
            new_cols_d[CIGAR].append(cigar)
            new_cols_d[SCORE].append(score)
            new_cols_d[ALIGN_SCORE].append(align_score)

        for col_name, col in new_cols_d.items():
            reads[col_name] = col

    @classmethod
    def _compute_needle_wunsch_alignment(cls, reference: DNASeq, read: DNASeq,
                                         aligner: Align.PairwiseAligner) -> Tuple[DNASeq, DNASeq, CigarPath, float]:
        """
        Compute needle wunsch alignment, cigar_path and score.
        :param read:  reads
        :param reference:  amplicon reference sequences
        :return: (alignment with ins, alignment with deletion, cigar path, score)
        """
        alignments = aligner.align(reference, read)
        [ref_with_ins, align, read_with_del, _] = format(alignments[0]).split("\n")
        cigar_path, match_bp = cls._compute_cigar_path_from_alignment(reference=ref_with_ins,
                                                                      read=read_with_del, align=align)

        # TODO - Consult with Zohar
        score = 100*(match_bp/len(reference))
        align_score = alignments[0].score

        return ref_with_ins, read_with_del, cigar_path, score, align_score

    @staticmethod
    def _compute_cigar_path_from_alignment(reference: DNASeq, read: DNASeq, align: str) -> Tuple[CigarPath, int]:
        """
        Function return cigar path from biopython alignment
        :param reference: reference with insertions
        :param read: read with deletions
        :param align: align path in biopython format
        :return: Tuple of cigar_path and number of matches.
        """

        cigar_path = []
        state: str = ""
        length = 0
        match_total = 0
        for ref_bp, align_bp, read_bp in zip(reference, align, read):
            # Insertion
            if (align_bp == "-") and (ref_bp == "-"):
                if (state != CIGAR_I) and (length != 0):
                    cigar_path.append("{}{}".format(length, state))
                    length = 1
                else:
                    length += 1
                state = CIGAR_I
            # Deletion
            elif (align_bp == "-") and (read_bp == "-"):
                if (state != CIGAR_D) and (length != 0):
                    cigar_path.append("{}{}".format(length, state))
                    length = 1
                else:
                    length += 1
                state = CIGAR_D
            # Match
            elif align_bp == "|":
                if (state != CIGAR_M) and (length != 0):
                    cigar_path.append("{}{}".format(length, state))
                    length = 1
                else:
                    length += 1
                state = CIGAR_M
                match_total += 1
            # Mismatch
            elif align_bp == "X":
                if (state != CIGAR_S) and (length != 0):
                    cigar_path.append("{}{}".format(length, state))
                    length = 1
                else:
                    length += 1
                state = CIGAR_S
            else:
                assert 0, "ERROR!!! shouldn't happen"  # TODO - change to an error

        # Push the last part of the path
        cigar_path.append("{}{}".format(length, state))

        return "".join(cigar_path), match_total

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
        for length, indel in re.findall(r'(\d+)([{}{}{}{}])'.format(CIGAR_D, CIGAR_I, CIGAR_S, CIGAR_M), cigar):
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

