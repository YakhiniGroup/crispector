import os
import subprocess
from exceptions import FastpRunTimeError, SgRNANotInReferenceSequence, CantOpenMergedFastqFile, UnknownAlignmentChar
from constants_and_types import AmpliconDf, ReadsDict, ExpType, ReadsDf, IndelType, Path, DNASeq, CigarPath, FASTP_DIR, \
    READ, ALIGNMENT_W_INS, ALIGNMENT_W_DEL, CIGAR, SCORE, ALIGN_SCORE, FREQ, INS_LEN, INS_POS, DEL_LEN, DEL_START, \
    DEL_END, SUB_CNT, SUB_POS, INDEL_COLS, COMPLEMENT, REFERENCE, SGRNA, SITE_NAME, CUT_SITE, CIGAR_D, CIGAR_I, \
    CIGAR_S, CIGAR_M, AlignedIndel
from utils import Logger, Configurator
import edlib      # TODO - add to install requirements, conda install python-edlib, consider conda install gcc before.
# TODO - no reason to install edlib with pip...
from typing import List, Tuple, Dict
import re
import pandas as pd  # TODO - add to install requirements
from Bio import Align  # TODO - add to install requirements, make sure conda install biopython=1.74!
from Bio.SubsMat import MatrixInfo
from collections import defaultdict
import json


class InputProcessing:
    """
    A container static class with all relevant function to process crispector input.
    """

    @staticmethod
    def fastp(in1: Path, in2: Path, fastp_options_string: Path, output_dir: Path, fastq_type: ExpType) -> \
        Tuple[Path, int]:
        """
        Wrapper for fastp SW.
        :param in1: read1 input
        :param in2: read2 input
        :param fastp_options_string: any additional fastp options
        :param output_dir: output path
        :param fastq_type: ExpType
        :return: path to merged fastq file and number of reads in the input
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

        # Get the number of reads in the input
        fastp_summary_path = os.path.join(fastp_output, "fastp.json")
        if os.path.isfile(fastp_summary_path):
            with open(fastp_summary_path) as json_file:
                summary = json.load(json_file)
                reads_in_input_num = summary["summary"]["before_filtering"]["total_reads"] // 2
        else:
            reads_in_input_num = -1

        logger.info("fastp for {} - Done.".format(fastq_type.name()))

        return merged_path, reads_in_input_num

    @classmethod
    def split_read_and_align(cls, merged_fastq: Path, ref_df: AmpliconDf, min_score: float, output: Path,
                             exp_type: ExpType, override_alignment: bool) -> Tuple[ReadsDict, int, int]:
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
        :return: A tuple of (dictionary with key=SITE_NAME and value=ReadsDf, merged reads num, aligned reads num)
        """

        logger = Logger.get_logger()
        amplicon_d = {k: v for (k, v) in zip(ref_df[SITE_NAME].values, ref_df[REFERENCE].values)}
        cut_site_d = {k: v for (k, v) in zip(ref_df[SITE_NAME].values, ref_df[CUT_SITE].values)}
        cfg = Configurator.get_cfg()
        window_size = cfg["window_size"]

        # TODO - remove when delete override_alignment
        if override_alignment:
            logger.debug("Alignment for {} - Skip alignment and read from input.".format(exp_type.name()))
            reads_df = pd.read_csv(merged_fastq)
            merged_reads_num, aligned_reads_num = reads_df[FREQ].sum(), reads_df[FREQ].sum()
        else:
            # Read file and convert to pandas
            logger.debug("Alignment for {} - Read merged fastq from {}".format(exp_type.name(), merged_fastq))
            reads = cls._parse_fastq_file(merged_fastq)
            reads_df: ReadsDf = pd.DataFrame(data=reads, columns=[READ])

            # Group identical reads together
            reads_df = reads_df.groupby(READ).size().to_frame(FREQ).reset_index()
            merged_reads_num = reads_df[FREQ].sum()

            # Match reads and amplicons
            logger.info("Alignment for {} - Start alignment for {} sites.".format(exp_type.name(), ref_df.shape[0]))
            logger.debug("Alignment for {} - Demultiplexing reads for {:,} sites.".format(exp_type.name(),
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
            cls._filter_low_score_reads(reads_df, output, exp_type, min_score, amplicon_d)
            aligned_reads_num = reads_df[FREQ].sum()
            logger.debug("Alignment for {} - Needleman-Wunsch alignment done.".format(exp_type.name()))

            # Shift modification into cut-site
            logger.debug("Alignment for {} - Start shift modifications into cut-site.".format(exp_type.name()))
            cls._shift_modifications_into_cut_site(reads_df, cut_site_d, window_size)
            logger.debug("Alignment for {} - Shift modifications into cut-site done.".format(exp_type.name()))

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
        return read_d, merged_reads_num, aligned_reads_num

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
    def _filter_low_score_reads(reads_df: ReadsDf, output: Path, exp_type: ExpType, min_score: float, amplicon_d: Dict):
        """
        - Filter low alignment score reads
        - Store all unaligned reads in a fasta format.
        :param reads_df: aligned reads df.
        :param output: output_path
        :param exp_type:
        :param min_score: user min alignment score (0-100)
        :param amplicon_d:
        :return:
        """
        logger = Logger.get_logger()

        # Find all indexes with lower score than the threshold
        unaligned_indexes = []
        for site in amplicon_d.keys():
            site_df = reads_df.loc[reads_df[SITE_NAME] == site]
            max_score = site_df[ALIGN_SCORE].max() # TODO - decide if to use alignment score or regular score, delete the other
            score_threshold = (min_score/100)*max_score
            unaligned_indexes += list(site_df.loc[site_df[ALIGN_SCORE] < score_threshold].index)

        total_reads_num = reads_df[FREQ].sum()
        unaligned_df = reads_df.loc[unaligned_indexes]
        unaligned_reads_num = unaligned_df[FREQ].sum()

        logger.info("Alignment for {} - {:,} reads weren't aligned ({:.2f}% of all reads)".format(exp_type.name(),
                    unaligned_reads_num, unaligned_reads_num/total_reads_num))

        with open(os.path.join(output, "{}_unaligned_reads.fasta".format(exp_type.name())), 'w') as file:
            for index, row in unaligned_df.iterrows():
                file.write("> unaligned read index {} - read has {} copies in the original fastq file and it's most "
                           "similar to site {}\n".format(index, row[FREQ], row[SITE_NAME]))
                file.write("{}\n".format(row[READ]))

        reads_df.drop(unaligned_df.index, inplace=True)

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
    def _compute_needle_wunsch_alignment(cls, reference: DNASeq, read: DNASeq, aligner: Align.PairwiseAligner) \
        -> Tuple[DNASeq, DNASeq, CigarPath, float, float]:
        """
        Compute needle wunsch alignment, cigar_path and score.
        :param read:  reads
        :param reference:  amplicon reference sequences
        :param : aligner: type Align.PairwiseAligner
        :return: (alignment with ins, alignment with deletion, cigar path, score)
        """
        alignments = aligner.align(reference, read)
        [ref_with_ins, align, read_with_del, _] = format(alignments[0]).split("\n")
        cigar_path = cls._compute_cigar_path_from_alignment(reference=ref_with_ins, read=read_with_del)

        # TODO - Consult with Zohar and delete one score type
        match_bp  = len(list(filter(cls.delete_non_match, align)))  # Can be taken from _compute_cigar_path_from_alignment(), but probably we don't need it
        score = 100*(match_bp/len(reference))
        align_score = alignments[0].score

        return ref_with_ins, read_with_del, cigar_path, score, align_score

    # TODO - delete this function!!!!
    @staticmethod
    def delete_non_match(base):
        return base == "|"

    @staticmethod
    def _compute_cigar_path_from_alignment(reference: DNASeq, read: DNASeq) -> CigarPath:
        """
        Function return cigar path from biopython alignment
        :param reference: reference with insertions
        :param read: read with deletions
        :param align: align path in biopython format
        :return: cigar_path
        """
        # TODO - Change to work without align
        cigar_path = []
        state: str = ""
        length = 0
        for ref_bp, read_bp in zip(reference, read):
            # Insertion
            if ref_bp == "-":
                if (state != CIGAR_I) and (length != 0):
                    cigar_path.append("{}{}".format(length, state))
                    length = 1
                else:
                    length += 1
                state = CIGAR_I
            # Deletion
            elif read_bp == "-":
                if (state != CIGAR_D) and (length != 0):
                    cigar_path.append("{}{}".format(length, state))
                    length = 1
                else:
                    length += 1
                state = CIGAR_D
            # Match
            elif ref_bp == read_bp:
                if (state != CIGAR_M) and (length != 0):
                    cigar_path.append("{}{}".format(length, state))
                    length = 1
                else:
                    length += 1
                state = CIGAR_M
            # Mismatch
            else:
                if (state != CIGAR_S) and (length != 0):
                    cigar_path.append("{}{}".format(length, state))
                    length = 1
                else:
                    length += 1
                state = CIGAR_S

        # Push the last part of the path
        cigar_path.append("{}{}".format(length, state))

        return "".join(cigar_path)

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
    def _shift_indel_from_left(read_a: DNASeq, read_b: DNASeq, length: int, align_idx: int, alignment_cut_site: int) \
        -> Tuple[DNASeq, bool]:
        """
        Shift indel from left to right if possible. only read_a is changed.
        For insertion - read_a is the reference and read_b is the read.
        For deletion - read_a is the read and read_b is the reference.
        Shift is possible in two scenarios:
        1. The "pushed" base in read_a is identical to his match in read_b (read_a[end_idx] == read_b[start_idx])
        2. The "pushed" base in read_a is already marked as a mismatch
        ((read_a[end_idx] != read_b[end_idx]) and (read_b[end_idx] != "-"))
        :param read_a: DNASeq
        :param read_b: DNASeq
        :param length: indel length
        :param align_idx: indel position
        :param alignment_cut_site: cut-site position
        :return: A tuple of shifted read_a and flag is read_a was changes
        """

        shift = 0
        start_idx = align_idx
        end_idx = align_idx + length  # Consider to replace the base at end_idx
        while end_idx < alignment_cut_site:
            if (read_a[end_idx] == read_b[start_idx]) or \
               ((read_a[end_idx] != read_b[end_idx]) and (read_b[end_idx] != "-")):
                shift += 1
            # Can't shift any further
            else:
                break
            start_idx += 1
            end_idx += 1

        if shift == 0:
            return read_a, False

        # Shift indels left
        shifted_bases = read_a[align_idx + length:align_idx + length + shift]
        indels = read_a[align_idx:align_idx + length]
        read_a = read_a[0:align_idx] + shifted_bases + indels + read_a[align_idx + length + shift:]

        return read_a, True

    @staticmethod
    def _shift_indel_from_right(read_a: DNASeq, read_b: DNASeq, length: int, align_idx: int, alignment_cut_site: int) \
        -> Tuple[DNASeq, bool]:
        """
        Shift indel from right to left if possible. only read_a is changed.
        For insertion - read_a is the reference and read_b is the read.
        For deletion - read_a is the read and read_b is the reference.
        Shift is possible in two scenarios:
        1. The "pushed" base in read_a is identical to his match in read_b (read_a[end_idx] == read_b[start_idx])
        2. The "pushed" base in read_a is already marked as a mismatch
        ((read_a[end_idx] != read_b[end_idx]) and (read_b[end_idx] != "-"))
        :param read_a: DNASeq
        :param read_b: DNASeq
        :param length: indel length
        :param align_idx: indel position
        :param alignment_cut_site: cut-site position
        :return: A tuple of shifted read_a and flag is read_a was changes
        """
        shift = 0
        start_idx = align_idx - 1  # Consider to replace base at start_idx
        end_idx = align_idx + length - 1
        while start_idx >= alignment_cut_site:
            if (read_a[start_idx] == read_b[end_idx]) or \
               ((read_a[start_idx] != read_b[start_idx]) and (read_b[start_idx] != "-")):
                shift += 1
            # Can't shift any further
            else:
                break
            start_idx -= 1
            end_idx -= 1

        if shift == 0:
            return read_a, False

        # Shift bp's in the read for deletions
        shifted_bases = read_a[align_idx - shift:align_idx]
        indels = read_a[align_idx:align_idx + length]
        read_a = read_a[0:align_idx - shift] + indels + shifted_bases + read_a[align_idx + length:]

        return read_a, True

    # TODO - keep only for debug?
    @classmethod
    def compute_alignment_score_from_cigar(cls, cigar):
        # Get alignment config
        cfg = Configurator.get_cfg()
        align_cfg = cfg["alignment"]
        match_score = align_cfg["match_score"]
        mismatch_score = align_cfg["mismatch_score"]
        open_gap_score = align_cfg["open_gap_score"]
        extend_gap_score = align_cfg["extend_gap_score"]
        score = 0

        # Find closest indels from left and from right
        for length, indel_type in cls.parse_cigar(cigar):
            if indel_type == IndelType.MATCH:
                score += length * match_score
            elif indel_type in [IndelType.INS, IndelType.DEL]:
                score += open_gap_score + (length - 1) * extend_gap_score
            else:
                score += length * mismatch_score

        return score

    @classmethod
    def _find_closest_indels_to_cut_site(cls, cigar: CigarPath, cut_site: int) -> Tuple[AlignedIndel, AlignedIndel, int]:
        """
        Find closest indels left to the cut-site and right to the cut-site.
        :param cigar: cigar path of the alignment
        :param cut_site: cut-site in reference coordinates
        :return: most_left indel, most right indel and aligned cut-site. return None if indel is already on cut-site
        or if there is no such indel.
        """
        most_left: AlignedIndel = None  # closest indel from left to the cut-site
        most_right: AlignedIndel = None  # closest indel from right to the cut-site
        pos_idx = 0  # position index for the original reference (with no indels)
        align_idx = 0  # position index for the alignment (with indels)
        aligned_cut_site = -1  # cut-site in align_idx coordinates
        search_finished = False  # Search finished flag. will be used to jump to return value

        # Find closest indels from left and from right
        for length, indel_type in cls.parse_cigar(cigar):

            # Deletions
            if indel_type == IndelType.DEL:
                # Deletion is left to the cut-site
                if pos_idx <= cut_site:
                    if (pos_idx + length) < cut_site:
                        most_left = (indel_type, length, align_idx)
                    # Deletion is already on the cut-site
                    else:
                        most_left = None
                        search_finished = True
                # Deletion is right to the cut-site
                elif most_right is None:
                    most_right = (indel_type, length, align_idx)
                    search_finished = True

            # Insertions
            elif indel_type == IndelType.INS:
                # Insertion is already on the cut-site
                if pos_idx == cut_site:
                    most_left = None
                    search_finished = True
                # Insertion is left to the cut-site
                elif pos_idx < cut_site:
                    most_left = (indel_type, length, align_idx)
                # Insertion is right to the cut-site
                elif most_right is None:
                    most_right = (indel_type, length, align_idx)
                    search_finished = True

            # Update indexes and compute cut-site in alignment coordinates
            if indel_type != IndelType.INS:
                # store cut-site position
                if cut_site in range(pos_idx, pos_idx + length + 1):
                    aligned_cut_site = align_idx + (cut_site - pos_idx)
                pos_idx += length
            align_idx += length

            if search_finished:
                break

        return most_left, most_right, aligned_cut_site

    @classmethod
    def _shift_modifications_into_cut_site(cls, reads: ReadsDf, cut_site_d: Dict[str, int], window_size: int):
        """
        Shift deletions and insertions with a region of 2*window_size into the cut-site
        :param reads: ReadsDf - All reads, already aligned with a cigar path.
        :param cut_site_d: cut-site dict
        :param window_size: config window_size
        :return: no return value. Change reads inplace.
        """
        update_idx = []  # Changed read indexes list
        reference_l = []  # Changed references (read with insertions) list
        read_l = []  # Changed reads (read with deletions) list
        cigar_l = []  # Changed cigar path list

        for row_idx, row in reads.iterrows():
            cut_site = cut_site_d[row[SITE_NAME]]  # cut-site in pos_idx coordinates
            cigar = row[CIGAR]
            reference = row[ALIGNMENT_W_INS]
            read = row[ALIGNMENT_W_DEL]
            changed_right = False
            changed_left = False

            # Find closest indels left and right to the cut-site and cut-site in alignment coordinates
            most_left, most_right, aligned_cut_site = cls._find_closest_indels_to_cut_site(cigar, cut_site)

            # Shift most left modification to the cut-site
            if most_left is not None:
                indel_type, length, align_idx = most_left  # most_left is type AlignedIndel
                # shift indels with a region of twice the window size from the cut-site
                if align_idx + length + (2*window_size) < aligned_cut_site:
                    changed_left = False
                elif indel_type == IndelType.DEL:
                    read, changed_left = cls._shift_indel_from_left(read, reference, length, align_idx,
                                                                    aligned_cut_site)
                else:
                    reference, changed_left = cls._shift_indel_from_left(reference, read, length, align_idx,
                                                                         aligned_cut_site)

            # Shift most right modification to the cut-site
            if most_right is not None:
                indel_type, length, align_idx = most_right  # most_left is type AlignedIndel
                # shift indels with a region of twice the window size from the cut-site
                if align_idx > aligned_cut_site + (2*window_size):
                    changed_right = False
                elif indel_type == IndelType.DEL:
                    read, changed_right = cls._shift_indel_from_right(read, reference, length, align_idx,
                                                                      aligned_cut_site)
                else:
                    reference, changed_right = cls._shift_indel_from_right(reference, read, length, align_idx,
                                                                           aligned_cut_site)
            # Mark read if it was changed
            if changed_right or changed_left:
                # Compute new cigar_path
                cigar = cls._compute_cigar_path_from_alignment(reference, read)
                # TODO - delete this assertion from final version
                original_score = cls.compute_alignment_score_from_cigar(row[CIGAR])
                new_score = cls.compute_alignment_score_from_cigar(cigar)
                assert new_score == original_score, "Bad modification shift! reference={}, read={}, cut-site={}".format(
                    row[ALIGNMENT_W_INS], row[ALIGNMENT_W_DEL], aligned_cut_site)

                update_idx.append(row_idx)
                reference_l.append(reference)
                read_l.append(read)
                cigar_l.append(cigar)

        # Update with all changed reads
        updated_reads_df = pd.DataFrame({ALIGNMENT_W_INS: reference_l, ALIGNMENT_W_DEL: read_l,
                                         CIGAR: cigar_l}, index=update_idx)
        reads.update(updated_reads_df)
