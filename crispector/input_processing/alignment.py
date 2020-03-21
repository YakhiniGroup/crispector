import gzip
import os
from utils.exceptions import AlignerSubstitutionDoesntExist
from utils.constants_and_types import ReadsDf, IndelType, Path, DNASeq, CigarPath, \
    READ, ALIGNMENT_W_INS, ALIGNMENT_W_DEL, CIGAR, ALIGN_SCORE, FREQ, INS_LEN, INS_POS, DEL_LEN, DEL_START, \
    DEL_END, SUB_CNT, SUB_POS, INDEL_COLS, CIGAR_D, CIGAR_I, \
    CIGAR_S, CIGAR_M, AlignedIndel, DEL_BASE, INS_BASE, SUB_BASE, REVERSED, CIGAR_LEN, CIGAR_LEN_THRESHOLD, \
    ALIGN_CUT_SITE, ALIGNMENT_HUMAN, FILTERED_PATH, ExpType
from input_processing.utils import reverse_complement, parse_cigar
from utils.logger import LoggerWrapper
from utils.configurator import Configurator
from typing import List, Tuple, Dict
import pandas as pd
from Bio import Align
from Bio.SubsMat import MatrixInfo
from collections import defaultdict


class Alignment:
    """
    All crispector alignment functionality - Needle-Wunsch, Shifting modification and so on.
    """
    def __init__(self, align_cfg: Dict, min_score: float, min_read_length_without_primers: int,
                 window_size: int):
        self._min_score = min_score
        self._window_size = window_size

        # Set logger
        logger = LoggerWrapper.get_logger()
        self._logger = logger

        self._min_primer_dimer_thresh = min_read_length_without_primers

        # Create Aligner
        # Init biopython aligner
        self._aligner = Align.PairwiseAligner()
        self._aligner.match_score = align_cfg["match_score"]
        self._aligner.mismatch_score = align_cfg["mismatch_score"]
        self._aligner.open_gap_score = align_cfg["open_gap_score"]
        self._aligner.extend_gap_score = align_cfg["extend_gap_score"]
        if align_cfg["substitution_matrix"] != "":
            if align_cfg["substitution_matrix"] in MatrixInfo.__dict__:
                self._aligner.substitution_matrix = MatrixInfo.__dict__[align_cfg["substitution_matrix"]]
                self._logger.warning("Shifting indels to cut-site isn't available for alignment with difference score"
                                     "substitution matrix. Contact package owner if this feature is required")
            else:
                raise AlignerSubstitutionDoesntExist(align_cfg["substitution_matrix"])

    def align_reads(self, reads_df: ReadsDf, reference: DNASeq, cut_site: int, primers_len: int,
                    output: Path, exp_name: str, exp_type: ExpType) -> ReadsDf:
        """
        - Align each read to his reference and filter noisy alignments.
        - Function add columns to reads_df in place.
        :param reads_df: all reads with a left and right matches for each site
        :param reference: reference sequence to align
        :param cut_site: cut_site position
        :param primers_len: the length of both primers together
        :param output: output path for filtered reads
        :param exp_name: experiment name
        :param exp_type:
        :return: reads_df with new columns & filtered reads (ReadDf type)
        """

        if reads_df.shape[0] == 0:
            return reads_df

        # Align reads to their amplicon
        self._logger.debug("Alignment for {} - Start Needleman-Wunsch alignment for all reads.".format(exp_name))

        self._align_reads_to_amplicon(reads_df, reference)

        # Filter reads with low alignment score
        self._filter_low_score_reads(reads_df, primers_len, output, exp_name, exp_type)

        self._logger.debug("Alignment for {} - Needleman-Wunsch alignment done.".format(exp_name))

        # Shift modification into cut-site
        self._logger.debug("Alignment for {} - Start shift modifications into cut-site.".format(exp_name))
        self._shift_modifications_into_cut_site(reads_df, cut_site)
        self._logger.debug("Alignment for {} - Shift modifications into cut-site done.".format(exp_name))

        # Split read_df to all the different sites
        reads_df = reads_df.sort_values(by=[FREQ], ascending=False).reset_index(drop=True)

        # Add indels columns to reads df
        self._add_indels_columns(reads_df, cut_site)

        # Remove unnecessary columns
        reads_df.drop(columns=[REVERSED], inplace=True)
        self._logger.info("Alignment for {} - Done.".format(exp_name))

        return reads_df

    def needle_wunsch_align(self, reference: DNASeq, read: DNASeq) -> Tuple[DNASeq, DNASeq, CigarPath, int, float]:
        """
        Compute needle wunsch alignment, cigar_path and score.
        :param read:  reads
        :param reference:  amplicon reference sequences
        :param : aligner: type Align.PairwiseAligner
        :return: (alignment with ins, alignment with deletion, cigar path, score)
        """
        alignments = self._aligner.align(reference, read)
        [ref_with_ins, _, read_with_del, _] = format(alignments[0]).split("\n")
        cigar_path, cigar_len = self._compute_cigar_path_from_alignment(reference=ref_with_ins, read=read_with_del)
        align_score = alignments[0].score

        return ref_with_ins, read_with_del, cigar_path, cigar_len, align_score

    def match_by_full_alignment(self, read: DNASeq, references: List[DNASeq], names: List[str],
                                reverse_l: List[bool], ref_score_l) -> Tuple[str, bool, float]:
        """
        Find read with highest alignment score.
        Time consuming function.
        :param read: DNASeq
        :param references: list of all possible references
        :param names: reference names
        :param reverse_l: indicate if references is reversed or not
        :param ref_score_l: List of alignment scores
        :return: site_name
        """
        alignments = self._aligner.align(references[0], read)
        max_name = names[0]
        max_score = alignments[0].score / ref_score_l[0]
        max_rev = reverse_l[0]

        for reference, name, ref_score, rev in zip(references[1:], names[1:], ref_score_l[1:], reverse_l[1:]):
            alignments = self._aligner.align(reference, read)
            score = alignments[0].score / ref_score
            if score > max_score:
                max_score = score
                max_name = name
                max_rev = rev

        return max_name, max_rev, max_score

    def needle_wunsch_align_debug(self, reads: ReadsDf, reference: DNASeq):
        """
        Used only for debug filtered reads.
        - Align all reads to their amplicon.
        - Compute cigar path.
        - Compute score
        :param reads: all reads
        :param reference: reference sequence
        :return:
        """

        new_cols_d = defaultdict(list)
        for row_idx, row in reads.iterrows():
            alignments = self._aligner.align(reference, row[READ])
            [ref_with_ins, align_human, read_with_del, _] = format(alignments[0]).split("\n")
            cigar, cigar_len = self._compute_cigar_path_from_alignment(reference=ref_with_ins, read=read_with_del)
            align_score = alignments[0].score

            new_cols_d[ALIGNMENT_W_INS].append(ref_with_ins)
            new_cols_d[ALIGNMENT_W_DEL].append(read_with_del)
            new_cols_d[CIGAR_LEN].append(cigar_len)
            new_cols_d[CIGAR].append(cigar)
            new_cols_d[ALIGN_SCORE].append(align_score)
            new_cols_d[ALIGNMENT_HUMAN].append(align_human)

        for col_name, col in new_cols_d.items():
            reads[col_name] = col


    #-------------------------------#
    ######### Private methods #######
    #-------------------------------#
    def _filter_low_score_reads(self, reads_df: ReadsDf, primers_len: int, output: Path, exp_name: str,
                                exp_type: ExpType):
        """
        - Filter low alignment score reads
        - Store all unaligned reads in a fasta format.
        :param reads_df: aligned reads df.
        :param primers_len: the length of both primers together
        :param output: output path for filtered reads
        :param exp_name: experiment name
        :param exp_type:
        :return:
        """
        # Find all indexes with lower score than the threshold
        unaligned_indexes = []
        max_score = reads_df[ALIGN_SCORE].max()
        score_threshold = (self._min_score/100)*max_score

        # Filter low alignment score
        low_score_df = reads_df.loc[(reads_df[ALIGN_SCORE] < score_threshold) & (reads_df[CIGAR_LEN] > CIGAR_LEN_THRESHOLD)]

        # Filter PRIMER-DIMER affect
        min_len = primers_len + self._min_primer_dimer_thresh
        short_reads_df = reads_df.loc[reads_df[READ].str.len() < min_len]

        unaligned_indexes += list(low_score_df.index) + list(short_reads_df.index)
        total_reads_num = reads_df[FREQ].sum()
        unaligned_df = reads_df.loc[unaligned_indexes]
        unaligned_reads_num = unaligned_df[FREQ].sum()

        self._logger.info("Alignment for {} - {:,} reads were filtered out ({:.2f}% of all reads)".format(exp_name,
                          unaligned_reads_num, 100*unaligned_reads_num/total_reads_num))

        file = gzip.open(os.path.join(output, FILTERED_PATH[exp_type]), 'wt')
        for _, row in unaligned_df.iterrows():
            file.write("> filtered read with {} copies in the original fastq file.".format(row[FREQ]))
            file.write("{}\n".format(row[READ]))
        file.close()

        # remove unaligned reads from reads_df
        reads_df.drop(unaligned_df.index, inplace=True)

    def _align_reads_to_amplicon(self, reads: ReadsDf, reference: DNASeq):
        """
        - Align all reads to their amplicon.
        - Compute cigar path.
        - Compute score
        :param reads: all reads
        :param reference - reference sequence
        :return:
        """

        # If reads were demultiplexed in house - reads have an indication if they reversed or not.
        # Otherwise - compute both directions of the alignment.

        rev_flag = REVERSED in reads.columns
        rev_reference = reverse_complement(reference)

        # Reverse all the reversed reads
        if rev_flag:
            rev_reads = reads.loc[reads[REVERSED], READ]
            reads.loc[rev_reads.index, READ] = rev_reads.apply(reverse_complement)

        new_cols_d = defaultdict(list)
        for row_idx, row in reads.iterrows():
            ref_w_ins, read_w_del, cigar, c_len, score = self.needle_wunsch_align(
                reference=reference, read=row[READ])
            # compute both directions of alignment
            if not rev_flag:
                rev_ref_w_ins, rev_read_w_del, rev_cigar, rev_c_len, rev_score = self.needle_wunsch_align(
                    reference=rev_reference, read=row[READ])
                if rev_score > score:
                    ref_w_ins, read_w_del, cigar, c_len, score = rev_ref_w_ins, rev_read_w_del, rev_cigar, rev_c_len, rev_score
                    new_cols_d[REVERSED].append(True)
                else:
                    new_cols_d[REVERSED].append(False)

            new_cols_d[ALIGNMENT_W_INS].append(ref_w_ins)
            new_cols_d[ALIGNMENT_W_DEL].append(read_w_del)
            new_cols_d[CIGAR].append(cigar)
            new_cols_d[CIGAR_LEN].append(c_len)
            new_cols_d[ALIGN_SCORE].append(score)

        for col_name, col in new_cols_d.items():
            reads[col_name] = col

        # Reverse all the reversed reads if REVERSED information wasn't exist
        if not rev_flag:
            rev_reads = reads.loc[reads[REVERSED], READ]
            reads.loc[rev_reads.index, READ] = rev_reads.apply(reverse_complement)

    @staticmethod
    def _compute_cigar_path_from_alignment(reference: DNASeq, read: DNASeq) -> Tuple[CigarPath, int]:
        """
        Function return cigar path from biopython alignment
        :param reference: reference with insertions
        :param read: read with deletions
        :return: cigar_path, cigar len
        """
        cigar_path = []
        state: str = ""
        length = 0
        cigar_length = 0
        for ref_bp, read_bp in zip(reference, read):
            # Insertion
            if ref_bp == "-":
                if (state != CIGAR_I) and (length != 0):
                    cigar_path.append("{}{}".format(length, state))
                    length = 1
                    cigar_length += 1
                else:
                    length += 1
                state = CIGAR_I
            # Deletion
            elif read_bp == "-":
                if (state != CIGAR_D) and (length != 0):
                    cigar_path.append("{}{}".format(length, state))
                    length = 1
                    cigar_length += 1
                else:
                    length += 1
                state = CIGAR_D
            # Match
            elif ref_bp == read_bp:
                if (state != CIGAR_M) and (length != 0):
                    cigar_path.append("{}{}".format(length, state))
                    length = 1
                    cigar_length += 1
                else:
                    length += 1
                state = CIGAR_M
            # Mismatch
            else:
                if (state != CIGAR_S) and (length != 0):
                    cigar_path.append("{}{}".format(length, state))
                    length = 1
                    cigar_length += 1
                else:
                    length += 1
                state = CIGAR_S

        # Push the last part of the path
        cigar_path.append("{}{}".format(length, state))
        return "".join(cigar_path), cigar_length

    def _add_indels_columns(self, reads: ReadsDf, cut_site: int):
        """
        Add columns with info about indels position and lengths in the qualification window.
        :param reads: The site aggregated reads
        :param cut_site: The site cut-site
        :return:
        """
        start_idx = cut_site - self._window_size  # Start index to include indel
        end_idx = cut_site + self._window_size  # end index to include indel
        new_col_d = dict()

        for col in INDEL_COLS:
            new_col_d[col] = reads.shape[0] * [None]

        for row_idx, row in reads.iterrows():
            pos_idx = 0  # position index for the original reference (with no indels)
            align_idx = 0  # position index for the alignment (with indels)
            reference = row[ALIGNMENT_W_INS]
            read = row[ALIGNMENT_W_DEL]
            for length, indel_type in parse_cigar(row[CIGAR]):
                # If outside the qualification window then move to the next read
                if pos_idx > end_idx:
                    break

                # Deletions
                elif indel_type == IndelType.DEL:
                    if (pos_idx + length > start_idx) and (pos_idx < end_idx):
                        # First deletion
                        if new_col_d[DEL_LEN][row_idx] is None:
                            new_col_d[DEL_LEN][row_idx] = str(length)
                            new_col_d[DEL_START][row_idx] = str(pos_idx)
                            new_col_d[DEL_END][row_idx] = str(pos_idx + length - 1)
                            new_col_d[DEL_BASE][row_idx] = reference[align_idx:align_idx+length]
                        else:
                            new_col_d[DEL_LEN][row_idx] += ", {}".format(length)
                            new_col_d[DEL_START][row_idx] += ", {}".format(pos_idx)
                            new_col_d[DEL_END][row_idx] += ", {}".format(pos_idx + length - 1)
                            new_col_d[DEL_BASE][row_idx] += ", {}".format(reference[align_idx:align_idx+length])

                # Substations
                elif indel_type == IndelType.SUB:
                    if (pos_idx + length > start_idx) and (pos_idx < end_idx):
                        # First snp
                        if new_col_d[SUB_CNT][row_idx] is None:
                            new_col_d[SUB_CNT][row_idx] = int(length)
                            new_col_d[SUB_POS][row_idx] = str(pos_idx)
                            new_col_d[SUB_POS][row_idx] += "".join([", {}".format(pos_idx+i) for i in range(1, length)])
                            new_col_d[SUB_BASE][row_idx] = read[align_idx]
                            new_col_d[SUB_BASE][row_idx] += "".join([", {}".format(read[align_idx+i]) for i in range(1, length)])
                        else:
                            new_col_d[SUB_CNT][row_idx] += int(length)
                            new_col_d[SUB_POS][row_idx] += "".join([", {}".format(pos_idx+i) for i in range(length)])
                            new_col_d[SUB_BASE][row_idx] += "".join([", {}".format(read[align_idx+i]) for i in range(1, length)])

                # Insertions
                elif indel_type == IndelType.INS:
                    if pos_idx >= start_idx:
                        # First Insertion
                        if new_col_d[INS_LEN][row_idx] is None:
                            new_col_d[INS_LEN][row_idx] = str(length)
                            new_col_d[INS_POS][row_idx] = str(pos_idx)
                            new_col_d[INS_BASE][row_idx] = read[align_idx:align_idx+length]
                        else:
                            new_col_d[INS_LEN][row_idx] += ", {}".format(length)
                            new_col_d[INS_POS][row_idx] += ", {}".format(pos_idx)
                            new_col_d[INS_BASE][row_idx] += ", {}".format(read[align_idx:align_idx+length])

                # update indexes and store aligned_cut_site
                if indel_type != IndelType.INS:
                    # store cut-site position
                    if cut_site in range(pos_idx, pos_idx + length + 1):
                        new_col_d[ALIGN_CUT_SITE][row_idx] = align_idx + (cut_site - pos_idx)
                    pos_idx += length
                align_idx += length

        # Add the new columns
        for col in INDEL_COLS:
            reads[col] = new_col_d[col]

    @staticmethod
    def compute_alignment_score_from_cigar(cigar):
        """
        Debug function. compute alignment score from cigar path.
        :param cigar:
        :return:
        """
        # Get alignment config
        cfg = Configurator.get_cfg()
        align_cfg = cfg["alignment"]
        match_score = align_cfg["match_score"]
        mismatch_score = align_cfg["mismatch_score"]
        open_gap_score = align_cfg["open_gap_score"]
        extend_gap_score = align_cfg["extend_gap_score"]
        score = 0

        # Find closest indels from left and from right
        for length, indel_type in parse_cigar(cigar):
            if indel_type == IndelType.MATCH:
                score += length * match_score
            elif indel_type in [IndelType.INS, IndelType.DEL]:
                score += open_gap_score + (length - 1) * extend_gap_score
            else:
                score += length * mismatch_score

        return score

    ######### Indel Modification #########
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

    @staticmethod
    def _find_closest_indels_to_cut_site(cigar: CigarPath, cut_site: int) -> Tuple[AlignedIndel, AlignedIndel, int]:
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
        for length, indel_type in parse_cigar(cigar):

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

    def _shift_modifications_into_cut_site(self, reads: ReadsDf, cut_site: int):
        """
        Shift deletions and insertions with a region of 2*window_size into the cut-site
        :param reads: ReadsDf - All reads, already aligned with a cigar path.
        :param cut_site: cure_site position
        :return: no return value. Change reads inplace.
        """
        update_idx = []  # Changed read indexes list
        reference_l = []  # Changed references (read with insertions) list
        read_l = []  # Changed reads (read with deletions) list
        cigar_l = []  # Changed cigar path list

        for row_idx, row in reads.iterrows():
            cigar = row[CIGAR]
            reference = row[ALIGNMENT_W_INS]
            read = row[ALIGNMENT_W_DEL]
            changed_right = False
            changed_left = False

            # Find closest indels left and right to the cut-site and cut-site in alignment coordinates
            most_left, most_right, aligned_cut_site = self._find_closest_indels_to_cut_site(cigar, cut_site)

            # Shift most left modification to the cut-site
            if most_left is not None:
                indel_type, length, align_idx = most_left  # most_left is type AlignedIndel
                # shift indels with a region of twice the window size from the cut-site
                if align_idx + length + (2*self._window_size) < aligned_cut_site:
                    changed_left = False
                elif indel_type == IndelType.DEL:
                    read, changed_left = self._shift_indel_from_left(read, reference, length, align_idx,
                                                                    aligned_cut_site)
                else:
                    reference, changed_left = self._shift_indel_from_left(reference, read, length, align_idx,
                                                                         aligned_cut_site)

            # Shift most right modification to the cut-site
            if most_right is not None:
                indel_type, length, align_idx = most_right  # most_left is type AlignedIndel
                # shift indels with a region of twice the window size from the cut-site
                if align_idx > aligned_cut_site + (2*self._window_size):
                    changed_right = False
                elif indel_type == IndelType.DEL:
                    read, changed_right = self._shift_indel_from_right(read, reference, length, align_idx,
                                                                      aligned_cut_site)
                else:
                    reference, changed_right = self._shift_indel_from_right(reference, read, length, align_idx,
                                                                           aligned_cut_site)
            # Mark read if it was changed
            if changed_right or changed_left:
                # Compute new cigar_path
                cigar, _ = self._compute_cigar_path_from_alignment(reference, read)
                update_idx.append(row_idx)
                reference_l.append(reference)
                read_l.append(read)
                cigar_l.append(cigar)

        # Update with all changed reads
        updated_reads_df = pd.DataFrame({ALIGNMENT_W_INS: reference_l, ALIGNMENT_W_DEL: read_l,
                                         CIGAR: cigar_l}, index=update_idx)
        reads.update(updated_reads_df)
