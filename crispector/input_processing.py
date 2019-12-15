import os
import subprocess
from exceptions import FastpRunTimeError, SgRNANotInReferenceSequence, CantOpenMergedFastqFile, \
    AlignerSubstitutionDoesntExist
from constants_and_types import AmpliconDf, ReadsDict, ExpType, ReadsDf, IndelType, Path, DNASeq, CigarPath, FASTP_DIR, \
    READ, ALIGNMENT_W_INS, ALIGNMENT_W_DEL, CIGAR, ALIGN_SCORE, FREQ, INS_LEN, INS_POS, DEL_LEN, DEL_START, \
    DEL_END, SUB_CNT, SUB_POS, INDEL_COLS, COMPLEMENT, REFERENCE, SGRNA, SITE_NAME, CUT_SITE, CIGAR_D, CIGAR_I, \
    CIGAR_S, CIGAR_M, AlignedIndel, DEL_BASE, INS_BASE, SUB_BASE, REVERSED, L_SITE, L_REV, R_SITE, R_REV, \
    L_READ, R_READ, PRIMER_LEN, TransDf, TRANS_NAME, BAD_AMPLICON_THRESHOLD, CIGAR_LEN, CIGAR_LEN_THRESHOLD, MAX_SCORE, \
    F_PRIMER, R_PRIMER, MIN_PRIMER_DIMER_THRESH, SGRNA_REVERSED, CS_SHIFT_L, CS_SHIFT_L, CS_SHIFT_R
from utils import Logger, Configurator
from typing import List, Tuple, Dict
import re
import pandas as pd
from Bio import Align
from Bio.SubsMat import MatrixInfo
from collections import defaultdict
import json
import edlib

class InputProcessing:
    """
    A helper class with all relevant functions to process crispector input.
    """
    def __init__(self, ref_df: AmpliconDf, output: Path, min_alignment_score: float, min_read_length: int,
                 max_error_on_primer: int):
        """
        :param ref_df: AmpliconDf type
        :param output: output path
        :param min_alignment_score: user min alignment score (0-100)
        :param min_read_length:  minimum read length
        :return:
        """
        self._ref_df = ref_df
        self._output = output
        self._min_score = min_alignment_score
        self._min_read_length = min_read_length
        self._max_error_on_primer = max_error_on_primer
        # Set logger
        logger = Logger.get_logger()
        self._logger = logger

        # Get config
        self._cfg = Configurator.get_cfg()

        # Create Aligner
        align_cfg = self._cfg["alignment"]

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

        # Add max_score column to ref_df
        max_score_list = []
        for _, row in self._ref_df.iterrows():
            reference = row[REFERENCE]
            _, _, _, _, max_score = self._compute_needle_wunsch_alignment(reference=reference, read=reference)
            max_score_list.append(max_score)
        self._ref_df[MAX_SCORE] = max_score_list

        # Add primers where values are None
        self._ref_df.loc[self._ref_df[F_PRIMER].isna(), F_PRIMER] = self._ref_df.loc[self._ref_df[F_PRIMER].isna(),
                                                                                     REFERENCE].str[0:PRIMER_LEN]

        self._ref_df.loc[self._ref_df[R_PRIMER].isna(), R_PRIMER] = self._ref_df.loc[self._ref_df[R_PRIMER].isna(),
                                                                                     REFERENCE].apply(self.reverse_complement).str[0:PRIMER_LEN]

        # Reads numbers for statistics
        self._input_n = dict()
        self._input_n[ExpType.TX] = -1
        self._input_n[ExpType.MOCK] = -1

        self._merged_n = dict()
        self._merged_n[ExpType.TX] = -1
        self._merged_n[ExpType.MOCK] = -1

        self._aligned_n = dict()
        self._aligned_n[ExpType.TX] = -1
        self._aligned_n[ExpType.MOCK] = -1

    # -------------------------------#
    ######### Public methods #########
    # -------------------------------#
    def convert_sgRNA_to_cut_site_position(self, cut_site_position: int):
        """
        :param cut_site_position: position relative to the PAM
        :return:
        """
        self._ref_df[[CUT_SITE, SGRNA_REVERSED]] = self._ref_df.apply(lambda row: self._get_expected_cut_site(row[REFERENCE],
                                                   row[SGRNA], cut_site_position, row[SITE_NAME]), axis=1, result_type='expand')


    def detect_ambiguous_cut_site(self, cut_site_pos: int, ambiguous_cut_site_detection: bool):
        """
        # Check if there is an alternative PAM
        :param cut_site_pos: position relative to the PAM
        :param ambiguous_cut_site_detection: Flag
        :return:
        """
        self._ref_df[CS_SHIFT_L] = False
        self._ref_df[CS_SHIFT_R] = False

        if not ambiguous_cut_site_detection:
            return

        for idx, row in self._ref_df.iterrows():
            reverse = row[SGRNA_REVERSED]
            ref = self.reverse_complement(row[REFERENCE]) if reverse else row[REFERENCE]
            cut_site = len(row[REFERENCE]) - row[CUT_SITE] if reverse else row[CUT_SITE]
            PAM_start = cut_site-cut_site_pos

            # cut-site can possibly shift to the right
            if ref[PAM_start+2:PAM_start+4] == "GG":
                if not reverse:
                    self._ref_df.at[idx, CS_SHIFT_R] = True
                else: # right is left for 3'->5' reference
                    self._ref_df.at[idx, CS_SHIFT_L] = True

            # cut-site can possibly shift to the left
            if ref[PAM_start:PAM_start+2] == "GG":
                if not reverse:
                    self._ref_df.at[idx, CS_SHIFT_L] = True
                else:  # left is right for 3'->5' reference
                    self._ref_df.at[idx, CS_SHIFT_R] = True

    def fastp(self, in1: Path, in2: Path, options_string: str, threads: int, exp_type: ExpType) -> Path:
        """
        Wrapper for fastp SW.
        :param in1: read1 input
        :param in2: read2 input
        :param options_string: any additional fastp options
        :param threads: number of threads for fastp
        :param exp_type: ExpType
        :return: path to merged fastq file
        """
        # Create output folder
        fastp_output = os.path.join(self._output, FASTP_DIR[exp_type])
        if not os.path.exists(fastp_output):
            os.makedirs(fastp_output)

        merged_path = os.path.join(fastp_output, "merged_reads.fastq")

        command = ["fastp", "-i", in1, "-I", in2, "-o", os.path.join(fastp_output, "r1_filtered_reads.fastq"),
                   "-O", os.path.join(fastp_output, "r2_filtered_reads.fastq"), "-m", "--merged_out", merged_path,
                   "-j", os.path.join(fastp_output, "fastp.json"), "-h", os.path.join(fastp_output, "fastp.html"),
                   "-w {}".format(threads), "--length_required {}".format(self._min_read_length),
                   options_string, ">> {} 2>&1".format(Logger.get_log_path())]

        command = " ".join(command)

        self._logger.debug("fastp for {} - Command {}".format(exp_type.name(), command))
        self._logger.info("fastp for {} - Run (may take a few minutes).".format(exp_type.name()))
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
                merged_reads_num = summary["summary"]["after_filtering"]["total_reads"]
        else:
            reads_in_input_num = -1
            merged_reads_num =-1

        self._logger.info("fastp for {} - Done.".format(exp_type.name()))
        self._input_n[exp_type] = reads_in_input_num
        self._merged_n[exp_type] = merged_reads_num

        return merged_path

    def demultiplex_reads(self, merged_fastq: Path, exp_type: ExpType) -> ReadsDf:
        """
        Demultiplex reads using edit distance on primers.
        Assign each read with a left primer match and right primer match
        :param merged_fastq: fastp merged fastq file
        :param exp_type: ExpType
        :return: ReadsDf
        """
        reads = self._parse_fastq_file(merged_fastq)
        reads_df = pd.DataFrame(data=reads, columns=[READ])

        # Prepare primers for match
        references = self._ref_df[REFERENCE].str
        left_primers = list(references[:PRIMER_LEN]) + list(references[-PRIMER_LEN:].apply(self.reverse_complement))
        right_primers = list(references[-PRIMER_LEN:]) + list(references[:PRIMER_LEN].apply(self.reverse_complement))
        primers_names = 2 * list(self._ref_df[SITE_NAME])
        primers_rev = self._ref_df.shape[0] * [False] + self._ref_df.shape[0] * [True]

        # Group identical reads together
        reads_df = reads_df.groupby(READ).size().to_frame(FREQ).reset_index()

        # Create temporary columns for read start and end
        reads_df[L_READ] = reads_df[READ].str[0:PRIMER_LEN]
        reads_df[R_READ] = reads_df[READ].str[-PRIMER_LEN:]

        self._logger.info("Demultiplexing for {} - Start Demultiplexing {:,} reads - May take a few minutes"
                    .format(exp_type.name(), reads_df[FREQ].sum()))

        # Find a match for left and right parts of the read
        l_match = self._compute_read_primer_matching(reads_df[L_READ].unique(), left_primers, primers_rev,
                                                     primers_names, self._max_error_on_primer)
        r_match = self._compute_read_primer_matching(reads_df[R_READ].unique(), right_primers, primers_rev,
                                                     primers_names, self._max_error_on_primer)

        reads_df[[L_SITE, L_REV]] = reads_df.apply((lambda row: l_match[row[L_READ]]), axis=1, result_type='expand')
        reads_df[[R_SITE, R_REV]] = reads_df.apply((lambda row: r_match[row[R_READ]]), axis=1, result_type='expand')

        self._logger.info("Demultiplexing for {} - Done".format(exp_type.name()))

        return reads_df

    def align_reads(self, reads_df: ReadsDf, exp_type: ExpType, allow_trans: bool, debug: bool) -> Tuple[ReadsDict, TransDf]:
        """
        - Align each read to his reference.
        - find possible translocations.
        - Filter noisy alignment.
        :param reads_df: all reads with a left and right matches for each site
        :param exp_type: ExpType
        :param allow_trans: flag to allow translocations check
        :return: A tuple of (dictionary with key=SITE_NAME and value=ReadsDf, translocations DF)
        """
        unmatched_df = reads_df.loc[(reads_df[L_SITE] != reads_df[R_SITE]) | (reads_df[L_REV] != reads_df[R_REV])].copy()
        reads_df.drop(index=unmatched_df.index, inplace=True)

        # Find translocations or match ambiguous reads
        trans_df, re_matched_df = self._find_translocations_or_new_matching(unmatched_df, allow_trans)

        # Add re matched reads to reads pool
        reads_df = pd.concat([reads_df, re_matched_df])
        reads_df[SITE_NAME] = reads_df[L_SITE]
        reads_df[REVERSED] = reads_df[L_REV]
        reads_df.drop(columns=[L_SITE, L_REV, R_SITE, R_REV, R_READ, L_READ], inplace=True)
        reads_df.reset_index(drop=True, inplace=True)

        # Detect bad amplicons (search for high frequency reads without a site match)
        self._detect_bad_amplicons(unmatched_df)

        # Reverse all the reversed reads
        rev_reads = reads_df.loc[reads_df[REVERSED], READ]
        reads_df.loc[rev_reads.index, READ] = rev_reads.apply(self.reverse_complement)

        # Align reads to their amplicon
        self._logger.debug("Alignment for {} - Start Needleman-Wunsch alignment for all reads.".format(exp_type.name()))

        self._align_reads_to_amplicon(reads_df)

        # Filter reads with low alignment score
        unaligned_df = self._filter_low_score_reads(unmatched_df, reads_df, exp_type)

        aligned_reads_num = reads_df[FREQ].sum()
        self._logger.debug("Alignment for {} - Needleman-Wunsch alignment done.".format(exp_type.name()))

        # Shift modification into cut-site
        self._logger.debug("Alignment for {} - Start shift modifications into cut-site.".format(exp_type.name()))
        self._shift_modifications_into_cut_site(reads_df)
        self._logger.debug("Alignment for {} - Shift modifications into cut-site done.".format(exp_type.name()))

        # Split read_df to all the different sites
        read_d: ReadsDict = dict()
        for site in self._ref_df.index:
            read_d[site] = reads_df.loc[reads_df[SITE_NAME] == site].sort_values(by=[FREQ],
                                                                                 ascending=False).reset_index(drop=True)
            # Add indels columns to reads df
            self._add_indels_columns(read_d[site], self._ref_df.loc[site, CUT_SITE])

        self._logger.info("Alignment for {} - Done.".format(exp_type.name(), self._ref_df.shape[0]))
        self._aligned_n[exp_type] = aligned_reads_num

        # TODO - delete!!!!!!!
        if debug:
            filtered_df = pd.concat([unaligned_df, unmatched_df])
            filtered_df.reset_index(inplace=True, drop=True)

            # prepare reference list and names
            self._logger.debug("Alignment for {} - debug start".format(exp_type.name()))

            references = list(self._ref_df[REFERENCE].values)
            names = list(self._ref_df[SITE_NAME].values)
            references += list(self._ref_df[REFERENCE].apply(self.reverse_complement))
            names += names
            ref_rev = self._ref_df.shape[0]*[False] + self._ref_df.shape[0]*[True]
            ref_scores = list(self._ref_df[MAX_SCORE].values) + list(self._ref_df[MAX_SCORE].values)

            # full demultiplexing
            filtered_df[['presumed_site','is_rev','score_normal']] = filtered_df.apply((lambda row: self.debug_match(row[READ],references, names, ref_rev, ref_scores)), axis=1, result_type='expand')

            # Reverse all the reversed reads
            rev_reads = filtered_df.loc[filtered_df['is_rev'], READ]
            filtered_df.loc[rev_reads.index, READ] = rev_reads.apply(self.reverse_complement)

            # full alignment
            self.debug_align(filtered_df)
            filtered_df.drop(inplace=True, columns=['left_primer_read','left_site_name','left_site_reversed',
                                                    'right_primer_read', 'right_site_reversed', 'right_site_name',
                                                    'is_rev', 'reversed_reads', 'site_name'])

            # draw alignment
            filtered_df.reset_index(inplace=True, drop=True)
            filtered_df['read_id'] = filtered_df.index
            filtered_df = pd.concat([filtered_df, filtered_df, filtered_df]).reset_index(drop=True)
            filtered_df.sort_values(by=['read_id'], ascending=False)

            filtered_df['left_alignment'] = None
            filtered_df['right_alignment'] = None
            filtered_df[CUT_SITE] = None
            filtered_df['indel_on_cutsite'] = False
            for read_id in filtered_df['read_id'].unique():
                filtered_read_df = filtered_df.loc[filtered_df['read_id'] == read_id]
                ref_cut_site = self._ref_df.at[filtered_read_df['presumed_site'].values[0], CUT_SITE]
                filtered_df.at[filtered_read_df.index, CUT_SITE] = ref_cut_site
                # Find alignment
                pos_idx = 0
                align_idx = 0
                aligned_cut_site = 0
                for length, indel_type in InputProcessing.parse_cigar(filtered_read_df[CIGAR].values[0]):
                    if indel_type != IndelType.INS:
                        if ref_cut_site in range(pos_idx, pos_idx + length + 1):
                            aligned_cut_site = align_idx + (ref_cut_site - pos_idx)
                        pos_idx += length
                    align_idx += length

                # Draw alignment
                for idx, (read_idx, row) in enumerate(filtered_read_df.iterrows()):
                    # draw ref
                    if idx == 0:
                        filtered_df.at[read_idx, 'left_alignment'] = row[ALIGNMENT_W_INS][0:aligned_cut_site]
                        filtered_df.at[read_idx, 'right_alignment'] = row[ALIGNMENT_W_INS][aligned_cut_site:]
                        if (row['alignment_human_readable'][aligned_cut_site-5:aligned_cut_site] != "|||||") or\
                           (row['alignment_human_readable'][aligned_cut_site:aligned_cut_site+5]!= "|||||"):
                            filtered_df.at[read_idx, 'indel_on_cutsite'] = True
                    # draw lines
                    elif idx == 1:
                        filtered_df.at[read_idx, 'left_alignment'] = row['alignment_human_readable'][0:aligned_cut_site]
                        filtered_df.at[read_idx, 'right_alignment'] = row['alignment_human_readable'][aligned_cut_site:]
                        filtered_df.at[read_idx, FREQ] = 0
                        if (row['alignment_human_readable'][aligned_cut_site-5:aligned_cut_site] != "|||||") or\
                           (row['alignment_human_readable'][aligned_cut_site:aligned_cut_site+5]!= "|||||"):
                            filtered_df.at[read_idx, 'indel_on_cutsite'] = True
                    # draw read
                    else:
                        filtered_df.at[read_idx, 'left_alignment'] = row[ALIGNMENT_W_DEL][0:aligned_cut_site]
                        filtered_df.at[read_idx, 'right_alignment'] = row[ALIGNMENT_W_DEL][aligned_cut_site:]
                        filtered_df.at[read_idx, FREQ] = 0
                        if (row['alignment_human_readable'][aligned_cut_site-5:aligned_cut_site] != "|||||") or\
                           (row['alignment_human_readable'][aligned_cut_site:aligned_cut_site+5]!= "|||||"):
                            filtered_df.at[read_idx, 'indel_on_cutsite'] = True
            filtered_df = filtered_df.reindex(columns=['read_id','frequency','presumed_site','cigar_length', 'cigar_path',
                                                       'score_normal', 'alignment_score', 'indel_on_cutsite', 'left_alignment',
                                                       'right_alignment','cut-site','read','alignment_human_readable',
                                                       'alignment_w_del','alignment_w_del'])
            filtered_df.to_csv(os.path.join(self._output, "{}_filtered.csv".format(exp_type.name())), index=False)
            self._logger.debug("Alignment for {} - debug done".format(exp_type.name()))

        return read_d, trans_df

    # TODO - delete
    def debug_match(self, read: DNASeq, references: List[DNASeq], names: List[str],
                    reverse_l: List[bool], ref_score_l):
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

    # TODO - delete
    def debug_align(self, reads: ReadsDf):
        """
        - Align all reads to their amplicon.
        - Compute cigar path.
        - Compute score
        :param reads: all reads
        :return:
        """

        new_cols_d = defaultdict(list)
        for row_idx, row in reads.iterrows():
            alignments = self._aligner.align(self._ref_df[REFERENCE][row['presumed_site']], row[READ])
            [ref_with_ins, align_human, read_with_del, _] = format(alignments[0]).split("\n")
            cigar, cigar_len = self._compute_cigar_path_from_alignment(reference=ref_with_ins, read=read_with_del)
            align_score = alignments[0].score

            new_cols_d[ALIGNMENT_W_INS].append(ref_with_ins)
            new_cols_d[ALIGNMENT_W_DEL].append(read_with_del)
            new_cols_d[CIGAR_LEN].append(cigar_len)
            new_cols_d[CIGAR].append(cigar)

            new_cols_d[ALIGN_SCORE].append(align_score)
            new_cols_d['alignment_human_readable'].append(align_human)

        for col_name, col in new_cols_d.items():
            reads[col_name] = col

    #-------------------------------#
    ######### Private methods #######
    #-------------------------------#

    ######### Demultiplex ###########

    def _compute_read_primer_matching(self, reads: List[DNASeq], primers: List[DNASeq], primers_revered: List[bool],
                                      primers_names: List[str], max_edit_distance: int) -> Dict[DNASeq, Tuple[str, bool]]:
        """
        Create a dictionary with a match between each possible read (key) and a primer (val).
        :param reads:  partial (left or right) reads
        :param primers: List of primers sequences
        :param primers_revered: List that indicates if primers reversed or not
        :param primers_names: primers names
        :param max_edit_distance: Max edit distance to account as a match
        :return: Dict[read, Tuple[site_name, reversed_flag]]
        """
        match = dict()
        for read in reads:
            match[read] = self._match_by_edit_distance(read, primers, primers_revered, primers_names, max_edit_distance)

        return match

    @staticmethod
    def _match_by_edit_distance(read: DNASeq, primers: List[DNASeq], primers_revered: List[bool],
                                primers_names: List[str], max_edit_distance: int) -> Tuple[str, bool]:
        """
        Find for every read the most similar sequence (primer).
        If the minimum edit distance is above DNASeq, no primer is matched.
        :param read:  reads
        :param primers: List of primers sequences
        :param primers_revered: List that indicates if primers reversed or not
        :param primers_names: primers names
        :param max_edit_distance: Max edit distance to account as a match
        :return: site_name, reversed flag
        """
        min_name = None
        min_dist = max_edit_distance + 1
        min_reversed = False

        for primer, reverse, name  in zip(primers, primers_revered, primers_names):
            d = edlib.align(read, primer, k=min_dist)['editDistance']
            if d < min_dist and d != -1:
                min_dist = d
                min_name = name
                min_reversed = reverse

        return min_name, min_reversed

    def _match_by_full_alignment(self, read: DNASeq, references: List[DNASeq], names: List[str]) -> str:
        """
        Find read with highest alignment score.
        Time consuming function.
        :param read: DNASeq
        :return: site_name
        """
        alignments = self._aligner.align(references[0], read)
        max_name = names[0]
        max_score = alignments[0].score

        for reference, name  in zip(references[1:], names[1:]):
            alignments = self._aligner.align(reference, read)
            score = alignments[0].score
            if score > max_score:
                max_score = score
                max_name = name

        return max_name

    ######### Translocations ########

    def _get_translocation_reference(self, l_name: str, l_rev: bool, r_name: str, r_rev: bool) -> Tuple[DNASeq, int, float]:
        """
        Compute possible translocation for between left site and right site
        :param l_name: left site match
        :param l_rev: left site reversed flag
        :param r_name: right site match
        :param r_rev: right site reversed flag
        :return: reference sequence, expected_cut_site and max alignment score
        """
        # Left part of the reference
        l_ref = self._ref_df[REFERENCE][l_name]
        l_cut_site = self._ref_df[CUT_SITE][l_name]
        l_ref = self.reverse_complement(l_ref)[l_cut_site:] if l_rev else l_ref[:l_cut_site]

        # Right part of the reference
        r_ref = self._ref_df[REFERENCE][r_name]
        r_cut_site = self._ref_df[CUT_SITE][r_name]
        r_ref = self.reverse_complement(r_ref)[:r_cut_site] if r_rev else r_ref[r_cut_site:]

        reference = l_ref + r_ref
        _, _, _, _, max_score = self._compute_needle_wunsch_alignment(reference=reference, read=reference)

        return reference, l_cut_site, max_score

    def _find_translocations_or_new_matching(self, unmatched_df: ReadsDf, allow_trans: bool) -> Tuple[TransDf, ReadsDf]:
        """
        For every non None match, find the best alignment betwwen: left site, right site and translocation between
        left & right.
        :param unmatched_df: unmatched reads
        :param allow_trans: If False, consider translocation as unmatched read.
        :return:
        """
        trans_d = defaultdict(list)
        trans_idx_list = [] # list of matched translocation
        re_matched_idx_list = [] # list of re matched indexes (initial bad demultiplexing)
        # trans_ref dict - Key is translocation ID and value is reference, cut_site & max alignment score
        trans_ref_d: Dict[Tuple[str, bool, str, bool], Tuple[DNASeq, int, float]] = dict()

        possible_trans_df = unmatched_df.loc[unmatched_df[L_SITE].notnull() & unmatched_df[R_SITE].notnull()]

        def get_trans_name(name, rev):
            direction = '5\'_to_3\'' if not rev else '3\'_to_5\''
            return name + "_" + direction

        if allow_trans:
            self._logger.info("Search possible Translocations for {} reads. May take a few minutes".format(
                unmatched_df.shape[0]))

        for idx, row in possible_trans_df.iterrows():
            # Compute translocation reference, cut-site & max alignment score
            l_name, l_rev, r_name, r_rev = row[[L_SITE, L_REV, R_SITE, R_REV]]
            if (l_name, l_rev, r_name, r_rev) in trans_ref_d:
                reference, cut_site, max_score = trans_ref_d[(l_name, l_rev, r_name, r_rev)]
            else:
                reference, cut_site, max_score = self._get_translocation_reference(l_name, l_rev, r_name, r_rev)
                trans_ref_d[(l_name, l_rev, r_name, r_rev)] = reference, cut_site, max_score

            # Compute alignment
            ref_w_ins, read_w_del, cigar, cigar_len, align_score = \
                self._compute_needle_wunsch_alignment(reference=reference, read=row[READ])
            normalized_score = align_score / max_score

            # Compute left site and right site alignment score - to verify that the highest is the translocation
            l_ref = self.reverse_complement(self._ref_df[REFERENCE][l_name]) if l_rev else self._ref_df[REFERENCE][l_name]
            _, _, _, _, l_align_score = self._compute_needle_wunsch_alignment(reference=l_ref, read=row[READ])
            l_normalized_score = l_align_score / self._ref_df[MAX_SCORE][l_name]

            r_ref = self.reverse_complement(self._ref_df[REFERENCE][r_name]) if r_rev else self._ref_df[REFERENCE][r_name]
            _, _, _, _, r_align_score = self._compute_needle_wunsch_alignment(reference=r_ref, read=row[READ])
            r_normalized_score = r_align_score / self._ref_df[MAX_SCORE][r_name]

            # Check if left or right sites have better alignment than the translocation
            if (r_normalized_score > normalized_score) or (l_normalized_score > normalized_score):
                if r_normalized_score > l_normalized_score:
                    unmatched_df.at[idx, L_SITE] = row[R_SITE]
                    unmatched_df.at[idx, L_REV] = row[R_REV]
                else:
                    unmatched_df.at[idx, R_SITE] = row[L_SITE]
                    unmatched_df.at[idx, R_REV] = row[L_REV]
                re_matched_idx_list.append(idx)
            # If read has high quality alignment, consider it as translocation
            # TODO - change min score?
            elif (align_score > (self._min_score / 100) * max_score) or (cigar_len < CIGAR_LEN_THRESHOLD):
                trans_idx_list.append(idx)
                trans_d[READ].append(row[READ])
                trans_d[FREQ].append(row[FREQ])
                trans_d[TRANS_NAME].append(get_trans_name(l_name, l_rev) + "_" + get_trans_name(r_name, r_rev))
                trans_d[REFERENCE].append(reference)
                trans_d[CUT_SITE].append(cut_site)
                trans_d[ALIGNMENT_W_INS].append(ref_w_ins)
                trans_d[ALIGNMENT_W_DEL].append(read_w_del)
                trans_d[CIGAR].append(cigar)
                trans_d[ALIGN_SCORE].append(align_score)
                trans_d[R_SITE].append(l_name)
                trans_d[L_SITE].append(r_name)
                trans_d[CIGAR_LEN].append(cigar_len)
                # TODO - for debug only
                trans_d['max_score'].append(max_score)
                trans_d['normalized_score'].append(normalized_score)


        # Convert translocation dictionary to translocation pandas
        if allow_trans:
            self._logger.info("Search possible Translocations - Done.")
            trans_df = pd.DataFrame.from_dict(trans_d, orient='columns')
            # TODO - for debug only, change to by FREQ
            if trans_df.shape[0] > 0:
                trans_df.sort_values(by=['normalized_score'], ascending=False, inplace=True)

            # Remove translocation from unmatched reads
            unmatched_df.drop(index=trans_idx_list, inplace=True)
        else:
            trans_df = pd.DataFrame()

        # Create new df for re matched reads
        re_matched_df = unmatched_df.copy()
        unmatched_df.drop(index=re_matched_idx_list, inplace=True)
        re_matched_df.drop(index=unmatched_df.index, inplace=True)

        return trans_df, re_matched_df

    ######### Alignment #############

    def _detect_bad_amplicons(self, unmatched_df: ReadsDf):
        """
        - detect unmatched reads with large frequency and warn user from bad amplicon.
        :param unmatched_df:
        :return:
        """
        # prepare reference list and names
        references = list(self._ref_df[REFERENCE].values)
        names = list(self._ref_df[SITE_NAME].values)
        references += list(self._ref_df[REFERENCE].apply(self.reverse_complement))
        names += ["{}_reversed".format(name) for name in names]

        high_freq_df = unmatched_df.loc[unmatched_df[FREQ] > BAD_AMPLICON_THRESHOLD].copy()
        high_freq_df.sort_values(by=[FREQ], ascending=False, inplace=True)

        for _, row in high_freq_df.iterrows():
            site_name = self._match_by_full_alignment(row[READ], references, names)
            self._logger.warning("The following read appears {}, but it doesn't match any site!"
                                 "This read is most similar to {}. Please check amplicons correctness"
                                 "read={}".format(row[FREQ], site_name, row[READ]))

    def _filter_low_score_reads(self, unmatched_df: ReadsDf, reads_df: ReadsDf, exp_type: ExpType):
        """
        - Filter low alignment score reads
        - Store all unaligned reads in a fasta format.
        :param reads_df: aligned reads df.
        :param exp_type:
        :return:
        """
        # Find all indexes with lower score than the threshold
        unaligned_indexes = []
        for site in self._ref_df[SITE_NAME].keys():
            site_df = reads_df.loc[reads_df[SITE_NAME] == site]
            max_score = site_df[ALIGN_SCORE].max()
            score_threshold = (self._min_score/100)*max_score
            # Filter low alignment score
            low_score_df = site_df.loc[(site_df[ALIGN_SCORE] < score_threshold) & (site_df[CIGAR_LEN] > CIGAR_LEN_THRESHOLD)]
            # Filter PRIMER-DIMER affect
            min_len = len(self._ref_df[F_PRIMER][site]) + len(self._ref_df[R_PRIMER][site]) + MIN_PRIMER_DIMER_THRESH
            short_reads_df = site_df.loc[site_df[READ].str.len() < min_len]
            unaligned_indexes += list(low_score_df.index) + list(short_reads_df.index)

        total_reads_num = reads_df[FREQ].sum() + unmatched_df[FREQ].sum()
        unaligned_df = reads_df.loc[unaligned_indexes]
        unaligned_reads_num = unaligned_df[FREQ].sum() + unmatched_df[FREQ].sum()

        self._logger.info("Alignment for {} - {:,} reads weren't aligned ({:.2f}% of all reads)".format(exp_type.name(),
                          unaligned_reads_num, 100*unaligned_reads_num/total_reads_num))

        with open(os.path.join(self._output, "{}_unaligned_reads.fasta".format(exp_type.name())), 'w') as file:
            for _, row in unmatched_df.iterrows():
                file.write("> unaligned read with {} copies in the original fastq file.\n".format(row[FREQ]))
                file.write("{}\n".format(row[READ]))

            for _, row in unaligned_df.iterrows():
                file.write("> unaligned read with {} copies in the original fastq file. Read most "
                           "similar to site {}\n".format(row[FREQ], row[SITE_NAME]))
                file.write("{}\n".format(row[READ]))

        # remove unaligned reads from reads_df
        reads_df.drop(unaligned_df.index, inplace=True)
        return unaligned_df # TODO - remove

    def _align_reads_to_amplicon(self, reads: ReadsDf):
        """
        - Align all reads to their amplicon.
        - Compute cigar path.
        - Compute score
        :param reads: all reads
        :return:
        """

        new_cols_d = defaultdict(list)
        for row_idx, row in reads.iterrows():
            ref_w_ins, read_w_del, cigar, cigar_len, align_score = self._compute_needle_wunsch_alignment(
                reference=self._ref_df[REFERENCE][row[SITE_NAME]], read=row[READ])

            new_cols_d[ALIGNMENT_W_INS].append(ref_w_ins)
            new_cols_d[ALIGNMENT_W_DEL].append(read_w_del)
            new_cols_d[CIGAR_LEN].append(cigar_len)
            new_cols_d[CIGAR].append(cigar)

            new_cols_d[ALIGN_SCORE].append(align_score)

        for col_name, col in new_cols_d.items():
            reads[col_name] = col

    def _compute_needle_wunsch_alignment(self, reference: DNASeq, read: DNASeq) -> Tuple[DNASeq, DNASeq, CigarPath, int, float]:
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
        window_size = self._cfg["window_size"]
        start_idx = cut_site - window_size  # Start index to include indel
        end_idx = cut_site + window_size  # end index to include indel
        new_col_d = dict()

        for col in INDEL_COLS:
            new_col_d[col] = reads.shape[0] * [None]

        for row_idx, row in reads.iterrows():
            pos_idx = 0  # position index for the original reference (with no indels)
            align_idx = 0  # position index for the alignment (with indels)
            reference = row[ALIGNMENT_W_INS]
            read = row[ALIGNMENT_W_DEL]
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
                            new_col_d[DEL_BASE][row_idx] = reference[align_idx:align_idx+length]
                        else:
                            new_col_d[DEL_LEN][row_idx] += ", {}".format(length)
                            new_col_d[DEL_START][row_idx] += ", {}".format(pos_idx)
                            new_col_d[DEL_END][row_idx] += ", {}".format(pos_idx + length - 1)
                            new_col_d[DEL_BASE][row_idx] += ", {}".format(reference[align_idx:align_idx+length])

                    pos_idx += length

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

                    pos_idx += length

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

                align_idx += length

        # Add the cut-site
        reads[CUT_SITE] = reads.shape[0] * [cut_site]

        # Add the new
        for col in INDEL_COLS:
            reads[col] = new_col_d[col]

    # TODO - keep only for debug?
    def compute_alignment_score_from_cigar(self, cigar):
        # Get alignment config
        cfg = Configurator.get_cfg()
        align_cfg = cfg["alignment"]
        match_score = align_cfg["match_score"]
        mismatch_score = align_cfg["mismatch_score"]
        open_gap_score = align_cfg["open_gap_score"]
        extend_gap_score = align_cfg["extend_gap_score"]
        score = 0

        # Find closest indels from left and from right
        for length, indel_type in self.parse_cigar(cigar):
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

    def _find_closest_indels_to_cut_site(self, cigar: CigarPath, cut_site: int) -> Tuple[AlignedIndel, AlignedIndel, int]:
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
        for length, indel_type in self.parse_cigar(cigar):

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

    def _shift_modifications_into_cut_site(self, reads: ReadsDf):
        """
        Shift deletions and insertions with a region of 2*window_size into the cut-site
        :param reads: ReadsDf - All reads, already aligned with a cigar path.
        :return: no return value. Change reads inplace.
        """
        window_size = self._cfg["window_size"]
        update_idx = []  # Changed read indexes list
        reference_l = []  # Changed references (read with insertions) list
        read_l = []  # Changed reads (read with deletions) list
        cigar_l = []  # Changed cigar path list

        for row_idx, row in reads.iterrows():
            cut_site = self._ref_df.loc[row[SITE_NAME], CUT_SITE]  # cut-site in pos_idx coordinates
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
                if align_idx + length + (2*window_size) < aligned_cut_site:
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
                if align_idx > aligned_cut_site + (2*window_size):
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
                # TODO - delete this assertion from final version
                original_score = self.compute_alignment_score_from_cigar(row[CIGAR])
                new_score = self.compute_alignment_score_from_cigar(cigar)
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

    ######### General ###############

    @staticmethod
    def reverse_complement(seq: DNASeq) -> DNASeq:
        return "".join(COMPLEMENT.get(base, base) for base in reversed(seq))

    def _get_expected_cut_site(self, reference: DNASeq, sgRNA: DNASeq, cut_site_position: int, site_name: str = '') \
        -> Tuple[int, bool]:
        """
        Find sgRNA (or the reverse complement) inside the reference and return the expected cut-site.
        The cut-site is LEFT to returned index
        :param reference: reference sequence
        :param sgRNA: sgRNA sequence
        :param cut_site_position: position relative to the PAM
        :param site_name: site name
        :return: expected cut-site, reversed sgRNA or not
        """
        reverse = False
        sgRNA_start_idx = reference.find(sgRNA)
        if sgRNA_start_idx == -1:
            sgRNA_start_idx = reference.find(self.reverse_complement(sgRNA))
            if sgRNA_start_idx == -1:
                raise SgRNANotInReferenceSequence(site_name)
            else:
                cut_site = sgRNA_start_idx - cut_site_position
                reverse = True
        else:
            cut_site = sgRNA_start_idx + len(sgRNA) + cut_site_position

        return cut_site, reverse

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
    def parse_cigar_with_adjacent_indels(cigar: str) -> List[Tuple[int, int, IndelType]]:
        """
        Function returns a list of Tuple[indel length, IndelType] where adjacent indels are aggregated
        into indel type = Indel.
        :param cigar:
        :return: list of Tuple[indel length, indel length without insertion, IndelType]
        """
        indel_list = []
        prev_indel = IndelType.MATCH
        prev_length = 0
        prev_length_wo_ins = 0  # used to understand with positions are relevant for this modification
        for length, indel in re.findall(r'(\d+)([{}{}{}{}])'.format(CIGAR_D, CIGAR_I, CIGAR_S, CIGAR_M), cigar):
            indel = IndelType.from_cigar(indel)
            length = int(length)
            length_wo_ins = int(length) if indel != IndelType.INS else 0
            if (indel != IndelType.MATCH) and (prev_indel != IndelType.MATCH):
                indel = IndelType.MIXED
                length += prev_length
                length_wo_ins += prev_length_wo_ins
                indel_list[-1] = (length, length_wo_ins, indel)
            else:
                indel_list.append((length, length_wo_ins, indel))

            # update prev
            prev_indel = indel
            prev_length = length
            prev_length_wo_ins = length_wo_ins

        return indel_list

    ######### Getters #############
    def read_numbers(self, exp_type: ExpType) -> Tuple[int, int, int]:
        """
        return a tuple of number of input, merged & aligned
        :param exp_type:
        :return:
        """
        return self._input_n[exp_type], self._merged_n[exp_type], self._aligned_n[exp_type]
