from utils.constants_and_types import COMPLEMENT, IndelType, Path, DNASeq, CIGAR_D, CIGAR_I, CIGAR_S, CIGAR_M
from typing import List, Tuple
import re
from utils.exceptions import CantOpenMergedFastqFile


def reverse_complement(seq: DNASeq) -> DNASeq:
    return "".join(COMPLEMENT.get(base, base) for base in reversed(seq))


def parse_fastq_file(file_name: Path) -> List[DNASeq]:
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
