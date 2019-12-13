#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `crispector` package."""

import pytest
from Bio import Align

from click.testing import CliRunner

from constants_and_types import CIGAR, SITE_NAME, ALIGNMENT_W_INS, ALIGNMENT_W_DEL
from crispector import crispector
from crispector import cli
from input_processing import InputProcessing
from utils import Configurator
from fixtures import rag1_amplicons

@pytest.fixture
def response():
    """Sample pytest fixture.

    See more at: http://doc.pytest.org/en/latest/fixture.html
    """
    # import requests
    # return requests.get('https://github.com/audreyr/cookiecutter-pypackage')


def test_content(response):
    """Sample pytest test function with the pytest fixture as an argument."""
    # from bs4 import BeautifulSoup
    # assert 'GitHub' in BeautifulSoup(response.content).title.string


def test_command_line_interface():
    """Test the CLI."""
    runner = CliRunner()
    result = runner.invoke(cli.main)
    assert result.exit_code == 0
    assert 'crispector.cli.main' in result.output
    help_result = runner.invoke(cli.main, ['--help'])
    assert help_result.exit_code == 0
    assert '--help  Show this message and exit.' in help_result.output


def test_read_aligment(rag1_amplicons):
    import time
    amplicons, amplicons_names, max_amplicon_len = rag1_amplicons
    aligner = Align.PairwiseAligner()
    aligner.match_score = 5
    aligner.mismatch_score = -4
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    read = "CCCAGCAAACAATTGAGTGACCAAATAAATGAATGACTCTGTGGCTAGCACAGGACGCGCCAGGCAGCACATGGAGGGTCTGTAGTAA"
    start_time = time.time()
    InputProcessing._match_by_edit_distance(read, amplicons, amplicons_names, max_amplicon_len)
    print(time.time()-start_time)

    start_time = time.time()
    InputProcessing._match_by_alignment_distance(read, amplicons, amplicons_names, aligner)
    print(time.time() - start_time)


def test_modification_change():
    import pandas as pd
    import numpy as np
    Configurator.set_cfg_path("/home/toky/Ido/crispector/crispector/config/default_config.yml")
    # ref = "TGAGAGACCAGTGCTGAGGGAGTGGGAGAGGGCGCTGTGCCAGGCAAGCAAGAGGTGGCTGTGGAGAGGGCAGTGCGGGGATGGAGCTGACCACGTTGGGGAGTCCCAGGGTGGGTGGGAAGAATGTAAGTTTGGAAATGATAGTATCCAGGTCAGTTGTCAGGGAG"
    # read = "TGAGAGACCAGTGCTGAGGGAGTGGGAGAGGGCGCTGTGCCAGGCAAGCAAGAGGTGGCTGTGGAGAGGGCAGTGCGGGGATGGAGCTGACCACGTTGGGGAGTCCCA----GGGTGGGAAGAATGTAAGTTTGGAAATGATAGTATCCAGGTCAGTTGTCAGGGAT"
    ref = "AAGCAGCCCCGAGATAATCCAAGGATCTCACCTCTTTAGATTCTGAGCGGAGTCCCGAGGTGTGTGGGCTGGCTTGCCCCGGAGGGAGGCAGGTCGCGAAT"
    # ref = "CCCAGCAAACAATTGAGTGACCAAATAAATGAATGACTCTGTGGCTAGCACAGGACGCGCCAGGCAGCACCCCTCCCACCCTTGGGGCTCAGGGGTAGCCCGTGTGAACCCTGACTCCTACCCACCCACCGCCCCCCTCACTTGCTGGGACCCAGCCATGGGGTTGAACAGGGAGCGGCCAGAGAAGGAGATCACAGGGATT"
    print(len(ref))
    read = "CCCAGCAAACAATTGAGTGACCAAATAAATGAATGACTCTGTGGCTAGCACAGGACGCGCCAGGCAGCACATGGAGGGTCTGTAGTAA"
    print(len(read))

    cut_site = 113
    df = pd.DataFrame()
    cigar = InputProcessing._compute_cigar_path_from_alignment(ref, read)
    df[ALIGNMENT_W_INS] = [np.nan, ref]
    df[ALIGNMENT_W_DEL] = [read, read]
    df[CIGAR] = [cigar, cigar]
    df[SITE_NAME] = ["site_1", "site_1"]
    df = df.dropna()
    cut_site_d = {"site_1": cut_site}
    print("before")
    print("reference\t= {}".format(ref))
    print("read\t\t= {}".format(read))
    print("cigar\t\t= {}".format(cigar))
    InputProcessing._shift_modifications_into_cut_site(df, cut_site_d, 20)
    print("after")
    print("reference\t= {}".format(df[ALIGNMENT_W_INS]))
    print("read\t\t= {}".format(df[ALIGNMENT_W_DEL]))
    print("cigar\t\t= {}".format(df[CIGAR]))
    check = 0
