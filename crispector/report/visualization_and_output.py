from crispector.utils.constants_and_types import SITE_NAME, IndelType, Path, FREQ, IS_EDIT, TX_READ_NUM, \
    MOCK_READ_NUM, TX_EDIT, EDIT_PERCENT, CI_LOW, CI_HIGH, READ_LEN_SIDE, ALIGN_CUT_SITE, ALIGNMENT_W_INS, \
    ALIGNMENT_W_DEL, POS_IDX_E, POS_IDX_S, INDEL_TYPE, ExpType, ON_TARGET, IsEdit, C_TX, C_MOCK, \
    SUMMARY_RESULTS_TITLES, OFF_TARGET_COLOR, ON_TARGET_COLOR, OUTPUT_DIR, DISCARDED_SITES, \
    EDITING_ACTIVITY, PLOT_PATH, TITLE, W, H, PDF_PATH, PAGE_TITLE, READING_STATS, MAPPING_STATS, MAPPING_PER_SITE, \
    FASTP_DIR, FASTP_TX_PATH, FASTP_MOCK_PATH, RESULT_TABLE, TAB_DATA, HTML_SITES, LOG_PATH, TransDf, AlgResultDf, \
    TransResultDf, TRANS_FDR, SITE_A, SITE_B, TX_TRANS_READ, TRANSLOCATIONS, TX_TRANS_PATH, MOCK_TRANS_PATH, \
    TRANS_RES_TAB, TRANS_HEATMAP_TAB, TRANS_RESULTS_TITLES, EDIT_SECTION, MOD_SECTION, CLS_RES_SECTION, CLS_RES_INS, \
    CLS_RES_DEL, CLS_RES_MIX, MOD_DIST, EDIT_DIST, EDIT_SIZE_DIST, READ_SECTION, READ_EDIT, READ_MOCK_ALL, READ_TX_ALL, \
    FILTERED_PATH, READ_TX_FILTER, READ_MOCK_FILTER, HTML_SITES, HTML_SITES_NAME_LIST, REPORT_PATH, LOGO_PATH, \
    EDIT_TEXT, UNBALANCED_READ_WARNING, UNMATCHED_PATH, UNMATCHED_TX_PATH, UNMATCHED_MOCK_PATH
import math
import os
import warnings
from typing import List, Tuple, Dict
from crispector.input_processing.input_processing import InputProcessing
from crispector.modifications.modification_types import ModificationTypes
from crispector.algorithm.core_algorithm import CoreAlgorithm
from crispector.modifications.modification_tables import ModificationTables
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import seaborn as sns #
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
from crispector.utils.logger import LoggerWrapper
from copy import deepcopy
from matplotlib.collections import QuadMesh
import shutil


def create_site_output(algorithm: CoreAlgorithm, modifications: ModificationTypes, mod_table: ModificationTables,
                       html_param_d: Dict, site_result: Dict, site_name: str, output: Path):

    html_param_d[HTML_SITES][HTML_SITES_NAME_LIST] += [site_name]
    html_param_d[HTML_SITES][site_name] = dict() # HTML parameters dict
    html_d = html_param_d[HTML_SITES][site_name]
    html_d[TITLE] = "{}".format(site_name)
    html_d[REPORT_PATH] = os.path.join(OUTPUT_DIR, site_name, "report.html")
    base_path = ""

    cut_site = algorithm.cut_site
    win_size = algorithm.win_size

    # Plot editing activity
    plot_site_editing_activity(algorithm, site_result, site_name, output, html_d, base_path)

    # Plot mutation distribution
    html_d[MOD_SECTION] = dict()
    html_d[MOD_SECTION][TITLE] = "Indels Level Information"
    plot_distribution_of_edit_events(mod_table, cut_site, win_size, output, html_d, base_path)
    plot_distribution_of_all_modifications(mod_table, cut_site, win_size, output, html_d, base_path)
    plot_distribution_of_edit_event_sizes(mod_table, output, html_d, base_path)

    # plot modification table
    html_d[CLS_RES_SECTION] = dict()
    html_d[CLS_RES_SECTION][TITLE] = "Classifier Results for Every Indel Type and Reference Position"
    del_table_idx = [i for i, x in enumerate(modifications.types) if x == IndelType.DEL]
    ins_table_idx = [i for i, x in enumerate(modifications.types) if x == IndelType.INS]
    mix_table_idx = [i for i, x in enumerate(modifications.types) if x in [IndelType.MIXED, IndelType.SUB]]
    plot_modification_tables(mod_table, modifications, algorithm.edited, algorithm.tables_offset,
                             del_table_idx, output, html_d, base_path, "del")
    plot_modification_tables(mod_table, modifications, algorithm.edited, algorithm.tables_offset,
                             ins_table_idx, output, html_d, base_path, "ins")
    plot_modification_tables(mod_table, modifications, algorithm.edited, algorithm.tables_offset,
                             mix_table_idx, output, html_d, base_path, "mix_and_sub")

    # Create edited read table
    html_d[READ_SECTION] = dict()
    html_d[READ_SECTION][TITLE] = "Reads Level Information"
    plot_edited_reads_to_table(mod_table, cut_site, output, html_d, base_path)

    # Dump .csv file with all reads
    mod_table.tx_reads.to_csv(os.path.join(output, "treatment_aligned_reads.csv.gz"), index=False, compression='gzip')
    mod_table.mock_reads.to_csv(os.path.join(output, "mock_aligned_reads.csv.gz"), index=False, compression='gzip')
    html_d[READ_SECTION][READ_TX_ALL] = os.path.join(base_path, "treatment_aligned_reads.csv.gz")
    html_d[READ_SECTION][READ_MOCK_ALL] = os.path.join(base_path, "mock_aligned_reads.csv.gz")

    # Add fastp links
    if os.path.exists(os.path.join(output, FASTP_DIR[ExpType.TX])):
        html_d[READ_SECTION][FASTP_TX_PATH] = os.path.join(base_path, "{}/fastp.html".format(FASTP_DIR[ExpType.TX]))
    else:
        html_d[READ_SECTION][FASTP_TX_PATH] = ""

    if os.path.exists(os.path.join(output, FASTP_DIR[ExpType.MOCK])):
        html_d[READ_SECTION][FASTP_MOCK_PATH] = os.path.join(base_path, "{}/fastp.html".format(FASTP_DIR[ExpType.MOCK]))
    else:
        html_d[READ_SECTION][FASTP_MOCK_PATH] = ""

    # Add filtered reads link
    if os.path.exists(os.path.join(output, FILTERED_PATH[ExpType.TX])):
        html_d[READ_SECTION][READ_TX_FILTER] = os.path.join(base_path, FILTERED_PATH[ExpType.TX])
    else:
        html_d[READ_SECTION][READ_TX_FILTER] = ""

    if os.path.exists(os.path.join(output, FILTERED_PATH[ExpType.MOCK])):
        html_d[READ_SECTION][READ_MOCK_FILTER] = os.path.join(base_path, FILTERED_PATH[ExpType.MOCK])
    else:
        html_d[READ_SECTION][READ_MOCK_FILTER] = ""


def create_experiment_output(result_df: AlgResultDf, tx_trans_df: TransDf, mock_trans_df: TransDf,
                             trans_result_df: TransResultDf, input_processing: InputProcessing, min_num_of_reads: int,
                             confidence_interval: float, editing_threshold: float, translocation_p_value: float,
                             output: Path, donor: bool):

    html_d = dict()  # html parameters dict

    # Dump summary results
    summary_result_to_excel(result_df, confidence_interval, output)

    # Create bar plot for editing activity
    html_d[EDIT_SECTION] = dict()
    html_d[EDIT_SECTION][TITLE] = "Editing Activity"
    plot_editing_activity(result_df, confidence_interval, editing_threshold, html_d, output)

    # Dump reads statistics
    tx_input_n, tx_merged_n, tx_aligned_n = input_processing.read_numbers(ExpType.TX)
    mock_input_n, mock_merged_n, mock_aligned_n = input_processing.read_numbers(ExpType.MOCK)
    html_d[READING_STATS] = dict()
    create_reads_statistics_report(result_df, tx_input_n, tx_merged_n, tx_aligned_n,
                                   mock_input_n, mock_merged_n, mock_aligned_n, html_d, output)

    # Create a text file with all discarded sites
    warnings_and_discarded_sites_text(result_df, min_num_of_reads, donor, html_d)

    html_d[TRANSLOCATIONS] = dict()
    html_d[TRANSLOCATIONS][TITLE] = "Translocations"
    # Dump all translocations reads
    if tx_trans_df.shape[0] > 0:
        tx_trans_df.to_csv(os.path.join(output, "tx_reads_with_primer_inconsistency.csv"), index=False)
        html_d[TRANSLOCATIONS][TX_TRANS_PATH] = os.path.join(OUTPUT_DIR, "tx_reads_with_primer_inconsistency.csv")
    else:
        html_d[TRANSLOCATIONS][TX_TRANS_PATH] = ""
    if mock_trans_df.shape[0] > 0:
        mock_trans_df.to_csv(os.path.join(output, "mock_reads_with_primer_inconsistency.csv"), index=False)
        html_d[TRANSLOCATIONS][MOCK_TRANS_PATH] = os.path.join(OUTPUT_DIR, "mock_reads_with_primer_inconsistency.csv")
    else:
        html_d[TRANSLOCATIONS][MOCK_TRANS_PATH] = ""

    # Save translocations results
    trans_result_df.to_csv(os.path.join(output, "translocations_results.csv"), index=False)
    if trans_result_df.shape[0] > 0:
        html_d[TRANSLOCATIONS][TRANS_RES_TAB] = dict()
        html_d[TRANSLOCATIONS][TRANS_RES_TAB][TITLE] = "Translocations Results Sorted by FDR Value"
        html_d[TRANSLOCATIONS][TRANS_RES_TAB][TAB_DATA] = dict()
        html_d[TRANSLOCATIONS][TRANS_RES_TAB][TAB_DATA][0] = list(trans_result_df.columns.values)
        for row_idx, row in trans_result_df.iterrows():
            html_d[TRANSLOCATIONS][TRANS_RES_TAB][TAB_DATA][row_idx+1] = list(row.values)
    else:
        html_d[TRANSLOCATIONS][TRANS_RES_TAB] = ""

    # Translocations Heatmap
    plot_translocations_heatmap(result_df, trans_result_df, translocation_p_value, html_d, output)

    # Add summary results to page
    html_d[EDIT_SECTION][RESULT_TABLE] = dict()
    html_d[EDIT_SECTION][RESULT_TABLE][TITLE] = "Results Table"
    html_d[EDIT_SECTION][RESULT_TABLE][TAB_DATA] = dict()
    html_d[EDIT_SECTION][RESULT_TABLE][TAB_DATA][0] = list(result_df.columns.values)
    for row_idx, row in result_df.iterrows():
        html_d[EDIT_SECTION][RESULT_TABLE][TAB_DATA][row_idx+1] = list(row.values)

    # Will be filled for each site
    html_d[HTML_SITES] = dict()
    html_d[HTML_SITES][HTML_SITES_NAME_LIST] = []

    # Add fastp links
    if os.path.exists(os.path.join(output, FASTP_DIR[ExpType.TX])):
        html_d[READING_STATS][FASTP_TX_PATH] = os.path.join(OUTPUT_DIR, "{}/fastp.html".format(FASTP_DIR[ExpType.TX]))
    else:
        html_d[READING_STATS][FASTP_TX_PATH] = ""

    if os.path.exists(os.path.join(output, FASTP_DIR[ExpType.MOCK])):
        html_d[READING_STATS][FASTP_MOCK_PATH] = os.path.join(OUTPUT_DIR, "{}/fastp.html".format(FASTP_DIR[ExpType.MOCK]))
    else:
        html_d[READING_STATS][FASTP_MOCK_PATH] = ""

    # Add unmacthed reads links
    if os.path.exists(os.path.join(output, UNMATCHED_PATH[ExpType.TX])):
        html_d[READING_STATS][UNMATCHED_TX_PATH] = os.path.join(OUTPUT_DIR, UNMATCHED_PATH[ExpType.TX])
    else:
        html_d[READING_STATS][UNMATCHED_TX_PATH] = ""

    if os.path.exists(os.path.join(output, UNMATCHED_PATH[ExpType.MOCK])):
        html_d[READING_STATS][UNMATCHED_MOCK_PATH] = os.path.join(OUTPUT_DIR, UNMATCHED_PATH[ExpType.MOCK])
    else:
        html_d[READING_STATS][UNMATCHED_MOCK_PATH] = ""

    html_d[LOG_PATH] = os.path.join(OUTPUT_DIR, LoggerWrapper.logger_name)

    # copy logo to user directory
    user_path = os.path.join(output, "crispector_logo.jpg")
    package_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'html_templates/crispector_logo.jpg')
    shutil.copy2(package_path, user_path)
    html_d[LOGO_PATH] = os.path.join(OUTPUT_DIR, "crispector_logo.jpg")

    return html_d

#####----------------------#####
#####----modifications-----#####
#####----------------------#####
def plot_modification_tables(mod_table: ModificationTables, modifications: ModificationTypes, edit_table: IsEdit,
                             table_offset: int, table_indexes: List[int], output: Path, html_d: Dict, base_output: Path,
                             name_suffix = ""):
    """
    Plot all modification tables around cut-site.
    Also display edit events.
    :param mod_table: ModificationTables
    :param modifications: ModificationTypes
    :param edit_table: edit as marked by crispector algorithm
    :param table_offset: reference offset (in bp) to the start of the qualification window
    :param output: output path
    :param table_indexes - list indexes for table plot
    :param name_suffix - prefix for file name
    :param figsize
    :return:
    """
    # Set font
    mpl.rcParams.update(mpl.rcParamsDefault)
    mpl.rcParams['font.size'] = 9
    mpl.rcParams['xtick.labelsize'] = 11
    mpl.rcParams['ytick.labelsize'] = 11
    mpl.rcParams['axes.labelsize'] = 11
    mpl.rcParams['legend.fontsize'] = 11
    mpl.rcParams['axes.titlesize'] = 11
    indel_size = 8
    bar_width = 0.4
    dpi = 300
    mock_color = '#3690c0'  # blue
    tx_color = '#f16a13'  # orange
    edit_color = '#dcdcdc'  # light grey
    cut_site_color = "red"
    grid_color = 'grey'
    # Create axes
    fig_w, fig_h = (8, 4.5)
    fig, axes = plt.subplots(nrows=len(table_indexes), ncols=1, figsize=(fig_w, fig_h), sharex=True)

    for axes_idx, table_idx in enumerate(table_indexes):

        # Insertions have a slightly different positions marking.
        is_ins = modifications.types[table_idx] == IndelType.INS

        # Positions and bars positions
        edit = edit_table[table_idx]
        positions = np.arange(len(edit)) - 0.5 * is_ins
        bar_ind = positions + 0.5
        cut_site = len(positions) / 2

        tx = mod_table.tables[table_idx][C_TX, table_offset:table_offset + len(positions)]
        mock = mod_table.tables[table_idx][C_MOCK, table_offset:table_offset + len(positions)]
        y_max = max(max(tx.max(), mock.max()), 1)

        # Set green color for Edit events
        for pos_idx, pos in enumerate(positions):
            if edit[pos_idx]:
                axes[axes_idx].axvspan(pos, pos + 1, facecolor=edit_color)

        # Set bp "grid"
        grid_positions = positions if is_ins else np.arange(len(edit) + 1)
        axes[axes_idx].vlines(grid_positions, ymin=0, ymax=y_max, color=grid_color, lw=0.5)
        axes[axes_idx].spines['bottom'].set_color(grid_color)
        axes[axes_idx].spines['top'].set_color(grid_color)
        axes[axes_idx].spines['right'].set_color(grid_color)
        axes[axes_idx].spines['left'].set_color(grid_color)

        # Set cut-site in red (Insertions cut site is a full bp cell)
        if is_ins:
            axes[axes_idx].vlines([cut_site - 1, cut_site], ymin=0, ymax=y_max, color=cut_site_color, lw=1)
            axes[axes_idx].hlines([0, y_max], xmin=cut_site - 1, xmax=cut_site, color=cut_site_color, lw=2)
        else:
            axes[axes_idx].axvline(cut_site, ymin=0, ymax=y_max, color=cut_site_color, lw=1)

        # Create bar plot
        axes[axes_idx].bar(bar_ind - bar_width / 2, tx, width=bar_width, color=tx_color,
                           label="Tx indels (out of {:,} Reads)".format(mod_table.n_reads_tx))
        axes[axes_idx].bar(bar_ind + bar_width / 2, mock, width=bar_width, color=mock_color,
                           label="Mock indels (out of {:,} Reads)".format(mod_table.n_reads_mock))

        # Set x, y lim & ticks
        axes[axes_idx].set_xlim(min(positions) - 0.5, max(positions) + 1.5)
        axes[axes_idx].set_xticks([])
        axes[axes_idx].set_yticks([])
        axes[axes_idx].set_ylim(0, y_max)

        # Set y_label - Use text due to alignment issues
        axes[axes_idx].text(x=positions[0] - 1, y=y_max / 2,
                             s=modifications.plot_name_at_idx(table_idx), ha="right", va="center")

        # Create legend at the bottom of the plot
        if axes_idx == len(table_indexes) - 1:
            axes[axes_idx].bar([0], [0], color=edit_color, label="Edit event", edgecolor='grey')
            axes[axes_idx].plot([], [], color=cut_site_color, label="Cut-Site")
            handles, labels = axes[axes_idx].get_legend_handles_labels()
            order = [1, 2, 3, 0]
            axes[axes_idx].legend([handles[idx] for idx in order], [labels[idx] for idx in order],
                                   bbox_to_anchor=(0.75, 0))

        # Add the reference sequence in the position of the title
        if axes_idx == 0:
            ref = mod_table.amplicon[table_offset:table_offset + len(positions)]
            offset = 0.5 if name_suffix == "ins" else 0
            for pos_idx, ii in enumerate(bar_ind):
                axes[axes_idx].text(ii+offset, y_max, ref[pos_idx], ha="center", va="bottom", weight='bold', size=15)
            # red cute-site
            axes[axes_idx].text(cut_site-offset, 1.1 * y_max, "|", ha="center", va="bottom", weight='bold', size=10,
                                color=cut_site_color)

        # Add bars values (numbers) as text
        for pos_idx, bar_i in enumerate(bar_ind):
            if tx[pos_idx]:
                axes[axes_idx].text(bar_i - bar_width / 2, 0.05 * y_max, "{:,}".format(tx[pos_idx]), rotation=90,
                                    size=indel_size, ha="center", va="bottom")
            if mock[pos_idx]:
                axes[axes_idx].text(bar_i + bar_width / 2, 0.05 * y_max, "{:,}".format(mock[pos_idx]), rotation=90,
                                    size=indel_size, ha="center", va="bottom")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        fig.savefig(os.path.join(output, 'classifier_results_by_position_{}.png'.format(name_suffix)),
                    bbox_inches='tight', dpi=dpi)
        plt.subplots_adjust(left=0.14, bottom=0.25, right=None, top=0.9)
        fig.savefig(os.path.join(output, 'classifier_results_by_position_{}.svg'.format(name_suffix)),
                    box_inches='tight')

        plt.close(fig)

    if name_suffix == "all":
        return

    if name_suffix == "ins":
        tab_name = CLS_RES_INS
        title = "Insertions"
    elif name_suffix == "del":
        tab_name = CLS_RES_DEL
        title = "Deletions"
    else:
        tab_name = CLS_RES_MIX
        title = "Mixed & Substitutions"

    html_d[CLS_RES_SECTION][tab_name] = dict()
    html_d[CLS_RES_SECTION][tab_name][PLOT_PATH] = os.path.join(base_output,'classifier_results_by'
                                                                            '_position_{}.png'.format(name_suffix))
    html_d[CLS_RES_SECTION][tab_name][PDF_PATH] = os.path.join(base_output, 'classifier_results_by'
                                                                            '_position_{}.svg'.format(name_suffix))
    html_d[CLS_RES_SECTION][tab_name][TITLE] = title
    html_d[CLS_RES_SECTION][tab_name][W] = dpi * fig_w
    html_d[CLS_RES_SECTION][tab_name][H] = dpi * fig_h


#####----------------------#####
#####---------Site---------#####
#####----------------------#####
def plot_distribution_of_all_modifications(tables: ModificationTables, cut_site: int, win_size: int,
                                           output: Path, html_d: Dict, base_path: Path):

    tx_dist_d = dict()
    mock_dist_d = dict()
    amplicon_length = len(tables.amplicon)
    for indel_type in [IndelType.DEL, IndelType.SUB, IndelType.INS]:
        tx_dist_d[indel_type] = np.zeros(amplicon_length + 1, dtype=np.int)
        mock_dist_d[indel_type] = np.zeros(amplicon_length + 1, dtype=np.int)

    # Tx - Aggregate modifications
    for row_idx, row in tables.tx_dist.iterrows():
        mod_range = range(row[POS_IDX_S], row[POS_IDX_E])
        tx_dist_d[row[INDEL_TYPE]][mod_range] += row[FREQ]

    # Mock - Aggregate modifications
    for row_idx, row in tables.mock_dist.iterrows():
        mod_range = range(row[POS_IDX_S], row[POS_IDX_E])
        mock_dist_d[row[INDEL_TYPE]][mod_range] += row[FREQ]

    # Find plots heights
    tx_max_indels = np.array([tx_dist_d[IndelType.INS], tx_dist_d[IndelType.DEL], tx_dist_d[IndelType.SUB]]).max()
    mock_max_indels = np.array(
        [mock_dist_d[IndelType.INS], mock_dist_d[IndelType.DEL], mock_dist_d[IndelType.SUB]]).max()
    max_indels = max(tx_max_indels, mock_max_indels)

    # Set font
    mpl.rcParams.update(mpl.rcParamsDefault)
    mpl.rcParams['font.size'] = 9
    mpl.rcParams['axes.labelsize'] = 11
    dpi = 300

    # Create axes
    fig_w, fig_h = 8, 4.5
    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(fig_w, fig_h), sharex=True)
    plt.subplots_adjust(wspace=0, hspace=0)

    # Plot distributions
    for exp_type in ExpType:
        dist_d = tx_dist_d if exp_type == ExpType.TX else mock_dist_d
        ax = axes[0] if exp_type == ExpType.TX else axes[1]
        positions = list(range(amplicon_length + 1))
        ax.axvline(x=cut_site, linestyle='-', color="red", label='Expected cut-site', linewidth=1)
        ax.axvline(x=cut_site - win_size, linestyle='--', color='k', label='Quantification window', linewidth=0.5)
        ax.axvline(x=cut_site + win_size, linestyle='--', color='k', linewidth=0.5)
        ax.plot(positions, dist_d[IndelType.INS], color=IndelType.INS.color, label="Insertions", linewidth=1.5, alpha=0.9)
        ax.plot(positions, dist_d[IndelType.SUB], color=IndelType.SUB.color, label="Substitutions", linewidth=1.5, alpha=0.9)
        ax.plot(positions, dist_d[IndelType.DEL], color=IndelType.DEL.color, label="Deletions", linewidth=1.5, alpha=0.9)
        ax.set_ylim(bottom=0, top=max(int(1.1 * max_indels), 10))
        ax.set_xlim(left=0, right=amplicon_length)
        ax.set_ylabel("{}\nIndel count".format(exp_type.name))

    plt.xlabel("Position")

    # Create legend, title and x-label
    axes[0].legend()

    # Remove Tx first y axe label (overlap with mock)
    plt.setp(axes[0].get_yticklabels()[0], visible=False)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        fig.savefig(os.path.join(output, 'distribution_of_all_modifications.png'), bbox_inches='tight', dpi=dpi)
        fig.savefig(os.path.join(output, 'distribution_of_all_modifications.svg'), pad_inches = 1, box_inches='tight')
        plt.close(fig)

    html_d[MOD_SECTION][MOD_DIST] = dict()
    html_d[MOD_SECTION][MOD_DIST][PLOT_PATH] = os.path.join(base_path, 'distribution_of_all_modifications.png')
    html_d[MOD_SECTION][MOD_DIST][PDF_PATH] = os.path.join(base_path, 'distribution_of_all_modifications.svg')
    html_d[MOD_SECTION][MOD_DIST][TITLE] = "All Indels Distribution"
    html_d[MOD_SECTION][MOD_DIST][W] = dpi*fig_w
    html_d[MOD_SECTION][MOD_DIST][H] = dpi*fig_h

def plot_distribution_of_edit_events(tables: ModificationTables, cut_site: int, win_size: int,
                                     output: Path, html_d: Dict, base_path: Path):
    dist_d = dict()
    amplicon_length = len(tables.amplicon)
    for indel_type in [IndelType.DEL, IndelType.SUB, IndelType.INS]:
        dist_d[indel_type] = np.zeros(amplicon_length + 1, dtype=np.int)

    # Aggregate modifications
    edit_mod = tables.tx_dist.loc[tables.tx_dist[IS_EDIT]]
    for row_idx, row in edit_mod.iterrows():
        mod_range = range(row[POS_IDX_S], row[POS_IDX_E])
        dist_d[row[INDEL_TYPE]][mod_range] += row[FREQ]

    # Set font
    mpl.rcParams.update(mpl.rcParamsDefault)
    mpl.rcParams['font.size'] = 9
    mpl.rcParams['axes.labelsize'] = 11
    dpi = 300
    fig_w, fig_h = 8, 4.5

    # Create axes
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(fig_w, fig_h))

    positions = list(range(amplicon_length + 1))
    ax.axvline(x=cut_site, linestyle='-', color="red", label='Expected cut-site', linewidth=1)
    ax.axvline(x=cut_site - win_size, linestyle='--', color='k', label='Quantification window', linewidth=0.5)
    ax.axvline(x=cut_site + win_size, linestyle='--', color='k', linewidth=0.5)

    for indel_type in [IndelType.INS, IndelType.SUB, IndelType.DEL]:
        ax.plot(positions, dist_d[indel_type], color=indel_type.color, label=indel_type.name, linewidth=1.5, alpha=0.9)

    # Create legend, axes limit and labels
    ax.legend()
    max_indels = np.array([dist_d[IndelType.INS], dist_d[IndelType.DEL], dist_d[IndelType.SUB]]).max()
    ax.set_ylim(bottom=0, top=max(int(1.1 * max_indels), 10))
    ax.set_xlim(left=0, right=amplicon_length)
    ax.set_xlabel("Position")
    ax.set_ylabel("Indel count")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        fig.savefig(os.path.join(output, 'distribution_of_edit_events.png'), bbox_inches='tight', dpi=dpi)
        fig.savefig(os.path.join(output, 'distribution_of_edit_events.svg'), box_inches='tight')
        plt.close(fig)

    html_d[MOD_SECTION][EDIT_DIST] = dict()
    html_d[MOD_SECTION][EDIT_DIST][PLOT_PATH] = os.path.join(base_path, 'distribution_of_edit_events.png')
    html_d[MOD_SECTION][EDIT_DIST][PDF_PATH] = os.path.join(base_path, 'distribution_of_edit_events.svg')
    html_d[MOD_SECTION][EDIT_DIST][TITLE] = "Edited Events Distribution"
    html_d[MOD_SECTION][EDIT_DIST][W] = dpi*fig_w
    html_d[MOD_SECTION][EDIT_DIST][H] = dpi*fig_h


def plot_distribution_of_edit_event_sizes(tables: ModificationTables, output: Path, html_d: Dict, base_path: Path):
    # Set font
    mpl.rcParams.update(mpl.rcParamsDefault)
    mpl.rcParams['font.size'] = 9
    mpl.rcParams['xtick.labelsize'] = 9
    mpl.rcParams['ytick.labelsize'] = 11
    mpl.rcParams['axes.labelsize'] = 11
    mpl.rcParams['axes.titlesize'] = 11
    mpl.rcParams['legend.fontsize'] = 11
    bar_width = 1
    dpi = 300
    amplicon_length = len(tables.amplicon)

    # Create bin borders & names:
    bins = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 20, 30, 40, 50, 100, amplicon_length]
    bins_names = []
    for l, r in zip(bins[:-1], bins[1:]):
        if r - l == 1:
            bins_names.append(str(r))
        else:
            bins_names.append("{}-{}".format(l + 1, r))
    bin_num = len(bins) - 1

    # Get indel size values
    values = {}
    edited_mod = tables.tx_dist.loc[tables.tx_dist[IS_EDIT]]
    for indel in [IndelType.DEL, IndelType.INS, IndelType.SUB]:
        indel_dist = edited_mod.loc[edited_mod['IndelType'] == indel]
        df_bins = pd.cut(indel_dist['indel_length'], bins)
        values[indel] = indel_dist.groupby(df_bins)[FREQ].agg(['sum'])['sum'].values

    # Create the X position of bars
    r = list(range(0, bin_num))

    # Create figure
    fig_w, fig_h = 8, 3
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(fig_w, fig_h), constrained_layout=True)

    # Create bar plot
    for idx, indel in enumerate([IndelType.DEL, IndelType.INS, IndelType.SUB]):
        axes[idx].bar(r, values[indel], label=indel.name, color=indel.color, edgecolor='black', alpha=0.9,
                      width=bar_width)
        axes[idx].set_xticks(r)
        axes[idx].set_xticklabels(bins_names, rotation='vertical')
        axes[idx].spines['top'].set_visible(False)
        axes[idx].spines['right'].set_visible(False)

        # Labels, ticks and title
        axes[idx].legend(loc='upper right')

    axes[0].set_ylabel('Number of modifications')
    axes[1].set_xlabel('Edit event length')

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        fig.savefig(os.path.join(output, 'distribution_of_edit_events_size.png'), bbox_inches='tight', dpi=dpi)
        fig.savefig(os.path.join(output, 'distribution_of_edit_events_size.svg'), pad_inches = 1, box_inches='tight')
        plt.close(fig)

    html_d[MOD_SECTION][EDIT_SIZE_DIST] = dict()
    html_d[MOD_SECTION][EDIT_SIZE_DIST][PLOT_PATH] = os.path.join(base_path, 'distribution_of_edit_events_size.png')
    html_d[MOD_SECTION][EDIT_SIZE_DIST][PDF_PATH] = os.path.join(base_path, 'distribution_of_edit_events_size.svg')
    html_d[MOD_SECTION][EDIT_SIZE_DIST][TITLE] = "Edit Event Size Distribution"
    html_d[MOD_SECTION][EDIT_SIZE_DIST][W] = dpi*fig_w
    html_d[MOD_SECTION][EDIT_SIZE_DIST][H] = dpi*fig_h


def plot_edited_reads_to_table(tables: ModificationTables, cut_site: int, output: Path, html_d: Dict,
                               base_path: Path):
    """
    Create an HTML page with all the edited reads (aligned to reference)
    """
    length = READ_LEN_SIDE

    # Columns variables
    column_order = list(range(1, length + 1)) + list(range(length + 2, 2 * length + 4)) + [length + 1]
    column_width = (2 * length + 3) * [40]
    column_width[-1] = 5  # cut-site marking
    column_width[-2] = 160  # Frequency
    column_width[-3] = 160  # number of reads

    # Header
    reference_l = [base for base in tables.amplicon[cut_site - length:cut_site + length]]
    header = []
    header += ["<b>{}<b>".format(bp) for bp in reference_l]
    header += ["<b>Number Of Reads<b>", "<b>Frequency (%)<b>", ""]
    header_colors = color_edit_background(reference_l, reference_l) + ['#c8d4e3', '#c8d4e3', 'black']

    # Prepare edit reads - Cut the relevant part and aggregate
    edit_reads = tables.tx_reads.loc[tables.tx_reads[IS_EDIT]]

    partial_reads = pd.DataFrame(columns=[FREQ, ALIGNMENT_W_INS, ALIGNMENT_W_DEL])
    if edit_reads.shape[0] > 0:
        partial_reads[FREQ] = edit_reads[FREQ].copy()
        partial_reads[ALIGNMENT_W_INS] = edit_reads.apply(lambda row: get_read_around_cut_site(row[ALIGNMENT_W_INS],
                                                                                               row[ALIGN_CUT_SITE],
                                                                                               length), axis=1)
        partial_reads[ALIGNMENT_W_DEL] = edit_reads.apply(lambda row: get_read_around_cut_site(row[ALIGNMENT_W_DEL],
                                                                                               row[ALIGN_CUT_SITE],
                                                                                               length), axis=1)
        partial_reads = partial_reads.groupby([ALIGNMENT_W_INS, ALIGNMENT_W_DEL])[FREQ].sum().to_frame(
            FREQ).reset_index()
        partial_reads = partial_reads.sort_values(by=[FREQ], ascending=False).reset_index(drop=True)

    # Cells values & colors
    values = np.zeros(shape=(partial_reads.shape[0], 2 * length + 3), dtype='<U20')
    colors = np.zeros(shape=(partial_reads.shape[0], 2 * length + 3), dtype='<U20')
    for row_idx, (_, row) in enumerate(partial_reads.iterrows()):
        read_w_ins = list(row[ALIGNMENT_W_INS])
        read_w_del = list(row[ALIGNMENT_W_DEL])
        # cell fill color
        colors[row_idx] = np.array(color_edit_background(read_w_ins, read_w_del) + ['lightgrey', 'lightgrey', 'black'])
        values[row_idx] = read_w_del + ["<b>{:,}<b>".format(row[FREQ]),
                                        "<b>{:.4f}<b>".format(100 * row[FREQ] / tables.n_reads_tx), ""]

    # convert to list
    values = [list(values[:, idx]) for idx in range(values.shape[1])]
    colors = [list(colors[:, idx]) for idx in range(colors.shape[1])]

    # Make the plot
    fig = make_subplots(
        rows=2, cols=1,
        row_heights=[0.99, 0.01],
        vertical_spacing=0.03,
        specs=[[{"type": "table"}],
               [{"type": "scatter"}]]
    )
    # Add the table
    fig.add_trace(
        go.Table(
            columnorder=column_order,
            columnwidth=column_width,
            header=dict(
                values=header,
                font=dict(size=14),
                fill=dict(color=header_colors),
                align="center"
            ),
            cells=dict(
                values=values,
                align="center",
                fill=dict(color=colors),
                font=dict(size=14, color='black'))
        ),
        row=1, col=1
    )

    # Add Scatter plots for legend
    fig.add_trace(
        go.Scatter(
            x=[0],
            y=[0],
            mode="lines",
            name="Insertions",
            line=dict(color='indianred')
        ),
        row=2, col=1
    )
    fig.add_trace(
        go.Scatter(
            x=[0],
            y=[0],
            mode="lines",
            name="Deletions",
            line=dict(color='lightgrey')
        ),
        row=2, col=1
    )
    fig.add_trace(
        go.Scatter(
            x=[0],
            y=[0],
            mode="lines",
            name="Substitutions",
            line=dict(color='royalblue')
        ),
        row=2, col=1
    )
    fig.add_trace(
        go.Scatter(
            x=[0],
            y=[0],
            mode="lines",
            name="Cut Site",
            line=dict(color='Black')
        ),
        row=2, col=1
    )

    fig.update_layout(
        title_text="Edited reads around the cut-site",
        title_x=0.5,
        titlefont=dict(size=30, color='black', family='Arial, sans-serif'),
        showlegend=True,
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        legend=go.layout.Legend(
            traceorder="normal",
            font=dict(
                family="sans-serif",
                size=20,
                color="black"
            ),
            bordercolor="Black",
            borderwidth=2
        )
    )
    fig.update_xaxes(color='rgba(0,0,0,0)')
    fig.update_yaxes(color='rgba(0,0,0,0)')

    # Dump to static HTML
    fig.write_html(os.path.join(output, 'edited_reads_table.html'))
    html_d[READ_SECTION][READ_EDIT] = os.path.join(base_path, 'edited_reads_table.html')


def plot_site_editing_activity(algorithm: CoreAlgorithm, result_d: Dict, site_name: str, output: Path, html_d: Dict,
                               base_path: Path):
    # Set font
    mpl.rcParams.update(mpl.rcParamsDefault)
    sns.set(style="whitegrid")
    mpl.rcParams['font.size'] = 9
    mpl.rcParams['xtick.labelsize'] = 11
    mpl.rcParams['ytick.labelsize'] = 11
    mpl.rcParams['axes.labelsize'] = 11
    mpl.rcParams['axes.titlesize'] = 11
    dpi = 300
    bar_width = 0.8

    if algorithm.is_on_target:
        bar_color = ON_TARGET_COLOR
    else:
        bar_color = OFF_TARGET_COLOR

    # Define fix and axes
    fig_w, fig_h = 4, 3
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(fig_w, fig_h))
    plt.subplots_adjust(left=0.2)

    # Get bar data
    editing = result_d[EDIT_PERCENT]
    CI_high = result_d[CI_HIGH] - editing
    CI_low = editing - result_d[CI_LOW]

    # plot bar
    ax.bar([0], [result_d[EDIT_PERCENT]], color=bar_color, width=bar_width,
           yerr=[[CI_low], [CI_high]], align='center', ecolor='black', capsize=7)

    # Set labels
    ax.set_ylabel("Editing Activity (%)")

    # Set scale and lim
    y_lim = max(min(1.5 * (editing + CI_high), 100), 0.1)
    ax.set_ylim(0, y_lim)
    ax.set_xlim(-1.2, 1.2)

    # Text below each bar plot + y ticks
    ax.set_xticks([0])
    ax.set_xticklabels([site_name])

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        fig.savefig(os.path.join(output, 'site_editing_activity.png'), bbox_inches='tight', dpi=dpi)
        fig.savefig(os.path.join(output, 'site_editing_activity.svg'), box_inches='tight')
        plt.close(fig)

    edit_text = "<p style=\"text-align: center\"> Number of edited reads - {:,} (out of {:,} reads).<p>".format(result_d[TX_EDIT], result_d[TX_READ_NUM])
    edit_text += "<p style=\"text-align: center\"> Editing activity - {:.2f}%, CI=({:.2f}% - {:.2f}%).<p>".format(editing, result_d[CI_LOW], result_d[CI_HIGH])
    html_d[EDITING_ACTIVITY] = dict()
    html_d[EDITING_ACTIVITY][PLOT_PATH] = os.path.join(base_path, 'site_editing_activity.png')
    html_d[EDITING_ACTIVITY][PDF_PATH] = os.path.join(base_path, 'site_editing_activity.svg')
    html_d[EDITING_ACTIVITY][TITLE] = "Editing Activity"
    html_d[EDITING_ACTIVITY][W] = dpi*fig_w
    html_d[EDITING_ACTIVITY][H] = dpi*fig_h
    html_d[EDITING_ACTIVITY][EDIT_TEXT] = edit_text


#####----------------------#####
#####------Experiment------#####
#####----------------------#####
def plot_editing_activity(result_df: AlgResultDf, confidence_interval: float, editing_threshold: float, html_d: Dict,
                          output: Path):

    # Set font
    mpl.rcParams.update(mpl.rcParamsDefault)
    mpl.rcParams['font.size'] = 9
    mpl.rcParams['xtick.labelsize'] = 9
    mpl.rcParams['ytick.labelsize'] = 9
    mpl.rcParams['axes.labelsize'] = 11
    mpl.rcParams['axes.titlesize'] = 11
    mpl.rcParams['legend.fontsize'] = 11
    editing_bar_text_size = 8
    dpi = 300
    title = ""

    # Filter all low editing activity sites
    edit_df = result_df.dropna()
    edit_df = edit_df.loc[edit_df[CI_HIGH] >= editing_threshold]

    # Sort experiments
    edit_df = edit_df.sort_values(by=EDIT_PERCENT, ascending=False)

    # Create axes
    max_bars = 20
    bar_num = edit_df.shape[0]
    plot_num = math.ceil(bar_num / max_bars)
    # set dynamic bar_width - according to the number of bars
    fig_w = 4 + 4 * (min(bar_num, max_bars) / 20)
    fig_h = plot_num * 3
    bar_width = 0.9
    fig, axes = plt.subplots(nrows=plot_num, ncols=1, figsize=(fig_w, fig_h), tight_layout=True)

    # Create bars and bar names
    editing = edit_df[EDIT_PERCENT].values
    site_names = edit_df[SITE_NAME].values
    CI_high = edit_df[CI_HIGH].values - editing
    CI_low = editing - edit_df[CI_LOW].values
    on_target = edit_df[ON_TARGET].values

    # Create bar plot
    for idx in range(plot_num):
        if plot_num == 1:
            axes = [axes]

        plt_editing = editing[max_bars * idx:min(max_bars * (idx + 1), len(editing))]
        plt_site_names = site_names[max_bars * idx:min(max_bars * (idx + 1), len(site_names))]
        plt_CI_high = CI_high[max_bars * idx:min(max_bars * (idx + 1), len(CI_high))]
        plt_CI_low = CI_low[max_bars * idx:min(max_bars * (idx + 1), len(CI_low))]
        plt_on_target = on_target[max_bars * idx:min(max_bars * (idx + 1), len(on_target))]

        # The X position of bars
        number_of_bars = len(plt_editing)
        bar_pos = list(range(1, number_of_bars + 1))

        # Bar plot
        capsize = 2
        bar_plot = axes[idx].bar(bar_pos, plt_editing, width=bar_width, color=OFF_TARGET_COLOR,
                                 yerr=[plt_CI_low, plt_CI_high], align='center', ecolor='black',
                                 capsize=capsize, error_kw={"elinewidth": 0.75})
        for site_idx, is_site_on_target in enumerate(plt_on_target):
            if is_site_on_target:
                bar_plot[site_idx].set_color(ON_TARGET_COLOR)

        # Add horizontal line
        axes[idx].axhline(y=editing_threshold, linewidth=1, color='k', linestyle="--", alpha=0.5)

        # Set labels
        axes[idx].set_xlabel("Site Name")
        axes[idx].set_ylabel("Editing Activity (%)")

        # Set scale and lim
        y_lim = 1e-2
        axes[idx].set_ylim(1e-2, 100)
        axes[idx].set_yscale('log')
        if plot_num > 1:
            axes[idx].set_xlim(0, max_bars + 1)
        else:
            axes[idx].set_xlim(0, min(number_of_bars + 3, max_bars + 1))

        # Text below each bar plot + y ticks
        axes[idx].set_xticks([r + 1 for r in range(number_of_bars)])
        axes[idx].set_xticklabels(plt_site_names, rotation='vertical')

        def format_func(value, tick_number):
            if value == 0.001:
                return "0.001"
            elif value == 0.01:
                return "0.01"
            elif value == 0.1:
                return "0.1"
            elif value == 1.0:
                return "1"
            elif value == 10.0:
                return "10"
            elif value == 100.0:
                return "100"
            elif value == 1000.0:
                return "1000"

        axes[idx].yaxis.set_major_formatter(plt.FuncFormatter(format_func))
        axes[idx].spines['top'].set_visible(False)
        axes[idx].spines['right'].set_visible(False)

        # Text on the top of each bar plot
        text_height = 1.1 * (plt_editing + plt_CI_high)
        for text_idx in range(number_of_bars):
            axes[idx].text(x=bar_pos[text_idx], y=text_height[text_idx],
                           s="{:.2f}".format(plt_editing[text_idx]), ha='center', va='bottom',
                           size=editing_bar_text_size)

        if idx == 0:
            title = "Editing Activity with {} % CI above {}%".format(100 * confidence_interval, editing_threshold)
            # axes[idx].set_title(title, weight='bold', family='serif')

            # Add legend
            axes[idx].bar([0], [y_lim], color=ON_TARGET_COLOR, label="On-Target")
            axes[idx].bar([0], [y_lim], color=OFF_TARGET_COLOR, label="Off-Target")
            axes[idx].legend(loc='upper right')

    if edit_df.shape[0] > 0:
        fig.savefig(os.path.join(output, 'editing_activity.png'), box_inches='tight', dpi=dpi)
        fig.savefig(os.path.join(output, 'editing_activity.svg'), pad_inches = 1, box_inches='tight')

    html_d[EDIT_SECTION][EDITING_ACTIVITY] = dict()
    html_d[EDIT_SECTION][EDITING_ACTIVITY][PLOT_PATH] = os.path.join(OUTPUT_DIR, 'editing_activity.png')
    html_d[EDIT_SECTION][EDITING_ACTIVITY][PDF_PATH] = os.path.join(OUTPUT_DIR, 'editing_activity.svg')
    html_d[EDIT_SECTION][EDITING_ACTIVITY][TITLE] = title
    html_d[EDIT_SECTION][EDITING_ACTIVITY][W] = "{}%".format(min(95, int(50 + 45 * (bar_num / 20))))
    html_d[EDIT_SECTION][EDITING_ACTIVITY][H] = dpi*fig_h

    plt.close(fig)


def create_reads_statistics_report(result_df: AlgResultDf, tx_in: int, tx_merged: int, tx_aligned: int,
                                   mock_in: int, mock_merged: int, mock_aligned: int, html_d: Dict, output: Path):
    html_d[READING_STATS] = dict()
    html_d[READING_STATS][TITLE] = "Reading Statistics"
    # Set font
    mpl.rcParams.update(mpl.rcParamsDefault)
    sns.set(style="whitegrid")
    mpl.rcParams['font.size'] = 9
    mpl.rcParams['ytick.labelsize'] = 11
    mpl.rcParams['xtick.labelsize'] = 11
    mpl.rcParams['axes.labelsize'] = 11
    mpl.rcParams['axes.titlesize'] = 11
    mpl.rcParams['legend.fontsize'] = 11
    dpi = 300
    fig_w, fig_h = 8, 4.6
    bar_width = 0.4
    mock_color = '#3690c0'  # blue
    mock_color_lighter = '#72b1d2'
    mock_color_lightest = '#9ac7df'
    tx_color = '#f16a13'  # orange
    tx_color_lighter = '#f59659'  # orange
    tx_color_lightest = '#f7b488'  # orange

    # set input == merged if input information isn't available
    if tx_in == 0:
        if tx_merged == 0:
            tx_in, tx_merged = 1, 1
        else:
            tx_in = tx_merged
    if mock_in == 0:
        if mock_merged == 0:
            mock_in, mock_merged = 1, 1
        else:
            mock_in = mock_merged

    # Create mapping statistics plot
    fig = plt.figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes([0, 0, 1, 1])

    bar_ind = np.arange(6)
    bars = np.array([tx_in, tx_merged, tx_aligned, mock_in, mock_merged, mock_aligned])
    colors = [tx_color, tx_color_lighter, tx_color_lightest, mock_color, mock_color_lighter, mock_color_lightest]

    # Create bar plot
    ax.bar(bar_ind, bars, width=bar_width, color=colors)

    # Add numbers above bars
    text_height = 1.01 * bars
    ax.text(x=bar_ind[0], y=text_height[0],
                 s="{:,}".format(tx_in), ha='center', va='bottom')
    ax.text(x=bar_ind[1], y=text_height[1],
                 s="{:,}\n({:.2f}%)".format(tx_merged, 100 * tx_merged / tx_in), ha='center', va='bottom')
    ax.text(x=bar_ind[2], y=text_height[2],
                 s="{:,}\n({:.2f}%)".format(tx_aligned, 100 * tx_aligned / tx_in), ha='center', va='bottom')
    ax.text(x=bar_ind[3], y=text_height[3],
                 s="{:,}".format(mock_in), ha='center', va='bottom')
    ax.text(x=bar_ind[4], y=text_height[4],
                 s="{:,}\n({:.2f}%)".format(mock_merged, 100 * mock_merged / mock_in), ha='center', va='bottom')
    ax.text(x=bar_ind[5], y=text_height[5],
                 s="{:,}\n({:.2f}%)".format(mock_aligned, 100 * mock_aligned / mock_in), ha='center', va='bottom')

    # Set x, y lim & ticks and title
    ax.set_xlim(min(bar_ind) - 0.5, max(bar_ind) + 0.5)
    ax.set_xticks(bar_ind)
    ax.set_xticklabels(['Treatment\nInput', 'Treatment\nMerged', 'Treatment\nAligned',
                             'Mock\nInput', 'Mock\nMerged', 'Mock\nAligned'])
    ax.set_ylim(0, 1.2 * np.max(bars))
    ax.set_ylabel("Number Of Reads")
    title = "Mapping Statistics"
    # ax.set_title(title, weight='bold', family='serif')

    html_d[READING_STATS][MAPPING_STATS] = dict()
    html_d[READING_STATS][MAPPING_STATS][PLOT_PATH] = os.path.join(OUTPUT_DIR, 'mapping_statistics.png')
    html_d[READING_STATS][MAPPING_STATS][PDF_PATH] = os.path.join(OUTPUT_DIR, 'mapping_statistics.svg')
    html_d[READING_STATS][MAPPING_STATS][TITLE] = title
    html_d[READING_STATS][MAPPING_STATS][W] = dpi*fig_w
    html_d[READING_STATS][MAPPING_STATS][H] = dpi*fig_h

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        fig.savefig(os.path.join(output, 'mapping_statistics.png'), bbox_inches='tight', dpi=dpi)
        fig.savefig(os.path.join(output, 'mapping_statistics.svg'), bbox_inches='tight', pad_inches=1)
        plt.close(fig)

    # Create reads box_plot
    fig = plt.figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes([0, 0, 1, 1])

    max_number_of_reads = max(result_df[TX_READ_NUM].max(), result_df[MOCK_READ_NUM].max())
    sns.scatterplot(x=TX_READ_NUM, y=MOCK_READ_NUM, data=result_df, ax=ax, s=115)
    plt.plot([0, max_number_of_reads], [0, max_number_of_reads], c='k') # balance line

    unbalanced_df = result_df.loc[(result_df[TX_READ_NUM] > UNBALANCED_READ_WARNING * result_df[MOCK_READ_NUM]) |
                                  (result_df[MOCK_READ_NUM] > UNBALANCED_READ_WARNING * result_df[TX_READ_NUM])]
    if unbalanced_df.shape[0] > 0:
        sns.scatterplot(x=TX_READ_NUM, y=MOCK_READ_NUM, data=unbalanced_df, ax=ax, s=115, color='r',
                        label="Sites with unbalanced read numbers")
        plt.legend()
    ax.set_ylabel("Number of reads in Mock")
    ax.set_xlabel("Number of reads in Treatment")
    title = "Number of Aligned Reads per Site"

    logger = LoggerWrapper.get_logger()
    if unbalanced_df.shape[0] > 0:
        logger.warning("Experiment has sites with unbalanced number reads between Treatment and Mock."
                       " These sites editing activity estimation can be inaccurate.")
    for _, row in unbalanced_df.iterrows():
        logger.warning("Site {} has highly unbalanced number of reads: Treatment - {:,}, Mock - {:,}.".format(row[SITE_NAME],
                                                                                                              row[TX_READ_NUM],
                                                                                                              row[MOCK_READ_NUM]))

    html_d[READING_STATS][MAPPING_PER_SITE] = dict()
    html_d[READING_STATS][MAPPING_PER_SITE][PLOT_PATH] = os.path.join(OUTPUT_DIR, 'number_of_aligned_reads_per_site.png')
    html_d[READING_STATS][MAPPING_PER_SITE][PDF_PATH] = os.path.join(OUTPUT_DIR, 'number_of_aligned_reads_per_site.svg')
    html_d[READING_STATS][MAPPING_PER_SITE][TITLE] = title
    html_d[READING_STATS][MAPPING_PER_SITE][W] = dpi*fig_w
    html_d[READING_STATS][MAPPING_PER_SITE][H] = dpi*fig_h

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        fig.savefig(os.path.join(output, 'number_of_aligned_reads_per_site.png'), bbox_inches='tight', dpi=dpi)
        fig.savefig(os.path.join(output, 'number_of_aligned_reads_per_site.svg'), bbox_inches='tight', pad_inches=1)
        plt.close(fig)


def plot_translocations_heatmap(result_df: pd.DataFrame, trans_result_df: TransResultDf, translocation_p_value: float,
                                html_d: Dict, output: Path):

    # Set font
    mpl.rcParams.update(mpl.rcParamsDefault)
    mpl.rcParams['font.size'] = 9
    mpl.rcParams['xtick.labelsize'] = 9
    mpl.rcParams['ytick.labelsize'] = 9
    mpl.rcParams['axes.labelsize'] = 11
    mpl.rcParams['axes.titlesize'] = 11
    mpl.rcParams['legend.fontsize'] = 11

    # Prepare Heat Map
    trans_df = trans_result_df.loc[trans_result_df[TRANS_FDR] < translocation_p_value]
    if trans_df.shape[0] == 0:  # create only if there are translocations
        html_d[TRANSLOCATIONS][TRANS_HEATMAP_TAB] = ""
        return

    trans_copy_df = deepcopy(trans_df)
    trans_copy_df[SITE_A] = trans_df[SITE_B]
    trans_copy_df[SITE_B] = trans_df[SITE_A]
    trans_df = pd.concat([trans_df, trans_copy_df])
    heat_df = trans_df.pivot(SITE_A, SITE_B, TX_TRANS_READ)
    heat_df[heat_df.isna()] = 0
    heat_df = heat_df.astype(int)
    max_trans_n = heat_df.max().max()
    active_site_n = heat_df.shape[0]
    heat_df['Total'] = np.zeros(active_site_n, dtype=np.int)
    heat_df['NHEJ Activity (%)'] = np.zeros(active_site_n, dtype=np.int)
    total_trans = heat_df.sum()

    # Prepare Heatmap annotations
    annot_df = deepcopy(heat_df)
    annot_df['Total'] = total_trans
    annot_df = annot_df.astype(str)
    editing = []
    for site in heat_df.index:
        editing.append("{:.1f}".format(result_df.loc[result_df[SITE_NAME] == site, EDIT_PERCENT].values[0]))
    annot_df['NHEJ Activity (%)'] = editing

    dpi = 300
    grid_kws = {"height_ratios": (.9, .05), "hspace": .1}

    # create diagonal mask
    mask = np.ones_like(heat_df)
    mask[np.triu_indices_from(mask)] = False
    mask[np.diag_indices(active_site_n)] = True

    # Create figure
    fig_w = 3 + 5*(active_site_n/20)
    fig_h = 2.5 + 5*(active_site_n/20)
    fig, (ax, cbar_ax) = plt.subplots(2, gridspec_kw=grid_kws, figsize=(fig_w, fig_h))

    # Create heat map
    heat_df.fillna(value=0, inplace=True)
    annot_df.fillna(value='0', inplace=True)
    sns.heatmap(heat_df, annot=annot_df.values, linewidths=0.05, linecolor='xkcd:light grey', fmt='s',
                mask=mask, cmap="Reds", vmax=max_trans_n, annot_kws={"size": 9},
                cbar_ax=cbar_ax, cbar_kws={"orientation": "horizontal"}, ax=ax)

    # make colors of the last two columns white
    quadmesh = ax.findobj(QuadMesh)[0]
    facecolors = quadmesh.get_facecolors()
    last_col_index = list(range(active_site_n, active_site_n * (active_site_n + 2), active_site_n + 2))
    last_col_index += list(range(active_site_n + 1, active_site_n * (active_site_n + 2), active_site_n + 2))
    facecolors[last_col_index] = np.array([1, 1, 1, 1])
    quadmesh.set_facecolors = facecolors

    # set marked area color
    # ax.set_facecolor('xkcd:light grey')

    # move x ticks and label to the top
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')

    # remove labels
    ax.set_ylabel("")
    ax.set_xlabel("")

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=-45, ha="right", rotation_mode="anchor", size=11)
    plt.setp(ax.get_yticklabels(), rotation=0, ha="right", rotation_mode="anchor", size=11)
    title = "Translocation Reads Heatmap for FDR < {:.2e}".format(translocation_p_value)

    # Add to final report
    html_d[TRANSLOCATIONS][TRANS_HEATMAP_TAB] = dict()
    html_d[TRANSLOCATIONS][TRANS_HEATMAP_TAB][PLOT_PATH] = os.path.join(OUTPUT_DIR, 'translocations_heatmap.png')
    html_d[TRANSLOCATIONS][TRANS_HEATMAP_TAB][PDF_PATH] = os.path.join(OUTPUT_DIR, 'translocations_heatmap.svg')
    html_d[TRANSLOCATIONS][TRANS_HEATMAP_TAB][TITLE] = title
    html_d[TRANSLOCATIONS][TRANS_HEATMAP_TAB][W] = "{}%".format(min(100, int(40 + 55 * (active_site_n/20))))
    html_d[TRANSLOCATIONS][TRANS_HEATMAP_TAB][H] = dpi*fig_h

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        fig.savefig(os.path.join(output, "translocations_heatmap.png"), bbox_inches='tight', dpi=dpi)
        fig.savefig(os.path.join(output, "translocations_heatmap.svg"), pad_inches=1, box_inches='tight')
        plt.close(fig)


def summary_result_to_excel(summary_result_df: pd.DataFrame, confidence_interval: float, output: str):
    """
    Dump results summary to excel
    :param summary_result_df
    :param confidence_interval
    :param output
    :return:
    """
    df = summary_result_df.copy()
    df.loc[df[ON_TARGET] == True, ON_TARGET] = "on-target"
    df.loc[df[ON_TARGET] == False, ON_TARGET] = "off-target"
    df.rename(columns={CI_LOW: 'CI low ({}%)'.format(100*confidence_interval),
                       CI_HIGH: 'CI high ({}%)'.format(100*confidence_interval)}, inplace=True)

    df.to_csv(os.path.join(output, "results_summary.csv"), index=False)


def warnings_and_discarded_sites_text(summary_result_df: pd.DataFrame, min_num_of_reads: int, donor: bool, html_d: Dict):
    # Print information on discarded reads
    if donor:
        summary_result_df = summary_result_df.loc[~summary_result_df[ON_TARGET]]

    discarded_df = summary_result_df.loc[summary_result_df[EDIT_PERCENT].isna()]

    html_d[READING_STATS][DISCARDED_SITES] = ""
    # Display discarded sites
    if discarded_df.shape[0] > 0:
        html_d[READING_STATS][DISCARDED_SITES] += "<p style=\"text-align: left\"> {} sites were discarded due to" \
                                                  " low number of reads (below {:,}):<p>".format(discarded_df.shape[0],
                                                                                                 min_num_of_reads)

        for row_idx, row in discarded_df.iterrows():
            html_d[READING_STATS][DISCARDED_SITES] += "<p style=\"text-align: left\"> {} - Treatment reads - {:,}. " \
                                                      "Mock reads - {:,}.<p>".format(row[SITE_NAME], row[TX_READ_NUM],
                                                                                     row[MOCK_READ_NUM])
    # Display all warnings
    logger = LoggerWrapper.get_logger()
    for warn_msg in logger.warning_msg_l:
        for idx in range(0, len(warn_msg), 100):
            current_msg = warn_msg[idx:min(idx + 100, len(warn_msg))]
            html_d[READING_STATS][DISCARDED_SITES] += "<p style=\"text-align: left\">" + current_msg + "<p>"

    # No warnings msg:
    if html_d[READING_STATS][DISCARDED_SITES] == "":
        html_d[READING_STATS][DISCARDED_SITES] = "No warnings"


# Edit read table utils
def get_read_around_cut_site(read, cut_site, length):
    return read[cut_site-length:cut_site+length]


def color_edit_background(read_w_ins: List[str], read_w_del: List[str]):
    colors = []
    for base_w_ins, base_w_del in zip(read_w_ins, read_w_del):
        if base_w_ins == "-":
            colors.append('indianred')
        elif base_w_del == "-":
            colors.append('lightgrey')
        elif base_w_ins != base_w_del:
            colors.append('royalblue')
        elif base_w_del == "A":
            colors.append('#8cb6d6') # light blue
        elif base_w_del == "T":
            colors.append('#ffba85') # light orange
        elif base_w_del == "G":
            colors.append('#cae6ca') # light green
        elif base_w_del == "C":
            colors.append('#c5aeda') # light purple

    return colors