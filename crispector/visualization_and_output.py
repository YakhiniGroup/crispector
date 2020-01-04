from constants_and_types import SITE_NAME, IndelType, Path, FREQ, IS_EDIT, TX_READ_NUM, \
    MOCK_READ_NUM, TX_EDIT, EDIT_PERCENT, CI_LOW, CI_HIGH, READ_LEN_SIDE, ALIGN_CUT_SITE, ALIGNMENT_W_INS, \
    ALIGNMENT_W_DEL, POS_IDX_E, POS_IDX_S, INDEL_TYPE, ExpType, ON_TARGET, IsEdit, C_TX, C_MOCK, \
    SUMMARY_RESULTS_TITLES, OFF_TARGET_COLOR, ON_TARGET_COLOR, OUTPUT_DIR, DISCARDED_SITES, \
    EDITING_ACTIVITY, PLOT_PATH, TITLE, W, H, PDF_PATH, PAGE_TITLE, READING_STATS, MAPPING_STATS, MAPPING_PER_SITE, \
    FASTP_DIR, FASTP_TX_PATH, FASTP_MOCK_PATH, RESULT_TABLE, TAB_DATA, HTML_SITE_NAMES, LOG_PATH, TransDf, AlgResultDf, \
    TransResultDf
import math
import os
import warnings
from typing import List, Tuple, Dict
from input_processing import InputProcessing
from modification_types import ModificationTypes
from core_algorithm import CoreAlgorithm
from modification_tables import ModificationTables
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import seaborn as sns #
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
from utils import Logger


def create_site_output(algorithm: CoreAlgorithm, modifications: ModificationTypes, mod_table: ModificationTables,
                       site_result: Dict, site_name: str, experiment_name: str, output: Path):
    base_path = os.path.join(OUTPUT_DIR, site_name)
    # TODO - Don't forget to add PDF's and remove titles.

    title_name = "{} - {}".format(experiment_name, site_name)

    cut_site = algorithm.cut_site
    win_size = algorithm.win_size

    # plot modification table
    # TODO - decide which plots are needed
    # TODO - needed? tables.dump_tables(self._edit, self._tables_offset, self._output)
    all_table_idx = list(range(modifications.size))
    del_table_idx = [i for i, x in enumerate(modifications.types) if x == IndelType.DEL]
    ins_table_idx = [i for i, x in enumerate(modifications.types) if x == IndelType.INS]
    mix_table_idx = [i for i, x in enumerate(modifications.types) if x in [IndelType.MIXED, IndelType.SUB]]
    plot_modification_tables(mod_table, modifications, algorithm.edited, algorithm.tables_offset, title_name,
                             all_table_idx, output, "all", figsize=(14, 20))
    plot_modification_tables(mod_table, modifications, algorithm.edited, algorithm.tables_offset, title_name,
                             del_table_idx, output, "del")
    plot_modification_tables(mod_table, modifications, algorithm.edited, algorithm.tables_offset, title_name,
                             ins_table_idx, output, "ins")
    plot_modification_tables(mod_table, modifications, algorithm.edited, algorithm.tables_offset, title_name,
                             mix_table_idx, output, "mix_and_sub")

    # Create edited read table
    plot_edited_reads_to_table(mod_table, cut_site, output)

    # Plot mutation distribution
    plot_distribution_of_edit_events(mod_table, cut_site, win_size, title_name, output)
    plot_distribution_of_all_modifications(mod_table, cut_site, win_size, title_name, output)
    # plot_distribution_of_edit_event_sizes(mod_table, title_name, output) TODO - which one is better
    plot_distribution_of_edit_event_sizes_3_plots(mod_table, title_name, output)

    # Plot editing activity
    plot_site_editing_activity(algorithm, site_result, site_name, title_name, output)

    # Dump .csv file with all reads
    mod_table.tx_reads.to_csv(os.path.join(output, "treatment_aligned_reads.csv.gz"), index=False, compression='gzip')
    mod_table.mock_reads.to_csv(os.path.join(output, "mock_aligned_reads.csv.gz"), index=False, compression='gzip')


def create_experiment_output(result_df: AlgResultDf, tx_trans_df: TransDf, mock_trans_df: TransDf,
                             trans_result_df: TransResultDf, input_processing: InputProcessing, min_num_of_reads: int,
                             confidence_interval: float, editing_threshold: float, experiment_name: str, output: Path):

    html_d = dict() # html parameters dict
    html_d[PAGE_TITLE] = experiment_name
    html_d[READING_STATS] = dict()

    # Dump summary results
    summary_result_to_excel(result_df, confidence_interval, output)

    # Create bar plot for editing activity
    plot_editing_activity(result_df, confidence_interval, editing_threshold, html_d, output)

    # Dump reads statistics
    tx_input_n, tx_merged_n, tx_aligned_n = input_processing.read_numbers(ExpType.TX)
    mock_input_n, mock_merged_n, mock_aligned_n = input_processing.read_numbers(ExpType.MOCK)
    create_reads_statistics_report(result_df, tx_input_n, tx_merged_n, tx_aligned_n,
                                   mock_input_n, mock_merged_n, mock_aligned_n, html_d, output)

    # Create a text file with all discarded sites
    discarded_sites_text(result_df, min_num_of_reads, html_d, output)

    # TODO - add all translocation outputs....
    tx_trans_df.to_csv(os.path.join(output, "tx_translocations_reads.csv"), index=False)
    mock_trans_df.to_csv(os.path.join(output, "mock_translocations_reads.csv"), index=False)
    trans_result_df.to_csv(os.path.join(output, "translocations_results.csv"), index=False)

    # Add summary results to page
    html_d[RESULT_TABLE] = dict()
    html_d[RESULT_TABLE][TITLE] = "Results Table"
    html_d[RESULT_TABLE][TAB_DATA] = dict()
    for col in SUMMARY_RESULTS_TITLES:
        html_d[RESULT_TABLE][TAB_DATA][col] = list(result_df[col].values)

    html_d[HTML_SITE_NAMES] = list(result_df.loc[result_df[EDIT_PERCENT].isna(), SITE_NAME].values)

    # Add fastp links
    if os.path.exists(os.path.join(output, FASTP_DIR[ExpType.TX])):
        html_d[READING_STATS][FASTP_TX_PATH] = os.path.join(OUTPUT_DIR, "{}/fastp.html".format(FASTP_DIR[ExpType.TX]))
    if os.path.exists(os.path.join(output, FASTP_DIR[ExpType.MOCK])):
        html_d[READING_STATS][FASTP_MOCK_PATH] = os.path.join(OUTPUT_DIR, "{}/fastp.html".format(FASTP_DIR[ExpType.MOCK]))

    html_d[LOG_PATH] = os.path.join(OUTPUT_DIR, Logger.logger_name)

    return html_d

#####----------------------#####
#####----modifications-----#####
#####----------------------#####
def plot_modification_tables(mod_table: ModificationTables, modifications: ModificationTypes, edit_table: IsEdit,
                             table_offset: int, experiment_name: str, table_indexes: List[int], output: Path,
                             name_suffix = "", figsize: Tuple = (16, 11)):
    """
    Plot all modification tables around cut-site.
    Also display edit events.
    :param mod_table: ModificationTables
    :param modifications: ModificationTypes
    :param edit_table: edit as marked by crispector algorithm
    :param table_offset: reference offset (in bp) to the start of the qualification window
    :param output: output path
    :param experiment_name
    :param table_indexes - list indexes for table plot
    :param name_suffix - prefix for file name
    :param figsize
    :return:
    """
    # Set font
    mpl.rcParams.update(mpl.rcParamsDefault)
    mpl.rcParams['font.size'] = 22
    mpl.rcParams['ytick.labelsize'] = 12
    indel_size = 12
    bar_width = 0.4
    mock_color = '#3690c0'  # blue
    tx_color = '#f16a13'  # orange
    edit_color = '#dcdcdc'  # light grey
    cut_site_color = "red"
    grid_color = 'grey'
    # Create axes
    fig, axes = plt.subplots(nrows=len(table_indexes), ncols=1, figsize=figsize, sharex=True)

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
        axes[axes_idx].vlines(grid_positions, ymin=0, ymax=y_max, color=grid_color, lw=1)
        axes[axes_idx].spines['bottom'].set_color(grid_color)
        axes[axes_idx].spines['top'].set_color(grid_color)
        axes[axes_idx].spines['right'].set_color(grid_color)
        axes[axes_idx].spines['left'].set_color(grid_color)

        # Set cut-site in red (Insertions cut site is a full bp cell)
        if is_ins:
            axes[axes_idx].vlines([cut_site - 1, cut_site], ymin=0, ymax=y_max, color=cut_site_color, lw=2)
            axes[axes_idx].hlines([0, y_max], xmin=cut_site - 1, xmax=cut_site, color=cut_site_color, lw=4)
        else:
            axes[axes_idx].axvline(cut_site, ymin=0, ymax=y_max, color=cut_site_color, lw=2)

        # Create bar plot
        axes[axes_idx].bar(bar_ind - bar_width / 2, tx, width=bar_width, color=tx_color,
                           label="Tx ({:,} Reads)".format(mod_table.n_reads_tx))
        axes[axes_idx].bar(bar_ind + bar_width / 2, mock, width=bar_width, color=mock_color,
                           label="Mock ({:,} Reads)".format(mod_table.n_reads_mock))

        # Set x, y lim & ticks
        axes[axes_idx].set_xlim(min(positions) - 0.5, max(positions) + 1.5)
        axes[axes_idx].set_xticks([])
        axes[axes_idx].set_yticks([])
        axes[axes_idx].set_ylim(0, y_max)

        # Set y_label - Use text due to alignment issues
        axes[axes_idx].text(x=positions[0] - 5 - 0.5 * (not is_ins), y=y_max / 2,
                             s=modifications.plot_name_at_idx(table_idx), ha="left", va="center")

        # Create legend in the middle of the plot
        if axes_idx == ((len(table_indexes) // 2) - 1):
            axes[axes_idx].bar([0], [0], color=edit_color, label="Edit event", edgecolor='grey')
            axes[axes_idx].plot([], [], color=cut_site_color, label="Cut-Site")
            handles, labels = axes[axes_idx].get_legend_handles_labels()
            order = [1, 2, 3, 0]
            axes[axes_idx].legend([handles[idx] for idx in order], [labels[idx] for idx in order],
                                   bbox_to_anchor=(1.48, 0))

        # Add the reference sequence in the position of the title
        if axes_idx == 0:
            axes[axes_idx].text(0.5, 0.95,
                                "CRISPECTOR classifier results, by position\n{}".format(experiment_name),
                                ha="center", va="bottom", weight='bold', size=35,
                                transform=fig.transFigure, family='serif')
            ref = mod_table.amplicon[table_offset:table_offset + len(positions)]
            for pos_idx, ii in enumerate(bar_ind):
                axes[axes_idx].text(ii, y_max, ref[pos_idx], ha="center", va="bottom", weight='bold', size=30)
                axes[axes_idx].text(ii, y_max, ref[pos_idx], ha="center", va="bottom", weight='bold', size=30)
            # red cute-site
            axes[axes_idx].text(cut_site, 1.1 * y_max, "|", ha="center", va="bottom", weight='bold', size=20,
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
                    bbox_inches='tight', dpi=200)
        plt.close(fig)

def dump_tables(mod_table: ModificationTables, modifications: ModificationTypes,
                edit_table: IsEdit, table_offset: int, win_size: int, output: Path):
    """
    Dump all modification tables around cut-site to a csv.
    Also display edit events.
    :param mod_table: ModificationTables
    :param modifications: ModificationTypes
    :param edit_table: edit as marked by crispector algorithm
    :param table_offset: reference offset (in bp) to the start of the qualification window
    :param output: output path
    :return:
    """
    ref = mod_table.amplicon[table_offset:table_offset + 2*win_size]
    tables_dict = dict()
    new_cols = []
    tables_df = pd.DataFrame()

    for table_idx in range(modifications.size):
        # Gather data
        is_ins = modifications.types[table_idx] == IndelType.INS
        edit = edit_table[table_idx]
        pos_len = len(edit)  # length of positions
        tx = mod_table.tables[table_idx][C_TX, table_offset:table_offset + pos_len]
        mock = mod_table.tables[table_idx][C_MOCK, table_offset:table_offset + pos_len]
        cut_site = np.zeros(pos_len, dtype=np.int)
        cut_site[(pos_len // 2) - 1 * (not is_ins):(pos_len // 2) + 1] = 1

        # Insert to dictionary
        table_name = modifications.name_at_idx(table_idx)
        tables_dict[table_name + " Tx"] = np.insert(tx, 0, mod_table.n_reads_tx)
        tables_dict[table_name + " Mock"] = np.insert(mock, 0, mod_table.n_reads_mock)
        tables_dict[table_name + " Is edited"] = [None] + edit
        tables_dict[table_name + " Is cut site"] = [None] + list(cut_site)

        tables_df = pd.DataFrame.from_dict(tables_dict, orient='index')
        new_cols = ['number of reads'] + [c for c in ref] + ['insertion at the end']

    tables_df.set_axis(new_cols, axis=1, inplace=True)
    tables_df.to_csv(os.path.join(output, "modification_tables.csv"))

#####----------------------#####
#####---------Site---------#####
#####----------------------#####
def plot_distribution_of_all_modifications(tables: ModificationTables, cut_site: int, win_size: int, experiment_name: str,
                                           output: Path):
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
    mpl.rcParams['font.size'] = 18
    mpl.rcParams['axes.labelsize'] = 22
    mpl.rcParams['axes.titlesize'] = 24

    # Create axes
    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(16, 9), sharex=True)
    plt.subplots_adjust(wspace=0, hspace=0)

    # Plot distributions
    for exp_type in ExpType:
        dist_d = tx_dist_d if exp_type == ExpType.TX else mock_dist_d
        ax = axes[0] if exp_type == ExpType.TX else axes[1]
        positions = list(range(amplicon_length + 1))
        ax.axvline(x=cut_site, linestyle='-', color="red", label='Expected cut-site', linewidth=2)
        ax.axvline(x=cut_site - win_size, linestyle='--', color='k', label='Quantification window', linewidth=1)
        ax.axvline(x=cut_site + win_size, linestyle='--', color='k', linewidth=1)
        ax.plot(positions, dist_d[IndelType.INS], color=IndelType.INS.color, label="Insertions", linewidth=3, alpha=0.9)
        ax.plot(positions, dist_d[IndelType.SUB], color=IndelType.SUB.color, label="Substitutions", linewidth=3, alpha=0.9)
        ax.plot(positions, dist_d[IndelType.DEL], color=IndelType.DEL.color, label="Deletions", linewidth=3, alpha=0.9)
        ax.set_ylim(bottom=0, top=max(int(1.1 * max_indels), 10))
        ax.set_xlim(left=0, right=amplicon_length)
        ax.set_ylabel("{}\nIndel count".format(exp_type.name()))

    # Create legend, title and x-label
    axes[0].legend()
    axes[0].set_title("Distribution Of All Modifications\n{}".format(experiment_name),
                      weight='bold', family='serif')

    # Remove Tx first y axe label (overlap with mock
    plt.setp(axes[0].get_yticklabels()[0], visible=False)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        fig.savefig(os.path.join(output, 'distribution_of_all_modifications.png'),
                    bbox_inches='tight', dpi=100)
        plt.close(fig)


def plot_distribution_of_edit_events(tables: ModificationTables, cut_site: int, win_size: int, experiment_name: str,
                                     output: Path):
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
    mpl.rcParams['font.size'] = 22
    mpl.rcParams['axes.labelsize'] = 30
    mpl.rcParams['axes.titlesize'] = 26

    # Create axes
    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_axes([0, 0, 1, 1])

    positions = list(range(amplicon_length + 1))
    ax.axvline(x=cut_site, linestyle='-', color="red", label='Expected cut-site', linewidth=2)
    ax.axvline(x=cut_site - win_size, linestyle='--', color='k',
               label='Quantification window', linewidth=1)
    ax.axvline(x=cut_site + win_size, linestyle='--', color='k', linewidth=1)

    for indel_type in [IndelType.INS, IndelType.SUB, IndelType.DEL]:
        ax.plot(positions, dist_d[indel_type], color=indel_type.color, label=indel_type.name, linewidth=3, alpha=0.9)

    # Create legend, axes limit and labels
    ax.legend()
    ax.set_title("Distribution Of Edited Events\n{}".format(experiment_name),
                 weight='bold', family='serif')
    max_indels = np.array([dist_d[IndelType.INS], dist_d[IndelType.DEL], dist_d[IndelType.SUB]]).max()
    ax.set_ylim(bottom=0, top=max(int(1.1 * max_indels), 10))
    ax.set_xlim(left=0, right=amplicon_length)
    ax.set_xlabel("Position")
    ax.set_ylabel("Indel count")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        fig.savefig(os.path.join(output, 'distribution_of_edit_events.png'),
                    bbox_inches='tight', dpi=100)
        plt.close(fig)


def plot_edited_reads_to_table(tables: ModificationTables, cut_site: int, output: Path):
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


def plot_site_editing_activity(algorithm: CoreAlgorithm, result_d: Dict, site_name: str, experiment_name: str,
                               output: Path):
    # Set font
    mpl.rcParams.update(mpl.rcParamsDefault)
    sns.set(style="whitegrid")
    mpl.rcParams['font.size'] = 14
    mpl.rcParams['xtick.labelsize'] = 14
    mpl.rcParams['ytick.labelsize'] = 14
    mpl.rcParams['axes.labelsize'] = 16
    mpl.rcParams['axes.titlesize'] = 16

    bar_width = 0.8

    if algorithm.is_on_target:
        bar_color = ON_TARGET_COLOR
    else:
        bar_color = OFF_TARGET_COLOR

    # Define fix and axes
    fig = plt.figure(figsize=(4, 4))
    ax = fig.add_axes([0, 0, 1, 1])

    # Get bar data
    editing = result_d[EDIT_PERCENT]
    CI_high = result_d[CI_HIGH] - editing
    CI_low = editing - result_d[CI_LOW]

    # plot bar
    ax.bar([0], [result_d[EDIT_PERCENT]], color=bar_color, width=bar_width,
           yerr=[[CI_low], [CI_high]], align='center', ecolor='black', capsize=15)

    # Set labels
    ax.set_ylabel("Editing Activity (%)")

    # Set scale and lim
    y_lim = max(min(1.2 * (editing + CI_high), 100), 0.1)
    ax.set_ylim(0, y_lim)
    ax.set_xlim(-0.8, 0.8)

    # Text below each bar plot + y ticks
    ax.set_xticks([0])
    ax.set_xticklabels([site_name])

    ax.set_title("Editing Activity with {} % CI\n{}".format(100 * algorithm.confidence, experiment_name),
                 weight='bold', family='serif')

    ax.text(x=-0.1, y=-0.3, s="Number of edited reads\nEditing activity",
            ha='left', va='bottom', transform=fig.transFigure, family='serif')
    ax.text(x=0.55, y=-0.3, s="- {:,} (out of {:,} reads).\n"
                              "- {:.2f}%, CI=({:.2f}%$-${:.2f}%).".format(result_d[TX_EDIT], result_d[TX_READ_NUM],
                                                                          editing, result_d[CI_LOW],
                                                                          result_d[CI_HIGH]),
            ha='left', va='bottom', transform=fig.transFigure, family='serif')

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        fig.savefig(os.path.join(output, 'site_editing_activity.png'),
                    bbox_inches='tight', dpi=200)
        plt.close(fig)


def plot_distribution_of_edit_event_sizes(tables: ModificationTables, experiment_name: str, output: Path):
    # Set font
    mpl.rcParams.update(mpl.rcParamsDefault)
    mpl.rcParams['font.size'] = 20
    mpl.rcParams['xtick.labelsize'] = 20
    mpl.rcParams['ytick.labelsize'] = 20
    mpl.rcParams['axes.labelsize'] = 24
    mpl.rcParams['axes.titlesize'] = 26
    mpl.rcParams['legend.fontsize'] = 24
    bar_width = 1

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
    r = {}
    for idx, indel in enumerate([IndelType.DEL, IndelType.INS, IndelType.SUB]):
        r[indel] = list(range(0 + idx, 4 * bin_num + idx, 4))

    # Create figure
    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_axes([0, 0, 1, 1])

    # Create bar plot
    for indel in [IndelType.DEL, IndelType.INS, IndelType.SUB]:
        ax.bar(r[indel], values[indel], label=indel.name, color=indel.color, edgecolor='black', alpha=0.9,
               width=bar_width)

    # Labels, ticks and title
    plt.xlabel('Edit event length')
    plt.ylabel('Number of modifications')
    plt.legend(loc='upper right')
    ax.set_xticks(r[IndelType.INS])
    ax.set_xticklabels(bins_names, rotation='vertical')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_title("Edit Event Size Distribution\n{}".format(experiment_name), weight='bold', family='serif')

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        fig.savefig(os.path.join(output, 'distribution_of_edit_events_size.png'),
                    bbox_inches='tight', dpi=100)
        plt.close(fig)


def plot_distribution_of_edit_event_sizes_3_plots(tables: ModificationTables, experiment_name: str, output: Path):
    # Set font
    mpl.rcParams.update(mpl.rcParamsDefault)
    mpl.rcParams['font.size'] = 16
    mpl.rcParams['xtick.labelsize'] = 16
    mpl.rcParams['ytick.labelsize'] = 16
    mpl.rcParams['axes.labelsize'] = 18
    mpl.rcParams['axes.titlesize'] = 26
    mpl.rcParams['legend.fontsize'] = 18
    bar_width = 1

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
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(16, 6), constrained_layout=True)

    # Create bar plot
    for idx, indel in enumerate([IndelType.DEL, IndelType.INS, IndelType.SUB]):
        axes[idx].bar(r, values[indel], label=indel.name, color=indel.color, edgecolor='black', alpha=0.9,
                      width=bar_width)
        axes[idx].set_xticks(r)
        axes[idx].set_xticklabels(bins_names, rotation='vertical')
        axes[idx].spines['top'].set_visible(False)
        axes[idx].spines['right'].set_visible(False)

        # Labels, ticks and title
        axes[idx].set_xlabel('Edit event length')
        axes[idx].set_ylabel('Number of modifications')
        axes[idx].legend(loc='upper right')

    fig.suptitle("Edit Event Size Distribution\n{}".format(experiment_name), weight='bold', family='serif')

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        fig.savefig(os.path.join(output, 'distribution_of_edit_events_size_3_plots.png'),
                    bbox_inches='tight', dpi=100)
        plt.close(fig)

#####----------------------#####
#####------Experiment------#####
#####----------------------#####
def plot_editing_activity(result_df: AlgResultDf, confidence_interval: float, editing_threshold: float, html_d: Dict,
                          output: Path):

    # Set font
    mpl.rcParams.update(mpl.rcParamsDefault)
    mpl.rcParams['font.size'] = 20
    mpl.rcParams['xtick.labelsize'] = 20
    mpl.rcParams['ytick.labelsize'] = 20
    mpl.rcParams['axes.labelsize'] = 24
    mpl.rcParams['axes.titlesize'] = 26
    mpl.rcParams['legend.fontsize'] = 24
    editing_bar_text_size = 18
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
    fig_w, fig_h = 20, plot_num * 6
    bar_width = 0.9 if plot_num > 1 else 0.9 * (0.5 + 0.5 * bar_num / max_bars)
    fig, axes = plt.subplots(nrows=plot_num, ncols=1, figsize=(fig_w, fig_h), constrained_layout=True)

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
        dynamic_capsize = 10 + 10 * (1 - number_of_bars / max_bars)
        bar_plot = axes[idx].bar(bar_pos, plt_editing, width=bar_width, color=OFF_TARGET_COLOR,
                                 yerr=[plt_CI_low, plt_CI_high], align='center', ecolor='black',
                                 capsize=dynamic_capsize)
        for site_idx, is_site_on_target in enumerate(plt_on_target):
            if is_site_on_target:
                bar_plot[site_idx].set_color(ON_TARGET_COLOR)

        # Add horizontal line
        axes[idx].axhline(y=editing_threshold, linewidth=1, color='k', linestyle="--", alpha=0.75)

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
            axes[idx].set_xlim(0, number_of_bars + 1)

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

    fig.savefig(os.path.join(output, 'editing_activity.png'), box_inches='tight', dpi=dpi)
    fig.savefig(os.path.join(output, 'editing_activity.pdf'), pad_inches = 1, box_inches='tight')

    html_d[EDITING_ACTIVITY] = dict()
    html_d[EDITING_ACTIVITY][PLOT_PATH] = os.path.join(OUTPUT_DIR, 'editing_activity.png')
    html_d[EDITING_ACTIVITY][PDF_PATH] = os.path.join(OUTPUT_DIR, 'editing_activity.pdf')
    html_d[EDITING_ACTIVITY][TITLE] = title
    html_d[EDITING_ACTIVITY][W] = dpi*fig_w
    html_d[EDITING_ACTIVITY][H] = dpi*fig_h

    plt.close(fig)


def create_reads_statistics_report(result_df: AlgResultDf, tx_in: int, tx_merged: int, tx_aligned: int,
                                   mock_in: int, mock_merged: int, mock_aligned: int, html_d: Dict, output: Path):
    html_d[READING_STATS] = dict()
    html_d[READING_STATS][TITLE] = "Reading Statistics"
    # Set font
    mpl.rcParams.update(mpl.rcParamsDefault)
    sns.set(style="whitegrid")
    mpl.rcParams['font.size'] = 18
    mpl.rcParams['ytick.labelsize'] = 16
    mpl.rcParams['xtick.labelsize'] = 18
    mpl.rcParams['axes.labelsize'] = 20
    mpl.rcParams['axes.titlesize'] = 22
    dpi = 200
    fig_w, fig_h = 14, 8
    bar_width = 0.4
    mock_color = '#3690c0'  # blue
    mock_color_lighter = '#72b1d2'
    mock_color_lightest = '#9ac7df'
    tx_color = '#f16a13'  # orange
    tx_color_lighter = '#f59659'  # orange
    tx_color_lightest = '#f7b488'  # orange

    # set input == merged if input information isn't available
    if tx_in == -1:
        tx_in = tx_merged
    if mock_in == -1:
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
    html_d[READING_STATS][MAPPING_STATS][PDF_PATH] = os.path.join(OUTPUT_DIR, 'mapping_statistics.pdf')
    html_d[READING_STATS][MAPPING_STATS][TITLE] = title
    html_d[READING_STATS][MAPPING_STATS][W] = dpi*fig_w
    html_d[READING_STATS][MAPPING_STATS][H] = dpi*fig_h

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        fig.savefig(os.path.join(output, 'mapping_statistics.png'), bbox_inches='tight', dpi=dpi)
        fig.savefig(os.path.join(output, 'mapping_statistics.pdf'), bbox_inches='tight', pad_inches=1)
        plt.close(fig)

    # Create reads box_plot
    fig = plt.figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes([0, 0, 1, 1])

    bplot = sns.boxplot(x=["Treatment", "Mock"],
                        y=[result_df[TX_READ_NUM], result_df[MOCK_READ_NUM]],
                        linewidth=2.5, ax=ax)
    txbox = bplot.artists[0]
    txbox.set_facecolor(tx_color)
    mockbox = bplot.artists[1]
    mockbox.set_facecolor(mock_color)

    # Set x, y lim & ticks and title
    ax.set_xlim(-1, 2)
    ax.set_ylabel("Number Of Reads")
    title = "Number Of Aligned Reads Per Site"
    # ax.set_title(,title weight='bold', family='serif')

    html_d[READING_STATS][MAPPING_PER_SITE] = dict()
    html_d[READING_STATS][MAPPING_PER_SITE][PLOT_PATH] = os.path.join(OUTPUT_DIR, 'number_of_aligned_reads_per_site.png')
    html_d[READING_STATS][MAPPING_PER_SITE][PDF_PATH] = os.path.join(OUTPUT_DIR, 'number_of_aligned_reads_per_site.pdf')
    html_d[READING_STATS][MAPPING_PER_SITE][TITLE] = title
    html_d[READING_STATS][MAPPING_PER_SITE][W] = dpi*fig_w
    html_d[READING_STATS][MAPPING_PER_SITE][H] = dpi*fig_h

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        fig.savefig(os.path.join(output, 'number_of_aligned_reads_per_site.png'), bbox_inches='tight', dpi=dpi)
        fig.savefig(os.path.join(output, 'number_of_aligned_reads_per_site.pdf'), bbox_inches='tight', pad_inches=1)
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

    df.to_excel(os.path.join(output, "results_summary.xlsx"), index=False, float_format='%.4f')


def discarded_sites_text(summary_result_df: pd.DataFrame, min_num_of_reads: int, html_d: Dict, output: str):
    # Print information on discarded reads
    discarded_df = summary_result_df.loc[summary_result_df[EDIT_PERCENT].isna()]

    with open(os.path.join(output, "discarded_sites.txt"), 'w') as file:

        if discarded_df.shape[0] > 0:
            opening_line = "{} sites were discarded due to low number of reads (below {:,}):\n".format(
                discarded_df.shape[0], min_num_of_reads)
            file.write(opening_line)
            site_lines = []
            for row_idx, row in discarded_df.iterrows():
                site_lines.append("{} - Treatment reads - {:,}. Mock reads - {:,}.\n".format(row[SITE_NAME],
                                                                                             row[TX_READ_NUM],
                                                                                             row[MOCK_READ_NUM]))
                file.write(site_lines[-1])

            html_d[READING_STATS][DISCARDED_SITES] = opening_line + "".join(site_lines)

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