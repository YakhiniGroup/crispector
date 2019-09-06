import logging
import os
import yaml # TODO - add to project requirements - This is conda install pyyaml and not "yaml"!!!
from exceptions import ConfiguratorIsCalledBeforeInitConfigPath
from constants_and_types import Path, AlgResult, CI_HIGH, EDIT_PERCENT, SITE_NAME, CI_LOW, ON_TARGET
from matplotlib import pyplot as plt
import matplotlib as mpl
import math

class Logger:
    """
    Singleton logger based on logging package.
    Dump all messages both to shell and crispector_main.log.
    """
    _configured = False
    _OUTPUT_DIR = None
    _logger_level = logging.DEBUG

    @classmethod
    def get_logger(cls):
        # Create a custom logger
        logger = logging.getLogger("default")
        logger.level = cls._logger_level

        if not cls._configured:
            cls._configured = True

            # Create handlers
            c_handler = logging.StreamHandler()
            logger_path = os.path.join(cls._OUTPUT_DIR, 'crispector_main.log')
            if os.path.exists(logger_path):
                os.remove(logger_path)
            f_handler = logging.FileHandler(logger_path)

            f_handler.setLevel(cls._logger_level)
            c_handler.setLevel(cls._logger_level)

            # Create formatters and add it to handlers
            f_format = logging.Formatter('%(asctime)s %(levelname)s\t %(message)s')
            c_format = logging.Formatter('%(asctime)s %(levelname)s\t %(message)s')
            f_handler.setFormatter(f_format)
            c_handler.setFormatter(c_format)

            # Add handlers to the logger
            logger.addHandler(f_handler)
            logger.addHandler(c_handler)
        return logger

    @classmethod
    def set_log_path(cls, path: Path):
        cls._OUTPUT_DIR = path

    @classmethod
    def set_logger_level(cls, mode):
        cls._logger_level = mode


class Configurator:
    """
    Singleton YAML configurator based on yaml package.
    """

    _config_file = None
    _CONFIG_PATH = None

    @classmethod
    def set_cfg_path(cls, path: Path):
        if path is None:
            path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'config/default_config.yml')
        cls._CONFIG_PATH = path

    @classmethod
    def get_cfg(cls):
        if cls._config_file is None:
            if cls._CONFIG_PATH is None:
                raise ConfiguratorIsCalledBeforeInitConfigPath()

            # Read YAML file
            with open(cls._CONFIG_PATH, 'r') as stream:
                cls._config_file = yaml.safe_load(stream)

        return cls._config_file


def plot_editing_activity(result_df: AlgResult, confidence_interval: float, editing_threshold: float, output: Path):
    # Set font
    mpl.rcParams.update(mpl.rcParamsDefault)
    mpl.rcParams['font.size'] = 20
    mpl.rcParams['xtick.labelsize'] = 20
    mpl.rcParams['ytick.labelsize'] = 20
    mpl.rcParams['axes.labelsize'] = 24
    mpl.rcParams['axes.titlesize'] = 26
    mpl.rcParams['legend.fontsize'] = 24
    editing_bar_text_size = 18

    off_target_color = '#d0743c'
    on_target_color = '#7b6888'
    # Filter all low editing activity sites
    result_df = result_df.dropna()
    edit_df = result_df.loc[result_df[CI_HIGH] >= editing_threshold]

    # Sort experiments
    edit_df = edit_df.sort_values(by=EDIT_PERCENT, ascending=False)

    # Create axes
    max_bars = 20
    bar_num = edit_df.shape[0]
    plot_num = math.ceil(bar_num / max_bars)
    # set dynamic bar_width - according to the number of bars
    bar_width = 0.9 if plot_num > 1 else 0.9 * (0.5 + 0.5 * bar_num / max_bars)
    fig, axes = plt.subplots(nrows=plot_num, ncols=1, figsize=(20, plot_num * 6), constrained_layout=True)

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
        bar_plot = axes[idx].bar(bar_pos, plt_editing, width=bar_width, color=off_target_color,
                                 yerr=[plt_CI_low, plt_CI_high], align='center', ecolor='black',
                                 capsize=dynamic_capsize)
        for site_idx, is_site_on_target in enumerate(plt_on_target):
            if is_site_on_target:
                bar_plot[site_idx].set_color(on_target_color)

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
            axes[idx].text(x=bar_pos[text_idx] - 0.45 + bar_width / 2., y=text_height[text_idx],
                           s="{:.2f}".format(plt_editing[text_idx]), ha='center', va='bottom',
                           size=editing_bar_text_size)

        if idx == 0:
            axes[idx].set_title(r"Editing Activity with {} % CI above {}%".format(confidence_interval,
                                                                                  editing_threshold), weight='bold')
            # Add legend
            axes[idx].bar([0], [y_lim], color=on_target_color, label="On-Target")
            axes[idx].bar([0], [y_lim], color=off_target_color, label="Off-Target")
            axes[idx].legend(loc='upper right')

    fig.savefig(os.path.join(output, 'editing_activity.png'), box_inches='tight', dpi=300)
    plt.close(fig)

