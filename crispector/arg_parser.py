"""Console script for crispector."""
import sys
import click
from crispector_main import run
import os


@click.command(context_settings=dict(max_content_width=120))
@click.option('--tx_in1', '-t_r1', type=click.Path(exists=True),
              help="Tx read 1 input path or Tx merged FASTQ file")
@click.option('--tx_in2', '-t_r2', type=click.Path(exists=True),
              help="Tx read 2 input path, if FASTQ files aren't merged [OPTIONAL]")
@click.option('--mock_in1', '-m_r1', type=click.Path(exists=True),
              help="Mock read 1 input path or mock merged FASTQ file")
@click.option('--mock_in2', '-m_r2', type=click.Path(exists=True),
              help="Mock read read 2 input path, if FASTQ files aren't merged [OPTIONAL] ")
@click.option("--experiment_config", '-c', type=click.Path(exists=True), required=True,
              help="A CSV (Comma Separated Values‏) file with the experiment data. Table has 11 columns: "
                   "SiteName, AmpliconReference, sgRNA, OnTarget, ForwardPrimer, ReversePrimer,"
                   "TxInput1Path TxInput2Path, MockInput1Path, MockInput2Path, DonorReference. "
                   "The first 4 columns are required, the rest are optional. "
                   "Header should be specified by the above order. "
                   "Please check the README on GitHub further details and examples.")
@click.option('--report_output', '-o', type=click.Path(), default="CRISPECTOR", help="Path to output folder")
@click.option("--cut_site_position", type=click.INT, default=-3, show_default=True,
              help="Expected cut-site position with respect to the 3' end of the provided sgRNA sequence. "
                   "Note, the sgRNA sequence must be entered without the PAM.")
@click.option("--crispector_config", type=click.Path(),
              help="Path to crispector configuration in YAML format."
                   "See \"Advanced usage\" section in README on GitHub for further.")
@click.option('--fastp_options_string', type=click.STRING, default="-w 2",
              help="Try \"fastp --help\" for more details")
@click.option("--min_num_of_reads", type=click.INT, default=500, show_default=True,
              help="Minimum number of reads (per locus site) to evaluate edit events")
@click.option("--min_read_length_without_primers", type=click.INT, default=10, show_default=True,
              help="Filter out any read shorter than min_read_length_without_primers + length of forward and "
                   "reverse primers. This threshold filters primer-dimmer effect reads.")
@click.option("--max_edit_distance_on_primers", type=click.INT, default=8, show_default=True,
              help="Maximum edit distance to consider a read prefix (or suffix) as a match for a primer.")
@click.option('--amplicon_min_score', type=click.FloatRange(min=0, max=100), default=30, show_default=True,
              help="Minimum normalized alignment score to consider a read alignment as valid. "
                   "Normalized alignment score is defined as the Needleman-Wunch alignment score divided by "
                   "the maximum possible score. Below this alignment threshold, reads are discarded.")
@click.option('--translocation_amplicon_min_score', type=click.FloatRange(min=0, max=100), default=80, show_default=True,
              help="Minimum alignment score to consider a read with primer inconsistency as a possible translocation. "
                   "Should be higher than --amplicon_min_score, because translocations reads are noisier."
                   "Score is normalized between 0 (not even one bp match) to 100 (read is identical to")
@click.option("--min_editing_activity", type=click.FloatRange(), default=0.1, show_default=True,
              help="Minimum editing activity (%). Sites with editing activity lower than the minimum, "
                   "will be discarded from the translocation detection.")
@click.option("--translocation_p_value", type=click.FloatRange(), default=0.05, show_default=True,
              help="Translocations statistical significance level. This threshold is applied on the corrected p_value,"
                   "FDR (false discovery rate).")
@click.option('--disable_translocations', is_flag=True, default=False, show_default=True,
              help="Disable translocations detection")
@click.option("--override_noise_estimation", is_flag=True, default=False, show_default=True,
              help="Override noise estimation with default q parameter from crispector_config file. "
                   "It's advisable to set this flag for experiment with a low number of off-target sites (<5). "
                   "q is defined as the probability of an indel to occur through an edit event. Check CRISPECTOR "
                   "paper for more details.")
@click.option("--confidence_interval", type=click.FloatRange(min=0, max=1), default=0.95, show_default=True,
              help="Confidence interval for the evaluated editing activity")
@click.option('--enable_substitutions', is_flag=True, default=False, show_default=True,
              help="Enable substitutions events for the quantification of edit events")
@click.option("--suppress_site_output",  is_flag=True, default=False, show_default=True,
              help="Do not create plots for sites (save memory and runtime)")
@click.option('--keep_intermediate_files', is_flag=True, default=False, show_default=True, required=True,
              help="Keep intermediate files for debug purposes")
@click.option('--verbose', is_flag=True, default=False, show_default=True,
              help="Higher verbosity")
def cli(**kwargs):
    """Accurate estimation of off-target editing activity from comparative NGS data"""
    kwargs["command_used"] = ' '.join(sys.argv)
    if kwargs["report_output"] is None:
        kwargs["report_output"] = os.path.abspath(os.getcwd())

    # Run crispector
    sys.exit(run(**kwargs))


if __name__ == "__main__":
    sys.exit(cli())
