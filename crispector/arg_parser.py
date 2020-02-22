"""Console script for crispector."""
import sys
import click
from crispector_main import run

@click.command()
@click.option('--tx_in1', '-t_r1', type=click.Path(exists=True), help="Treatment read 1 input path or treatment merged FASTQ file")
@click.option('--tx_in2', '-t_r2', type=click.Path(exists=True), help="Treatment read 2 input path, if FASTQ's files aren't merged [OPTIONAL] ")
@click.option('--mock_in1', '-m_r1', type=click.Path(exists=True), help="Mock read 1 input path or mock merged FASTQ file")
@click.option('--mock_in2', '-m_r2', type=click.Path(exists=True), help="Mock read read 2 input path, if FASTQ's files aren't merged [OPTIONAL] ")
@click.option("--experiment_config", '-c', type=click.Path(exists=True), required=True,
              help="A CSV (Comma-separated values‏) file with the experiment data. Table has 11 columns: "
                   "SiteName, AmpliconReference, sgRNA, OnTarget, ForwardPrimer, ReversePrimer, TxInput1Path"
                   "TxInput2Path, MockInput1Path, MockInput2Path, DonorReference.\n"
              "The first 4 columns are required, the rest are optional. Header should be specified by the above order."
              "Please check the README file for further details and examples.") #TODO - sgRNA only 5'->3'
@click.option('--report_output', '-o', type=click.Path(), default="CRISPECTOR", show_default=True, required=True,
              help="Path to output folder (string)")
@click.option("--cut_site_position", type=click.INT, default=-3, show_default=True,
              help="Cut-site position with respect to the 3' end of the provided sgRNA sequence. Note, the sgRNA sequence must be entered without the PAM.")
@click.option("--crispector_config", type=click.Path(),
              help="Path YAML configuration file. See README on GitHub (####) for more details.") #TODO - add description
@click.option('--fastp_options_string', type=click.STRING, default="-w 2", help="Try \"fastp --help\" for more details")
@click.option("--min_num_of_reads", type=click.INT, default=500, show_default=True,
              help="Minimum number of reads (per site) to evaluate edit events")
@click.option("--min_read_length", type=click.INT, default=40, show_default=True,
              help="Filter out any read shorter than min read length")
@click.option('--amplicon_min_score', type=click.FloatRange(min=0, max=100), default=30, show_default=True,
              help="Minimum alignment score to consider a read alignment as valid."
                   "Score is normalized between 0 (not even one bp match) to 100 (read is identical to"
                   "reference). Below this alignment threshold, reads are discarded."
                   "This is useful for filtering erroneous reads that do not align to their target amplicon.")
@click.option('--translocation_amplicon_min_score', type=click.FloatRange(min=0, max=100), default=80, show_default=True,
              help="Minimum alignment score to consider a read alignment as valid. Should be higher than --amplicon_min_score"
                   "because translocations reads are noisier."
                   "Score is normalized between 0 (not even one bp match) to 100 (read is identical to"
                   "reference). Below this alignment threshold, reads are discarded."
                   "This is useful for filtering erroneous reads that do not align to their translocation amplicon.")
@click.option("--override_binomial_p", is_flag=True, default=False, show_default=True,
              help="Override binomial coin estimation with default value from config file. It's advisable to set"
              "this flag for low number of Off-Target sites (<5)")
@click.option("--confidence_interval", type=click.FloatRange(min=0, max=1), default=0.95, show_default=True,
              help="Confidence interval for the evaluated editing activity")
@click.option("--editing_threshold", type=click.FloatRange(), default=0.1, show_default=True,
              help="The editing activity threshold (%). Below this threshold, editing activity won't be"
                   "displayed at final plots and possible translocation will be filtered out.")
@click.option("--translocation_p_value", type=click.FloatRange(), default=0.05, show_default=True,
              help="Translocations statistical significance level. This threshold is applied on the corrected p_value,"
                   "FDR (false discovery rate).")
@click.option("--suppress_site_output",  is_flag=True, default=False, show_default=True,
              help="Do not dump plots and tables for sites")
@click.option('--disable_translocations', is_flag=True, default=False, show_default=True,
              help="Disable translocations search")
@click.option('--enable_substitutions', is_flag=True, default=False, show_default=True, help="Enable substitutions"
              "events for the quantification of edit events")
@click.option('--keep_intermediate_files', is_flag=True, default=False, show_default=True, required=True,
              help="Keep intermediate files for debug purposes")
@click.option('--verbose', is_flag=True, default=False, show_default=True, help="Higher verbosity")
def main(**kwargs):
    """CRISPECTOR - Console script"""

    # Run crispector
    run(**kwargs)

if __name__ == "__main__":
    sys.exit(main())
