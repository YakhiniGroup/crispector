# -*- coding: utf-8 -*-

"""Console script for crispector."""
import sys
import click
from crispector_main import run
import os

# TODO - add a config from user + window_size + effective qualification window
@click.command()
@click.option('--tx_in1', '-t_r1', type=click.Path(exists=True), required=True,
              help="Treatment read 1 input path (string)")
@click.option('--tx_in2', '-t_r2', type=click.Path(exists=True), help="Treatment read 2 input path (string)")
@click.option('--mock_in1', '-m_r1', type=click.Path(exists=True), required=True,
              help="Mock read 1 input path (string)")
@click.option('--mock_in2', '-m_r2', type=click.Path(exists=True), help="Mock read 2 input path (string)")
@click.option("--amplicons_csv", '-a', type=click.Path(exists=True), required=True,
              help="A CSV (Comma-separated values‏) file with all amplicon sequences. Table has 4 columns:\
              amplicon_name (string) - A unique identifier.\
              amplicon_sequence (4 letter string) - The amplicon reference sequence.\
              sgRNA_sequence (4 letter string) - The sgRNA sequence without the PAM. Should be inside the\
              amplicon_sequence.\
              on_target (Bool) - Indicate if the site is on-targer or off-target.\
              All fields are required. No header should be specified.\
              Please check the README file for further details and examples.")
@click.option('--report_output', '-o', type=click.Path(), default="CRISPECTOR", show_default=True,
              help="Output folder path (string)")
@click.option('--fastp_options_string', type=click.STRING, default="", help="Try \"fastp --help\" for more details")
# TODO - Add --length_required 40 to filter anything under 40 (15 is the default). see fastp TODO
@click.option('--fastp_threads', type=click.INT, default=2, help="fastp worker thread number")
@click.option('--override_fastp', is_flag=True, default=False, show_default=True, # TODO -Change default to 1
              help="Override fastp and require merged fastq files (pre-processing is necessary).\
                    Set paths to merged fastq files at --tx_in1 and --mock_in1.\
                    Can't be used with --fastp_options_string")
@click.option('--verbose', is_flag=True, default=False, show_default=True, help="Higher verbosity")
@click.option('--keep_fastp_output', is_flag=True, default=False, show_default=True, help="Keep fastp output directory")
@click.option("--min_num_of_reads", type=click.INT, default=500, show_default=True,
              help="Minimum number of reads (per site) to evaluate edit events")
@click.option("--cut_site_position", type=click.INT, default=-3, show_default=True,
              help="Cut-site position relative to PAM (minus sign for upstream)")
@click.option("--config", '-c', type=click.Path(),
              help="Path YAML configuration file. See README on GitHub (####) for more details.")
@click.option("--override_binomial_p", is_flag=True, default=False, show_default=True,
              help="Override binomial coin estimation with default value from config file. It's advisable to set"
              "this flag for low number of sites (< ############)") # TODO - fix this description
@click.option("--confidence_interval", type=click.FloatRange(min=0, max=1), default=0.95, show_default=True,
              help="Confidence interval for the evaluated editing activity")
@click.option("--editing_threshold", type=click.FloatRange(), default=0.1, show_default=True,
              help="The editing activity threshold (%). Below this threshold, editing activity won't be"
                   "displayed at final plots and possible translocation will be filtered out.")
@click.option("--translocation_p_value", type=click.FloatRange(), default=0.05, show_default=True,
              help="Translocations statistical significance level. This threshold is applied on the corrected p_value,"
                   "FDR (false discovery rate).")
@click.option("--suppress_site_output",  is_flag=True, default=False, show_default=True,
              help="Do not dump plots and reads for all sites")
# TODO - this values and description should be changed.
@click.option('--amplicon_min_score', type=click.FloatRange(min=0, max=100), default=30, show_default=True,
              help="Minimum alignment score to consider a read alignment to a specific amplicon reference sequence."
                   "Score is normalized between 0 (not even one bp match) to 100 (the read is identical to"
                   "the reference). Below this alignment threshold, reads are discarded."
                   "This is useful for filtering erroneous reads that do not align to any target amplicon.")
@click.option('--translocation_amplicon_min_score', type=click.FloatRange(min=0, max=100), default=80, show_default=True,
              help="Alignment minimum score for translocations reads. Default value is higher due to higher level of" 
                   "Noise in translocations reads.")
@click.option('--experiment_name', type=str, default=" ", show_default=True,
              help="Experiment name as will be reported in CRISPECTOR plots")
@click.option('--allow_translocations', is_flag=True, default=True, show_default=True,
              help="Add something") # TODO - add description
@click.option("--min_read_length", type=click.INT, default=40, show_default=True,
              help="Filter out any read shorter than min read length")
@click.option("--max_error_on_primer", type=click.INT, default=8, show_default=True,
              help="") # TODO - change to  max_edit distance on primer
@click.option('--enable_substitutions', is_flag=True, default=False, show_default=True, help="Enable substitutions"
              "events for the quantification of edit events")
@click.option('--ambiguous_cut_site_detection', is_flag=True, default=True, show_default=True,
              help="Detect ambiguous cut-site (e.g. PAM is GGG, so cut-site can be shift one base to the left."
              "Set False if this isn't a CAS9 experiment")
@click.option('--debug', is_flag=True, default=False, show_default=True,
              help="Delete...")
@click.option('--alignment_input',type=click.Path(), required=False, help="") # TODO - delete
@click.option('--table_input', type=click.Path(), required=False, help="") # TODO - delete
def main(**kwargs):
    """CRISPECTOR - Console script"""
    override_fastp = kwargs["override_fastp"]
    tx_in2 = kwargs["tx_in2"]
    mock_in2 = kwargs["mock_in2"]
    report_output = kwargs["report_output"]

    # Input verification
    if override_fastp:
        if (tx_in2 is not None) or (mock_in2 is not None):
            raise click.BadOptionUsage(override_fastp,
                                       "--tx_in2 and --mock_in2 can't be set when override_fastp is used!")
    else:
        if tx_in2 is None:
            raise click.BadOptionUsage(tx_in2, "--tx_in2 is missing!")

        if mock_in2 is None:
            raise click.BadOptionUsage(mock_in2, "--mock_in2 is missing!")

    # Create output folder
    if not os.path.exists(report_output):
        os.makedirs(report_output)

    # Run crispector
    run(**kwargs)

if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
