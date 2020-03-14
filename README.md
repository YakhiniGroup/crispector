<img src="https://github.com/iamit87/crispector/blob/master/crispector/report/html_templates/crispector_logo.jpg" height="100" />

#
CRISPECTOR is a software package designed to support the detection, evaluation, and quantification of on and off-target genome editing activity. CRISPECTOR accepts FASTQ files resulting from running treatment vs. mock experiments followed by multiplex-PCR and NGS. The tool analyzes the NGS input and applies statistical modeling to determine and quantify NHEJ edit activity at every interrogated locus as well as adverse translocation activity in all relevant pairs of loci.

Briefly,  CRISPECTOR assigns each read in the treatment and mock FASTQ files to a specific locus of interest or a putative translocation. Then, a Bayesian inference classifier accurately estimates the NHEJ editing activity, and a hypergeometric test is performed to detect translocation reads.

 **CRISPECTOR Workflow:**  
  <img src="https://github.com/iamit87/crispector/blob/master/CRISPECTOR_workflow.png" />
  
* [Installation](#Installation)
* [Usage](#Usage)
* [citation](#examples-of-report)

# Installation


# Usage

CRISPECTOR is designed to run on two multiplex-PCR experiments, treatment and mock. In default, CRISPECTOR demultiplexes the treatment and mock reads. Namely, CRISPECTOR assigns each read to its locus. 
There is also an option to run CRISPECTOR on an already demultiplexed data (where each locus has a FASTQ file contains all its reads).  
Please note that in both modes, adapters need to be trimmed in a pre-processing step. 

## Usage - Default mode
CRISPECTOR requires three parameters: 
1. Treatment input sequences in the form of FASTQ files. Given by the  `-t_r1`  and  `-t_r2`.  If the input is already pair-end merged or is a single-end, then omit `-t_r2`. FASTQ files can be gzip-compressed.
2. Mock input sequences in the form of FASTQ files. Given by the  `-m_r1`  and  `-m_r2`.  If the input is already pair-end merged, then omit `-m_r2`.  FASTQ files can be gzip-compressed.
3.  An experiment config file ( `-c`). The experiment description in a CSV (Comma Separated Values‏) format. The table has 11 columns:
	-  **SiteName** [REQUIRED] - an identifier for the reference locus. 
	- **AmpliconReference**  [REQUIRED] - amplicon sequence used for the experiment (5'->3').
	- **sgRNA** [REQUIRED] - sgRNA sequence for each locus site. sequence should be supplied without the PAM sequence and without insertions or deletions.
	-  **OnTarget**  [REQUIRED] - a boolian indicating if the site is on-target (`True`) or off-target (`False`). 
	- **ForwardPrimer** [Optional] - forward primers as were used in the experiment. If not supplied, primers are inffered from the amplicon reference.
	- 	**ReversePrimer** [Optional] - reverse primers as were used in the experiment. If not supplied, primers are inffered from the amplicon reference.
	- **TxInput1Path** - Leave empty for experiment with multiplexed input. 
	- **TxInput2Path** - Leave empty for experiment with multiplexed input. 
	- **MockInput1Path** - Leave empty for experiment with multiplexed input. 
	- **MockInput2Path** - Leave empty for experiment with multiplexed input. 
	- **DonorReference** - If experiment is designed with HDR, then insert amplicon sequence in the on-target row. Note that editing HDR activity isn't evaluatied by CRISPECTOR. 

	Where data is not available, leave the cell empty. If an entire collum is empty, leave the entire collum empty. 


**Command:**
```
crispector -t_r1 tx_R1.fq.gz -t_r2 tx_R2.fq.gz -m_r1 mock_R1.fq.gz -mock_r2 m_R2.fq.gz -c exp_config.csv
```
**Exmaple:**
You can download data and configuration for EMX1 experiemt (performed with [rhAmpSeq](https://eu.idtdna.com/pages/products/next-generation-sequencing/amplicon-sequencing?utm_source=google&utm_medium=cpc&utm_campaign=ga_rhampseq&utm_content=ad_group_rhampseq&gclid=Cj0KCQjw3qzzBRDnARIsAECmryqo5fO62fqk95a4PfkqES-9G07br5kdtTpjJInnYFjqYw2OxYI2gRwaAmTQEALw_wcB)) . Experiment was designed with one on-target site and 10 off-target sites. The FASTQ files contain the first 50,000 reads of the expriment. 
The compressed data can be found [here](https://github.com/iamit87/crispector/raw/master/example/EMX1_11_sites_50k_reads.zip), and it contain the following files: 
- Files in   EMX1_config.csv
- EMX1_tx_R1.fq.gz
- EMX1_tx_R2.fq.gz
- EMX1_mock_R1.fq.gz
- EMX1_mock_R2.fq.gz

```
crispector -t_r1 EMX1_tx_R1.fq.gz -t_r2 EMX1_tx_R2.fq.gz -m_r1 EMX1_mock_R1.fq.gz -m_r2 EMX1_mock_R2.fq.gz -c EMX1_config.csv
``` 
## Usage - Reads are pre-demultiplexed 

## All options

## Advanced usage

## CRISPECTOR output
CRISPECTOR generates an HTML-based report to support user interpretation and further analysis of the outcomes. The report contains plots and tables. 
For example: multiple off-target editing activity including statistical confidence indicators and translocation results. The report also contains read statistics, such as the number of assigned and aligned reads in each processing step in the different loci. In addition, plots are generated for each individual locus: distribution of edit events, classifier results for each modification type and reference position, and alignments of all edited reads, in a graphical format. Furthermore, all aligned reads in treatment and mock, whether assigned to a specific locus or determined as primer inconsistent (and therefore maybe representing a translocation event), are saved in a set of CSV files.
An example report can be found in ##########.




