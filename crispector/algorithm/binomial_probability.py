from utils.constants_and_types import Pr, IndelType, C_TX, C_MOCK, ON_TARGET, CUT_SITE, SITE_NAME, R_PRIMER, F_PRIMER, \
    BINOM_PERCENTILE, BINOM_AVERAGE_P, AmpliconDf
from utils.logger import LoggerWrapper
from utils.configurator import Configurator
from typing import Dict, List
from modifications.modification_tables import ModificationTables
from modifications.modification_types import ModificationTypes
import numpy as np


def compute_binom_p(tables: Dict[str ,ModificationTables], modifications: ModificationTypes,
                    override_coin: bool, ref_df: AmpliconDf) -> Dict[str, List[Pr]]:
    """
    Compute binomial probability for each site.
    :param tables: Modification tables
    :param modifications: Modification Types
    :param override_coin:
    :param ref_df:
    :param donor: is donor experiment flag
    :return: Dict[SITE_NAME, LIST of probability for each modification type]
    """

    logger = LoggerWrapper.get_logger()
    cfg = Configurator.get_cfg()
    binom_p_d = dict()
    default_q = cfg["NHEJ_inference"]["default_q"]  # the probability of an indel to occur through an edit event
    default_q_warning_text = "Will use default q parameter (= {}) from config file for the Bayesian NHEJ inference " \
                             "(check CRISPECTOR paper for more details).".format(default_q)
    # Use default Binomial probability when can't estimate signal from on target.
    if ref_df[ON_TARGET].sum() > 1:
        override_coin = True
        logger.warning("Experiment has Multiple on-target sites. {}".format(default_q_warning_text))
    elif ref_df[ON_TARGET].sum() == 0:
        override_coin = True
        logger.warning("Experiment doesn't have on-target site. {}".format(default_q_warning_text))
    else:
        on_site_name = ref_df.loc[ref_df[ON_TARGET], SITE_NAME].values[0]
        if on_site_name not in tables:
            override_coin = True
            logger.warning("On-target site was discarded from evaluation. {}".format(default_q_warning_text))
        elif ref_df.shape[0] == 1:
            override_coin = True
            logger.warning("Experiment doesn't have off-target sites. {}".format(default_q_warning_text))

    # Use default_p if override_coin is True
    if override_coin:
        for site_name in tables.keys():
            binom_p_d[site_name] = modifications.size*[default_q]
        return binom_p_d

    n_on = tables[on_site_name].n_reads_tx
    p_d = dict() # probability dict, for all indels and sites
    for indel_type in [IndelType.DEL, IndelType.INS]:
        e_on = 0
        site_noise_d = dict() # Keys - sites, values - mock noise vector (for each position)
        p_d[indel_type] = dict()

        # Aggregate all tables for current indel_type
        for table_idx in range(modifications.size):
            if indel_type != modifications.types[table_idx]:
                continue

            # Concatenate alpha from all sites
            for site_name, site_table in tables.items():
                table = site_table.tables[table_idx]

                # Calculate e_on - number of modifications ("edits") on the cut-site
                if site_name == on_site_name:
                    cut_site = ref_df.loc[site_name, CUT_SITE]
                    # For insertion - edit on the cut-site
                    if indel_type == IndelType.INS:
                        e_on += site_table.tables[table_idx][C_TX, cut_site]
                    # For deletion - max between left and right to the cut-site
                    else:
                        e_on += np.max(site_table.tables[table_idx][C_TX, cut_site:cut_site+1])

                    # Don't compute alpha from on-target due to possible contamination
                    continue

                # Compute alpha - percentage of "noise" on each position outside the primers
                f_primer = len(ref_df.loc[site_name, F_PRIMER])
                r_primer = len(ref_df.loc[site_name, R_PRIMER])
                if site_name in site_noise_d:
                    site_noise_d[site_name] += table[C_MOCK, f_primer:-r_primer] / site_table.n_reads_mock
                else:
                    site_noise_d[site_name] = table[C_MOCK, f_primer:-r_primer] / site_table.n_reads_mock

        # Estimate experiment noise vector from all Concatenated site noise values
        experiment_noise = np.zeros(len(site_noise_d))
        for site_idx, site_noise in enumerate(site_noise_d.values()):
            site_noise = site_noise[site_noise > 0]
            if len(site_noise) != 0:
                # Percentile on different position in the same site
                experiment_noise[site_idx] = np.percentile(site_noise, BINOM_PERCENTILE)

        # Estimate Binomial probability for each site
        for site_name, site_table in tables.items():
            alpha = np.percentile(experiment_noise, BINOM_PERCENTILE)  # Percentile on different sites
            n_mock = site_table.n_reads_mock
            n_tx = site_table.n_reads_tx
            binom_p = e_on / (e_on + (alpha * n_on * n_mock / n_tx))
            p_d[indel_type][site_name] = (binom_p + BINOM_AVERAGE_P) / 2 # Average for stability

    # Convert p_d to crispector algorithm format
    for site_name in tables.keys():
        binom_p_d[site_name] = []
        for indel_type in modifications.types:
            if indel_type in [IndelType.DEL, IndelType.INS]:
                binom_p_d[site_name].append(p_d[indel_type][site_name])
            # For mixed indels - average insertions and deletions
            elif indel_type == IndelType.MIXED:
                binom_p_d[site_name].append((p_d[IndelType.DEL][site_name] + p_d[IndelType.INS][site_name]) / 2)
            # For substitutions - use the default coin
            elif indel_type == IndelType.SUB:
                binom_p_d[site_name].append(default_q)

    return binom_p_d

