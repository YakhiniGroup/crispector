from collections import defaultdict

from constants_and_types import Pr, IndelType, C_TX, C_MOCK, ON_TARGET, CUT_SITE, SITE_NAME, R_PRIMER, F_PRIMER
from utils import Configurator, Logger
from typing import Dict, List
from modification_tables import ModificationTables
from modification_types import ModificationTypes
import numpy as np
import pandas as pd
import os

def compute_binom_p(tables: Dict[str ,ModificationTables], modifications: ModificationTypes,
                    override_coin: bool, ref_df, output) -> Dict[str, List[Pr]]:
    """
    Compute binomial probability for each site.
    :param tables: Modification tables
    :param modifications: Modification Types
    :param override_coin:
    :param ref_df:
    :param output:
    :return: Dict[SITE_NAME, LIST of probability for each modification type]
    """
    cfg = Configurator.get_cfg()
    binom_p_d = dict()
    default_p = cfg["default_binom_p"]
    average_p = cfg["average_binom_p"]

    # TODO - Delete this section
    coin_d = defaultdict(list)
    percents = [50, 90, 95, 99]
    on_site_name = ref_df.loc[ref_df[ON_TARGET], SITE_NAME].values[0]
    # TODO - mixed, snps..

    # Use default_p if override_coin is True is on_target has low number of reads
    if override_coin or (not on_site_name in  tables):
        for site_name in tables.keys():
            binom_p_d[site_name] = modifications.size*[default_p]
        return binom_p_d

    n_on = tables[on_site_name].n_reads_tx

    p_d = dict() # probability dict
    for indel_type in [IndelType.DEL, IndelType.INS]:
        e_on = 0
        alpha_d = dict()
        p_d[indel_type] = dict()
        for table_idx in range(modifications.size):
            if indel_type != modifications.types[table_idx]:
                continue

            for site_name, site_table in tables.items():
                table = site_table.tables[table_idx]
                cut_site = ref_df.loc[site_name, CUT_SITE]
                f_primer = len(ref_df.loc[site_name, F_PRIMER])
                r_primer = len(ref_df.loc[site_name, R_PRIMER])

                # Calculate e_on - number of modifications ("edits") on the cut-site
                if site_name == on_site_name:
                    # For insertion - edit on the cut-site
                    if indel_type == IndelType.INS:
                        e_on += site_table.tables[table_idx][C_TX, cut_site]
                    # For deletion - max between left and right to the cut-site
                    else:
                        e_on += np.max(site_table.tables[table_idx][C_TX, cut_site:cut_site+1])

                    # Don't compute alpha from on-target due to possible contamination
                    continue

                # Compute alpha - percentage of "noise" on each position
                if site_name in alpha_d:
                    alpha_d[site_name] += table[C_MOCK, f_primer:-r_primer] / site_table.n_reads_mock
                else:
                    alpha_d[site_name] = table[C_MOCK, f_primer:-r_primer] / site_table.n_reads_mock

        # Estimate alpha from all alpha values
        alpha_est_d = dict()
        for percent_a in percents:
            alpha_est_d[percent_a] = np.zeros(len(alpha_d))
            for idx, alpha_v in enumerate(alpha_d.values()):
                alpha_v = alpha_v[alpha_v > 0]
                if len(alpha_v) != 0:
                    alpha_est_d[percent_a][idx] = np.percentile(alpha_v, percent_a)


        # Estimate Binomial probability for each site
        for site_name, site_table in tables.items():
            coin_d["indel_type"].append(indel_type.name)
            coin_d["e_on"].append(e_on)
            coin_d["n_on"].append(n_on)
            coin_d["site_name"].append(site_name)
            for percent_a in percents:
                for percent_b in percents:
                    a_mock = np.percentile(alpha_est_d[percent_a], percent_b)
                    n_mock = site_table.n_reads_mock
                    n_tx = site_table.n_reads_tx
                    binom_p = e_on / (e_on + (a_mock * n_on * n_mock / n_tx))
                    binom_p = (binom_p + average_p) / 2
                    coin_d["p_{}_{}".format(percent_a, percent_b)].append(binom_p)
                    if (percent_a == 95) and (percent_b == 95):
                        p_d[indel_type][site_name] = binom_p

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
                binom_p_d[site_name].append(default_p)

    coin_df = pd.DataFrame.from_dict(coin_d, orient='columns')
    coin_df.to_csv(os.path.join(output, "coin_results.csv"), index=False)
    return binom_p_d

