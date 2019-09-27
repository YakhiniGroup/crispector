from constants_and_types import Pr, IndelType, C_TX, C_MOCK, ON_TARGET, SITE_NAME, CUT_SITE
from utils import Configurator, Logger
from typing import Dict, List
from modification_tables import ModificationTables
from modification_types import ModificationTypes
import numpy as np
import pandas as pd
import os

# TODO - upgrade function.
# TODO - decide if need to override values here with other user parameters.
def compute_binom_p(tables: Dict[str ,ModificationTables], modifications: ModificationTypes,
                    override_coin: bool, ref_df, output) -> List[Pr]:

    logger = Logger.get_logger()
    binom_p_l = []
    if override_coin:
        cfg = Configurator.get_cfg()
        p = cfg["default_binom_p"]
        binom_p_l = modifications.size*[p]
        return binom_p_l

    coin_d = dict()
    percents = [50, 60, 70, 80, 90, 95]

    for table_idx in range(modifications.size):
        alpha = dict()
        for percent in percents:
            alpha[percent] = []

        for site_name, site_table in tables.items():
            # TODO - support in multi on-target sites
            table = site_table.tables[table_idx]
            on_target = ref_df.loc[site_name, ON_TARGET]
            cut_site = ref_df.loc[site_name, CUT_SITE]
            if on_target:
                n_on = site_table.n_reads_tx
                # For insertion - edit on t he cut-site
                if modifications.types[table_idx] == IndelType.INS:
                    e_on = site_table.tables[table_idx][C_TX, cut_site]
                # For deletion/substitutions - max between left and right to the cut-site
                else:
                    e_on = np.max(site_table.tables[table_idx][C_TX, cut_site:cut_site+1])
            else:
                # mock_indel_v = table[C_MOCK, :][table[C_MOCK, :] > 0]
                alpha_v = table[C_MOCK, 10:-10] / site_table.n_reads_mock
                alpha_v = alpha_v[alpha_v > 0]
                for percent in percents:
                    if len(alpha_v) != 0:
                        alpha[percent].append(np.percentile(alpha_v, percent))

        for percent in percents:
            alpha[percent] = np.array(alpha[percent])

        # alpha_w_zeros_v = np.concatenate(alpha_list, axis=0)
        # alpha_v = alpha_w_zeros_v[alpha_w_zeros_v > 0]
        # alpha_opt = np.array([np.mean(alpha_w_zeros_v), np.mean(alpha_v), np.percentile(alpha_v, 50), np.percentile(alpha_v, 90)])

        # binom_p_l.append(p_v[1])
        name = modifications.name_at_idx(table_idx)
        coin_d[name] = dict()
        coin_d[name]["e_on"] = e_on
        coin_d[name]["n_on"] = n_on
        for percent_a in percents:
            for percent_b in percents:
                a_mock = np.percentile(alpha[percent_a], percent_b)
                coin_d[name]["p_{}_{}".format(percent_a, percent_b)] = e_on / (e_on + n_on * a_mock)

        coin_df = pd.DataFrame.from_dict(coin_d, orient='index')
        coin_df["IndelType"] = coin_df.index
        coin_df.to_csv(os.path.join(output, "coin_results.csv"), index=False)
    # TODO - delete this code
    cfg = Configurator.get_cfg()
    p = cfg["default_binom_p"]
    binom_p_l = modifications.size*[p]
    return binom_p_l

