from crispector_constants import C_MOCK, ON_TARGET, SITE_NAME, C_TX, CUT_SITE
from crispector_types import Pr, IndelType
from crispector_utils import Configurator, Logger
from typing import Dict, List
from modification_tables import ModificationTables
from modification_types import ModificationTypes
import numpy as np


# TODO - upgrade function.
# TODO - decide if need to override values here with other user parameters.
def compute_binom_p(tables: Dict[str ,ModificationTables], modifications: ModificationTypes,
                    override_coin: bool, ref_df) -> List[Pr]:

    logger = Logger.get_logger()
    binom_p_l = []
    if override_coin:
        cfg = Configurator.get_cfg()
        p = cfg["default_binom_p"]
        binom_p_l = modifications.size*[p]
        return binom_p_l

    for table_idx in range(modifications.size):
        alpha_list = []
        for site_name, site_table in tables.items():
            # TODO - support in multi on-target sites
            table = site_table.tables[table_idx]
            on_target = ref_df.loc[ref_df[SITE_NAME] == site_name, ON_TARGET].values[0]
            cut_site = ref_df.loc[ref_df[SITE_NAME] == site_name, CUT_SITE].values[0]
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
                alpha_list.append(table[C_MOCK, 10:-10] / site_table.n_reads_mock)

        alpha_w_zeros_v = np.concatenate(alpha_list, axis=0)
        alpha_v = alpha_w_zeros_v[alpha_w_zeros_v > 0]
        alpha_opt = np.array([np.mean(alpha_w_zeros_v), np.mean(alpha_v), np.percentile(alpha_v, 50), np.percentile(alpha_v, 90)])
        p_v = e_on / (e_on + n_on*alpha_opt)

        logger.debug("{},{},{},{:.4f},{:.4f},{:.4f},{:.4f},{:.4f},{:.4f},{:.4f},{:.4f}".format(table_idx, e_on, n_on,
                     alpha_opt[0], alpha_opt[1], alpha_opt[2], alpha_opt[3], p_v[0], p_v[1], p_v[2], p_v[3]))
        binom_p_l.append(p_v[1])

    return binom_p_l
