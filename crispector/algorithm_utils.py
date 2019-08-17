from crispector_types import Pr
from crispector_utils import Configurator
from typing import Dict, List
from modification_tables import ModificationTables
from modification_types import ModificationTypes


# TODO - upgrade function.
# TODO - decide if need to override values here with other user parameters.
def compute_binom_p(tables: Dict[str ,ModificationTables], modifications: ModificationTypes,
                    override_coin: bool) -> List[Pr]:
    binom_p_l = []
    if override_coin:
        cfg = Configurator.get_cfg()
        p = cfg["default_binom_p"]
        binom_p_l = modifications.size*[p]

    return binom_p_l
