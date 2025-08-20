# MIT License
# Simple materials library for Meep (approximate dispersive params or n,k tables)
# For RF and IR regimes; extend as needed.

import numpy as np

def constant_eps(n_real):
    return n_real**2

RF_MATERIALS = {
    "FR4": {"eps": 4.2, "tan_delta": 0.02},
    "RO4350B": {"eps": 3.66, "tan_delta": 0.0037},
    "AIR": {"eps": 1.0, "tan_delta": 0.0},
}

IR_MATERIALS = {
    "PDMS": {"n": 1.41},
    "SiO2": {"n": 1.45},
    "TiO2": {"n": 2.40},
    "AIR": {"n": 1.0},
}

def get_rf_eps(material_name):
    return RF_MATERIALS[material_name]["eps"]

def get_ir_n(material_name):
    return IR_MATERIALS[material_name]["n"]
