# Constants for this analysis
# Today: 2026-02-04

import uproot
import matplotlib.pyplot as plt
from makedf.constants import *
from os import path
import pandas as pd

from pyanalib.split_df_helpers import *
from pyanalib.pandas_helpers import *
from analysis_village.cc1pi.Constants import CTE as CTE
from analysis_village.cc1pi.GraphUtils.Utils import *
plot = False


DETECTOR = "SBND_nohighyz"
EPSILON = 1e-6 # for clipping distributions at bin ranges
PROTON_MASS = 0.938272
NEUTRON_MASS = 0.939565
MUON_MASS = 0.105658
PION_MASS = 0.139570

MASS_A = 22*NEUTRON_MASS + 18*PROTON_MASS - 0.34381
BE = 0.0295
MASS_Ap = MASS_A - NEUTRON_MASS + BE

# --- for integrated flux calculation ---
RHO = 1.3836  #g/cm3, liquid Ar density
N_A = 6.02214076e23 # Avogadro’s number
M_AR = 39.95 # g, molar mass of argon

plot = True
#data_tot_pot = 4.516288e+18
data_tot_pot = 4.4343664e+18
data_gates = 948132
integrated_fv_flux = 0.000151169785920044047/1e4

def get_xsec_unit():
    print("data_tot_pot: %.3e" %(data_tot_pot))

    integrated_flux = data_tot_pot * integrated_fv_flux
    print("Integrated flux: %.3e" % integrated_flux)

    
    V_SBND = (
        185 * 380 * 440 + 
        185 * 380 * 240 + 
        185 * 290 * 200
    )# cm3, the active volume of the detector 
    NTARGETS = RHO * V_SBND * N_A / M_AR
    
    print("# of targets: ", NTARGETS)

    
    flux_64 = np.float64(integrated_flux)
    targets_64 = np.float64(NTARGETS)
    
    denom = flux_64 * targets_64
    
    xsec_unit = 1. / denom
    # TODO: fix scalar overflow error in python v3.10+
    if xsec_unit == 0:
        print("XSEC_UNIT is 0, setting to 1e-38")
        xsec_unit = 1e-38
        
    print("xsec unit: ", xsec_unit)
    return xsec_unit
    
XSEC_UNIT = get_xsec_unit()