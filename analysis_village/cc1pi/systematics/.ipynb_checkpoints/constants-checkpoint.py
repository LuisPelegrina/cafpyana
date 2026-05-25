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

'''
def get_xsec_unit():
    # ==== xsec unit calculation ====
    # TODO: z-dependence?
    # flux file, units: /m^2/10^6 POT 
    # 50 MeV bins
    fluxfile = "/exp/sbnd/data/users/munjung/flux/sbnd_original_flux.root"
    flux = uproot.open(fluxfile)
    numu_flux = flux["flux_sbnd_numu"].to_numpy()
    bin_edges = numu_flux[1]
    flux_vals = numu_flux[0]

    if plot:
        fig, ax = plt.subplots()
        plt.hist(bin_edges[:-1], bins=bin_edges, weights=flux_vals, histtype="step", linewidth=2, color="C0")
        plt.xlim(0, 3)
        plt.xlabel("Neutrino Energy [GeV]")
        plt.ylabel("Flux [/m$^{2}$/10$^{6}$ POT]")
        plt.title("SBND $\\nu_\\mu$ Flux")
        plt.savefig("sbnd-flux.pdf", bbox_inches='tight')

    # get integrated flux
    # TODO: refactor this to be always consistent with what's used in the ana notebooks
    _mc_tot_pot = data_tot_pot
    print("mc_tot_pot: %.3e" %(_mc_tot_pot))

    integrated_flux = _mc_tot_pot * flux_vals.sum() / (1e4  * 1e6) # to cm2 # to POT
    print("Integrated flux: %.3e" % integrated_flux)

    V_SBND = 380 * 380 * 470 # cm3, the active volume of the detector 
    NTARGETS = RHO * V_SBND * N_A / M_AR
    print("# of targets: ", NTARGETS)

    xsec_unit = 1 / (integrated_flux * NTARGETS)
    # TODO: fix scalar overflow error in python v3.10+
    if xsec_unit == 0:
        print("XSEC_UNIT is 0, setting to 1e-38")
        xsec_unit = 1e-38
    print("xsec unit: ", xsec_unit)
    return xsec_unit
'''

plot = True
data_tot_pot = 4.516288e+18
data_gates = 948132
def get_xsec_unit():
    # ==== xsec unit calculation ====
    # TODO: z-dependence?
    # flux file, units: /m^2/10^6 POT 
    # 50 MeV bins
    fluxfile = "/exp/sbnd/data/users/munjung/flux/sbnd_original_flux.root"
    flux = uproot.open(fluxfile)
    numu_flux = flux["flux_sbnd_numu"].to_numpy()
    bin_edges = numu_flux[1]
    flux_vals = numu_flux[0]

    if plot:
        fig, ax = plt.subplots()
        plt.hist(bin_edges[:-1], bins=bin_edges, weights=flux_vals, histtype="step", linewidth=2, color="C0")
        plt.xlim(0, 3)
        plt.xlabel("Neutrino Energy [GeV]")
        plt.ylabel("Flux [/m$^{2}$/10$^{6}$ POT]")
        plt.title("SBND $\\nu_\\mu$ Flux")
        plt.savefig("sbnd-flux.pdf", bbox_inches='tight')

    # get integrated flux
    # TODO: refactor this to be always consistent with what's used in the ana notebooks
    _mc_tot_pot = data_tot_pot
    print("mc_tot_pot: %.3e" %(_mc_tot_pot))

    integrated_flux = _mc_tot_pot * flux_vals.sum() / (1e4  * 1e6) # to cm2 # to POT
    print("Integrated flux: %.3e" % integrated_flux)
    print("Integrated flux: %.3e" % flux_vals.sum())
   
    V_SBND = (
        (400 - 2*CTE.min_distance_to_wall_x_y)
        *(400 - 2*CTE.min_distance_to_wall_x_y)
        *(500 - CTE.min_distance_to_first_z_wall - CTE.min_distance_to_last_z_wall)
    )# cm3, the active volume of the detector 
    NTARGETS = RHO * V_SBND * N_A / M_AR
    
    print("# of targets: ", NTARGETS)

    print(integrated_flux)
    
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