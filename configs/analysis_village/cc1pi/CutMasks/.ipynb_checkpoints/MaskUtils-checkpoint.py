
import numpy as np
import pandas as pd
   

    
from analysis_village.cc1pi.Constants import CTE as CTE
group_levels = ['__ntuple', 'entry', 'rec.slc..index']
def build_event_cumulative_masks(evt_df, sideband = "shower"):
    # Define which BDT mask to use
    cut_sequence = [
        ("cosmic",      evt_df.slc.cut.obvious_cosmic),
        ("t0",          evt_df.slc.cut.t0),
        ("FV",          evt_df.slc.cut.inside_FV),
        ("nu_score",    evt_df.slc.nu_score > CTE.min_nu_score),
        ("track",       evt_df.slc.cut.track),
        ("chi2",        evt_df.slc.cut.MIP_candidates),
        ("shower",      evt_df.slc.cut.shower),
        ("angle",       evt_df.slc.cut.angle),
        ("proton_BDT",  evt_df.slc.cut.proton_BDT),
        #("containment", evt_df.slc.cut.containment),
        ("containment", (evt_df.slc.cut_var.n_exiting_pfps == 0) & (evt_df.slc.cut.no_high_yz == True)),
        ("TPC_containment",  evt_df.slc.cut.TPC_containment ),
        ("michel",      evt_df.slc.cut.michel),
        ("extra_pion",  evt_df.slc.cut.extra_pion),
        ("energy",      evt_df.slc.cut.energy),
    ]

    
    if sideband == "proton": 
        cut_sequence = [
            ("cosmic",      evt_df.slc.cut.obvious_cosmic),
            ("t0",          evt_df.slc.cut.t0),
            ("FV",          evt_df.slc.cut.inside_FV),
            ("nu_score",    evt_df.slc.nu_score > CTE.min_nu_score),
            ("track",       evt_df.slc.cut.track),
            ("chi2",        evt_df.slc.cut.MIP_candidates),
            ("shower",      evt_df.slc.cut.shower),
            ("angle",       evt_df.slc.cut.angle),
            ("proton_BDT",  evt_df.slc.cut.proton_BDT_sideband),
            ("containment", (evt_df.slc.cut_var.n_exiting_pfps == 0) & (evt_df.slc.cut.no_high_yz == True)),
            ("TPC_containment",  evt_df.slc.cut.TPC_containment),
            ("michel",      evt_df.slc.cut.michel),
            ("extra_pion",  evt_df.slc.cut.extra_pion),
            ("energy",      evt_df.slc.cut.energy),
            #("TPC_containment",  evt_df.slc.cut.TPC_containment),
        ]
    elif sideband == "two_pions": 
        cut_sequence = [
            ("cosmic",      evt_df.slc.cut.obvious_cosmic),
            ("t0",          evt_df.slc.cut.t0),
            ("FV",          evt_df.slc.cut.inside_FV),
            ("nu_score",    evt_df.slc.nu_score > CTE.min_nu_score),
            ("track",       evt_df.slc.cut.track),
            ("chi2",        evt_df.slc.cut_var.n_MIP_candidates > 2),
            ("shower",      evt_df.slc.cut.shower),
            ("angle",       evt_df.slc.cut.angle),
            ("proton_BDT",  evt_df.slc.cut.proton_BDT_2pi),
            ("containment", (evt_df.slc.cut_var.n_exiting_pfps == 0) & (evt_df.slc.cut.no_high_yz == True)),
            ("TPC_containment",  evt_df.slc.cut.TPC_containment),
            ("michel",      evt_df.slc.cut.michel),
            ("energy",      evt_df.slc.cut.energy),
            #("TPC_containment",  evt_df.slc.cut.TPC_containment),
        ]

    # 2️⃣ Build cumulative masks
    cumulative_masks = {}
    current_mask = None

    for name, mask in cut_sequence:

        if current_mask is None:
            current_mask = mask.copy()
        else:
            current_mask = current_mask & mask

        cumulative_masks[name] = current_mask.copy()

    return cumulative_masks

GROUP_LEVELS = ['__ntuple', 'entry', 'rec.slc..index']
WGT_COL = ('slc', 'wgt', '', '', '', '')

def get_n_evt(df, mask=None, use_weight=True):
    if WGT_COL not in df.columns:
        raise ValueError("Weight column not found")

    wgt = df[WGT_COL] if mask is None else df.loc[mask, WGT_COL]
    grouped = wgt.groupby(level=GROUP_LEVELS)

    if not use_weight:
        return grouped.ngroups

    return grouped.first().sum()


'''
def get_n_evt(df, use_weight=True):
    # Use .codes on each level directly — no index copy
    level_arrays = [df.index.get_level_values(i) for i in range(3)]
    
    if not use_weight:
        # Combine levels into a single structured array for unique counting
        # Much cheaper than creating a new MultiIndex
        combined = list(zip(*level_arrays))  # avoids full MultiIndex rebuild
        return len(set(combined))

    wgt_col = ('slc', 'wgt', '', '', '', '')
    if wgt_col not in df.columns:
        raise ValueError("Weight column not found")

    # Work on a minimal 2-column frame: avoid copying the whole df
    slim = df[[wgt_col]].copy(deep=False)  # shallow copy, no data duplication
    slim.index = pd.MultiIndex.from_arrays(level_arrays)  # drop unused levels in-place
    
    # ~first per event: use groupby on the slimmed index
    return slim[wgt_col].groupby(level=[0, 1, 2]).first().sum()
'''