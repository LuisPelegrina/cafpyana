
import numpy as np
import pandas as pd
   

group_levels = ['__ntuple', 'entry', 'rec.slc..index']
    
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
        ("containment", evt_df.slc.cut_var.n_exiting_pfps == 0),
        ("TPC_containment",  evt_df.slc.cut.TPC_containment),
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
            ("containment", evt_df.slc.cut_var.n_exiting_pfps == 0),
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
            ("containment", evt_df.slc.cut_var.n_exiting_pfps == 0),
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


def get_n_evt(df, use_weight=True):
    # event identifier = first two index levels
    evt_index = df.index.droplevel(list(df.index.names[3:]))

    if not use_weight:
        return evt_index.nunique()

    wgt_col = ('slc','wgt','','','','')

    if wgt_col not in df.columns:
        raise ValueError("Weight column not found")

    # select the weight column first, then group
    weights = df[wgt_col].groupby(evt_index).first()

    return weights.sum()