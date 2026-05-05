
import numpy as np
import pandas as pd

def InFV_strict(df):
    xmax = 190.
    zmin = 10.
    zmax = 450.
    ymax_highz = 100.
    pass_xz = (np.abs(df.slc.vertex.x) < xmax) & (df.slc.vertex.z > zmin) & (df.slc.vertex.z < zmax)
    pass_y = ((df.slc.vertex.z < 250) & (np.abs(df.slc.vertex.y) < 190.)) | ((df.slc.vertex.z > 250) & (df.slc.vertex.y > -190.) & (df.slc.vertex.y < ymax_highz))
    return pass_xz & pass_y



def cathode_crossing_pfp_mask(df):
    xmin = -CTE.min_distance_to_consider_contained
    xmax = CTE.min_distance_to_consider_contained
    
    crossing_cathode = (df.pfp.trk.start.x > xmin) & (df.pfp.trk.start.x < xmax) 
    return crossing_cathode

    
def TPC_containment_mask(df, group_levels):
    valid_df = df[df[('pfp','trk','len','','','')] > 0]
    
    vertex_x = valid_df[('slc','vertex','x','','','')].groupby(level=group_levels).first()
    
    vertex_sign = np.sign(vertex_x)
    row_vertex_sign = valid_df.index.droplevel('rec.slc.reco.pfp..index').map(vertex_sign)
    
    start_bad = np.sign(valid_df[('pfp','trk','start','x','','')]) != row_vertex_sign
    end_bad   = np.sign(valid_df[('pfp','trk','end','x','','')])   != row_vertex_sign

    violating_df = valid_df[start_bad | end_bad]
    violating_slices = violating_df.groupby(level=group_levels).size()
    invalid_slices = violating_slices[violating_slices > 0].index
     
    mask = pd.Series(
        ~df.index.droplevel('rec.slc.reco.pfp..index').isin(invalid_slices),
        index=df.index
    )
    
    cathode_crossing_df = valid_df[cathode_crossing_pfp_mask(valid_df)]
    cathode_crossing_counts = cathode_crossing_df.groupby(level=group_levels).size()
    cathode_crossing_invalid_slices = cathode_crossing_counts[cathode_crossing_counts > 0].index 
    cathode_crossing_mask = pd.Series(~df.index.droplevel('rec.slc.reco.pfp..index').isin(cathode_crossing_invalid_slices), index=df.index)

    return mask & cathode_crossing_mask
    

group_levels = ['__ntuple', 'entry', 'rec.slc..index']
    
from analysis_village.cc1pi.Constants import CTE as CTE
group_levels = ['__ntuple', 'entry', 'rec.slc..index']
def build_event_cumulative_masks(evt_df, sideband = "shower"):
    # Define which BDT mask to use
    
    cut_sequence = [
        ("cosmic",      evt_df.slc.cut.obvious_cosmic),
        ("t0",          evt_df.slc.cut.t0),
        ("FV",          InFV_strict(evt_df)),
        ("nu_score",    evt_df.slc.nu_score > CTE.min_nu_score),
        ("track",       evt_df.slc.cut.track),
        ("chi2",        evt_df.slc.cut.MIP_candidates),
        ("shower",      evt_df.slc.cut.shower),
        ("angle",       evt_df.slc.cut.angle),
        ("proton_BDT",  evt_df.slc.cut.proton_BDT),
        #("containment", evt_df.slc.cut.containment),
        ("containment", evt_df.slc.cut_var.n_exiting_pfps == 0),
        ("TPC_containment",  TPC_containment_mask(evt_df,group_levels)),
        ("michel",      evt_df.slc.cut.michel),
        ("extra_pion",  evt_df.slc.cut.extra_pion),
        ("energy",      evt_df.slc.cut.energy),
    ]
    

    if sideband == "proton": 
        cut_sequence = [
            ("cosmic",      evt_df.slc.cut.obvious_cosmic),
            ("t0",          evt_df.slc.cut.t0),
            ("FV",          InFV_strict(evt_df)),
            ("nu_score",    evt_df.slc.nu_score > CTE.min_nu_score),
            ("track",       evt_df.slc.cut.track),
            ("chi2",        evt_df.slc.cut.MIP_candidates),
            ("shower",      evt_df.slc.cut.shower),
            ("angle",       evt_df.slc.cut.angle),
            ("proton_BDT",  evt_df.slc.cut.proton_BDT_sideband),
            ("containment", evt_df.slc.cut_var.n_exiting_pfps == 0),
            ("TPC_containment",  TPC_containment_mask(evt_df,group_levels)),
            ("michel",      evt_df.slc.cut.michel),
            ("extra_pion",  evt_df.slc.cut.extra_pion),
            ("energy",      evt_df.slc.cut.energy),
            #("TPC_containment",  evt_df.slc.cut.TPC_containment),
        ]
    elif sideband == "two_pions": 
        cut_sequence = [
            ("cosmic",      evt_df.slc.cut.obvious_cosmic),
            ("t0",          evt_df.slc.cut.t0),
            ("FV",          InFV_strict(evt_df)),
            ("nu_score",    evt_df.slc.nu_score > CTE.min_nu_score),
            ("track",       evt_df.slc.cut.track),
            ("chi2",        evt_df.slc.cut_var.n_MIP_candidates > 2),
            ("shower",      evt_df.slc.cut.shower),
            ("angle",       evt_df.slc.cut.angle),
            ("proton_BDT",  evt_df.slc.cut.proton_BDT_2pi),
            ("containment", evt_df.slc.cut_var.n_exiting_pfps == 0),
            ("TPC_containment",  TPC_containment_mask(evt_df,group_levels)),
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