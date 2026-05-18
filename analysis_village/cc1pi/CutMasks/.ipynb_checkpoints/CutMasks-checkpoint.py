import sys
import numpy as np
import math
import uproot as uproot
import pickle
import pandas as pd


from analysis_village.cc1pi.Constants import CTE

# Absolute path to cafpyana directory
print('Importing cafpyana utils...')
cafpyana_root = "/home/lpelegri/cafpyana"
# Add CAFpyana to the Python search path -- this will allow you to import cafpyana modules
sys.path.insert(0, cafpyana_root)

def is_obvious_cosmic_cut_mask(df):
    return df.slc.is_clear_cosmic == 0
    
def t0_cut_mask(df):
    return (df.slc.barycenterFM.flashTime > -3) & (df.slc.barycenterFM.flashTime < 5) & (df.slc.barycenterFM.score > CTE.min_bc_score)

def is_inside_FV_cut_mask(df):
    xmin = -200 + CTE.min_distance_to_wall_x_y
    xmax = 200 - CTE.min_distance_to_wall_x_y
    ymin = -200 + CTE.min_distance_to_wall_x_y
    ymax = 200 - CTE.min_distance_to_wall_x_y
    zmin = CTE.min_distance_to_first_z_wall
    zmax = 500 - CTE.min_distance_to_last_z_wall

    pass_fv = (df.slc.vertex.x > xmin) & (df.slc.vertex.x < xmax) & (df.slc.vertex.y > ymin) & (df.slc.vertex.y < ymax) & (df.slc.vertex.z < zmax) & (df.slc.vertex.z > zmin)
    return pass_fv
    

    
def nu_score_cut_mask(df):
    return df.slc.nu_score > CTE.min_nu_score

def is_pandora_primary_mask(df):
    is_primary_mask = (df.pfp.parent_is_primary == True) 
    return is_primary_mask

def is_analysis_primary_mask(df):
    is_primary_mask = (df.pfp.parent_is_primary == True) & (df.pfp.dist_to_vertex < CTE.max_primary_distance_to_vertex)
    return is_primary_mask
    
def is_primary_track_mask(df):
    is_primary_mask = (df.pfp.parent_is_primary == True) & (df.pfp.dist_to_vertex < CTE.max_primary_distance_to_vertex)
    is_track_mask = (df.pfp.trk.len > CTE.min_track_lenght) & (df.pfp.trackScore > CTE.min_track_score)
    return is_primary_mask & is_track_mask

    
def is_primary_shower_mask(df):
    is_pandora_primary_mask = (df.pfp.parent_is_primary == True)
    is_shower_mask = (df.pfp.trackScore >= 0) & (df.pfp.trackScore < CTE.max_shower_track_score)
    return is_pandora_primary_mask & is_shower_mask


def track_cut_mask(df, group_levels):
    track_df = df[is_primary_track_mask(df)]

    # Count how many pfps per slice
    track_counts = track_df.groupby(level=group_levels).size()
 
    # Get only slices with at least 2 pfps
    valid_slices = track_counts[track_counts > 1].index

    # Apply the mask to original DataFrame
    final_mask = pd.Series(df.index.droplevel('rec.slc.reco.pfp..index').isin(valid_slices), index=df.index)

    return final_mask

def shower_cut_mask(df, group_levels):
    is_pandora_primary_mask = (df.pfp.parent_is_primary == True)
    is_shower_mask = (df.pfp.trackScore >= 0) & (df.pfp.trackScore < CTE.max_shower_track_score)
    energy_mask = (df.pfp.shw.bestplane_energy > CTE.min_shower_ke)
    shower_df = df[is_pandora_primary_mask & is_shower_mask & energy_mask] 
        
    # Count how many pfps per slice
    shower_counts = shower_df.groupby(level=group_levels).size()
    
    # Get only slices with at least 2 pfps
    non_valid_slices = shower_counts[(shower_counts > 0)].index

    # Apply the mask to original DataFrame
    final_mask = pd.Series(~df.index.droplevel('rec.slc.reco.pfp..index').isin(non_valid_slices), index=df.index)
   
    return final_mask

def is_MIP_candidate_mask(df):
    is_primary_mask = (df.pfp.parent_is_primary == True) & (df.pfp.dist_to_vertex < CTE.max_primary_distance_to_vertex)
    is_track_mask = (df.pfp.trk.len > CTE.min_track_lenght) & (df.pfp.trackScore > CTE.min_track_score)
    track_mask = is_primary_mask & is_track_mask

    len_mask = df.pfp.trk.len > CTE.MIP_candidate_min_TL
    chi2_mask = (df.pfp.trk.chi2pid.best.chi2_muon < CTE.MIP_candidate_max_muon_score) & (df.pfp.trk.chi2pid.best.chi2_proton > CTE.MIP_candidate_min_proton_score)
    return track_mask & chi2_mask & len_mask

def chi2_cut_mask(df, group_levels):
    chi2_df = df[is_MIP_candidate_mask(df)]
    
    # Count how many pfps per slice
    MIP_counts = chi2_df.groupby(level=group_levels).size()
 
    # Get only slices with at least 2 pfps
    valid_slices = MIP_counts[MIP_counts == 2].index

    # Apply the mask to original DataFrame
    final_mask = pd.Series(df.index.droplevel('rec.slc.reco.pfp..index').isin(valid_slices), index=df.index)

    return final_mask

def angle_cut_mask(df):
    return df.slc.measure_var.max_angle_between_candidates < CTE.max_angle_between_candidates

def exiting_pfp_mask(df):
    xmin = -200 + CTE.min_distance_to_consider_contained
    xmax = 200 - CTE.min_distance_to_consider_contained
    ymin = -200 + CTE.min_distance_to_consider_contained
    ymax = 200 - CTE.min_distance_to_consider_contained
    zmin = CTE.min_distance_to_consider_contained
    zmax = 500 - CTE.min_distance_to_consider_contained
    
    not_in_fv_start = (df.pfp.trk.start.x < xmin) | (df.pfp.trk.start.x > xmax) | (df.pfp.trk.start.y < ymin) | (df.pfp.trk.start.y > ymax) | (df.pfp.trk.start.z < zmin) | (df.pfp.trk.start.z  > zmax)
    not_in_fv_end = (df.pfp.trk.end.x < xmin) | (df.pfp.trk.end.x > xmax) | (df.pfp.trk.end.y < ymin) | (df.pfp.trk.end.y > ymax) | (df.pfp.trk.end.z < zmin) | (df.pfp.trk.end.z  > zmax)
    return not_in_fv_start | not_in_fv_end

def exiting_z_pfp_mask(df):
    return (df.pfp.trk.end.z < CTE.min_z_for_escaping_p) | (df.pfp.trk.start.z < CTE.min_z_for_escaping_p)

def containment_z_cut_mask(df, group_levels):
    exiting_df = df[exiting_pfp_mask(df)]
    exiting_z_df = df[exiting_z_pfp_mask(df)]

    exiting_counts = exiting_df.groupby(level=group_levels).size()
    exiting_z_counts = exiting_z_df.groupby(level=group_levels).size()
    
    exiting_invalid_slices = exiting_counts[exiting_counts > 1].index
    exiting_z_invalid_slices = exiting_z_counts[exiting_z_counts > 0].index
    
    exiting_mask = pd.Series(~df.index.droplevel('rec.slc.reco.pfp..index').isin(exiting_invalid_slices), index=df.index)
    exiting_z_mask = pd.Series(~df.index.droplevel('rec.slc.reco.pfp..index').isin(exiting_z_invalid_slices), index=df.index)
 
    return exiting_mask & exiting_z_mask


def containment_cut_mask(df, group_levels):
    exiting_df = df[exiting_pfp_mask(df)]
    exiting_counts = exiting_df.groupby(level=group_levels).size()
    
    exiting_invalid_slices = exiting_counts[exiting_counts > 0].index
 
    exiting_mask = pd.Series(~df.index.droplevel('rec.slc.reco.pfp..index').isin(exiting_invalid_slices), index=df.index)
   
    return exiting_mask

 
def michel_cut_mask_old(df, group_levels, max_mean_dEdx = CTE.michel_max_dEdx_mean, max_TS = CTE.michel_max_track_score, max_ke = CTE.michel_max_visible_energy):

    michel_df = df[(is_MIP_candidate_mask(df)) 
        & (df.pfp.trk.mean_dEdx < max_mean_dEdx)
        & (df.pfp.trackScore < max_TS)
        & (df.pfp.trk.calo.best.ke < max_ke)
        & (df.pfp.max_daughter_hits == 0)]

    # Count how many pfps per slice
    candidate_counts = michel_df.groupby(level=group_levels).size()
 
    # Get only slices with at least 2 pfps
    invalid_slices = candidate_counts[candidate_counts > 0].index

    # Apply the mask to original DataFrame
    final_mask = pd.Series(~df.index.droplevel('rec.slc.reco.pfp..index').isin(invalid_slices), index=df.index)

    return final_mask


def michel_cut_mask(df, group_levels):

    michel_df = df[(is_MIP_candidate_mask(df)) 
        & (df.pfp.trackScore < CTE.michel_max_track_score)
        & (df.pfp.trk.calo.best.ke < CTE.michel_max_visible_energy)
        & (df.pfp.max_daughter_hits == 0)]

    # Count how many pfps per slice
    candidate_counts = michel_df.groupby(level=group_levels).size()
 
    # Get only slices with at least 2 pfps
    invalid_slices = candidate_counts[candidate_counts > 0].index

    # Apply the mask to original DataFrame
    final_mask = pd.Series(~df.index.droplevel('rec.slc.reco.pfp..index').isin(invalid_slices), index=df.index)

    return final_mask
def extra_pion_cut_mask(df, group_levels):
    chi2_mask = (df.pfp.trk.chi2pid.best.chi2_muon < CTE.MIP_candidate_max_muon_score) & (df.pfp.trk.chi2pid.best.chi2_proton > CTE.MIP_candidate_min_proton_score)
    relaxed_MIP_df = df[(is_primary_track_mask(df)) & (chi2_mask)]
    
    # Count how many pfps per slice
    candidate_counts = relaxed_MIP_df.groupby(level=group_levels).size()
 
    # Get only slices with at least 2 pfps
    valid_slices = candidate_counts[candidate_counts == 2].index

    # Apply the mask to original DataFrame
    final_mask = pd.Series(df.index.droplevel('rec.slc.reco.pfp..index').isin(valid_slices), index=df.index)

    return final_mask
    
# Create a column filled with NaN first
def proton_BDT_cut_mask(df, group_levels):
    BDT_proton_df = df[(is_MIP_candidate_mask(df)) & (df.pfp.trk.bdt_proton_score < CTE.BDT_proton_max_score)]
    
    # Count how many pfps per slice
    candidate_counts = BDT_proton_df.groupby(level=group_levels).size()
 
    # Get only slices with at least 2 pfps
    invalid_slices = candidate_counts[candidate_counts > 0].index

    # Apply the mask to original DataFrame
    final_mask = pd.Series(~df.index.droplevel('rec.slc.reco.pfp..index').isin(invalid_slices), index=df.index)

    return final_mask
    
def proton_BDT_sideband_mask(df, group_levels):
    BDT_proton_df = df[(is_MIP_candidate_mask(df)) & (df.pfp.trk.bdt_proton_score < CTE.BDT_proton_max_score_sideband)]
    
    # Count how many pfps per slice
    candidate_counts = BDT_proton_df.groupby(level=group_levels).size()
 
    # Get only slices with at least 2 pfps
    valid_slices = candidate_counts[candidate_counts == 1].index

    # Apply the mask to original DataFrame
    final_mask = pd.Series(df.index.droplevel('rec.slc.reco.pfp..index').isin(valid_slices), index=df.index)

    return final_mask

def proton_BDT_cut_mask_2pi(df, group_levels):
    BDT_proton_df = df[(is_MIP_candidate_mask(df)) & (df.pfp.trk.bdt_proton_score < CTE.BDT_proton_max_score_sideband_pion)]
    
    # Count how many pfps per sliceVenga, mucho ani
    candidate_counts = BDT_proton_df.groupby(level=group_levels).size()
 
    # Get only slices with at least 2 pfps
    invalid_slices = candidate_counts[candidate_counts > 0].index

    # Apply the mask to original DataFrame
    final_mask = pd.Series(~df.index.droplevel('rec.slc.reco.pfp..index').isin(invalid_slices), index=df.index)

    return final_mask

def InFV_strict(df):
    # Access the series directly to avoid creating extra variable references
    x = df.slc.vertex.x
    y = df.slc.vertex.y
    z = df.slc.vertex.z

    # Combine everything into a single evaluation path. 
    # Python will evaluate this much faster and clear out temporary arrays.
    contained = (np.abs(x) > 5) & (np.abs(x) < 190) & (
        ((z > 10)  & (z < 250) & (np.abs(y) < 190)) |
        ((z > 250) & (z < 450) & (y > -190) & (y < 100) & (x < 0)) |
        ((z > 250) & (z < 450) & (y > -190) & (y < 190) & (x > 0))
    )
    
    return contained


def cathode_crossing_pfp_mask(df):
    xmin = -CTE.min_distance_to_consider_contained
    xmax = CTE.min_distance_to_consider_contained
    
    crossing_cathode = ((df.pfp.trk.start.x > xmin) & (df.pfp.trk.start.x < xmax))|((df.pfp.trk.end.x > xmin) & (df.pfp.trk.end.x < xmax) )
    return crossing_cathode

def ends_in_high_y_high_z(df):
    x = df.pfp.trk.end.x
    y = df.pfp.trk.end.y
    z = df.pfp.trk.end.z
    in_high_y_high_z = (z > 250) & (y > 100) & (x < 0)
    
    return in_high_y_high_z


def not_in_high_y_high_z_containment_mask(df, group_levels):
    high_yz_df = df[ends_in_high_y_high_z(df)]
    
    # Count how many pfps per slice
    candidate_counts = high_yz_df.groupby(level=group_levels).size()
 
    # Get only slices with at least 2 pfps
    invalid_slices = candidate_counts[candidate_counts > 0].index

    # Apply the mask to original DataFrame
    final_mask = pd.Series(~df.index.droplevel('rec.slc.reco.pfp..index').isin(invalid_slices), index=df.index)

    return final_mask
    
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