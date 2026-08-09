from makedf.makedf import *
from pyanalib.pandas_helpers import *
from makedf.util import *
from makedf.geniesyst import *
from tqdm import tqdm
from analysis_village.cc1pi.Constants import CTE
from analysis_village.cc1pi.CutMasks import CutMasks
from analysis_village.cc1pi.dEdxCleaning import dEdxCleaning
from analysis_village.cc1pi.makedf.make_cc1pidf import *
import lmfit
import uproot as uproot
from ROOT import TMVA
import os
import ROOT
import subprocess
from ROOT import std
import pickle
from itertools import combinations

cols = [
        #Pfp truth info
        ('pfp', 'trk', 'truth', 'p', 'p_type', ''),
        ('pfp', 'trk', 'truth', 'p', 'pdg', ''),
        ('pfp', 'trk', 'truth', 'p', 'end_process', ''),
        ('pfp', 'trk', 'truth', 'genp', 'mag', ''),
        ('pfp', 'trk', 'rangeP', 'p_pion', '', ''),
        ('pfp', 'trk', 'rangeP', 'p_muon', '', ''),
        ('pfp', 'trk', 'chi2_exp_pol', '', '', ''),
        ('pfp', 'trk', 'len', '', '', ''),
        ('pfp', 'trk', 'chi2pid', 'best', 'chi2_pion', ''),
        ('pfp', 'trk', 'chi2pid', 'best', 'chi2_proton', ''),
        ('pfp', 'trk', 'start', 'x', '', ''),
        ('pfp', 'trk', 'start', 'y', '', ''),
        ('pfp', 'trk', 'start', 'z', '', ''),
        ('pfp', 'trk', 'end', 'x', '', ''),
        ('pfp', 'trk', 'end', 'y', '', ''),
        ('pfp', 'trk', 'end', 'z', '', ''),
        ('pfp', 'nhit0', '', '', '', ''),
        ('pfp', 'nhit1', '', '', '', ''),
        ('pfp', 'nhit2', '', '', '', ''),
    
    ]

def get_pion_candidate_id(group):
    """
    Applies the same muon/pion separation logic as get_mu_pi_vars,
    but only returns the ('pfp', 'id') of the selected pion candidate.
    Returns -999 if no valid pion candidate is found.
    """
    INVALID = pd.Series({'pion_candidate_id': -999})

    mip_df = group[CutMasks.is_MIP_candidate_mask(group)]
    
    # Safety check
    if len(mip_df) != 2:
        return INVALID

    bdt_score_col = ('pfp', 'trk', 'bdt_muon_pion_score', '', '', '')
    len_col       = ('pfp', 'trk', 'len', '', '', '')

    exiting_mask = CutMasks.exiting_pfp_mask(mip_df)

    # ✅ FIXED: Explicitly evaluate as a scalar boolean using .any() or int comparison
    if exiting_mask.any():
        muon_row = mip_df.loc[exiting_mask].iloc[[0]]
    else:
        muon_row = mip_df.sort_values(bdt_score_col).iloc[[-1]]

    remaining_mips = mip_df.drop(muon_row.index)
    pion_row = remaining_mips.sort_values(len_col).iloc[[-1]]

    # Extract ID safely
    pion_id = pion_row[('pfp', 'id', '', '', '', '')].iloc[0]
    return pd.Series({'pion_candidate_id': pion_id})


def make_pion_selection_df(f, updatecalo = None, select_stopping = False):
    
    pandora_df = make_pandora_df(f, trkScoreCut = False, trkDistCut= -1, cutClearCosmic = True, requireFiducial=False, updatecalo=updatecalo)
    if pandora_df.empty:
       idx = pd.MultiIndex(
           levels=[[], [], []],
           codes=[[], [], []],
           names=['entry', 'rec.slc..index', 'rec.slc.reco.pfp..index']
       )
       empty_df = pd.DataFrame(index=idx, columns=pd.MultiIndex.from_tuples(cols))
       return empty_df
        
    cc1pi_shwbranches = [shwbranch + 'bestplane_energy']
    shw_df = loadbranches(f["recTree"], cc1pi_shwbranches)
    shw_df = shw_df.rec.slc.reco

    # create mcs df
    mcs_df = loadbranches(f["recTree"], [trkmcsbranches[1]]).rec.slc.reco.pfp.trk.mcsP

    hit0_df = make_trkhitdf_plane0(f, updatecalo = updatecalo).sort_values('rr', ascending=True)
    hit1_df = make_trkhitdf_plane1(f, updatecalo = updatecalo).sort_values('rr', ascending=True)
    hit2_df = make_trkhitdf_plane2(f, updatecalo = updatecalo).sort_values('rr', ascending=True)
        
    hit_dfs = [hit0_df, hit1_df, hit2_df]
    hit_names = ['nhit0', 'nhit1', 'nhit2']  

    best_hit_df = get_best_hit_df(pandora_df, hit_dfs)
    pandora_df = add_nhit_column(pandora_df, hit_dfs, hit_names, 1000)
    pandora_df = multicol_merge(pandora_df, shw_df, left_index=True, right_index=True, how="left", validate="one_to_one")  
    
    pandora_df = pandora_df[CutMasks.is_obvious_cosmic_cut_mask(pandora_df)]
    pandora_df[('slc', 'vertex_inside_fv', '', '', '', '')] = CutMasks.is_inside_FV_cut_mask(pandora_df)
                 
    pandora_df = add_best_ke_column(pandora_df)
    pandora_df = add_p_type_column(pandora_df)  
    pandora_df = add_best_chi2_columns(pandora_df, update_calo = updatecalo)
    pandora_df = add_max_daughter_hits_column(pandora_df)

    # Remove non-primary particles
    pandora_df = pandora_df[pandora_df.pfp.parent_is_primary == True]
    
    fixed_hit_df = (
        best_hit_df
        .groupby(level=group_levels, group_keys=False)
        .apply(dEdxCleaning.get_fix_hit_df, plot = False)
    )
    fixed_hit_df['dE'] = fixed_hit_df.dedx * fixed_hit_df.pitch
    
    pandora_df = add_mean_dedx_column(pandora_df, fixed_hit_df)  
    pandora_df = add_chi2_exp_pol_column(pandora_df, fixed_hit_df)
    pandora_df = add_frac50_column(pandora_df, fixed_hit_df)
    pandora_df = add_scatter_angle_ratio_column(pandora_df, mcs_df)
    pandora_df = add_max_angle_between_candidates_column(pandora_df)
    pandora_df[('slc', 'measure_var', 'angle_between_candidates', '', '', '')] = pandora_df[('slc', 'measure_var', 'max_angle_between_candidates', '', '', '')]
     
    pandora_df[('pfp', 'is_exiting', '', '', '', '')] = CutMasks.exiting_pfp_mask(pandora_df)
    pandora_df[('pfp', 'is_exiting_z', '', '', '', '')] = CutMasks.exiting_z_pfp_mask(pandora_df)
    
    # Add BDT columns
    proton_BDT_input_columns = [
        ('pfp','trk','chi2pid','best','chi2_muon',''),
        ('pfp','trk','chi2pid','best','chi2_proton',''),
        ('pfp','trk','chi2_exp_pol','','',''),
        ('pfp','trk','frac50','','','')
    ]
    col_BDTG_score_proton = ('pfp','trk','bdt_proton_score','','','')
    pandora_df = add_bdt_score(
        pandora_df,
        base_running_path + "/../BDTs/bdtg_model_proton.pkl",
        proton_BDT_input_columns,
        col_BDTG_score_proton
    )

    muon_pion_BDT_input_columns = [
        ('pfp', 'trk', 'chi2pid', 'best', 'chi2_muon', ''),
        ('pfp', 'trk', 'chi2pid', 'best', 'chi2_proton', ''),
        ('pfp', 'trk', 'chi2_exp_pol', '', '', ''),
        ('pfp', 'scatter_angle_ratio', '', '', '', ''),
        ('pfp', 'max_daughter_hits', '', '', '', '')
    ]
    col_BDTG_score_muon_pion = ('pfp','trk','bdt_muon_pion_score','','','')
    pandora_df = add_bdt_score(
        pandora_df,
        base_running_path + "/../BDTs/bdtg_model_muon_pion.pkl",
        muon_pion_BDT_input_columns,
        col_BDTG_score_muon_pion
    )

    # CC1Pi selection cuts
    pandora_df = pandora_df[CutMasks.is_obvious_cosmic_cut_mask(pandora_df)]
    pandora_df = pandora_df[CutMasks.is_inside_FV_cut_mask(pandora_df)]
    pandora_df = pandora_df[CutMasks.t0_cut_mask(pandora_df)]
    pandora_df = pandora_df[CutMasks.nu_score_cut_mask(pandora_df)]
    pandora_df = pandora_df[CutMasks.track_cut_mask(pandora_df, group_levels)]
    pandora_df = pandora_df[CutMasks.shower_cut_mask(pandora_df, group_levels)]
    pandora_df = pandora_df[CutMasks.chi2_cut_mask(pandora_df, group_levels)]
    pandora_df = pandora_df[CutMasks.angle_cut_mask(pandora_df)]
    pandora_df = pandora_df[CutMasks.proton_BDT_cut_mask(pandora_df, group_levels)]
    pandora_df = pandora_df[CutMasks.containment_z_cut_mask(pandora_df, group_levels)]
    pandora_df = pandora_df[CutMasks.michel_cut_mask(pandora_df, group_levels)]
    pandora_df = pandora_df[CutMasks.extra_pion_cut_mask(pandora_df, group_levels)]

    # Select MIP candidates
    pandora_df = pandora_df[CutMasks.is_MIP_candidate_mask(pandora_df)]
    if pandora_df.empty:
        return pd.DataFrame(index=pandora_df.index, columns=pd.MultiIndex.from_tuples(cols))


    
    #Add true p 
    genp_x = pandora_df[('pfp', 'trk', 'truth', 'p', 'genp', 'x')]  # check 'p' level naming
    genp_y = pandora_df[('pfp', 'trk', 'truth', 'p', 'genp', 'y')]
    genp_z = pandora_df[('pfp', 'trk', 'truth', 'p', 'genp', 'z')]
    mag = np.sqrt(genp_x**2 + genp_y**2 + genp_z**2)
    pandora_df.loc[:, ('pfp', 'trk', 'truth', 'genp', 'mag', '')] = mag

    
    # --- Apply pion candidate selection per slice ---
    pion_candidates = pandora_df.groupby(
        level=['entry', 'rec.slc..index'], group_keys=False
    ).apply(get_pion_candidate_id)

    # Convert output series/dataframe into a predictable single-level frame for filtering
    if isinstance(pion_candidates, pd.Series):
        pion_candidates = pion_candidates.to_frame()
    if isinstance(pion_candidates.columns, pd.MultiIndex):
        pion_candidates.columns = pion_candidates.columns.get_level_values(0)

    # Filter out invalid entries (-999)
    pion_candidates = pion_candidates[pion_candidates['pion_candidate_id'] != -999]
    if pion_candidates.empty:
        return pd.DataFrame(index=pandora_df.index, columns=pd.MultiIndex.from_tuples(cols))

    # --- Index-based filtering instead of merge to avoid MultiIndex column corruption ---
    # Match the pfp id column directly against the candidate map
    pfp_ids = pandora_df[('pfp', 'id', '', '', '', '')]
    slice_indices = pandora_df.index.droplevel('rec.slc.reco.pfp..index')
    
    # Map the chosen pion_candidate_id back onto pandora_df
    target_pion_ids = slice_indices.map(pion_candidates['pion_candidate_id'])
    
    # Filter pandora_df directly
    pandora_pion_df = pandora_df[pfp_ids == target_pion_ids].copy()
    
    # Veto high-y/high-z region and cathode region
    pandora_pion_df = pandora_pion_df[~CutMasks.ends_in_high_y_high_z(pandora_pion_df)]
    pandora_pion_df = pandora_pion_df[~CutMasks.starts_in_high_y_high_z(pandora_pion_df)]
    pandora_pion_df = pandora_pion_df[~CutMasks.exiting_pfp_mask(pandora_pion_df)]
    
    x_start = pandora_pion_df[('pfp', 'trk', 'start', 'x', '', '')].abs()
    x_end   = pandora_pion_df[('pfp', 'trk', 'end', 'x', '', '')].abs()
    
    pandora_pion_df = pandora_pion_df[(x_start > 10) & (x_end > 10)]


    # Select stopping candidates if requested
    if select_stopping:
        pandora_pion_df = pandora_pion_df[pandora_pion_df.pfp.trk.chi2pid.best.chi2_pion < 6]
        pandora_pion_df = pandora_pion_df[pandora_pion_df.pfp.trk.chi2_exp_pol > 0.25]
        pandora_pion_df = pandora_pion_df[pandora_pion_df.pfp.trk.chi2_exp_pol < 0.75]
        pandora_pion_df = pandora_pion_df[pandora_pion_df.pfp.trk.chi2pid.best.chi2_proton < 130]
        pandora_pion_df = pandora_pion_df[pandora_pion_df.pfp.trk.len < 100]
    
    min_df = pandora_pion_df[cols].copy()
    min_df = min_df[min_df[('pfp', 'trk', 'len', '', '', '')] > 0]
    return min_df



def get_muon_candidate_id(group):
    """
    Finds MIP candidates in the group and returns the ('pfp', 'id') 
    of the track with the maximum length as the muon candidate.
    Returns -999 if no valid MIP candidates are present.
    """
    INVALID = pd.Series({'muon_candidate_id': -999})

    mip_df = group[CutMasks.is_MIP_candidate_mask(group)]

    # Safety check: ensure at least one MIP candidate exists
    if mip_df.empty:
        return INVALID

    len_col = ('pfp', 'trk', 'len', '', '', '')
    id_col  = ('pfp', 'id', '', '', '', '')

    # Sort by track length ascending and select the last row (longest track)
    muon_row = mip_df.sort_values(len_col).iloc[[-1]]

    # Extract ID safely
    muon_id = muon_row[id_col].iloc[0]
    
    return pd.Series({'muon_candidate_id': muon_id})
    

def make_muon_selection_df(f, updatecalo = None):
    
    pandora_df = make_pandora_df(f, trkScoreCut = False, trkDistCut= -1, cutClearCosmic = True, requireFiducial=False, updatecalo=updatecalo)
    if pandora_df.empty:
       idx = pd.MultiIndex(
           levels=[[], [], []],
           codes=[[], [], []],
           names=['entry', 'rec.slc..index', 'rec.slc.reco.pfp..index']
       )
       empty_df = pd.DataFrame(index=idx, columns=pd.MultiIndex.from_tuples(cols))
       return empty_df
        
    hit0_df = make_trkhitdf_plane0(f, updatecalo = updatecalo).sort_values('rr', ascending=True)
    hit1_df = make_trkhitdf_plane1(f, updatecalo = updatecalo).sort_values('rr', ascending=True)
    hit2_df = make_trkhitdf_plane2(f, updatecalo = updatecalo).sort_values('rr', ascending=True)
        
    hit_dfs = [hit0_df, hit1_df, hit2_df]
    hit_names = ['nhit0', 'nhit1', 'nhit2']  
    pandora_df = add_nhit_column(pandora_df, hit_dfs, hit_names, 1000)

    pandora_df = add_p_type_column(pandora_df)
    pandora_df = add_best_chi2_columns(pandora_df, update_calo = updatecalo)
    
    best_hit_df = get_best_hit_df(pandora_df, hit_dfs)
    fixed_hit_df = (
        best_hit_df
        .groupby(level=group_levels, group_keys=False)
        .apply(dEdxCleaning.get_fix_hit_df, plot = False)
    )
    fixed_hit_df['dE'] = fixed_hit_df.dedx * fixed_hit_df.pitch
    
    pandora_df = add_chi2_exp_pol_column(pandora_df, fixed_hit_df)
    
    pandora_df = add_n_exiting_pfps_column(pandora_df)
    
    pandora_df = pandora_df[CutMasks.is_obvious_cosmic_cut_mask(pandora_df)]
    pandora_df = pandora_df[CutMasks.is_inside_FV_cut_mask(pandora_df)]
    pandora_df = pandora_df[CutMasks.t0_cut_mask(pandora_df)]
    pandora_df = pandora_df[CutMasks.nu_score_cut_mask(pandora_df)]
    pandora_df = pandora_df[CutMasks.track_cut_mask(pandora_df, group_levels)]
    pandora_df = pandora_df[(pandora_df.slc.cut_var.n_exiting_pfps == 0)]
    
    # Select MIP candidates
    pandora_df = pandora_df[CutMasks.is_MIP_candidate_mask(pandora_df)]
    if pandora_df.empty:
        return pd.DataFrame(index=pandora_df.index, columns=pd.MultiIndex.from_tuples(cols))

    #Add true p 
    genp_x = pandora_df[('pfp', 'trk', 'truth', 'p', 'genp', 'x')]  # check 'p' level naming
    genp_y = pandora_df[('pfp', 'trk', 'truth', 'p', 'genp', 'y')]
    genp_z = pandora_df[('pfp', 'trk', 'truth', 'p', 'genp', 'z')]
    mag = np.sqrt(genp_x**2 + genp_y**2 + genp_z**2)
    pandora_df.loc[:, ('pfp', 'trk', 'truth', 'genp', 'mag', '')] = mag

    
    # --- Apply muon candidate selection per slice --- (Just take the longest pfp)
    muon_candidates = pandora_df.groupby(
        level=['entry', 'rec.slc..index'], group_keys=False
    ).apply(get_muon_candidate_id)

    # Convert output series/dataframe into a predictable single-level frame for filtering
    if isinstance(muon_candidates, pd.Series):
        muon_candidates = muon_candidates.to_frame()
    if isinstance(muon_candidates.columns, pd.MultiIndex):
        muon_candidates.columns = muon_candidates.columns.get_level_values(0)

    # Filter out invalid entries (-999)
    muon_candidates = muon_candidates[muon_candidates['muon_candidate_id'] != -999]
    if muon_candidates.empty:
        return pd.DataFrame(index=pandora_df.index, columns=pd.MultiIndex.from_tuples(cols))

    # --- Index-based filtering instead of merge to avoid MultiIndex column corruption ---
    # Match the pfp id column directly against the candidate map
    pfp_ids = pandora_df[('pfp', 'id', '', '', '', '')]
    slice_indices = pandora_df.index.droplevel('rec.slc.reco.pfp..index')
    
    # Map the chosen muon_candidate_id back onto pandora_df
    target_muon_ids = slice_indices.map(muon_candidates['muon_candidate_id'])
    
    # Filter pandora_df directly
    pandora_muon_df = pandora_df[pfp_ids == target_muon_ids].copy()
    
    # Veto high-y/high-z region and cathode region
    pandora_muon_df = pandora_muon_df[~CutMasks.ends_in_high_y_high_z(pandora_muon_df)]
    pandora_muon_df = pandora_muon_df[~CutMasks.starts_in_high_y_high_z(pandora_muon_df)]
    pandora_muon_df = pandora_muon_df[~CutMasks.exiting_pfp_mask(pandora_muon_df)]
    x_start = pandora_muon_df[('pfp', 'trk', 'start', 'x', '', '')].abs()
    x_end   = pandora_muon_df[('pfp', 'trk', 'end', 'x', '', '')].abs()
    pandora_muon_df = pandora_muon_df[(x_start > 10) & (x_end > 10)]

    # Extra conditions for making sure its a muon
    pandora_muon_df = pandora_muon_df[pandora_muon_df.pfp.trk.chi2pid.best.chi2_muon < 6]
    pandora_muon_df = pandora_muon_df[pandora_muon_df.pfp.trk.len > 50]

    min_df = pandora_muon_df[cols].copy()
    min_df = min_df[min_df[('pfp', 'trk', 'len', '', '', '')] > 0]
    return min_df


def make_trkhitdf_selection_df(f, plane=0, pdg=211, select_stopping=False, updatecalo=None):
    hit_df = make_trkhitdf(f, plane, updatecalo=updatecalo)

    if pdg == 211:
        selection_df = make_pion_selection_df(
            f, updatecalo=updatecalo, select_stopping=select_stopping
        )
    elif pdg == 13:
        selection_df = make_muon_selection_df(
            f, updatecalo=updatecalo
        )
    else:
        print("ERROR PDG NOT IN SELECTION TYPES")
        # Return an empty DataFrame with hit_df's index structure if available
        return hit_df.iloc[0:0] if 'hit_df' in locals() and not hit_df.empty else pd.DataFrame()

    # --- Early return with preserved MultiIndex structure ---
    if selection_df.empty or hit_df.empty:
        # iloc[0:0] keeps 0 rows while preserving full MultiIndex & column dtypes
        return hit_df.iloc[0:0]

    # Shared levels between hit_df and selection_df
    common_levels = ['entry', 'rec.slc..index', 'rec.slc.reco.pfp..index']

    # Safely align level ordering by extracting common columns explicitly
    hit_common_tuples = set(
        hit_df.index.to_frame()[common_levels].itertuples(index=False, name=None)
    )
    sel_common_tuples = set(
        selection_df.index.to_frame()[common_levels].itertuples(index=False, name=None)
    )

    # Intersection of shared tuples
    valid_tuples = hit_common_tuples.intersection(sel_common_tuples)

    if not valid_tuples:
        return hit_df.iloc[0:0]

    # Build mask ensuring level order matches common_levels
    hit_tuples_list = list(
        hit_df.index.to_frame()[common_levels].itertuples(index=False, name=None)
    )
    mask = [tup in valid_tuples for tup in hit_tuples_list]

    filtered_hit_df = hit_df[mask]
    
    return filtered_hit_df




def make_pion_selection_all_df(f):
    return make_pion_selection_df(f, updatecalo = None)

def make_trkhitdf_plane0_pion_selection_all(f):
    return make_trkhitdf_selection_df(f, plane=0, pdg=211, select_stopping=False, updatecalo=None)
    
def make_trkhitdf_plane1_pion_selection_all(f):
    return make_trkhitdf_selection_df(f, plane = 1, pdg = 211, select_stopping = False, updatecalo=None)
    
def make_trkhitdf_plane2_pion_selection_all(f):
    return make_trkhitdf_selection_df(f, plane = 2, pdg = 211, select_stopping = False, updatecalo=None)


    
def make_pion_selection_stopping_df(f):
    return make_pion_selection_df(f, updatecalo = None, select_stopping = True)

def make_trkhitdf_plane0_pion_selection_stopping(f):
    return make_trkhitdf_selection_df(f, plane = 0, pdg = 211, select_stopping = True, updatecalo=None)

def make_trkhitdf_plane1_pion_selection_stopping(f):
    return make_trkhitdf_selection_df(f, plane = 1, pdg = 211, select_stopping = True, updatecalo=None)
    
def make_trkhitdf_plane2_pion_selection_stopping(f):
    return make_trkhitdf_selection_df(f, plane = 2, pdg = 211, select_stopping = True, updatecalo=None)


    

def make_pion_selection_stopping_update_calo_cv_df(f):
    return make_pion_selection_df(f, updatecalo = "cv", select_stopping = True)

def make_trkhitdf_plane0_pion_selection_stopping_update_calo_cv(f):
    return make_trkhitdf_selection_df(f, plane = 0, pdg = 211, select_stopping = True, updatecalo="cv")

def make_trkhitdf_plane1_pion_selection_stopping_update_calo_cv(f):
    return make_trkhitdf_selection_df(f, plane = 1, pdg = 211, select_stopping = True, updatecalo="cv")
    
def make_trkhitdf_plane2_pion_selection_stopping_update_calo_cv(f):
    return make_trkhitdf_selection_df(f, plane = 2, pdg = 211, select_stopping = True, updatecalo="cv")






    

def make_muon_selection_no_calo_df(f):
    return make_muon_selection_df(f, updatecalo = None)

def make_trkhitdf_plane0_muon_selection(f):
    return make_trkhitdf_selection_df(f, plane=0, pdg=13, select_stopping=False, updatecalo=None)
    
def make_trkhitdf_plane1_muon_selection(f):
    return make_trkhitdf_selection_df(f, plane = 1, pdg = 13, select_stopping = False, updatecalo=None)
    
def make_trkhitdf_plane2_muon_selection(f):
    return make_trkhitdf_selection_df(f, plane = 2, pdg = 13, select_stopping = False, updatecalo=None)