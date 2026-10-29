from makedf.makedf import *
from pyanalib.pandas_helpers import *
from makedf.util import *
from makedf.geniesyst import *
from tqdm import tqdm
from analysis_village.cc1pi.Constants import CTE
from analysis_village.cc1pi.CutMasks import CutMasks
from analysis_village.cc1pi.dEdxCleaning import dEdxCleaning
import lmfit
import uproot as uproot
#from ROOT import TMVA
import os
import ROOT
import subprocess
from ROOT import std
import pickle
from itertools import combinations


base_running_path = os.path.dirname(os.path.abspath(__file__))
base_path = base_running_path + "/../TLExtensionMethod"
 
try:
    root_lib_dir = subprocess.check_output(['root-config', '--libdir'], text=True).strip()
except:
    root_lib_dir = "/cvmfs/larsoft.opensciencegrid.org/products/root/v6_28_12/Linux64bit+3.10-2.17-e26-p3915-prof/lib"


# 2. Update Environment
os.environ['LD_LIBRARY_PATH'] = f"{root_lib_dir}:{os.environ.get('LD_LIBRARY_PATH', '')}"
ROOT.gSystem.AddDynamicPath(root_lib_dir)

# 3. Load the "Big" ROOT dependency blocks
# Loading these covers almost all physics class requirements (Hist, Geom, Graf, etc.)
root_libs = ["libRint","libROOTGpadv7","libMathMore","libCore", "libRIO", "libNet", "libHist", "libGraf", "libGraf3d", "libGpad", "libTree", "libMathCore", "libThread", "libMatrix", "libGeom","libROOTHist"]
for lib in root_libs:
    ROOT.gSystem.Load(lib)

def load_custom_class(class_name):
    source_file = os.path.join(base_path, f"{class_name}.cpp")
    header_file = os.path.join(base_path, f"{class_name}.h")
    so_file = os.path.join(base_path, f"{class_name}_cpp.so")

    status = ROOT.gSystem.CompileMacro(source_file, "kOf")
    if status >= 0:
        if ROOT.gSystem.Load(so_file) < 0:
            print(f"⚠️  Compiled {class_name} but failed to load {so_file}")
            return False
        ROOT.gInterpreter.Declare(f'#include "{header_file}"')  # keep this!
        print(f"📦 Successfully compiled and linked: {class_name}")
        return True
    else:
        print(f"❌ Failed to compile: {class_name}")
        return False

if load_custom_class("PhysdEdx"):
    if load_custom_class("Hypfit"):
        try:
            ROOT.gInterpreter.ProcessLine(f'#include "{os.path.join(base_path, "Hypfit.h")}"')
            h_fit = ROOT.Hypfit()
            print("🚀 Success! Hypfit object initialized and linked.")
        except Exception as e:
            print(f"❌ Error during instantiation: {e}")
            ROOT.gInterpreter.ProcessLine("Hypfit* h_fit_ptr = new Hypfit();")
            h_fit = ROOT.h_fit_ptr
            print("🚀 Success! Hypfit object initialized via Global Pointer fallback.")
  
# Keep your original calls (3-arg for MC, 4-arg with True for data)
_CONV_TF1_MAPS = {
    (211, False): h_fit.get_conv_function_map(211, "og", 500),
    (211, True):  h_fit.get_conv_function_map(211, "og", 500, True),
    (13,  False): h_fit.get_conv_function_map(13,  "og", 500),
    (13,  True):  h_fit.get_conv_function_map(13,  "og", 500, True),
}
 
 
def select_conv_tf1_map(target_pdg, use_data):
    key = (abs(int(target_pdg)), bool(use_data))
    try:
        return _CONV_TF1_MAPS[key]
    except KeyError:
        valid = sorted({pdg for pdg, _ in _CONV_TF1_MAPS})
        raise ValueError(f"No conv TF1 map for pdg={target_pdg}. Available: {valid}")


 
 


         
# Silence ROOT
#ROOT.gSystem.RedirectOutput(os.devnull, "w")
#ROOT.gROOT.SetBatch(True)

group_levels = ['entry', 'rec.slc..index']
pfp_levels = ['entry','rec.slc..index','rec.slc.reco.pfp..index']



def TruthInFV(data):
    x_region = (np.abs(data.x) > 5) & (np.abs(data.x) < 190)
    z_region1 = (data.z > 10)  & (data.z < 250) & (np.abs(data.y) < 190)
    z_region2 = (data.z > 250) & (data.z < 450) & (data.y > -190) & (data.y < 100) & (data.x < 0)
    z_region3 = (data.z > 250) & (data.z < 450) & (data.y > -190) & (data.y < 190) & (data.x > 0)
    
    contained = x_region & (z_region1 | z_region2 | z_region3)
    return contained

    
def IsNu(df):
    is_numu = abs(df.pdg) == 14
    is_nue = abs(df.pdg) == 12
    return is_numu | is_nue   

def isCC1Pi(df): # definition
    is_1pi1mu = (df.nmu_P_100MeV_3000MeV == 1) & (df.npi_P_130MeV_2000MeV == 1) & (df.npi_P_85MeV_10000MeV == 1)
    is_NpiNmuNnNp = df.nprim - df.nmu - df.npi - df.np - df.nn == 0

    # Initialize full theta mask (False by default)
    is_theta = pd.Series(False, index=df.index)

    # Only compute angles where needed
    df_sel = df.loc[is_1pi1mu]
    
    if len(df_sel) > 0:
        cpi_vec = df_sel.loc[:, ('cpi','genp',['x','y','z'])].to_numpy()
        mu_vec  = df_sel.loc[:, ('mu','genp',['x','y','z'])].to_numpy()
     

        mu_mag  = np.linalg.norm(mu_vec, axis=1)
        cpi_mag = np.linalg.norm(cpi_vec, axis=1)
        dot     = np.sum(mu_vec * cpi_vec, axis=1)

        cos_theta = dot / np.clip(mu_mag * cpi_mag, 1e-12, None)
        theta     = np.arccos(np.clip(cos_theta, -1.0, 1.0))
            
        # Assign back using the SAME index subset
        is_theta.loc[df_sel.index] = theta < CTE.max_angle_between_candidates

    is_mu_p =  df.mu.totp < 1
    #return is_1pi1mu & is_NpiNmuNnNp & is_theta & is_mu_contained
    return is_1pi1mu & is_NpiNmuNnNp & is_theta & is_mu_p

def add_nu_categ_column(df, is_truth_df = False):
    if(is_truth_df):
        truth_df = df
    else:
        truth_df = df.slc.truth # Make a copy to safely assign

    is_inside_fv = TruthInFV(truth_df.position)
    is_nu = IsNu(truth_df)
    is_signal = isCC1Pi(truth_df)
    is_cc = truth_df.iscc.astype(bool)
    is_nu_mu_cc = is_cc & (abs(truth_df.pdg) == 14)

    nu_categ = pd.Series("none", index=truth_df.index, dtype="object")
    # Apply categories
    nu_categ[~is_nu] = "cosmic"
    nu_categ[is_nu & ~is_inside_fv] = "out_AV_nu"
    nu_categ[is_nu & is_inside_fv & ~is_cc] = "NC"
    nu_categ[is_nu & is_inside_fv & is_cc & (abs(truth_df.pdg) == 12)] = "CC_e"
    nu_categ[is_nu & is_inside_fv & is_nu_mu_cc & (truth_df.npi_P_85MeV_10000MeV == 0) & (truth_df.np_P_325MeV_10000MeV == 1)] = "CC_mu_0pi_1p"
    nu_categ[is_nu & is_inside_fv & is_nu_mu_cc & (truth_df.npi_P_85MeV_10000MeV == 0) & (truth_df.np_P_325MeV_10000MeV > 1)] = "CC_mu_0pi_2p"
    nu_categ[is_nu & is_inside_fv & is_nu_mu_cc & (truth_df.npi_P_85MeV_10000MeV == 0) & (truth_df.np_P_325MeV_10000MeV == 0)] = "CC_mu_0pi_0p"
    nu_categ[is_nu & is_inside_fv & is_nu_mu_cc & (truth_df.npi_P_85MeV_10000MeV > 1)] = "CC_mu_2pi"
    nu_categ[is_nu & is_inside_fv & is_nu_mu_cc & is_signal] = "CC1pi"
    nu_categ[is_nu & is_inside_fv & is_nu_mu_cc & ~is_signal & (truth_df.npi_P_85MeV_10000MeV == 1)] = "other_CC1pi"
    
    
    if(is_truth_df):
        df['nu_categ'] = nu_categ
    else:
        df[('slc','truth', 'nu_categ', '', '','')] = nu_categ 
    return df

   
def add_nu_categ_proton_reduced_column(df, is_truth_df = False):
    if(is_truth_df):
        truth_df = df
    else:
        truth_df = df.truth # Make a copy to safely assign

    is_inside_fv = TruthInFV(truth_df.position)
    is_nu = IsNu(truth_df)
    is_signal = isCC1Pi(truth_df)
    is_cc = truth_df.iscc.astype(bool)
    is_nu_mu_cc = is_cc & (abs(truth_df.pdg) == 14)

    nu_categ = pd.Series("none", index=truth_df.index, dtype="object")
    # Apply categories
    nu_categ[~is_nu] = "cosmic"
    nu_categ[is_nu & ~is_inside_fv] = "out_AV_nu"
    nu_categ[is_nu & is_inside_fv & ~is_signal] = "other_nu"
    nu_categ[is_nu & is_inside_fv & is_nu_mu_cc & (truth_df.npi_P_85MeV_10000MeV == 0)] = "CC_mu_0pi"
    nu_categ[is_nu & is_inside_fv & is_nu_mu_cc & (truth_df.npi_P_85MeV_10000MeV > 1)] = "CC_mu_2pi"
    nu_categ[is_nu & is_inside_fv & is_nu_mu_cc & is_signal & (truth_df.np_P_325MeV_10000MeV == 0)] = "0p_CC1Pi"
    nu_categ[is_nu & is_inside_fv & is_nu_mu_cc & is_signal & (truth_df.np_P_325MeV_10000MeV == 1)] = "1p_CC1Pi"
    nu_categ[is_nu & is_inside_fv & is_nu_mu_cc & is_signal & (truth_df.np_P_325MeV_10000MeV > 1)] = "plus2p_CC1Pi"
    
    if(is_truth_df):
        df['nu_categ_proton_reduced'] = nu_categ
    else:
        df[('truth', 'nu_categ_proton_reduced', '', '','','')] = nu_categ 
    return df

def add_genie_categ_column(df, is_truth_df=False):
    truth_df = df if is_truth_df else df.truth
 
    is_inside_fv = TruthInFV(truth_df.position)
    is_nu = IsNu(truth_df)
    is_cc = truth_df.iscc.astype(bool)
    is_nu_mu_cc = is_cc & (abs(truth_df.pdg) == 14)
    base = is_nu & is_inside_fv
 
    genie_categ = pd.Series("other", index=truth_df.index, dtype="object")
    genie_categ[~is_nu] = "cosmic"
    genie_categ[is_nu & ~is_inside_fv] = "out_AV_nu"
    genie_categ[base & ~is_cc & (abs(truth_df.pdg) == 14)] = "nu_mu_NC"
    for mode, name in [(0, "nu_mu_CC_QE"), (10, "nu_mu_CC_MEC"), (1, "nu_mu_CC_Res"), (2, "nu_mu_CC_Dis")]:
        genie_categ[base & is_nu_mu_cc & (df.genie_mode == mode)] = name
 
    if is_truth_df:
        df['genie_categ'] = genie_categ
    else:
        df[('truth', 'genie_categ', '', '', '', '')] = genie_categ   # was: nu_categ (NameError)
    return df
    
def add_max_angle_between_candidates_column(df, group_levels=['entry', 'rec.slc..index']):
    # 1. Initialize result series at the slice level
    angle_series = pd.Series(np.nan, index=df.index)
    
    # 2. Filter for MIP and Proton candidates
    mip_df = df[CutMasks.is_MIP_candidate_mask(df)]
    proton_df = df[CutMasks.is_primary_proton_mask(df)]
    
    # 3. Group both dataframes by the slice levels
    mip_groups = mip_df.groupby(level=group_levels)
    proton_groups = proton_df.groupby(level=group_levels)
    
    # Get the union of all unique slice indices present in either dataframe
    all_slice_indices = set(mip_groups.groups.keys()).union(proton_groups.groups.keys())
    
    # Define the column coordinates for the direction vectors and track length
    dir_cols = [
        ('pfp', 'trk', 'dir', 'x', '', ''),
        ('pfp', 'trk', 'dir', 'y', '', ''),
        ('pfp', 'trk', 'dir', 'z', '', '')
    ]
    len_col = ('pfp', 'trk', 'len', '', '', '')  # Standard column label for track length
    
    # 4. Iterate through every unique slice index
    for idx in all_slice_indices:
        # Retrieve the mip and proton subgroups if they exist for this slice index
        mip_group = mip_groups.get_group(idx) if idx in mip_groups.groups else pd.DataFrame()
        proton_group = proton_groups.get_group(idx) if idx in proton_groups.groups else pd.DataFrame()
        
        # Condition 1: Must have at least 1 MIP candidate
        if len(mip_group) < 1:
            continue
            
        # Condition 2: Total candidate particles must be at least 2
        total_len = len(mip_group) + len(proton_group)
        if total_len < 2:
            continue
        max_angle = 0.0
        
        # Case A: 2 or more MIPs -> Only use MIPs
        if len(mip_group) >= 2:
            dirs = mip_group.loc[:, dir_cols].astype(float).values
            
            # Check every unique pair among MIPs
            for dir_1, dir_2 in combinations(dirs, 2):
                mag1 = np.linalg.norm(dir_1)
                mag2 = np.linalg.norm(dir_2)
                if mag1 > 0 and mag2 > 0:
                    dot = np.dot(dir_1, dir_2)
                    angle = np.arccos(np.clip(dot / (mag1 * mag2), -1.0, 1.0))
                    if angle > max_angle:
                        max_angle = angle

        # Case B: Exactly 1 MIP -> Combine with the longest particle in the proton group
        elif len(mip_group) == 1 and len(proton_group) >= 1:
            # Extract single MIP direction vector
            mip_dir = mip_group.loc[:, dir_cols].astype(float).values[0]
            
            # Find the row index of the longest proton track in this slice group
            longest_proton_idx = proton_group[len_col].astype(float).idxmax()
            
            # Extract that specific proton's direction vector
            proton_dir = proton_group.loc[longest_proton_idx, dir_cols].astype(float).values
            
            # Compute the single angle between the 1 MIP and the longest Proton
            mag1 = np.linalg.norm(mip_dir)
            mag2 = np.linalg.norm(proton_dir)
            if mag1 > 0 and mag2 > 0:
                dot = np.dot(mip_dir, proton_dir)
                max_angle = np.arccos(np.clip(dot / (mag1 * mag2), -1.0, 1.0))
                
        else:
            # Safety catch-all (e.g., 1 MIP and 0 Protons, which total_len < 2 already filtered)
            continue
        
        # 7. Broadcast the calculated max angle back to this entire slice index
        angle_series.loc[idx] = max_angle
        
    # Add the column to the original dataframe using your multi-index structure
    df[('slc', 'measure_var', 'max_angle_between_candidates', '', '', '')] = angle_series

    return df
    
def add_n_protons_column(df):
    # 1. Calculate your masks and size
    is_track_mask = (df.pfp.trk.len > CTE.min_track_lenght) & (df.pfp.trackScore > CTE.min_track_score)
    is_primary_mask = (df.pfp.parent_is_primary == True) & (df.pfp.dist_to_vertex < CTE.max_primary_distance_to_vertex)
    
    # 2. Calculate num_protons (this is at the Slice level)
    num_protons = df[is_primary_mask & is_track_mask].groupby(level=group_levels).size() - 2
    
    # 3. Broadcast and join to the original 4-level index
    # We use reindex to "stretch" the 3-level data across the 4-level structure
    target_key = ('slc', 'measure_var', 'num_protons', '', '', '')
    
    df.loc[:,target_key] = df.index.droplevel('rec.slc.reco.pfp..index').map(num_protons)
    
    # 3. Fill NaNs
    df.loc[:,target_key] = df[target_key].fillna(-2)
    df.loc[df.slc.measure_var.num_protons > 2, target_key] = 2.
    return df
    


def add_n_true_proton_column(df): 
    df[('true_var', 'num_protons','')] = df[('np_P_325MeV_10000MeV','','')]
    df.loc[df.true_var.num_protons >= 2,('true_var', 'num_protons','')] = 2
    return df

def add_best_chi2_columns(df, update_calo = None):
    chi2_vars = [
        'chi2_muon',
        'chi2_pion',
        'chi2_proton',
    ]

    best_plane = df.pfp.trk.bestplane
    conditions = [
        best_plane == 2,
        best_plane == 1,
        best_plane == 0,
    ]
    
    for var in chi2_vars:
        if update_calo is None:
            choices = [
                df[('pfp','trk', 'chi2pid', 'I2', var, '')],
                df[('pfp','trk', 'chi2pid', 'I1', var, '')],
                df[('pfp','trk', 'chi2pid', 'I0', var, '')],
            ]
        else:
            choices = [
                df[('pfp','trk', 'chi2pid', 'I2', var + '_new', '')],
                df[('pfp','trk', 'chi2pid', 'I1', var + '_new', '')],
                df[('pfp','trk', 'chi2pid', 'I0', var + '_new', '')],
            ]
        
        df[('pfp','trk', 'chi2pid', 'best', var, '')] = np.select(
            conditions,
            choices,
            default=np.nan
        )
    return df


def add_nhit_column(slc_df, hit_dfs, hit_names, dedx_max=1000):
    if len(hit_dfs) != len(hit_names):
        raise ValueError("hit_dfs and hit_names must have the same length")
 
    cols = []
    for hit_df, name in zip(hit_dfs, hit_names):
        col = ('pfp', name, '', '', '', '')
        counts = hit_df[hit_df.dedx < dedx_max].groupby(level=pfp_levels).size()
        slc_df[col] = counts.reindex(slc_df.index, fill_value=0)
        cols.append(col)
 
    slc_df[('pfp', 'nhit', '', '', '', '')] = slc_df[cols].sum(axis=1)
    return slc_df
    

def add_max_daughter_hits_column(df):
    # Column names
    group_cols = ['entry','rec.slc..index']
    id_col = ('pfp','id','','','','')
    nhit_col = ('pfp','nhit','','','','')
    parent_col = ('pfp','parent','','','','')
    new_col = ('pfp','max_daughter_hits','','','','')

    # 1️⃣ Make a copy of nhit and parent info with a simple DataFrame
    tmp = df[[nhit_col, parent_col]].copy()

    # 2️⃣ Reset the index to make grouping easier
    tmp = tmp.reset_index()

    # 3️⃣ Group by parent
    # Group by entry, rec.slc..index, parent_id
    max_nhit_per_parent = (
        tmp.groupby(['entry','rec.slc..index', parent_col])
           [[nhit_col]]
           .max()
    )
    
    # 4️⃣ Rename index level for clarity
    max_nhit_per_parent.index = max_nhit_per_parent.index.rename(['entry','rec.slc..index','pfp_id'])
    last_level = max_nhit_per_parent.index.names[-1]
    max_nhit_per_parent = max_nhit_per_parent[max_nhit_per_parent.index.get_level_values(last_level) != -1]
    max_nhit_per_parent = max_nhit_per_parent.rename(columns={'nhit': 'max_daughter_hits'},level=1)

    df = df.reset_index()
    max_nhit_per_parent = max_nhit_per_parent.reset_index()

    df = df.merge(
        max_nhit_per_parent,
        how='left',
        left_on=group_cols + [id_col],
        right_on=group_cols + ['pfp_id']
    )

    df = df.set_index(['entry','rec.slc..index','rec.slc.reco.pfp..index'])
    df = df.drop(columns=['pfp_id'])
    df[new_col] = df[('pfp','max_daughter_hits','','','','')].fillna(0).astype(int)
    return df


def get_best_hit_df(slc_df, hit_dfs):
    for df in hit_dfs:
        df.index = df.index.set_names('rec.slc.reco.pfp.trk.calo.points..index', level=-1)
 
    best_plane = slc_df.pfp.trk.bestplane
    best_hit_df = (
        hit_dfs[0].where(best_plane == 0)
        .combine_first(hit_dfs[1].where(best_plane == 1))
        .combine_first(hit_dfs[2].where(best_plane == 2))
    )
    return best_hit_df.sort_values('rr', ascending=True)

def get_reduced_hit_df(hit_df):
    reduced_best_hit_df = hit_df[hit_df.rr < 30].copy()
    return reduced_best_hit_df

    
def compute_frac50(group):
    group = group.sort_values('rr')
    total_E = group['dE'].sum()
    max_rr = group['rr'].max()
    cumE = group['dE'].cumsum()

    # find first point where cumE >= 50% total
    mask = cumE >= 0.5 * total_E
    if not mask.any():
        return -1

    if max_rr == 0:
        return -1
        
    rr_50 = group.loc[mask, 'rr'].iloc[0]
    return rr_50 / max_rr


def add_frac50_column(df, best_hit_df):
    # Make sure we're not modifying a view
    df = df.copy()

    frac50 = (
        best_hit_df
        .groupby(level=pfp_levels)
        .apply(compute_frac50)
    )

    df[('pfp','trk','frac50','','','')] = frac50.reindex(df.index)

    return df

def chi2_exp_over_pol_lmfit(group):
    x = group['rr'].values
    y = group['dedx'].values

    #print(x,y)
    if len(x) < 3:
        return np.nan  # Not enough points to fit

    # --- Constant model ---
    model_const = lmfit.models.ConstantModel(prefix="c0_")
    result_const = model_const.fit(y, x=x)  # Constant model has no free parameters besides c0_ 

    # --- Exp + constant model ---
    model_exp = lmfit.Model(
        lambda x, p1, p2:p1*np.exp(p2*x),
        name="Const+Exp"
    )

    params = model_exp.make_params(
        p1=3.5,
        p2=-0.05
    )
    # Create parameters with initial guesses
    params['p1'].set(min=0.0)   # p1 > 0
    
    # Fit using params
    result_exp = model_exp.fit(y, params, x=x)
    # Return ratio of chi-squareds 
    return (result_exp.chisqr / (result_const.chisqr)) *(result_const.nfree / result_exp.nfree) 

def chi2_exp_over_pol_lmfit_3var(group):
    x = group['rr'].values
    y = group['dedx'].values

    #print(x,y)
    if len(x) < 3:
        return np.nan  # Not enough points to fit

    # --- Constant model ---
    model_const = lmfit.models.ConstantModel(prefix="c0_")
    result_const = model_const.fit(y, x=x)  # Constant model has no free parameters besides c0_ 

    # --- Exp + constant model ---
    model_exp = lmfit.Model(
        lambda x, p0, p1, p2: p0 + p1*np.exp(p2*x),
        name="Const+Exp"
    )

    params = model_exp.make_params(
        p0=0,
        p1=3.5,
        p2=-0.05
    )
    # Create parameters with initial guesses
    params['p0'].set(min=0.0)   # p0 > 0
    params['p1'].set(min=0.0)   # p1 > 0
    
    # Fit using params
    result_exp = model_exp.fit(y, params, x=x)
    # Return ratio of chi-squareds 
    if (result_const.chisqr != 0) & (result_exp.nfree != 0):
        return (result_exp.chisqr / (result_const.chisqr)) *(result_const.nfree / result_exp.nfree) 
    else:
        return -1


def add_chi2_exp_pol_column(df, best_hit_df):
    # 1️⃣ Ensure we're not modifying a view
    df = df.copy()

    # 2️⃣ Compute per-PFP value
    chi2_exp_pol_col = (
        best_hit_df
        .groupby(level=pfp_levels)
        .apply(chi2_exp_over_pol_lmfit)
    )

    # 3️⃣ Safe, aligned assignment
    df[('pfp','trk','chi2_exp_pol','','','')] = chi2_exp_pol_col.reindex(df.index)

    return df

def add_chi2_exp_pol_column_3var_fit(df, best_hit_df):
    # 1️⃣ Ensure we're not modifying a view
    df = df.copy()

    # 2️⃣ Compute per-PFP value
    chi2_exp_pol_col = (
        best_hit_df
        .groupby(level=pfp_levels)
        .apply(chi2_exp_over_pol_lmfit_3var)
    )

    # 3️⃣ Safe, aligned assignment
    df[('pfp','trk','chi2_exp_pol_3var','','','')] = chi2_exp_pol_col.reindex(df.index)

    return df

'''
def add_BDT_TMVA_proton_column(df, xmlWeightsFile, output_column):   
    # Initialize TMVA
    TMVA.Tools.Instance()
    reader_proton = ROOT.TMVA.Reader("!Color:!Silent")

    # Use NumPy float32 arrays of size 1 as pointers
    chi2_p        = np.zeros(1, dtype=np.float32)
    chi2_mu       = np.zeros(1, dtype=np.float32)
    chi2_exp_pol  = np.zeros(1, dtype=np.float32)
    frac_50       = np.zeros(1, dtype=np.float32)

    # Add variables to the reader (pass the array as pointer)
    reader_proton.AddVariable("chi2p", chi2_p)
    reader_proton.AddVariable("chi2mu", chi2_mu)
    reader_proton.AddVariable("chi2_exp_pol0", chi2_exp_pol)
    reader_proton.AddVariable("fraction_50_drift_percent", frac_50)

    # Book the MVA
    reader_proton.BookMVA(ROOT.TString("BDT"), ROOT.TString(xmlWeightsFile))

    # Fill values in one go
    bdt_scores = np.zeros(len(df), dtype=np.float32)

    for i, (_, row) in enumerate(tqdm(df.iterrows(), total=len(df), desc="Processing BDT scores")):
        chi2_p[0]       = row[('pfp','trk','chi2pid','best','chi2_proton','')]
        chi2_mu[0]      = row[('pfp','trk','chi2pid','best','chi2_muon','')]
        chi2_exp_pol[0] = row[('pfp','trk','chi2_exp_pol','','','')]
        frac_50[0]      = row[('pfp','trk','frac50','','','')]

        vals = [chi2_p[0], chi2_mu[0], chi2_exp_pol[0], frac_50[0]]

        if not np.all(np.isfinite(vals)):
            bdt_scores[i] = -1
            continue
        
        bdt_scores[i] = reader_proton.EvaluateMVA(ROOT.TString("BDT"))
    # Assign the full array in **one step**
    df[output_column] = bdt_scores
    return(df)
    
'''
# Unpack back if you want the same variable names
def add_best_ke_column(df):
    best_plane = df.pfp.trk.bestplane
    conditions = [
        best_plane == 2,
        best_plane == 1,
        best_plane == 0,
    ]

    choices = [
        df[('pfp','trk', 'calo', 'I2', 'ke', '')],
        df[('pfp','trk', 'calo', 'I1', 'ke', '')],
        df[('pfp','trk', 'calo', 'I0', 'ke', '')],
    ]

    df.loc[:, ('pfp','trk', 'calo', 'best', 'ke', '')] = np.select(
            conditions,
            choices,
            default=np.nan
    )
    return df

def compute_mean_dedx(group):
    return group['dedx'].mean()

def add_mean_dedx_column(df, best_hit_df):
    mean_dedx = (
        best_hit_df
        .groupby(level=pfp_levels)
        .apply(compute_mean_dedx)
    )
    df.loc[:, ('pfp','trk', 'mean_dEdx', '', '', '')] = mean_dedx
    return df

def add_p_type_column(df):        
    # ----------------------    
    pdg_col = ('pfp', 'trk', 'truth', 'p', 'pdg', '')
    endp_col = ('pfp', 'trk', 'truth', 'p', 'end_process', '')
    ptype_col = ('pfp', 'trk', 'truth', 'p', 'p_type', '')
    
    pdg = df[pdg_col].abs()
    endp = df[endp_col]
    df.loc[:,ptype_col] = 'other'

    #if(muon_row.pfp.trk.truth.p.pdg.iloc[0] > -2147483648):
    # Muon
    df.loc[pdg == 13, ptype_col] = 'muon'

    # Inelastic pion
    inelastic_pion_mask = (
        (pdg == 211) &
        (endp.isin([9, 10]))
    )
    df.loc[inelastic_pion_mask, ptype_col] = 'inelastic pion'

    # Stopping pion
    stopping_pion_mask = (
        (pdg == 211) &
        (endp.isin([3, 45]))
    )
    df.loc[stopping_pion_mask, ptype_col] = 'stopping pion'

    # Proton
    df.loc[pdg == 2212, ptype_col] = 'proton'
    shower_mask = pdg.isin([ 22, 11])
    df.loc[shower_mask, ptype_col] = 'shower'
    
    return df

def add_scatter_angle_ratio_column(df,mcs_df):
    # Output column (matches your dataframe structure)
    out_col = ('pfp', 'scatter_angle_ratio', '', '', '', '')
    mcs_pos = mcs_df[mcs_df['seg_scatter_angles'] >= 0]
    
    scatter_angle_ratio = (
        mcs_pos
        .groupby(level=pfp_levels)['seg_scatter_angles']
        .apply(lambda x: x.max() / x.sum() if x.sum() > 0 else np.nan)
    )

    df[out_col] = scatter_angle_ratio
    return df

def add_bdt_score(evt_df, bdt_model_file, input_cols, output_col):
    with open(bdt_model_file, "rb") as f:
        bdt = pickle.load(f)

    mask = np.ones(len(evt_df), dtype=bool)
    for col in input_cols:
        mask &= evt_df[col].notna()
        # Only use >= 0 if your specific physics vars (like chi2) can't be negative
        mask &= (evt_df[col] >= 0) 

    X_BDT = evt_df.loc[mask, input_cols].values

    # CRITICAL: Match the logic of your overtraining plot
    if hasattr(bdt, "decision_function"):
        scores = bdt.decision_function(X_BDT)
        default_val = -999.0 # Raw scores can be negative, so use a clear outlier
    else:
        scores = bdt.predict_proba(X_BDT)[:, 1]
        default_val = -1.0

    evt_df.loc[:, output_col] = default_val
    evt_df.loc[mask, output_col] = scores

    return evt_df

'''
def add_TMVA_BDT_muon_pion_column(df, xmlWeightsFile_mp, output_column):  
    # Initialize TMVA

    TMVA.Tools.Instance()
    reader_muon_pion = ROOT.TMVA.Reader("!Color:!Silent")

    # Use NumPy float32 arrays of size 1 as pointers
    chi2_p_mp        = np.zeros(1, dtype=np.float32)
    chi2_mu_mp       = np.zeros(1, dtype=np.float32)
    frac_50_mp       = np.zeros(1, dtype=np.float32)
    track_mcs_scatter_max_ratio = np.zeros(1, dtype=np.float32) 
    max_daughter_hits = np.zeros(1, dtype=np.float32) 

    # Add variables to the reader (pass the array as pointer)
    reader_muon_pion.AddVariable("chi2p", chi2_p_mp)
    reader_muon_pion.AddVariable("chi2mu", chi2_mu_mp)
    reader_muon_pion.AddVariable("fraction_50_drift_percent", frac_50_mp)
    reader_muon_pion.AddVariable("track_mcs_scatter_max_ratio", track_mcs_scatter_max_ratio)
    reader_muon_pion.AddVariable("max_daughter_hits", max_daughter_hits)

    # Book the MVA
    reader_muon_pion.BookMVA(ROOT.TString("BDT"), ROOT.TString(xmlWeightsFile_mp))

    # Fill values in one go
    bdt_scores_muon_pion = np.zeros(len(df), dtype=np.float32)

    for i, (_, row) in enumerate(tqdm(df.iterrows(), total=len(df), desc="Processing BDT scores")):
        chi2_p_mp[0]       = row[('pfp','trk','chi2pid','best','chi2_proton','')]
        chi2_mu_mp[0]      = row[('pfp','trk','chi2pid','best','chi2_muon','')]
        frac_50_mp[0]      = row[('pfp','trk','frac50','','','')]
        track_mcs_scatter_max_ratio[0]   = row[('pfp', 'scatter_angle_ratio','','','','')]
        max_daughter_hits[0]      = row[('pfp','max_daughter_hits','','','','')]

        vals = [chi2_p_mp[0], chi2_mu_mp[0], frac_50_mp[0],track_mcs_scatter_max_ratio[0], max_daughter_hits[0]]

        if not np.all(np.isfinite(vals)):
            bdt_scores_muon_pion[i] = -1
            continue
        
        bdt_scores_muon_pion[i] = reader_muon_pion.EvaluateMVA(ROOT.TString("BDT"))

    # Assign the full array in **one step**
    df[output_column] = bdt_scores_muon_pion
    return df
'''

_PFP_LEVEL = 'rec.slc.reco.pfp..index'
 
 
def _add_slice_count(df, mask, key, offset=0, fill=0):
    counts = df[mask].groupby(level=group_levels).size() + offset
    df[key] = df.index.droplevel(_PFP_LEVEL).map(counts)
    df[key] = df[key].fillna(fill)
    return df
 
 
def add_n_primary_tracks_column(df):
    return _add_slice_count(df, CutMasks.is_primary_track_mask(df),
                            ('slc', 'cut_var', 'n_prim_tracks', '', '', ''))
 
 
def add_n_primary_showers_column(df):
    mask = ((df.pfp.parent_is_primary == True)
            & (df.pfp.trackScore >= 0) & (df.pfp.trackScore < CTE.max_shower_track_score)
            & (df.pfp.shw.bestplane_energy > CTE.min_shower_ke))
    return _add_slice_count(df, mask, ('slc', 'cut_var', 'n_prim_showers', '', '', ''))
 
 
def add_n_primary_MIP_column(df):
    return _add_slice_count(df, CutMasks.is_MIP_candidate_mask(df),
                            ('slc', 'cut_var', 'n_MIP_candidates', '', '', ''))
 
 
def add_n_primary_MIP_proton_column(df):
    mask = CutMasks.is_MIP_candidate_mask(df) & (df.pfp.trk.bdt_proton_score < CTE.BDT_proton_max_score)
    return _add_slice_count(df, mask, ('slc', 'cut_var', 'n_MIP_candidates_proton', '', '', ''))
 
 
def add_n_exiting_pfps_column(df):
    return _add_slice_count(df, CutMasks.exiting_pfp_mask(df),
                            ('slc', 'cut_var', 'n_exiting_pfps', '', '', ''))
 
 
def add_n_exiting_z_pfps_column(df):
    return _add_slice_count(df, CutMasks.exiting_z_pfp_mask(df),
                            ('slc', 'cut_var', 'n_exiting_z_pfps', '', '', ''))
 
 
def add_n_MIP_candidate_michel_column(df):
    mask = (CutMasks.is_MIP_candidate_mask(df)
            & (df.pfp.trackScore < CTE.michel_max_track_score)
            & (df.pfp.trk.calo.best.ke < CTE.michel_max_visible_energy)
            & (df.pfp.max_daughter_hits == 0))
    return _add_slice_count(df, mask, ('slc', 'cut_var', 'n_MIP_candidate_michel', '', '', ''))
 
 
def add_n_extra_pions_column(df):
    chi2 = df.pfp.trk.chi2pid.best
    mask = (CutMasks.is_primary_track_mask(df)
            & (chi2.chi2_muon < CTE.MIP_candidate_max_muon_score)
            & (chi2.chi2_proton > CTE.MIP_candidate_min_proton_score))
    # NOTE: fill=0 means "no relaxed MIPs" -> 0 but "1 relaxed MIP" -> -1.
    # Check this is what your cut expects (num_protons uses fill=-2).
    return _add_slice_count(df, mask, ('slc', 'cut_var', 'n_extra_pions', '', '', ''),
                            offset=-2, fill=-2)
 
 
def add_n_protons_column(df):
    is_track = (df.pfp.trk.len > CTE.min_track_lenght) & (df.pfp.trackScore > CTE.min_track_score)
    is_primary = (df.pfp.parent_is_primary == True) & (df.pfp.dist_to_vertex < CTE.max_primary_distance_to_vertex)
    key = ('slc', 'measure_var', 'num_protons', '', '', '')
    df = _add_slice_count(df, is_primary & is_track, key, offset=-2, fill=-2)
    df.loc[df[key] > 2, key] = 2.
    return df
 

def calculate_hypfit_p(track_id, best_plane, hit_df, target_pdg=211, cleaning="all", use_data=False):
    """Return the Hypfit momentum (GeV) for one track, or -1.0 on any failure."""
    try:
        plane = int(best_plane.iloc[0] if hasattr(best_plane, "iloc") else best_plane)
        if plane not in (0, 1, 2):
            return -1.0
        pdg = int(target_pdg)

        track_hits = hit_df.loc[track_id]
        arrays = [np.ascontiguousarray(track_hits[c].to_numpy(dtype=np.float64))
                  for c in ("rr", "dedx", "pitch")]
        finite = np.logical_and.reduce([np.isfinite(a) for a in arrays])
        rr, dedx, pitch = (a[finite] for a in arrays)
        if len(rr) == 0:
            return -1.0

        v_rr = ROOT.std.vector('double')(rr)
        v_dedx = ROOT.std.vector('double')(dedx)
        v_pitch = ROOT.std.vector('double')(pitch)
        bad_hits = h_fit.get_hits_to_ignore(v_rr, v_dedx, str(cleaning))

        conv_map = select_conv_tf1_map(pdg, use_data)
        TLE_p = h_fit.NormLikelihood_w_tf1_vec(
            v_dedx, v_rr, v_pitch, bad_hits, True, pdg, conv_map[plane])

        if TLE_p > 0:
            return h_fit.map_PhysdEdx.at(target_pdg).KEtoMomentum(TLE_p) / 1000.0
        return TLE_p 

    
    except KeyError:
        return -1.0
    except Exception as e:
        print(f"Error for track {track_id}: {e}")
        return -1.0
        
    finally:
        # 1. Clear the vectors
        if 'v_rr' in locals():
            v_rr.clear()
            v_dedx.clear()
            v_pitch.clear()
        
        # 2. Force ROOT to cleanup any temporary objects/TGraphs/TF1s 
        # created during the C++ call
        ROOT.gDirectory.Clear()
        
        # 3. Specifically clear the global list of functions to prevent 
        # accumulation of TF1 objects created in C++
        ROOT.gROOT.GetListOfFunctions().Clear()



def _mask_invalid_tki(df):
    """Set TKI variables to -999 where the muon/pion momenta are not valid."""
    mv = df.slc.measure_var
    valid = (mv.p_mu_x > -900) & (mv.p_pi_x > -900) & (mv.TLE_p_pi > 0)
    for name in ('delta_pT', 'delta_alpha_T', 'delta_phi_T'):
        df.loc[~valid, ('slc', 'measure_var', name, '', '', '')] = -999.0
    return df
 
 
def _angle_between(a, b):
    """Row-wise angle between two (N,3) arrays; -999 where either norm is 0."""
    norm = np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1)
    out = np.full(len(a), -999.0)
    ok = norm > 0
    dot = np.einsum('ij,ij->i', a[ok], b[ok])
    out[ok] = np.arccos(np.clip(dot / norm[ok], -1, 1))
    return out
 
 
def add_transverse_vars_column(df):
    len_col = ('pfp', 'trk', 'len', '', '', '')
    slice_index = df.index.droplevel(_PFP_LEVEL)
 
    # --- Proton-like primary tracks: non-MIP primary tracks ---
    proton_mask = ~CutMasks.is_MIP_candidate_mask(df) & CutMasks.is_primary_track_mask(df)
 
    # --- Remove the pion from the proton sum (fix for double counting) ---
    # In get_mu_pi_vars, when a slice has exactly 1 MIP, the pion is the longest
    # track of is_primary_proton_mask. That track is already counted in p_pi, so
    # it must not be added again to the proton momentum.
    n_mip = CutMasks.is_MIP_candidate_mask(df).groupby(level=group_levels).transform('sum')
    pion_cands = df[(n_mip == 1) & CutMasks.is_primary_proton_mask(df)]
    if len(pion_cands) > 0:
        pion_idx = pion_cands[len_col].astype(float).groupby(level=group_levels).idxmax()
        proton_mask &= ~df.index.isin(pion_idx.tolist())
 
    proton_df = df[proton_mask]
 
    # --- Sum proton momentum components per slice ---
    p_col = ('pfp', 'trk', 'rangeP', 'p_proton', '', '')
    comps = {}
    for axis in 'xyz':
        comps[axis] = proton_df[p_col] * proton_df[('pfp', 'trk', 'dir', axis, '', '')]
    proton_sum = (
        pd.DataFrame(comps)
        .groupby(level=group_levels).sum()
        .reindex(slice_index.unique(), fill_value=0)
        .reindex(slice_index)
    )
 
    mv = df.slc.measure_var
    for axis in 'xyz':
        df[('slc', f'p_proton_total_{axis}', '', '', '', '')] = proton_sum[axis].values
        df[('slc', f'p_hadron_{axis}', '', '', '', '')] = (
            mv[f'p_pi_{axis}'].values + proton_sum[axis].values
        )
 
    # --- 2D (x, y, 0) vectors for muon and hadron ---
    zeros = np.zeros(len(df))
    p_mu_2d = np.stack([mv.p_mu_x.values, mv.p_mu_y.values, zeros], axis=1)
    p_had_2d = np.stack([df[('slc', 'p_hadron_x', '', '', '', '')].values,
                         df[('slc', 'p_hadron_y', '', '', '', '')].values, zeros], axis=1)
 
    delta_pT_vec = p_mu_2d + p_had_2d
    df[('slc', 'measure_var', 'delta_pT', '', '', '')] = np.linalg.norm(delta_pT_vec, axis=1)
    df[('slc', 'measure_var', 'delta_alpha_T', '', '', '')] = _angle_between(-p_mu_2d, delta_pT_vec)
    df[('slc', 'measure_var', 'delta_phi_T', '', '', '')] = _angle_between(-p_mu_2d, p_had_2d)
 
    return _mask_invalid_tki(df)
 
    
def get_mu_pi_vars(group, best_hit_df, use_data = False):
    """
    Selects a Muon and a Pion from a group of PFPs based on:
    - If >= 2 MIPs: Muon and Pion are both selected from the MIP group.
    - If == 1 MIP: That MIP is the Muon; the Pion is the longest track in the Proton group.
    """
    
    # --- 1. SELECTION SPLIT AND SAFETY CHECKS ---
    mip_df = group[CutMasks.is_MIP_candidate_mask(group)]
    proton_df = group[CutMasks.is_primary_proton_mask(group)]
    
    # Safety condition: Total group must be >= 2 AND there must be at least 1 MIP
    if len(group) < 2 or len(mip_df) < 1:
        return pd.Series({
            'p_mu_x': -999.0, 'p_mu_y': -999.0, 'p_mu_z': -999.0,
            'p_pi_x': -999.0, 'p_pi_y': -999.0, 'p_pi_z': -999.0,
            'reco_p_mu': -999.0, 'cos_theta_mu': -999.0,
            'mu_chi2_proton': -999.0, 'mu_chi2_mu': -999.0, 'mu_chi2_exp_pol': -999.0,
            'mu_scatter_angle_ratio': -999.0, 'mu_max_daughter_hits': -999.0,
            'mu_frac50': -999.0, 'mu_bdt_score_proton': -999.0, 'mu_bdt_score_muon_pion': -999.0,
            'range_p_pi': -999.0, 'TLE_p_pi': -999.0, 'cos_theta_pi': -999.0,
            'pi_chi2_proton': -999.0, 'pi_chi2_mu': -999.0, 'pi_chi2_exp_pol': -999.0,
            'pi_scatter_angle_ratio': -999.0, 'pi_max_daughter_hits': -999.0,
            'pi_frac50': -999.0, 'pi_bdt_score_proton': -999.0, 'pi_bdt_score_muon_pion': -999.0,
            'muon_contained': False,
            'mu_true_pdg': -1, 'mu_true_p_type': "none", 'mu_true_end_process': -1,
            'pi_true_pdg': -1, 'pi_true_p_type': "none", 'pi_true_end_process': -1,
            'pi_true_p': -999.0, 'mu_true_p': -999.0,
            'pi_true_costheta': -999.0, 'mu_true_costheta': -999.0             
        })

    # Column name shortcuts for readability
    bdt_score_col = ('pfp', 'trk', 'bdt_muon_pion_score', '', '', '')
    len_col       = ('pfp', 'trk', 'len', '', '', '')

    # --- 2. BRANCH SELECTION LOGIC ---
    
    # CASE A: 2 or more MIP candidates -> Perform selection ONLY using the MIP candidates
    if len(mip_df) >= 2:
        exiting_mask = CutMasks.exiting_pfp_mask(mip_df)
        
        # Select Muon from within MIP candidates
        if exiting_mask.sum() >= 1:
            muon_row = mip_df.loc[exiting_mask].iloc[[0]]
            muon_contained = False
            p_mu = muon_row.pfp.trk.mcsP.fwdP_muon.iloc[0]
        else:
            mip_sorted = mip_df.sort_values(bdt_score_col)
            muon_row = mip_sorted.iloc[[-1]]
            muon_contained = True
            p_mu = muon_row.pfp.trk.rangeP.p_muon.iloc[0]
            
        # Select Pion from remaining MIP candidates
        remaining_mips = mip_df.drop(muon_row.index)
        pion_row = remaining_mips.sort_values(len_col).iloc[[-1]]

    # CASE B: Exactly 1 MIP candidate -> That MIP is the Muon; Pion is the longest Proton
    else:
        # If there are no proton tracks to pick a pion from, return empty series
        if len(proton_df) < 1:
            return pd.Series({
                'p_mu_x': -999.0, 'p_mu_y': -999.0, 'p_mu_z': -999.0,
                'p_pi_x': -999.0, 'p_pi_y': -999.0, 'p_pi_z': -999.0,
                'reco_p_mu': -999.0, 'cos_theta_mu': -999.0,
                'mu_chi2_proton': -999.0, 'mu_chi2_mu': -999.0, 'mu_chi2_exp_pol': -999.0,
                'mu_scatter_angle_ratio': -999.0, 'mu_max_daughter_hits': -999.0,
                'mu_frac50': -999.0, 'mu_bdt_score_proton': -999.0, 'mu_bdt_score_muon_pion': -999.0,
                'range_p_pi': -999.0, 'TLE_p_pi': -999.0, 'cos_theta_pi': -999.0,
                'pi_chi2_proton': -999.0, 'pi_chi2_mu': -999.0, 'pi_chi2_exp_pol': -999.0,
                'pi_scatter_angle_ratio': -999.0, 'pi_max_daughter_hits': -999.0,
                'pi_frac50': -999.0, 'pi_bdt_score_proton': -999.0, 'pi_bdt_score_muon_pion': -999.0,
                'muon_contained': False,
                'mu_true_pdg': -1, 'mu_true_p_type': "none", 'mu_true_end_process': -1,
                'pi_true_pdg': -1, 'pi_true_p_type': "none", 'pi_true_end_process': -1,
                'pi_true_p': -999.0, 'mu_true_p': -999.0,
                'pi_true_costheta': -999.0, 'mu_true_costheta': -999.0             
            })
          
        # The single MIP candidate is the Muon
        muon_row = mip_df.iloc[[0]]
        exiting_mask = CutMasks.exiting_pfp_mask(muon_row)
        
        if exiting_mask.sum() >= 1:
            muon_contained = False
            p_mu = muon_row.pfp.trk.mcsP.fwdP_muon.iloc[0]
        else:
            muon_contained = True
            p_mu = muon_row.pfp.trk.rangeP.p_muon.iloc[0]
            
        # The Pion is the longest track inside the proton dataframe
        pion_row = proton_df.sort_values(len_col).iloc[[-1]]
        
    # --- 4. VARIABLE EXTRACTION ---
    # Create helper aliases to keep lines short
    
    # Muon kinematics
    cos_theta_mu = muon_row.pfp.trk.dir.z.iloc[0]
    p_mu_x = p_mu * muon_row.pfp.trk.dir.x.iloc[0]
    p_mu_y = p_mu * muon_row.pfp.trk.dir.y.iloc[0]
    p_mu_z = p_mu * muon_row.pfp.trk.dir.z.iloc[0]

    # Pion kinematics
    p_pi_range = pion_row.pfp.trk.rangeP.p_pion.iloc[0]
    cos_theta_pi = pion_row. pfp.trk.dir.z.iloc[0]
    pion_id = pion_row.index[0] # (ntuple, entry, slc, pfp)
    best_plane = pion_row.pfp.trk.bestplane.iloc[0]

    # Calculate Hypfit (TLE) Momentum
    p_pi_TLE = calculate_hypfit_p(track_id=pion_id, best_plane=best_plane, 
                                  hit_df=best_hit_df, target_pdg=211, cleaning="all", use_data = use_data)
    
    p_pi_x = p_pi_TLE * pion_row.pfp.trk.dir.x.iloc[0]
    p_pi_y = p_pi_TLE * pion_row.pfp.trk.dir.y.iloc[0]
    p_pi_z = p_pi_TLE * pion_row.pfp.trk.dir.z.iloc[0]

    # --- 5. TRUTH MATCHING ---
    # Check if the muon has a valid truth association
    if (muon_row.pfp.trk.truth.p.pdg > -2147483648).any():
        mu_true_pdg = muon_row.pfp.trk.truth.p.pdg.iloc[0]
        mu_true_p_type =muon_row.pfp.trk.truth.p.p_type.iloc[0]
        mu_true_end_process = muon_row.pfp.trk.truth.p.end_process.iloc[0]
        
        pi_true_pdg = pion_row.pfp.trk.truth.p.pdg.iloc[0]
        pi_true_p_type = pion_row.pfp.trk.truth.p.p_type.iloc[0]
        pi_true_end_process = pion_row.pfp.trk.truth.p.end_process.iloc[0]
        
        mu_true_p = magdf(muon_row.pfp.trk.truth.p.genp).iloc[0]
        pi_true_p = magdf(pion_row.pfp.trk.truth.p.genp).iloc[0]
        
        mu_true_costheta = muon_row.pfp.trk.truth.p.genp.z.iloc[0] / mu_true_p
        pi_true_costheta = pion_row.pfp.trk.truth.p.genp.z.iloc[0] / pi_true_p
    else:
        mu_true_pdg = pi_true_pdg = -1
        mu_true_p_type = pi_true_p_type = "none"
        mu_true_end_process = pi_true_end_process = -1
        mu_true_p = pi_true_p = -1
        mu_true_costheta = pi_true_costheta = -1

    # --- 6. RETURN RESULTS ---
    return pd.Series({
        'p_mu_x': p_mu_x, 'p_mu_y': p_mu_y, 'p_mu_z': p_mu_z,
        'p_pi_x': p_pi_x, 'p_pi_y': p_pi_y, 'p_pi_z': p_pi_z,
        
        'reco_p_mu': p_mu,
        'cos_theta_mu': cos_theta_mu,
        'mu_chi2_proton': muon_row.pfp.trk.chi2pid.best.chi2_proton.iloc[0],
        'mu_chi2_mu': muon_row.pfp.trk.chi2pid.best.chi2_muon.iloc[0],
        'mu_chi2_exp_pol': muon_row.pfp.trk.chi2_exp_pol.iloc[0],
        'mu_scatter_angle_ratio': muon_row.pfp.scatter_angle_ratio.iloc[0],
        'mu_max_daughter_hits': muon_row.pfp.max_daughter_hits.iloc[0],
        'mu_frac50': muon_row.pfp.trk.frac50.iloc[0],
        'mu_bdt_score_proton': muon_row.pfp.trk.bdt_proton_score.iloc[0],
        'mu_bdt_score_muon_pion': muon_row.pfp.trk.bdt_muon_pion_score.iloc[0],

        'range_p_pi': p_pi_range,
        'TLE_p_pi': p_pi_TLE,
        'cos_theta_pi': cos_theta_pi,        
        'pi_chi2_proton': pion_row.pfp.trk.chi2pid.best.chi2_proton.iloc[0],
        'pi_chi2_mu': pion_row.pfp.trk.chi2pid.best.chi2_muon.iloc[0],
        'pi_chi2_exp_pol': pion_row.pfp.trk.chi2_exp_pol.iloc[0],
        'pi_scatter_angle_ratio': pion_row.pfp.scatter_angle_ratio.iloc[0],
        'pi_max_daughter_hits': pion_row.pfp.max_daughter_hits.iloc[0],
        'pi_frac50': pion_row.pfp.trk.frac50.iloc[0],
        'pi_bdt_score_proton': pion_row.pfp.trk.bdt_proton_score.iloc[0],
        'pi_bdt_score_muon_pion': pion_row.pfp.trk.bdt_muon_pion_score.iloc[0],

        'muon_contained': muon_contained,

        'mu_true_pdg': mu_true_pdg,
        'mu_true_p_type': mu_true_p_type,
        'mu_true_end_process': mu_true_end_process,
        'pi_true_pdg': pi_true_pdg,
        'pi_true_p_type': pi_true_p_type,
        'pi_true_end_process': pi_true_end_process,
        
        'pi_true_p': pi_true_p,
        'mu_true_p': mu_true_p,
        'pi_true_costheta': pi_true_costheta,
        'mu_true_costheta': mu_true_costheta             
    })
    
import numpy as np

def add_transverse_observables_mc(df):
    p_p_x = df[('ptot','pmag','')] * df[('ptot','dir','x')]
    p_p_y = df[('ptot','pmag','')] * df[('ptot','dir','y')]
    p_pi_x = df[('pitot','pmag','')] * df[('pitot','dir','x')]
    p_pi_y = df[('pitot','pmag','')] * df[('pitot','dir','y')]
    p_had_x = p_p_x + p_pi_x
    p_had_y = p_p_y + p_pi_y

    p_mu_x = df[('mu', 'genp', 'x')]
    p_mu_y = df[('mu', 'genp', 'y')]

    # --- Build 2D stacked vectors (x, y, 0) ---
    p_mu_2d = np.stack([p_mu_x, p_mu_y, np.zeros(len(df))], axis=1)
    p_had_2d = np.stack([p_had_x, p_had_y, np.zeros(len(df))], axis=1)

    # --- Delta pT vector and magnitude ---
    delta_pT_vec = p_mu_2d + p_had_2d
    df[('true_var', 'delta_pT','')] = np.linalg.norm(delta_pT_vec, axis=1)

    # --- Delta alphaT ---
    dot_alpha = np.einsum('ij,ij->i', -p_mu_2d, delta_pT_vec)
    norm_alpha = np.linalg.norm(p_mu_2d, axis=1) * np.linalg.norm(delta_pT_vec, axis=1)
    safe_alpha = norm_alpha > 0
    delta_alpha_T = np.full(len(df), -999.0)
    delta_alpha_T[safe_alpha] = np.arccos(np.clip(dot_alpha[safe_alpha] / norm_alpha[safe_alpha], -1, 1))
    df[('true_var', 'delta_alpha_T','')] = delta_alpha_T
    # --- Delta phiT ---
    dot_phi = np.einsum('ij,ij->i', -p_mu_2d, p_had_2d)
    norm_phi = np.linalg.norm(p_mu_2d, axis=1) * np.linalg.norm(p_had_2d, axis=1)
    safe_phi = norm_phi > 0
    delta_phi_T = np.full(len(df), -999.0)
    delta_phi_T[safe_phi] = np.arccos(np.clip(dot_phi[safe_phi] / norm_phi[safe_phi], -1, 1))
    df[('true_var', 'delta_phi_T','')] = delta_phi_T
    
    return df



def make_cc1pinudf(f, include_weights=True, multisim_nuniv=100, genie_multisim_nuniv=100, wgt_types=["bnb","genie","g4"], slim=True, genie_systematics=None):
    
    nudf = make_mcnudf(f, include_weights=include_weights, multisim_nuniv=multisim_nuniv, genie_multisim_nuniv=genie_multisim_nuniv, wgt_types= wgt_types, slim=slim, genie_systematics=genie_systematics)
    
    mu_p_series = magdf(nudf.mu.genp)
    cpi_p_series = magdf(nudf.cpi.genp)

    mu_cos_theta_series = nudf.mu.genp.z / mu_p_series
    cpi_cos_theta_series = nudf.cpi.genp.z / cpi_p_series

    p_cpi = nudf.loc[:,[('cpi','genp','x'),('cpi','genp','y'),('cpi','genp','z')]].to_numpy(dtype=float)
    p_mu = nudf.loc[:, [('mu','genp','x'),('mu','genp','y'),('mu','genp','z')]].to_numpy(dtype=float)
    dot = np.einsum('ij,ij->i', p_cpi, p_mu)
    mag = np.linalg.norm(p_cpi, axis=1) * np.linalg.norm(p_mu, axis=1)
    
    mu_pi_angle_series = np.full(len(nudf), np.nan)
    valid = mag > 0
    mu_pi_angle_series[valid] = np.arccos(np.clip(dot[valid] / mag[valid], -1.0, 1.0))

    this_nudf = pd.DataFrame({
        'true_p_mu': mu_p_series,
        'true_p_pi': cpi_p_series,
        'true_cos_theta_mu': mu_cos_theta_series,
        'true_cos_theta_pi': cpi_cos_theta_series,
        'true_mu_pi_angle': mu_pi_angle_series,
    })
    this_nudf.columns = pd.MultiIndex.from_tuples([('true_var',col,'') for col in this_nudf.columns])        
    cols_to_add = this_nudf.columns.difference(nudf.columns)
    nudf = nudf.join(this_nudf[cols_to_add])
    
    #add nu categories
    nudf = add_nu_categ_column(nudf, True)
    nudf = add_nu_categ_proton_reduced_column(nudf, True)
    nudf = add_genie_categ_column(nudf, True)
    nudf = add_n_true_proton_column(nudf)

    #add TKI
    nudf = add_transverse_observables_mc(nudf)
    
    return shrink_df(nudf)


def make_cc1pinudf_expanded_syst(f):
    return make_cc1pinudf(f, include_weights=True, multisim_nuniv=100, genie_multisim_nuniv=100, wgt_types=["bnb","genie","g4"], slim=False, genie_systematics=None)

def make_cc1pinudf_Ar23p_expanded_syst(f):
    return make_cc1pinudf(f, include_weights=True, multisim_nuniv=100, genie_multisim_nuniv=100, wgt_types=["bnb","genie","g4"], slim=False, genie_systematics=ar23p_genie_systematics + regen_systematics_no_duplicates)

def make_cc1pinudf_Ar23p(f):
    return make_cc1pinudf(f, include_weights=True, multisim_nuniv=100, genie_multisim_nuniv=100, wgt_types=["bnb","genie","g4"], slim=True, genie_systematics=ar23p_genie_systematics + regen_systematics_no_duplicates)

    
def make_cc1pinudf_no_syst(f):
    return make_cc1pinudf(f, include_weights=False, multisim_nuniv=100, genie_multisim_nuniv=100, wgt_types=["bnb","genie","g4"], slim=False, genie_systematics=None)
    
def make_mcs_df(f):
    mcsdf = loadbranches(f["recTree"], [trkmcsbranches[1]]).rec.slc.reco.pfp.trk.mcsP
    return mcsdf


cols = [
        ('slc', 'tmatch', 'idx', '', '', ''),
        ('slc', 'self', '', '', '', ''),
        
        #Cosmic cut optimization
        ('slc', 'nu_score', '', '', '', ''),
        ('slc', 'barycenterFM', 'score', '', '', ''),
        ('slc', 'barycenterFM', 'flashTime', '', '', ''),
        
        #Vertex position        
        ('slc', 'vertex', 'x', '', '', ''),
        ('slc', 'vertex', 'y', '', '', ''),
        ('slc', 'vertex', 'z', '', '', ''),

        #Pfp truth info
        ('pfp', 'trk', 'truth', 'p', 'p_type', ''),
        ('pfp', 'trk', 'truth', 'p', 'pdg', ''),
        ('pfp', 'trk', 'truth', 'p', 'end_process', ''),
        
        #('pfp', 'trk', 'truth', 'p', 'genp', 'x'),
        #('pfp', 'trk', 'truth', 'p', 'genp', 'y'),
        #('pfp', 'trk', 'truth', 'p', 'genp', 'z'),

        #track optimization
        ('pfp', 'trk', 'len', '', '', ''),
        ('pfp', 'trackScore', '', '', '', ''),
        ('pfp', 'dist_to_vertex', '', '', '', ''),
        ('pfp', 'parent_is_primary', '', '', '', ''),
        
        #shower optimization
        ('pfp', 'shw', 'bestplane_energy', '', '', ''),
        
        #chi2 optimization
        ('pfp', 'trk', 'chi2pid', 'best', 'chi2_muon', ''),
        ('pfp', 'trk', 'chi2pid', 'best', 'chi2_proton', ''),
        
        #michel cut
        ('pfp', 'trk', 'calo', 'best', 'ke', ''),
        ('pfp', 'trk', 'mean_dEdx', '', '', ''),
        ('pfp', 'max_daughter_hits', '', '', '', ''),

        #angle cut
        ('slc', 'measure_var', 'max_angle_between_candidates', '', '', ''),
        ('slc', 'measure_var', 'angle_between_candidates', '', '', ''),
    
        #BDT vars
        ('pfp', 'trk', 'chi2_exp_pol', '', '', ''),
        #('pfp', 'trk', 'chi2_exp_pol_3var', '', '', ''),
        ('pfp', 'trk', 'frac50', '', '', ''),
        ('pfp', 'scatter_angle_ratio', '', '', '', ''),

        #pfp containment masks
        ('pfp', 'is_exiting', '', '', '', ''),
        ('pfp', 'is_exiting_z', '', '', '', ''),
    
        #BDT columns
        ('pfp','trk','bdt_muon_pion_score','','',''),
        ('pfp','trk','bdt_proton_score','','',''),
        
        #cuts boolean
        ('slc', 'cut', 'obvious_cosmic', '', '', ''),
        ('slc', 'cut', 't0', '', '', ''),
        ('slc', 'cut', 'inside_FV', '', '', ''),
        ('slc', 'cut', 'nu_score', '', '', ''),
        ('slc', 'cut', 'track', '', '', ''),
        ('slc', 'cut', 'shower', '', '', ''),
        ('slc', 'cut', 'MIP_candidates', '', '', ''),
        ('slc', 'cut', 'angle', '', '', ''),
        ('slc', 'cut', 'proton_BDT', '', '', ''),
        ('slc', 'cut', 'proton_BDT_sideband', '', '', ''),
        ('slc', 'cut', 'containment', '', '', ''),
        ('slc', 'cut', 'michel', '', '', ''),
        ('slc', 'cut', 'extra_pion', '', '', ''),
        ('slc','cut','energy','','',''),
        ('slc', 'cut', 'TPC_containment', '', '', ''),
        ('slc','cut','no_high_yz','','',''),

        #Measure variables
        ('slc', 'measure_var', 'num_protons', '', '', ''),

        ('slc','measure_var','reco_p_mu','','',''),
        ('slc','measure_var','reco_cos_theta_mu','','',''),
        ('slc','measure_var','mu_chi2_proton','','',''),
        ('slc','measure_var','mu_chi2_mu','','',''),
        ('slc','measure_var','mu_chi2_exp_pol','','',''),
        ('slc','measure_var','mu_scatter_angle_ratio','','',''),
        ('slc','measure_var','mu_max_daughter_hits','','',''),
        ('slc','measure_var','mu_frac50','','',''),
        ('slc','measure_var','mu_bdt_score_proton','','',''),
        ('slc','measure_var','mu_bdt_score_muon_pion','','',''),

        ('slc','measure_var','range_p_pi','','',''),
        ('slc','measure_var','TLE_p_pi','','',''),
        ('slc','measure_var','reco_cos_theta_pi','','',''),      
        ('slc','measure_var','pi_chi2_proton','','',''),
        ('slc','measure_var','pi_chi2_mu','','',''),
        ('slc','measure_var','pi_chi2_exp_pol','','',''),
        ('slc','measure_var','pi_scatter_angle_ratio','','',''),
        ('slc','measure_var','pi_max_daughter_hits','','',''),
        ('slc','measure_var','pi_frac50','','',''),
        ('slc','measure_var','pi_bdt_score_proton','','',''),
        ('slc','measure_var','pi_bdt_score_muon_pion','','',''),
    
        ('slc','measure_var','muon_contained','','',''),
        
        ('slc','measure_var','mu_true_pdg','','',''),
        ('slc','measure_var','mu_true_p_type','','',''),
        ('slc','measure_var','mu_true_end_process','','',''),
        ('slc','measure_var','pi_true_pdg','','',''),
        ('slc','measure_var','pi_true_p_type','','',''),
        ('slc','measure_var','pi_true_end_process','','',''),

        ('slc','measure_var','pi_true_p','','',''),
        ('slc','measure_var','mu_true_p','','',''),
        ('slc','measure_var','pi_true_costheta','','',''),
        ('slc','measure_var','mu_true_costheta','','',''),   

        #Cut vars
        ('slc','cut_var','n_prim_tracks','','',''), 
        ('slc','cut_var','n_prim_showers','','',''),
        ('slc','cut_var','n_MIP_candidates','','',''),
        ('slc','cut_var','n_MIP_candidates_proton','','',''),
        ('slc','cut_var','n_exiting_pfps','','',''),
        ('slc','cut_var','n_exiting_z_pfps','','',''),
        ('slc','cut_var','n_MIP_candidate_michel','','',''),
        ('slc','cut_var','n_extra_pions','','',''),
 
        ('slc','measure_var','delta_pT','','',''),
        ('slc','measure_var','delta_alpha_T','','',''),
        ('slc','measure_var','delta_phi_T','','',''),

        ('pfp', 'trk', 'dir', 'x', '', ''),
        ('pfp', 'trk', 'dir', 'y', '', ''),
        ('pfp', 'trk', 'dir', 'z', '', ''),
    
        ('pfp', 'trk', 'start', 'x', '', ''),
        ('pfp', 'trk', 'start', 'y', '', ''),
        ('pfp', 'trk', 'start', 'z', '', ''),

        ('pfp', 'trk', 'end', 'x', '', ''),
        ('pfp', 'trk', 'end', 'y', '', ''),
        ('pfp', 'trk', 'end', 'z', '', ''),
    ]
# Requires the module-level helpers defined above get_mu_pi_vars:
#   EXPECTED_COLS, _MU_PI_SCHEMA, _empty_mu_pi_vars, _enforce_schema
# (see transverse_and_schema.py)


# =====================================================================
# Schema helpers (module level, above get_mu_pi_vars)
# =====================================================================
_MU_PI_FLOAT_COLS = [
    'p_mu_x', 'p_mu_y', 'p_mu_z', 'p_pi_x', 'p_pi_y', 'p_pi_z',
    'reco_p_mu', 'cos_theta_mu',
    'mu_chi2_proton', 'mu_chi2_mu', 'mu_chi2_exp_pol', 'mu_scatter_angle_ratio',
    'mu_max_daughter_hits', 'mu_frac50', 'mu_bdt_score_proton', 'mu_bdt_score_muon_pion',
    'range_p_pi', 'TLE_p_pi', 'cos_theta_pi',
    'pi_chi2_proton', 'pi_chi2_mu', 'pi_chi2_exp_pol', 'pi_scatter_angle_ratio',
    'pi_max_daughter_hits', 'pi_frac50', 'pi_bdt_score_proton', 'pi_bdt_score_muon_pion',
    'pi_true_p', 'mu_true_p', 'pi_true_costheta', 'mu_true_costheta',
]
 
_MU_PI_DEFAULTS = {
    **{c: -999.0 for c in _MU_PI_FLOAT_COLS},
    'muon_contained': False,
    'mu_true_pdg': -1, 'mu_true_end_process': -1,
    'pi_true_pdg': -1, 'pi_true_end_process': -1,
    'mu_true_p_type': "none", 'pi_true_p_type': "none",
}
 
EXPECTED_COLS = list(_MU_PI_DEFAULTS)
 
_MU_PI_SCHEMA = {
    **{c: 'float32' for c in _MU_PI_FLOAT_COLS},
    'muon_contained': 'bool',
    'mu_true_pdg': 'int32', 'mu_true_end_process': 'int32',
    'pi_true_pdg': 'int32', 'pi_true_end_process': 'int32',
    'mu_true_p_type': 'object', 'pi_true_p_type': 'object',
}
 
 
def _empty_mu_pi_vars():
    return pd.Series(_MU_PI_DEFAULTS)

def _col_tokens(label):
    """Set of string tokens of a column label (works for tuple labels)."""
    if isinstance(label, tuple):
        return {str(t) for t in label}
    return {str(label)}


def _float_is_safe(s):
    """True if casting the float64 Series to float32 loses nothing important:
    no overflow, and no large integer-valued data (ids, counters) stored as float."""
    a = s.to_numpy()
    fin = a[np.isfinite(a)]
    if fin.size == 0:
        return True
    amax = np.abs(fin).max()
    if amax >= _F32_MAX:
        return False
    if amax >= _F32_EXACT_INT and np.all(fin == np.round(fin)):
        return False
    return True

 
def _enforce_schema(slcdf):
    for col, dtype in _MU_PI_SCHEMA.items():
        if dtype == 'bool':
            slcdf[col] = slcdf[col].fillna(False).astype(bool)
        elif dtype == 'int32':
            slcdf[col] = pd.to_numeric(slcdf[col], errors='coerce').fillna(-9999).astype('int32')
        elif dtype == 'object':
            slcdf[col] = slcdf[col].fillna("none").astype('object')
        else:
            slcdf[col] = pd.to_numeric(slcdf[col], errors='coerce').astype('float32')
    return slcdf

SHRINK_SKIP_COLS = {"run", "subrun", "evt", "event", "pot", "spill", "ntuple", "entry"}

_F32_MAX = 3.0e38
_F32_EXACT_INT = 2 ** 24      # float32 represents integers exactly only up to here


def _col_tokens(label):
    """Set of string tokens of a column label (works for tuple labels)."""
    if isinstance(label, tuple):
        return {str(t) for t in label}
    return {str(label)}

import gc
import pyarrow as pa
import pandas as pd



def shrink_df(df, floats=True, ints=True, skip_cols=None, verbose=False, **_ignored):
    """In place, column by column:
       float64 -> float32  (only if safe: no overflow, no large integer-valued ids stored as float)
       int64   -> int32    (only if all values fit in int32)
    Everything else (objects, bools, categoricals, already-small dtypes) is untouched.
    The index is never modified. Do NOT use on the hdr DataFrame (POT sums need float64)."""
    if df is None or df.empty:
        return df

    skip = set(SHRINK_SKIP_COLS) | set(skip_cols or ())
    before = df.memory_usage(deep=True).sum() if verbose else 0
    i32 = np.iinfo(np.int32)

    for c in list(df.columns):
        if _col_tokens(c) & skip:
            continue
        s = df[c]
        if not isinstance(s, pd.Series):          # duplicated column label -> leave alone
            continue
        dt = s.dtype

        if floats and dt == np.float64:
            if _float_is_safe(s):
                df[c] = s.astype(np.float32)
        elif (ints and pd.api.types.is_integer_dtype(dt) and dt != np.bool_
              and np.dtype(dt).itemsize > 4 and len(s) > 0):
            if s.min() >= i32.min and s.max() <= i32.max:
                df[c] = s.astype(np.int32)

    if verbose:
        after = df.memory_usage(deep=True).sum()
        print(f"    shrink_df: {before / 1e9:.3f} GB -> {after / 1e9:.3f} GB")
    return df

import gc
import pandas as pd


def make_cc1pi_finaldf(f, updatecalo=None, use_recalo = False):

    pandora_df = make_pandora_df(f, trkScoreCut=False, trkDistCut=-1, cutClearCosmic=True,
                                 requireFiducial=False, updatecalo=updatecalo)
    if pandora_df.empty:
        # Define the 3-level empty MultiIndex
        idx = pd.MultiIndex(
            levels=[[], [], []],
            codes=[[], [], []],
            names=['entry', 'rec.slc..index', 'rec.slc.reco.pfp..index']
        )
        return pd.DataFrame(index=idx, columns=pd.MultiIndex.from_tuples(cols))

    # Data vs MC: data has no truth matching. Computed on the whole file,
    # BEFORE any slice filtering.
    tm = pandora_df.slc.tmatch
    is_data = bool(tm.idx.eq(-999).all() and tm.eff.isna().all() and tm.pur.isna().all())

    cc1pi_shwbranches = [shwbranch + 'bestplane_energy']
    shw_df = loadbranches(f["recTree"], cc1pi_shwbranches)
    shw_df = shw_df.rec.slc.reco

    # create mcs df
    mcs_df = loadbranches(f["recTree"], [trkmcsbranches[1]]).rec.slc.reco.pfp.trk.mcsP

    hit0_df = make_trkhitdf_plane0(f, updatecalo=updatecalo)
    hit1_df = make_trkhitdf_plane1(f, updatecalo=updatecalo)
    hit2_df = make_trkhitdf_plane2(f, updatecalo=updatecalo)

    hit0_df = hit0_df.sort_values('rr', ascending=True)
    hit1_df = hit1_df.sort_values('rr', ascending=True)
    hit2_df = hit2_df.sort_values('rr', ascending=True)

    hit_dfs = [hit0_df, hit1_df, hit2_df]
    hit_names = ['nhit0', 'nhit1', 'nhit2']

    best_hit_df = get_best_hit_df(pandora_df, hit_dfs)
    pandora_df = add_nhit_column(pandora_df, hit_dfs, hit_names, 1000)
    pandora_df = multicol_merge(pandora_df, shw_df, left_index=True, right_index=True,
                                how="left", validate="one_to_one")

    pandora_df = pandora_df[CutMasks.is_obvious_cosmic_cut_mask(pandora_df)]
    pandora_df[('slc', 'vertex_inside_fv', '', '', '', '')] = CutMasks.InFV_strict(pandora_df)

    pandora_df = add_best_ke_column(pandora_df)
    pandora_df = add_p_type_column(pandora_df)
    pandora_df = add_best_chi2_columns(pandora_df, update_calo=updatecalo)
    pandora_df = add_max_daughter_hits_column(pandora_df)

    # Remove non primary particles (not useful for analysis)
    pandora_df = pandora_df[pandora_df.pfp.parent_is_primary == True]
    pandora_df = add_n_protons_column(pandora_df)

    fixed_hit_df = (
        best_hit_df
        .groupby(level=pfp_levels, group_keys=False)
        .apply(dEdxCleaning.get_fix_hit_df, plot=False)
    )
    fixed_hit_df['dE'] = fixed_hit_df.dedx * fixed_hit_df.pitch

    pandora_df = add_mean_dedx_column(pandora_df, fixed_hit_df)
    pandora_df = add_chi2_exp_pol_column(pandora_df, fixed_hit_df)
    pandora_df = add_frac50_column(pandora_df, fixed_hit_df)
    pandora_df = add_scatter_angle_ratio_column(pandora_df, mcs_df)
    pandora_df = add_max_angle_between_candidates_column(pandora_df)
    pandora_df[('slc', 'measure_var', 'angle_between_candidates', '', '', '')] = \
        pandora_df[('slc', 'measure_var', 'max_angle_between_candidates', '', '', '')]

    pandora_df[('pfp', 'is_exiting', '', '', '', '')] = CutMasks.exiting_pfp_mask(pandora_df)
    pandora_df[('pfp', 'is_exiting_z', '', '', '', '')] = CutMasks.exiting_z_pfp_mask(pandora_df)

    # add proton BDT column
    proton_BDT_input_columns = [
        ('pfp', 'trk', 'chi2pid', 'best', 'chi2_muon', ''),
        ('pfp', 'trk', 'chi2pid', 'best', 'chi2_proton', ''),
        ('pfp', 'trk', 'chi2_exp_pol', '', '', ''),
        ('pfp', 'trk', 'frac50', '', '', '')
    ]
    col_BDTG_score_proton = ('pfp', 'trk', 'bdt_proton_score', '', '', '')

    BDT_proton_path = "/../BDTs/bdtg_model_proton_gen1_og.pkl"
    if use_recalo:
        BDT_proton_path = "/../BDTs/bdtg_model_proton_gen1_recalo.pkl"

        
    pandora_df = add_bdt_score(
        pandora_df,
        base_running_path + BDT_proton_path,
        proton_BDT_input_columns,
        col_BDTG_score_proton
    )

    # add muon pion BDT column
    muon_pion_BDT_input_columns = [
        ('pfp', 'trk', 'chi2pid', 'best', 'chi2_muon', ''),
        ('pfp', 'trk', 'chi2pid', 'best', 'chi2_proton', ''),
        ('pfp', 'trk', 'chi2_exp_pol', '', '', ''),
        ('pfp', 'scatter_angle_ratio', '', '', '', ''),
        ('pfp', 'max_daughter_hits', '', '', '', '')
    ]

    
    BDT_muon_pion_path = "/../BDTs/bdtg_model_muon_pion_gen1_og.pkl"
    if use_recalo:
        BDT_muon_pion_path = "/../BDTs/bdtg_model_muon_pion_gen1_recalo.pkl"
        
    col_BDTG_score_muon_pion = ('pfp', 'trk', 'bdt_muon_pion_score', '', '', '')
    pandora_df = add_bdt_score(
        pandora_df,
        base_running_path + BDT_muon_pion_path,
        muon_pion_BDT_input_columns,
        col_BDTG_score_muon_pion
    )

    # Passing cuts or not boolean
    pandora_df[('slc', 'cut', 'obvious_cosmic', '', '', '')] = CutMasks.is_obvious_cosmic_cut_mask(pandora_df)
    pandora_df[('slc', 'cut', 't0', '', '', '')] = CutMasks.t0_cut_mask(pandora_df)
    pandora_df[('slc', 'cut', 'inside_FV', '', '', '')] = CutMasks.InFV_strict(pandora_df)
    pandora_df[('slc', 'cut', 'nu_score', '', '', '')] = CutMasks.nu_score_cut_mask(pandora_df)
    pandora_df[('slc', 'cut', 'track', '', '', '')] = CutMasks.track_cut_mask(pandora_df, group_levels)
    pandora_df[('slc', 'cut', 'shower', '', '', '')] = CutMasks.shower_cut_mask(pandora_df, group_levels)
    pandora_df[('slc', 'cut', 'MIP_candidates', '', '', '')] = CutMasks.chi2_cut_mask(pandora_df, group_levels)
    pandora_df[('slc', 'cut', 'angle', '', '', '')] = CutMasks.angle_cut_mask(pandora_df)
    pandora_df[('slc', 'cut', 'containment', '', '', '')] = CutMasks.containment_z_cut_mask(pandora_df, group_levels)
    pandora_df[('slc', 'cut', 'michel', '', '', '')] = CutMasks.michel_cut_mask(pandora_df, group_levels)
    pandora_df[('slc', 'cut', 'extra_pion', '', '', '')] = CutMasks.extra_pion_cut_mask(pandora_df, group_levels)
    pandora_df[('slc', 'cut', 'proton_BDT', '', '', '')] = CutMasks.proton_BDT_cut_mask(pandora_df, group_levels, use_recalo)
    pandora_df[('slc', 'cut', 'proton_BDT_sideband', '', '', '')] = CutMasks.proton_BDT_sideband_mask(pandora_df, group_levels)
    pandora_df[('slc', 'cut', 'proton_BDT_2pi', '', '', '')] = CutMasks.proton_BDT_cut_mask_2pi(pandora_df, group_levels)
    pandora_df[('slc', 'cut', 'TPC_containment', '', '', '')] = CutMasks.TPC_containment_mask(pandora_df, group_levels)
    pandora_df[('slc', 'cut', 'no_high_yz', '', '', '')] = CutMasks.not_in_high_y_high_z_containment_mask(pandora_df, group_levels)

    pandora_df = add_n_primary_tracks_column(pandora_df)
    pandora_df = add_n_primary_showers_column(pandora_df)
    pandora_df = add_n_primary_MIP_column(pandora_df)
    pandora_df = add_n_primary_MIP_proton_column(pandora_df)
    pandora_df = add_n_exiting_pfps_column(pandora_df)
    pandora_df = add_n_exiting_z_pfps_column(pandora_df)
    pandora_df = add_n_MIP_candidate_michel_column(pandora_df)
    pandora_df = add_n_extra_pions_column(pandora_df)

    # Candidate slices for the muon/pion kinematics
    candidate_mask = (pandora_df.slc.cut.inside_FV & pandora_df.slc.cut.t0
                      & CutMasks.long_primary_track_cut_mask(pandora_df, group_levels)
                      & (pandora_df.slc.cut_var.n_MIP_candidates >= 1)
                      & pandora_df.slc.cut.containment)
    candidate_df = pandora_df[candidate_mask]

    muon_pion_vars = (
        candidate_df
        .groupby(level=group_levels, group_keys=False)
        .apply(get_mu_pi_vars, best_hit_df=best_hit_df, use_data=is_data)
    )
    

    # Build slcdf aligned to all slices in pandora_df
    # (every slice gets a row, even if missing in muon_pion_vars)
    all_slices_idx = pandora_df.index.droplevel('rec.slc.reco.pfp..index').unique()

    if muon_pion_vars.empty:
        slcdf = pd.DataFrame(index=all_slices_idx, columns=EXPECTED_COLS)
    else:
        if isinstance(muon_pion_vars, pd.Series):
            # A single group can come back as a Series
            if not isinstance(muon_pion_vars.index, pd.MultiIndex):
                muon_pion_vars = muon_pion_vars.to_frame().T
            else:
                muon_pion_vars = muon_pion_vars.unstack()

        slcdf = muon_pion_vars.reindex(columns=EXPECTED_COLS)
        slcdf = slcdf.reindex(all_slices_idx)

    # Enforce schema and dtypes (replaces the old for-loop over `schema`)
    slcdf = _enforce_schema(slcdf)

    slcdf = slcdf.rename(columns={
        'cos_theta_mu': 'reco_cos_theta_mu',
        'cos_theta_pi': 'reco_cos_theta_pi'
    })

    slcdf.columns = pd.MultiIndex.from_tuples(
        [('slc', 'measure_var', col, '', '', '') for col in slcdf.columns]
    )
    

    # Join back
    pandora_df = pandora_df.join(slcdf, on=['entry', 'rec.slc..index'])

    pandora_df = add_transverse_vars_column(pandora_df)

    mv = pandora_df.slc.measure_var
    pandora_df[('slc', 'cut', 'energy', '', '', '')] = (
        (mv.reco_p_mu > 0.1) & (mv.reco_p_mu < 1)
        & (mv.TLE_p_pi > 0.13) & (mv.TLE_p_pi < 2)
    )

    min_df = shrink_df(pandora_df[cols].copy())
    min_df = min_df[min_df.pfp.trk.len > 0]
    return min_df
    


def make_cc1pi_final_df_slim(f):
    pandora_df = make_cc1pi_finaldf(f)

    cols_to_keep_reco = [
        ('slc', 'self', '', '', '', ''),
        ('slc', 'tmatch', 'idx', '', '', ''),
        ('slc', 'nu_score', '', '', '', ''),
        ('slc', 'cut', 'obvious_cosmic', '', '', ''),
        ('slc', 'cut', 't0', '', '', ''),
        ('slc', 'cut', 'inside_FV', '', '', ''),
        ('slc', 'cut', 'nu_score', '', '', ''),
        ('slc', 'cut', 'track', '', '', ''),
        ('slc', 'cut', 'shower', '', '', ''),
        ('slc', 'cut', 'MIP_candidates', '', '', ''),
        ('slc', 'cut', 'angle', '', '', ''),
        ('slc', 'cut', 'proton_BDT', '', '', ''),
        ('slc', 'cut', 'proton_BDT_sideband', '', '', ''),
        ('slc', 'cut', 'containment', '', '', ''),
        ('slc', 'cut', 'TPC_containment', '', '', ''),
        ('slc','cut','no_high_yz','','',''),
        ('slc', 'cut', 'michel', '', '', ''),
        ('slc', 'cut', 'extra_pion', '', '', ''),
        ('slc','cut','energy','','',''),
        ('slc', 'measure_var', 'angle_between_candidates', '', '', ''),
        ('slc', 'measure_var', 'num_protons', '', '', ''),
        ('slc','measure_var','reco_p_mu','','',''),
        ('slc','measure_var','reco_cos_theta_mu','','',''),
        ('slc','measure_var','TLE_p_pi','','',''),
        ('slc','measure_var','reco_cos_theta_pi','','',''), 
        ('slc','cut_var','n_MIP_candidates','','',''),
        ('slc','measure_var','delta_pT','','',''),
        ('slc','measure_var','delta_alpha_T','','',''),
        ('slc','measure_var','delta_phi_T','','','')
    ]
    pandora_df = pandora_df[cols_to_keep_reco]
    
    slc_df = (
            pandora_df
            .groupby(['entry', 'rec.slc..index'])
            .first()
        )
    slc_df = shrink_df(slc_df.sort_index())
    
    return slc_df


    
def make_cc1pi_final_df_recalo_update_cv(f):
    return make_cc1pi_finaldf(f, updatecalo = "cv", use_recalo = True)
    

def make_cc1pi_final_df_recalo_ccal_p(f):
    return make_cc1pi_finaldf(f, updatecalo = "ccal_p")

def make_cc1pi_final_df_recalo_ccal_m(f):
    return make_cc1pi_finaldf(f, updatecalo = "ccal_m")

def make_cc1pi_final_df_recalo_alpha_p(f):
    return make_cc1pi_finaldf(f, updatecalo = "alpha_p")

def make_cc1pi_final_df_recalo_alpha_m(f):
    return make_cc1pi_finaldf(f, updatecalo = "alpha_m")

def make_cc1pi_final_df_recalo_beta_p(f):
    return make_cc1pi_finaldf(f, updatecalo = "beta_p")

def make_cc1pi_final_df_recalo_beta_m(f):
    return make_cc1pi_finaldf(f, updatecalo = "beta_m")

def make_cc1pi_final_df_recalo_r_p(f):
    return make_cc1pi_finaldf(f, updatecalo = "R_p")

def make_cc1pi_final_df_recalo_r_m(f):
    return make_cc1pi_finaldf(f, updatecalo = "R_m")

def make_trkhitdf_plane2_ccal_p(f):
    return make_trkhitdf_plane2(f, updatecalo = "ccal_p")

def make_trkhitdf_plane2_ccal_m(f):
    return make_trkhitdf_plane2(f, updatecalo = "ccal_m")