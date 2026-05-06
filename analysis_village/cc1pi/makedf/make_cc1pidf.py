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
from ROOT import TMVA
import os
import ROOT
import subprocess
from ROOT import std
import pickle
from itertools import combinations


# 1. Setup Paths
base_running_path = os.path.dirname(os.path.abspath(__file__))
base_path = base_running_path + "/../TLExtensionMethod"

try:
    root_lib_dir = subprocess.check_output(['root-config', '--libdir'], text=True).strip()
except:
    root_lib_dir = "/cvmfs/larsoft.opensciencegrid.org/products/root/v6_28_12/Linux64bit+3.10-2.17-e26-p3915-prof/lib"

# 2. Update Environment
os.environ['LD_LIBRARY_PATH'] = f"{root_lib_dir}:{os.environ.get('LD_LIBRARY_PATH', '')}"
ROOT.gSystem.AddDynamicPath(root_lib_dir)

# 3. Add Include Path so classes can find each other's headers
ROOT.gInterpreter.AddIncludePath(base_path)

# 4. Explicit Compilation and Loading Function
def load_custom_class(class_name):
    source_file = os.path.join(base_path, f"{class_name}.cpp")
    header_file = os.path.join(base_path, f"{class_name}.h")
    so_file = os.path.join(base_path, f"{class_name}_cpp.so")
    
    # Compile with 'k' (keep), 'O' (optimize), 'f' (force)
    # Force is useful to ensure the symbol table is rebuilt correctly
    status = ROOT.gSystem.CompileMacro(source_file, "kOf")
    
    if status >= 0:
        # CRITICAL: Manually load the shared library to resolve symbols
        if ROOT.gSystem.Load(so_file) < 0:
            print(f"⚠️  Compiled {class_name} but failed to load {so_file}")
            return False
            
        # Declare the header to the interpreter to help dictionary lookup
        ROOT.gInterpreter.Declare(f'#include "{header_file}"')
        print(f"📦 Successfully compiled and linked: {class_name}")
        return True
    else:
        print(f"❌ Failed to compile: {class_name}")
        return False

# 5. Execute Sequence
# We must load PhysdEdx first as Hypfit depends on it
if load_custom_class("PhysdEdx"):
    if load_custom_class("Hypfit"):
        try:
            # Re-declaring just to be safe before instantiation
            ROOT.gInterpreter.ProcessLine(f'#include "{os.path.join(base_path, "Hypfit.h")}"')
            
            # Instantiate
            h_fit = ROOT.Hypfit()
            print("🚀 Success! Hypfit object initialized and linked.")
        except Exception as e:
            print(f"❌ Error during instantiation: {e}")
            # Final fallback: access via C++ global pointer if Python attribute fails
            ROOT.gInterpreter.ProcessLine("Hypfit* h_fit_ptr = new Hypfit();")
            h_fit = ROOT.h_fit_ptr
            print("🚀 Success! Hypfit object initialized via Global Pointer fallback.")
            
# Silence ROOT
#ROOT.gSystem.RedirectOutput(os.devnull, "w")
#ROOT.gROOT.SetBatch(True)

group_levels = ['entry', 'rec.slc..index']
pfp_levels = ['entry','rec.slc..index','rec.slc.reco.pfp..index']

def TruthInFV(data):
    xmax = 190.
    zmin = 10.
    zmax = 450.
    ymax_highz = 100.
    pass_xz = (np.abs(data.x) < xmax) & (data.z > zmin) & (data.z < zmax)
    pass_y = ((data.z < 250) & (np.abs(data.y) < 190.)) | ((data.z > 250) & (data.y > -190.) & (data.y < ymax_highz))
    return pass_xz & pass_y
    
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
    
    #return is_1pi1mu & is_NpiNmuNnNp & is_theta & is_mu_contained
    return is_1pi1mu & is_NpiNmuNnNp & is_theta

def add_max_angle_between_candidates_column(df):
    # 1. Initialize result series at the slice level
    # We use nan so that slices with < 2 particles naturally fail the cut
    angle_series = pd.Series(np.nan, index=df.index)
    
    # 2. Filter for MIP candidates only
    mip_df = df[CutMasks.is_MIP_candidate_mask(df)]
    
    for idx, group in mip_df.groupby(level=group_levels):
        # We need at least 2 particles to have an angle
        if len(group) < 2:
            continue

        # Extract direction vectors for ALL candidates in the slice
        # Shape: (N_particles, 3)
        dirs = group.loc[:, [
            ('pfp', 'trk', 'dir', 'x', '', ''),
            ('pfp', 'trk', 'dir', 'y', '', ''),
            ('pfp', 'trk', 'dir', 'z', '', '')
        ]].astype(float).values

        max_angle = 0.0
        
        # 4. Check every unique pair of particles
        for dir_1, dir_2 in combinations(dirs, 2):
            mag1 = np.linalg.norm(dir_1)
            mag2 = np.linalg.norm(dir_2)
            
            if mag1 > 0 and mag2 > 0:
                dot = np.dot(dir_1, dir_2)
                # Compute angle and update max if this pair is wider
                angle = np.arccos(np.clip(dot / (mag1 * mag2), -1.0, 1.0))
                if angle > max_angle:
                    max_angle = angle

        # 5. Broadcast back to the full dataframe for all PFP entries in this slice
        # This ensures every particle in the slice 'knows' the max angle of the group
        angle_series.loc[idx] = max_angle
            
    # Add column to the original dataframe
    df[('slc', 'measure_var', 'max_angle_between_candidates', '', '', '')] = angle_series

    return df

def add_min_angle_between_candidates_column(df):
    # 1. Initialize result series at the slice level
    # We use nan so that slices with < 2 particles naturally fail the cut
    angle_series = pd.Series(np.nan, index=df.index)
    
    # 2. Filter for MIP candidates only
    mip_df = df[CutMasks.is_MIP_candidate_mask(df)]
    
    for idx, group in mip_df.groupby(level=group_levels):
        # We need at least 2 particles to have an angle
        if len(group) < 2:
            continue

        # Extract direction vectors for ALL candidates in the slice
        # Shape: (N_particles, 3)
        dirs = group.loc[:, [
            ('pfp', 'trk', 'dir', 'x', '', ''),
            ('pfp', 'trk', 'dir', 'y', '', ''),
            ('pfp', 'trk', 'dir', 'z', '', '')
        ]].astype(float).values

        max_angle = 3000.0
        
        # 4. Check every unique pair of particles
        for dir_1, dir_2 in combinations(dirs, 2):
            mag1 = np.linalg.norm(dir_1)
            mag2 = np.linalg.norm(dir_2)
            
            if mag1 > 0 and mag2 > 0:
                dot = np.dot(dir_1, dir_2)
                # Compute angle and update max if this pair is wider
                angle = np.arccos(np.clip(dot / (mag1 * mag2), -1.0, 1.0))
                if angle < max_angle:
                    max_angle = angle

        # 5. Broadcast back to the full dataframe for all PFP entries in this slice
        # This ensures every particle in the slice 'knows' the max angle of the group
        angle_series.loc[idx] = max_angle
            
    # Add column to the original dataframe
    df[('slc', 'measure_var', 'min_angle_between_candidates', '', '', '')] = angle_series

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
    
def add_nu_categ_column(df, is_truth_df = False):
    if(is_truth_df):
        truth_df = df
    else:
        truth_df = df.slc.truth # Make a copy to safely assign

    is_inside_fv = TruthInFV(truth_df.position)
    is_nu = IsNu(truth_df)
    is_signal = isCC1Pi(truth_df)
    is_cc = truth_df.iscc
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


def add_n_true_proton_column(df): 
    df[('true_var', 'num_protons','')] = df[('np_P_325MeV_10000MeV','','')]
    df.loc[df.true_var.num_protons >= 2,('true_var', 'num_protons','')] = 2
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

def add_genie_categ_column(df, is_truth_df = False):
    if(is_truth_df):
        truth_df = df
    else:
        truth_df = df.truth # Make a copy to safely assign

    is_inside_fv = TruthInFV(truth_df.position)
    is_nu = IsNu(truth_df)
    is_cc = truth_df.iscc.astype(bool)
    is_nu_mu_cc = is_cc & (abs(truth_df.pdg) == 14)

    genie_categ = pd.Series("other", index=truth_df.index, dtype="object")
    # Apply categories
    genie_categ[~is_nu] = "cosmic"
    genie_categ[is_nu & ~is_inside_fv] = "out_AV_nu"
    genie_categ[is_nu & is_inside_fv & ~is_cc & (abs(truth_df.pdg) == 14)] = "nu_mu_NC" 
    genie_categ[is_nu & is_inside_fv & is_nu_mu_cc & (df.genie_mode == 0)] = "nu_mu_CC_QE" 
    genie_categ[is_nu & is_inside_fv & is_nu_mu_cc & (df.genie_mode == 10)] = "nu_mu_CC_MEC" 
    genie_categ[is_nu & is_inside_fv & is_nu_mu_cc & (df.genie_mode == 1)] = "nu_mu_CC_Res" 
    genie_categ[is_nu & is_inside_fv & is_nu_mu_cc & (df.genie_mode == 2)] = "nu_mu_CC_Dis" 
    
    if(is_truth_df):
        df['genie_categ'] = genie_categ
    else:
        df[('truth', 'genie_categ', '', '','','')] = nu_categ 
    return df

def add_best_chi2_columns(df, update_calo = None):
    chi2_vars = [
        'chi2_muon',
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

def add_nhit_column(
    slc_df,
    hit_dfs,
    hit_names,
    dedx_max=1000
):

    if len(hit_dfs) != len(hit_names):
        raise ValueError("hit_dfs and hit_names must have the same length")

    for df, name in zip(hit_dfs, hit_names):
        out_col = ('pfp', name, '', '', '', '')

        # Count hits per particle
        hit_counts = (
            df[df.dedx < dedx_max]
            .groupby(level=[
                'entry',
                'rec.slc..index',
                'rec.slc.reco.pfp..index'
            ])
            .size()
            .rename(out_col)
        )

        # Map counts into slc_df (0 if no hits)
        slc_df.loc[:, out_col] = hit_counts.reindex(
            slc_df.index,
            fill_value=0
        )

    # Total nhit column
    nhit_cols = [('pfp', name, '', '', '', '') for name in hit_names]
    total_col = ('pfp', 'nhit', '', '', '', '')

    slc_df.loc[:, total_col] = slc_df[nhit_cols[0]] + slc_df[nhit_cols[1]] + slc_df[nhit_cols[2]]

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
        df.index = df.index.set_names(
            'rec.slc.reco.pfp.trk.calo.points..index',
            level=-1
        )
    
    for df in hit_dfs:
        level_name = f"rec.slc.reco.pfp.trk.calo.points..index"
        mask = df.index.droplevel(level_name).isin(slc_df.index)
        df = df[mask]
  
    # Unpack back if you want the same variable names
    best_plane = slc_df.pfp.trk.bestplane

    best_hit_df = (
        hit_dfs[0].where(best_plane == 0)
        .combine_first(hit_dfs[1].where(best_plane == 1))
        .combine_first(hit_dfs[2].where(best_plane == 2))
    )
    
    best_hit_df = best_hit_df.sort_values('rr', ascending=True)
    return best_hit_df

def get_reduced_hit_df(hit_df):
    reduced_best_hit_df = hit_df[hit_df.rr < 30].copy()
    return reduced_best_hit_df

def compute_frac50(group):
    total_E = group['dE'].sum()
    max_rr = group['rr'].max()
    cumE = group['dE'].cumsum()

    # find first point where cumE >= 50% total
    mask = cumE >= 0.5 * total_E
    if not mask.any():
        return 1

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

def add_n_primary_tracks_column(df):
    track_df = df[CutMasks.is_primary_track_mask(df)]
    track_counts = track_df.groupby(level=group_levels).size()

    target_key = ('slc','cut_var','n_prim_tracks','','','')
    
    # Expand counts to full MultiIndex (broadcast to pfp level)
    df.loc[:,target_key] = df.index.droplevel('rec.slc.reco.pfp..index').map(track_counts)
    
    # Replace NaN (slices with 0 tracks) with 0
    df.loc[:,target_key] = df[target_key].fillna(0)
    
    return df

def add_n_primary_showers_column(df):
    is_pandora_primary_mask = (df.pfp.parent_is_primary == True)
    is_shower_mask = (df.pfp.trackScore >= 0) & (df.pfp.trackScore < CTE.max_shower_track_score)
    energy_mask = (df.pfp.shw.bestplane_energy > CTE.min_shower_ke)
    shower_df = df[is_pandora_primary_mask & is_shower_mask & energy_mask] 
        
    # Count how many pfps per slice
    shower_counts = shower_df.groupby(level=group_levels).size()

    target_key = ('slc','cut_var','n_prim_showers','','','')
    
    # Expand counts to full MultiIndex (broadcast to pfp level)
    df.loc[:,target_key] = df.index.droplevel('rec.slc.reco.pfp..index').map(shower_counts)
    
    # Replace NaN (slices with 0 tracks) with 0
    df.loc[:,target_key] = df[target_key].fillna(0)
    
    return df
    
def add_n_primary_MIP_column(df):
    MIP_df = df[CutMasks.is_MIP_candidate_mask(df)]
    MIP_counts = MIP_df.groupby(level=group_levels).size()

    target_key = ('slc','cut_var','n_MIP_candidates','','','')
    
    # Expand counts to full MultiIndex (broadcast to pfp level)
    df.loc[:,target_key] = df.index.droplevel('rec.slc.reco.pfp..index').map(MIP_counts)
    
    # Replace NaN (slices with 0 tracks) with 0
    df.loc[:,target_key] = df[target_key].fillna(0)
    
    return df

def add_n_primary_MIP_proton_column(df):
    BDT_proton_df = df[(CutMasks.is_MIP_candidate_mask(df)) & (df.pfp.trk.bdt_proton_score < CTE.BDT_proton_max_score)]
    BDT_proton_counts = BDT_proton_df.groupby(level=group_levels).size()

    target_key = ('slc','cut_var','n_MIP_candidates_proton','','','')
    
    # Expand counts to full MultiIndex (broadcast to pfp level)
    df.loc[:,target_key] = df.index.droplevel('rec.slc.reco.pfp..index').map(BDT_proton_counts)
    
    # Replace NaN (slices with 0 tracks) with 0
    df.loc[:,target_key] = df[target_key].fillna(0)
    
    return df

def add_n_exiting_pfps_column(df):
    exiting_pfp_df = df[(CutMasks.exiting_pfp_mask(df))]
    exiting_pfp_counts = exiting_pfp_df.groupby(level=group_levels).size()

    target_key = ('slc','cut_var','n_exiting_pfps','','','')
    
    # Expand counts to full MultiIndex (broadcast to pfp level)
    df.loc[:,target_key] = df.index.droplevel('rec.slc.reco.pfp..index').map(exiting_pfp_counts)
    
    # Replace NaN (slices with 0 tracks) with 0
    df.loc[:,target_key] = df[target_key].fillna(0)
    
    return df


def add_n_exiting_z_pfps_column(df):
    exiting_z_pfp_df = df[(CutMasks.exiting_z_pfp_mask(df))]
    exiting_z_pfp_counts = exiting_z_pfp_df.groupby(level=group_levels).size()

    target_key = ('slc','cut_var','n_exiting_z_pfps','','','')
    
    # Expand counts to full MultiIndex (broadcast to pfp level)
    df.loc[:,target_key] = df.index.droplevel('rec.slc.reco.pfp..index').map(exiting_z_pfp_counts)
    
    # Replace NaN (slices with 0 tracks) with 0
    df.loc[:,target_key] = df[target_key].fillna(0)
    
    return df

def add_n_MIP_candidate_michel_column(df):
    michel_df = df[(CutMasks.is_MIP_candidate_mask(df)) 
        & (df.pfp.trackScore < CTE.michel_max_track_score)
        & (df.pfp.trk.calo.best.ke < CTE.michel_max_visible_energy)
        & (df.pfp.max_daughter_hits == 0)]
    michel_pfp_counts = michel_df.groupby(level=group_levels).size()

    target_key = ('slc','cut_var','n_MIP_candidate_michel','','','')
    
    # Expand counts to full MultiIndex (broadcast to pfp level)
    df.loc[:,target_key] = df.index.droplevel('rec.slc.reco.pfp..index').map(michel_pfp_counts)
    
    # Replace NaN (slices with 0 tracks) with 0
    df.loc[:,target_key] = df[target_key].fillna(0)
    
    return df


def add_n_extra_pions_column(df):
    chi2_mask = (df.pfp.trk.chi2pid.best.chi2_muon < CTE.MIP_candidate_max_muon_score) & (df.pfp.trk.chi2pid.best.chi2_proton > CTE.MIP_candidate_min_proton_score)
    relaxed_MIP_df = df[(CutMasks.is_primary_track_mask(df)) & (chi2_mask)]
    
    # Count how many pfps per slice
    candidate_counts = relaxed_MIP_df.groupby(level=group_levels).size() - 2
 
    target_key = ('slc','cut_var','n_extra_pions','','','')
    
    # Expand counts to full MultiIndex (broadcast to pfp level)
    df.loc[:,target_key] = df.index.droplevel('rec.slc.reco.pfp..index').map(candidate_counts)
    
    # Replace NaN (slices with 0 tracks) with 0
    df.loc[:,target_key] = df[target_key].fillna(0)
    
    return df
  
def calculate_hypfit_p(track_id, best_plane, hit_df, target_pdg=211, cleaning="none"):
    try:
        # 1. Extract and sort hits
        track_hits = hit_df.loc[track_id]
        # 2. Strict type casting for C++ 
        # Fix: Handle cases where best_plane is a Series or a scalar
        if hasattr(best_plane, "iloc"):
            c_plane = int(best_plane.iloc[0])
        else:
            c_plane = int(best_plane)
            
        c_pdg = int(target_pdg)
        c_cleaning = str(cleaning)

        # Convert to contiguous float64 numpy arrays
        rr    = np.ascontiguousarray(track_hits['rr'].to_numpy(dtype=np.float64))
        dedx  = np.ascontiguousarray(track_hits['dedx'].to_numpy(dtype=np.float64))
        pitch = np.ascontiguousarray(track_hits['pitch'].to_numpy(dtype=np.float64))

        # Build mask: keep only entries where ALL values are finite
        mask = np.isfinite(rr) & np.isfinite(dedx) & np.isfinite(pitch)
        
        # Apply mask consistently to all arrays
        rr    = rr[mask]
        dedx  = dedx[mask]
        pitch = pitch[mask]
        
        # Now create ROOT vectors
        v_rr    = ROOT.std.vector('double')(rr)
        v_dedx  = ROOT.std.vector('double')(dedx)
        v_pitch = ROOT.std.vector('double')(pitch)

        reco_p = -10000
        reco_p = h_fit.GetTLExtensionP(c_pdg, v_rr, v_dedx, v_pitch, c_plane, c_cleaning)
        # Explicitly clear/delete to help PyROOT
        if 'v_rr' in locals():
            v_rr.clear()
            v_dedx.clear()
            v_pitch.clear()
            # Explicitly delete the Python reference to the C++ object
            del v_rr, v_dedx, v_pitch
            
        return reco_p
        
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

'''
def get_mu_pi_vars(group):
    # group has 2 rows (2 PFPs)

    exiting = CutMasks.exiting_pfp_mask(group)

    # Case 1️⃣: exactly one exiting PFP
    if exiting.sum() == 1:
        p_mu = group.loc[exiting].pfp.trk.mcsP.fwdP_muon.iloc[0]
        p_pi = group.loc[~exiting].pfp.trk.rangeP.p_pion.iloc[0]
        cos_theta_mu =  group.loc[exiting].pfp.trk.dir.z.iloc[0]
        cos_theta_pi =  group.loc[~exiting].pfp.trk.dir.z.iloc[0]
        muon_contained = False

    # Case 2️⃣: no exiting PFP
    else:
        group_sorted = group.sort_values(
            ('pfp','trk','bdt_muon_pion_score','','','')
        )

        p_mu = group_sorted.pfp.trk.rangeP.p_muon.iloc[1]
        p_pi = group_sorted.pfp.trk.rangeP.p_pion.iloc[0]
        cos_theta_mu =  group_sorted.pfp.trk.dir.z.iloc[1]
        cos_theta_pi =  group_sorted.pfp.trk.dir.z.iloc[0]
        muon_contained = True

    return pd.Series({
        'reco_p_mu': p_mu,
        'reco_p_pi': p_pi,
        'cos_theta_mu': cos_theta_mu,
        'cos_theta_pi': cos_theta_pi,
        'muon_contained': muon_contained,
    })
'''


def build_p3d(df, p_col):
    return np.vstack([
        df[p_col] * df[('pfp','trk','dir','x','','')],
        df[p_col] * df[('pfp','trk','dir','y','','')],
        df[p_col] * df[('pfp','trk','dir','z','','')]
    ]).T


def add_transverse_vars_column(df):
    
    # --- Select primary MIP protons ---
    proton_mask = ~CutMasks.is_MIP_candidate_mask(df) & CutMasks.is_primary_track_mask(df)
    proton_df = df[proton_mask].copy()

    # --- Build 3D proton vectors ---
    p_col = ('pfp','trk','rangeP','p_proton','','')
    proton_df[('slc','measure_var','p_px','','','')] = proton_df[p_col] * proton_df[('pfp','trk','dir','x','','')]
    proton_df[('slc','measure_var','p_py','','','')] = proton_df[p_col] * proton_df[('pfp','trk','dir','y','','')]
    proton_df[('slc','measure_var','p_pz','','','')] = proton_df[p_col] * proton_df[('pfp','trk','dir','z','','')]

    # --- Sum proton components per slice ---
    proton_sum = proton_df.groupby(group_levels)[[('slc','measure_var','p_px','','',''),('slc','measure_var','p_py','','',''),('slc','measure_var','p_pz','','','')]].sum()
    slice_index = df.index.droplevel('rec.slc.reco.pfp..index').unique()
    proton_sum = proton_sum.reindex(slice_index, fill_value=0)

    df[('slc','p_proton_total_x','','','','')] = proton_sum[('slc','measure_var','p_px','','','')].reindex(df.index.droplevel('rec.slc.reco.pfp..index')).values
    df[('slc','p_proton_total_y','','','','')] = proton_sum[('slc','measure_var','p_py','','','')].reindex(df.index.droplevel('rec.slc.reco.pfp..index')).values
    df[('slc','p_proton_total_z','','','','')] = proton_sum[('slc','measure_var','p_pz','','','')].reindex(df.index.droplevel('rec.slc.reco.pfp..index')).values

    # --- Sum proton + pion component-wise for hadron ---

    df[('slc','p_hadron_x','','','','')] = df[('slc','measure_var','p_pi_x','','','')] + df[('slc','p_proton_total_x','','','','')]
    df[('slc','p_hadron_y','','','','')] = df[('slc','measure_var','p_pi_y','','','')] + df[('slc','p_proton_total_y','','','','')]
    df[('slc','p_hadron_z','','','','')] = df[('slc','measure_var','p_pi_z','','','')] + df[('slc','p_proton_total_z','','','','')]
    
    # --- Build 3D arrays for muon and hadron ---
    p_mu_2d = np.stack([
        df[('slc','measure_var','p_mu_x','','','')],
        df[('slc','measure_var','p_mu_y','','','')],
        np.zeros(len(df))
    ], axis=1)
    
    p_had_2d = np.stack([
        df[('slc','p_hadron_x','','','','')],
        df[('slc','p_hadron_y','','','','')],
        np.zeros(len(df))
    ], axis=1)

  
    
    # --- Delta PT vector and magnitude ---
    delta_pT_vec = p_mu_2d + p_had_2d
    df[('slc','measure_var','delta_pT','','','')] = np.linalg.norm(delta_pT_vec, axis=1)

    # --- Delta alphaT ---
    dot_alpha = np.einsum('ij,ij->i', -p_mu_2d, delta_pT_vec)
    norm_alpha = np.linalg.norm(p_mu_2d, axis=1) * np.linalg.norm(delta_pT_vec, axis=1)
    safe = norm_alpha > 0
    delta_alpha_T = np.full(len(df), -999.0)
    delta_alpha_T[safe] = np.arccos(np.clip(dot_alpha[safe] / norm_alpha[safe], -1, 1))
    df[('slc','measure_var','delta_alpha_T','','','')] = delta_alpha_T

    # --- Delta phiT ---
    dot_phi = np.einsum('ij,ij->i', -p_mu_2d, p_had_2d)
    norm_phi = np.linalg.norm(p_mu_2d, axis=1) * np.linalg.norm(p_had_2d, axis=1)
    safe_phi = norm_phi > 0
    delta_phi_T = np.full(len(df), -999.0)
    delta_phi_T[safe_phi] = np.arccos(np.clip(dot_phi[safe_phi] / norm_phi[safe_phi], -1, 1))
    df[('slc','measure_var','delta_phi_T','','','')] = delta_phi_T

    return df


'''
def get_mu_pi_vars(group, best_hit_df):

    # group has 2 rows (2 PFPs)
    exiting = CutMasks.exiting_pfp_mask(group)

    # Case 1: exactly one exiting PFP
    if exiting.sum() == 1:
        muon_row = group.loc[exiting]
        pion_row = group.loc[~exiting]
        
        
        # NEW: Calculate p_pi using Hypfit
        
        p_mu = muon_row.pfp.trk.mcsP.fwdP_muon.iloc[0]
        cos_theta_mu = muon_row.pfp.trk.dir.z.iloc[0]
        mu_chi2_proton = muon_row.pfp.trk.chi2pid.best.chi2_proton.iloc[0]
        mu_chi2_mu = muon_row.pfp.trk.chi2pid.best.chi2_muon.iloc[0]
        mu_chi2_exp_pol = muon_row.pfp.trk.chi2_exp_pol.iloc[0]
        mu_scatter_angle_ratio = muon_row.pfp.scatter_angle_ratio.iloc[0]
        mu_max_daughter_hits = muon_row.pfp.max_daughter_hits.iloc[0]
        mu_frac50 = muon_row.pfp.trk.frac50.iloc[0]
        mu_bdt_score_proton = muon_row.pfp.trk.bdt_proton_score.iloc[0]
        mu_bdt_score_muon_pion = muon_row.pfp.trk.bdt_muon_pion_score.iloc[0]
        
        # split components
        p_mu_x = p_mu * muon_row.pfp.trk.dir.x.iloc[0]
        p_mu_y = p_mu * muon_row.pfp.trk.dir.y.iloc[0]
        p_mu_z = p_mu * muon_row.pfp.trk.dir.z.iloc[0]

        
        
        p_pi_range = pion_row.pfp.trk.rangeP.p_pion.iloc[0]
        cos_theta_pi = pion_row.pfp.trk.dir.z.iloc[0]
        pion_id = pion_row.index[0] # This is the (ntuple, entry, slc, pfp) tuple
        best_plane = pion_row.pfp.trk.bestplane.iloc[0]
        #best_plane = 2
        p_pi_TLE = calculate_hypfit_p(track_id=pion_id, best_plane=best_plane, hit_df=best_hit_df, target_pdg=211, cleaning="all")
        #p_pi_TLE = 0.5
        pi_chi2_proton = pion_row.pfp.trk.chi2pid.best.chi2_proton.iloc[0]
        pi_chi2_mu = pion_row.pfp.trk.chi2pid.best.chi2_muon.iloc[0]
        pi_chi2_exp_pol = pion_row.pfp.trk.chi2_exp_pol.iloc[0]
        pi_scatter_angle_ratio = pion_row.pfp.scatter_angle_ratio.iloc[0]
        pi_max_daughter_hits = pion_row.pfp.max_daughter_hits.iloc[0]
        pi_frac50 = pion_row.pfp.trk.frac50.iloc[0]
        pi_bdt_score_proton = pion_row.pfp.trk.bdt_proton_score.iloc[0]
        pi_bdt_score_muon_pion = pion_row.pfp.trk.bdt_muon_pion_score.iloc[0]
       
        # split components
        p_pi_x = p_pi_TLE * pion_row.pfp.trk.dir.x.iloc[0]
        p_pi_y = p_pi_TLE * pion_row.pfp.trk.dir.y.iloc[0]
        p_pi_z = p_pi_TLE * pion_row.pfp.trk.dir.z.iloc[0]


        
        muon_contained = False
        
        if(muon_row.pfp.trk.truth.p.pdg.iloc[0] > -2147483648):
            mu_true_pdg = muon_row.pfp.trk.truth.p.pdg.iloc[0]   
            mu_true_p_type = muon_row.pfp.trk.truth.p.p_type.iloc[0]
            mu_true_end_process = muon_row.pfp.trk.truth.p.end_process.iloc[0]
            pi_true_pdg = pion_row.pfp.trk.truth.p.pdg.iloc[0]
            pi_true_p_type = pion_row.pfp.trk.truth.p.p_type.iloc[0]
            pi_true_end_process = pion_row.pfp.trk.truth.p.end_process.iloc[0]
            mu_true_p = magdf(muon_row.pfp.trk.truth.p.genp).iloc[0]
            pi_true_p = magdf(pion_row.pfp.trk.truth.p.genp).iloc[0]          
            mu_true_costheta = muon_row.pfp.trk.truth.p.genp.z.iloc[0]/mu_true_p
            pi_true_costheta = pion_row.pfp.trk.truth.p.genp.z.iloc[0]/pi_true_p
        else:
            mu_true_pdg = -1   
            mu_true_p_type = "none"
            mu_true_end_process = -1
            pi_true_pdg = -1
            pi_true_p_type = "none"
            pi_true_end_process = -1
            mu_true_p = -1
            pi_true_p = -1         
            mu_true_costheta = -1
            pi_true_costheta = -1
        
        
    # Case 2: no exiting PFP
    else:
        group_sorted = group.sort_values(('pfp','trk','bdt_muon_pion_score','','',''))
        
        p_mu = group_sorted.pfp.trk.rangeP.p_muon.iloc[1]
        cos_theta_mu = group_sorted.pfp.trk.dir.z.iloc[1]
        mu_chi2_proton = group_sorted.pfp.trk.chi2pid.best.chi2_proton.iloc[1]
        mu_chi2_mu = group_sorted.pfp.trk.chi2pid.best.chi2_muon.iloc[1]
        mu_chi2_exp_pol = group_sorted.pfp.trk.chi2_exp_pol.iloc[1]
        mu_scatter_angle_ratio = group_sorted.pfp.scatter_angle_ratio.iloc[1]
        mu_max_daughter_hits = group_sorted.pfp.max_daughter_hits.iloc[1]
        mu_frac50 = group_sorted.pfp.trk.frac50.iloc[1]
        mu_bdt_score_proton = group_sorted.pfp.trk.bdt_proton_score.iloc[1]
        mu_bdt_score_muon_pion = group_sorted.pfp.trk.bdt_muon_pion_score.iloc[1]
        
        # split components
        p_mu_x = p_mu * group_sorted.pfp.trk.dir.x.iloc[1]
        p_mu_y = p_mu * group_sorted.pfp.trk.dir.y.iloc[1]
        p_mu_z = p_mu * group_sorted.pfp.trk.dir.z.iloc[1]

            
        p_pi_range = group_sorted.pfp.trk.rangeP.p_pion.iloc[0]
        cos_theta_pi = group_sorted.pfp.trk.dir.z.iloc[0]
        pion_id = group_sorted.index[0] # This is the (ntuple, entry, slc, pfp) tuple
        best_plane = group_sorted.pfp.trk.bestplane.iloc[0]
        #best_plane = 2
        p_pi_TLE = calculate_hypfit_p(track_id=pion_id, best_plane=best_plane, hit_df=best_hit_df, target_pdg=211, cleaning="all")
        #p_pi_TLE = 0.5
        pi_chi2_proton = group_sorted.pfp.trk.chi2pid.best.chi2_proton.iloc[0]
        pi_chi2_mu = group_sorted.pfp.trk.chi2pid.best.chi2_muon.iloc[0]
        pi_chi2_exp_pol = group_sorted.pfp.trk.chi2_exp_pol.iloc[0]
        pi_scatter_angle_ratio = group_sorted.pfp.scatter_angle_ratio.iloc[0]
        pi_max_daughter_hits = group_sorted.pfp.max_daughter_hits.iloc[0]
        pi_frac50 = group_sorted.pfp.trk.frac50.iloc[0]
        pi_bdt_score_proton = group_sorted.pfp.trk.bdt_proton_score.iloc[0]
        pi_bdt_score_muon_pion = group_sorted.pfp.trk.bdt_muon_pion_score.iloc[0]
       
        # split components
        p_pi_x = p_pi_TLE * group_sorted.pfp.trk.dir.x.iloc[0]
        p_pi_y = p_pi_TLE * group_sorted.pfp.trk.dir.y.iloc[0]
        p_pi_z = p_pi_TLE * group_sorted.pfp.trk.dir.z.iloc[0]
        
        muon_contained = True

        if(group_sorted.pfp.trk.truth.p.pdg.iloc[1] > -2147483648):
            mu_true_pdg = group_sorted.pfp.trk.truth.p.pdg.iloc[1]
            mu_true_p_type = group_sorted.pfp.trk.truth.p.p_type.iloc[1]
            mu_true_end_process = group_sorted.pfp.trk.truth.p.end_process.iloc[1]
            pi_true_pdg = group_sorted.pfp.trk.truth.p.pdg.iloc[0]
            pi_true_p_type = group_sorted.pfp.trk.truth.p.p_type.iloc[0]
            pi_true_end_process = group_sorted.pfp.trk.truth.p.end_process.iloc[0]
            mu_true_p = magdf(group_sorted.pfp.trk.truth.p.genp).iloc[1]
            pi_true_p = magdf(group_sorted.pfp.trk.truth.p.genp).iloc[0]
            mu_true_costheta = group_sorted.pfp.trk.truth.p.genp.z.iloc[1]/mu_true_p
            pi_true_costheta = group_sorted.pfp.trk.truth.p.genp.z.iloc[0]/pi_true_p

        else:
            mu_true_pdg = -1   
            mu_true_p_type = "none"
            mu_true_end_process = -1
            pi_true_pdg = -1
            pi_true_p_type = "none"
            pi_true_end_process = -1
            mu_true_p = -1
            pi_true_p = -1         
            mu_true_costheta = -1
            pi_true_costheta = -1
        
    return pd.Series({
        'p_mu_x': p_mu_x,
        'p_mu_y': p_mu_y,
        'p_mu_z': p_mu_z,
        'p_pi_x': p_pi_x,
        'p_pi_y': p_pi_y,
        'p_pi_z': p_pi_z,
        
        'reco_p_mu': p_mu,
        'cos_theta_mu': cos_theta_mu,
        'mu_chi2_proton': mu_chi2_proton,
        'mu_chi2_mu': mu_chi2_mu,
        'mu_chi2_exp_pol': mu_chi2_exp_pol,
        'mu_scatter_angle_ratio': mu_scatter_angle_ratio,
        'mu_max_daughter_hits': mu_max_daughter_hits,
        'mu_frac50': mu_frac50,
        'mu_bdt_score_proton': mu_bdt_score_proton,
        'mu_bdt_score_muon_pion': mu_bdt_score_muon_pion,

        'range_p_pi': p_pi_range,
        'TLE_p_pi': p_pi_TLE,
        'cos_theta_pi': cos_theta_pi,        
        'pi_chi2_proton': pi_chi2_proton,
        'pi_chi2_mu': pi_chi2_mu,
        'pi_chi2_exp_pol': pi_chi2_exp_pol,
        'pi_scatter_angle_ratio': pi_scatter_angle_ratio,
        'pi_max_daughter_hits': pi_max_daughter_hits,
        'pi_frac50': pi_frac50,
        'pi_bdt_score_proton': pi_bdt_score_proton,
        'pi_bdt_score_muon_pion': pi_bdt_score_muon_pion,

        'muon_contained': muon_contained,

        'mu_true_pdg': mu_true_pdg,
        'mu_true_p_type': mu_true_p_type,
        'mu_true_end_process': mu_true_end_process,
        'pi_true_pdg': pi_true_pdg,
        'pi_true_p_type': pi_true_p_type,
        'pi_true_end_process': pi_true_end_process,
        
        'pi_true_p':pi_true_p,
        'mu_true_p':mu_true_p,
        'pi_true_costheta':pi_true_costheta,
        'mu_true_costheta':mu_true_costheta             
    })
    
'''
def get_mu_pi_vars(group, best_hit_df):
    """
    Selects a Muon and a Pion from a group of PFPs based on:
    1. Muon: Exiting track preferred; otherwise highest BDT muon/pion score.
    2. Pion: Longest remaining track (by range) after muon is removed.
    """
    
    # --- 1. SAFETY CHECK ---
    # If less than 2 PFPs, we cannot have a Muon + Pion pair.
    if len(group) < 2:
        return pd.Series({  # <--- MUST BE pd.Series
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
    # --- 2. MUON SELECTION ---
    exiting_mask = CutMasks.exiting_pfp_mask(group)
    
    if exiting_mask.sum() >= 1:
        # Case A: At least one exiting track. Pick the first exiting one as the muon.
        muon_row = group.loc[exiting_mask].iloc[[0]]
        muon_contained = False
        # For exiting muons, we usually use MCS momentum
        p_mu = muon_row.pfp.trk.mcsP.fwdP_muon.iloc[0]
    else:
        # Case B: All tracks contained. Pick the one with the highest BDT muon/pion score.
        # We sort by score and take the last one (highest).
        group_sorted = group.sort_values(('pfp','trk','bdt_muon_pion_score','','',''))
        muon_row = group_sorted.iloc[[-1]]
        muon_contained = True
        # For contained muons, we use Range momentum
        p_mu = muon_row.pfp.trk.rangeP.p_muon.iloc[0]
    
    # --- 3. PION SELECTION ---
    # The pion is the longest remaining track. 
    # We drop the muon's index to ensure we don't pick the same track twice.
    remaining_pfps = group.drop(muon_row.index)
    
    # Sort remaining by rangeP.p_pion (proxy for length) and take the highest
    pion_row = remaining_pfps.sort_values(('pfp','trk','len','','','')).iloc[[-1]]
    
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
                                  hit_df=best_hit_df, target_pdg=211, cleaning="all")
    
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
    
    return nudf


def make_cc1pinudf_expanded_syst(f):
    return make_cc1pinudf(f, include_weights=True, multisim_nuniv=100, genie_multisim_nuniv=100, wgt_types=["bnb","genie","g4"], slim=False, genie_systematics=None)

def make_cc1pinudf_Ar23p_expanded_syst(f):
    return make_cc1pinudf(f, include_weights=True, multisim_nuniv=100, genie_multisim_nuniv=100, wgt_types=["bnb","genie","g4"], slim=False, genie_systematics=ar23p_genie_systematics + regen_systematics)

def make_cc1pinudf_Ar23p(f):
    return make_cc1pinudf(f, include_weights=True, multisim_nuniv=100, genie_multisim_nuniv=100, wgt_types=["bnb","genie","g4"], slim=True, genie_systematics=ar23p_genie_systematics + regen_systematics)

    
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
        ('slc', 'measure_var', 'min_angle_between_candidates', '', '', ''),
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

def make_cc1pi_finaldf(f, updatecalo = None):
    
    pandora_df  =  make_pandora_df(f, trkScoreCut = False, trkDistCut= -1, cutClearCosmic = True, requireFiducial=False, updatecalo=updatecalo)
    if pandora_df.empty:
       # Define the 3-level empty MultiIndex
        idx = pd.MultiIndex(
            levels=[[], [], []],  # Three levels
            codes=[[], [], []],   # Three levels
            names=['entry', 'rec.slc..index', 'rec.slc.reco.pfp..index']
        )
        
        # Create the DataFrame using the columns from your 'cols' list
        empty_df = pd.DataFrame(index=idx, columns=pd.MultiIndex.from_tuples(cols))
        
        return empty_df
        
    cc1pi_shwbranches = [
        shwbranch + 'bestplane_energy'
    ]
    shw_df = loadbranches(f["recTree"], cc1pi_shwbranches)
    shw_df = shw_df.rec.slc.reco

    #create mcs df
    mcs_df = loadbranches(f["recTree"], [trkmcsbranches[1]]).rec.slc.reco.pfp.trk.mcsP

    hit0_df = make_trkhitdf_plane0(f, updatecalo = updatecalo)
    hit1_df = make_trkhitdf_plane1(f, updatecalo = updatecalo)
    hit2_df = make_trkhitdf_plane2(f, updatecalo = updatecalo)
    
    hit0_df = hit0_df.sort_values('rr', ascending=True)
    hit1_df = hit1_df.sort_values('rr', ascending=True)
    hit2_df = hit2_df.sort_values('rr', ascending=True)
        
    hit_dfs = [hit0_df, hit1_df, hit2_df]
    hit_names = ['nhit0', 'nhit1', 'nhit2']  

    
    best_hit_df = get_best_hit_df(pandora_df, hit_dfs)
    pandora_df = add_nhit_column(pandora_df, hit_dfs, hit_names,1000)
    pandora_df = multicol_merge(pandora_df, shw_df, left_index=True, right_index=True, how="left", validate="one_to_one")  
    
    pandora_df = pandora_df[CutMasks.is_obvious_cosmic_cut_mask(pandora_df)]
    pandora_df[('slc', 'vertex_inside_fv', '', '', '', '')] = CutMasks.is_inside_FV_cut_mask(pandora_df)
    
    #pandora_df = pandora_df[pandora_df.slc.vertex_inside_fv == True]

    #take only primary particles for optimization of cosmic, track, shower,Michel, MIP and Angle cut
                
    pandora_df = add_best_ke_column(pandora_df)
    pandora_df = add_p_type_column(pandora_df)  
    pandora_df = add_best_chi2_columns(pandora_df, update_calo = updatecalo)
    pandora_df = add_max_daughter_hits_column(pandora_df)

    #Remove non primary particles (not usefull for analysis)
    pandora_df = pandora_df[pandora_df.pfp.parent_is_primary == True]
    pandora_df = add_n_protons_column(pandora_df)
    
    
    fixed_hit_df = (
        best_hit_df
        .groupby(level=group_levels, group_keys=False)
        .apply(dEdxCleaning.get_fix_hit_df, plot = False)
    )
    fixed_hit_df['dE'] = fixed_hit_df.dedx * fixed_hit_df.pitch
    
    pandora_df = add_mean_dedx_column(pandora_df,fixed_hit_df)  
    pandora_df = add_chi2_exp_pol_column(pandora_df, fixed_hit_df)
    #pandora_df = add_chi2_exp_pol_column_3var_fit(pandora_df, fixed_hit_df)
    pandora_df = add_frac50_column(pandora_df, fixed_hit_df)
    pandora_df = add_scatter_angle_ratio_column(pandora_df,mcs_df)
    pandora_df = add_max_angle_between_candidates_column(pandora_df)
    pandora_df = add_min_angle_between_candidates_column(pandora_df)
    pandora_df[('slc', 'measure_var', 'angle_between_candidates', '', '', '')] = pandora_df[('slc', 'measure_var', 'max_angle_between_candidates', '', '', '')]
    
    
    pandora_df[('pfp', 'is_exiting', '', '', '', '')] = CutMasks.exiting_pfp_mask(pandora_df)
    pandora_df[('pfp', 'is_exiting_z', '', '', '', '')] = CutMasks.exiting_z_pfp_mask(pandora_df)
    
    #add proton BDT column
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

    #add muon pion BDT column
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

    #Passing cuts or not boolean
    pandora_df[('slc', 'cut', 'obvious_cosmic', '', '', '')] = CutMasks.is_obvious_cosmic_cut_mask(pandora_df)
    pandora_df[('slc', 'cut', 't0', '', '', '')] = CutMasks.t0_cut_mask(pandora_df)
    pandora_df[('slc', 'cut', 'inside_FV', '', '', '')] = CutMasks.is_inside_FV_cut_mask(pandora_df)
    pandora_df[('slc', 'cut', 'nu_score', '', '', '')] = CutMasks.nu_score_cut_mask(pandora_df)
    pandora_df[('slc', 'cut', 'track', '', '', '')] = CutMasks.track_cut_mask(pandora_df, group_levels)
    pandora_df[('slc', 'cut', 'shower', '', '', '')] = CutMasks.shower_cut_mask(pandora_df, group_levels)
    pandora_df[('slc', 'cut', 'MIP_candidates', '', '', '')] = CutMasks.chi2_cut_mask(pandora_df, group_levels)
    pandora_df[('slc', 'cut', 'angle', '', '', '')] = CutMasks.angle_cut_mask(pandora_df)
    pandora_df[('slc', 'cut', 'containment', '', '', '')] = CutMasks.containment_z_cut_mask(pandora_df, group_levels)
    pandora_df[('slc', 'cut', 'michel', '', '', '')] = CutMasks.michel_cut_mask(pandora_df, group_levels)
    pandora_df[('slc', 'cut', 'extra_pion', '', '', '')] = CutMasks.extra_pion_cut_mask(pandora_df, group_levels)
    pandora_df[('slc', 'cut', 'proton_BDT', '', '', '')] = CutMasks.proton_BDT_cut_mask(pandora_df, group_levels)
    pandora_df[('slc', 'cut', 'proton_BDT_sideband', '', '', '')] = CutMasks.proton_BDT_sideband_mask(pandora_df, group_levels)

    
    #candidate_df = pandora_df[pandora_df.slc.cut.MIP_candidates & CutMasks.is_MIP_candidate_mask(pandora_df) & pandora_df.slc.cut.containment]
    candidate_df = pandora_df[pandora_df.slc.cut.inside_FV & pandora_df.slc.cut.t0 & pandora_df.slc.cut.track & CutMasks.is_MIP_candidate_mask(pandora_df) & pandora_df.slc.cut.containment]
 
    pandora_df = add_n_primary_tracks_column(pandora_df)
    pandora_df = add_n_primary_showers_column(pandora_df)
    pandora_df = add_n_primary_MIP_column(pandora_df)
    pandora_df = add_n_primary_MIP_proton_column(pandora_df)
    pandora_df = add_n_exiting_pfps_column(pandora_df)
    pandora_df = add_n_exiting_z_pfps_column(pandora_df)
    pandora_df = add_n_MIP_candidate_michel_column(pandora_df)
    pandora_df = add_n_extra_pions_column(pandora_df)

    EXPECTED_COLS = [
        'p_mu_x',
        'p_mu_y',
        'p_mu_z',
        'p_pi_x',
        'p_pi_y',
        'p_pi_z',
        
        'reco_p_mu',
        'cos_theta_mu',
        'mu_chi2_proton',
        'mu_chi2_mu',
        'mu_chi2_exp_pol',
        'mu_scatter_angle_ratio',
        'mu_max_daughter_hits',
        'mu_frac50',
        'mu_bdt_score_proton',
        'mu_bdt_score_muon_pion',

        'range_p_pi',
        'TLE_p_pi',
        'cos_theta_pi',        
        'pi_chi2_proton',
        'pi_chi2_mu',
        'pi_chi2_exp_pol',
        'pi_scatter_angle_ratio',
        'pi_max_daughter_hits',
        'pi_frac50',
        'pi_bdt_score_proton',
        'pi_bdt_score_muon_pion',

        'muon_contained',

        'mu_true_pdg',
        'mu_true_p_type',
        'mu_true_end_process',
        'pi_true_pdg',
        'pi_true_p_type',
        'pi_true_end_process',
        
        'pi_true_p',
        'mu_true_p',
        'pi_true_costheta',
        'mu_true_costheta'
    ]
    
    # 1. Define schema: default float32, override specific types
    schema = {col: 'float32' for col in EXPECTED_COLS}
    schema.update({
        'muon_contained': 'bool',      # nullable bool
        'mu_true_pdg': 'int32',           # nullable int
        'mu_true_end_process': 'int32',
        'pi_true_pdg': 'int32',
        'pi_true_end_process': 'int32',
        'mu_true_p_type': 'object',
        'pi_true_p_type': 'object'
    })

    
    #process only with plane 2
    # 2. Apply function → Series → DataFrame
    '''
    muon_pion_vars = (
        candidate_df
            .groupby(level=group_levels, group_keys=False)
            .apply(get_mu_pi_vars, best_hit_df=hit2_df)
    )
    '''
    muon_pion_vars = (
        candidate_df
            .groupby(level=group_levels, group_keys=False)
            .apply(get_mu_pi_vars, best_hit_df=best_hit_df)
    )
    
    # 3. Build slcdf aligned to all slices in pandora_df
    #    Ensures every slice has a row, even if missing in muon_pion_vars
    all_slices_idx = pandora_df.index.droplevel('rec.slc.reco.pfp..index').unique()
    
    if muon_pion_vars.empty:
        slcdf = pd.DataFrame(index=all_slices_idx, columns=EXPECTED_COLS)
    else:
        if isinstance(muon_pion_vars, pd.Series):
            # If it's a Series but has the Expected Columns as its index, 
            # we need to unstack it or convert it.
            # But usually, if apply returns a Series per group, 
            # the result is already a DF. 
            # If it's a single group, it might be a Series.
            if not isinstance(muon_pion_vars.index, pd.MultiIndex):
                 muon_pion_vars = muon_pion_vars.to_frame().T
            else:
                 # It's a Series with a MultiIndex (likely from a single group)
                 # We want the columns to be the inner-most labels
                 muon_pion_vars = muon_pion_vars.unstack()
                
        # Keep only expected columns
        slcdf = muon_pion_vars.reindex(columns=EXPECTED_COLS)
        # Reindex to include all slices
        slcdf = slcdf.reindex(all_slices_idx)
    
    # 4. Enforce schema and nullable dtypes
    '''
    for col, dtype in schema.items():
        # Replace missing values for non-nullable types
        if dtype == 'bool':
            slcdf[col] = slcdf[col].fillna(False).astype(bool)
        elif dtype == 'int32':
            slcdf[col] = slcdf[col].fillna(-9999).astype('int32')
        elif dtype == 'object':
            slcdf[col] = slcdf[col].fillna("none").astype('object')
        else:  # float32 or object
            slcdf[col] = slcdf[col].astype(dtype)
    '''
    for col, dtype in schema.items():
        if dtype == 'bool':
            slcdf[col] = slcdf[col].fillna(False).astype(bool)
        elif dtype == 'int32':
            # Use to_numeric first to safely handle any Series-like objects in the cells
            slcdf[col] = pd.to_numeric(slcdf[col], errors='coerce').fillna(-9999).astype('int32')
        elif dtype == 'object':
            slcdf[col] = slcdf[col].fillna("none").astype('object')
        else:  # float32 
            # Use to_numeric here as well to kill the FutureWarning
            slcdf[col] = pd.to_numeric(slcdf[col], errors='coerce').astype('float32')
            
    # 5. Optional: rename columns
    slcdf = slcdf.rename(columns={
        'cos_theta_mu': 'reco_cos_theta_mu',
        'cos_theta_pi': 'reco_cos_theta_pi'
    })
    
    # 6. Apply MultiIndex columns
    slcdf.columns = pd.MultiIndex.from_tuples(
        [('slc', 'measure_var', col, '', '', '') for col in slcdf.columns]
    )

    # 7. Join back safely
    pandora_df = pandora_df.join(
        slcdf,
        on=['entry', 'rec.slc..index']
    )

    pandora_df = add_transverse_vars_column(pandora_df)
    
    pandora_df[('slc', 'cut', 'energy', '', '', '')] = (pandora_df.slc.measure_var.reco_p_mu > 0.1) & (pandora_df.slc.measure_var.reco_p_mu < 3) & (pandora_df.slc.measure_var.TLE_p_pi > 0.13) & (pandora_df.slc.measure_var.TLE_p_pi < 2)   
    
    
    min_df = pandora_df[cols].copy()
    min_df = min_df[min_df.pfp.trk.len > 0]
    return min_df

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