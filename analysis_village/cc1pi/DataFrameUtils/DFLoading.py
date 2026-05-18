import pandas as pd
import pyanalib.pandas_helpers as ph
import pyanalib.split_df_helpers as splh
import pyanalib.stat_helpers as sh
from analysis_village.cc1pi.DataFrameUtils import DFCleaning
from analysis_village.cc1pi.CutMasks import CutMasks
from analysis_village.cc1pi.Constants import CTE

import numpy as np


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
 


def concat_shift_first_index(df1: pd.DataFrame, df2: pd.DataFrame) -> pd.DataFrame:
    # Safety checks
    if not isinstance(df1.index, pd.MultiIndex) or not isinstance(df2.index, pd.MultiIndex):
        raise ValueError("Both DataFrames must have MultiIndex")

    if df1.index.nlevels != df2.index.nlevels:
        raise ValueError("MultiIndex levels must match")

    # Get max of first level in df1
    max_first_level = df1.index.get_level_values(0).max()

    # Compute shift
    shift = max_first_level + 1

    # Extract df2 index as DataFrame
    index_df2 = df2.index.to_frame(index=False)

    # Shift first level
    index_df2.iloc[:, 0] = index_df2.iloc[:, 0] + shift

    # Rebuild MultiIndex
    new_index_df2 = pd.MultiIndex.from_frame(index_df2, names=df2.index.names)

    # Assign new index
    df2_shifted = df2.copy()
    df2_shifted.index = new_index_df2

    # Concatenate
    combined = pd.concat([df1, df2_shifted])

    return combined


def load_df(file, keys2load, n_max_concat = 100, filter_df = True, reprocess_df = True, reprocess_truth = True):
    
    print(f"keys in {file}")
    splh.print_keys(file)
    
    ## Check split multiplicity
    print("n_split: %d" %splh.get_n_split(file))
    
    ## Define keys to load
    print('dataframes')
     ## for big files, each key could have more than one split
    df = splh.load_dfs(file, keys2load, n_max_concat)
    print('loaded!')

    print(df['hdr']['pot'].sum())

    if reprocess_df:
        if "cc1pi" in keys2load:
            print("Changing CC1pi")
            df['cc1pi'][('slc', 'cut', 'proton_BDT_2pi', '', '', '')] = CutMasks.proton_BDT_cut_mask_2pi(df['cc1pi'], ['__ntuple', 'entry', 'rec.slc..index'])
            print("Adding proton BDT sidebans")
            df['cc1pi'][('slc', 'cut', 'proton_BDT_sideband', '', '', '')] = CutMasks.proton_BDT_sideband_mask(df['cc1pi'], ['__ntuple', 'entry', 'rec.slc..index'])
            print("Adding TPC containment")
            df['cc1pi'][('slc', 'cut', 'TPC_containment', '', '', '')] = CutMasks.TPC_containment_mask(df['cc1pi'], ['__ntuple', 'entry', 'rec.slc..index'])
            print("Adding FV")
            df['cc1pi'][('slc', 'cut', 'inside_FV', '', '', '')] = CutMasks.InFV_strict(df['cc1pi'])
            print("Adding high y z ")
            df['cc1pi'][('slc', 'cut', 'no_high_yz', '', '', '')] = CutMasks.not_in_high_y_high_z_containment_mask(df['cc1pi'], ['__ntuple', 'entry', 'rec.slc..index'])
            print("Finish")
    
    df['cc1pi'][('slc', 'cut', 'energy', '', '', '')] = (df['cc1pi'].slc.measure_var.reco_p_mu > 0.1) & (df['cc1pi'].slc.measure_var.reco_p_mu < 1) & (df['cc1pi'].slc.measure_var.TLE_p_pi > 0.13) & (df['cc1pi'].slc.measure_var.TLE_p_pi < 2)   
    print("Change energy")
    
    if reprocess_truth:   
        if "nudf" in keys2load:
            print("Reprocessing truth")
            df['nudf'] = df['nudf'].loc[~df['nudf'].index.duplicated(keep='first')]
            print("CHANGING nudf")
            df['nudf'] = add_nu_categ_column(df['nudf'], True)
            df['nudf'] = add_nu_categ_proton_reduced_column(df['nudf'], True)
            df['nudf'] = add_genie_categ_column(df['nudf'], True)
     
    
    if filter_df:
        #Perform duplication validation
        print(f"duplication for {file}")
        DFCleaning.find_duplicate_run_evt_combinations(df['hdr'])
        DFCleaning.plot_duplicate_run_subrun_evt_distribution(df["hdr"], file)
        
        ### Filter the hdr DataFrame first, then filter other DataFrames by matching with the hdr DataFrame
        df["hdr"] = DFCleaning.filter_unique_events(df["hdr"])
        DFCleaning.find_duplicate_run_evt_combinations(df["hdr"])
        #filter the rest of dataframe keys
        for key in keys2load:
            if key == "hdr":
                continue
            df[key] = DFCleaning.filter_using_hdr(df[key], df["hdr"])
    return df


def perform_truth_matching(mc_evt_df,mc_nu_df):
    new_columns = []
    for c in mc_nu_df.columns:
        new_columns.append(('truth',) + c + ('',) + ('',))  # prepend 'truth'
    mc_nu_df.columns = pd.MultiIndex.from_tuples(new_columns)

    matchdf = ph.multicol_merge(mc_evt_df.reset_index(), mc_nu_df.reset_index(),
                                left_on=[("__ntuple", "","","","",""),("entry", "","","","",""), ("slc", "tmatch","idx","","","")],
                                right_on=[("__ntuple", "","","","",""),("entry", "","","","",""), ("rec.mc.nu..index", "","","","","")], 
                                how="left") ## -- save all sllices
    #Reindex so it is again "__ntuple","entry", "slice_id"
    matchdf = matchdf.set_index(mc_evt_df.index.names, verify_integrity=True)
    #Remove "rec.mc.nu..index"
    matchdf = matchdf.drop(columns=[('rec.mc.nu..index','','','','','')])
    matchdf.loc[:, ('truth', 'nu_categ','','','','')] = (
        matchdf.loc[:, ('truth', 'nu_categ','','','','')].fillna('cosmic')
    )

    matchdf.loc[:, ('truth', 'genie_categ','','','','')] = (
        matchdf.loc[:, ('truth', 'genie_categ','','','','')].fillna('cosmic')
    )
    
    matchdf.loc[:, ('truth', 'nu_categ_proton_reduced','','','','')] = (
        matchdf.loc[:, ('truth', 'nu_categ_proton_reduced','','','','')].fillna('cosmic')
    )

    return matchdf

'''
import gc
import pandas as pd
from tqdm import tqdm  # Use standard tqdm instead of tqdm.auto

def perform_truth_matching_low_memmory(mc_evt_df, mc_nu_df, ntuple_chunk_size=10):
    # 1. Ensure both are sorted for fast slicing
    mc_evt_df = mc_evt_df.sort_index()
    mc_nu_df = mc_nu_df.sort_index()
    
    unique_ntuples = mc_evt_df.index.get_level_values("__ntuple").unique()
    num_ntuples = len(unique_ntuples)
    matched_chunks = []

    # 2. Setup the Progress Bar
    # Use tqdm(range(...)) to track the steps through the ntuple list
    pbar = tqdm(range(0, num_ntuples, ntuple_chunk_size), desc="Truth Matching")

    for i in pbar:
        # Select a range of ntuples for this batch
        batch_ntuples = unique_ntuples[i : i + ntuple_chunk_size]
        
        # Use isin() or .loc with the explicit list for safer slicing
        evt_chunk = mc_evt_df.loc[batch_ntuples]
        nu_chunk = mc_nu_df.loc[batch_ntuples]
        
        # Merge the batch
        matched_batch = perform_truth_matching(evt_chunk, nu_chunk)
        matched_chunks.append(matched_batch)

        # Update the progress bar suffix with current ntuple count
        current_count = min(i + ntuple_chunk_size, num_ntuples)
        pbar.set_postfix({"ntuples": f"{current_count}/{num_ntuples}"})
        
        # Memory Management
        del evt_chunk
        del nu_chunk
        gc.collect()

    print("\nFinalizing concatenation...")
    return pd.concat(matched_chunks)
'''

import gc
import pandas as pd
from tqdm import tqdm

def perform_truth_matching_low_memmory(mc_evt_df, mc_nu_df, ntuple_chunk_size=10):
    mc_evt_df = mc_evt_df.sort_index()
    mc_nu_df = mc_nu_df.sort_index()
    
    unique_ntuples = mc_evt_df.index.get_level_values("__ntuple").unique()
    num_ntuples = len(unique_ntuples)

    # 1. We use a generator function to yield batches
    # This prevents the creation of a massive intermediate list object
    def batch_generator():
        pbar = tqdm(range(0, num_ntuples, ntuple_chunk_size), desc="Matching Batches")
        for i in pbar:
            batch_ntuples = unique_ntuples[i : i + ntuple_chunk_size]
            
            # Use .loc with a slice if the index is sorted (much faster than .isin)
            start, end = batch_ntuples[0], batch_ntuples[-1]
            evt_chunk = mc_evt_df.loc[start:end].copy()
            nu_chunk = mc_nu_df.loc[start:end].copy()
            
            matched_batch = perform_truth_matching(evt_chunk, nu_chunk)
            
            if not matched_batch.empty:
                yield matched_batch

            # Explicitly clear chunk memory before the next yield
            del evt_chunk
            del nu_chunk
            gc.collect()

    print("Starting Truth Matching...")
    
    # 2. Use pd.concat directly on the generator
    # This is more memory efficient than building a list first
    final_df = pd.concat(batch_generator(), copy=False)
    
    # 3. Final cleanup of the generator and intermediate data
    gc.collect()
    
    print("Finalizing Result...")
    return final_df