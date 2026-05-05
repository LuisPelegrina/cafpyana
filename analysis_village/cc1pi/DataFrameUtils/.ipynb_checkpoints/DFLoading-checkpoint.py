import pandas as pd
import pyanalib.pandas_helpers as ph
import pyanalib.split_df_helpers as splh
import pyanalib.stat_helpers as sh
from analysis_village.cc1pi.DataFrameUtils import DFCleaning
from analysis_village.cc1pi.CutMasks import CutMasks
from analysis_village.cc1pi.Constants import CTE

import numpy as np




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


def load_df(file, keys2load, n_max_concat = 100, filter_df = True):
    
    print(f"keys in {file}")
    splh.print_keys(file)
    
    ## Check split multiplicity
    print("n_split: %d" %splh.get_n_split(file))
    
    ## Define keys to load
    print('dataframes')
     ## for big files, each key could have more than one split
    df = splh.load_dfs(file, keys2load, n_max_concat)
    print('loaded!')
    
    if "cc1pi" in keys2load:
        df['cc1pi'][('slc', 'cut', 'proton_BDT_2pi', '', '', '')] = CutMasks.proton_BDT_cut_mask_2pi(df['cc1pi'], ['__ntuple', 'entry', 'rec.slc..index'])
        #df['cc1pi'][('slc', 'cut', 'TPC_containment', '', '', '')] = CutMasks.TPC_containment_mask(df['cc1pi'], ['__ntuple', 'entry', 'rec.slc..index'])

    '''
    if "nudf" in keys2load:
        print("CHANGING nudf")
        df['nudf'] = add_nu_categ_column(df['nudf'], True)
        df['nudf'] = add_nu_categ_proton_reduced_column(df['nudf'], True)
        df['nudf'] = add_genie_categ_column(df['nudf'], True)
    '''
    
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