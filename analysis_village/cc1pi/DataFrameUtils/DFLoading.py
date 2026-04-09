import pandas as pd
import pyanalib.pandas_helpers as ph
import pyanalib.split_df_helpers as splh
import pyanalib.stat_helpers as sh
from analysis_village.cc1pi.DataFrameUtils import DFCleaning


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

    if filter_df:
        #Perform duplication validation
        print(f"duplication for {file}")
        #DFCleaning.find_duplicate_run_evt_combinations(df['hdr'])
        #DFCleaning.plot_duplicate_run_subrun_evt_distribution(df["hdr"], file)
        
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