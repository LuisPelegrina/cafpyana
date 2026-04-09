import sys
import os
from typing import List, Dict
import pandas as pd

# Absolute path to cafpyana directory
cafpyana_root = "/home/lpelegri/cafpyana"
# Add CAFpyana to the Python search path -- this will allow you to import cafpyana modules
sys.path.insert(0, cafpyana_root)

# import this repo's classes
import pyanalib.pandas_helpers as ph
import pyanalib.split_df_helpers as splh
import pyanalib.stat_helpers as sh

def filename_to_dataframe(filenames: List[str], keys2load: List[str], n_max_concat: int = 2) -> Dict[str, pd.DataFrame]:
    """
    Given a list with the names of the files, returns the corresponding MultiIndex DataFrames.
    Parameters:
    -----------
        - filenames: List of filenames (strings)
        - n_max_concat: max number of concatenations per key
        - keys2load: List of keys to load from all files
    -----------
    Returns:
    -----------
        - Dict[str, pd.DataFrame]
        A dictionary where each key is a filename and each value is the
        corresponding MultiIndex DataFrame.
    -----------
    """

    dataframes_MultiIndex: Dict[str, pd.DataFrame] = {}

    for name in filenames:

        ## Check keys in each file
        print(f"Keys in file {name}:")
        splh.print_keys(name)

        ## Check split multiplicity
        print(f"File {name} n_split: {splh.get_n_split(name)}")

        _df = splh.load_dfs(name, keys2load, n_max_concat)

        dataframes_MultiIndex[name] = _df
        # dataframes_flat.append( _df.columns = ['_'.join([str(c) for c in col if c]) for col in _df.columns] )
        # NEED TO BE ANOTHER FUNCTION

    return dataframes_MultiIndex

####
# I WILL HAVE A FUNCTION TO FLATTEN THE DATAFRAMES HERE
####