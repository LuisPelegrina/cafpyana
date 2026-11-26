import sys
import numpy as np
import math
import uproot as uproot
import pickle
import pandas as pd

from analysis_village.cc1pi.GraphUtils.Utils import *
from analysis_village.cc1pi.Constants import CTE as CTE

def print_purity(df, nu_column=('truth', 'nu_categ', '', '', '', '')):
    # Extract columns and handle weights
    nu = df[nu_column]
    weights = df[('slc', 'wgt', '', '', '', '')]

    # Build per-slice dataframe
    # Using groupby and first() to ensure we count unique slices correctly
    tmp = pd.DataFrame({
        "nu": nu,
        "wgt": weights
    }).groupby(['__ntuple', 'entry', 'rec.slc..index']).first()

    # Calculate weighted counts for every category present in the data
    weighted_counts = tmp.groupby("nu")["wgt"].sum()
    total_weighted = weighted_counts.sum()

    print(f"Total weighted slices: {total_weighted:.3f}")
    print("Purity (weighted / total weighted)\n")

    # Iterate through all categories actually present in the weighted_counts index
    # This replaces the hardcoded 'final_states' list
    for category in weighted_counts.index:
        wcount = weighted_counts[category]
        purity = wcount / total_weighted if total_weighted > 0 else 0.0

        # Optional: Use your bkg_name_nice_map for cleaner printing
        # If the category isn't in your map, it defaults to the raw name
        pretty_name = bkg_name_nice_map.get(category, category)

        print(
            f"{pretty_name}: "
            f"purity: {purity*100:.4f}%  "
            f"counts: {wcount:.3f}"
        )

def print_category_metrics(initial_df, selected_df, target_categ="CC1pi"):
    """
    Calculates Efficiency, Purity, and Eff*Pur for a specific neutrino category.
    """
    
    def get_nu_counts(df):
        if df.empty:
            return {}
        # Identify unique slices and count occurrences of each neutrino category
        return (df.truth.nu_categ
                .groupby(['__ntuple', 'entry', 'rec.slc..index'])
                .first()
                .value_counts()
                .to_dict())

    # 1. Get counts for the denominator (initial) and numerator (selected)
    initial_counts = get_nu_counts(initial_df)
    selected_counts = get_nu_counts(selected_df)
    
    # 2. Extract specific counts for the signal category
    total_target_initial = initial_counts.get(target_categ, 0)
    total_target_selected = selected_counts.get(target_categ, 0)
    
    # 3. Total number of all slices passing the cut
    total_slices_selected = sum(selected_counts.values())
    
    # 4. Calculate Metrics (using decimals for the product calculation)
    eff = (total_target_selected / total_target_initial) if total_target_initial > 0 else 0
    pur = (total_target_selected / total_slices_selected) if total_slices_selected > 0 else 0
    eff_times_pur = eff * pur
    
    # 5. Output
    print(f"--- Metrics for Category: {target_categ} ---")
    print(f"Initial Signal Events:    {total_target_initial}")
    print(f"Selected Signal Events:   {total_target_selected}")
    print(f"Total Selected (S+B):     {total_slices_selected}")
    print(f"------------------------------------------")
    print(f"Efficiency:     {eff*100:.2f}%")
    print(f"Purity:         {pur*100:.2f}%")
    print(f"Eff * Pur:      {eff_times_pur:.4f}") 
    print(f"------------------------------------------\n")
    
    return {"eff": eff, "pur": pur, "eff_pur": eff_times_pur}
