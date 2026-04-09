import sys
import numpy as np
import math
import uproot as uproot
import pickle
import pandas as pd
import matplotlib.pyplot as plt

MIN_SIZE_FOR_FIX = 5

MAX_DEDX = 1000

DISJOINED_THRESHOLD = 10
MIN_SEGMENT_SIZE = 10

BRAGG_PEAK_AREA_RR = 3
MIN_MAX_FOR_BRAGG_CORRECTION = 4

MAX_RESIDUAL_RANGE_CUT = 30

VERTEX_AREA_RR = 5
TRUNCATED_MEAN_RR = 2

def plot_group(group, rr_col="rr", dedx_col="dedx", title=None):
    """
    Plot rr vs dEdx for a single group.

    Parameters
    ----------
    group : pd.DataFrame
        One calorimetry group
    rr_col : str
        Residual range column name
    dedx_col : str
        dE/dx column name
    title : str, optional
        Plot title
    """

    df = group.sort_values(rr_col)

    plt.figure(figsize=(7, 5))
    plt.scatter(df[rr_col], df[dedx_col], s=35)

    plt.xlabel("Residual range [cm]")
    plt.ylabel("dE/dx")

    if title:
        plt.title(title)

    plt.grid(True)
    plt.tight_layout()
    plt.show()

def plot_segments(group, segments, title=None):
    plt.figure(figsize=(7, 5))

    # Plot full group in light gray (context)
    plt.scatter(
        group["rr"],
        group["dedx"],
        alpha=0.3,
        label="All hits"
    )

    # Plot each segment with a different color
    for i, seg in enumerate(segments):
        plt.scatter(
            seg["rr"],
            seg["dedx"],
            label=f"Segment {i}",
            s=30
        )

    plt.xlabel("Residual range [cm]")
    plt.ylabel("dE/dx")
    if title:
        plt.title(title)

    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def get_disjoined_segments(df):
    segments = []
    current = []

    rr = df["rr"].to_numpy()
    new_segment = np.abs(np.diff(rr, prepend=rr[0])) > DISJOINED_THRESHOLD
    segment_id = np.cumsum(new_segment)

    return [seg for _, seg in df.groupby(segment_id)]

def get_best_segment(segments):
    candidates = []

    for i, seg in enumerate(segments):
        if len(seg) < MIN_SEGMENT_SIZE:
            continue
        min_rr = seg["rr"].min()
        candidates.append((min_rr, i))

    if not candidates:
        return segments[0]

    _, idx = min(candidates, key=lambda x: x[0])
    chosen = segments[idx].copy()

    chosen["rr"] -= chosen["rr"].min()
    return chosen

def correct_area_close_to_vertex(df):
    max_rr = df["rr"].max()
    if max_rr < VERTEX_AREA_RR:
        return df

    mask = (df["rr"] < max_rr - TRUNCATED_MEAN_RR) & (df["rr"] > TRUNCATED_MEAN_RR)
    ref = df.loc[mask, "dedx"]
    mean = ref.mean()
    std = ref.std(ddof=0)

    # Start with everything kept
    keep_mask = np.ones(len(df), dtype=bool)
    
    # Condition: Eliminate high dEdx points close to vertex
    keep_mask &= ~((df["rr"] > max_rr - VERTEX_AREA_RR) & (np.abs(df["dedx"] - mean) > 5 * std)
    )
    # Condition: Eliminate low dEdx points close to vertex (more restrictive)
    keep_mask &= ~(
        (df["rr"] > max_rr - VERTEX_AREA_RR) & (df["dedx"] < mean) & (np.abs(df["dedx"] - mean) > 2 * std)
    )

    return df.loc[keep_mask]

def correct_bragg_peak(df):
    max_rr = df["rr"].max()
    if max_rr < BRAGG_PEAK_AREA_RR:
        return df

    # Bragg peak region
    bragg_area_df = df[df["rr"] <= BRAGG_PEAK_AREA_RR]
    if bragg_area_df.empty:
        return df

        MIN_MAX_FOR_BRAGG_CORRECTION

    max_bragg = bragg_area_df["dedx"].max()
    if max_bragg < MIN_MAX_FOR_BRAGG_CORRECTION:
        return df
    
    rr_max = bragg_area_df.loc[bragg_area_df["dedx"].idxmax(), "rr"]
    drop_mask = (
        (df["rr"] < rr_max) & (df["dedx"] < max_bragg)
    )
    corrected = df.loc[~drop_mask].copy()
    
    # Re-zero residual range
    corrected["rr"] -= corrected["rr"].min()

    return corrected


def get_fix_hit_df(group, plot = False):
    
    if len(group) < MIN_SIZE_FOR_FIX:
        return group
        
    fixed = group[group["dedx"] <= MAX_DEDX]
    segments = get_disjoined_segments(fixed)
    if(plot):
        plot_segments(fixed, segments, title="Segments")
    
    fixed = get_best_segment(segments)
    if(plot):
        plot_group(fixed, title= "Best Segment")
    
    fixed = correct_bragg_peak(fixed)   
    if(plot):
        plot_group(fixed, title= "Bragg peak fix")
    
    if fixed["rr"].max() > MAX_RESIDUAL_RANGE_CUT:
        fixed = fixed[fixed["rr"] <= MAX_RESIDUAL_RANGE_CUT]
    if(plot):
        plot_group(fixed, title= "RR trimmed")
        
    fixed = correct_area_close_to_vertex(fixed)
    if(plot):
        plot_group(fixed, title= "Vertex Area Corrected")

    return fixed