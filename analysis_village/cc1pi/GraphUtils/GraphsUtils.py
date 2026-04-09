import pandas as pd
import numpy as np
from matplotlib.colors import to_rgba
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib.lines as mlines
from matplotlib.colors import LinearSegmentedColormap
from dataclasses import dataclass, field

# --------------------------------------------------
# Config class (columns, bins, colors, labels)
# --------------------------------------------------

# Apply ROOT-like global styles
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['DejaVu Sans', 'Arial', 'Helvetica', 'sans-serif'],
    'font.weight': 'normal',
    'axes.labelweight': 'normal',
    'axes.titleweight': 'normal',
    'axes.linewidth': 2,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'axes.labelsize': 16,
    'legend.frameon': True,
    'mathtext.default': 'it',
})


@dataclass
class HistogramConfig:
    data_column: tuple               # Column with the data to plot
    bins: np.ndarray = np.linspace(-0.5, 0.5, 51)
    color_map: dict = field(default_factory=dict)
    xlabel: str = 'Score'
    ylabel: str = 'Entries'
    title: str = None
    
COLORS = [
    "#0072B2",  # strong blue
    "#d62728",  # red
    "#E69F00",  # orange
    "#009E73",  # bluish green
    "#D55E00",  # vermillion
    "#CC79A7",  # reddish purple
    "#F0E442",  # yellow
    "#000000",  # black
    "#56B4E9",  # sky blue
    "#8C564B",  # brown (distinct from orange)
    "#7F7F7F",  # neutral gray

    # --- Replaced ones (more separated hues) ---
    "#A6761D",  # warm brown (distinct from orange)
    "#F4C2C2",  # vivid magenta (distinct from purple)
    "#98FB98",  # lime green (distinct from teal/green)
    "#7570B3",  # slate violet (cool toned, not deep purple)
    "#A6CEE3",  # pale blue (very different brightness)
]

from matplotlib.colors import LinearSegmentedColormap
import numpy as np

# 1. RGB colors in 0-255 format
sunset_rgb = [
    (247, 242, 229),     # dark purple/black
    (236, 226, 192), # light pink (this will be your second color)
    (227, 158, 62),  # orange
    (217, 84, 34),   # deep orange
    (157, 39, 46),   # red
    (30, 9, 31),     # dark purple/black
]

# Normalize to 0-1
sunset_colors = [(r/255, g/255, b/255) for r, g, b in sunset_rgb]

# 2. Add White at the very beginning
# We use 'nodes' to force the white to occupy only the very bottom
# 0.0 is White, 0.0001 starts the light pink
colors_with_white = [(1, 1, 1)] + sunset_colors
nodes = [0.0, 0.0001] + list(np.linspace(0.25, 1.0, len(sunset_colors)-1))

sunset_cmap = LinearSegmentedColormap.from_list(
    "kSunset_white_min", 
    list(zip(nodes, colors_with_white)), 
    N=256
)

# 3. Optional: Ensure empty/NaN bins also show as white
sunset_cmap.set_bad(color='white')


def plot_hist_with_outline(
    ax,
    scores,
    bins=50,
    range=None,
    weights=None,
    color="blue",
    alpha_fill=0.3,
    linewidth_fill=0,
    linewidth_outline=2.0,
    label=None,
    stacked=False
):

    # Filled histogram
    ax.hist(
        scores,
        bins=bins,
        range=range,
        weights=weights,
        histtype="stepfilled",
        color=color,
        alpha=alpha_fill,
        linewidth=linewidth_fill,
        label=label,
        stacked=stacked
    )

    # Outline histogram on top
    ax.hist(
        scores,
        bins=bins,
        range=range,
        weights=weights,
        histtype="step",
        color=color,
        alpha=1.0,
        linewidth=linewidth_outline,
        stacked=stacked
    )


def make_hist_legend_from_ax(ax, alpha_fill=0.3, linewidth=2.0):
    handles, labels = ax.get_legend_handles_labels()
    legend_handles = []

    # Loop over all handles and pick only the filled histograms (Patch)
    for h, l in zip(handles, labels):
        # Some histograms may be Line2D (outline), skip them
        if isinstance(h, Patch):
            # Extract color from the Patch
            facecolor = to_rgba(h.get_facecolor(), alpha_fill)
            edgecolor = to_rgba(h.get_facecolor(), 1.0)
            proxy = Patch(facecolor=facecolor, edgecolor=edgecolor, linewidth=linewidth, label=l)
            legend_handles.append(proxy)

    return legend_handles

# --------------------------------------------------
# Plotting function
# --------------------------------------------------
def plot_stacked_histogram(
    df,
    config: HistogramConfig,
    type_column: tuple,
    first_per_slice: bool = False,
):
  

    slice_levels = ['__ntuple', 'entry', 'rec.slc..index']

    # Optionally reduce to first PFP per slice
    if first_per_slice:
        df = (
            df
            .groupby(level=slice_levels, sort=False)
            .first()
        )

    # Extract series
    data = df[config.data_column]
    types = df[type_column]

    # Unique type values (preserve order of appearance)
    type_values = list(types.dropna().unique())

    # Build stack data
    stack_data = {
        t: data[types == t].dropna()
        for t in type_values
    }

    # Colors (list-based, stable)
    colors = [COLORS[i % len(COLORS)] for i in range(len(type_values))]

    # Create figure
    fig, ax = plt.subplots()

    # Filled stack
    ax.hist(
        [stack_data[t] for t in type_values],
        bins=config.bins,
        stacked=True,
        histtype='stepfilled',
        color=colors,
        alpha=0.1,
        linewidth=0,
    )

    # Outline stack
    ax.hist(
        [stack_data[t] for t in type_values],
        bins=config.bins,
        stacked=True,
        histtype='step',
        color=colors,
        linewidth=2.0,
    )

    # Legend
    legend_handles = [
        Patch(
            facecolor=to_rgba(c, 0.1),
            edgecolor=to_rgba(c, 1.0),
            linewidth=2.0,
            label=t
        )
        for c, t in zip(colors, type_values)
    ]

    ax.set_xlabel(config.xlabel)
    ax.set_ylabel(config.ylabel)
    if config.title:
        ax.set_title(config.title)

    ax.legend(handles=legend_handles, loc='upper right')

    plt.show()