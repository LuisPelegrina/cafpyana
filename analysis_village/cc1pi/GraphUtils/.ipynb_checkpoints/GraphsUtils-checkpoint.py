import pandas as pd
import numpy as np
from matplotlib.colors import to_rgba
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib.lines as mlines
from matplotlib.colors import LinearSegmentedColormap
from dataclasses import dataclass, field

# ----------------import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D 
from matplotlib.patches import Patch
from matplotlib.colors import to_rgba
import matplotlib.gridspec as gridspec
from matplotlib.legend_handler import HandlerTuple, HandlerBase

import numpy as np
from scipy.stats import chi2 as chi2_dist

from analysis_village.cc1pi.GraphUtils.Utils import *

def get_stat_covariance_matrix(cv_contents, sum_w2):
    """
    cv_contents: array of bin contents (sum of weights)
    sum_w2: array of the sum of the squares of the weights per bin
    """
    cv_contents = np.asarray(cv_contents)
    sum_w2 = np.asarray(sum_w2)
    n_bins = len(cv_contents)

    # 1. Variance for weighted Poisson is Sum(W^2)
    cov = np.diag(sum_w2)

    # 2. Fractional covariance: Var / (Content^2) = Sum(W^2) / (Sum W)^2
    with np.errstate(divide='ignore', invalid='ignore'):
        # This is the squared fractional error
        frac_variance = np.where(cv_contents > 0, sum_w2 / (cv_contents**2), 0.0)
        cov_frac = np.diag(frac_variance)

    # 3. Correlation matrix
    corr = np.eye(n_bins)

    return {
        "cov": cov,
        "cov_frac": cov_frac,
        "corr": corr,
    }
    
from scipy.stats import chi2 as chi2_dist

def get_chi2(data, mc, cov, n_params=0):
    data = np.atleast_1d(data)
    mc   = np.atleast_1d(mc)

    if cov.shape != (len(data), len(data)):
        raise ValueError(f"Covariance shape {cov.shape} doesn't match data length {len(data)}")

    mask = (mc > 0) & np.isfinite(data) & np.isfinite(mc)

    data_f = data[mask]
    mc_f   = mc[mask]
    cov_f  = cov[np.ix_(mask, mask)]

    n_bins = mask.sum()
    ndof   = n_bins - n_params

    if ndof <= 0:
        return np.nan, ndof, np.nan

    if np.linalg.cond(cov_f) > 1e12:
        return np.nan, ndof, np.nan

    delta = mc_f - data_f
    try:
        chi2_val = delta @ np.linalg.solve(cov_f, delta)
    except np.linalg.LinAlgError:
        return np.nan, ndof, np.nan

    reduced_chi2 = chi2_val / ndof
    pval         = chi2_dist.sf(chi2_val, ndof)

    return chi2_val, ndof, pval
    
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
def plot_stacked_histogram_with_ratio(
    mc_df,
    data_df,
    config,
    cov_frac_matrix=None,
    cov_matrix=None,
    title: str = None,
    weight_column: tuple = None,
    data_pot: float = None,
    normalize: bool = False,
    show_stats: bool = True,
    symmetric_ratio: bool = False,
    divide_by_bin_width: bool = False
):
    # 1. Pre-processing and Scaling
    slice_levels = ['__ntuple', 'entry', 'rec.slc..index']
    if config.first_per_slice:
        mc_df   = mc_df.groupby(level=slice_levels, sort=False).first()
        data_df = data_df.groupby(level=slice_levels, sort=False).first()

    mc_data    = mc_df[config.var_evt_reco_col]
    mc_types   = mc_df[config.truth_column]
    mc_weights = mc_df[weight_column] if weight_column in mc_df.columns else pd.Series(1.0, index=mc_df.index)

    # --- Palette Selection ---
    present_categories = set(mc_types.dropna().unique())
    palettes = [category_colors, category_colors_pfp, proton_distinction_category_colors, genie_category_colors]
    chosen_map, max_overlap = category_colors, -1
    for p in palettes:
        overlap = len(present_categories.intersection(p.keys()))
        if overlap > max_overlap:
            max_overlap = overlap
            chosen_map  = p

    # 2. Sorting Categories
    category_totals = [(t, mc_weights[mc_types == t].sum()) for t in mc_types.dropna().unique()]
    signal_keys     = ["CC1pi"]
    signals     = sorted([x for x in category_totals if     x[0] in signal_keys], key=lambda x: x[1], reverse=True)
    backgrounds = sorted([x for x in category_totals if not x[0] in signal_keys], key=lambda x: x[1], reverse=True)

    sorted_types   = [x[0] for x in backgrounds + signals]
    stack_data_mc  = [mc_data[mc_types == t].dropna() for t in sorted_types]
    stack_weights  = [mc_weights[mc_types == t].loc[mc_data[mc_types == t].dropna().index] for t in sorted_types]

    # 3. Setup Figure
    fig      = plt.figure(figsize=(10, 8))
    gs       = gridspec.GridSpec(2, 1, height_ratios=[4, 1], hspace=0.07)
    ax_top   = fig.add_subplot(gs[0])
    ax_ratio = fig.add_subplot(gs[1], sharex=ax_top)

    # 4. TOP PLOT — clip if needed
    max_bin_edge      = config.bins[-1]
    stack_mc_clipped  = (
        [np.clip(d, config.bins[0], max_bin_edge) for d in stack_data_mc]
        if config.clip else stack_data_mc
    )

    colors        = [chosen_map.get(t, "#7f7f7f") for t in sorted_types]
    all_mc_weights = pd.concat(stack_weights)

    mc_sum, bins = np.histogram(
        pd.concat(stack_mc_clipped),
        bins=config.bins,
        weights=all_mc_weights
    )
    bin_centers = (bins[:-1] + bins[1:]) / 2
    bin_widths  = np.diff(bins)

    # --- MC UNCERTAINTY ---
    if cov_frac_matrix is not None:
        frac_err = np.sqrt(np.diag(cov_frac_matrix))
        mc_error = frac_err * mc_sum
    else:
        mc_sum_w2, _ = np.histogram(
            pd.concat(stack_mc_clipped),
            bins=config.bins,
            weights=all_mc_weights ** 2
        )
        mc_error         = np.sqrt(mc_sum_w2)
        stat_cov         = get_stat_covariance_matrix(mc_sum, mc_sum_w2)
        cov_matrix       = stat_cov["cov"]
        cov_frac_matrix  = stat_cov["cov_frac"]

    if normalize:
        norm_fact = mc_sum.sum()
        if norm_fact > 0:
            mc_sum    /= norm_fact
            mc_error  /= norm_fact
            stack_weights = [w / norm_fact for w in stack_weights]

    # --- DATA ---
    cols_to_check = [config.var_evt_reco_col]
    if weight_column in data_df.columns:
        cols_to_check.append(weight_column)
    valid_df     = data_df.dropna(subset=cols_to_check)
    data         = valid_df[config.var_evt_reco_col]
    data_weights = valid_df[weight_column] if weight_column in valid_df.columns else np.ones(len(data))

    data_clipped = np.clip(data, config.bins[0], max_bin_edge) if config.clip else data

    data_counts, _ = np.histogram(data_clipped, bins=bins, weights=data_weights)
    data_sum_w2, _ = np.histogram(data_clipped, bins=bins, weights=data_weights ** 2)

    # FIX 1: use data_sum_w2 consistently; never overwrite with sqrt(counts)
    if normalize:
        denom            = len(data) if len(data) > 0 else 1
        data_plot_counts = data_counts / denom
        data_errors      = np.sqrt(data_sum_w2) / denom
    else:
        data_plot_counts = data_counts.astype(float)
        data_errors      = np.sqrt(data_sum_w2)   # FIX: was sqrt(data_counts)

    # 6. STATISTICS
    
    # FIX 2: pass data_sum_w2 as second argument, not data_plot_counts
    ret_stats_data_rate = get_stat_covariance_matrix(data_plot_counts, data_sum_w2)
    #print(cov_matrix)
    #print("_----")
    #print(ret_stats_data_rate["cov"])
    chi2_val, ndof, p_val = get_chi2(data_counts, mc_sum, cov_matrix + ret_stats_data_rate["cov"])

    if divide_by_bin_width:
        mc_sum           /= bin_widths
        mc_error         /= bin_widths
        data_plot_counts /= bin_widths
        data_errors      /= bin_widths
        new_stack_weights = []
        for i, d in enumerate(stack_mc_clipped):
            bin_indices = np.clip(np.digitize(d, bins) - 1, 0, len(bin_widths) - 1)
            new_stack_weights.append(stack_weights[i] / bin_widths[bin_indices])
        stack_weights = new_stack_weights

    ax_top.hist(stack_mc_clipped, bins=bins, stacked=True, weights=stack_weights,
                histtype='stepfilled', color=colors, alpha=0.3)
    ax_top.hist(stack_mc_clipped, bins=bins, stacked=True, weights=stack_weights,
                histtype='step', color=colors, linewidth=2)
    ax_top.bar(bin_centers, 2 * mc_error, bottom=mc_sum - mc_error, width=bin_widths,
               edgecolor='grey', facecolor='grey', alpha=0.2, linewidth=0)
    ax_top.bar(bin_centers, 2 * mc_error, bottom=mc_sum - mc_error, width=bin_widths,
               edgecolor='grey', facecolor='none', hatch='////', alpha=0.5, linewidth=0)
    ax_top.errorbar(bin_centers, data_plot_counts, yerr=data_errors, xerr=bin_widths / 2,
                    fmt='ko', markersize=6, zorder=10, capsize=0)

    # 5. RATIO PLOT
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio       = np.divide(data_plot_counts, mc_sum, out=np.zeros_like(data_plot_counts), where=mc_sum != 0)
        ratio_error = np.divide(data_errors,       mc_sum, out=np.zeros_like(data_errors),       where=mc_sum != 0)
        mc_rel_error = np.divide(mc_error,          mc_sum, out=np.zeros_like(mc_error),           where=mc_sum != 0)

    valid_ratio = mc_sum > 0
    if np.any(valid_ratio) and symmetric_ratio:
        data_extrema  = np.maximum(
            np.abs((ratio[valid_ratio] + ratio_error[valid_ratio]) - 1),
            np.abs((ratio[valid_ratio] - ratio_error[valid_ratio]) - 1)
        )
        mc_extrema    = mc_rel_error[valid_ratio]
        max_deviation = np.max(np.maximum(data_extrema, mc_extrema))
        y_padding     = max(max_deviation * 1.4, 0.1)
        ax_ratio.set_ylim(1 - y_padding, 1 + y_padding)
    elif np.any(valid_ratio):
        max_deviation = np.max(np.abs(ratio[valid_ratio] - 1) + ratio_error[valid_ratio])
        y_padding     = max(min(max_deviation * 1.4, 0.75), 0.1)
        ax_ratio.set_ylim(1 - y_padding, 1 + y_padding)
    else:
        ax_ratio.set_ylim(0.5, 1.5)

    ax_ratio.bar(bin_centers, 2 * mc_rel_error, bottom=1 - mc_rel_error, width=bin_widths,
                 edgecolor='grey', facecolor='grey', alpha=0.2, linewidth=0)
    ax_ratio.bar(bin_centers, 2 * mc_rel_error, bottom=1 - mc_rel_error, width=bin_widths,
                 edgecolor='grey', facecolor='none', hatch='////', alpha=0.4, linewidth=0)
    ax_ratio.errorbar(bin_centers, ratio, yerr=ratio_error, xerr=bin_widths / 2,
                      fmt='ko', markersize=6, capsize=0)
    ax_ratio.axhline(1.0, color='#d62728', linestyle='--', linewidth=2)

    # 5.5 VERTICAL CUT LINE
    # FIX 3: guard against None cut_value
    cut_val  = getattr(config, 'cut_value', None)
    draw_cut = False
    if cut_val is not None:
        draw_cut = (
            any(v != -999 for v in cut_val)
            if isinstance(cut_val, (list, np.ndarray))
            else (cut_val != -999)
        )

    cut_hand = Line2D([0], [0], color='black', linestyle='--', linewidth=2, label='Selection Cut')
    if draw_cut:
        cuts_to_draw = cut_val if isinstance(cut_val, (list, np.ndarray)) else [cut_val]
        for ax in [ax_top, ax_ratio]:
            for val in cuts_to_draw:
                if val != -999:
                    ax.axvline(val, color='black', linestyle='--', linewidth=2, zorder=2)

    # 7. LEGEND
    total_data_counts = len(data_clipped)
    mc_hand  = [Patch(facecolor=to_rgba(chosen_map.get(t, "#7f7f7f"), 0.3),
                      edgecolor=chosen_map.get(t, "#7f7f7f"),
                      label=bkg_name_nice_map.get(t, t)) for t in sorted_types]
    err_label = 'MC Total Error' if show_stats else 'MC Stat. Error'
    err_hand  = Patch(edgecolor='grey', facecolor='none', hatch='////', alpha=0.5, label=err_label)
    dat_hand  = Line2D([0], [0], color='black', marker='o', linestyle='', label='Data', markersize=8)

    all_handles = mc_hand + [err_hand, dat_hand]
    all_labels  = [h.get_label() for h in mc_hand] + [err_label, 'Data']

    if draw_cut:
        all_handles.append(cut_hand)
        all_labels.append(cut_hand.get_label())

    if show_stats:
        # FIX 4: guard against ndof == 0 in the chi2 string
        if ndof > 0 and chi2_val is not None:
            chi2_str = f"$\\chi^{{2}}$ / ndf: {chi2_val:.2f} / {ndof} = {chi2_val/ndof:.3f}"
        else:
            chi2_str = f"$\\chi^{{2}}$: N/A"
        p_value_str      = f"$p_{{\\mathrm{{value}}}}$ = {p_val:.3f}" if p_val is not None else "$p$: N/A"
        data_str         = f'$N_{{\\mathrm{{Data}}}}$ = {total_data_counts}'
        chi2_handle      = Patch(color='none', label=chi2_str)
        p_value_handle   = Patch(color='none', label=p_value_str)
        N_data_evts_handle = Patch(color='none', label=data_str)
        all_handles += [chi2_handle, p_value_handle, N_data_evts_handle]
        all_labels  += [chi2_str, p_value_str, data_str]

    n_cols = (len(all_handles) + 3) // 4
    leg = ax_top.legend(
        handles=all_handles, labels=all_labels,
        loc='upper ' + config.stats_horizontal_alignment, ncol=n_cols,
        fontsize=12, framealpha=1.0, edgecolor='black', fancybox=False,
        borderaxespad=1, columnspacing=1.5, handlelength=1.5, handletextpad=0.5,
    )
    leg.get_frame().set_linewidth(1.5)

    plt.draw()
    if show_stats:
        texts = leg.get_texts()
        for t in texts[-3:]:
            t.set_position((-28, 0))

    # 8. FINAL STYLING
    ax_top.set_xlim(bins[0], bins[-1])
    ax_top.set_ylim(0, ax_top.get_ylim()[1] * 1.4)

    ylabel = config.ylabel
    if divide_by_bin_width:
        ylabel += " / bin width"

    ax_top.set_ylabel(f"{ylabel} (POT = {data_pot:.2e})")
    ax_ratio.set_ylabel("Data/MC")
    ax_ratio.set_xlabel(config.xlabel, fontsize=20)
    ax_ratio.tick_params(axis='x', which='both', direction='inout', length=6)
    ax_top.set_title("")

    plt.setp(ax_top.get_xticklabels(), visible=False)
    fig.align_ylabels([ax_top, ax_ratio])
    plt.subplots_adjust(top=0.92, bottom=0.12, left=0.12, right=0.95, hspace=0.07)

    plt.show()
    return fig, p_val