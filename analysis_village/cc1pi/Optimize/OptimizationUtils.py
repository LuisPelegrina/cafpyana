
from analysis_village.cc1pi.GraphUtils.GraphsUtils import *
from sklearn.metrics import roc_curve
import numpy as np
import pandas as pd
from itertools import product
from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D # This is the fix

def optimize_cut_eff_pur(
    signal_df,
    bkg_df,
    column,
    cut_type=">",
    xlabel="Score",
    title="Cut optimization",
    signal_name=r"$\nu$",
    bkg_name="Background",
    xlim=None,
    nbins=50,
    legend_loc="upper right",
    cut_unit="",
    normalize_hist=False,
    min_pur=0.0  # <-- new parameter
):
    """
    Optimize a 1D cut by maximizing efficiency × purity,
    optionally requiring a minimum purity threshold.

    Parameters
    ----------
    min_pur : float
        Minimum purity threshold to consider a cut.
    """

    # -------------------------
    # Extract scores
    # -------------------------
    sig_scores = signal_df[column].dropna()
    bkg_scores = bkg_df[column].dropna()

    y_true = np.concatenate([
        np.ones(len(sig_scores)),
        np.zeros(len(bkg_scores)),
    ])

    scores = np.concatenate([
        sig_scores.values,
        bkg_scores.values,
    ])

    # -------------------------
    # Handle cut direction
    # -------------------------
    if cut_type == "<":
        scores_for_roc = -scores
        cut_label = f"{xlabel} < cut"
    elif cut_type == ">":
        scores_for_roc = scores
        cut_label = f"{xlabel} > cut"
    else:
        raise ValueError("cut_type must be '>' or '<'")

    # -------------------------
    # ROC
    # -------------------------
    fpr, tpr, thresholds = roc_curve(y_true, scores_for_roc)
    cuts = -thresholds if cut_type == "<" else thresholds

    valid = np.isfinite(cuts)
    cuts = cuts[valid]
    tpr = tpr[valid]
    fpr = fpr[valid]

    # -------------------------
    # Efficiency & purity
    # -------------------------
    n_signal = len(sig_scores)
    n_bkg = len(bkg_scores)

    n_signal_pass = tpr * n_signal
    n_bkg_pass = fpr * n_bkg

    eff = tpr
    pur = n_signal_pass / (n_signal_pass + n_bkg_pass)

    df = pd.DataFrame({
        "cut": cuts,
        "eff": eff,
        "pur": pur,
        "eff_x_pur": eff * pur,
    })

    # -------------------------
    # Apply min_pur threshold
    # -------------------------
    df = df[df["pur"] >= min_pur]

    if df.empty:
        print(f"No cuts satisfy min_pur >= {min_pur}")
        return df, None

    best = df.loc[
    df[
        (df["pur"] >= min_pur)
    ]["eff_x_pur"].idxmax()
    ]



    print("\nBest cut:")
    print(f"  Cut logic   : {cut_label}")
    print(f"  Cut value   : {best['cut']:.4f} {cut_unit}")
    print(f"  Efficiency  : {best['eff']:.4f}")
    print(f"  Purity      : {best['pur']:.4f}")
    print(f"  Eff × Pur   : {best['eff_x_pur']:.4f}")

    # -------------------------
    # Plot
    # -------------------------
    fig, ax1 = plt.subplots(figsize=(9, 6))

    # Histogram range
    hist_range = xlim if xlim is not None else (min(scores), max(scores))

    # Determine weights
    if normalize_hist:
        sig_weights = np.ones_like(sig_scores) / len(sig_scores)
        bkg_weights = np.ones_like(bkg_scores) / len(bkg_scores)
        ylabel = "Normalized entries"
    else:
        sig_weights = np.ones_like(sig_scores)  # adapt as needed
        bkg_weights = np.ones_like(bkg_scores)  # adapt as needed
        ylabel = "Events"

    # Histograms
    plot_hist_with_outline(
        ax1,
        scores=sig_scores,
        bins=nbins,
        range=hist_range,
        weights=sig_weights,
        color=COLORS[0],
        alpha_fill=0.3,
        linewidth_outline=2.0,
        label=signal_name
    )
    
    plot_hist_with_outline(
        ax1,
        scores=bkg_scores,
        bins=nbins,
        range=hist_range,
        weights=bkg_weights,
        color=COLORS[1],
        alpha_fill=0.3,
        linewidth_outline=2.0,
        label=bkg_name
    )

    # Best cut line
    ax1.axvline(
        best["cut"],
        color="black",
        linestyle="--",
        linewidth=2,
        label=f"Best cut = {best['cut']:.3f} {cut_unit}" if cut_unit else f"Best cut = {best['cut']:.3f}"
    )

    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)

    # Metrics on secondary axis
    ax2 = ax1.twinx()
    ax2.plot(df["cut"], df["eff"], label="Efficiency", color=COLORS[0], linestyle="--", linewidth=2)
    ax2.plot(df["cut"], df["pur"], label="Purity", color=COLORS[1], linestyle=":", linewidth=3)
    ax2.plot(df["cut"], df["eff_x_pur"], label="Eff × Pur", color="black", linewidth=2)
    ax2.set_ylabel("Metric value")
    ax2.set_ylim(0, 1.05)

    # Build combined legend automatically
    hist_legend_handles = make_hist_legend_from_ax(ax1, alpha_fill=0.3, linewidth=2.0)
    cut_line_handle = Line2D([0], [0], color="black", linestyle="--", linewidth=2,
                             label=f"Best cut = {best['cut']:.3f} {cut_unit}" if cut_unit else f"Best cut = {best['cut']:.3f}")
    lines2, labels2 = ax2.get_legend_handles_labels()
    all_handles = hist_legend_handles + [cut_line_handle] + lines2
    all_labels = [h.get_label() for h in hist_legend_handles] + [cut_line_handle.get_label()] + labels2
    ax1.legend(handles=all_handles, labels=all_labels, loc=legend_loc)

    if xlim is not None:
        ax1.set_xlim(xlim)

    plt.title(title)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    
    return df, best, fig

def manual_cut_optimizer(
    signal_df,
    bkg_df,
    cut_columns,
    cut_sign,
    cut_min,
    cut_max,
    n_steps=50,
    min_eff=0.0,
    min_pur=0.0,  # <-- new parameter
):
    """
    Manual grid-scan cut optimization with per-variable binning.

    Parameters
    ----------
    n_steps : int or list[int]
        Number of scan points per variable.
        If int, the same value is used for all variables.
    min_pur : float
        Minimum purity threshold to consider a cut combination.
    """

    n_vars = len(cut_columns)

    # -------------------------
    # Normalize n_steps
    # -------------------------
    if isinstance(n_steps, int):
        n_steps = [n_steps] * n_vars
    elif len(n_steps) != n_vars:
        raise ValueError("n_steps must be int or have same length as cut_columns")

    # -------------------------
    # Drop NaNs
    # -------------------------
    sig = signal_df[cut_columns].dropna()
    bkg = bkg_df[cut_columns].dropna()

    n_sig = len(sig)
    n_bkg = len(bkg)

    # -------------------------
    # Build cut grids
    # -------------------------
    grids = [
        np.linspace(cmin, cmax, steps)
        for cmin, cmax, steps in zip(cut_min, cut_max, n_steps)
    ]

    cut_labels = [f"cut_{i}" for i in range(n_vars)]
    total_combinations = np.prod([len(g) for g in grids])

    results = []

    # -------------------------
    # Grid scan with progress bar
    # -------------------------
    for cuts in tqdm(
        product(*grids),
        total=total_combinations,
        desc="Scanning cut space",
        unit="cuts"
    ):
        sig_mask = np.ones(n_sig, dtype=bool)
        bkg_mask = np.ones(n_bkg, dtype=bool)

        for col, sign, cut in zip(cut_columns, cut_sign, cuts):
            if sign == ">":
                sig_mask &= sig[col].values > cut
                bkg_mask &= bkg[col].values > cut
            elif sign == "<":
                sig_mask &= sig[col].values < cut
                bkg_mask &= bkg[col].values < cut
            else:
                raise ValueError("cut_sign must be '<' or '>'")

        sig_pass = sig_mask.sum()
        bkg_pass = bkg_mask.sum()

        eff = sig_pass / n_sig if n_sig else 0
        pur = sig_pass / (sig_pass + bkg_pass) if (sig_pass + bkg_pass) else 0
     
        results.append(
            dict(
                **{label: val for label, val in zip(cut_labels, cuts)},
                eff=eff,
                pur=pur,
                eff_x_pur=eff * pur
            )
        )

    results_df = pd.DataFrame(results)

    if results_df.empty:
        print("No cut combinations passed the min_eff and min_pur thresholds.")
        return None, results_df

    best_row = results_df.loc[
    results_df[
        (results_df["eff"] >= min_eff) &
        (results_df["pur"] >= min_pur)
    ]["eff_x_pur"].idxmax()
    ]
    
    best = {
        "cuts": best_row[cut_labels].values,
        "eff": best_row["eff"],
        "pur": best_row["pur"],
        "eff_x_pur": best_row["eff_x_pur"],
    }

    # -------------------------
    # Print summary
    # -------------------------
    print("\nBest cut combination:")
    for i, val in enumerate(best["cuts"]):
        print(f"  cut_{i}: {val:.3f} ({cut_sign[i]})")

    print(f"  Efficiency : {best['eff']:.4f}")
    print(f"  Purity     : {best['pur']:.4f}")
    print(f"  Eff × Pur  : {best['eff_x_pur']:.4f}")

    return best, results_df

def plot_2d_cut_individual_manual(
    signal_df,
    bkg_df,
    col1,
    col2,
    best,
    cutnum1,
    cutnum2,
    save_dir=None,      # New parameter
    xlabel="Variable 1",
    ylabel="Variable 2",
    signal_label="Signal",
    bkg_label="Background",
    xlim=None,
    ylim=None,
    bins=50,
    cmap="plasma"
):
    """
    Produce and save 2D plots (Signal, Background, and Combined) 
    using manual cut optimization results.
    """
    import os
    
    # Create directory if it doesn't exist
    if save_dir and not os.path.exists(save_dir):
        os.makedirs(save_dir)

    cut1_best = best["cuts"][cutnum1]
    cut2_best = best["cuts"][cutnum2]
    
    # -------------------------
    # Drop NaNs
    # -------------------------
    sig = signal_df.loc[:, [col1, col2]].dropna()
    bkg = bkg_df.loc[:, [col1, col2]].dropna()
    hist_range = [xlim, ylim] if xlim and ylim else None

    # Helper function to avoid repeating plot logic
    def finalize_and_save(filename):
        plt.colorbar(label="Density")
        plt.axvline(cut1_best, color="black", linestyle="--")
        plt.axhline(cut2_best, color="black", linestyle="--")
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.tight_layout()
        if save_dir:
            path = os.path.join(save_dir, filename)
            plt.savefig(path, format='pdf', dpi=300)
            print(f"Saved: {path}")
        plt.show()

    # -------------------------
    # 1️⃣ Signal 2D histogram
    # -------------------------
    plt.figure(figsize=(7, 6))
    plt.hist2d(sig.loc[:, col1], sig.loc[:, col2], bins=bins, range=hist_range, cmap=cmap)
    plt.title(signal_label)
    finalize_and_save("signal.pdf")

    # -------------------------
    # 2️⃣ Background 2D histogram
    # -------------------------
    plt.figure(figsize=(7, 6))
    plt.hist2d(bkg.loc[:, col1], bkg.loc[:, col2], bins=bins, range=hist_range, cmap=cmap)
    plt.title(bkg_label)
    finalize_and_save("bkg.pdf")

    # -------------------------
    # 3️⃣ Signal + Background
    # -------------------------
    all_x = np.concatenate([sig.loc[:, col1], bkg.loc[:, col1]])
    all_y = np.concatenate([sig.loc[:, col2], bkg.loc[:, col2]])

    plt.figure(figsize=(7, 6))
    plt.hist2d(all_x, all_y, bins=bins, range=hist_range, cmap=cmap)
    plt.title(signal_label + " + " + bkg_label)
    finalize_and_save("signal_bkg.pdf")

def plot_2d_cut_metric_heatmap(
    signal_df,
    bkg_df,
    cutnum1,
    cutnum2,
    best,
    results_df,
    xlabel="Variable 1",
    ylabel="Variable 2",
    signal_label="Signal",
    bkg_label="Background",
    xlim=None,
    ylim=None,
    bins=50,
    cmap="plasma"   # ROOT-like (similar to kSunset)
):

     # -------------------------
    # Best cuts
    # -------------------------
    cut1_best = best["cuts"][cutnum1]
    cut2_best = best["cuts"][cutnum2]

    
    # -------------------------
    # 5️⃣ Metric heatmaps
    # -------------------------
    label_1 = "cut_" + str(cutnum1)
    label_2 = "cut_" + str(cutnum2)
    cut1_vals = np.sort(results_df[label_1].unique())
    cut2_vals = np.sort(results_df[label_2].unique())

    def make_grid(metric):
        grid = np.full((len(cut2_vals), len(cut1_vals)), np.nan)
        for _, r in results_df.iterrows():
            i = np.where(cut1_vals == r[label_1])[0][0]
            j = np.where(cut2_vals == r[label_2])[0][0]
            grid[j, i] = r[metric]
        return grid

    for metric, title in [
        ("eff", "Efficiency"),
        ("pur", "Purity"),
        ("eff_x_pur", "Eff × Pur"),
    ]:
        grid = make_grid(metric)

        plt.figure(figsize=(7, 6))
        plt.imshow(
            grid,
            origin="lower",
            extent=[cut1_vals.min(), cut1_vals.max(),
                    cut2_vals.min(), cut2_vals.max()],
            aspect="auto",
            cmap=cmap.reversed()
        )
        plt.colorbar(label=title)
        plt.axvline(cut1_best, color="black", linestyle="--", linewidth=2)
        plt.axhline(cut2_best, color="black", linestyle="--", linewidth=2)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(f"{title} for all cut combinations")
        plt.tight_layout()
        plt.show()

def optimize_cut_accuracy(
    signal_df,
    bkg_df,
    column,
    cut_type=">",
    xlabel="Score",
    title="Cut optimization (Accuracy)",
    signal_name=r"$\nu$",
    bkg_name="Background",
    xlim=None,
    nbins=50,
    legend_loc="upper right",
    cut_unit="",
    normalize_hist=False
):
    """
    Optimize a 1D cut by maximizing Accuracy.
    Accuracy = (Signal_Passing + Bkg_Failing) / Total_Events
    """

    # -------------------------
    # Extract scores
    # -------------------------
    sig_scores = signal_df[column].dropna()
    bkg_scores = bkg_df[column].dropna()

    y_true = np.concatenate([
        np.ones(len(sig_scores)),
        np.zeros(len(bkg_scores)),
    ])

    scores = np.concatenate([
        sig_scores.values,
        bkg_scores.values,
    ])

    # -------------------------
    # Handle cut direction
    # -------------------------
    if cut_type == "<":
        scores_for_roc = -scores
        cut_label = f"{xlabel} < cut"
    elif cut_type == ">":
        scores_for_roc = scores
        cut_label = f"{xlabel} > cut"
    else:
        raise ValueError("cut_type must be '>' or '<'")

    # -------------------------
    # ROC and Statistics
    # -------------------------
    fpr, tpr, thresholds = roc_curve(y_true, scores_for_roc)
    cuts = -thresholds if cut_type == "<" else thresholds

    valid = np.isfinite(cuts)
    cuts = cuts[valid]
    tpr = tpr[valid]
    fpr = fpr[valid]

    n_signal = len(sig_scores)
    n_bkg = len(bkg_scores)
    n_total = n_signal + n_bkg

    # -------------------------
    # Accuracy Calculation
    # -------------------------
    # Accuracy = (True Positives + True Negatives) / Total
    accuracy = (tpr * n_signal + (1 - fpr) * n_bkg) / n_total

    df = pd.DataFrame({
        "cut": cuts,
        "accuracy": accuracy,
        "eff": tpr,
    })

    best = df.loc[df["accuracy"].idxmax()]

    print("\nBest cut (Accuracy Optimization):")
    print(f"  Cut logic   : {cut_label}")
    print(f"  Cut value   : {best['cut']:.4f} {cut_unit}")
    print(f"  Accuracy    : {best['accuracy']:.4f}")
    print(f"  Efficiency  : {best['eff']:.4f}")

    # -------------------------
    # Plot
    # -------------------------
    fig, ax1 = plt.subplots(figsize=(9, 6))
    hist_range = xlim if xlim is not None else (min(scores), max(scores))

    if normalize_hist:
        sig_weights = np.ones_like(sig_scores) / len(sig_scores)
        bkg_weights = np.ones_like(bkg_scores) / len(bkg_scores)
        ylabel = "Normalized entries"
    else:
        sig_weights = np.ones_like(sig_scores)
        bkg_weights = np.ones_like(bkg_scores)
        ylabel = "Events"

    # Histograms
    plot_hist_with_outline(ax1, scores=sig_scores, bins=nbins, range=hist_range,
                           weights=sig_weights, color=COLORS[0], alpha_fill=0.3, label=signal_name)
    
    plot_hist_with_outline(ax1, scores=bkg_scores, bins=nbins, range=hist_range,
                           weights=bkg_weights, color=COLORS[1], alpha_fill=0.3, label=bkg_name)

    # Best cut vertical line
    cut_line = ax1.axvline(best["cut"], color="black", linestyle="--", linewidth=2,
                label=f"Best cut = {best['cut']:.3f}")

    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)

    # Accuracy on secondary axis (Changed to Green)
    ax2 = ax1.twinx()
    acc_line, = ax2.plot(df["cut"], df["accuracy"], label="Accuracy", color="green", linewidth=3)
    ax2.set_ylabel("Accuracy")
    ax2.set_ylim(0, 1.05)

    # Combine legends from both axes
    hist_legend_handles = make_hist_legend_from_ax(ax1, alpha_fill=0.3, linewidth=2.0)
    all_handles = hist_legend_handles + [cut_line, acc_line]
    all_labels = [h.get_label() for h in all_handles]
    
    ax1.legend(handles=all_handles, labels=all_labels, loc=legend_loc)

    if xlim is not None:
        ax1.set_xlim(xlim)

    plt.title(title)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()

    return df, best, fig