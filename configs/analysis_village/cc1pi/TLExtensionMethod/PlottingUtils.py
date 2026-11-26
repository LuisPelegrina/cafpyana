
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import ROOT
from ROOT import TH1D, TF1
from analysis_village.cc1pi.TLExtensionMethod.GaussianFactorFittingUtils import *
from analysis_village.cc1pi.TLExtensionMethod.ConvolutionTestUtils import *

def plot_slice_diagnostic_convolution_test(
    df,
    plane,
    tpc,
    rr_min,
    rr_max,
    hfit=None,
    pdg=13,
    mass=None,
    dedx_col="dedx",
    nbins=50,
    hist_xmin=0.0,
    hist_xmax=20.0,
    first_stage_range=(0.0, 10.0),
    conv_tf1_map=None,
    conv_shifted_tf1_map=None,
    langau_map=None,
    conv_max_rr=1000,
    x_limits=None,
    plot_only_theory=False,
):
    # If x_limits is specified, use it to redefine the histogram boundaries
    if x_limits is not None:
        hist_xmin, hist_xmax = x_limits

    sl = df if tpc == -1 else df[df["tpc"] == tpc]
    sl = sl[(sl["rr"] >= rr_min) & (sl["rr"] < rr_max) & (sl["pitch"] <= 2)]
    rr_center = 0.5 * (rr_min + rr_max)

    i_rr_lookup = int(np.floor(rr_center))
    i_rr_lookup = max(0, min(i_rr_lookup, conv_max_rr - 1))
    rr_aligned = i_rr_lookup + 0.5

    # If plot_only_theory is True, disable all convolution/map fits
    if plot_only_theory:
        conv_tf1_map = None
        conv_shifted_tf1_map = None
        langau_map = None

    fig, ax = plt.subplots(figsize=(8, 5.5))
    counts, edges, _ = ax.hist(
        sl[dedx_col].dropna(),
        bins=nbins,
        range=(hist_xmin, hist_xmax),
        histtype="stepfilled",
        alpha=0.3,
        color="steelblue",
        edgecolor="navy",
        label=f"Measured PDF (N={len(sl)})",
    )
    bin_width = (hist_xmax - hist_xmin) / nbins
    x_eval = np.linspace(hist_xmin, hist_xmax, 400)
    data_max = counts.max() if len(counts) else 0.0
    curve_max = data_max

    # Track curves for visible range scaling
    y_theo, y_conv, y_conv_shifted, y_langau_map = None, None, None, None

    hist = th1_from_series(
        sl[dedx_col],
        f"diag_p{plane}_t{tpc}",
        "",
        nbins,
        hist_xmin,
        hist_xmax,
    )

    if hfit is not None:
        mean_pitch = 0.32
        pdf = build_theoretical_pdf(hfit, pdg, rr_aligned, mean_pitch, mass=mass)
        norm = len(sl) * bin_width
        y_theo = np.array([pdf.Eval(xv) * norm for xv in x_eval])
        if y_theo.max() > 0 and data_max > 0:
            y_theo = y_theo / y_theo.max() * data_max
        ax.plot(
            x_eval,
            y_theo,
            color="#785EF0",  # Colorblind-friendly purple
            linestyle=":",
            lw=2.5,
            label="Theoretical PDF",
        )
        curve_max = max(curve_max, y_theo.max())

    if conv_tf1_map is not None:
        f_conv = get_conv_tf1_for_slice(conv_tf1_map, plane, rr_center, conv_max_rr)
        if f_conv is not None:
            y_conv = np.array([f_conv.Eval(xv) for xv in x_eval])
            if y_conv.max() > 0 and data_max > 0:
                y_conv = y_conv / y_conv.max() * data_max
            ax.plot(
                x_eval,
                y_conv,
                color="crimson",
                linestyle="-.",
                lw=2.0,
                label="Convolution",
            )
            curve_max = max(curve_max, y_conv.max())

    if conv_shifted_tf1_map is not None:
        f_conv_shifted = get_conv_tf1_for_slice(
            conv_shifted_tf1_map, plane, rr_center, conv_max_rr
        )
        if f_conv_shifted is not None:
            y_conv_shifted = np.array([f_conv_shifted.Eval(xv) for xv in x_eval])
            if y_conv_shifted.max() > 0 and data_max > 0:
                y_conv_shifted = y_conv_shifted / y_conv_shifted.max() * data_max
            ax.plot(
                x_eval,
                y_conv_shifted,
                color="teal",
                linestyle="-.",
                lw=2.0,
                label="Convolution (shifted)",
            )
            curve_max = max(curve_max, y_conv_shifted.max())

    if langau_map is not None:
        f_langau_map = get_conv_tf1_for_slice(langau_map, plane, rr_center, conv_max_rr)
        if f_langau_map is not None:
            y_langau_map = np.array([f_langau_map.Eval(xv) for xv in x_eval])
            if y_langau_map.max() > 0 and data_max > 0:
                y_langau_map = y_langau_map / y_langau_map.max() * data_max
            ax.plot(
                x_eval,
                y_langau_map,
                color="darkorange",
                linestyle="-.",
                lw=2.0,
                label="LanGau approximation",
            )
            curve_max = max(curve_max, y_langau_map.max())

    # Set axes limits cleanly based on current bin boundaries
    ax.set_xlim(hist_xmin, hist_xmax)
    ax.set_ylim(0, curve_max * 1.15 if curve_max > 0 else 1.0)

    ax.set_xlabel(r"$dE/dx$ [MeV/cm]")
    ax.set_ylabel(f"hits / {bin_width:.2f} MeV/cm")
    
    tpc_str = "TPCs Combined" if tpc == -1 else f"tpc {tpc}"
    ax.set_title(f"plane {plane}, {tpc_str}, {rr_min:g} <= rr < {rr_max:g} cm")

    ax.legend(fontsize=9)
    fig.tight_layout()
    plt.show()
    return fig



def plot_shift_diagnostics(diag_df, plane, tpc, pitch=None, out_prefix="shift_diag", fit_diff_hist_theory_curve=True):
    """Single-panel plot of shift diagnostics vs rr_center."""
    if pitch is not None:
        print(f"theoretical MPV computed at pitch = {pitch:.2f} cm")

    fig, ax = plt.subplots(figsize=(8, 5.5))

    ax.axhline(0, color="gray", lw=1, linestyle=":")
    ax.plot(
        diag_df["rr_center"], diag_df["dx_shift"], "o-",
        color="darkorange", ms=4,
        label=r"conv(theory, $\sigma_G$) $-$ theory  (mode-shift from smearing)",
    )
    ax.plot(
        diag_df["rr_center"], diag_df["diff_hist_theory"], "s-",
        color="mediumpurple", ms=4,
        label=r"hist. max $-$ theory  (data peak vs. unsmeared theory)",
    )

    if fit_diff_hist_theory_curve:
        popt, perr = fit_diff_hist_theory(diag_df)
        if popt is not None:
            a, b, c = popt
            a_e, b_e, c_e = perr
            print(f"diff_hist_theory fit: a={a:.4f}+/-{a_e:.4f}, "
                  f"b={b:.4f}+/-{b_e:.4f}, c={c:.4f}+/-{c_e:.4f}")
            xs_fit = np.linspace(diag_df["rr_center"].min(), diag_df["rr_center"].max(), 200)
            ys_fit = exp_decay_plateau(xs_fit, *popt)
            ax.plot(
                xs_fit, ys_fit, "--", color="mediumpurple", lw=1.8, alpha=0.8,
                label=rf"fit: ${a:.3f} + {b:.3f}\,e^{{-rr/{c:.3f}}}$",
            )
        else:
            print("diff_hist_theory fit failed or too few points -- skipping overlay")

    ax.set_xlabel("rr [cm]")
    ax.set_ylabel(r"$\Delta$ MPV [MeV/cm]")
    tpc_label = "combined" if tpc == -1 else tpc
    title = f"plane {plane}, tpc {tpc_label}"
    if pitch is not None:
        title += f", pitch={pitch:.2f} cm"
    ax.set_title(title)
    ax.legend(fontsize=9)
    ax.grid(True, linestyle=":", alpha=0.5)

    fig.tight_layout()
    fig.savefig(f"{out_prefix}_p{plane}_t{tpc}.pdf")
    plt.show()
    return fig

    
def plot_mpv_vs_pitch(
    hfit,
    pdg,
    rr_center,
    mass=None,
    pitch_min=0.1,
    pitch_max=2.0,
    n_pitch=50,
):
    """Sweep pitch through build_theoretical_pdf and plot the resulting
    theoretical MPV (peak of the unconvolved PDF) as a function of pitch,
    at fixed rr_center. Useful for seeing how sensitive the theory MPV is
    to the pitch value -- e.g. checking whether the hardcoded pitch=0.32
    in get_conv_function_map is close enough to the data's actual pitch
    range to not matter, or whether it's a real source of MPV mismatch.
    """
    pitches = np.linspace(pitch_min, pitch_max, n_pitch)
    mpvs = []
    for p in pitches:
        pdf = build_theoretical_pdf(hfit, pdg, rr_center, p, mass=mass)
        mpvs.append(pdf.GetMaximumX())
    mpvs = np.array(mpvs)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(pitches, mpvs, marker="o", markersize=3, lw=1.5, color="darkorange")
    ax.axvline(0.32, color="gray", linestyle="--", lw=1,
               label="hardcoded C++ pitch (0.32)")
    ax.set_xlabel("pitch [cm]")
    ax.set_ylabel("theory MPV [MeV/cm]")
    ax.set_title(f"Theoretical MPV vs pitch (rr={rr_center:g} cm, pdg={pdg})")
    ax.legend(fontsize=9)
    fig.tight_layout()
    plt.show()
    return fig, pitches, mpvs


def plot_slice_diagnostic_w_pitch(
    df,
    plane,
    tpc,
    rr_min,
    rr_max,
    hfit=None,
    pdg=13,
    mass=None,
    dedx_col="dedx",
    nbins=100,
    hist_xmin=0.0,
    hist_xmax=20.0,
    pitch_min=0.31,
    pitch_max=0.33,
    fit_params=None,
    first_stage_range=(0.0, 10.0),
    conv_tf1_map=None,
    conv_shifted_tf1_map=None,
    langau_map=None,
    conv_max_rr=1000,
    x_limits=None,
    plot_only_theory=False,
):
    # If x_limits is specified, use it to redefine the histogram boundaries
    if x_limits is not None:
        hist_xmin, hist_xmax = x_limits

    # Filter by TPC, Residual Range (rr), and specific Pitch window
    sl = df if tpc == -1 else df[df["tpc"] == tpc]
    sl = sl[
        (sl["rr"] >= rr_min)
        & (sl["rr"] < rr_max)
        & (sl["pitch"] >= pitch_min)
        & (sl["pitch"] <= pitch_max)
    ]
    rr_center = 0.5 * (rr_min + rr_max)

    i_rr_lookup = int(np.floor(rr_center))
    i_rr_lookup = max(0, min(i_rr_lookup, conv_max_rr - 1))
    rr_aligned = i_rr_lookup + 0.5

    # If plot_only_theory is True, disable all convolution/map fits
    if plot_only_theory:
        conv_tf1_map = None
        conv_shifted_tf1_map = None
        langau_map = None

    fig, ax = plt.subplots(figsize=(8, 5.5))
    counts, edges, _ = ax.hist(
        sl[dedx_col].dropna(),
        bins=nbins,
        range=(hist_xmin, hist_xmax),
        histtype="stepfilled",
        alpha=0.3,
        color="steelblue",
        edgecolor="navy",
        label=f"Measured PDF (N={len(sl)})",
    )
    bin_width = (hist_xmax - hist_xmin) / nbins
    x_eval = np.linspace(hist_xmin, hist_xmax, 400)
    data_max = counts.max() if len(counts) else 0.0
    curve_max = data_max

    # Track curves for visible range scaling
    y_theo, y_conv, y_conv_shifted, y_langau_map = None, None, None, None

    hist = th1_from_series(
        sl[dedx_col],
        f"diag_p{plane}_t{tpc}",
        "",
        nbins,
        hist_xmin,
        hist_xmax,
    )

    if hfit is not None:
        mean_pitch = 0.5 * (pitch_min + pitch_max)
        
        # Pass fit_params to build_theoretical_pdf if provided
        pdf_kwargs = {"mass": mass}
        if fit_params is not None:
            pdf_kwargs["fit_params"] = fit_params

        pdf = build_theoretical_pdf(hfit, pdg, rr_aligned, mean_pitch, **pdf_kwargs)
        norm = len(sl) * bin_width
        y_theo = np.array([pdf.Eval(xv) * norm for xv in x_eval])
        if y_theo.max() > 0 and data_max > 0:
            y_theo = y_theo / y_theo.max() * data_max
        ax.plot(
            x_eval,
            y_theo,
            color="#785EF0",  # Colorblind-friendly purple
            linestyle=":",
            lw=2.5,
            label="Theoretical PDF",
        )
        curve_max = max(curve_max, y_theo.max())

    if conv_tf1_map is not None:
        f_conv = get_conv_tf1_for_slice(conv_tf1_map, plane, rr_center, conv_max_rr)
        if f_conv is not None:
            y_conv = np.array([f_conv.Eval(xv) for xv in x_eval])
            if y_conv.max() > 0 and data_max > 0:
                y_conv = y_conv / y_conv.max() * data_max
            ax.plot(
                x_eval,
                y_conv,
                color="crimson",
                linestyle="-.",
                lw=2.0,
                label="Convolution",
            )
            curve_max = max(curve_max, y_conv.max())

    if conv_shifted_tf1_map is not None:
        f_conv_shifted = get_conv_tf1_for_slice(
            conv_shifted_tf1_map, plane, rr_center, conv_max_rr
        )
        if f_conv_shifted is not None:
            y_conv_shifted = np.array([f_conv_shifted.Eval(xv) for xv in x_eval])
            if y_conv_shifted.max() > 0 and data_max > 0:
                y_conv_shifted = y_conv_shifted / y_conv_shifted.max() * data_max
            ax.plot(
                x_eval,
                y_conv_shifted,
                color="teal",
                linestyle="-.",
                lw=2.0,
                label="Convolution (shifted)",
            )
            curve_max = max(curve_max, y_conv_shifted.max())

    if langau_map is not None:
        f_langau_map = get_conv_tf1_for_slice(langau_map, plane, rr_center, conv_max_rr)
        if f_langau_map is not None:
            y_langau_map = np.array([f_langau_map.Eval(xv) for xv in x_eval])
            if y_langau_map.max() > 0 and data_max > 0:
                y_langau_map = y_langau_map / y_langau_map.max() * data_max
            ax.plot(
                x_eval,
                y_langau_map,
                color="darkorange",
                linestyle="-.",
                lw=2.0,
                label="LanGau approximation",
            )
            curve_max = max(curve_max, y_langau_map.max())

    # Set axes limits cleanly based on current bin boundaries
    ax.set_xlim(hist_xmin, hist_xmax)
    ax.set_ylim(0, curve_max * 1.15 if curve_max > 0 else 1.0)

    ax.set_xlabel(r"$dE/dx$ [MeV/cm]")
    ax.set_ylabel(f"hits / {bin_width:.2f} MeV/cm")

    tpc_str = "TPCs Combined" if tpc == -1 else f"tpc {tpc}"
    ax.set_title(
        f"plane {plane}, {tpc_str}, {rr_min:g} <= rr < {rr_max:g} cm, "
        f"{pitch_min:g} <= pitch <= {pitch_max:g} cm"
    )

    ax.legend(fontsize=9)
    fig.tight_layout()
    plt.show()
    return fig


def plot_slice_diagnostic_convolution_test_mc_vs_data(
    mc_df,
    data_df,
    plane,
    tpc,
    rr_min,
    rr_max,
    conv_tf1_map=None,             # MC, unshifted
    conv_shifted_tf1_map=None,     # MC, shifted
    conv_tf1_map_data=None,        # data, unshifted
    conv_shifted_tf1_map_data=None,# data, shifted
    dedx_col="dedx",
    nbins=100,
    hist_xmin=0.0,
    hist_xmax=20.0,
    conv_max_rr=500,
    x_limits=None,
    target_height=1.0,
):
    """Overlay MC (blue) and data (red) dE/dx slice histograms for the same
    rr window, both rescaled to share the same peak height (target_height),
    each with its unshifted (dashed) and shifted (solid) convolution TF1
    curves drawn in the matching color and normalized to that same height."""

    def get_slice(df):
        sl = df if tpc == -1 else df[df["tpc"] == tpc]
        return sl[(sl["rr"] >= rr_min) & (sl["rr"] < rr_max) & (sl["pitch"] <= 2)]

    sl_mc = get_slice(mc_df)
    sl_data = get_slice(data_df)

    rr_center = 0.5 * (rr_min + rr_max)
    i_rr_lookup = int(np.floor(rr_center))
    i_rr_lookup = max(0, min(i_rr_lookup, conv_max_rr - 1))
    rr_aligned = i_rr_lookup + 0.5

    fig, ax = plt.subplots(figsize=(8, 5.5))
    bin_edges = np.linspace(hist_xmin, hist_xmax, nbins + 1)
    bin_width = (hist_xmax - hist_xmin) / nbins
    x_eval = np.linspace(hist_xmin, hist_xmax, 400)

    mc_vals = sl_mc[dedx_col].dropna().values
    data_vals = sl_data[dedx_col].dropna().values

    # --- Precompute raw counts to find each histogram's own peak ---
    raw_counts_mc, _ = np.histogram(mc_vals, bins=bin_edges)
    raw_counts_data, _ = np.histogram(data_vals, bins=bin_edges)
    raw_peak_mc = raw_counts_mc.max() if len(raw_counts_mc) else 0.0
    raw_peak_data = raw_counts_data.max() if len(raw_counts_data) else 0.0

    # --- Per-entry weights so each histogram's peak lands at target_height ---
    w_mc = np.full(len(mc_vals), target_height / raw_peak_mc) if raw_peak_mc > 0 else np.ones(len(mc_vals))
    w_data = np.full(len(data_vals), target_height / raw_peak_data) if raw_peak_data > 0 else np.ones(len(data_vals))

    # --- Histograms (now both peak at target_height) ---
    counts_mc, edges, _ = ax.hist(
        mc_vals,
        bins=bin_edges,
        weights=w_mc,
        histtype="stepfilled",
        alpha=0.25,
        color="steelblue",
        edgecolor="navy",
        label=f"MC (N={len(sl_mc)})",
    )
    counts_data, _, _ = ax.hist(
        data_vals,
        bins=bin_edges,
        weights=w_data,
        histtype="stepfilled",
        alpha=0.25,
        color="salmon",
        edgecolor="darkred",
        label=f"Data (N={len(sl_data)})",
    )

    curve_max = target_height
    all_y_curves = []  # for zoom rescaling later

    def draw_conv_curve(tf1_map, color, linestyle, label_prefix):
        if tf1_map is None:
            return None
        f_conv = get_conv_tf1_for_slice(tf1_map, plane, rr_center, conv_max_rr)
        if f_conv is None:
            return None
        y = np.array([f_conv.Eval(xv) for xv in x_eval])
        if y.max() > 0:
            y = y / y.max() * target_height
        mpv = f_conv.GetMaximumX(0.0, 20.0)
        ax.plot(
            x_eval, y,
            color=color, linestyle=linestyle, lw=2.0,
            label=f"{label_prefix} (MPV={mpv:.2f})",
        )
        all_y_curves.append(y)
        return y

    # MC: unshifted dashed, shifted solid, blue
    draw_conv_curve(conv_tf1_map, "steelblue", "--", "MC conv (unshifted)")
    draw_conv_curve(conv_shifted_tf1_map, "steelblue", "-", "MC conv (shifted)")

    # Data: unshifted dashed, shifted solid, red
    draw_conv_curve(conv_tf1_map_data, "crimson", "--", "Data conv (unshifted)")
    draw_conv_curve(conv_shifted_tf1_map_data, "crimson", "-", "Data conv (shifted)")

    curve_max = max([curve_max] + [y.max() for y in all_y_curves]) if all_y_curves else curve_max

    # --- Zoom limits & dynamic y-rescaling ---
    bin_centers = 0.5 * (edges[:-1] + edges[1:])
    if x_limits is not None:
        ax.set_xlim(x_limits)
        mask = (x_eval >= x_limits[0]) & (x_eval <= x_limits[1])
        hist_mask = (bin_centers >= x_limits[0]) & (bin_centers <= x_limits[1])

        visible_max = 0.0
        if len(counts_mc) and np.any(hist_mask):
            visible_max = max(visible_max, counts_mc[hist_mask].max())
        if len(counts_data) and np.any(hist_mask):
            visible_max = max(visible_max, counts_data[hist_mask].max())
        for y_arr in all_y_curves:
            if len(y_arr[mask]) > 0:
                visible_max = max(visible_max, y_arr[mask].max())

        ax.set_ylim(0, visible_max * 1.15 if visible_max > 0 else 1.0)
    else:
        ax.set_ylim(0, curve_max * 1.15 if curve_max > 0 else 1.0)

    ax.set_xlabel(r"$dE/dx$ [MeV/cm]")
    ax.set_ylabel("normalized (peak = 1)")
    tpc_str = "TPCs Combined" if tpc == -1 else f"tpc {tpc}"
    ax.set_title(f"plane {plane}, {tpc_str}, {rr_min:g} <= rr < {rr_max:g} cm")
    ax.legend(fontsize=9)
    fig.tight_layout()
    plt.show()
    return fig