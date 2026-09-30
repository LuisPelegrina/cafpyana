import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import ROOT
from ROOT import TF1, TH1D
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy.signal import fftconvolve, savgol_filter
from scipy.stats import gaussian_kde

from analysis_village.cc1pi.TLExtensionMethod.GaussianFactorFittingUtils import *


def get_conv_tf1_for_slice(conv_map, plane, rr_center, max_rr):
    """Look up the TF1* in a conv_tf1_map / conv_shifted_tf1_map for the
    residual-range bin closest to rr_center. Bins are 1 cm wide, centered
    at i_rr + 0.5 for i_rr = 0 .. max_rr-1 (see get_conv_function_map).
    Returns None if plane is missing or the index falls outside the map.
    """
    if plane not in conv_map:
        return None
    i_rr = int(np.floor(rr_center))  # rr_center == i_rr + 0.5 for an aligned bin
    if i_rr < 0 or i_rr >= max_rr or i_rr >= len(conv_map[plane]):
        return None
    return conv_map[plane][i_rr]


class InterpFunc:
    """Helper wrapper to expose an .Eval(x) method for 1D arrays,
    allowing interpolated arrays/histograms to be passed to robust_max_x.
    """

    def __init__(self, x, y):
        self._fn = interp1d(x, y, kind="linear", bounds_error=False, fill_value=0.0)

    def Eval(self, x):
        return float(self._fn(x))


# ===========================================================================
# 9. Shift diagnostics: (a) the mode-shift dx between the theoretical PDF's
#    MPV and its Gaussian-smeared convolution's MPV, and (b) the offset
#    between the raw histogram's peak and the theoretical MPV -- both
#    plotted vs rr_center, one point per rr slice.
# ===========================================================================


def gaussian_kernel(x, sigma):
    """Zero-mean 1D Gaussian kernel."""
    return np.exp(-0.5 * (x / sigma) ** 2) / (sigma * np.sqrt(2.0 * np.pi))


def robust_hist_max_kde(values, xmin, xmax, n_grid=2000, bw_method=None):
    """Alternative to robust_hist_max: find the peak (mode) via a Gaussian
    KDE of the raw sample instead of a smoothed histogram. Less prone to
    the mode-shift bias Savitzky-Golay smoothing introduces on skewed
    (Landau-like) peaks, since it doesn't rely on a locally-symmetric
    polynomial fit.

    bw_method: passed to scipy.stats.gaussian_kde (e.g. a float scale
    factor, or 'scott'/'silverman'). None uses scipy's default (Scott's
    rule) -- worth tuning if the peak still looks over/under-smoothed.
    """
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    if len(vals) < 5:
        return np.nan

    kde = gaussian_kde(vals, bw_method=bw_method)
    x_grid = np.linspace(xmin, xmax, n_grid)
    density = kde(x_grid)
    i_max = int(np.argmax(density))

    # Parabolic refinement around the grid maximum, same idea as
    # robust_hist_max's sub-bin interpolation
    if 0 < i_max < len(density) - 1:
        y0, y1, y2 = density[i_max - 1], density[i_max], density[i_max + 1]
        denom = y0 - 2.0 * y1 + y2
        step = x_grid[1] - x_grid[0]
        delta = 0.5 * (y0 - y2) / denom if denom != 0 else 0.0
        delta = np.clip(delta, -1.0, 1.0)
        return x_grid[i_max] + delta * step
    return x_grid[i_max]


def convolved_pdf_curve(x_grid, pdf_vals, dx, sigma):
    """Theory PDF (x) zero-mean Gaussian(sigma), evaluated back on x_grid
    itself -- no amplitude scaling, since only the peak *location* is
    needed here, not the height.
    """
    x_centered = x_grid - x_grid[len(x_grid) // 2]
    kernel = gaussian_kernel(x_centered, sigma)
    kernel = kernel / (kernel.sum() * dx)
    conv = fftconvolve(pdf_vals, kernel, mode="same") * dx
    return conv


def conv_mpv_from_grid(x_grid, conv_vals, xmin=0.0, xmax=10.0):
    """Grid-search peak location of a convolved curve, restricted to
    [xmin, xmax] (mirrors robust_max_x_py's search window).
    """
    mask = (x_grid >= xmin) & (x_grid <= xmax)
    idx = np.argmax(conv_vals[mask])
    return x_grid[mask][idx]


def compute_shift_diagnostics(
    df,
    plane,
    tpc,
    hfit,
    pdg,
    fit_params,
    dedx_col="dedx",
    rr_min=3.0,
    rr_max=40.0,
    rr_bin_width=1.0,
    pitch=0.55,
    mass=None,
    min_entries=200,
    hist_nbins=1000,
    hist_xmin=0.0,
    hist_xmax=10.0,
    savgol_window=51,
    savgol_poly=3,
    peak_method="savgol",  # "savgol" or "kde"
    kde_bw_method=None,
    verbose=False,
):
    popt, _ = fit_params.get((plane, tpc), (None, None))
    if popt is None:
        raise ValueError(
            f"No sigma_G(MPV) power-law fit available for (plane={plane}, tpc={tpc}); "
            "run analyze(...)/analyze_theoretical(...) first."
        )

    sl_all = df if tpc == -1 else df[df["tpc"] == tpc]
    edges = make_rr_edges(rr_min, rr_max, rr_bin_width)
    rows = []

    for lo, hi in zip(edges[:-1], edges[1:]):
        sl = sl_all[(sl_all["rr"] >= lo) & (sl_all["rr"] < hi)]
        if len(sl) < min_entries:
            continue
        rr_center = 0.5 * (lo + hi)

        pdf = build_theoretical_pdf(hfit, pdg, rr_center, pitch, mass=mass)
        theory_mpv = robust_max_x_py(pdf, 0.0, 10.0, 2000)
        x_grid, pdf_vals, dx = build_pdf_grid(pdf)

        # Predicted smearing sigma at this slice's theory MPV, and the
        # resulting convolution's MPV
        predicted_sigma_g = power_law(theory_mpv, *popt)
        conv_vals = convolved_pdf_curve(x_grid, pdf_vals, dx, predicted_sigma_g)
        conv_mpv = conv_mpv_from_grid(x_grid, conv_vals, xmin=0.0, xmax=10.0)
        dx_shift = conv_mpv - theory_mpv  # was theory_mpv - conv_mpv

        if peak_method == "kde":
            hist_max = robust_hist_max_kde(
                sl[dedx_col].to_numpy(),
                hist_xmin,
                hist_xmax,
                bw_method=kde_bw_method,
            )
        else:
            hist_max, _, _, _ = robust_hist_max(
                sl[dedx_col].to_numpy(),
                hist_xmin,
                hist_xmax,
                nbins=hist_nbins,
                savgol_window=savgol_window,
                savgol_poly=savgol_poly,
            )
        diff_hist_theory = hist_max - theory_mpv

        if verbose:
            print(
                f"  rr={rr_center:.2f}: theory_mpv={theory_mpv:.3f}, "
                f"sigma_G={predicted_sigma_g:.3f}, conv_mpv={conv_mpv:.3f}, "
                f"dx_shift={dx_shift:.3f}, hist_max={hist_max:.3f}, "
                f"diff_hist_theory={diff_hist_theory:.3f}"
            )

        rows.append(
            dict(
                plane=plane,
                tpc=tpc,
                rr_center=rr_center,
                n_hits=len(sl),
                theory_mpv=theory_mpv,
                predicted_sigma_g=predicted_sigma_g,
                conv_mpv=conv_mpv,
                dx_shift=dx_shift,
                hist_max=hist_max,
                diff_hist_theory=diff_hist_theory,
            )
        )

    return pd.DataFrame(rows)


def fit_diff_hist_theory(
    diag_df, x_col="rr_center", y_col="diff_hist_theory", p0=None
):
    """Fits exp_decay_plateau to diff_hist_theory vs rr_center.
    Returns (popt, perr) or (None, None) if the fit fails or there
    aren't enough points.
    """
    xs = diag_df[x_col].to_numpy()
    ys = diag_df[y_col].to_numpy()
    mask = np.isfinite(xs) & np.isfinite(ys)
    xs, ys = xs[mask], ys[mask]

    if len(xs) < 4:
        return None, None

    if p0 is None:
        # Rough starting guess: plateau ~ value at largest rr,
        # amplitude ~ (value at smallest rr) - plateau, decay length
        # ~ a third of the rr range
        order = np.argsort(xs)
        a0 = ys[order][-1]
        b0 = ys[order][0] - a0
        c0 = max((xs.max() - xs.min()) / 3.0, 1e-3)
        p0 = [a0, b0, c0]

    try:
        popt, pcov = curve_fit(exp_decay_plateau, xs, ys, p0=p0, maxfev=20000)
        perr = np.sqrt(np.diag(pcov))
        return popt, perr
    except Exception:
        return None, None


def robust_hist_max(values, xmin, xmax, nbins=1000, savgol_window=51, savgol_poly=3):
    """Find the peak (mode) of a sample via a finely-binned histogram,
    smoothed with a Savitzky-Golay filter, then refined to sub-bin
    precision with parabolic interpolation around the smoothed maximum.

    Raw argmax on an O(1000)-bin histogram is dominated by Poisson noise
    bin-to-bin; smoothing first (rather than just coarsening the binning)
    keeps the underlying peak shape while suppressing that noise, and the
    parabolic refinement recovers precision the smoothing softens.
    """
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    counts, edges = np.histogram(vals, bins=nbins, range=(xmin, xmax))
    centers = 0.5 * (edges[:-1] + edges[1:])
    bin_width = centers[1] - centers[0]

    # Savgol window must be odd and <= number of bins, and > polyorder
    window = min(savgol_window, len(counts))
    if window % 2 == 0:
        window -= 1
    window = max(
        window, savgol_poly + 2 if (savgol_poly + 2) % 2 else savgol_poly + 3
    )

    if window >= len(counts) or window < 5:
        # Not enough bins to smooth meaningfully -- fall back to raw argmax
        i_max = int(np.argmax(counts))
        return centers[i_max], centers, counts, counts.astype(float)

    smoothed = savgol_filter(
        counts.astype(float), window_length=window, polyorder=savgol_poly
    )
    i_max = int(np.argmax(smoothed))

    # Parabolic (3-point) interpolation around the smoothed peak for
    # sub-bin precision -- fits a parabola through (i-1, i, i+1) and
    # returns the vertex position rather than just the discrete bin center.
    if 0 < i_max < len(smoothed) - 1:
        y0, y1, y2 = smoothed[i_max - 1], smoothed[i_max], smoothed[i_max + 1]
        denom = y0 - 2.0 * y1 + y2
        delta = 0.5 * (y0 - y2) / denom if denom != 0 else 0.0
        delta = np.clip(delta, -1.0, 1.0)
        x_max = centers[i_max] + delta * bin_width
    else:
        x_max = centers[i_max]

    return x_max, centers, counts, smoothed