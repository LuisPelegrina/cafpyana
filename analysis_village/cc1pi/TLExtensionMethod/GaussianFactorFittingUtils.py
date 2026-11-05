"""
langau_rr_fit.py
=================
Python/PyROOT port of the ROOT macro `FitPlotWholeDataset.C`.

Instead of building TH2D(rr, dEdx) histograms from a ROOT file and slicing
them bin-by-bin, this version works directly on per-hit pandas dataframes
(e.g. hit_dfs = [mc_bnb_hit0_df, mc_bnb_hit1_df, mc_bnb_hit2_df], one
dataframe per wire plane, each carrying a 'tpc' column for TPC0/TPC1).
For each (plane, tpc) pair the hits are sliced in residual range (rr),
and the dE/dx distribution of each slice is fit with the same two-stage
Landau-convoluted-with-Gaussian ("Langau") fit as the original macro.
The resulting Gaussian-smearing sigma is then fit vs. MPV with the same
power law used by `f_pol2` in the macro.

Requirements: numpy, pandas, matplotlib, scipy, and a PyROOT build of ROOT
(needed for TH1D/TF1 and TMath::Landau/Gaus, which is the actual fitting
engine -- reimplementing the Landau/Gaussian convolution fit in pure numpy
would give a different, unvalidated minimizer).

-----------------------------------------------------------------------
Mapping from the original macro to this module
-----------------------------------------------------------------------
FitPlotWholeDataset.C construct              -> here
--------------------------------------------  ------------------------------
langaufun                                     -> langau_cpp (declared via
                                                  ROOT.gInterpreter.Declare)
langaufit                                     -> _langau_fit_once
broad fit, then MPV-centered refit            -> langau_fit_two_stage
TH2D "slice_xbin_%d" loop over ix             -> fit_rr_slices (slices the
                                                  hit-level dataframe by rr
                                                  directly, no TH2D needed)
gr_MPV_gaussian_smearing[i] + f_pol2[i] fit   -> per-(plane, tpc) DataFrame
                                                  + fit_sigma_vs_mpv
canvas c3 (all 6 planes/tpc overlay)          -> plot_all_planes
canvas c4 (3 subplots, tpc0 vs tpc1)          -> plot_by_plane
PhysdEdx / Hypfit "theoretical" MPV           -> load_physics_classes +
                                                  make_theoretical_mpv_func
                                                  (see note below)

Note on the "theoretical" MPV
------------------------------
The original macro plots sigma_G against a *theoretical* MPV coming from
`hfit.map_PhysdEdx[PDG]` (Bethe-Bloch mean dE/dx, Landau's xi, Wmax, and a
KE-from-range spline) rather than the fitted MPV of the slice. Since
Hypfit.h/.cpp and PhysdEdx.h/.cpp are available next to this script, this
module loads the *real* classes straight into the ROOT interpreter (the
same way the macro did with `#include "PhysdEdx.h"`/`"PhysdEdx.cpp"`) via
`load_physics_classes`, then reproduces the macro's
KEFromRangeSpline/Landau_xi/Get_Wmax/meandEdx/dEdx_PDF_fuction block in
`make_theoretical_mpv_func`, rather than approximating the physics.
If those files aren't available in a given run, `theoretical_mpv_func`
can simply be left as `None` and the *fitted* MPV of each slice is used
instead (mirrors the macro's commented-out "Doing it with reco" branch).
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import ROOT
from ROOT import TH1D, TF1

# ===========================================================================
# 1. Landau (x) Gaussian convolution -- ported 1:1 from `langaufun`
#    Parameter order matches the original macro:
#       p[0] = Width   (Landau scale)
#       p[1] = MPV     (Landau most probable value)
#       p[2] = Area
#       p[3] = GSigma  (Gaussian smearing sigma)
# ===========================================================================
ROOT.gInterpreter.Declare(r"""
#include "TMath.h"

double langau_cpp(double *x, double *p) {
    const double width = p[0];
    const double mpv   = p[1];
    const double area  = p[2];
    const double sigG  = p[3];

    if (width <= 0 || sigG <= 0) return 0.0;

    const double invsq2pi = 0.398942280401;
    const double np = 500.0;   // integration slices, same as the macro
    const double sc = 5.0;     // convolution extends +/- sc * sigG
    const double xx = x[0];

    double xlow = xx - sc * sigG;
    double xupp = xx + sc * sigG;
    double step = (xupp - xlow) / np;

    double sum = 0.0;
    for (double i = 1.0; i <= np / 2.0; i += 1.0) {
        double t1 = xlow + (i - 0.5) * step;
        double fland1 = TMath::Landau(t1, mpv, width) / width;
        sum += fland1 * TMath::Gaus(xx, t1, sigG);

        double t2 = xupp - (i - 0.5) * step;
        double fland2 = TMath::Landau(t2, mpv, width) / width;
        sum += fland2 * TMath::Gaus(xx, t2, sigG);
    }
    return area * step * sum * invsq2pi / sigG;
}
""")

# Slice fits. The fit models and fit functions (power_law, exp_decay_plateau, exp_decay_linear,
# power_law_linear, SHIFT_MODELS, SHIFT_MODEL, eval_shift_fit, shift_fit_params,
# fit_sigma_vs_mpv, fit_shift_vs_rr) come from the "Fit models and fit functions" cell and are
# no longer redefined here, so the order in which the two cells are run does not matter.
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import ROOT
from ROOT import TF1, TH1D
from scipy.optimize import curve_fit
from scipy.signal import fftconvolve


def load_physics_classes():
    base_dir = (
        "/home/lpelegri/cafpyana/analysis_village/cc1pi/TLExtensionMethod"
    )

    # Compile and load libraries dynamically
    res1 = ROOT.gSystem.CompileMacro(
        os.path.join(base_dir, "PhysdEdx.cpp"), "k"
    )
    res2 = ROOT.gSystem.CompileMacro(os.path.join(base_dir, "Hypfit.cpp"), "k")

    if res1 != 1 or res2 != 1:
        raise RuntimeError(
            "C++ compilation failed! Check stdout for syntax or redefinition errors."
        )

    return ROOT.Hypfit()


def robust_max_x(f, xmin, xmax, n_scan=10000):
    """Python port of Hypfit::robust_max_x -- identical grid-scan +
    parabolic-refinement algorithm, so comparisons against C++-derived
    peaks (PDF_max, f0_max_x) aren't confounded by using a different
    optimizer (ROOT's GetMaximumX uses Brent's method internally, which
    can converge to a slightly different x than a grid+parabola search
    even on the same function/domain)."""
    step = (xmax - xmin) / (n_scan - 1)
    ys = np.array([f.Eval(xmin + i * step) for i in range(n_scan)])
    best_i = int(np.argmax(ys))
    best_x = xmin + best_i * step
    if best_i <= 0 or best_i >= n_scan - 1:
        return best_x
    y0, y1, y2 = ys[best_i - 1], ys[best_i], ys[best_i + 1]
    denom = y0 - 2.0 * y1 + y2
    if denom == 0:
        return best_x
    delta = 0.5 * (y0 - y2) / denom
    return best_x + delta * step


def th1_from_series(values, name, title="", nbins=100, xmin=0.0, xmax=20.0):
    """Build a TH1D (with Sumw2) from a 1D array-like of dE/dx values."""
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    counts, _edges = np.histogram(vals, bins=nbins, range=(xmin, xmax))
    h = TH1D(name, title, nbins, xmin, xmax)
    h.Sumw2()
    for i, c in enumerate(counts, start=1):
        h.SetBinContent(i, float(c))
        h.SetBinError(i, np.sqrt(c) if c > 0 else 1.0)
    return h


# ===========================================================================
# 3b. Load theoretical-MPV functions
# ===========================================================================
PDG_MASS = {
    13: 105.6583755,
    -13: 105.6583755,  # muon
    211: 139.57039,
    -211: 139.57039,  # charged pion
    2212: 938.27208943,  # proton
}


def build_theoretical_pdf(hfit, pdg, rr_center, mean_pitch, mass=None):
    """Reproduces the macro's PDF construction logic and returns the TF1 object."""
    if mass is None:
        mass = PDG_MASS.get(abs(pdg), 105.6583755)

    phys = hfit.map_PhysdEdx[pdg]
    this_KE = phys.KEFromRangeSpline(rr_center)
    gamma = (this_KE / mass) + 1.0
    beta2 = 1.0 - 1.0 / (gamma * gamma)
    this_xi = phys.Landau_xi(this_KE, mean_pitch)
    this_Wmax = phys.Get_Wmax(this_KE)
    this_kappa = this_xi / this_Wmax
    this_dEdx_BB = phys.meandEdx(this_KE)

    pdf = ROOT.TF1("", ROOT.PhysdEdx.dEdx_PDF_function, -10.0, 20.0, 5)
    pdf.SetParameters(this_kappa, beta2, this_xi, this_dEdx_BB, mean_pitch)
    return pdf


# ===========================================================================
# 4. Slice a per-hit dataframe in residual range and fit each slice.
# ===========================================================================
def make_rr_edges(rr_min, rr_max, bin_width):
    return np.arange(rr_min, rr_max + bin_width, bin_width)


# ===========================================================================
# 7. Theoretical-PDF (x) zero-mean-Gaussian smearing-only fit
#
#    Instead of fitting a free 4-parameter Langau (Width, MPV, Area, GSigma)
#    to each rr slice, this fixes the underlying physics PDF entirely
#    (pitch = 0.32, rr = slice center, via build_theoretical_pdf) and only
#    lets a Gaussian smearing sigma (mean fixed at 0) float, plus an
#    amplitude to match the histogram's raw counts. The convolution is done
#    numerically with scipy.signal.fftconvolve instead of TF1Convolution --
#    more robust/easier to inspect than the ROOT FFT convolution, and no
#    ROOT-side TF1 fit machinery is needed for this part.
# ===========================================================================
def build_pdf_grid(pdf, xmin=-5.0, xmax=25.0, n_points=6001):
    """Evaluate the theoretical TF1 PDF once on a fixed fine grid.
    Reused across all sigma trials during the fit -- only the Gaussian
    kernel and the convolution are recomputed per curve_fit iteration."""
    x_grid = np.linspace(xmin, xmax, n_points)
    pdf_vals = np.array([pdf.Eval(x) for x in x_grid])
    dx = x_grid[1] - x_grid[0]
    return x_grid, pdf_vals, dx


def gaussian_kernel(x, sigma):
    """Zero-mean Gaussian sampled on x, normalized to unit area by
    the caller (not here, since spacing dx isn't known inside this fn)."""
    if sigma <= 0:
        kernel = np.zeros_like(x)
        kernel[np.argmin(np.abs(x))] = 1.0
        return kernel
    return np.exp(-0.5 * (x / sigma) ** 2)


def convolved_theoretical_pdf(
    x_query, sigma, shift, amplitude, x_grid, pdf_vals, dx
):
    """3-parameter model: sigma, shift, and amplitude are fit.

    The theoretical PDF is convolved with a unit-area Gaussian kernel of width `sigma`,
    shifted horizontally by `shift`, and scaled by an overall `amplitude` factor.
    """
    x_centered = x_grid - x_grid[len(x_grid) // 2]
    kernel = gaussian_kernel(x_centered, sigma)
    kernel = kernel / (kernel.sum() * dx)  # Normalize Gaussian kernel to unit area
    conv = fftconvolve(pdf_vals, kernel, mode="same") * dx

    # Shift and scale by the free amplitude parameter
    return amplitude * np.interp(x_query - shift, x_grid, conv)


def fit_theoretical_conv_slice(
    bin_centers,
    bin_counts,
    bin_errs,
    x_grid,
    pdf_vals,
    dx,
    mpv_theory,
    fit_range=None,
    sigma0=0.2,
    shift0=0.0,
    amp0=None,
    sigma_bounds=(1e-4, 5.0),
    shift_bounds=(-3.0, 3.0),
    amp_bounds=(0.0, np.inf),
):
    """Fits (sigma, shift, amplitude) using scipy curve_fit.

    Also returns chi2 and ndf of the fit (over the fitted bins) and the
    correlation coefficient between sigma and shift.
    """
    if fit_range is None:
        fit_range = (0.6 * mpv_theory, 1.4 * mpv_theory)
    mask = (
        (bin_centers >= fit_range[0])
        & (bin_centers <= fit_range[1])
        & (bin_counts > 0)
    )
    if mask.sum() < 5:
        return None
    xs, ys, yerr = bin_centers[mask], bin_counts[mask], bin_errs[mask]
    yerr = np.where(yerr > 0, yerr, 1.0)

    # Estimate initial amplitude guess if not provided (data peak / theory peak)
    if amp0 is None:
        peak_pdf = pdf_vals.max() if pdf_vals.max() > 0 else 1.0
        amp0 = ys.max() / peak_pdf

    def model(x, sigma, shift, amplitude):
        return convolved_theoretical_pdf(
            x, sigma, shift, amplitude, x_grid, pdf_vals, dx
        )

    try:
        popt, pcov = curve_fit(
            model,
            xs,
            ys,
            p0=[sigma0, shift0, amp0],
            sigma=yerr,
            absolute_sigma=True,
            bounds=(
                [sigma_bounds[0], shift_bounds[0], amp_bounds[0]],
                [sigma_bounds[1], shift_bounds[1], amp_bounds[1]],
            ),
            maxfev=20000,
        )
    except Exception:
        return None

    perr = np.sqrt(np.diag(pcov))

    resid = (ys - model(xs, *popt)) / yerr
    chi2 = float(np.sum(resid**2))
    ndf = len(xs) - len(popt)

    denom = np.sqrt(pcov[0, 0] * pcov[1, 1])
    corr_sigma_shift = (
        pcov[0, 1] / denom if np.isfinite(denom) and denom > 0 else np.nan
    )

    return dict(
        sigma=popt[0],
        sigma_err=perr[0],
        shift=popt[1],
        shift_err=perr[1],
        amplitude=popt[2],
        amplitude_err=perr[2],
        chi2=chi2,
        ndf=ndf,
        corr_sigma_shift=corr_sigma_shift,
        fit_range=fit_range,
    )


def fit_rr_slices_theoretical(
    df,
    plane,
    tpc,
    hfit,
    pdg,
    dedx_col="dedx",
    rr_min=3.0,
    rr_max=40.0,
    rr_bin_width=1.0,
    hist_nbins=150,
    hist_xmin=0.0,
    hist_xmax=20.0,
    min_entries=30,
    pitch=0.32,
    mass=None,
    sigma0=0.2,
    sigma_bounds=(1e-4, 5.0),
    max_sigma_err=0.2,
    verbose=False,
):
    edges = make_rr_edges(rr_min, rr_max, rr_bin_width)
    rows = []

    for lo, hi in zip(edges[:-1], edges[1:]):
        sl = df[(df["rr"] >= lo) & (df["rr"] < hi)]
        if len(sl) < min_entries:
            continue
        rr_center = 0.5 * (lo + hi)

        pdf = build_theoretical_pdf(hfit, pdg, rr_center, pitch, mass=mass)
        mpv_theory = robust_max_x_py(pdf, 0.0, 10.0, 2000)
        x_grid, pdf_vals, dx = build_pdf_grid(pdf)

        counts, bin_edges = np.histogram(
            sl[dedx_col].to_numpy(),
            bins=hist_nbins,
            range=(hist_xmin, hist_xmax),
        )
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        bin_errs = np.where(counts > 0, np.sqrt(counts), 1.0)

        # No seeding from a data-side peak finder here -- shift0=0.0 lets
        # curve_fit locate the shift purely from the fit itself, with no
        # correction toward a pre-measured data peak.
        fit = fit_theoretical_conv_slice(
            bin_centers,
            counts.astype(float),
            bin_errs,
            x_grid,
            pdf_vals,
            dx,
            mpv_theory,
            sigma0=sigma0,
            shift0=0.0,
            sigma_bounds=sigma_bounds,
        )

        if fit is None or fit["sigma_err"] > max_sigma_err:
            continue

        rows.append(
            dict(
                plane=plane,
                tpc=tpc,
                rr_center=rr_center,
                n_hits=len(sl),
                mpv_x=mpv_theory,
                mpv_theory=mpv_theory,
                gsigma=fit["sigma"],
                gsigma_err=fit["sigma_err"],
                shift=fit["shift"],
                shift_err=fit["shift_err"],
                chi2=fit["chi2"],
                ndf=fit["ndf"],
                corr_sigma_shift=fit["corr_sigma_shift"],
                fit_range_lo=fit["fit_range"][0],
                fit_range_hi=fit["fit_range"][1],
            )
        )

    return pd.DataFrame(rows)


# ===========================================================================
# 5. Langau parameters vs residual range (power law from the fit-models cell)
# ===========================================================================
def fit_param_curve(res_df, param_col, err_col, rr_lo=5.0, rr_hi=None, model="power_law"):
    """Fits a decay curve (power law by default) to a parameter vs rr_center."""
    if len(res_df) < 3:
        return None, None

    mask = res_df["rr_center"] >= rr_lo
    if rr_hi is not None:
        mask &= res_df["rr_center"] <= rr_hi

    xs = res_df.loc[mask, "rr_center"].values
    ys = res_df.loc[mask, param_col].values
    yerrs = res_df.loc[mask, err_col].values

    valid = np.isfinite(xs) & np.isfinite(ys) & np.isfinite(yerrs) & (yerrs > 0) & (xs > 0)
    if np.sum(valid) < 3:
        return None, None

    x_val, y_val, y_err = xs[valid], ys[valid], yerrs[valid]

    if model == "power_law":
        func = power_law
        c_guess = np.percentile(y_val, 10)
        a_guess = max((np.percentile(y_val, 90) - c_guess) * (x_val.min() ** 1.0), 1e-3)
        p0 = [a_guess, 1.0, c_guess]
        bounds = ([0, 0, -np.inf], [np.inf, np.inf, np.inf])
    elif model == "double_exp":
        func = double_exp
        c_guess = np.percentile(y_val, 10)
        span = np.percentile(y_val, 90) - c_guess
        p0 = [span * 0.7, 0.3, span * 0.3, 0.02, c_guess]
        bounds = ([0, 0, 0, 0, -np.inf], [np.inf, np.inf, np.inf, np.inf, np.inf])
    elif model == "inverse":
        func = inverse_law
        c_guess = np.percentile(y_val, 10)
        a_guess = (np.percentile(y_val, 90) - c_guess) * x_val.min()
        p0 = [a_guess, c_guess]
        bounds = ([0, -np.inf], [np.inf, np.inf])
    else:
        raise ValueError(f"Unknown model: {model!r}")

    try:
        popt, pcov = curve_fit(
            func, x_val, y_val, sigma=y_err, absolute_sigma=True,
            p0=p0, bounds=bounds, maxfev=5000,
        )
        perr = np.sqrt(np.diag(pcov))
        return popt, perr
    except Exception:
        return None, None


PARAM_SPECS = [
    ("mpv_reco", "mpv_reco_err", "LanGau MPV", "MPV [MeV/cm]"),
    ("gsigma", "gsigma_err", r"Gaussian Smearing $\sigma_G$", r"$\sigma_G$ [MeV/cm]"),
    ("width", "width_err", r"Landau Scale Width $\xi$", r"$\xi$ [MeV/cm]"),
    ("area", "area_err", "Fit Area", "Area"),
]


def analyze(
    hit_dfs,
    particle="muon",
    dedx_col="dedx",
    rr_max_by_particle=None,
    theoretical_mpv_func=None,
    out_prefix="langau_rr",
    make_summary_plots=True,
    verbose=True,
    fit_model="power_law",
    fit_rr_lo=5.0,
):
    """Analyzes hit DataFrames and performs power-law fits for MPV, gsigma, and width vs RR.

    Returns:
        all_results (dict): {(plane, tpc): DataFrame_per_slice}
        fit_params (dict):  {(plane, tpc): {param_name: (popt, perr)}}
        combined (pd.DataFrame): Master DataFrame of all slices, planes, and fits.
    """
    if rr_max_by_particle is None:
        rr_max_by_particle = {"muon": 80.0, "pion": 60.0, "proton": 60.0}
    rr_max = rr_max_by_particle.get(particle, 80.0)
    all_results = {}
    fit_params = {}
    for plane, df in enumerate(hit_dfs):
        for tpc in (0, 1, -1):
            if verbose:
                tpc_label = "combined" if tpc == -1 else tpc
                print(f"Plane {plane}, TPC {tpc_label}")
            sub = df if tpc == -1 else df[df["tpc"] == tpc]
            res = fit_rr_slices(
                sub,
                plane,
                tpc,
                dedx_col=dedx_col,
                rr_max=rr_max,
                theoretical_mpv_func=theoretical_mpv_func,
                verbose=verbose,
            )
            # Store full results (including the Bragg-peak-adjacent slices)
            all_results[(plane, tpc)] = res
            # Filter out slices too close to the Bragg peak ONLY for the curve-fitting step
            res_fit = res[res["rr_center"] >= fit_rr_lo] if "rr_center" in res.columns else res
            # Fit power law (default) to MPV, GSigma, and Width vs Residual Range
            plane_fits = {}
            if len(res_fit) >= 3:
                for p_col, p_err, _, _ in PARAM_SPECS:
                    if p_col == "area":
                        continue  # Skip fit for Area as requested
                    popt, perr = fit_param_curve(
                        res_fit, p_col, p_err, rr_lo=fit_rr_lo, model=fit_model
                    )
                    plane_fits[p_col] = (popt, perr)
            fit_params[(plane, tpc)] = plane_fits
    non_empty = [
        r.assign(plane=p, tpc=t)
        for (p, t), r in all_results.items()
        if len(r)
    ]
    combined = (
        pd.concat(non_empty, ignore_index=True)
        if non_empty
        else pd.DataFrame()
    )
    if len(combined):
        combined.to_hdf(
            f"{out_prefix}_slices.h5",
            key="fits",
            mode="w",
            format="table",
            complib="blosc",
            complevel=9,
        )
    return all_results, fit_params, combined


def fit_rr_slices(
    df,
    plane,
    tpc,
    dedx_col="dedx",
    rr_min=2.0,
    rr_max=40.0,
    rr_bin_width=1.0,
    hist_nbins=150,
    hist_xmin=0,
    hist_xmax=20.0,
    min_entries=30,
    first_stage_range=(0.5, 20.0),
    theoretical_mpv_func=None,
    max_gsigma_err=0.2,
    min_gsigma_err=0.0001,
    verbose=False,
):
    """Fits individual RR slices and collects MPV, Landau Width, Area, and Gaussian Smearing parameters."""
    edges = make_rr_edges(rr_min, rr_max, rr_bin_width)
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        sl = df[(df["rr"] >= lo) & (df["rr"] < hi)]
        if len(sl) < min_entries:
            continue

        rr_center = 0.5 * (lo + hi)
        hname = f"slice_p{plane}_t{tpc}_rr{rr_center:.2f}"
        hist = th1_from_series(
            sl[dedx_col], hname, "", hist_nbins, hist_xmin, hist_xmax
        )

        fit = langau_fit_two_stage(hist, first_stage_range=first_stage_range)
        if fit is None:
            if verbose:
                print(
                    f"   [plane {plane}, tpc {tpc}] rr={rr_center:.2f}: fit rejected"
                )
            continue

        width, mpv_reco, area, gsigma = fit["pars"]
        w_e, mpv_e, area_e, gsigma_e = fit["errs"]

        if gsigma_e > max_gsigma_err or not np.isfinite(gsigma_e) or gsigma_e < min_gsigma_err:
            if verbose:
                print(
                    f"   [plane {plane}, tpc {tpc}] rr={rr_center:.2f}: dropped due to gsigma_err={gsigma_e:.3f}"
                )
            continue

        if theoretical_mpv_func is not None:
            mean_pitch = sl["pitch"].mean()
            mpv_x = theoretical_mpv_func(rr_center, mean_pitch)
            mpv_x_err = 0.0
        else:
            mpv_x = mpv_reco
            mpv_x_err = mpv_e

        rows.append(
            dict(
                plane=plane,
                tpc=tpc,
                rr_center=rr_center,
                n_hits=len(sl),
                mpv_x=mpv_x,
                mpv_x_err=mpv_x_err,
                mpv_reco=mpv_reco,
                mpv_reco_err=mpv_e,
                width=width,
                width_err=w_e,
                area=area,
                area_err=area_e,
                gsigma=gsigma,
                gsigma_err=gsigma_e,
                chi2=fit["chi2"],
                ndf=fit["ndf"],
                status=fit["status"],
            )
        )

    return pd.DataFrame(rows)


def robust_max_x_py(tf1, xmin, xmax, n_points=2000):
    """Python port of the C++ robust_max_x: scan a grid rather than trust
    TF1::GetMaximumX(), which can get stuck depending on start point."""
    xs = np.linspace(xmin, xmax, n_points)
    vals = np.array([tf1.Eval(x) for x in xs])
    return xs[np.argmax(vals)]


def analyze_theoretical(
    hit_dfs,
    hfit,
    pdg,
    particle="muon",
    dedx_col="dedx",
    rr_max_by_particle=None,
    pitch=0.32,
    mass=None,
    sigma0=0.2,
    sigma_bounds=(1e-4, 5.0),
    max_sigma_err=0.2,
    out_prefix="langau_rr_theoretical",
    verbose=True,
    sigma_fit_min=0.0001,   # drop points with gsigma_err below this
    sigma_fit_max=0.05,     # drop points with gsigma_err above this
    shift_model=None,       # key of SHIFT_MODELS; None -> SHIFT_MODEL
):
    """Theoretical-PDF slice fits for every (plane, tpc), then sigma_G vs MPV
    (fit_sigma_vs_mpv) and shift vs residual range (fit_shift_vs_rr with shift_model).

    Returns:
        all_results (dict):      {(plane, tpc): DataFrame_per_slice}
        sigma_fit_params (dict): {(plane, tpc): (popt, perr)}
        shift_fit_params (dict): {(plane, tpc): (popt, perr)}, all with model shift_model
    """
    if shift_model is None:
        shift_model = SHIFT_MODEL
    if rr_max_by_particle is None:
        rr_max_by_particle = {"muon": 250.0, "pion": 80.0, "proton": 60.0}
    rr_max = rr_max_by_particle.get(particle, 80.0)
    all_results = {}
    sigma_fit_params = {}
    shift_fit_params = {}
    for plane, df in enumerate(hit_dfs):
        for tpc in (0, 1, -1):
            if verbose:
                tpc_label = "combined" if tpc == -1 else tpc
                print(f"[theoretical] Plane {plane}, TPC {tpc_label}")
            sub = df if tpc == -1 else df[df["tpc"] == tpc]
            res = fit_rr_slices_theoretical(
                sub,
                plane,
                tpc,
                hfit=hfit,
                pdg=pdg,
                dedx_col=dedx_col,
                rr_max=rr_max,
                pitch=pitch,
                mass=mass,
                sigma0=sigma0,
                sigma_bounds=sigma_bounds,
                max_sigma_err=max_sigma_err,
                verbose=verbose,
            )

            # Drop degenerate/out-of-range points ONCE, up front. The
            # filtered result -- not the raw `res` -- is what gets used
            # everywhere downstream: both fits, all_results, and the
            # final HDF5 export. Dropped points never reappear anywhere.
            if "gsigma_err" in res.columns:
                mask = (
                    (res["gsigma_err"] >= sigma_fit_min)
                    & (res["gsigma_err"] <= sigma_fit_max)
                )
                n_dropped = len(res) - mask.sum()
                if verbose and n_dropped > 0:
                    print(
                        f"  -> dropped {n_dropped} slice(s) with gsigma_err "
                        f"outside [{sigma_fit_min}, {sigma_fit_max}]"
                    )
                res = res[mask]

            all_results[(plane, tpc)] = res

            # Both fit functions check the number of points themselves and return
            # (None, None) if the fit is not possible or fails; the empty-result guard
            # only protects them from a DataFrame without columns.
            if len(res):
                sigma_fit_params[(plane, tpc)] = fit_sigma_vs_mpv(res)
                shift_fit_params[(plane, tpc)] = fit_shift_vs_rr(
                    res, model=shift_model, verbose=verbose
                )
            else:
                sigma_fit_params[(plane, tpc)] = (None, None)
                shift_fit_params[(plane, tpc)] = (None, None)

    non_empty = [
        r.assign(plane=p, tpc=t)
        for (p, t), r in all_results.items()
        if len(r)
    ]
    combined = (
        pd.concat(non_empty, ignore_index=True)
        if non_empty
        else pd.DataFrame()
    )
    if len(combined):
        combined.to_hdf(
            f"{out_prefix}_slices.h5",
            key="fits",
            mode="w",
            format="table",
            complib="blosc",
            complevel=9,
        )
    return all_results, sigma_fit_params, shift_fit_params


def _langau_fit_once(hist, frange, start, lo, hi, fname, fit_options="RBOSQNL"):
    """Enhanced single-stage Langau fit with ROOT object cleanup and boundary checks."""
    old = ROOT.gROOT.GetListOfFunctions().FindObject(fname)
    if old:
        old.Delete()

    f = ROOT.TF1(fname, ROOT.langau_cpp, frange[0], frange[1], 4)
    f.SetParameters(*start)
    f.SetParNames("Width", "MPV", "Area", "GSigma")

    for i in range(4):
        # Prevent lower bound >= upper bound ROOT crashes
        if lo[i] >= hi[i]:
            lo[i] = start[i] * 0.1 if start[i] > 0 else 0.001
            hi[i] = start[i] * 10.0 if start[i] > 0 else 10.0
        f.SetParLimits(i, lo[i], hi[i])

    fitres = hist.Fit(fname, fit_options)

    pars = [f.GetParameter(i) for i in range(4)]
    errs = [f.GetParError(i) for i in range(4)]
    chi2 = f.GetChisquare()
    ndf = f.GetNDF()
    status = fitres.CovMatrixStatus() if fitres else -1

    return f, pars, errs, chi2, ndf, status


def langau_fit_single_stage(hist, first_stage_range=(0.0, 20.0), max_par_err=1.0):
    """Single-stage Landau-Gaussian fit performed over the range [0.5, max_x],

    where `max_x` is the x-position of the maximum bin in the histogram.

    Rejects the fit if MPV or GSigma errors exceed `max_par_err`.

    Returns a dict(func, pars, errs, chi2, ndf, status), or None if the
    fit fails or parameter errors are too large.
    """
    max_x = hist.GetBinCenter(hist.GetMaximumBin())
    bin_width = hist.GetBinWidth(1)

    sv = [0.1, max_x, hist.Integral() * 0.05 * bin_width, 0.2]
    lo = [0.01 * v if v != 0 else 0.0 for v in sv]
    hi = [100.0 * v if v != 0 else 1.0 for v in sv]

    try:
        f, p, e, c, n, s = _langau_fit_once(
            hist, first_stage_range, sv, lo, hi, "langau_fit"
        )
    except Exception:
        return None

    # Error checking on MPV (p[1]) and GSigma (p[3])
    if e[3] > max_par_err or e[1] > max_par_err:
        return None

    return dict(func=f, pars=p, errs=e, chi2=c, ndf=n, status=s)


def langau_fit_two_stage(
    hist, first_stage_range=(0.0, 20.0), max_par_err=1.0
):
    """Two-stage Landau-Gaussian fit, ported 1:1 from the ROOT macro's
    inline double-fit logic (langaufit called twice with a narrowed window).

    Returns a dict with the Stage-2 function/parameters, or None if the
    Stage-2 MPV or GSigma error exceeds max_par_err.
    """
    # ROOT uses the raw max-bin center as the MPV start guess -- no clamping
    max_x = hist.GetBinCenter(hist.GetMaximumBin())

    bin_width = hist.GetBinWidth(1)
    integral = hist.Integral()
    fit_opts = "RBOSQN"

    # --- Stage 1 Fit Setup ---
    # sv[2] restored to match ROOT's "Integral * 0.05 * bin_width" exactly
    sv = [0.1, max_x, integral * 0.05 * bin_width, 0.2]
    # ROOT uses one uniform rule for all four parameters: [0.01x, 100x]
    lo = [0.01 * v for v in sv]
    hi = [100.0 * v for v in sv]

    try:
        f1, p1, e1, c1, n1, s1 = _langau_fit_once(
            hist,
            first_stage_range,
            sv,
            lo,
            hi,
            f"langau_stage1_{hist.GetName()}",
            fit_options=fit_opts,
        )
    except Exception:
        return None

    # ROOT has no acceptance gate here -- it always proceeds to Stage 2
    # regardless of Stage-1 MPV sign or parameter error.

    # --- Stage 2 Window ---
    # ROOT always uses a fixed 0.8x-1.5x MPV1 window, with no clamping
    # back into first_stage_range.
    mpv_stage1 = p1[1]
    second_stage_range = (mpv_stage1 * 0.6, mpv_stage1 * 1.4)

    # ROOT reuses the *original* Stage-1 bound arrays (lo, hi) unchanged for
    # Stage 2 -- only the start values are updated, to the Stage-1 result.
    try:
        f2, p2, e2, c2, n2, s2 = _langau_fit_once(
            hist,
            second_stage_range,
            p1,
            lo,
            hi,
            f"langau_stage2_{hist.GetName()}",
            fit_options=fit_opts,
        )
    except Exception:
        return None

    # ROOT's only post-fit cut: reject if the MPV or GSigma error > 1
    # (here generalized to max_par_err). No status/MPV-sign check, and
    # no fallback to the Stage-1 result -- a rejected Stage-2 fit is just
    # dropped, same as the macro's "continue".
    if e2[1] > max_par_err or e2[3] > max_par_err:
        return None

    return dict(
        func=f2,
        pars=p2,
        errs=e2,
        chi2=c2,
        ndf=n2,
        status=s2,
        fit_stage=2,
        stage1_func=f1,
        stage1_pars=p1,
        stage2_range=second_stage_range,
    )


def make_theoretical_mpv_func(hfit, pdg, mass=None):
    """Wrapper around `build_theoretical_pdf` that computes the scalar MPV."""

    def theoretical_mpv(rr_center, mean_pitch):
        pdf = build_theoretical_pdf(
            hfit, pdg, rr_center, mean_pitch, mass=mass
        )
        return pdf.GetMaximumX()

    return theoretical_mpv
# ===========================================================================
# 6. Driver + Data vs. MC Comparison Plotting
# ===========================================================================
PLANE_COLORS = {0: "tab:blue", 1: "tab:orange", 2: "tab:green"}



def power_law(x, a, b, c):
    return a + b * np.power(x, c)


def exp_decay_plateau(rr, a, b, c):
    """a + b * exp(-rr / c)"""
    return a + b * np.exp(-rr / c)


def exp_decay_linear(rr, a, b, c, d):
    """a + b * exp(-rr / c) + d * rr"""
    return a + b * np.exp(-rr / c) + d * rr


def power_law_linear(rr, a, b, p, d):
    """a + b * rr^p + d * rr"""
    return a + b * np.power(rr, p) + d * rr


# Shift models: name -> (function, parameter names)
SHIFT_MODELS = {
    "exp_plateau": (exp_decay_plateau, ("a", "b", "c")),
    "exp_linear": (exp_decay_linear, ("a", "b", "c", "d")),
    "power_linear": (power_law_linear, ("a", "b", "p", "d")),
}
SHIFT_MODEL = "power_linear"   # shift model of every shift fit, a key of SHIFT_MODELS


def eval_shift_fit(x, popt, model, rr_floor=None):
    """Shift model `model` at residual range x. With rr_floor, the model is evaluated
    at max(x, rr_floor), so it stays constant below rr_floor."""
    x = np.asarray(x, dtype=float)
    if rr_floor is not None:
        x = np.maximum(x, rr_floor)
    return SHIFT_MODELS[model][0](x, *popt)


def shift_fit_params(popt, perr, model):
    """{"model": model, name: value, name_err: error} for the shift-fit parameters."""
    out = {"model": model}
    for name, val, err in zip(SHIFT_MODELS[model][1], popt, perr):
        out[name], out[name + "_err"] = float(val), float(err)
    return out


def fit_sigma_vs_mpv(res, p0=(0.05, 0.1, 2.0),
                     bounds=((0.0, 0.0, -20.0), (2.0, 100.0, 20.0))):
    """sigma_G = a + b * MPV^c, MPV = theoretical MPV of the slice. Uses only the
    slices flagged sigma_ok by run_slice_fits. Returns (popt, perr)."""
    v = res.dropna(subset=["mpv_theory", "gsigma", "gsigma_err"])
    v = v[v["gsigma_err"] > 0]
    if "sigma_ok" in v:
        v = v[v["sigma_ok"]]
    if len(v) < 4:
        return None, None
    try:
        popt, pcov = curve_fit(power_law, v["mpv_theory"].to_numpy(), v["gsigma"].to_numpy(),
                               p0=p0, sigma=v["gsigma_err"].to_numpy(), absolute_sigma=True,
                               bounds=bounds, maxfev=20000)
    except Exception as exc:
        print(f"  sigma fit failed: {exc}")
        return None, None
    return popt, np.sqrt(np.diag(pcov))


def fit_shift_vs_rr(points, model="exp_linear", c_max=50.0, p_bounds=(-4.0, -0.05),
                    verbose=False):
    """Shift vs residual range. model: "exp_linear" (a, b, c, d), "exp_plateau" (a, b, c)
    or "power_linear" (a, b, p, d). c_max bounds the decay length, so the exponential
    stays the fast component; p_bounds bounds the power-law exponent."""
    if model not in SHIFT_MODELS:
        raise ValueError(f"Unknown shift model: {model}")
    func, names = SHIFT_MODELS[model]
    v = points.dropna(subset=["rr_center", "shift", "shift_err"])
    v = v[v["shift_err"] > 0].sort_values("rr_center")
    if len(v) < len(names) + 2:
        return None, None
    x, y, yerr = (v[c].to_numpy() for c in ("rr_center", "shift", "shift_err"))

    c0, p_start = 5.0, -1.0
    if model in ("exp_linear", "power_linear"):
        tail = x >= x.min() + (x.max() - x.min()) / 3.0
        d0, a0 = np.polyfit(x[tail], y[tail], 1, w=1.0 / yerr[tail])
        excess0 = y[0] - (a0 + d0 * x[0])
    if model == "exp_linear":
        p0 = [a0, excess0 * np.exp(x[0] / c0), c0, d0]
        bounds = ([-np.inf, -np.inf, 0.1, -np.inf], [np.inf, np.inf, c_max, np.inf])
    elif model == "power_linear":
        p0 = [a0, excess0 * np.power(x[0], -p_start), p_start, d0]
        bounds = ([-np.inf, -np.inf, p_bounds[0], -np.inf],
                  [np.inf, np.inf, p_bounds[1], np.inf])
    else:
        p0 = [y[-1], y[0] - y[-1], c0]
        bounds = ([-np.inf, -np.inf, 0.1], [np.inf, np.inf, c_max])

    try:
        popt, pcov = curve_fit(func, x, y, p0=p0, sigma=yerr, absolute_sigma=True,
                               bounds=bounds, maxfev=20000)
    except Exception as exc:
        print(f"  shift fit failed: {exc}")
        return None, None
    if verbose:
        chi2 = np.sum(((y - func(x, *popt)) / yerr) ** 2)
        print(f"  shift fit ({model}): chi2/ndf = {chi2:.1f}/{len(x) - len(popt)}")
    return popt, np.sqrt(np.diag(pcov))



''''
# ===========================================================================
# 2. Histogram construction from a pandas Series (ported from
#    `th1_from_series`, simplified to fixed binning / no weights, since the
#    hit dataframes don't carry a weight column)
# ===========================================================================


# ===========================================================================
# 3. Two-stage Langau fit, ported from `langaufit` + the broad-then-refined
#    double-fit logic that lived inline in FitPlotWholeDataset()
# ===========================================================================
import ROOT
from ROOT import TF1








def fit_rr_slices(
    df,
    plane,
    tpc,
    dedx_col="dedx",
    rr_min=3.0,
    rr_max=40.0,
    rr_bin_width=1.0,
    hist_nbins=150,
    hist_xmin=0,
    hist_xmax=20.0,
    min_entries=30,
    first_stage_range=(0.5, 20.0),
    theoretical_mpv_func=None,
    max_gsigma_err=0.2,  # Added max error threshold
    min_gsigma_err=0.0001,  # Added max error threshold
    verbose=False,
):
    """Fits individual RR slices and collects MPV and Gaussian Smearing parameters."""
    edges = make_rr_edges(rr_min, rr_max, rr_bin_width)
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        sl = df[(df["rr"] >= lo) & (df["rr"] < hi)]
        if len(sl) < min_entries:
            continue

        rr_center = 0.5 * (lo + hi)
        hname = f"slice_p{plane}_t{tpc}_rr{rr_center:.2f}"
        hist = th1_from_series(
            sl[dedx_col], hname, "", hist_nbins, hist_xmin, hist_xmax
        )

        fit = langau_fit_two_stage(hist, first_stage_range=first_stage_range)
        #fit = langau_fit_single_stage(hist, first_stage_range=first_stage_range)
        if fit is None:
            if verbose:
                print(
                    f"   [plane {plane}, tpc {tpc}] rr={rr_center:.2f}: fit rejected"
                )
            continue

        width, mpv_reco, area, gsigma = fit["pars"]
        w_e, mpv_e, area_e, gsigma_e = fit["errs"]

        # Filter out slices where gsigma error exceeds threshold right at collection time
        if gsigma_e > max_gsigma_err or not np.isfinite(gsigma_e) or gsigma_e < min_gsigma_err:
            if verbose:
                print(
                    f"   [plane {plane}, tpc {tpc}] rr={rr_center:.2f}: dropped due to gsigma_err={gsigma_e:.3f} > {max_gsigma_err}"
                )
            continue

        if theoretical_mpv_func is not None:
            mean_pitch = sl["pitch"].mean()
            mpv_x = theoretical_mpv_func(rr_center, mean_pitch)
            mpv_x_err = 0.0
        else:
            mpv_x = mpv_reco
            mpv_x_err = mpv_e

        rows.append(
            dict(
                plane=plane,
                tpc=tpc,
                rr_center=rr_center,
                n_hits=len(sl),
                mpv_x=mpv_x,
                mpv_x_err=mpv_x_err,
                mpv_reco=mpv_reco,
                mpv_reco_err=mpv_e,
                width=width,
                width_err=w_e,
                area=area,
                area_err=area_e,
                gsigma=gsigma,
                gsigma_err=gsigma_e,
                chi2=fit["chi2"],
                ndf=fit["ndf"],
                status=fit["status"],
            )
        )

    return pd.DataFrame(rows)



def fit_sigma_vs_mpv(
    res_df,
    x_col="mpv_x",
    y_col="gsigma",
    yerr_col="gsigma_err",
    max_yerr=0.1,  # Added threshold parameter
    p0=(0.2, 0.1, 1.5),
    bounds=((0.0, 0.0, -10.0), (1e6, 1.0, 10.0)),
):
    """Fits power-law model to sigma_G vs MPV, excluding points with yerr > max_yerr."""
    # 1. Filter the DataFrame to remove points with yerr > max_yerr or invalid values
    mask = (res_df[yerr_col] <= max_yerr) & (res_df[yerr_col] > 0) & np.isfinite(res_df[yerr_col])
    filtered_df = res_df[mask]

    # Ensure we still have enough points to perform a 3-parameter fit
    if len(filtered_df) < 3:
        raise RuntimeError(
            f"Not enough valid points left for fitting after filtering (yerr <= {max_yerr}). "
            f"Only {len(filtered_df)} points remaining."
        )

    # 2. Extract arrays from filtered data
    x = filtered_df[x_col].to_numpy()
    y = filtered_df[y_col].to_numpy()
    yerr = filtered_df[yerr_col].to_numpy()

    # 3. Fallback handling for non-finite values if any remain
    fallback = np.nanmedian(yerr)
    yerr = np.where(
        np.isfinite(yerr), yerr, fallback if np.isfinite(fallback) else 1.0
    )

    # 4. Perform the curve fit
    popt, pcov = curve_fit(
        power_law,
        x,
        y,
        p0=p0,
        sigma=yerr,
        absolute_sigma=True,
        bounds=bounds,
        maxfev=20000,
    )
    perr = np.sqrt(np.diag(pcov))
    return popt, perr



def analyze(
    hit_dfs,
    particle="muon",
    dedx_col="dedx",
    rr_max_by_particle=None,
    theoretical_mpv_func=None,
    out_prefix="langau_rr",
    make_summary_plots=True,
    verbose=True,
):
    """Analyzes a set of hit DataFrames across planes and TPCs.

    In addition to the per-TPC fits (tpc = 0, 1), each plane is also fit
    with both TPCs combined, stored under tpc = -1.
    """
    if rr_max_by_particle is None:
        rr_max_by_particle = {"muon": 80.0, "pion": 40.0, "proton": 60.0}
    rr_max = rr_max_by_particle.get(particle, 80.0)

    all_results = {}
    fit_params = {}

    for plane, df in enumerate(hit_dfs):
        # tpc == -1 means "both TPCs combined": fit the full per-plane
        # dataframe instead of filtering down to a single TPC.
        for tpc in (0, 1, -1):
            if verbose:
                tpc_label = "combined" if tpc == -1 else tpc
                print(f"Plane {plane}, TPC {tpc_label}")

            sub = df if tpc == -1 else df[df["tpc"] == tpc]

            res = fit_rr_slices(
                sub,
                plane,
                tpc,
                dedx_col=dedx_col,
                rr_max=rr_max,
                theoretical_mpv_func=theoretical_mpv_func,
                verbose=verbose,
            )
            all_results[(plane, tpc)] = res

            popt = perr = None
            if len(res) >= 3:
                try:
                    popt, perr = fit_sigma_vs_mpv(res)
                except RuntimeError as exc:
                    if verbose:
                        print(f"  power-law fit failed: {exc}")
            fit_params[(plane, tpc)] = (popt, perr)

    non_empty = [
        r.assign(plane=p, tpc=t)
        for (p, t), r in all_results.items()
        if len(r)
    ]
    combined = (
        pd.concat(non_empty, ignore_index=True)
        if non_empty
        else pd.DataFrame()
    )

    if len(combined):
        combined.to_hdf(
            f"{out_prefix}_slices.h5",
            key="fits",
            mode="w",
            format="table",
            complib="blosc",
            complevel=9,
        )

    return all_results, fit_params
# ===========================================================================
# 8. Driver for the theoretical-PDF (x) smearing-only fit -- parallel to
#    `analyze`, but calls `fit_rr_slices_theoretical` instead of
#    `fit_rr_slices`, and uses `mpv_theory` (not `mpv_x`) as the x-axis
#    for the sigma-vs-MPV power-law fit, since here the theoretical MPV
#    *is* the x-axis by construction (pitch fixed at 0.32, no reco MPV).
# ===========================================================================
def analyze_theoretical(
    hit_dfs,
    hfit,
    pdg,
    particle="muon",
    dedx_col="dedx",
    rr_max_by_particle=None,
    pitch=0.32,
    mass=None,
    sigma0=0.2,
    sigma_bounds=(1e-4, 5.0),
    max_sigma_err=0.2,
    out_prefix="langau_rr_theoretical",
    make_summary_plots=True,
    verbose=True,
):
    """Analyzes a set of hit DataFrames across planes and TPCs, fitting the
    fixed theoretical PDF (x) zero-mean-Gaussian(sigma) to each rr slice
    instead of a free 4-parameter Langau.

    In addition to the per-TPC fits (tpc = 0, 1), each plane is also fit
    with both TPCs combined, stored under tpc = -1 -- same convention as
    `analyze`.
    """
    if rr_max_by_particle is None:
        rr_max_by_particle = {"muon": 80.0, "pion": 40.0, "proton": 60.0}
    rr_max = rr_max_by_particle.get(particle, 80.0)

    all_results = {}
    fit_params = {}

    for plane, df in enumerate(hit_dfs):
        # tpc == -1 means "both TPCs combined": fit the full per-plane
        # dataframe instead of filtering down to a single TPC.
        for tpc in (0, 1, -1):
            if verbose:
                tpc_label = "combined" if tpc == -1 else tpc
                print(f"[theoretical] Plane {plane}, TPC {tpc_label}")

            sub = df if tpc == -1 else df[df["tpc"] == tpc]

            res = fit_rr_slices_theoretical(
                sub,
                plane,
                tpc,
                hfit=hfit,
                pdg=pdg,
                dedx_col=dedx_col,
                rr_max=rr_max,
                pitch=pitch,
                mass=mass,
                sigma0=sigma0,
                sigma_bounds=sigma_bounds,
                max_sigma_err=max_sigma_err,
                verbose=verbose,
            )
            all_results[(plane, tpc)] = res

            popt = perr = None
            if len(res) >= 3:
                try:
                    popt, perr = fit_sigma_vs_mpv(
                        res, x_col="mpv_theory", y_col="gsigma", yerr_col="gsigma_err"
                    )
                except RuntimeError as exc:
                    if verbose:
                        print(f"  power-law fit failed: {exc}")
            fit_params[(plane, tpc)] = (popt, perr)

    non_empty = [
        r.assign(plane=p, tpc=t)
        for (p, t), r in all_results.items()
        if len(r)
    ]
    combined = (
        pd.concat(non_empty, ignore_index=True)
        if non_empty
        else pd.DataFrame()
    )

    if len(combined):
        combined.to_hdf(
            f"{out_prefix}_slices.h5",
            key="fits",
            mode="w",
            format="table",
            complib="blosc",
            complevel=9,
        )

    return all_results, fit_params







def convolved_theoretical_pdf(x_query, sigma, amplitude, x_grid, pdf_vals, dx):
    """model(x) = amplitude * (theoretical PDF (x) Gaussian(0, sigma))(x).

    x_grid/pdf_vals/dx are the *fixed* physics PDF (pitch=0.32, rr=slice
    center) -- only sigma and amplitude change during the fit.
    """
    x_centered = x_grid - x_grid[len(x_grid) // 2]
    kernel = gaussian_kernel(x_centered, sigma)
    kernel = kernel / (kernel.sum() * dx)  # == unit-area Gaussian
    conv = fftconvolve(pdf_vals, kernel, mode="same") * dx
    return amplitude * np.interp(x_query, x_grid, conv)


def fit_theoretical_conv_slice(
    bin_centers,
    bin_counts,
    bin_errs,
    x_grid,
    pdf_vals,
    dx,
    mpv_theory,
    fit_range=None,
    sigma0=0.2,
    sigma_bounds=(1e-4, 5.0),
):
    """curve_fit wrapper: only (sigma, amplitude) float. Fit range defaults
    to [0.7, 1.5] x mpv_theory, matching the Langau Stage-2 window."""
    if fit_range is None:
        fit_range = (0.7 * mpv_theory, 1.5 * mpv_theory)

    mask = (bin_centers >= fit_range[0]) & (bin_centers <= fit_range[1]) & (bin_counts > 0)
    if mask.sum() < 5:
        return None

    xs, ys, yerr = bin_centers[mask], bin_counts[mask], bin_errs[mask]
    yerr = np.where(yerr > 0, yerr, 1.0)

    # == Rough amplitude start: match data peak height to the (unsmeared)
    # == PDF peak height -- curve_fit refines both params from here.
    amp0 = ys.max() / max(pdf_vals.max(), 1e-9)

    def model(x, sigma, amplitude):
        return convolved_theoretical_pdf(x, sigma, amplitude, x_grid, pdf_vals, dx)

    try:
        popt, pcov = curve_fit(
            model, xs, ys, p0=[sigma0, amp0],
            sigma=yerr, absolute_sigma=True,
            bounds=([sigma_bounds[0], 0.0], [sigma_bounds[1], np.inf]),
            maxfev=20000,
        )
    except Exception:
        return None

    perr = np.sqrt(np.diag(pcov))
    sigma_fit, amp_fit = popt
    sigma_err, amp_err = perr
    return dict(
        sigma=sigma_fit, sigma_err=sigma_err,
        amplitude=amp_fit, amplitude_err=amp_err,
        fit_range=fit_range,
    )


def fit_rr_slices_theoretical(
    df,
    plane,
    tpc,
    hfit,
    pdg,
    dedx_col="dedx",
    rr_min=3.0,
    rr_max=40.0,
    rr_bin_width=1.0,
    hist_nbins=150,
    hist_xmin=0.0,
    hist_xmax=20.0,
    min_entries=30,
    pitch=0.32,  # == fixed, per the hardcoded C++ pitch / distribution mode
    mass=None,
    sigma0=0.2,
    sigma_bounds=(1e-4, 5.0),
    max_sigma_err=0.2,
    verbose=False,
):
    """Same rr-slicing loop as fit_rr_slices, but fits the fixed theoretical
    PDF (x) zero-mean-Gaussian(sigma) instead of a free 4-parameter Langau.

    Output columns are aligned with fit_rr_slices's output (mpv_x/mpv_x_err
    is the theoretical MPV here, mpv_x_err fixed at 0 since it isn't a
    fitted quantity) so downstream code (fit_sigma_vs_mpv, plot_all_planes,
    plot_by_plane) works unchanged on either analyze(...) or
    analyze_theoretical(...) output.
    """
    edges = make_rr_edges(rr_min, rr_max, rr_bin_width)
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        sl = df[(df["rr"] >= lo) & (df["rr"] < hi)]
        if len(sl) < min_entries:
            continue
        rr_center = 0.5 * (lo + hi)

        pdf = build_theoretical_pdf(hfit, pdg, rr_center, pitch, mass=mass)
        mpv_theory = robust_max_x_py(pdf, 0.0, 10.0, 2000)
        x_grid, pdf_vals, dx = build_pdf_grid(pdf)

        counts, bin_edges = np.histogram(
            sl[dedx_col].to_numpy(), bins=hist_nbins, range=(hist_xmin, hist_xmax)
        )
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        bin_errs = np.where(counts > 0, np.sqrt(counts), 1.0)

        fit = fit_theoretical_conv_slice(
            bin_centers, counts.astype(float), bin_errs,
            x_grid, pdf_vals, dx, mpv_theory,
            sigma0=sigma0, sigma_bounds=sigma_bounds,
        )
        if fit is None:
            if verbose:
                print(f"   [plane {plane}, tpc {tpc}] rr={rr_center:.2f}: fit failed")
            continue
        if fit["sigma_err"] > max_sigma_err or not np.isfinite(fit["sigma_err"]):
            if verbose:
                print(f"   [plane {plane}, tpc {tpc}] rr={rr_center:.2f}: "
                      f"dropped, sigma_err={fit['sigma_err']:.3f}")
            continue

        rows.append(dict(
            plane=plane,
            tpc=tpc,
            rr_center=rr_center,
            n_hits=len(sl),
            # == mpv_x/mpv_x_err: same names fit_rr_slices uses, so
            # == downstream code (fit_sigma_vs_mpv, plotting) is identical
            # == regardless of which fitter produced the dataframe.
            mpv_x=mpv_theory,
            mpv_x_err=0.0,
            # == kept too, for anything that wants the theoretical value by
            # == its more explicit name
            mpv_theory=mpv_theory,
            width=np.nan,       # == no free Landau width in this fit
            width_err=np.nan,
            area=fit["amplitude"],
            area_err=fit["amplitude_err"],
            gsigma=fit["sigma"],
            gsigma_err=fit["sigma_err"],
            chi2=np.nan,        # == not computed by curve_fit here
            ndf=np.nan,
            status=np.nan,
            fit_range_lo=fit["fit_range"][0],
            fit_range_hi=fit["fit_range"][1],
        ))

    return pd.DataFrame(rows)
'''