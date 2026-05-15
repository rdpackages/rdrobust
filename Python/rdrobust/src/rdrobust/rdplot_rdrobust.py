"""
Diagnostic plot for rdrobust objects.

Mirrors the R function plot.rdrobust(): binned scatter within the bandwidth
window, local-polynomial fit lines, pointwise CI bands, and an optional
"effect" forest-plot panel.
"""

import numpy as np
import pandas as pd
from plotnine import (
    ggplot, aes, geom_point, geom_line, geom_ribbon, geom_vline,
    geom_errorbarh, annotate, labs, theme_bw, theme, element_blank,
    element_text, element_line, scale_size_identity, coord_cartesian,
    scale_y_discrete
)
from plotnine.options import figure_size
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import norm


def plot_rdrobust(obj, y, x,
                  nbins=20,
                  ci=True,
                  show_effect=False,
                  title=None,
                  x_label="Running Variable",
                  y_label="Outcome",
                  x_lim=None,
                  y_lim=None,
                  col_l="#3B7DD8",
                  col_r="#D95F3B",
                  base_size=14,
                  figsize=(8, 5)):
    """
    Plot an rdrobust result object.

    Parameters
    ----------
    obj : rdrobust_output
        Result from rdrobust().
    y : array-like
        Outcome variable (same data used in rdrobust call).
    x : array-like
        Running variable (same data used in rdrobust call).
    nbins : int
        Number of bins per side (default 20).
    ci : bool
        Show pointwise CI bands for the polynomial fit (default True).
    show_effect : bool
        Show a second panel with the RD effect + CI (default False).
    title : str or None
        Plot title.
    x_label, y_label : str
        Axis labels.
    x_lim, y_lim : tuple or None
        Axis limits as (min, max).
    col_l, col_r : str
        Colors for left/right of cutoff.
    base_size : int
        Base font size.
    figsize : tuple
        Figure size in inches (width, height).

    Returns
    -------
    matplotlib.figure.Figure
    """
    x = np.asarray(x).ravel()
    y = np.asarray(y).ravel()

    # --- extract from rdrobust object ---
    c_val = obj.c
    h_l = obj.bws.iloc[0, 0]
    h_r = obj.bws.iloc[0, 1]
    coefs_l = np.asarray(obj.beta_p_l).ravel()
    coefs_r = np.asarray(obj.beta_p_r).ravel()
    V_l = np.asarray(obj.V_cl_l)
    V_r = np.asarray(obj.V_cl_r)
    kernel = obj.kernel
    tau_cl = obj.coef.iloc[0, 0]
    ci_rb = [obj.ci.iloc[2, 0], obj.ci.iloc[2, 1]]
    lev = obj.level
    z_crit = norm.ppf(1 - (1 - lev / 100) / 2)
    is_fuzzy = obj.rdmodel is not None and "fuzzy" in obj.rdmodel.lower()

    # --- kernel weight function ---
    def kern_w(u):
        u = np.abs(u)
        k = kernel.lower()[:3] if isinstance(kernel, str) else "tri"
        if k == "tri":
            return np.maximum(1 - u, 0)
        elif k == "epa":
            return np.maximum(0.75 * (1 - u**2), 0)
        elif k == "uni":
            return (u <= 1).astype(float) * 0.5
        return np.maximum(1 - u, 0)

    # --- filter to bandwidth window ---
    in_win = (x >= c_val - h_l) & (x <= c_val + h_r)
    x_w = x[in_win]
    y_w = y[in_win]
    side = np.where(x_w < c_val, "left", "right")

    kw = np.where(
        side == "left",
        kern_w((x_w - c_val) / h_l),
        kern_w((x_w - c_val) / h_r)
    )

    # --- binned means + average kernel weight ---
    def make_bins(xv, yv, wv, nb):
        brks = np.linspace(np.min(xv), np.max(xv), nb + 1)
        grp = np.digitize(xv, brks, right=False)
        grp = np.clip(grp, 1, nb)
        xmid = np.array([np.mean(xv[grp == i]) for i in range(1, nb + 1) if np.any(grp == i)])
        ymid = np.array([np.mean(yv[grp == i]) for i in range(1, nb + 1) if np.any(grp == i)])
        wavg = np.array([np.mean(wv[grp == i]) for i in range(1, nb + 1) if np.any(grp == i)])
        return xmid, ymid, wavg

    idx_l = side == "left"
    idx_r = side == "right"

    if np.sum(idx_l) > 0:
        xmid_l, ymid_l, wavg_l = make_bins(x_w[idx_l], y_w[idx_l], kw[idx_l], nbins)
    else:
        xmid_l = ymid_l = wavg_l = np.array([])

    if np.sum(idx_r) > 0:
        xmid_r, ymid_r, wavg_r = make_bins(x_w[idx_r], y_w[idx_r], kw[idx_r], nbins)
    else:
        xmid_r = ymid_r = wavg_r = np.array([])

    # rescale kernel weights to point size range [1.5, 5]
    all_w = np.concatenate([wavg_l, wavg_r])
    w_min = np.min(all_w) if len(all_w) > 0 else 0
    w_rng = np.max(all_w) - w_min if len(all_w) > 0 else 1
    if w_rng == 0:
        w_rng = 1

    def rescale_sz(w):
        return 1.5 + 3.5 * (w - w_min) / w_rng

    sz_l = rescale_sz(wavg_l)
    sz_r = rescale_sz(wavg_r)

    # --- polynomial fit + pointwise CI bands ---
    def poly_fit(xv, coefs, c0):
        xc = xv - c0
        mat = np.column_stack([xc**j for j in range(len(coefs))])
        return mat @ coefs

    def poly_se(xv, V, c0):
        xc = xv - c0
        mat = np.column_stack([xc**j for j in range(V.shape[0])])
        vv = np.sum((mat @ V) * mat, axis=1)
        return np.sqrt(np.maximum(vv, 0))

    xseq_l = np.linspace(c_val - h_l, c_val, 200)
    xseq_r = np.linspace(c_val, c_val + h_r, 200)
    yhat_l = poly_fit(xseq_l, coefs_l, c_val)
    yhat_r = poly_fit(xseq_r, coefs_r, c_val)

    y0_l = coefs_l[0]
    y0_r = coefs_r[0]

    # --- annotation label ---
    fuzzy_str = "\n(Fuzzy)" if is_fuzzy else ""
    lbl = (f"RD = {tau_cl:.3f}\n"
           f"{lev:.0f}% RBC CI: [{ci_rb[0]:.3f}, {ci_rb[1]:.3f}]"
           f"{fuzzy_str}")

    # --- build figure with matplotlib ---
    if show_effect:
        fig = plt.figure(figsize=(figsize[0], figsize[1] * 1.3))
        gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1], hspace=0.25)
        ax_rd = fig.add_subplot(gs[0])
        ax_eff = fig.add_subplot(gs[1])
    else:
        fig, ax_rd = plt.subplots(1, 1, figsize=figsize)
        ax_eff = None

    # --- RD scatter + fit plot ---
    # CI bands
    if ci:
        se_l = poly_se(xseq_l, V_l, c_val)
        se_r = poly_se(xseq_r, V_r, c_val)
        ax_rd.fill_between(xseq_l, yhat_l - z_crit * se_l,
                           yhat_l + z_crit * se_l,
                           color=col_l, alpha=0.12)
        ax_rd.fill_between(xseq_r, yhat_r - z_crit * se_r,
                           yhat_r + z_crit * se_r,
                           color=col_r, alpha=0.12)

    # Binned scatter (size = kernel weight)
    if len(xmid_l) > 0:
        ax_rd.scatter(xmid_l, ymid_l, s=sz_l * 15, c=col_l, alpha=0.85, zorder=3)
    if len(xmid_r) > 0:
        ax_rd.scatter(xmid_r, ymid_r, s=sz_r * 15, c=col_r, alpha=0.85, zorder=3)

    # Polynomial fits
    ax_rd.plot(xseq_l, yhat_l, color=col_l, linewidth=1.5)
    ax_rd.plot(xseq_r, yhat_r, color=col_r, linewidth=1.5)

    # Cutoff line
    ax_rd.axvline(c_val, linestyle="--", color="grey", linewidth=0.5)

    # Open circles at cutoff
    ax_rd.plot(c_val, y0_l, 'o', markersize=7, markerfacecolor='white',
               markeredgecolor=col_l, markeredgewidth=1.5, zorder=5)
    ax_rd.plot(c_val, y0_r, 'o', markersize=7, markerfacecolor='white',
               markeredgecolor=col_r, markeredgewidth=1.5, zorder=5)

    # Annotation
    ax_rd.annotate(lbl,
                   xy=(c_val + 0.04 * (h_l + h_r), ax_rd.get_ylim()[1] if ax_rd.get_ylim()[1] else y0_r),
                   fontsize=base_size * 0.7,
                   color="grey",
                   verticalalignment='top',
                   linespacing=1.3)

    # Style
    ax_rd.set_xlabel(x_label, fontsize=base_size)
    ax_rd.set_ylabel(y_label, fontsize=base_size)
    if title:
        ax_rd.set_title(title, fontsize=base_size + 2, fontweight='bold')
    ax_rd.tick_params(labelsize=base_size - 1)
    ax_rd.spines['top'].set_visible(False)
    ax_rd.spines['right'].set_visible(False)
    ax_rd.grid(True, axis='both', color='#EBEBEB', linewidth=0.4)

    if x_lim is not None:
        ax_rd.set_xlim(x_lim)
    if y_lim is not None:
        ax_rd.set_ylim(y_lim)

    # Re-annotate with correct y-position after axis limits are set
    # (clear previous annotation, re-add)
    for child in ax_rd.texts[:]:
        child.remove()
    ylims = ax_rd.get_ylim()
    ax_rd.annotate(lbl,
                   xy=(c_val + 0.04 * (h_l + h_r), ylims[1]),
                   fontsize=base_size * 0.7,
                   color="grey",
                   verticalalignment='top',
                   linespacing=1.3)

    # --- Effect panel ---
    if show_effect and ax_eff is not None:
        pv_rb = obj.pv.iloc[2, 0]
        if pv_rb < 0.01:
            stars = "***"
        elif pv_rb < 0.05:
            stars = "**"
        elif pv_rb < 0.10:
            stars = "*"
        else:
            stars = ""

        eff_lbl = f"{tau_cl:.3f}{stars}  [{ci_rb[0]:.3f}, {ci_rb[1]:.3f}]"
        eff_size = base_size - 4

        ax_eff.axvline(0, linestyle="--", color="grey", linewidth=0.4)
        ax_eff.errorbar(tau_cl, 0, xerr=[[tau_cl - ci_rb[0]], [ci_rb[1] - tau_cl]],
                        fmt='o', color='grey', markersize=5, linewidth=0.8,
                        capsize=0)
        ax_eff.annotate(eff_lbl, xy=(tau_cl, 0.25),
                        fontsize=eff_size * 0.8, color="grey",
                        ha='center')
        ax_eff.set_yticks([0])
        ax_eff.set_yticklabels(["RD Effect"], fontsize=eff_size)
        ax_eff.set_xlabel(f"Estimate ({lev:.0f}% RBC CI)", fontsize=eff_size)
        ax_eff.set_ylim(-0.6, 0.6)
        ax_eff.spines['top'].set_visible(False)
        ax_eff.spines['right'].set_visible(False)
        ax_eff.spines['left'].set_visible(False)
        ax_eff.tick_params(left=False, labelsize=eff_size - 1)
        ax_eff.grid(True, axis='x', color='#EBEBEB', linewidth=0.3)

        caption = "Point estimate: conventional. CI: robust bias-corrected. * p<0.10  ** p<0.05  *** p<0.01"
        fig.text(0.12, 0.01, caption, fontsize=eff_size - 2, color="grey")
        fig.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.08)
    else:
        fig.tight_layout()
    return fig
