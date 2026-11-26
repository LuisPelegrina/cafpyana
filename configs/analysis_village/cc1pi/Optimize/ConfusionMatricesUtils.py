
from analysis_village.cc1pi.GraphUtils.GraphsUtils import *
import numpy as np
def apply_nd_cuts(df, cols, cuts, cut_sign):
    """
    Apply N-dimensional cuts with per-cut direction.

    cut_sign: list of "<" or ">"
    """
    mask = np.ones(len(df), dtype=bool)

    for col, cut, sign in zip(cols, cuts, cut_sign):
        if sign == ">":
            mask &= df[col] > cut
        elif sign == "<":
            mask &= df[col] < cut
        else:
            raise ValueError(f"Invalid cut sign: {sign}")

    return mask

def build_confusion_matrix_from_cuts(
    signal_df,
    bkg_df,
    cols,
    cuts,
    cut_sign,
):
    sig_pass = apply_nd_cuts(signal_df,cols,cuts,cut_sign)
    bkg_pass = apply_nd_cuts(bkg_df,cols,cuts,cut_sign)

    cm = np.array([
        [np.sum(sig_pass) ,  np.sum(bkg_pass)],  
        [np.sum(~sig_pass),  np.sum(~bkg_pass)],
    ])
    '''
    cm = np.array([
        [1 ,  2],  
        [3,  4],
    ])
    '''

    return cm

def plot_confusion_matrix(
    cm,
    x_labels,   # True labels
    y_labels,   # Reco labels
    cmap=sunset_cmap,
):
    cm = np.asarray(cm, dtype=float)

 # -------------------------
    # Normalizations
    # -------------------------
    # Purity (row-normalized): Sum over True for a given Reco
    row_denom = cm.sum(axis=1, keepdims=True)
    cm_purity = np.divide(cm, row_denom, where=row_denom > 0) * 100

    # Efficiency (column-normalized): Sum over Reco for a given True
    col_denom = cm.sum(axis=0, keepdims=True)
    cm_eff = np.divide(cm, col_denom, where=col_denom > 0) * 100

    # Eff * Pur
    cm_effpur = cm_purity * cm_eff / 100.0

    # -------------------------
    # Plot
    # -------------------------
    fig, axes = plt.subplots(
        1, 2, figsize=(13, 6), sharey=True, constrained_layout=True
    )


    matrices = [
        (cm_purity, "Purity", True),
        (cm_eff, "Efficiency", True),
        #(cm_effpur, "Eff × Pur", False),
    ]

    for ax, (cm_plot, title, show_counts) in zip(axes, matrices):
        im = ax.imshow(cm_plot, cmap=cmap, vmin=0, vmax=100, origin="lower", aspect='auto')

        ax.set_xticks(np.arange(len(x_labels)))
        ax.set_yticks(np.arange(len(y_labels)))
        
        ax.set_xticklabels(x_labels, rotation=45)
        ax.set_yticklabels(y_labels)

        ax.set_xlabel("True")
        ax.set_title(title)

        for i in range(cm.shape[0]):    
            for j in range(cm.shape[1]): 
                value = cm_plot[i, j]
                txt = f"{int(cm[i,j])}\n{value:.1f}%" if show_counts else f"{value:.1f}%"
                text_color = "white" if value > 60 else "black"
                ax.text(j, i, txt, ha="center", va="center", color=text_color, fontsize=14)

    axes[0].set_ylabel("Reco")
    cbar = fig.colorbar(im, ax=axes.ravel().tolist())
    cbar.set_label("Percentage (%)")

    # --- CRITICAL FIX: Ensure this is indented inside the function ---
    return fig