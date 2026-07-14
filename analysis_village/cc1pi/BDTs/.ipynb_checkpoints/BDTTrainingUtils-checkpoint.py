import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc
from scipy.stats import ks_2samp
import pickle

variables_name_map = {
    "pfp_trk_chi2pid_best_chi2_muon": r"$\chi^2_{\mu}$",
    "pfp_trk_chi2pid_best_chi2_proton": r"$\chi^2_{p}$",
    "pfp_trk_chi2_exp_pol": r"$\chi^2_{pol_0} / \chi^2_{exp}$",
    "pfp_trk_frac50": r"RR frac. 50% E",
    "pfp_max_daughter_hits": r"daughter max hits",
    "pfp_scatter_angle_ratio": r" $\frac{MCS\ max\ scatter}{MCS\ total\ scatter}$",
}

from sklearn.base import clone
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, QuantileTransformer
from sklearn.decomposition import PCA
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import AdaBoostClassifier, GradientBoostingClassifier, BaggingClassifier, RandomForestClassifier
from xgboost import XGBClassifier
from sklearn.model_selection import train_test_split
from sklearn.datasets import make_classification

from sklearn.metrics import confusion_matrix

def bdt_quality_mask(df, columns):
    mask = np.ones(len(df), dtype=bool)

    for col in columns:
        mask &= df[col].notna()
        mask &= df[col] >= 0

    return mask



def create_models_for_columns(input_cols):
    model_vec = {}

    def wrap_in_pipeline(clf):
        return Pipeline([
            ('classifier', clf)
        ])
    
    # --- BDT (AdaBoost) ---
    bdt_clf = AdaBoostClassifier(
        estimator=DecisionTreeClassifier(max_depth=3, min_samples_leaf=0.025, criterion="gini"),
        n_estimators=850, learning_rate=0.5, random_state=40
    )
    model_vec["BDT"] = wrap_in_pipeline(bdt_clf)

    # --- BDTG (Gradient Boost) ---

    bdtg_clf = GradientBoostingClassifier(
        n_estimators=850, learning_rate=0.1, max_depth=3, 
        subsample=0.5, min_samples_leaf=0.025, random_state=40
    )

    model_vec["BDTG"] = wrap_in_pipeline(bdtg_clf)

    # --- XGBoost ---
    xgb_clf = XGBClassifier(
        n_estimators=850,
        learning_rate=0.1,
        max_depth=3,
        subsample=0.5,
        colsample_bytree=0.8,      # Similar to RF max_features
        tree_method='hist',        # Fast histogram-based method
        random_state=40,
        n_jobs=-1                  # Use all cores
    )
    model_vec["XGB"] = wrap_in_pipeline(xgb_clf)

    # --- BDTB (Bagging) ---
    bdtb_clf = BaggingClassifier(
        estimator=DecisionTreeClassifier(max_depth=3, min_samples_leaf=0.025, criterion="gini"),
        n_estimators=400, bootstrap=True, random_state=40, n_jobs=-1
    )
    model_vec["BDTB"] = wrap_in_pipeline(bdtb_clf)

    # --- Random Forest ---
    rf_clf = RandomForestClassifier(
        n_estimators=100, criterion='gini', max_depth=None, 
        max_features='sqrt', bootstrap=True, n_jobs=-1, random_state=40
    )
    model_vec["RF"] = wrap_in_pipeline(rf_clf)

    return model_vec


def train_all_models(signal_df, bkg_df, input_cols, model_vec, test_size=0.2, random_state=40):
    # 1. Build input arrays
    X_sig = signal_df[input_cols].values
    X_bkg = bkg_df[input_cols].values
    y_sig = np.ones(len(X_sig))
    y_bkg = np.zeros(len(X_bkg))
    X = np.vstack([X_sig, X_bkg])
    y = np.hstack([y_sig, y_bkg])
    # 2. Train/test split (Using stratify to maintain signal/bkg ratio)
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=test_size, stratify=y, random_state=random_state
    )
    n_train_sig = int(y_train.sum())
    n_train_bkg = int((y_train == 0).sum())
    n_test_sig  = int(y_test.sum())
    n_test_bkg  = int((y_test == 0).sum())

    print(f"Total samples   — Signal: {len(X_sig)},  Background: {len(X_bkg)}")
    print(f"Training samples — Signal: {n_train_sig}, Background: {n_train_bkg} (total: {len(X_train)})")
    print(f"Testing samples  — Signal: {n_test_sig},  Background: {n_test_bkg}  (total: {len(X_test)})")
    print(f"Signal fraction in train: {y_train.mean():.4f}")

    # 3. Train each model
    trained_models = {}
    for name, pipeline in model_vec.items():
        new_pipeline = clone(pipeline)
        new_pipeline.fit(X_train, y_train)
        trained_models[name] = new_pipeline
        print(f"Trained {name}")
        
    return trained_models, X_train, X_test, y_train, y_test
    

def get_scores(clf, X_train, X_test):
    """
    Safely get the BDT response/score for any sklearn classifier.
    Uses decision_function if available (TMVA-like), otherwise predict_proba.
    """
    # Try decision_function first (outputs values usually centered around 0)
    if hasattr(clf, "decision_function"):
        score_train = clf.decision_function(X_train)
        score_test  = clf.decision_function(X_test)
    # Fallback to predict_proba (outputs values between 0 and 1)
    else:
        score_train = clf.predict_proba(X_train)[:, 1]
        score_test  = clf.predict_proba(X_test)[:, 1]
        
    return score_train, score_test


def plot_roc(y_true, scores, label):
    fpr, tpr, _ = roc_curve(y_true, scores)
    roc_auc = auc(fpr, tpr)
    plt.plot(fpr, tpr, label=f"{label} (AUC = {roc_auc:.3f})")
    return roc_auc


def ks_test(signal_train, signal_test, bkg_train, bkg_test):
    ks_signal = ks_2samp(signal_train, signal_test).pvalue
    ks_bkg    = ks_2samp(bkg_train, bkg_test).pvalue
    return ks_signal, ks_bkg
    
def plot_transformed_importance(pipeline, original_columns):
    # 1. Access the steps from the pipeline
    pca = pipeline.named_steps['pca']
    bdt = pipeline.named_steps['classifier']
    
    # 2. Get BDT importances (for the Principal Components)
    bdt_importances = bdt.feature_importances_
    
    # 3. Get PCA Loadings (How much each original variable contributes to each PC)
    # pca.components_ has shape (n_components, n_features)
    pca_loadings = np.abs(pca.components_) 
    
    # 4. Project BDT importance back to original features
    # We multiply the importance of each PC by the weights of the original variables in that PC
    feat_importance_original = np.dot(bdt_importances, pca_loadings)
    
    # 5. Format names and plot
    feature_names = ["_".join([str(c) for c in col if c != ""]) for col in original_columns]
    
    df_final = pd.DataFrame({
        "original_feature": feature_names,
        "total_contribution": feat_importance_original
    }).sort_values(by="total_contribution", ascending=True)

    plt.figure(figsize=(8, 5))
    plt.barh(df_final['original_feature'], df_final['total_contribution'], color='lightcoral', edgecolor='black')
    plt.xlabel("Integrated Importance (BDT + PCA)")
    plt.title("Importance Projected onto Original Variables")
    plt.tight_layout()
    plt.show()

    return df_final



    
def plot_correlation_matrix(df, columns, title=None):
    """
    Plot correlation matrix for a list of MultiIndex columns in a DataFrame using plt.imshow.

    MultiIndex column names are merged into single strings for labeling.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing the data.
    columns : list of tuples
        MultiIndex columns to compute correlation on.
    title : str, optional
        Plot title.
    """
    # 1️⃣ Select valid rows (no NaN, >=0)
    mask = np.ones(len(df), dtype=bool)
    for col in columns:
        mask &= df[col].notna()
        mask &= df[col] >= 0

    df_vars = df.loc[mask, columns]

    # 2️⃣ Merge MultiIndex names for labeling
    feature_names = ["_".join([str(c) for c in col if c != ""]) for col in columns]
    df_vars.columns = feature_names

    # 3️⃣ Compute correlation matrix
    corr = df_vars.corr().values

    # 4️⃣ Plot with plt.imshow
    plt.figure(figsize=(6,5))
    im = plt.imshow(corr, cmap="coolwarm", vmin=-1, vmax=1)
    plt.colorbar(im, label="Correlation")

    # 5️⃣ Set ticks and labels
    plt.xticks(ticks=np.arange(len(feature_names)), labels=feature_names, rotation=45, ha="right")
    plt.yticks(ticks=np.arange(len(feature_names)), labels=feature_names)

    # 6️⃣ Annotate values
    for i in range(len(feature_names)):
        for j in range(len(feature_names)):
            plt.text(j, i, f"{corr[i,j]:.2f}", ha="center", va="center", color="black")

    plt.title(title or "Correlation matrix")
    plt.tight_layout()
    plt.show()

    return corr



def plot_response(sig_train, sig_test, bkg_train, bkg_test, title="BDT Response", bins=30,
                   signal_label="Signal", bkg_label="Background", save_path=None):
    """
    Plots the BDT response for Signal and Background, comparing Train and Test.
    Replicates the TMVA 'Overtraining Check' control plot.
    """
    # 1. Calculate KS Test Stats (compares Train vs Test distribution)
    ks_sig = ks_2samp(sig_train, sig_test)
    ks_bkg = ks_2samp(bkg_train, bkg_test)
    # 2. Define range for the x-axis (0 to 1 for predict_proba, or auto for decision_function)
    all_scores = np.concatenate([sig_train, sig_test, bkg_train, bkg_test])
    x_range = (np.min(all_scores), np.max(all_scores))
    
    fig = plt.figure(figsize=(8, 6))
    
# --- Helper to calculate density error bars for the 'Test' points ---
    def get_points_and_errors(data, bins, r):
        counts, bin_edges = np.histogram(data, bins=bins, range=r)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        bin_width = bin_edges[1] - bin_edges[0]
        # Normalize to density
        density = counts / (len(data) * bin_width)
        errors = np.sqrt(counts) / (len(data) * bin_width)
        return bin_centers, density, errors
    # --- Plot Training Distributions (Filled Histograms) ---
    plt.hist(sig_train, bins=bins, range=x_range, density=True, 
             alpha=0.2, color='blue', label=f'{signal_label} (Train)', histtype='stepfilled')
    plt.hist(sig_train, bins=bins, range=x_range, density=True, 
             color='blue', histtype='step', lw=1.5)
    plt.hist(bkg_train, bins=bins, range=x_range, density=True, 
             alpha=0.2, color='red', label=f'{bkg_label} (Train)', histtype='stepfilled')
    plt.hist(bkg_train, bins=bins, range=x_range, density=True, 
             color='red', histtype='step', lw=1.5)
    # --- Plot Testing Distributions (Points with Error Bars) ---
    # Signal Test
    x_s, y_s, err_s = get_points_and_errors(sig_test, bins, x_range)
    plt.errorbar(x_s, y_s, yerr=err_s, fmt='o', color='blue', 
                 label=f'{signal_label} (Test), KS p={ks_sig.pvalue:.3f}', markersize=4)
    # Background Test
    x_b, y_b, err_b = get_points_and_errors(bkg_test, bins, x_range)
    plt.errorbar(x_b, y_b, yerr=err_b, fmt='o', color='red', 
                 label=f'{bkg_label} (Test), KS p={ks_bkg.pvalue:.3f}', markersize=4)
    # --- Final Touches ---
    plt.xlabel(title)
    plt.ylabel("A.U.")
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2, frameon=False)
    plt.grid(alpha=0.2, linestyle='--')
    plt.tight_layout()

    if save_path is not None:
        fig.savefig(save_path, format='pdf', bbox_inches='tight')

    plt.show()




def prepare_muon_pion_cm(df, score_col = ('pfp', 'BDT_score_muon_pion', '', '', '', '')):
    """
    Groups by slice, assigns reco labels by BDT score, 
    and returns a raw confusion matrix for the plotting function.
    """
    SLICE_LEVELS = ["__ntuple", "entry", "rec.slc..index"]
    pdg_col = ('pfp', 'trk', 'truth', 'p', 'pdg', '')

    # 1. Identify rows for Muon (highest score) and Pion (lowest score) per slice
    grouped = df[score_col].groupby(level=SLICE_LEVELS)
    
    idx_reco_muon = grouped.idxmax()
    idx_reco_pion = grouped.idxmin()

    # 2. Extract true PDGs
    true_pdg_muon = df.loc[idx_reco_muon, pdg_col].abs().values
    true_pdg_pion = df.loc[idx_reco_pion, pdg_col].abs().values

    # 3. Create the 'Matched' (Truth) and 'Reco' arrays
    # Row index in CM will be RECO, Column index will be TRUTH
    y_true = np.concatenate([true_pdg_muon, true_pdg_pion])
    y_reco = np.concatenate([["Muon"]*len(true_pdg_muon), ["Pion"]*len(true_pdg_pion)])

    # Map PDG to string labels
    mapping = {13: "Muon", 211: "Pion"}
    y_true_labels = [mapping.get(code, "Other") for code in y_true]

    # 4. Generate the raw confusion matrix
    # Labels order defines the [0,0], [0,1] etc. positions
    labels = ["Muon", "Pion", "Other"]
    # confusion_matrix(y_true, y_pred) puts Truth on rows, Pred on cols. 
    # To match your plot (Reco on rows), we swap the arguments:
    cm_raw = confusion_matrix(y_reco, y_true_labels, labels=labels)
    
    return cm_raw, labels





def plot_correlation_matrices(signal_df, bkg_df, columns,
                               signal_title=r"$\mu/\pi$",
                               bkg_title="proton",
                               name_map=None, save_folder = None):
    """
    Plot side-by-side correlation matrices for signal and background DataFrames,
    for a list of MultiIndex columns, using plt.imshow.

    Parameters
    ----------
    signal_df : pd.DataFrame
        Signal DataFrame.
    bkg_df : pd.DataFrame
        Background DataFrame.
    columns : list of tuples
        MultiIndex columns to compute correlation on.
    signal_title : str, optional
        Title for the signal subplot.
    bkg_title : str, optional
        Title for the background subplot.
    name_map : dict, optional
        Mapping from raw merged column names to display labels (e.g. LaTeX).
        Falls back to the raw name if not found.

    Returns
    -------
    corr_signal, corr_bkg : np.ndarray
        Correlation matrices for signal and background.
    """
    name_map = name_map or {}

    def compute_corr(df, columns):
        mask = np.ones(len(df), dtype=bool)
        for col in columns:
            mask &= df[col].notna()
            mask &= df[col] >= 0
        df_vars = df.loc[mask, columns]
        feature_names = ["_".join([str(c) for c in col if c != ""]) for col in columns]
        df_vars.columns = feature_names
        return df_vars.corr().values, feature_names

    corr_signal, feature_names = compute_corr(signal_df, columns)
    corr_bkg, _ = compute_corr(bkg_df, columns)

    display_names = [name_map.get(name, name) for name in feature_names]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    for k, (ax, corr, title) in enumerate(zip(axes, [corr_signal, corr_bkg], [signal_title, bkg_title])):
        im = ax.imshow(corr, cmap="coolwarm", vmin=-1, vmax=1)
        ax.set_xticks(np.arange(len(display_names)))
        ax.set_xticklabels(display_names, rotation=45, ha="right")
        ax.set_yticks(np.arange(len(display_names)))
        if k == 0:
            ax.set_yticklabels(display_names)
        else:
            ax.set_yticklabels([])
        for i in range(len(display_names)):
            for j in range(len(display_names)):
                ax.text(j, i, f"{corr[i,j]:.2f}", ha="center", va="center", color="black")
        ax.set_title(title)

    fig.colorbar(im, ax=axes, label="Correlation", shrink=0.8)
    
    if save_folder is not None:
        fig.savefig(save_folder, format='pdf', bbox_inches='tight')
        
    plt.show()
        
    return corr_signal, corr_bkg


def plot_importance_no_transform(pipeline, original_columns, name_map=None, save_path=None):
    name_map = name_map or {}
    # 1. Access the classifier directly from the pipeline
    # We assume the last step is named 'classifier'
    bdt = pipeline.named_steps['classifier']
    
    # 2. Get BDT importances directly
    # Since there's no PCA, these align 1-to-1 with your input features
    feat_importance = bdt.feature_importances_
    
    # 3. Format names 
    # This handles the multi-index tuples you have in your LArSoft dataframes
    feature_names = ["_".join([str(c) for c in col if c != ""]) for col in original_columns]
    display_names = [name_map.get(name, name) for name in feature_names]
    
    # 4. Create DataFrame and sort
    df_final = pd.DataFrame({
        "original_feature": feature_names,
        "display_feature": display_names,
        "importance": feat_importance
    }).sort_values(by="importance", ascending=True)
    
    # 5. Plotting
    fig = plt.figure(figsize=(10, 6))
    plt.barh(df_final['display_feature'], df_final['importance'], 
             color='skyblue', edgecolor='black')
    
    plt.xlabel("Feature Importance")
    plt.ylabel("Variables")
    plt.grid(axis='x', linestyle='--', alpha=0.7)
    plt.tight_layout()
    
    if save_path is not None:
        fig.savefig(save_path, format='pdf', bbox_inches='tight')
        
    plt.show()
    
    return df_final