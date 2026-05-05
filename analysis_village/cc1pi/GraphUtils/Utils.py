bkg_name_nice_map = {
    "no_match": "Not matched",
    "out_AV_nu": "out-AV $\\nu$",
    "cosmic": "cosmic",
    "CC1pi": "$\\nu_{\\mu}$CC1$\\pi^{\\pm}$", # Check if your df uses 'CC1pi' or 'CC1Pi'
    "CC_mu_0pi_1p": "$\\nu_{\\mu}$CC0$\\pi^{\\pm}$1p",
    "CC_mu_0pi_2p": "$\\nu_{\\mu}$CC0$\\pi^{\\pm}2^{+}p$",
    "CC_mu_0pi_0p": "$\\nu_{\\mu}$CC0$\\pi^{\\pm}$0p",
    "CC_mu_0pi": "$\\nu_{\\mu}$CC0$\\pi^{\\pm}$",
    "CC_mu_2pi": "$\\nu_{\\mu}$CC2$^{+}\\pi^{\\pm}$",
    "NC": "NC",
    "other_nu": "Other $\\nu$ events",
    "other_CC1pi": "Excluded 1$\\mu$1$\\pi^{\\pm}$",
    "CC_e": "$\\nu_{e}$CC",
    "0p_CC1Pi": "$\\nu_{\\mu}$CC1$\\pi^{\\pm}$0p",
    "1p_CC1Pi": "$\\nu_{\\mu}$CC1$\\pi^{\\pm}$1p",
    "plus2p_CC1Pi": "$\\nu_{\\mu}$CC1$\\pi^{\\pm}$2p",
    "nu_mu_CC_Res": "$\\nu_{\\mu}$CC Res",
    "nu_mu_CC_QE": "$\\nu_{\\mu}$CC QE",
    "nu_mu_CC_Dis": "$\\nu_{\\mu}$CC DIS",
    "nu_mu_CC_MEC": "$\\nu_{\\mu}$CC MEC",
    "nu_mu_NC": "$\\nu_{\\mu}$NC",
    "bkg": "Background",
    "other": "other",
    
}

cut_name_nice_map = {
    "cosmic_rejection": "Cosmic rejection",
    "cosmic": "Clear cosmic",
    "t0": "Timing",
    "FV": "FV",
    "nu_score": "Nu Score Cut",
    "track": "2 Tracks",
    "shower": "No shower",
    "chi2": "2 MIP candidates",
    "angle": "Angle restriction",
    "proton_BDT": "HE proton removal",
    "containment": "No exiting particles",
    "michel": "Michel removal",
    "extra_pion": "Extra pion removal",
    "energy": "Energy",
    "MIP_refinement": "Michel & extra pion removal",
    "TPC_containment": "Local TPC Containment",
}

category_colors = {
    "CC1pi":        "#1f77b4", # Blue
    "CC_mu_0pi_1p": "#d62728", # Red
    "NC":           "#17becf", # Cyan-ish
    "CC_mu_0pi_0p": "#ffbb78", # Orange-Yellow
    "CC_mu_0pi_2p": "#2ca02c", # Green
    "out_AV_nu":    "#bcbd22", # Olive/Yellow-green
    "cosmic":       "#ff7f0e", # Orange
    "other_CC1pi":  "#7f7f7f", # Grey
    "CC_mu_2pi":    "#e377c2", # Pink
    "CC_e":         "#756bb1", # Light Green
}


genie_category_colors = { 
    "other": "#e377c2",
    "nu_mu_CC_Res": "#1f77b4",
    "nu_mu_CC_QE": "#d62728",
    "nu_mu_CC_Dis": "#2ca02c",
    "nu_mu_CC_MEC": "#7f7f7f",
    "nu_mu_NC": "#ffbb78",
    "out_AV_nu":    "#bcbd22", # Olive/Yellow-green
    "cosmic":       "#ff7f0e", # Orange
}


proton_distinction_category_colors = {
    "0p_CC1Pi":        "#1f77b4", # Blue
    "1p_CC1Pi": "#d62728", # Red
    "plus2p_CC1Pi": "#2ca02c", # Green
    "CC_mu_2pi":           "#17becf", # Cyan-ish
    "CC_mu_0pi": "#ffbb78", # Orange-Yellow
    "out_AV_nu":    "#bcbd22", # Olive/Yellow-green
    "cosmic":       "#ff7f0e", # Orange
    "other_nu":  "#7f7f7f", # Grey
}

category_colors_pfp = {
    "muon":        "#1f77b4", # Blue
    "stopping pion": "#d62728", # Red
    "inelastic pion": "#2ca02c", # Green
    "proton":       "#ff7f0e", # Orange
    "other": "#7f7f7f", # Orange-Yellow
    "shower":    "#e377c2", # Olive/Yellow-green
}


topology_list = ["CC1pi",       
    "CC_mu_0pi_0p",       
    "CC_mu_0pi_1p",      
    "CC_mu_0pi_2p",   
    "CC_mu_2pi",       
    "other_CC1pi",      
    "NC",               
    "out_AV_nu",        
    "cosmic",
    "CC_e"]    


topology_labels = [bkg_name_nice_map[topo] for topo in topology_list]

topology_color_dict = {
    "CC1pi":        "#1f77b4", # Blue
    "CC_mu_0pi_1p": "#d62728", # Red
    "NC":           "#17becf", # Cyan-ish
    "CC_mu_0pi_0p": "#ffbb78", # Orange-Yellow
    "CC_mu_0pi_2p": "#2ca02c", # Green
    "out_AV_nu":    "#bcbd22", # Olive/Yellow-green
    "cosmic":       "#ff7f0e", # Orange
    "other_CC1pi":  "#7f7f7f", # Grey
    "CC_mu_2pi":    "#e377c2", # Pink
    "CC_e":         "#756bb1", # Light Green
}

topology_colors = [topology_color_dict[topo] for topo in topology_list]

