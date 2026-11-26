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
    "cosmic": "Initial sample",
    "t0": "Timing",
    "FV": "FV",
    "nu_score": "Nu Score Cut",
    "track": "2 Tracks",
    "shower": "No shower",
    "chi2": "2 MIP candidates",
    "angle": "Angle restriction",
    "proton_BDT": "HE proton removal",
    "containment": "No particles exiting or in high-yz",
    "michel": "Michel removal",
    "extra_pion": "Extra pion removal",
    "energy": "Kinematic constraints",
    "MIP_refinement": "Michel & extra pion removal",
    "TPC_containment": "Local TPC Containment",
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




# Paul Tol's Muted 10-Color Palette (Left to Right from Image)
TOL_MUTED = [
    "#88CCEE",  # 1. Pale Cyan (Light sky blue)
    "#44AA99",  # 2. Teal (Mid-tone cool teal)
    "#117733",  # 3. Dark Green (Deep forest accent)
    "#332288",  # 4. Indigo (Deep purple-blue baseline)
    "#DDCC77",  # 5. Sand / Pale Yellow (Warm neutral)
    "#999933",  # 6. Olive Green (Muted yellow-green)
    "#CC6677",  # 7. Rose / Dusty Pink (Soft warm accent)
    "#882255",  # 8. Wine / Maroon (Deep reddish-purple)
    "#AA4499",  # 9. Purple / Amethyst (Vibrant muted purple)
    "#DDDDDD"   # 10. Light Gray (Subtle control group/baseline)
]

category_colors_old = {
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

category_colors = {
    "CC1pi":        "#0072B2", # 5. Blue
    "CC_mu_0pi_1p": "#D55E00", # 6. Vermilion
    "NC":           "#009E73", # 3. Bluish Green
    "CC_mu_0pi_0p":  "#CC79A7",  # 7. Reddish Purple
    "CC_mu_0pi_2p":   "#F5C5A3",  # 8.Peach
    "out_AV_nu":    "#6C6C6C",  # 10. Mid Gray ← new    
    "cosmic":       "#E69F00",  # 1. Orange
    "other_CC1pi": "#56B4E9",  # 2. Sky Blue
    "CC_mu_2pi":   "#994F00",  # 9. Brown ← new
    "CC_e":         "#F0E442",  # 4. Yellow 
}



topology_color_dict = category_colors
topology_colors = [topology_color_dict[topo] for topo in topology_list]




genie_category_colors = { 
    "other":   "#F0E442",  # 4. Yellow 
    "nu_mu_CC_Res":  "#0072B2", # 5. Blue
    "nu_mu_CC_QE": "#D55E00", # 6. Vermilion
    "nu_mu_CC_Dis": "#56B4E9",  # 2. Sky Blue
    "nu_mu_CC_MEC":"#CC79A7",  # 7. Reddish Purple
    "nu_mu_NC": "#009E73", # 3. Bluish Green
    "out_AV_nu":   "#6C6C6C",  # 10. Mid Gray ← new  
    "cosmic":   "#E69F00",  # 1. Orange
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


# Paul Tol's Muted 10-Color Palette (Left to Right from Image)
OKABE_ITO = [
    "#E69F00",  # 1. Orange
    "#56B4E9",  # 2. Sky Blue
    "#009E73",  # 3. Bluish Green
    "#F0E442",  # 4. Yellow
    "#0072B2",  # 5. Blue
    "#D55E00",  # 6. Vermilion
    "#CC79A7",  # 7. Reddish Purple
    "#F5C5A3",  # 8.Peach
    "#994F00",  # 9. Brown ← new
    "#6C6C6C",  # 10. Mid Gray ← new
]



category_colors_pfp = {
    "muon":        "#0072B2", # 5. Blue
    "stopping pion": "#D55E00",  # 6. Vermilion
    "inelastic pion": "#009E73", # 3. Bluish Green
    "proton":    "#E69F00",  # 1. Orange
    "other": "#6C6C6C",  # 10. Mid Gray ← new    
    "shower":    "#CC79A7",  # 7. Reddish Purple
}

category_colors_pfp_old = {
    "muon":        "#1f77b4", # Blue
    "stopping pion": "#d62728", # Red
    "inelastic pion": "#2ca02c", # Green
    "proton":       "#ff7f0e", # Orange
    "other": "#7f7f7f", # Orange-Yellow
    "shower":    "#e377c2", # Olive/Yellow-green
}



