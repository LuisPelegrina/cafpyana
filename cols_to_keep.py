# config_columns.py
import pandas as pd

def get_truth_slim_cols():
    # 1. Base physics columns
    cols = [
        ('E', '', '', '', ''),
        ('pdg', '', '', '', ''),
        ('iscc', '', '', '', ''),
        ('genie_mode', '', '', '', ''),
        ('genweight', '', '', '', ''),
        ('nprim', '', '', '', ''),
        ('nn', '', '', '', ''),
        ('np', '', '', '', ''),
        ('nmu', '', '', '', ''),
        ('npi', '', '', '', ''),
        ('nmu_P_100MeV_3000MeV', '', '', '', ''),
        ('npi_P_85MeV_10000MeV', '', '', '', ''),
        ('npi_P_130MeV_10000MeV', '', '', '', ''),
        ('npi_P_130MeV_2000MeV', '', '', '', ''),
        ('np_P_325MeV_10000MeV', '', '', '', ''),
        ('mu', 'genp', 'x', '', ''),
        ('mu', 'genp', 'y', '', ''),
        ('mu', 'genp', 'z', '', ''),
        ('cpi', 'genp', 'x', '', ''),
        ('cpi', 'genp', 'y', '', ''),
        ('cpi', 'genp', 'z', '', ''),
        ('mu', 'totp', '', '', ''),
        ('true_var', 'true_cos_theta_mu', '', '', ''),
        ('true_var', 'true_cos_theta_pi', '', '', ''),
        ('true_var', 'true_mu_pi_angle', '', '', ''),
        ('true_var', 'true_p_mu', '', '', ''),
        ('true_var', 'true_p_pi', '', '', ''),
        ('nu_categ', '', '', '', ''),
        ('nu_categ_proton_reduced', '', '', '', ''),
        ('genie_categ', '', '', '', ''),
        ('true_var', 'num_protons', '', '', ''),
        ('true_var', 'delta_pT', '', '', ''),
        ('true_var', 'delta_alpha_T', '', '', ''),
        ('true_var', 'delta_phi_T', '', '', '')
    ]

    # 2. Add Systematic Universes
    for label in ['Flux', 'GENIE', 'G4']:
        cols.extend([(label, f'univ_{i}', '', '', '') for i in range(100)])
        
    return cols

# Create the final variable
truth_cols_to_keep_slim = get_truth_slim_cols()

min_truth_cols_to_keep= [
        ('E', '', '', '', ''),
        ('pdg', '', '', '', ''),
        ('iscc', '', '', '', ''),
        ('genie_mode', '', '', '', ''),
        ('genweight', '', '', '', ''),
        ('nprim', '', '', '', ''),
        ('nn', '', '', '', ''),
        ('np', '', '', '', ''),
        ('nmu', '', '', '', ''),
        ('npi', '', '', '', ''),
        ('nmu_P_100MeV_3000MeV', '', '', '', ''),
        ('npi_P_85MeV_10000MeV', '', '', '', ''),
        ('npi_P_130MeV_10000MeV', '', '', '', ''),
        ('npi_P_130MeV_2000MeV', '', '', '', ''),
        ('np_P_325MeV_10000MeV', '', '', '', ''),
        ('mu', 'genp', 'x', '', ''),
        ('mu', 'genp', 'y', '', ''),
        ('mu', 'genp', 'z', '', ''),
        ('cpi', 'genp', 'x', '', ''),
        ('cpi', 'genp', 'y', '', ''),
        ('cpi', 'genp', 'z', '', ''),
        ('mu', 'totp', '', '', ''),
        ('true_var', 'true_cos_theta_mu', '', '', ''),
        ('true_var', 'true_cos_theta_pi', '', '', ''),
        ('true_var', 'true_mu_pi_angle', '', '', ''),
        ('true_var', 'true_p_mu', '', '', ''),
        ('true_var', 'true_p_pi', '', '', ''),
        ('nu_categ', '', '', '', ''),
        ('nu_categ_proton_reduced', '', '', '', ''),
        ('genie_categ', '', '', '', ''),
        ('true_var', 'num_protons', '', '', ''),
        ('true_var', 'delta_pT', '', '', ''),
        ('true_var', 'delta_alpha_T', '', '', ''),
        ('true_var', 'delta_phi_T', '', '', '')
    ]


reco_cols_to_keep = [
        ('slc', 'tmatch', 'idx', '', '', ''),
        ('slc', 'self', '', '', '', ''),
        
        #Cosmic cut optimization
        ('slc', 'nu_score', '', '', '', ''),
        ('slc', 'barycenterFM', 'score', '', '', ''),
        ('slc', 'barycenterFM', 'flashTime', '', '', ''),
        
        #Vertex position        
        ('slc', 'vertex', 'x', '', '', ''),
        ('slc', 'vertex', 'y', '', '', ''),
        ('slc', 'vertex', 'z', '', '', ''),

        #Pfp truth info
        ('pfp', 'trk', 'truth', 'p', 'p_type', ''),
        ('pfp', 'trk', 'truth', 'p', 'pdg', ''),
        ('pfp', 'trk', 'truth', 'p', 'end_process', ''),
        
        #track optimization
        ('pfp', 'trk', 'len', '', '', ''),
        ('pfp', 'trackScore', '', '', '', ''),
        ('pfp', 'dist_to_vertex', '', '', '', ''),
        ('pfp', 'parent_is_primary', '', '', '', ''),
        
        #shower optimization
        ('pfp', 'shw', 'bestplane_energy', '', '', ''),
        
        #chi2 optimization
        ('pfp', 'trk', 'chi2pid', 'best', 'chi2_muon', ''),
        ('pfp', 'trk', 'chi2pid', 'best', 'chi2_proton', ''),
        
        #michel cut
        ('pfp', 'trk', 'calo', 'best', 'ke', ''),
        ('pfp', 'trk', 'mean_dEdx', '', '', ''),
        ('pfp', 'max_daughter_hits', '', '', '', ''),

        #angle cut
        ('slc', 'measure_var', 'angle_between_candidates', '', '', ''),
    
        #BDT vars
        ('pfp', 'trk', 'chi2_exp_pol', '', '', ''),
        ('pfp', 'trk', 'frac50', '', '', ''),
        ('pfp', 'scatter_angle_ratio', '', '', '', ''),

        #pfp containment masks
        ('pfp', 'is_exiting', '', '', '', ''),
        
        #BDT columns
        ('pfp','trk','bdt_muon_pion_score','','',''),
        ('pfp','trk','bdt_proton_score','','',''),
        
        #cuts boolean
        ('slc', 'cut', 'obvious_cosmic', '', '', ''),
        ('slc', 'cut', 't0', '', '', ''),
        ('slc', 'cut', 'inside_FV', '', '', ''),
        ('slc', 'cut', 'nu_score', '', '', ''),
        ('slc', 'cut', 'track', '', '', ''),
        ('slc', 'cut', 'shower', '', '', ''),
        ('slc', 'cut', 'MIP_candidates', '', '', ''),
        ('slc', 'cut', 'angle', '', '', ''),
        ('slc', 'cut', 'proton_BDT', '', '', ''),
        ('slc', 'cut', 'proton_BDT_sideband', '', '', ''),
        ('slc', 'cut', 'proton_BDT_2pi', '', '', ''),
        ('slc', 'cut', 'TPC_containment', '', '', ''),
        ('slc', 'cut', 'containment', '', '', ''),
        ('slc', 'cut', 'michel', '', '', ''),
        ('slc', 'cut', 'extra_pion', '', '', ''),
        ('slc','cut','energy','','',''),
        ('slc','cut','no_high_yz','','',''),

        #Measure variables
        ('slc', 'measure_var', 'num_protons', '', '', ''),
        ('slc','measure_var','reco_p_mu','','',''),
        ('slc','measure_var','reco_cos_theta_mu','','',''),
        ('slc','measure_var','mu_chi2_proton','','',''),
        ('slc','measure_var','mu_chi2_mu','','',''),
        ('slc','measure_var','mu_chi2_exp_pol','','',''),
        ('slc','measure_var','mu_scatter_angle_ratio','','',''),
        ('slc','measure_var','mu_max_daughter_hits','','',''),
        ('slc','measure_var','mu_frac50','','',''),
        ('slc','measure_var','mu_bdt_score_proton','','',''),
        ('slc','measure_var','mu_bdt_score_muon_pion','','',''),
        ('slc','measure_var','range_p_pi','','',''),
        ('slc','measure_var','TLE_p_pi','','',''),
        ('slc','measure_var','reco_cos_theta_pi','','',''),      
        ('slc','measure_var','pi_chi2_proton','','',''),
        ('slc','measure_var','pi_chi2_mu','','',''),
        ('slc','measure_var','pi_chi2_exp_pol','','',''),
        ('slc','measure_var','pi_scatter_angle_ratio','','',''),
        ('slc','measure_var','pi_max_daughter_hits','','',''),
        ('slc','measure_var','pi_frac50','','',''),
        ('slc','measure_var','pi_bdt_score_proton','','',''),
        ('slc','measure_var','pi_bdt_score_muon_pion','','',''),
      
        ('slc','measure_var','mu_true_pdg','','',''),
        ('slc','measure_var','mu_true_p_type','','',''),
        ('slc','measure_var','mu_true_end_process','','',''),
        ('slc','measure_var','pi_true_pdg','','',''),
        ('slc','measure_var','pi_true_p_type','','',''),
        ('slc','measure_var','pi_true_end_process','','',''),
        ('slc','measure_var','pi_true_p','','',''),
        ('slc','measure_var','mu_true_p','','',''),
        ('slc','measure_var','pi_true_costheta','','',''),
        ('slc','measure_var','mu_true_costheta','','',''), 
        ('slc','measure_var','muon_contained','','',''),  

        #Cut vars
        ('slc','cut_var','n_prim_tracks','','',''), 
        ('slc','cut_var','n_prim_showers','','',''),
        ('slc','cut_var','n_MIP_candidates','','',''),
        ('slc','cut_var','n_MIP_candidates_proton','','',''),
        ('slc','cut_var','n_exiting_pfps','','',''),
        ('slc','cut_var','n_exiting_z_pfps','','',''),
        ('slc','cut_var','n_MIP_candidate_michel','','',''),
        ('slc','cut_var','n_extra_pions','','',''),
 
        ('slc','measure_var','delta_pT','','',''),
        ('slc','measure_var','delta_alpha_T','','',''),
        ('slc','measure_var','delta_phi_T','','',''),
    ]





min_reco_cols_to_keep = [
        ('slc', 'self', '', '', '', ''),
        ('slc', 'tmatch', 'idx', '', '', ''),
        ('slc', 'nu_score', '', '', '', ''),
        ('slc', 'cut', 'obvious_cosmic', '', '', ''),
        ('slc', 'cut', 't0', '', '', ''),
        ('slc', 'cut', 'inside_FV', '', '', ''),
        ('slc', 'cut', 'nu_score', '', '', ''),
        ('slc', 'cut', 'track', '', '', ''),
        ('slc', 'cut', 'shower', '', '', ''),
        ('slc', 'cut', 'MIP_candidates', '', '', ''),
        ('slc', 'cut', 'angle', '', '', ''),
        ('slc', 'cut', 'proton_BDT', '', '', ''),
        ('slc', 'cut', 'proton_BDT_sideband', '', '', ''),
        ('slc', 'cut', 'containment', '', '', ''),
        ('slc', 'cut', 'proton_BDT_2pi', '', '', ''),
        ('slc', 'cut', 'TPC_containment', '', '', ''),
        ('slc', 'cut', 'michel', '', '', ''),
        ('slc', 'cut', 'extra_pion', '', '', ''),
        ('slc','cut','energy','','',''),
        ('slc','cut','no_high_yz','','',''),
        ('slc','cut_var','n_MIP_candidates','','',''),
        ('slc','cut_var','n_exiting_pfps','','',''),
        ('slc', 'measure_var', 'angle_between_candidates', '', '', ''),
        ('slc', 'measure_var', 'num_protons', '', '', ''),
        ('slc','measure_var','reco_p_mu','','',''),
        ('slc','measure_var','reco_cos_theta_mu','','',''),
        ('slc','measure_var','TLE_p_pi','','',''),
        ('slc','measure_var','range_p_pi','','',''),
        ('slc','measure_var','reco_cos_theta_pi','','',''), 
        ('slc','measure_var','delta_pT','','',''),
        ('slc','measure_var','delta_alpha_T','','',''),
        ('slc','measure_var','delta_phi_T','','',''),
        ('slc','measure_var','mu_true_p_type','','',''),
        ('slc','measure_var','pi_true_p_type','','',''),
        ('slc','measure_var','mu_true_pdg','','',''),
        ('slc','measure_var','pi_true_pdg','','',''),
        ('slc','measure_var','pi_true_p','','',''),
        ('slc','measure_var','mu_true_p','','',''),
        ('slc','measure_var','pi_true_costheta','','',''),
        ('slc','measure_var','mu_true_costheta','','',''),
        ('slc','measure_var','muon_contained','','',''),
]




reco_cols_to_keep_for_data_mc_comparison = [
        ('slc', 'tmatch', 'idx', '', '', ''),
        ('slc', 'self', '', '', '', ''),
        
        ('slc', 'nu_score', '', '', '', ''),
    
        #Pfp truth info
        ('pfp', 'trk', 'truth', 'p', 'p_type', ''),
        ('pfp', 'trk', 'truth', 'p', 'pdg', ''),
        ('pfp', 'trk', 'truth', 'p', 'end_process', ''),
        
        #track optimization
        ('pfp', 'trk', 'len', '', '', ''),
        ('pfp', 'trackScore', '', '', '', ''),
        ('pfp', 'dist_to_vertex', '', '', '', ''),
        ('pfp', 'parent_is_primary', '', '', '', ''),
        
        #chi2 optimization
        ('pfp', 'trk', 'chi2pid', 'best', 'chi2_muon', ''),
        ('pfp', 'trk', 'chi2pid', 'best', 'chi2_proton', ''),
        
        #angle cut
        ('slc', 'measure_var', 'angle_between_candidates', '', '', ''),
    
        #BDT vars
        ('pfp', 'trk', 'chi2_exp_pol', '', '', ''),
        ('pfp', 'trk', 'frac50', '', '', ''),
        ('pfp', 'scatter_angle_ratio', '', '', '', ''),

        #pfp containment masks
        ('pfp', 'is_exiting', '', '', '', ''),
        
        #BDT columns
        ('pfp','trk','bdt_muon_pion_score','','',''),
        ('pfp','trk','bdt_proton_score','','',''),
        
        #cuts boolean
        ('slc', 'cut', 'obvious_cosmic', '', '', ''),
        ('slc', 'cut', 't0', '', '', ''),
        ('slc', 'cut', 'inside_FV', '', '', ''),
        ('slc', 'cut', 'nu_score', '', '', ''),
        ('slc', 'cut', 'track', '', '', ''),
        ('slc', 'cut', 'shower', '', '', ''),
        ('slc', 'cut', 'MIP_candidates', '', '', ''),
        ('slc', 'cut', 'angle', '', '', ''),
        ('slc', 'cut', 'proton_BDT', '', '', ''),
        ('slc', 'cut', 'proton_BDT_sideband', '', '', ''),
        ('slc', 'cut', 'proton_BDT_2pi', '', '', ''),
        ('slc', 'cut', 'TPC_containment', '', '', ''),
        ('slc', 'cut', 'containment', '', '', ''),
        ('slc', 'cut', 'michel', '', '', ''),
        ('slc', 'cut', 'extra_pion', '', '', ''),
        ('slc','cut','energy','','',''),
        ('slc','cut','no_high_yz','','',''),

        #Measure variables
        ('slc', 'measure_var', 'num_protons', '', '', ''),
        ('slc','measure_var','reco_p_mu','','',''),
        ('slc','measure_var','reco_cos_theta_mu','','',''),
        ('slc','measure_var','range_p_pi','','',''),
        ('slc','measure_var','TLE_p_pi','','',''),
        ('slc','measure_var','reco_cos_theta_pi','','',''),
    
        ('slc','measure_var','mu_true_p_type','','',''),
        ('slc','measure_var','pi_true_p_type','','',''),
    
        #Cut vars
        ('slc','cut_var','n_prim_tracks','','',''), 
        ('slc','cut_var','n_prim_showers','','',''),
        ('slc','cut_var','n_MIP_candidates','','',''),
        ('slc','cut_var','n_MIP_candidates_proton','','',''),
        ('slc','cut_var','n_exiting_pfps','','',''),
        ('slc','cut_var','n_exiting_z_pfps','','',''),
        ('slc','cut_var','n_MIP_candidate_michel','','',''),
        ('slc','cut_var','n_extra_pions','','',''),
 
        ('slc','measure_var','delta_pT','','',''),
        ('slc','measure_var','delta_alpha_T','','',''),
        ('slc','measure_var','delta_phi_T','','',''),
    ]





reco_cols_to_keep_for_BDT_comparison = [
        ('slc', 'tmatch', 'idx', '', '', ''),
        ('slc', 'self', '', '', '', ''),
        
        ('slc', 'nu_score', '', '', '', ''),
    
        #Pfp truth info
        ('pfp', 'trk', 'truth', 'p', 'p_type', ''),
        ('pfp', 'trk', 'truth', 'p', 'pdg', ''),
        ('pfp', 'trk', 'truth', 'p', 'end_process', ''),
        
        #track optimization
        ('pfp', 'trk', 'len', '', '', ''),
        ('pfp', 'trackScore', '', '', '', ''),
        ('pfp', 'dist_to_vertex', '', '', '', ''),
        ('pfp', 'parent_is_primary', '', '', '', ''),
        
        #chi2 optimization
        ('pfp', 'trk', 'chi2pid', 'best', 'chi2_muon', ''),
        ('pfp', 'trk', 'chi2pid', 'best', 'chi2_proton', ''),
        
        #BDT vars
        ('pfp', 'trk', 'chi2_exp_pol', '', '', ''),
        ('pfp', 'trk', 'frac50', '', '', ''),
        ('pfp', 'scatter_angle_ratio', '', '', '', ''),

        #pfp containment masks
        ('pfp', 'is_exiting', '', '', '', ''),
        
        #BDT columns
        ('pfp','trk','bdt_muon_pion_score','','',''),
        ('pfp','trk','bdt_proton_score','','',''),

    
        ('pfp', 'trk', 'calo', 'best', 'ke', ''),
        ('pfp', 'trk', 'mean_dEdx', '', '', ''),
        ('pfp', 'max_daughter_hits', '', '', '', ''),
    
        #cuts boolean
        ('slc', 'cut', 'obvious_cosmic', '', '', ''),
        ('slc', 'cut', 't0', '', '', ''),
        ('slc', 'cut', 'inside_FV', '', '', ''),
        ('slc', 'cut', 'nu_score', '', '', ''),
        ('slc', 'cut', 'track', '', '', ''),
        ('slc', 'cut', 'shower', '', '', ''),
        ('slc', 'cut', 'MIP_candidates', '', '', ''),
        ('slc', 'cut', 'angle', '', '', ''),
        ('slc', 'cut', 'proton_BDT', '', '', ''),
        ('slc', 'cut', 'proton_BDT_sideband', '', '', ''),
        ('slc', 'cut', 'proton_BDT_2pi', '', '', ''),
        ('slc', 'cut', 'TPC_containment', '', '', ''),
        ('slc', 'cut', 'containment', '', '', ''),
        ('slc', 'cut', 'michel', '', '', ''),
        ('slc', 'cut', 'extra_pion', '', '', ''),
        ('slc','cut','energy','','',''),
        ('slc','cut','no_high_yz','','',''),
    
        ('slc','measure_var','reco_p_mu','','',''),
        ('slc','measure_var','TLE_p_pi','','',''),
    ]