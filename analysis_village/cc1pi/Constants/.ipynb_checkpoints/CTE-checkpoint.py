import numpy as np
# Geometry constants
# Detector geometry
y_min = -200.
y_max = +200.
x_min = -200.
x_max = +200.
z_min = 0.
z_max = 500.

# Masses of proton and muon (GeV)
mass_p = 0.9382720813
mass_mu = 0.1056583745
mass_pi = 0.13957018

# cuts constants

min_bc_score = 0.03
min_nu_score = 0.5

min_track_score = 0.5
min_track_lenght = 3
max_primary_distance_to_vertex = 10


min_distance_to_CPA = 0
min_distance_to_wall_x_y = 10
min_distance_to_first_z_wall = 10
min_distance_to_last_z_wall = 20

min_distance_to_consider_contained = 5
min_z_for_escaping_p = 5

max_shower_track_score = 0.5
min_shower_ke = 0.055

MIP_candidate_min_TL = 10
MIP_candidate_max_muon_score = 20
MIP_candidate_min_proton_score = 85

max_angle_between_candidates = 2.65

michel_max_track_score = 0.62
michel_max_visible_energy = 45
michel_max_dEdx_mean = 200

TMVA_BDT_proton_max_score = -0.07
BDT_proton_max_score = 0.62
clear_pion_BDT_score = 0


#final_states = ["out_AV_nu", "cosmic", "CC1pi", "CC_mu_0pi_1p", "CC_mu_0pi_2p", "CC_mu_0pi_0p",  "CC_mu_2pi", "NC", "other_CC1pi","CC_e"]
final_states = ["out_AV_nu", "cosmic", "CC1pi", "CC_mu_0pi_1p", "CC_mu_0pi_2p", "CC_mu_0pi_0p",  "CC_mu_2pi", "NC", "other_CC1pi","CC_e"]


extended_final_states = ["cosmic", "out_AV_nu", "0p_CC1Pi", "1p_CC1Pi", "plus2p_CC1Pi", "CC_mu_0pi", "CC_mu_2pi", "other_nu"]

# Colors for the plots in which I show different backgrounds
BCKG_COLORS = ['steelblue',  # Blue            
          'yellow',  # Yellow
          'moccasin',  # Orange 1
          'lime',  # Green
          'violet',  # Purple
          'peru',  # brown
          'pink',  # Pink
          'tomato',  # Red
          'orange',  # Orange 2
          'grey']  # Gray

# Colors for the plots in which I show different interaction modes
MODE_COLORS = ['steelblue',  # Blue 
               'lime',  # Green
               'tomato',  # Red    
               'moccasin',  # Orange
               'violet',  # Purple
               'pink',  # Pink
               'grey']  # Gray

MODE_LABELS = [r"QE", r"RES", r"MEC", r"DIS", r"MC Cosmics", "Others", "InTimeCosmics"]