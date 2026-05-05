import numpy as np
from dataclasses import dataclass, field
import pandas as pd
from analysis_village.cc1pi.CutMasks import CutMasks
from analysis_village.cc1pi.Constants import CTE as CTE

mask_dict = {
    "none": lambda df: pd.Series(True, index=df.index),
    "pandora_primary": lambda df: CutMasks.is_pandora_primary_mask(df),
    "analysis_primary": lambda df: CutMasks.is_analysis_primary_mask(df),
    "primary_track": lambda df: CutMasks.is_primary_track_mask(df),
    "primary_shower": lambda df: CutMasks.is_primary_shower_mask(df),
    "MIP_candidate": lambda df: CutMasks.is_MIP_candidate_mask(df),
    "contained_MIP_candidate": lambda df: CutMasks.is_MIP_candidate_mask(df) & (df.pfp.is_exiting == False),
    "exiting_MIP_candidate": lambda df: CutMasks.is_MIP_candidate_mask(df) & (df.pfp.is_exiting == True),
    "muon_exiting": lambda df: df.slc.measure_var.muon_contained == False,
    "muon_contained": lambda df: df.slc.measure_var.muon_contained == True,
    "0p": lambda df: df.slc.measure_var.num_protons == 0,
    "1p": lambda df: df.slc.measure_var.num_protons == 1,
    "2plusp": lambda df: df.slc.measure_var.num_protons > 1,
}


@dataclass
class FullHistogramConfig:
    file_name: str = None
    var_evt_reco_col: tuple = None
    truth_column: tuple = None
    first_per_slice: bool = False
    start_cut: str = "cosmic"
    end_cut: str = "energy"
    extra_mask: str = "none"
    bins: np.ndarray = field(default_factory=lambda: np.linspace(-0.5, 0.5, 51))
    xlabel: str = 'Score'
    ylabel: str = 'Entries'
    cut_value: any = field(default_factory=lambda: [-999.0])
    clip: bool = True
    stats_horizontal_alignment: str = 'right'

    # Will be automatically computed
    bin_centers: np.ndarray = field(init=False)

    def __post_init__(self):
        self.bin_centers = (self.bins[:-1] + self.bins[1:]) / 2.
            
slices_y_label = f'Candidate Slices'
pfps_y_label = f'PFPs'

nu_categ_column = ('truth','nu_categ','','','','')
pfp_truth_column = ('pfp','trk','truth','p','p_type','')
cuts = ["cosmic","t0","FV","nu_score","track","chi2","shower","angle","proton_BDT","containment","michel","extra_pion","energy"]



config_bc_flash_time = FullHistogramConfig(
    file_name = "bc_flash_time", 
    var_evt_reco_col=('slc','barycenterFM','flashTime','','',''),
    truth_column=nu_categ_column,
    first_per_slice = True,
    start_cut = "cosmic",
    end_cut = "FV",
    bins=np.linspace(-3, 5, 41),
    xlabel=r'BarycenterFM time [$\mu s$]',
    ylabel= slices_y_label,
    clip = False
)

config_bc_flash_score = FullHistogramConfig(
    file_name = "bc_flash_score", 
    var_evt_reco_col=('slc','barycenterFM','score','','',''),
    truth_column=nu_categ_column,
    first_per_slice = True,
    start_cut = "cosmic",
    end_cut = "FV",
    bins=np.linspace(0, 0.2, 41),
    xlabel=r'BarycenterFM score',
    ylabel= slices_y_label,
    cut_value = [CTE.min_bc_score],
    clip = False
)

config_nu_score = FullHistogramConfig(
    file_name = "nu_score", 
    var_evt_reco_col=('slc','nu_score','','','',''),
    truth_column=nu_categ_column,
    first_per_slice = True,
    start_cut = "cosmic",
    end_cut = "energy",
    bins=np.linspace(0, 1, 41),
    xlabel=r'Neutrino score',
    ylabel= slices_y_label,
    cut_value = [CTE.min_nu_score],
    clip = False
)

config_vtx_x = FullHistogramConfig(
    file_name = "vtx_x", 
    var_evt_reco_col=('slc','vertex','x','','',''),
    truth_column=nu_categ_column,
    first_per_slice = True,
    start_cut = "cosmic",
    end_cut = "energy",
    bins=np.linspace(-200, 200, 41),
    xlabel=r'Slice vertex x [cm]',
    cut_value = [-200 + CTE.min_distance_to_wall_x_y, 200 - CTE.min_distance_to_wall_x_y],
    ylabel= slices_y_label
)


config_vtx_y = FullHistogramConfig(
    file_name = "vtx_y", 
    var_evt_reco_col=('slc','vertex','y','','',''),
    truth_column=nu_categ_column,
    first_per_slice = True,
    start_cut = "cosmic",
    end_cut = "energy",
    bins=np.linspace(-200, 200, 41),
    xlabel=r'Slice vertex y [cm]',
    cut_value = [-200 + CTE.min_distance_to_wall_x_y, 200 - CTE.min_distance_to_wall_x_y],
    ylabel= slices_y_label
)

config_vtx_z = FullHistogramConfig(
    file_name = "vtx_z", 
    var_evt_reco_col=('slc','vertex','z','','',''),
    truth_column=nu_categ_column,
    first_per_slice = True,
    start_cut = "cosmic",
    end_cut = "energy",
    bins=np.linspace(0, 500, 41),
    xlabel=r'Slice vertex z [cm]',
    cut_value = [CTE.min_distance_to_first_z_wall, 500 - CTE.min_distance_to_last_z_wall],
    ylabel= slices_y_label
)

simple_vars_vec = [config_bc_flash_time, config_bc_flash_score, config_vtx_x, config_vtx_y, config_vtx_z,config_nu_score]
simple_vars_vec_no_flash_time = [config_bc_flash_score, config_vtx_x, config_vtx_y, config_vtx_z,config_nu_score]


####################

config_n_prim_tracks = FullHistogramConfig(
    file_name = "n_prim_tracks", 
    var_evt_reco_col=('slc','cut_var','n_prim_tracks','','',''),
    truth_column=nu_categ_column,
    first_per_slice = True,
    start_cut = "nu_score",
    end_cut = "shower",
    bins=np.linspace(0, 6, 7),
    xlabel=r'Number of primary tracks',
    ylabel= slices_y_label,
    cut_value = [2],
    clip = False
)


config_n_prim_showers = FullHistogramConfig(
    file_name = "n_prim_showers", 
    var_evt_reco_col=('slc','cut_var','n_prim_showers','','',''),
    truth_column=nu_categ_column,
    first_per_slice = True,
    start_cut = "chi2",
    end_cut = "shower",
    bins=np.linspace(0, 5, 6),
    xlabel=r'Number of primary showers',
    ylabel= slices_y_label,
    cut_value = [1],
    clip = False
)

config_n_MIP_candidates = FullHistogramConfig(
    file_name = "n_MIP_candidates", 
    var_evt_reco_col=('slc','cut_var','n_MIP_candidates','','',''),
    truth_column=nu_categ_column,
    first_per_slice = True,
    start_cut = "track",
    end_cut = "shower",
    bins=np.linspace(0, 5, 6),
    xlabel=r'Number of MIP candidates',
    ylabel= slices_y_label,
    cut_value = [2],
    clip = False
)

config_n_MIP_candidates_proton = FullHistogramConfig(
    file_name = "n_MIP_candidates_proton", 
    var_evt_reco_col=('slc','cut_var','n_MIP_candidates_proton','','',''),
    truth_column=nu_categ_column,
    first_per_slice = True,
    start_cut = "angle",
    end_cut = "containment",
    bins=np.linspace(0, 3, 4),
    xlabel=r'Number of proton like MIP candidates',
    ylabel= slices_y_label,
    cut_value = [2],
    clip = False
)

config_n_exiting_pfps = FullHistogramConfig(
    file_name = "n_exiting_pfps", 
    var_evt_reco_col=('slc','cut_var','n_exiting_pfps','','',''),
    truth_column=nu_categ_column,
    first_per_slice = True,
    start_cut = "proton_BDT",
    end_cut = "michel",
    bins=np.linspace(0, 5, 6),
    xlabel=r'Number of exiting pfps',
    ylabel= slices_y_label,
    cut_value = [2],
    clip = False
)

config_n_exiting_z_pfps = FullHistogramConfig(
    file_name = "n_exiting_z_pfps", 
    var_evt_reco_col=('slc','cut_var','n_exiting_z_pfps','','',''),
    truth_column=nu_categ_column,
    first_per_slice = True,
    start_cut = "proton_BDT",
    end_cut = "michel",
    bins=np.linspace(0, 5, 6),
    xlabel=r'Number of exiting pfps (z < 5 cm)',
    ylabel= slices_y_label,
    cut_value = [1],
    clip = False
)

config_n_MIP_candidate_michel = FullHistogramConfig(
    file_name = "n_MIP_candidate_michel", 
    var_evt_reco_col=('slc','cut_var','n_MIP_candidate_michel','','',''),
    truth_column=nu_categ_column,
    first_per_slice = True,
    start_cut = "containment",
    end_cut = "extra_pion",
    bins=np.linspace(0, 3, 4),
    xlabel=r'Number of Michel candidates',
    ylabel= slices_y_label,
    cut_value = [1],
    clip = False
)

config_n_extra_pions = FullHistogramConfig(
    file_name = "n_extra_pions", 
    var_evt_reco_col=('slc','cut_var','n_extra_pions','','',''),
    truth_column=nu_categ_column,
    first_per_slice = True,
    start_cut = "michel",
    end_cut = "energy",
    bins=np.linspace(0, 3, 4),
    xlabel=r'Number of extra pions',
    ylabel= slices_y_label,
    cut_value = [1],
    clip = False
)

n_pfps_config_vec = [
    config_n_prim_tracks,
    config_n_prim_showers,
    config_n_MIP_candidates,
    config_n_MIP_candidates_proton,
    config_n_exiting_pfps,
    config_n_exiting_z_pfps,
    config_n_MIP_candidate_michel,
    config_n_extra_pions
]


####################




#vtx distance
config_pandora_primary_vtx_distance = FullHistogramConfig(
    file_name = "pandora_primary_vtx_distance", 
    var_evt_reco_col=('pfp','dist_to_vertex','','','',''),
    truth_column=pfp_truth_column,
    first_per_slice = False,
    start_cut = "nu_score",
    end_cut = "nu_score",
    extra_mask = "pandora_primary",
    bins=np.linspace(0, 20, 41),
    xlabel='primary pfp distance to vtx [cm]',
    ylabel= pfps_y_label,
    cut_value = [CTE.max_primary_distance_to_vertex],
    clip = False
)

#track score
config_pandora_primary_track_score = FullHistogramConfig(
    file_name = "pandora_primary_track_score", 
    var_evt_reco_col=('pfp','trackScore','','','',''),
    truth_column=pfp_truth_column,
    first_per_slice = False,
    start_cut = "track",
    end_cut = "shower",
    extra_mask = "pandora_primary",
    bins=np.linspace(0, 1, 41),
    xlabel='primary pfp track score',
    ylabel= pfps_y_label,
    cut_value = [CTE.min_track_score],
    clip = False
)


config_analysis_primary_track_score = FullHistogramConfig(
    file_name = "analysis_primary_track_score", 
    var_evt_reco_col=('pfp','trackScore','','','',''),
    truth_column=pfp_truth_column,
    first_per_slice = False,
    start_cut = "nu_score",
    end_cut = "track",
    extra_mask = "analysis_primary",
    bins=np.linspace(0, 1, 41),
    xlabel='analysis primary track score (vtx dist inc.)',
    ylabel= pfps_y_label,
    cut_value = [CTE.min_track_score],
    clip = False
)

config_contained_MIP_candidate_track_score = FullHistogramConfig(
    file_name = "contained_MIP_candidate_track_score", 
    var_evt_reco_col=('pfp','trackScore','','','',''),
    truth_column=pfp_truth_column,
    first_per_slice = False,
    start_cut = "chi2",
    end_cut = "michel",
    extra_mask = "contained_MIP_candidate",
    bins=np.linspace(0, 1, 41),
    xlabel='contained MIP candidates track score',
    ylabel= pfps_y_label,
    cut_value = [CTE.michel_max_track_score],
    clip = False
)



#track KE
config_contained_MIP_candidate_KE = FullHistogramConfig(
    file_name = "contained_MIP_candidate_KE", 
    var_evt_reco_col=('pfp', 'trk', 'calo', 'best', 'ke', ''),
    truth_column=pfp_truth_column,
    first_per_slice = False,
    start_cut = "chi2",
    end_cut = "michel",
    extra_mask = "contained_MIP_candidate",
    bins=np.linspace(0, 100, 41),
    cut_value = [CTE.michel_max_visible_energy],
    xlabel='contained MIP candidates track KE [MeV]',
    ylabel= pfps_y_label,
    clip = False
)

#track len
config_analysis_primary_len = FullHistogramConfig(
    file_name = "analysis_primary_len", 
    var_evt_reco_col=('pfp','trk','len','','',''),
    truth_column=pfp_truth_column,
    first_per_slice = False,
    start_cut = "nu_score",
    end_cut = "track",
    extra_mask = "analysis_primary",
    bins=np.linspace(0, 100, 41),
    xlabel='primary pfp lenght [cm] (vtx dist inc.)',
    ylabel= pfps_y_label,
    cut_value = [CTE.min_track_lenght],
    clip = False
)

config_primary_track_len = FullHistogramConfig(
    file_name = "primary_track_len", 
    var_evt_reco_col=('pfp','trk','len','','',''),
    truth_column=pfp_truth_column,
    first_per_slice = False,
    start_cut = "nu_score",
    end_cut = "chi2",
    extra_mask = "primary_track",
    bins=np.linspace(0, 100, 41),
    xlabel='primary track lenght [cm]',
    ylabel= pfps_y_label,
    cut_value = [CTE.MIP_candidate_min_TL],
    clip = False
)


#Chi2mu
config_primary_track_chi2_mu = FullHistogramConfig(
    file_name = "primary_track_chi2_mu", 
    var_evt_reco_col=('pfp','trk','chi2pid','best','chi2_muon',''),
    truth_column=pfp_truth_column,
    first_per_slice = False,
    start_cut = "track",
    end_cut = "chi2",
    extra_mask = "primary_track",
    bins=np.linspace(0, 100, 41),
    xlabel=r'primary tracks $\chi^2_{\mu}$',
    ylabel= pfps_y_label,
    cut_value = [CTE.MIP_candidate_max_muon_score],
    clip = False
)

config_MIP_candidate_chi2_mu = FullHistogramConfig(
    file_name = "MIP_candidate_chi2_mu", 
    var_evt_reco_col=('pfp','trk','chi2pid','best','chi2_muon',''),
    truth_column=pfp_truth_column,
    first_per_slice = False,
    start_cut = "track",
    end_cut = "energy",
    extra_mask = "MIP_candidate",
    bins=np.linspace(0, 20, 41),
    xlabel=r'MIP candidates $\chi^2_{\mu}$',
    ylabel= pfps_y_label,
    cut_value = [CTE.MIP_candidate_max_muon_score],
    clip = False
)
config_contained_MIP_candidate_chi2_mu = FullHistogramConfig(
    file_name = "contained_MIP_candidate_chi2_mu", 
    var_evt_reco_col=('pfp','trk','chi2pid','best','chi2_muon',''),
    truth_column=pfp_truth_column,
    first_per_slice = False,
    start_cut = "track",
    end_cut = "energy",
    extra_mask = "contained_MIP_candidate",
    bins=np.linspace(0, 20, 41),
    xlabel=r'contained MIP candidates $\chi^2_{\mu}$',
    ylabel= pfps_y_label,
    cut_value = [CTE.MIP_candidate_max_muon_score],
    clip = False
)

config_exiting_MIP_candidate_chi2_mu = FullHistogramConfig(
    file_name = "exiting_MIP_candidate_chi2_mu", 
    var_evt_reco_col=('pfp','trk','chi2pid','best','chi2_muon',''),
    truth_column=pfp_truth_column,
    first_per_slice = False,
    start_cut = "track",
    end_cut = "energy",
    extra_mask = "exiting_MIP_candidate",
    bins=np.linspace(0, 20, 41),
    xlabel=r'exiting MIP candidates $\chi^2_{\mu}$',
    ylabel= pfps_y_label,
    cut_value = [CTE.MIP_candidate_max_muon_score],
    clip = False
)
    
#Chi2p
config_primary_track_chi2_p = FullHistogramConfig(
    file_name = "primary_track_chi2_p", 
    var_evt_reco_col=('pfp','trk','chi2pid','best','chi2_proton',''),
    truth_column=pfp_truth_column,
    first_per_slice = False,
    start_cut = "track",
    end_cut = "chi2",
    extra_mask = "primary_track",
    bins=np.linspace(0, 300, 41),
    xlabel=r'primary tracks $\chi^2_{p}$',
    ylabel= pfps_y_label,
    cut_value = [CTE.MIP_candidate_min_proton_score],
    clip = False
)

config_MIP_candidate_chi2_p = FullHistogramConfig(
    file_name = "MIP_candidate_chi2_p", 
    var_evt_reco_col=('pfp','trk','chi2pid','best','chi2_proton',''),
    truth_column=pfp_truth_column,
    first_per_slice = False,
    start_cut = "track",
    end_cut = "energy",
    extra_mask = "MIP_candidate",
    bins=np.linspace(0, 300, 41),
    xlabel=r'MIP candidates $\chi^2_{p}$',
    ylabel= pfps_y_label,
    cut_value = [CTE.MIP_candidate_min_proton_score],
    clip = False
)

config_contained_MIP_candidate_chi2_p = FullHistogramConfig(
    file_name = "contained_MIP_candidate_chi2_p", 
    var_evt_reco_col=('pfp','trk','chi2pid','best','chi2_proton',''),
    truth_column=pfp_truth_column,
    first_per_slice = False,
    start_cut = "track",
    end_cut = "energy",
    extra_mask = "contained_MIP_candidate",
    bins=np.linspace(0, 300, 41),
    xlabel=r'contained MIP candidates $\chi^2_{p}$',
    ylabel= pfps_y_label,
    cut_value = [CTE.MIP_candidate_min_proton_score],
    clip = False
)

config_exiting_MIP_candidate_chi2_p = FullHistogramConfig(
    file_name = "exiting_MIP_candidate_chi2_p", 
    var_evt_reco_col=('pfp','trk','chi2pid','best','chi2_proton',''),
    truth_column=pfp_truth_column,
    first_per_slice = False,
    start_cut = "track",
    end_cut = "energy",
    extra_mask = "exiting_MIP_candidate",
    bins=np.linspace(0, 300, 41),
    xlabel=r'exiting MIP candidates $\chi^2_{p}$',
    ylabel= pfps_y_label,
    cut_value = [CTE.MIP_candidate_min_proton_score],
    clip = False
)


#shower energy
config_primary_shower_energy = FullHistogramConfig(
    file_name = "primary_shower_energy", 
    var_evt_reco_col=('pfp','shw','bestplane_energy','','',''),
    truth_column=pfp_truth_column,
    first_per_slice = False,
    start_cut = "chi2",
    end_cut = "shower",
    extra_mask = "primary_shower",
    bins=np.linspace(0, 0.5, 41),
    xlabel='primary shower energy [GeV]',
    ylabel= pfps_y_label,
    cut_value = [CTE.min_shower_ke],
    clip = False
)



#scatter angle
config_MIP_candidate_scatter_angle = FullHistogramConfig(
    file_name = "MIP_candidate_scatter_angle_ratio", 
    var_evt_reco_col=('pfp','scatter_angle_ratio','','','',''),
    truth_column=pfp_truth_column,
    first_per_slice = False,
    start_cut = "chi2",
    end_cut = "energy",
    extra_mask = "MIP_candidate",
    bins=np.linspace(0, 1, 41),
    xlabel=r'MIP candidates $\frac{MCS max scatter}{MCS total scatter}$',
    ylabel= pfps_y_label
)

config_contained_MIP_candidate_scatter_angle = FullHistogramConfig(
    file_name = "contained_MIP_candidate_scatter_angle", 
    var_evt_reco_col=('pfp','scatter_angle_ratio','','','',''),
    truth_column=pfp_truth_column,
    first_per_slice = False,
    start_cut = "chi2",
    end_cut = "energy",
    extra_mask = "contained_MIP_candidate",
    bins=np.linspace(0, 1, 41),
    xlabel=r'Contained MIP candidates $\frac{MCS max scatter}{MCS total scatter}$',
    ylabel= pfps_y_label
)

config_exiting_MIP_candidate_scatter_angle = FullHistogramConfig(
    file_name = "exiting_MIP_candidate_scatter_angle", 
    var_evt_reco_col=('pfp','scatter_angle_ratio','','','',''),
    truth_column=pfp_truth_column,
    first_per_slice = False,
    start_cut = "chi2",
    end_cut = "energy",
    extra_mask = "exiting_MIP_candidate",
    bins=np.linspace(0, 1, 41),
    xlabel=r'Exiting MIP candidates $\frac{MCS max scatter}{MCS total scatter}$',
    ylabel= pfps_y_label
)


#daughter hits
config_MIP_candidate_max_daughter_hits = FullHistogramConfig(
    file_name = "MIP_candidate_max_daughter_hits", 
    var_evt_reco_col=('pfp','max_daughter_hits','','','',''),
    truth_column=pfp_truth_column,
    first_per_slice = False,
    start_cut = "chi2",
    end_cut = "energy",
    extra_mask = "MIP_candidate",
    bins=np.linspace(10, 300, 41),
    xlabel=r'MIP candidates daughter max hits$',
    ylabel= pfps_y_label
)

config_contained_MIP_candidate_max_daughter_hits = FullHistogramConfig(
    file_name = "contained_MIP_candidate_max_daughter_hits", 
    var_evt_reco_col=('pfp','max_daughter_hits','','','',''),
    truth_column=pfp_truth_column,
    first_per_slice = False,
    start_cut = "chi2",
    end_cut = "energy",
    extra_mask = "contained_MIP_candidate",
    bins=np.linspace(10, 300, 41),
    xlabel=r'Contained MIP candidates daughter max hits',
    ylabel= pfps_y_label
)

config_exiting_MIP_candidate_max_daughter_hits = FullHistogramConfig(
    file_name = "exiting_MIP_candidate_max_daughter_hits", 
    var_evt_reco_col=('pfp','max_daughter_hits','','','',''),
    truth_column=pfp_truth_column,
    first_per_slice = False,
    start_cut = "chi2",
    end_cut = "energy",
    extra_mask = "exiting_MIP_candidate",
    bins=np.linspace(10, 300, 41),
    xlabel=r'Exiting MIP candidates daughter max hits',
    ylabel= pfps_y_label
)


#frac 50
config_MIP_candidate_frac50 = FullHistogramConfig(
    file_name = "MIP_candidate_scatter_frac50", 
    var_evt_reco_col=('pfp','trk','frac50','','',''),
    truth_column=pfp_truth_column,
    first_per_slice = False,
    start_cut = "chi2",
    end_cut = "energy",
    extra_mask = "MIP_candidate",
    bins=np.linspace(0, 1, 41),
    xlabel=r'MIP candidates RR frac. with 50% energy',
    ylabel= pfps_y_label
)

config_contained_MIP_candidate_frac50 = FullHistogramConfig(
    file_name = "contained_MIP_candidate_frac50", 
    var_evt_reco_col=('pfp','trk','frac50','','',''),
    truth_column=pfp_truth_column,
    first_per_slice = False,
    start_cut = "chi2",
    end_cut = "energy",
    extra_mask = "contained_MIP_candidate",
    bins=np.linspace(0, 1, 41),
    xlabel=r'Contained MIP candidates RR frac. with 50% energy',
    ylabel= pfps_y_label
)

config_exiting_MIP_candidate_frac50 = FullHistogramConfig(
    file_name = "exiting_MIP_candidate_frac50", 
    var_evt_reco_col=('pfp','trk','frac50','','',''),
    truth_column=pfp_truth_column,
    first_per_slice = False,
    start_cut = "chi2",
    end_cut = "energy",
    extra_mask = "exiting_MIP_candidate",
    bins=np.linspace(0, 1, 41),
    xlabel=r'Exiting MIP candidates RR frac. with 50% energy',
    ylabel= pfps_y_label
)

#expol0
config_MIP_candidate_chi2_exp_pol = FullHistogramConfig(
    file_name = "MIP_candidate_scatter_chi2_exp_pol", 
    var_evt_reco_col=('pfp','trk','chi2_exp_pol','','',''),
    truth_column=pfp_truth_column,
    first_per_slice = False,
    start_cut = "chi2",
    end_cut = "energy",
    extra_mask = "MIP_candidate",
    bins=np.linspace(0, 1, 41),
    xlabel=r'MIP candidates $\chi^2_{pol_0} / \chi^2_{exp}$',
    ylabel= pfps_y_label
)

config_contained_MIP_candidate_chi2_exp_pol = FullHistogramConfig(
    file_name = "contained_MIP_candidate_chi2_exp_pol", 
    var_evt_reco_col=('pfp','trk','chi2_exp_pol','','',''),
    truth_column=pfp_truth_column,
    first_per_slice = False,
    start_cut = "chi2",
    end_cut = "energy",
    extra_mask = "contained_MIP_candidate",
    bins=np.linspace(0, 1, 41),
    xlabel=r'Contained MIP candidates $\chi^2_{pol_0} / \chi^2_{exp}$',
    ylabel= pfps_y_label
)

config_exiting_MIP_candidate_chi2_exp_pol = FullHistogramConfig(
    file_name = "exiting_MIP_candidate_chi2_exp_pol", 
    var_evt_reco_col=('pfp','trk','chi2_exp_pol','','',''),
    truth_column=pfp_truth_column,
    first_per_slice = False,
    start_cut = "chi2",
    end_cut = "energy",
    extra_mask = "exiting_MIP_candidate",
    bins=np.linspace(0, 1, 41),
    xlabel=r'Exiting MIP candidates $\chi^2_{pol_0} / \chi^2_{exp}$',
    ylabel= pfps_y_label
)

#bdt_score_proton
config_MIP_candidate_bdt_score_proton = FullHistogramConfig(
    file_name = "MIP_candidate_bdt_score_proton", 
    var_evt_reco_col=('pfp','trk','bdt_proton_score','','',''),
    truth_column=pfp_truth_column,
    first_per_slice = False,
    start_cut = "chi2",
    end_cut = "energy",
    extra_mask = "MIP_candidate",
    bins=np.linspace(-5, 9, 41),
    xlabel=r'MIP candidates proton BDT score',
    cut_value = [CTE.BDT_proton_max_score],
    ylabel= pfps_y_label
)

config_contained_MIP_candidate_bdt_score_proton = FullHistogramConfig(
    file_name = "contained_MIP_candidate_bdt_score_proton", 
    var_evt_reco_col=('pfp','trk','bdt_proton_score','','',''),
    truth_column=pfp_truth_column,
    first_per_slice = False,
    start_cut = "chi2",
    end_cut = "energy",
    extra_mask = "contained_MIP_candidate",
    bins=np.linspace(-5, 9, 41),
    xlabel=r'Contained MIP candidates proton BDT score',
    cut_value = [CTE.BDT_proton_max_score],
    ylabel= pfps_y_label
)

config_exiting_MIP_candidate_bdt_score_proton = FullHistogramConfig(
    file_name = "exiting_MIP_candidate_bdt_score_proton", 
    var_evt_reco_col=('pfp','trk','bdt_proton_score','','',''),
    truth_column=pfp_truth_column,
    first_per_slice = False,
    start_cut = "chi2",
    end_cut = "energy",
    extra_mask = "exiting_MIP_candidate",
    bins=np.linspace(-5, 9, 41),
    xlabel=r'Exiting MIP candidates proton BDT score',
    cut_value = [CTE.BDT_proton_max_score],
    ylabel= pfps_y_label
)

#bdt_score_muon_pion
config_MIP_candidate_bdt_score_muon_pion = FullHistogramConfig(
    file_name = "MIP_candidate_bdt_score_muon_pion", 
    var_evt_reco_col=('pfp','trk','bdt_muon_pion_score','','',''),
    truth_column=pfp_truth_column,
    first_per_slice = False,
    start_cut = "chi2",
    end_cut = "energy",
    extra_mask = "MIP_candidate",
    bins=np.linspace(-5, 7, 41),
    xlabel=r'MIP candidates $\mu/\pi$ separation BDT score',
    ylabel= pfps_y_label
)

config_contained_MIP_candidate_bdt_score_muon_pion = FullHistogramConfig(
    file_name = "contained_MIP_candidate_bdt_score_muon_pion", 
    var_evt_reco_col=('pfp','trk','bdt_muon_pion_score','','',''),
    truth_column=pfp_truth_column,
    first_per_slice = False,
    start_cut = "chi2",
    end_cut = "energy",
    extra_mask = "contained_MIP_candidate",
    bins=np.linspace(-5, 7, 41),
    xlabel=r'Contained MIP candidates $\mu/\pi$ separation BDT score',
    ylabel= pfps_y_label
)

config_exiting_MIP_candidate_bdt_score_muon_pion = FullHistogramConfig(
    file_name = "exiting_MIP_candidate_bdt_score_muon_pion", 
    var_evt_reco_col=('pfp','trk','bdt_muon_pion_score','','',''),
    truth_column=pfp_truth_column,
    first_per_slice = False,
    start_cut = "chi2",
    end_cut = "energy",
    extra_mask = "exiting_MIP_candidate",
    bins=np.linspace(-5, 7, 41),
    xlabel=r'Exiting MIP candidates $\mu/\pi$ separation BDT score',
    ylabel= pfps_y_label
)

particle_vars_vec = [
    config_pandora_primary_vtx_distance,
    config_pandora_primary_track_score,
    config_analysis_primary_track_score ,
    config_contained_MIP_candidate_track_score,
    config_contained_MIP_candidate_KE,
    config_analysis_primary_len,
    config_primary_track_len,
    config_primary_track_chi2_mu,
    config_MIP_candidate_chi2_mu,
    config_contained_MIP_candidate_chi2_mu,
    config_exiting_MIP_candidate_chi2_mu,
    config_primary_track_chi2_p,
    config_MIP_candidate_chi2_p,
    config_contained_MIP_candidate_chi2_p,
    config_exiting_MIP_candidate_chi2_p,
    config_primary_shower_energy,
    config_MIP_candidate_scatter_angle,
    config_contained_MIP_candidate_scatter_angle,
    config_exiting_MIP_candidate_scatter_angle,
    config_MIP_candidate_max_daughter_hits,
    config_contained_MIP_candidate_max_daughter_hits,
    config_exiting_MIP_candidate_max_daughter_hits,
    config_MIP_candidate_frac50,
    config_contained_MIP_candidate_frac50,
    config_exiting_MIP_candidate_frac50,
    config_MIP_candidate_chi2_exp_pol,
    config_contained_MIP_candidate_chi2_exp_pol,
    config_exiting_MIP_candidate_chi2_exp_pol,
    config_MIP_candidate_bdt_score_proton,
    config_contained_MIP_candidate_bdt_score_proton,
    config_exiting_MIP_candidate_bdt_score_proton,
    config_MIP_candidate_bdt_score_muon_pion,
    config_contained_MIP_candidate_bdt_score_muon_pion,
    config_exiting_MIP_candidate_bdt_score_muon_pion

]



#####################
config_angle_between_candidates = FullHistogramConfig(
    file_name = "angle_between_candidates",      
    var_evt_reco_col=('slc','measure_var','angle_between_candidates','','',''),
    truth_column = nu_categ_column,
    first_per_slice = True,
    start_cut = "chi2",
    end_cut = "energy",
    bins=np.linspace(0, np.pi, 41),
    xlabel=r'Angle between MIP candidates [rad]',
    cut_value = [CTE.max_angle_between_candidates],
    ylabel=slices_y_label
)

config_num_protons = FullHistogramConfig(
    file_name = "num_protons",      
    var_evt_reco_col=('slc','measure_var','num_protons','','',''),
    truth_column = nu_categ_column,
    first_per_slice = True,
    start_cut = "energy",
    end_cut = "energy",
    bins=np.linspace(0, 3, 4),
    xlabel=r'Number of protons',
    ylabel=slices_y_label
)

pion_p_type_column = ('slc','measure_var','pi_true_p_type','','','')
muon_p_type_column = ('slc','measure_var','mu_true_p_type','','','')

#Pion Measure vars
config_p_pi = FullHistogramConfig(
    file_name = "TLE_p_pi",      
    var_evt_reco_col=('slc','measure_var','TLE_p_pi','','',''),
    truth_column = pion_p_type_column,
    first_per_slice = True,
    start_cut = "energy",
    end_cut = "energy",
    bins=np.linspace(0.13, 0.8, 41),
    xlabel=r'Pion candidate P [GeV]',
    ylabel=slices_y_label
)

config_chi2_proton_pi = FullHistogramConfig(
    file_name = "chi2_proton_pi",      
    var_evt_reco_col=('slc','measure_var','pi_chi2_proton','','',''),
    truth_column = pion_p_type_column,
    first_per_slice = True,
    start_cut = "energy",
    end_cut = "energy",
    bins=np.linspace(0, 300, 41),
    xlabel=r'Pion candidate $\chi^2_{p}$',
    ylabel=slices_y_label
)

config_chi2_mu_pi = FullHistogramConfig(
    file_name = "chi2_mu_pi",      
    var_evt_reco_col=('slc','measure_var','pi_chi2_mu','','',''),
    truth_column = pion_p_type_column,
    first_per_slice = True,
    start_cut = "energy",
    end_cut = "energy",
    bins=np.linspace(0, 20, 41),
    xlabel=r'Pion candidate $\chi^2_{\mu}$',
    ylabel=slices_y_label
)

config_cos_theta_pi = FullHistogramConfig(
    file_name = "cos_theta_pi",      
    var_evt_reco_col=('slc','measure_var','reco_cos_theta_pi','','',''),
    truth_column = pion_p_type_column,
    first_per_slice = True,
    start_cut = "energy",
    end_cut = "energy",
    bins=np.linspace(-1, 1, 41),
    xlabel=r'Pion candidate $cos_{\theta_z}$',
    ylabel=slices_y_label,
    stats_horizontal_alignment = 'left'
)

config_pi_chi2_exp_pol = FullHistogramConfig(
    file_name = "pi_chi2_exp_pol",      
    var_evt_reco_col=('slc','measure_var','pi_chi2_exp_pol','','',''),
    truth_column = pion_p_type_column,
    first_per_slice = True,
    start_cut = "energy",
    end_cut = "energy",
    bins=np.linspace(0, 1, 41),
    xlabel=r'Pion candidate $\chi^2_{pol_0} / \chi^2_{exp}$',
    ylabel=slices_y_label,
    stats_horizontal_alignment = 'left'
)

config_pi_scatter_angle_ratio = FullHistogramConfig(
    file_name = "pi_scatter_angle_ratio",      
    var_evt_reco_col=('slc','measure_var','pi_scatter_angle_ratio','','',''),
    truth_column = pion_p_type_column,
    first_per_slice = True,
    start_cut = "energy",
    end_cut = "energy",
    bins=np.linspace(0, 1, 41),
    xlabel=r'Pion candidate $\frac{MCS max scatter}{MCS total scatter}$',
    ylabel=slices_y_label
)

config_pi_max_daughter_hits = FullHistogramConfig(
    file_name = "pi_max_daughter_hits",      
    var_evt_reco_col=('slc','measure_var','pi_max_daughter_hits','','',''),
    truth_column = pion_p_type_column,
    first_per_slice = True,
    start_cut = "energy",
    end_cut = "energy",
    bins=np.linspace(1, 300, 41),
    xlabel=r'Pion candidate daughter max hits',
    ylabel=slices_y_label
)

config_pi_frac50 = FullHistogramConfig(
    file_name = "pi_frac50",      
    var_evt_reco_col=('slc','measure_var','pi_frac50','','',''),
    truth_column = pion_p_type_column,
    first_per_slice = True,
    start_cut = "energy",
    end_cut = "energy",
    bins=np.linspace(0, 1, 41),
    xlabel=r'Pion candidate RR frac. with 50% energy',
    ylabel=slices_y_label
)

config_pi_bdt_score_proton = FullHistogramConfig(
    file_name = "pi_bdt_score_proton",      
    var_evt_reco_col=('slc','measure_var','pi_bdt_score_proton','','',''),
    truth_column = pion_p_type_column,
    first_per_slice = True,
    start_cut = "energy",
    end_cut = "energy",
    bins=np.linspace(-5, 9, 41),
    xlabel=r'Pion candidate proton BDT score',
    ylabel=slices_y_label,
    stats_horizontal_alignment = 'left'
)

config_pi_bdt_score_muon_pion = FullHistogramConfig(
    file_name = "pi_bdt_score_muon_pion",      
    var_evt_reco_col=('slc','measure_var','pi_bdt_score_muon_pion','','',''),
    truth_column = pion_p_type_column,
    first_per_slice = True,
    start_cut = "energy",
    end_cut = "energy",
    bins=np.linspace(-5, 7, 41),
    xlabel=r'Pion candidate $\mu/\pi$ separation BDT score',
    ylabel=slices_y_label,
)




#Muon Measure vars
config_p_mu = FullHistogramConfig(
    file_name = "reco_p_mu",      
    var_evt_reco_col=('slc','measure_var','reco_p_mu','','',''),
    truth_column = muon_p_type_column,
    first_per_slice = True,
    start_cut = "energy",
    end_cut = "energy",
    bins=np.linspace(0.1,3, 41),
    xlabel=r'Muon candidate P [GeV]',
    ylabel=slices_y_label
)

config_p_mu_contained = FullHistogramConfig(
    file_name = "reco_p_mu_contained",      
    var_evt_reco_col=('slc','measure_var','reco_p_mu','','',''),
    truth_column = muon_p_type_column,
    first_per_slice = True,
    start_cut = "energy",
    end_cut = "energy",
    extra_mask = "muon_contained",
    bins=np.linspace(0.1,3, 41),
    xlabel=r'Contained muon candidate P [GeV]',
    ylabel=slices_y_label
)

config_p_mu_exiting = FullHistogramConfig(
    file_name = "reco_p_mu_exiting",      
    var_evt_reco_col=('slc','measure_var','reco_p_mu','','',''),
    truth_column = muon_p_type_column,
    first_per_slice = True,
    start_cut = "energy",
    end_cut = "energy",
    extra_mask = "muon_exiting",
    bins=np.linspace(0.1,3, 41),
    xlabel=r'Exiting muon candidate P [GeV]',
    ylabel=slices_y_label
)

config_chi2_proton_mu = FullHistogramConfig(
    file_name = "chi2_proton_mu",      
    var_evt_reco_col=('slc','measure_var','mu_chi2_proton','','',''),
    truth_column = muon_p_type_column,
    first_per_slice = True,
    start_cut = "energy",
    end_cut = "energy",
    bins=np.linspace(0, 300, 41),
    xlabel=r'Muon candidate $\chi^2_{p}$',
    ylabel=slices_y_label
)

config_chi2_mu_mu = FullHistogramConfig(
    file_name = "chi2_mu_mu",      
    var_evt_reco_col=('slc','measure_var','mu_chi2_mu','','',''),
    truth_column = muon_p_type_column,
    first_per_slice = True,
    start_cut = "energy",
    end_cut = "energy",
    bins=np.linspace(0, 20, 41),
    xlabel=r'Muon candidate $\chi^2_{\mu}$',
    ylabel=slices_y_label
)

config_cos_theta_mu = FullHistogramConfig(
    file_name = "cos_theta_mu",      
    var_evt_reco_col=('slc','measure_var','reco_cos_theta_mu','','',''),
    truth_column = muon_p_type_column,
    first_per_slice = True,
    start_cut = "energy",
    end_cut = "energy",
    bins=np.linspace(-1, 1, 41),
    xlabel=r'Muon candidate $cos_{\theta_z}$',
    ylabel=slices_y_label,
    stats_horizontal_alignment = 'left'
)

config_mu_chi2_exp_pol = FullHistogramConfig(
    file_name = "mu_chi2_exp_pol",      
    var_evt_reco_col=('slc','measure_var','mu_chi2_exp_pol','','',''),
    truth_column = muon_p_type_column,
    first_per_slice = True,
    start_cut = "energy",
    end_cut = "energy",
    bins=np.linspace(0, 1, 41),
    xlabel=r'Muon candidate $\chi^2_{pol_0} / \chi^2_{exp}$',
    ylabel=slices_y_label,
    stats_horizontal_alignment =  'left',
)

config_mu_scatter_angle_ratio = FullHistogramConfig(
    file_name = "mu_scatter_angle_ratio",      
    var_evt_reco_col=('slc','measure_var','mu_scatter_angle_ratio','','',''),
    truth_column = muon_p_type_column,
    first_per_slice = True,
    start_cut = "energy",
    end_cut = "energy",
    bins=np.linspace(0, 1, 41),
    xlabel=r'Muon candidate $\frac{MCS max scatter}{MCS total scatter}$',
    ylabel=slices_y_label
)

config_mu_max_daughter_hits = FullHistogramConfig(
    file_name = "mu_max_daughter_hits",      
    var_evt_reco_col=('slc','measure_var','mu_max_daughter_hits','','',''),
    truth_column = muon_p_type_column,
    first_per_slice = True,
    start_cut = "energy",
    end_cut = "energy",
    bins=np.linspace(1, 300, 41),
    xlabel=r'Muon candidate daughter max hits',
    ylabel=slices_y_label
)

config_mu_frac50 = FullHistogramConfig(
    file_name = "mu_frac50",      
    var_evt_reco_col=('slc','measure_var','mu_frac50','','',''),
    truth_column = muon_p_type_column,
    first_per_slice = True,
    start_cut = "energy",
    end_cut = "energy",
    bins=np.linspace(0, 1, 41),
    xlabel=r'Muon candidate RR frac. with 50% energy',
    ylabel=slices_y_label
)

config_mu_bdt_score_proton = FullHistogramConfig(
    file_name = "mu_bdt_score_proton",      
    var_evt_reco_col=('slc','measure_var','mu_bdt_score_proton','','',''),
    truth_column = muon_p_type_column,
    first_per_slice = True,
    start_cut = "energy",
    end_cut = "energy",
    bins=np.linspace(-5, 9, 41),
    xlabel=r'Muon candidate proton BDT score',
    ylabel=slices_y_label,
    stats_horizontal_alignment = 'left'
)

config_mu_bdt_score_muon_pion = FullHistogramConfig(
    file_name = "mu_bdt_score_muon_pion",      
    var_evt_reco_col=('slc','measure_var','mu_bdt_score_muon_pion','','',''),
    truth_column = muon_p_type_column,
    first_per_slice = True,
    start_cut = "energy",
    end_cut = "energy",
    bins=np.linspace(-5, 7, 41),
    xlabel=r'Muon candidate $\mu/\pi$ separation BDT score',
    ylabel=slices_y_label
)



measure_var_vec = [
    config_p_pi,
    config_chi2_proton_pi,
    config_chi2_mu_pi,
    config_cos_theta_pi,
    config_pi_chi2_exp_pol,
    config_pi_scatter_angle_ratio,
    config_pi_max_daughter_hits,
    config_pi_frac50,
    config_pi_bdt_score_proton,
    config_pi_bdt_score_muon_pion,
    config_p_mu,
    config_p_mu_contained,
    config_p_mu_exiting,
    config_chi2_proton_mu,
    config_chi2_mu_mu,
    config_cos_theta_mu,
    config_mu_chi2_exp_pol,
    config_mu_scatter_angle_ratio,
    config_mu_max_daughter_hits,
    config_mu_frac50,
    config_mu_bdt_score_proton,
    config_mu_bdt_score_muon_pion,
    
    config_angle_between_candidates,
    config_num_protons
]



###########

config_delta_alpha_T = FullHistogramConfig(
    file_name = "delta_alpha_T",      
    var_evt_reco_col=('slc','measure_var','delta_alpha_T','','',''),
    truth_column = nu_categ_column,
    first_per_slice = True,
    start_cut = "energy",
    end_cut = "energy",
    bins=np.linspace(0, np.pi, 41),
    xlabel=r'$\delta \alpha_T$ [rad]',
    ylabel=slices_y_label,
)


config_delta_pT = FullHistogramConfig(
    file_name = "delta_pT",      
    var_evt_reco_col=('slc','measure_var','delta_pT','','',''),
    truth_column = nu_categ_column,
    first_per_slice = True,
    start_cut = "energy",
    end_cut = "energy",
    bins=np.linspace(0, 1, 41),
    xlabel=r'$\delta p_T$ [GeV]',
    ylabel=slices_y_label
)

config_delta_phi_T = FullHistogramConfig(
    file_name = "delta_phi_T",      
    var_evt_reco_col=('slc','measure_var','delta_phi_T','','',''),
    truth_column = nu_categ_column,
    first_per_slice = True,
    start_cut = "energy",
    end_cut = "energy",
    bins=np.linspace(0, np.pi, 41),
    xlabel=r'$\delta \phi_T$ [rad]',
    ylabel=slices_y_label
)

TKI_vec = [config_delta_alpha_T, config_delta_pT, config_delta_phi_T]


config_delta_alpha_T_genie_categ = FullHistogramConfig(
    file_name = "delta_alpha_T_genie_categ",      
    var_evt_reco_col=('slc','measure_var','delta_alpha_T','','',''),
    truth_column = ('truth', 'genie_categ', '', '','',''),
    first_per_slice = True,
    start_cut = "energy",
    end_cut = "energy",
    bins=np.linspace(0, np.pi, 41),
    xlabel=r'$\delta \alpha_T$ [rad]',
    ylabel=slices_y_label,
)

config_delta_pT_genie_categ = FullHistogramConfig(
    file_name = "delta_pT_genie_categ",      
    var_evt_reco_col=('slc','measure_var','delta_pT','','',''),
    truth_column = ('truth', 'genie_categ', '', '','',''),
    first_per_slice = True,
    start_cut = "energy",
    end_cut = "energy",
    bins=np.linspace(0, 1, 41),
    xlabel=r'$\delta p_T$ [GeV]',
    ylabel=slices_y_label
)

config_delta_phi_T_genie_categ = FullHistogramConfig(
    file_name = "delta_phi_T_genie_categ",      
    var_evt_reco_col=('slc','measure_var','delta_phi_T','','',''),
    truth_column = ('truth', 'genie_categ', '', '','',''),
    first_per_slice = True,
    start_cut = "energy",
    end_cut = "energy",
    bins=np.linspace(0, np.pi, 41),
    xlabel=r'$\delta \phi_T$ [rad]',
    ylabel=slices_y_label
)


TKI_genie_categ_vec= [config_delta_alpha_T_genie_categ, config_delta_pT_genie_categ, config_delta_phi_T_genie_categ]






config_all_evts_final = FullHistogramConfig(
    file_name = "all_evts",      
    var_evt_reco_col=('slc','measure_var','reco_p_mu','','',''),
    truth_column = nu_categ_column,
    first_per_slice = True,
    start_cut = "energy",
    end_cut = "energy",
    bins=np.linspace(0., 1.,2),
    xlabel=r'All Events',
    ylabel=slices_y_label
)


config_p_mu_final = FullHistogramConfig(
    file_name = "muon_p",      
    var_evt_reco_col=('slc','measure_var','reco_p_mu','','',''),
    truth_column = nu_categ_column,
    first_per_slice = True,
    start_cut = "energy",
    end_cut = "energy",
    bins= np.array([0.1, 0.242, 0.335,0.467, 1]),
    xlabel=r'Muon candidate P [GeV]',
    ylabel=slices_y_label
)

config_cos_theta_mu_final = FullHistogramConfig(
    file_name = "muon_dir_z",      
    var_evt_reco_col=('slc','measure_var','reco_cos_theta_mu','','',''),
    truth_column = nu_categ_column,
    first_per_slice = True,
    start_cut = "energy",
    end_cut = "energy",
    bins= np.array([-1.,0.324, 0.667, 0.855, 1.]),
    xlabel=r'Muon candidate $cos_{\theta_z}$',
    ylabel=slices_y_label,
    stats_horizontal_alignment = 'left'
)

config_p_pi_final = FullHistogramConfig(
    file_name = "pion_p",      
    var_evt_reco_col=('slc','measure_var','TLE_p_pi','','',''),
    truth_column = nu_categ_column,
    first_per_slice = True,
    start_cut = "energy",
    end_cut = "energy",
    bins= np.array([0.13, 0.218, 0.296,0.415,1]),
    xlabel=r'Pion candidate P [GeV]',
    ylabel=slices_y_label
)

config_cos_theta_pi_final = FullHistogramConfig(
    file_name = "pion_dir_z",      
    var_evt_reco_col=('slc','measure_var','reco_cos_theta_pi','','',''),
    truth_column = nu_categ_column,
    first_per_slice = True,
    start_cut = "energy",
    end_cut = "energy",
    bins= np.array([-1.,0.131, 0.592, 0.830, 1.]),
    xlabel=r'Pion candidate $cos_{\theta_z}$',
    ylabel=slices_y_label,
    stats_horizontal_alignment = 'left'
)


config_delta_alpha_T_final = FullHistogramConfig(
    file_name = "delta_alpha_T",      
    var_evt_reco_col=('slc','measure_var','delta_alpha_T','','',''),
    truth_column = nu_categ_column,
    first_per_slice = True,
    start_cut = "energy",
    end_cut = "energy",
    bins= np.array([0., 1.637, 2.643, np.pi]),
    xlabel=r'$\delta \alpha_T$ [rad]',
    ylabel=slices_y_label,
)


config_delta_pT_final = FullHistogramConfig(
    file_name = "delta_pt",      
    var_evt_reco_col=('slc','measure_var','delta_pT','','',''),
    truth_column = nu_categ_column,
    first_per_slice = True,
    start_cut = "energy",
    end_cut = "energy",
    bins= np.array([0., 0.158, 0.372, 0.8]),
    xlabel=r'$\delta p_T$ [GeV]',
    ylabel=slices_y_label
)

config_delta_phi_T_final = FullHistogramConfig(
    file_name = "delta_phi_T",      
    var_evt_reco_col=('slc','measure_var','delta_phi_T','','',''),
    truth_column = nu_categ_column,
    first_per_slice = True,
    start_cut = "energy",
    end_cut = "energy",
    bins= np.array([0., 0.542, 1.96, np.pi]),
    xlabel=r'$\delta \phi_T$ [rad]',
    ylabel=slices_y_label
)


config_num_protons = FullHistogramConfig(
    file_name = "num_protons",      
    var_evt_reco_col=('slc','measure_var','num_protons','','',''),
    truth_column = nu_categ_column,
    first_per_slice = True,
    start_cut = "energy",
    end_cut = "energy",
    bins=np.linspace(0, 3, 4),
    xlabel=r'Number of protons',
    ylabel=slices_y_label
)

config_angle_between_candidates_final = FullHistogramConfig(
    file_name = "angle_between_candidates",      
    var_evt_reco_col=('slc','measure_var','angle_between_candidates','','',''),
    truth_column = nu_categ_column,
    first_per_slice = True,
    start_cut = "energy",
    end_cut = "energy",
    bins= np.array([0., 0.954,1.392, 1.848, 2.65]),
    xlabel=r'Angle between MIP candidates [rad]',
    ylabel=slices_y_label
)

final_var_configs = [config_all_evts_final, config_p_mu_final, config_cos_theta_mu_final, config_p_pi_final,
                     config_cos_theta_pi_final, config_delta_alpha_T_final, config_delta_pT_final,
                    config_delta_phi_T_final, config_num_protons, config_angle_between_candidates_final]
