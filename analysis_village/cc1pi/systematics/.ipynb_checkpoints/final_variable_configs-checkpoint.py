# variable to unfold
import numpy as np
import inspect

class VariableConfig:
    """
    A configurable class for setting up unfolding variable configurations.
    Choose a configuration using one of the provided class methods,
    or instantiate directly with custom parameters.
    """
    def __init__(self, var_save_name, var_plot_name, var_unit, bins, var_evt_reco_col, var_evt_truth_col, var_nu_col, xsec_label):
        self.var_save_name = var_save_name
        self.var_plot_name = var_plot_name
        self.var_unit = var_unit
        unit_str = rf" {var_unit}" if var_unit else ""
        self.var_labels = [
            rf"{var_plot_name}{unit_str}", 
            rf"{var_plot_name}$^{{\text{{reco}}}}${unit_str}", 
            rf"{var_plot_name}$^{{\text{{true}}}}${unit_str}"
        ]
        self.bins = bins
        self.bin_centers = (bins[:-1] + bins[1:]) / 2.
        self.var_evt_reco_col = var_evt_reco_col
        self.var_evt_truth_col = var_evt_truth_col
        self.var_nu_col = var_nu_col
        self.xsec_label = xsec_label


    @classmethod
    def all_evts(cls):
        return cls(
            var_save_name="all_evts",
            var_plot_name=r"All Events",
            var_unit="",
            bins=np.linspace(0., 1.,2),
            var_evt_reco_col=('slc','measure_var', 'reco_p_mu', ''),
            var_evt_truth_col=('truth', 'true_var', 'true_p_mu', ''),
            var_nu_col=('truth', 'true_var', 'true_p_mu', ''),
            xsec_label=""        
        )
        
    @classmethod
    def muon_momentum(cls):
        return cls(
            var_save_name="muon_p",
            var_plot_name=r"$P_\mu$",
            var_unit="[GeV/c]",
            bins= np.array([0.1, 0.242, 0.335,0.467, 1]),
            var_evt_reco_col=('slc','measure_var', 'reco_p_mu', ''),
            var_evt_truth_col=('truth', 'true_var', 'true_p_mu', ''),
            var_nu_col=('truth', 'true_var', 'true_p_mu', ''),
            xsec_label=r"$\frac{d\sigma}{dP_\mu}$ ($\mathrm{cm^2}$ [$\mathrm{GeV/c}$])"        
        )

    @classmethod
    def muon_direction(cls):
        return cls(
            var_save_name="muon_dir_z",
            var_plot_name=r"$cos(\theta_\mu)$",
            var_unit="",
            bins= np.array([-1.,0.324, 0.667, 0.855, 1.]),
            var_evt_reco_col=('slc','measure_var', 'reco_cos_theta_mu', ''),
            var_evt_truth_col=('truth', 'true_var', 'true_cos_theta_mu', ''),
            var_nu_col=('truth', 'true_var', 'true_cos_theta_mu', ''),
            xsec_label=r"$\frac{d\sigma}{d\cos(\theta_\mu)}$ ($\mathrm{cm}^2$)"
        )

    @classmethod
    def pion_momentum(cls):
        return cls(
            var_save_name="pion_p",
            var_plot_name=r"$P_\pi$",
            var_unit="[GeV/c]",
            bins= np.array([0.13, 0.218, 0.296,0.415,2]),
            var_evt_reco_col=('slc', 'measure_var', 'TLE_p_pi', ''),
            var_evt_truth_col=('truth', 'true_var', 'true_p_pi', ''),
            var_nu_col=('truth', 'true_var', 'true_p_pi', ''),
            xsec_label=r"$\frac{d\sigma}{dP_\pi}$ ($\mathrm{cm}^2$ / GeV/c)"
        )

    @classmethod
    def pion_direction(cls):
        return cls(
            var_save_name="pion_dir_z",
            var_plot_name=r"$cos(\theta_\pi)$",
            var_unit="",            
            bins= np.array([-1.,0.131, 0.592, 0.830, 1.]),
            var_evt_reco_col=('slc','measure_var', 'reco_cos_theta_pi', ''),
            var_evt_truth_col=('truth', 'true_var', 'true_cos_theta_pi', ''),
            var_nu_col=('truth', 'true_var', 'true_cos_theta_pi', ''),
            xsec_label=r"$\frac{d\sigma}{dcos(\theta_\pi)}$ ($\mathrm{cm}^2$)"
        )

    @classmethod
    def angle_between_candidates(cls):
        return cls(
            var_save_name="angle_between_candidates",
            var_plot_name=r"$\theta_{\pi,\mu}$",
            var_unit="[rad]", 
            bins= np.array([0., 0.954,1.392, 1.848, 2.65]),
            var_evt_reco_col=('slc', 'measure_var', 'angle_between_candidates', ''),
            var_evt_truth_col=('truth', 'true_var', 'true_mu_pi_angle', ''),
            var_nu_col=('truth', 'true_var', 'true_mu_pi_angle', ''),
            xsec_label=r"$\frac{d\sigma}{dcos(\theta_{\pi,\mu})}$ ($\mathrm{cm}^2$)"
        )
        
    @classmethod
    def num_protons(cls):
        return cls(
            var_save_name="num_protons",
            var_plot_name=r"# protons",
            var_unit="",
            bins=np.linspace(0, 3, 4),
            var_evt_reco_col=('slc', 'measure_var', 'num_protons', ''),
            var_evt_truth_col=('truth', 'true_var', 'num_protons', ''),
            var_nu_col=('truth', 'true_var', 'num_protons', ''),
            xsec_label=r"$\frac{d\sigma}{dn}$ ($\mathrm{cm}^2$)"
        )

        
    @classmethod
    def delta_pt(cls):
        return cls(
            var_save_name="delta_pt",
            var_plot_name=r"$\delta p_T$",
            var_unit="[GeV]",
            bins= np.array([0., 0.158, 0.372, 0.8]),
            var_evt_reco_col=('slc', 'measure_var', 'delta_pT', ''),
            var_evt_truth_col=('truth', 'true_var', 'delta_pT', ''),
            var_nu_col=('truth', 'true_var', 'delta_pT', ''),
            xsec_label=r"$\frac{d\sigma}{d \delta p_T}$ ($\mathrm{cm}^2$)"
        )
    @classmethod
    def delta_alpha_T(cls):
        return cls(
            var_save_name="delta_alpha_T",
            var_plot_name=r"$\delta \alpha_T$",
            var_unit="[rad]",
            bins= np.array([0., 1.637, 2.643, np.pi]),
            var_evt_reco_col=('slc', 'measure_var', 'delta_alpha_T', ''),
            var_evt_truth_col=('truth', 'true_var', 'delta_alpha_T', ''),
            var_nu_col=('truth', 'true_var', 'delta_alpha_T', ''),
            xsec_label=r'$\frac{d\sigma}{d\delta \alpha_T}$ ($\mathrm{cm}^2$)'
        )
    
    @classmethod
    def delta_phi_T(cls):
        return cls(
            var_save_name="delta_phi_T",
            var_plot_name=r"$\delta\phi_T$",
            var_unit="[rad]",
            bins= np.array([0., 0.542, 1.96, np.pi]),
            var_evt_reco_col=('slc', 'measure_var', 'delta_phi_T', ''),
            var_evt_truth_col=('truth', 'true_var', 'delta_phi_T', ''),
            var_nu_col=('truth', 'true_var', 'delta_phi_T', ''),
            xsec_label=r"$\frac{d\sigma}{d\delta \phi_T}$ ($\mathrm{cm}^2$)"
        )
    
    @classmethod
    def list_all_configs(cls, print_summary=True):
        members = inspect.getmembers(cls, predicate=inspect.isroutine)
        config_methods = [name for name, func in members
                              if getattr(func, "__self__", None) == cls
                              and not name.startswith("_")
                              and name not in ("list_all_configs")]
        print("Available VariableConfig options:")
        for name in config_methods:
            config = getattr(cls, name)()
            if print_summary:
                print(f"  {name:20}")
                print(f"    var_save_name   : {config.var_save_name}")
                print(f"    var_plot_name   : {config.var_plot_name}")
                print(f"    var_labels      : {config.var_labels}")
                print(f"    bins            : {np.array2string(config.bins, separator=', ')}")
                print(f"    bin_centers     : {np.array2string(config.bin_centers, separator=', ')}")
                print(f"    var_evt_reco_col: {config.var_evt_reco_col}")
                print(f"    var_evt_truth_col: {config.var_evt_truth_col}")
                print(f"    var_nu_col      : {config.var_nu_col}")
                print(f"    xsec_label      : {config.xsec_label}")
                print()
            else:
                print(f"  {name:20}")    
        return config_methods  