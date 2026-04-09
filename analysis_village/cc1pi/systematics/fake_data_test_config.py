import numpy as np
from scipy.special import erf

class FakeDataWeights:
    def __init__(self, mc_evt_df, mc_nu_df, var_config):
        self.mc_evt_df = mc_evt_df
        self.mc_nu_df = mc_nu_df
        self.var_config = var_config

    def get_weights(self, test_name, **kwargs):
        # Always start from ones (unless otherwise needed)
        # TODO: place normalization here?
        weights_fake_data = np.ones(len(self.mc_evt_df))
        weight_fakedata_signal_truth = np.ones(len(self.mc_nu_df[self.mc_nu_df.truth.nu_categ == "CC1pi"]))

        print(test_name)
        print(test_name == "mec_test")
        # MEC normalization
        if test_name.startswith("mec_test"):
            scale_factor = kwargs.get("scale_factor", 0.5)
            weights_fake_data[self.mc_evt_df.truth.genie_mode == 10] *= scale_factor
            weight_fakedata_signal_truth[self.mc_nu_df[self.mc_nu_df.truth.nu_categ == "CC1pi"].truth.genie_mode == 10] *= scale_factor
        # RES normalization
        elif test_name.startswith("res_test"):
            scale_factor = kwargs.get("scale_factor", 1.2)
            weights_fake_data[self.mc_evt_df.truth.genie_mode == 1] *= scale_factor
            weight_fakedata_signal_truth[self.mc_nu_df[self.mc_nu_df.truth.nu_categ == "CC1pi"].truth.genie_mode == 1] *= scale_factor
            
        # QE normalization
        elif test_name.startswith("qe_test"):
            scale_factor = kwargs.get("scale_factor", 1.5)
            weights_fake_data[self.mc_evt_df.truth.genie_mode == 0] *= scale_factor
            weight_fakedata_signal_truth[self.mc_nu_df[self.mc_nu_df.truth.nu_categ == "CC1pi"].truth.genie_mode == 0] *= scale_factor

        # np normalization
        elif test_name.startswith("1pi1p_test"):
            scale_factor = kwargs.get("scale_factor", 1.2)
            weights_fake_data[self.mc_evt_df.truth.nu_categ_proton_reduced == "1p_CC1Pi"] *= scale_factor   
            weight_fakedata_signal_truth[self.mc_nu_df[self.mc_nu_df.truth.nu_categ == "CC1pi"].truth.nu_categ_proton_reduced == "1p_CC1Pi"] *= scale_factor

        # sig normalization
        elif test_name.startswith("sig_test"):
            scale_factor = kwargs.get("scale_factor", 1.2)
            weights_fake_data[self.mc_evt_df.truth.nu_categ == "CC1pi"] *= scale_factor
            weight_fakedata_signal_truth[self.mc_nu_df[self.mc_nu_df.truth.nu_categ == "CC1pi"].truth.nu_categ == "CC1pi"] *= scale_factor

        # cos(theta) scale
        elif test_name.startswith("costh_weight_scale_"):
            scale_factor = kwargs.get("scale_factor", 0.7)
            weights_fake_data[self.mc_evt_df.truth.mu.dir.z > 0.8] *= scale_factor
            weight_fakedata_signal_truth[self.mc_nu_df[self.mc_nu_df.truth.nu_categ == "CC1pi"].truth.mu.dir.z > 0.8] *= scale_factor
 
        # Pion P tilt
        elif test_name.startswith("pion_P_weight_scale_"):
            scale_factor = kwargs.get("scale_factor", 0.3)
            col_x = ('truth','cpi','genp', 'x', '', '')
            col_y = ('truth','cpi','genp', 'y', '', '')
            col_z = ('truth','cpi','genp', 'z', '', '')
            col_tot = ('truth','cpi', 'totp', '', '','') # La nueva columna
            
            self.mc_evt_df[col_tot] = np.sqrt(
                self.mc_evt_df[col_x]**2 + 
                self.mc_evt_df[col_y]**2 + 
                self.mc_evt_df[col_z]**2
            )
            weights_fake_data[self.mc_evt_df.truth.cpi.totp < 0.4] *= scale_factor
            
        
            # Calculamos el momento total
            self.mc_nu_df[col_tot] = np.sqrt(
                self.mc_nu_df[col_x]**2 + 
                self.mc_nu_df[col_y]**2 + 
                self.mc_nu_df[col_z]**2
            )
            weight_fakedata_signal_truth[self.mc_nu_df[self.mc_nu_df.truth.nu_categ == "CC1pi"].truth.cpi.totp < 0.4] *= scale_factor
            

        elif test_name.startswith("pion_P_tilt_alpha_"):
            scale_factor = kwargs.get("scale_factor", 0.3)
            col_x = ('truth','cpi','genp', 'x', '', '')
            col_y = ('truth','cpi','genp', 'y', '', '')
            col_z = ('truth','cpi','genp', 'z', '', '')
            col_tot = ('truth','cpi', 'totp', '', '','') # La nueva columna
            
            self.mc_evt_df[col_tot] = np.sqrt(
                self.mc_evt_df[col_x]**2 + 
                self.mc_evt_df[col_y]**2 + 
                self.mc_evt_df[col_z]**2
            )
            
            P_pi_evt = self.mc_evt_df.truth.cpi.totp
            weights_fake_data = (np.ones(len(self.mc_evt_df)) + 
                     scale_factor * (P_pi_evt - P_pi_evt.mean()) / P_pi_evt.mean()
                    ).fillna(1)
           
            # Calculamos el momento total
            self.mc_nu_df[col_tot] = np.sqrt(
                self.mc_nu_df[col_x]**2 + 
                self.mc_nu_df[col_y]**2 + 
                self.mc_nu_df[col_z]**2
            )
            
            P_p_pi = self.mc_nu_df[self.mc_nu_df.truth.nu_categ == "CC1pi"].truth.cpi.totp
        
            weight_fakedata_signal_truth = np.ones(len(P_p_pi)) + scale_factor * (P_p_pi - P_p_pi.mean())/P_p_pi.mean()
            
        # Q2 tilt
        elif test_name.startswith("q2_test_alpha_"):
            scale_factor = kwargs.get("scale_factor", 0.3)
            Q2 = self.mc_evt_df.truth.Q2
            weights_fake_data *= np.ones(len(self.mc_evt_df)) + scale_factor * (Q2 - Q2.mean())/Q2.mean()
            weights_fake_data[np.isnan(weights_fake_data)] = 1
            assert np.isnan(weights_fake_data).sum() == 0
            Q2_nu = self.mc_nu_df[self.mc_nu_df.truth.nu_categ == "CC1pi"].truth.Q2
            weight_fakedata_signal_truth *= np.ones(len(Q2_nu)) + scale_factor * (Q2_nu - Q2_nu.mean())/Q2_nu.mean()


       
        #elif test_name.startswith("bump_"):
        #    bump_pos = kwargs.get("bump_pos", 0.6)
        #    bump_width = kwargs.get("bump_width", 0.0015)
        #    bump_height = kwargs.get("bump_height", 0.001)
        #    bump_height = bump_height*len(self.mc_evt_df)/len(self.var_config.bin_centers)
        #    bump_var_evt = self.mc_evt_df[self.var_config.var_evt_truth_col]
        #    weights_fake_data = np.ones(len(self.mc_evt_df)) + bump_height * np.exp(-0.5 * (bump_var_evt - bump_pos)**2 / bump_width**2)
        #    weights_fake_data[np.isnan(weights_fake_data)] = 1.
        #    bump_var_nu = self.mc_nu_df[self.mc_nu_df.truth.nu_categ == "CC1pi"][self.var_config.var_nu_col]
        #    weight_fakedata_signal_truth = np.ones(len(bump_var_nu)) + bump_height * np.exp(-0.5 * (bump_var_nu - bump_pos)**2 / bump_width**2)

        
        
        elif test_name.startswith("bump_"):
            # 1. Parameters
            bump_pos = kwargs.get("bump_pos", 0.6)
            bump_strength = kwargs.get("bump_strength", 0.1) # e.g., 0.1 for 10%
            bump_width = kwargs.get("bump_width", 0.1)
            
            # Histogram limits to handle truncation/leakage
            x_min, x_max = self.var_config.bins[0], self.var_config.bins[-1]
            sigma = (x_max - x_min) * bump_width
            
            # 2. Calculate the "Leakage" correction
            # This ensures that even if the bump is at the edge, the sum is correct
            fraction_inside = 0.5 * (erf((x_max - bump_pos) / (sigma * np.sqrt(2))) - 
                                     erf((x_min - bump_pos) / (sigma * np.sqrt(2))))
            
            # 3. Calculate target area
            wgt_col = ('slc', 'wgt', '', '', '', '')
            total_entries = self.mc_evt_df[wgt_col].sum()
            # The total weight we need to add is:
            target_added_weight = total_entries * bump_strength
            
            # 4. The Weighting Math
            # To ensure sum(extra_weight * original_weight) = target_added_weight,
            # we scale by the density of the points.
            bump_var_evt = self.mc_evt_df[self.var_config.var_evt_truth_col]
            
            # Basic Gaussian shape
            gauss_shape = np.exp(-0.5 * (bump_var_evt - bump_pos)**2 / sigma**2)
            
            # Normalization: target / (integral of shape)
            # We approximate the integral by the sum of (shape * original_weights)
            current_integral = (gauss_shape * self.mc_evt_df[wgt_col]).sum()
            
            if current_integral > 0:
                # This factor ensures that: sum(weights_fake_data * wgt_col) = total_entries + target_added_weight
                norm_factor = target_added_weight / current_integral
            else:
                norm_factor = 0
                
            weights_fake_data = 1.0 + (norm_factor * gauss_shape)
            weights_fake_data = np.nan_to_num(weights_fake_data, nan=1.0)

            # 5. Apply same logic to CC1pi Signal Truth
            nu_mask = self.mc_nu_df.truth.nu_categ == "CC1pi"
            total_entries_nu = self.mc_nu_df[nu_mask][wgt_col].sum()
            target_added_weight_nu = total_entries_nu * bump_strength
            
            bump_var_nu = self.mc_nu_df[nu_mask][self.var_config.var_nu_col]
            gauss_shape_nu = np.exp(-0.5 * (bump_var_nu - bump_pos)**2 / sigma**2)
            
            current_integral_nu = (gauss_shape_nu * self.mc_nu_df[nu_mask][wgt_col]).sum()
            
            if current_integral_nu > 0:
                norm_factor_nu = target_added_weight_nu / current_integral_nu
            else:
                norm_factor_nu = 0
                
            weight_fakedata_signal_truth = 1.0 + (norm_factor_nu * gauss_shape_nu)
            weight_fakedata_signal_truth = np.nan_to_num(weight_fakedata_signal_truth, nan=1.0)

       
        else:
            raise ValueError(f"Unknown test_name '{test_name}' provided.")
      

        # TODO
        # # Enhance tail
        # elif test_name.startswith("enhance_tail_alpha_"):
        #     alpha = kwargs.get("alpha", 0.7)
        #     weights_fake_data = np.ones(len(self.mc_evt_df))
        #     weights_fake_data[self.mc_evt_df.del_p > 0.25] = alpha
        #     weights_fake_data[np.isnan(weights_fake_data)] = 1.
        #     weight_fakedata_signal_truth = np.ones(len(self.mc_nu_df[self.mc_nu_df.nuint_categ == 1]))
        #     weight_fakedata_signal_truth[self.mc_nu_df[self.mc_nu_df.nuint_categ == 1].mc.del_p > 0.25] = alpha
        #     weight_fakedata_signal_truth[np.isnan(weight_fakedata_signal_truth)] = 1.


        return weights_fake_data, weight_fakedata_signal_truth