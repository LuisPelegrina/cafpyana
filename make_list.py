import subprocess
import os
import sys

sam_definitions = {
    
    "mc_MCP2025B_5e18_v10_06_00_09_prodgenie_corsika_proton_rockbox_sbnd_CV_caf_flat_caf_sbnd": 
"data_lists/5e18POT/mc/mc_MCP2025B_5e18_v10_06_00_09_prodgenie_corsika_proton_rockbox_sbnd_CV_caf_flat_caf_sbnd.list",

  "mc_MCP2025B_5e18_v10_06_00_09_prodcorsika_proton_intime_sbnd_CV_caf_flat_caf_sbnd": 
    "data_lists/5e18POT/mc/mc_in_time_cosmics_5e18.list",
    
    "mc_SBND2026A_gen1_prodgenie_corsika_proton_rockbox_sbnd_GIBUU_CV_v10_06_00_09_flatcaf_sbnd": 
"data_lists/GIBUU/mc_SBND2026A_gen1_prodgenie_corsika_proton_rockbox_sbnd_GIBUU_CV_v10_06_00_09_flatcaf_sbnd.list",
            
  "mc_MCP2025B_v10_06_00_09_prodgenie_corsika_proton_rockbox_lowenergydirt_sbnd_CV_caf_flat_caf_sbnd": 
"data_lists/1e20/CV/mc_MCP2025B_v10_06_00_09_prodgenie_corsika_proton_rockbox_lowenergydirt_sbnd_CV_caf_flat_caf_sbnd.list",

  
    "data_MCP2025C_Spring25_reprocess_FixedDev_bnblight_v10_06_00_09_flatcaf_sbnd": 
"data_lists/5e18POT/data/data_MCP2025C_Spring25_reprocess_fixed_Dev_bnblight_v10_06_00_09_flatcaf_sbnd.list",
        
    "data_MCP2025C_Spring25_reprocess_RollingDev_bnblight_v10_06_00_09_flatcaf_sbnd": 
"data_lists/5e18POT/data/data_MCP2025C_Spring25_reprocess_rolling_Dev_bnblight_v10_06_00_09_flatcaf_sbnd.list",
    
    "data_MCP2025C_Spring25_reprocess_Intime_offbeamlight_v10_06_00_09_flatcaf_sbnd": 
"data_lists/5e18POT/data/data_MCP2025C_Spring25_reprocess_Intime_offbeamlight_v10_06_00_09_flatcaf_sbnd.list",

    "mc_MCP2025B_1e20_10_prodgenie_corsika_proton_rockbox_sbnd_SystVar_CV_caf_flat_caf_sbnd": 
"data_lists/1e20/SystVar/mc_MCP2025B_1e20_10_prodgenie_corsika_proton_rockbox_sbnd_SystVar_CV_caf_flat_caf_sbnd.list",
        
    "mc_MCP2025B_1e20_10_prodgenie_corsika_proton_rockbox_sbnd_SystVar_0xSCE_caf_flat_caf_sbnd": 
"data_lists/1e20/SystVar/mc_MCP2025B_1e20_10_prodgenie_corsika_proton_rockbox_sbnd_SystVar_0xSCE_caf_flat_caf_sbnd.list",

    "mc_MCP2025B_1e20_10_prodgenie_corsika_proton_rockbox_sbnd_SystVar_2xSCE_caf_flat_caf_sbnd": 
"data_lists/1e20/SystVar/mc_MCP2025B_1e20_10_prodgenie_corsika_proton_rockbox_sbnd_SystVar_2xSCE_caf_flat_caf_sbnd.list",

    "mc_MCP2025B_1e20_10_prodgenie_corsika_proton_rockbox_sbnd_SystVar_PMTGainFluct_caf_flat_caf_sbnd": 
"data_lists/1e20/SystVar/mc_MCP2025B_1e20_10_prodgenie_corsika_proton_rockbox_sbnd_SystVar_PMTGainFluct_caf_flat_caf_sbnd.list",
    
    "mc_MCP2025B_1e20_10_prodgenie_corsika_proton_rockbox_sbnd_SystVar_PMTLowEff_caf_flat_caf_sbnd": 
"data_lists/1e20/SystVar/mc_MCP2025B_1e20_10_prodgenie_corsika_proton_rockbox_sbnd_SystVar_PMTLowEff_caf_flat_caf_sbnd.list",
        "mc_MCP2025B_1e20_10_prodgenie_corsika_proton_rockbox_sbnd_SystVar_PMTHighNoise_caf_flat_caf_sbnd": 
"data_lists/1e20/SystVar/mc_MCP2025B_1e20_10_prodgenie_corsika_proton_rockbox_sbnd_SystVar_PMTHighNoise_caf_flat_caf_sbnd.list",
    "mc_MCP2025B_1e20_10_prodgenie_corsika_proton_rockbox_sbnd_SystVar_WireMod_XThetaXW_caf_flat_caf_sbnd": 
"data_lists/1e20/SystVar/mc_MCP2025B_1e20_10_prodgenie_corsika_proton_rockbox_sbnd_SystVar_WireMod_XThetaXW_caf_flat_caf_sbnd.list",

    "mc_MCP2025B_1e20_10_prodgenie_corsika_proton_rockbox_sbnd_SystVar_WireMod_YZ_caf_flat_caf_sbnd": 
"data_lists/1e20/SystVar/mc_MCP2025B_1e20_10_prodgenie_corsika_proton_rockbox_sbnd_SystVar_WireMod_YZ_caf_flat_caf_sbnd.list",
        
}

def generate_lists(definitions_dict):
    script_path = os.path.abspath("data/sample_lists/create_list.py")
    
    # Create a copy of the current environment and add SAM_EXPERIMENT
    env = os.environ.copy()
    env["SAM_EXPERIMENT"] = "sbnd"

    for sam_def, out_path in definitions_dict.items():
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        print(f"Querying SAM for: {sam_def}...")
        
        try:
            with open(out_path, "w") as f:
                subprocess.run(
                    [sys.executable, script_path, sam_def], # -e is now redundant but fine
                    stdout=f,
                    stderr=subprocess.PIPE,
                    check=True,
                    text=True,
                    env=env # <--- Pass the environment here
                )
            print(f"Successfully created: {out_path}")
        except subprocess.CalledProcessError as e:
            print(f"Error processing {sam_def}:")
            print(e.stderr)

    input_files = [
            "data_lists/5e18POT/data/data_MCP2025C_Spring25_reprocess_fixed_Dev_bnblight_v10_06_00_09_flatcaf_sbnd.list",
            "data_lists/5e18POT/data/data_MCP2025C_Spring25_reprocess_rolling_Dev_bnblight_v10_06_00_09_flatcaf_sbnd.list"
        ]
        
    # Define the combined output file
    output_file = "data_lists/5e18POT/data/data_MCP2025C_Spring25_reprocess_Dev_bnblight_v10_06_00_09_flatcaf_sbnd.list"
    
    # Ensure the output directory exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    try:
        with open(output_file, 'w') as outfile:
            for fname in input_files:
                print(f"Reading {fname}...")
                with open(fname, 'r') as infile:
                    # Read lines, strip potential trailing empty lines, and write
                    for line in infile:
                        if line.strip():  # Only write non-empty lines
                            outfile.write(line.strip() + "\n")
                            
        print(f"\nSuccessfully combined files into:\n{output_file}")
    
    except FileNotFoundError as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    generate_lists(sam_definitions)