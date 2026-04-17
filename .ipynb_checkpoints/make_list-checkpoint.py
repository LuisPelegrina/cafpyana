import subprocess
import os
import sys

sam_definitions = {
    "mc_MCP2025B_5e18_v10_06_00_09_prodgenie_corsika_proton_rockbox_sbnd_CV_caf_flat_caf_sbnd": 
    "data_lists/5e18POT/mc/mc_MCP2025B_5e18_v10_06_00_09_prodgenie_corsika_proton_rockbox_sbnd_CV_caf_flat_caf_sbnd.list",
}

def generate_lists(definitions_dict):
    # Get the absolute path to the helper script
    # This assumes create_list.py is in data/sample_lists/ relative to where you run this
    script_path = os.path.abspath("data/sample_lists/create_list.py")

    for sam_def, out_path in definitions_dict.items():
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        print(f"Querying SAM for: {sam_def}...")
        
        try:
            with open(out_path, "w") as f:
                # Use sys.executable to ensure we use the SAME python environment
                # that is currently running this master script.
                subprocess.run(
                    [sys.executable, script_path, sam_def, "-e", "sbnd"],
                    stdout=f,
                    stderr=subprocess.PIPE, # Capture errors for debugging
                    check=True,
                    text=True
                )
            print(f"Successfully created: {out_path}")
        except subprocess.CalledProcessError as e:
            print(f"Error processing {sam_def}:")
            print(e.stderr) # This will tell you exactly why SAM failed

if __name__ == "__main__":
    generate_lists(sam_definitions)