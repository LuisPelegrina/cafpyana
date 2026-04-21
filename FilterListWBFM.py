import uproot
import os

input_list = "data_lists/AR23p.list"
output_list = "data_lists/AR23p_filtered.list"
target_branch = "rec.slc.barycenterFM.chi2"

def check_file_for_branch(file_path):
    try:
        # We only open the file and the tree header, not the data
        with uproot.open(file_path) as f:
            if "recTree" not in f:
                print(f"[-] Skip: {file_path} (No recTree found)")
                return False
            
            # tree.keys() only reads the metadata
            tree = f["recTree"]
            if target_branch in tree.keys():
                return True
            else:
                print(f"[-] Skip: {file_path} (Missing {target_branch})")
                return False
    except Exception as e:
        print(f"[!] Error opening {file_path}: {e}")
        return False

# 1. Read the current list
with open(input_list, "r") as f:
    all_files = [line.strip() for line in f if line.strip()]

print(f"Checking {len(all_files)} files...")

# 2. Filter
valid_files = []
for i, path in enumerate(all_files):
    if i % 10 == 0:
        print(f"Progress: {i}/{len(all_files)}")
        
    if check_file_for_branch(path):
        valid_files.append(path)

# 3. Save new list
with open(output_list, "w") as f:
    for path in valid_files:
        f.write(path + "\n")

print(f"\nDone! Filtered {len(all_files)} down to {len(valid_files)} files.")
print(f"New list saved to: {output_list}")