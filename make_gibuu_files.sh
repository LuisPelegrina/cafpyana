#!/bin/bash

# Configuration
INPUT_LIST="/home/lpelegri/cafpyana/data_lists/GIBUU/mc_SBND2026A_gen1_prodgenie_corsika_proton_rockbox_sbnd_GIBUU_CV_v10_06_00_09_flatcaf_sbnd.list"
CONFIG="analysis_village/cc1pi/configs/cc1pi_no_systematics_df.py"
SCRATCH_DIR="/scratch/7DayLifetime/lpelegrina/GiBUUFiles"
N_SPLITS=50

mkdir -p $SCRATCH_DIR

# 1. Split the master list into temporary sub-lists
# This creates files: xaa, xab, xac... in the scratch dir
total_lines=$(wc -l < "$INPUT_LIST")
lines_per_file=$(( (total_lines + N_SPLITS - 1) / N_SPLITS ))

echo "Splitting $INPUT_LIST ($total_lines lines) into chunks of $lines_per_file..."
split -l $lines_per_file -d --additional-suffix=.list "$INPUT_LIST" "$SCRATCH_DIR/sublist_"

# 2. Loop through the created lists and run the python command
for sublist in $SCRATCH_DIR/sublist_*.list; do
    # Extract the number from the filename (e.g., 00, 01)
    num=$(basename "$sublist" .list | sed 's/sublist_//')
    output_name="$SCRATCH_DIR/cc1pi_GIBUU_subfile$num"
    
    echo "----------------------------------------------------"
    echo "Processing Chunk $num: $sublist"
    echo "----------------------------------------------------"
    
    # Run the command
    python run_df_maker.py -c "$CONFIG" -l "$sublist" -o "$output_name"
    
    # Optional: remove the temp sublist after success
    if [ $? -eq 0 ]; then
        rm "$sublist"
    fi
done

echo "Bash processing complete."