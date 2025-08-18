#!/bin/bash
source /programs/sbgrid.shrc
# Check if the user provided an input folder as an argument
if [ -z "$1" ]; then
    echo "Usage: $0 <path_to_cif_dir>"
    exit 1
fi

# Assign the input folder to a variable
CIF_DIR="$1"
OUTPUT_FOLDER="$2"

# Check if the specified folder exists
if [ ! -d "$CIF_DIR" ]; then
    echo "Error: The specified file does not exist: $CIF_DIR"
    exit 1
fi

cd "$OUTPUT_FOLDER"
# Extract the filename without the extension
# Find all .cif files in the folder
cif_files=($(find "$CIF_DIR" -maxdepth 1 -type f -iname "*.cif"))

# Check if any .cif files are found
if [ ${#cif_files[@]} -eq 0 ]; then
    echo "No .cif files found in $CIF_DIR."
    exit 1
fi

# Loop through each found .cif file
for cif_file in "${cif_files[@]}"; do
    # Extract the filename without the extension
    base_name=$(basename "$cif_file" .cif)

    # Define the output log file (optional)
    output_log="${OUTPUT_FOLDER}/${base_name}_prep.log"
    output_pdb="${OUTPUT_FOLDER}/${base_name}_prep.pdb"
    if [ -f "$output_pdb" ]; then
        echo "Preped file already exists: $output_pdb. Skipping."
        continue
    fi
    # Run the poseviewer_interactions.py script
    echo "Processing file: $cif_file"
    $SCHRODINGER/utilities/prepwizard -noccd -assign_all_residues -rehtreat -epik_pH 7.4 -propka_pH 7.4 -noimpref "$cif_file" "$output_pdb" > "$output_log" 2>&1

    if [ $? -eq 0 ]; then
        echo "Successfully processed: $cif_file. Log saved to $output_log."
    else
        echo "Error processing: $cif_file. Check $output_log for details."
    fi
done


