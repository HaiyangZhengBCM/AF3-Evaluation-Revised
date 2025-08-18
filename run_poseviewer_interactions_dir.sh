#!/bin/bash
source /programs/sbgrid.shrc
# Check if the user provided an input folder as an argument
if [ -z "$1" ]; then
    echo "Usage: $0 <path_to_pdb_folder>"
    exit 1
fi

# Assign the input folder to a variable
PDB_FOLDER="$1"

cd "$PDB_FOLDER"

# Check if the specified folder exists
if [ ! -d "$PDB_FOLDER" ]; then
    echo "Error: The specified folder does not exist: $PDB_FOLDER"
    exit 1
fi

# Find all .pdb files in the folder
pdb_files=($(find "$PDB_FOLDER" -maxdepth 1 -type f -iname "*.pdb"))

# Check if any .pdb files are found
if [ ${#pdb_files[@]} -eq 0 ]; then
    echo "No .pdb files found in $PDB_FOLDER."
    exit 1
fi

# Loop through each found .pdb file
for pdb_file in "${pdb_files[@]}"; do
    # Extract the filename without the extension
    base_name=$(basename "$pdb_file" .pdb)

    # Define the output log file (optional)
    output_log="${base_name}_poseviewer.log"

    # Run the poseviewer_interactions.py script
    echo "Processing file: $pdb_file"
    $SCHRODINGER/run poseviewer_interactions.py "$pdb_file" > "$output_log" 2>&1

    if [ $? -eq 0 ]; then
        echo "Successfully processed: $pdb_file. Log saved to $output_log."
    else
        echo "Error processing: $pdb_file. Check $output_log for details."
    fi
done

