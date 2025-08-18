# AF3-Evaluation
The Revised Repository for AlphaFold3 in Drug Discovery: A Comprehensive Assessment of Capabilities, Limitations, and Applications.

# Key changes
1. Most of the workflows were moved to SBGrid's Schordinger. If you are using different versions of Schrodinger, you need to change the launching/activation command to utilize your shcrodinger. 
2. Removed analysis based on reviwers' comments 
3. You need to have access to Schrodinger


# example-analysis-schrodinger.ipynb - use environment_analysis-struct.yml
This notebook collects, prepares structures and conducts analysis (using mannually written align-batch-schrodinger-smililarity.py and poseviewer_interactions.py by run_poseviewer_interactions_dir.sh).

# align-batch-schrodinger-smililarity.py - SBGrid environment needed (or Schordinger environment if using a different Schrodinger)
This script uses many functions from schordinger to conduct benchmark

# CaSR_evaluation.ipynb - use environment_analysis-struct.yml
This notebook calculates the protein RMSD for CaSR prediction evaluation.

# input_compile.ipynb - use environment_AF3.yml
This notebook contains the AF3 input file compilation functions that we used. Example use has been included in the script.

# ternary-schrodinger.ipynb - use environment_analysis-struct.yml
This notebook uses similar functions with example-analysis-schrodinger.ipynb, but has additional functions for comparing ternary PPI interactions using protein_interaction_analysis.py from Schrodinger

# run_prep_dir_cif.sh and run_prep_dir_pdb.sh
These scripts deal with preparation of structures using Schordinger prepwizard

# run_poseviewer_interactions_dir.sh
This script run poseviewer_interactions.py for a given directory

# evaluation-datasets
Structures from PLINDER used in benchmarking of AF3. Separated into three datasets.

# Installation
To use scripts from this directory, SBGrid is needed for all the Schrodinger functions. Join SBGrid at https://sbgrid.org/join. For conda environment, build the environment with the corresponding .yml. Installation time depends on conda. The code should work on all modern hardwares.  