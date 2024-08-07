#!/usr/bin/bash
# This script will run 10 times foldX RepairPDB command on compute canada server.
# Runs are parallelise per PDB_DIR file

# Function to display usage instructions
usage() {
  echo "Usage: $(basename $0) [-h] [-i] [-o] -p PDB_DIR [-n NUM_PROCESSES] [--mem 1G] [--time 1-00:00:00]"
  echo "Options:"
  echo "  -h, --help    Show this help message and exit"
  echo "  -i, --in      Input directory (default: current directory/input)"
  echo "  -o, --out     Output directory (default: current directory/repair)"
  echo "  -p, --pdb     PDB directory (default: current directory/pdb)"
  echo "  -n, --np      Number of process to request for SLURM (default: 20)"
  echo "      --mem     Memory request from SLURM (default: 10G)"
  echo "      --time    Time limit (default: 05:00:00)"
  exit 1
}

# Parse command line arguments using getopt
VALID_ARGS=$(getopt -o hi:o:p:n: --long help,in:,out:,pdb:,np:,mem:,time: -- "$@")
if [[ $? != 0 ]]; then
    exit 1;
fi

eval set -- "$VALID_ARGS"

# Default values for variables
INPUT_DIR=$PWD/input
OUTPUT_DIR=$PWD/repair
PDB_DIR=$PWD/pdb
NUM_PROCESSES=20
TIME=05:00:00
MEM=""

# Loop to handle the arguments passed to the script
while true; do
  case "$1" in
    -h | --help ) usage; shift ;;
    -i | --in ) INPUT_DIR="$2"; shift 2;;
    -o | --out ) OUTPUT_DIR="$2"; shift 2;;
    -p | --pdb ) PDB_DIR="$2"; shift 2;;
    -n | --np ) NUM_PROCESSES="$2"; shift 2;;
         --mem ) MEM="$2"; shift 2;;
	 --time ) TIME="$2"; shift 2;;
    -- ) shift; break;;
    * ) break ;;
  esac
done

# Set default memory if not provided
if [ -z "$MEM" ]; then
  MEM=$(($NUM_PROCESSES/2))G # Default memory is half the number of processes in GB
fi

cat > run_repair.sbatch << EOF
#!/bin/bash

#SBATCH -D $PWD
#SBATCH -J repairPDB
#SBATCH -o repairPDB.out
#SBATCH --account=def-clandry
#SBATCH --ntasks-per-node=$NUM_PROCESSES
#SBATCH --mem=$MEM
#SBATCH --time=$TIME
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alicia.pageau.1@ulaval.ca

module load python/3.10.2
virtualenv --no-download ENV_mutatex
source ENV_mutatex/bin/activate
pip install --no-index --upgrade pip
pip install --no-index -r /home/alpag18/ENV_mutatex/requirements.txt
pip install /home/alpag18/ENV_mutatex/MutateX-0.8-py3-none-any.whl


# Function to run FoldX RepairPDB multiple times on a given PDB file
run_foldx_multiple() {
    cd $PWD
    prot=\$(basename \$1 .pdb)
    
    # Create a directory for the repaired structure outputs
    mkdir -p $OUTPUT_DIR/\${prot}_repair
    
    # Copy the original PDB file to the repair directory
    cp $PDB_DIR/\${prot}.pdb $OUTPUT_DIR/\${prot}_repair
    cd $OUTPUT_DIR/\${prot}_repair
    
    # Loop to run FoldX RepairPDB command 10 times
    for j in {1..10}; do
        # Run FoldX RepairPDB command with specific parameters
        foldx5 --command=RepairPDB --pdb=\${prot}.pdb \
            --ionStrength=0.05 \
            --temperature=298 \
            --pH=7 \
            --water=-PREDICT \
            --vdwDesign=2 \
            --pdbHydrogens=false \
	    --pdbWaters=true \
	    --complexWithDNA=false > \${prot}_Repair_\${j}.log
    
        wait
        
        # Create a subdirectory for each repair run
        mkdir \${prot}_repair_\${j}
        
        # Move the output files to the respective subdirectory
        mv \${prot}_Repair* \${prot}_repair_\${j}
        
        # Remove the previous PDB file
        rm \${prot}.pdb
        
        # Copy the newly repaired PDB file for the next iteration
        cp \${prot}_repair_\${j}/\${prot}_Repair.pdb ./\${prot}.pdb
    done
}

# Export the function so it can be used by GNU Parallel
export -f run_foldx_multiple

# List all PDB files in the PDB directory and run the `run_foldx_multiple` function in parallel
ls $PDB_DIR | grep '\.pdb$' | parallel run_foldx_multiple {}

EOF

# Submit the SBATCH job script
sbatch run_repair.sbatch
