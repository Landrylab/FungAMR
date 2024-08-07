#!/usr/bin/bash
# This script will run a foldX scan of the whole protein with mutateX without repair on compute canada server.
# FoldX compute ddG of stabillity

# Your input directory should have :
#     1) PDB file (Must not have any small molecule not recognized by foldx)
#     2) Mutations list file (or positions list + residue types for saturation mutagenesis)
#     4) Configuration run file (default : flexddg_talaris2014)
#     5) Configuration settings file (default : valeria_rosettampi)
# Output will be in the flexddg folder created in your running directory


# Function to display usage instructions
usage() {
  echo "Usage: $(basename $0) [-h] [-i INPUT_DIR] [-o OUTPUT_DIR] [-p PDB_DIR] [--pdb-list PDB_LIST] [--nruns NRUNS] [--clean CLEAN] [-n NUM_PROCESSES] [--mem MEM] [--time TIME]"
  echo "Options:"
  echo "  -h, --help           Show this help message and exit"
  echo "  -i, --input          Input directory (default: current directory/input)"
  echo "  -o, --output         Output directory (default: current directory)"
  echo "  -p, --pdb-dir        PDB directory (default: current directory/repair)"
  echo "      --pdb-list       Text file with PDB ids list to compute"
  echo "      --nruns          Number of foldx runs (default: 5)"
  echo "      --clean          {partial, deep, none} (default: deep)"
  echo "  -n, --np             Number of process to request from SLURM (default: 64)"
  echo "      --mem            Memory request from SLURM (default: 32G)"
  echo "      --time           Time limit (default: 1-00:00:00)"
  exit 1
}

# Parse command line arguments using getopt
VALID_ARGS=$(getopt -o hi:o:p:n: --long help,input:,output:,pdb-dir:,pdb-list:,nruns:,poslist:,clean:,np:,mem:,time: -- "$@")
if [[ $? != 0 ]]; then
    exit 1;
fi

eval set -- "$VALID_ARGS"

# Default values for variables
INPUT_DIR=$PWD/input
OUTPUT_DIR=$PWD
PDB=$PWD/repair
CLEAN=deep
NRUNS=5
NP=64
MEM=32G
TIME=1-00:00:00

# Loop to handle the arguments passed to the script
while true; do
  case "$1" in
    -h | --help ) usage; shift;;
    -i | --input ) INPUT_DIR="$2"; shift 2;;
    -o | --output ) OUTPUT_DIR="$2"; shift 2;;
    -p | --pdb-dir ) PDB="$2"; shift 2;;
    --pdb-list ) PDB_LIST="$2"; shift 2;;
    --np ) NP="$2"; shift 2;;
    --nruns ) NRUNS="$2"; shift 2;;
    --clean ) CLEAN="$2"; shift 2;;
    --mem ) MEM="$2"; shift 2;;
    --time ) TIME="$2"; shift 2;;
    -- ) shift; break;;
    * ) break ;;
  esac
done

cat > mutatex.sbatch << EOF
#!/bin/bash

#SBATCH -D $PWD
#SBATCH -J mutateX
#SBATCH -o mutateX.out
#SBATCH --account=def-clandry
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=$NP
#SBATCH --mem=$MEM
#SBATCH --time=$TIME
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alicia.pageau.1@ulaval.ca

module load StdEnv/2023
module load python/3.10.2
virtualenv --no-download ENV_mutatex
source ENV_mutatex/bin/activate
pip install --no-index --upgrade pip
pip install --no-index -r /home/alpag18/ENV_mutatex/requirements.txt
pip install /home/alpag18/ENV_mutatex/MutateX-0.8-py3-none-any.whl
 
# Check if the PDB_LIST file exists and is readable
if [ -r $PDB_LIST ]; then
    # Read each line from the file into the array 'ids'
    while IFS= read -r line; do
        # Append the line (PDB ID) to the array
        ids+=("$line")
    done < $PDB_LIST
else
    # If the PDB_LIST file doesn't exist or isn't readable, set 'ids' to an empty array
    ids=()
fi

cd $OUTPUT_DIR

# Loop over all PDB files found in the subdirectories of the PDB directory
for pdb_file in \$(find $PDB/*/ -type f -name '*.pdb'); do
    # Extract the PDB ID from the filename (basename without the .pdb extension)
    pdb_id=\$(basename "\$pdb_file" .pdb)

    # Check if the priority list is empty or if the PDB ID is in the priority list
    if [ \${#ids[@]} -eq 0 ] || [[ " \${ids[@]} " =~ " \$pdb_id " ]]; then

    # Run mutateX saturation mutagenesis scan
	mutatex \$pdb_file \
	--np $NP \
	-f suite5 \
	-M $INPUT_DIR/mutate_runfile_template.txt \
	--foldx-log \
	--skip-repair \
	--clean $CLEAN
    fi
done
EOF

# Submit the SBATCH job script
sbatch mutatex.sbatch
