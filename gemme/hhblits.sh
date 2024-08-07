#!/usr/bin/bash
# This script is structured to optimize the number of sequences obtained by hhblits in the MSA by gradually
# relaxing the search criteria if the desired number of sequences is not achieved in the initial runs.

# Load conda and mamba environment scripts
source /home/aliciapageau/mambaforge/etc/profile.d/conda.sh
source /home/aliciapageau/mambaforge/etc/profile.d/mamba.sh

# Activate the HH-suite3 environment
mamba activate HH-suite3

# Set the input and output directories from the script arguments
INPUT_DIR=$1
OUTPUT_DIR=$2

# Change to the output directory, or exit if it fails
cd "${OUTPUT_DIR}" || exit 1

# Loop through all .fasta files in the INPUT_DIR
for i in `ls ${INPUT_DIR}*.fasta`
do
    # Run hhblits with 2 iterations for the first run
    hhblits -i "$i" -oa3m "${OUTPUT_DIR}/$(basename "$i" .fasta).a3m" \
    -d /home/aliciapageau/hhsuiteDB/UniRef30_2022_02 -cpu 24 -qid 20 -id 98 -cov 80 -e 0.0001 -n 2 -v 0

    # Check the number of sequences in the resulting MSA (Multiple Sequence Alignment)
    num_sequences=$(grep -c '^>' "${OUTPUT_DIR}/$(basename "$i" .fasta).a3m")

    # If the number of sequences is less than 200, run hhblits again with different parameters
    if [ "$num_sequences" -lt 200 ]; then
        hhblits -i "${OUTPUT_DIR}/$(basename "$i" .fasta).a3m" -oa3m "${OUTPUT_DIR}/$(basename "$i" .fasta)_v2.a3m" \
        -d /home/aliciapageau/hhsuiteDB/UniRef30_2022_02 -cpu 24 -qid 18 -id 99 -cov 70 -e 0.001 -n 1 -v 0

        # Check the number of sequences in the new MSA
        num_sequences=$(grep -c '^>' "${OUTPUT_DIR}/$(basename "$i" .fasta)_v2.a3m")
        
        # If the number of sequences is still less than 200, run hhblits again
        if [ "$num_sequences" -lt 200 ]; then
            hhblits -i "${OUTPUT_DIR}/$(basename "$i" .fasta)_v2.a3m" -oa3m "${OUTPUT_DIR}/$(basename "$i" .fasta)_v3.a3m" \
            -d /home/aliciapageau/hhsuiteDB/UniRef30_2022_02 -cpu 24 -qid 15 -id 99 -cov 60 -e 0.01 -n 1 -v 0
            
            # Check the number of sequences in the new MSA
            num_sequences=$(grep -c '^>' "${OUTPUT_DIR}/$(basename "$i" .fasta)_v3.a3m")
            
            # If the number of sequences is still less than 200, run hhblits one last time
            if [ "$num_sequences" -lt 200 ]; then
                hhblits -i "${OUTPUT_DIR}/$(basename "$i" .fasta)_v3.a3m" -oa3m "${OUTPUT_DIR}/$(basename "$i" .fasta)_v4.a3m" \
                -d /home/aliciapageau/hhsuiteDB/UniRef30_2022_02 -cpu 24 -qid 12 -id 100 -cov 50 -e 0.1 -n 1 -v 0
            fi
        fi
    fi
done

# Deactivate the HH-suite3 environment
mamba deactivate

