#!/usr/bin/bash
# This script automates the process of running GEMME analysis on multiple .a3m files by utilizing Docker to ensure a consistent environment.

# Loop through all .a3m files in the specified directory
for file in $1/*.a3m; do
    # Create a directory with the same name as the .a3m file (without the extension)
    mkdir $(basename $file .a3m)
    
    # Run a Docker container with the 'elodielaine/gemme:gemme' image
    docker run -ti --rm --mount type=bind,source=$PWD,target=/project elodielaine/gemme:gemme bash -c "\
        # Convert the alignment to another format
        convertAli.sh $(basename $file .a3m) && \
        # Change directory to the newly created directory
        cd $(basename $file .a3m) && \
        # Run the GEMME analysis on the .fasta file
        python2.7 \$GEMME_PATH/gemme.py ../$(basename $file .a3m).fasta -r input -f ../$(basename $file .a3m).fasta &&\
        # Exit the Docker container
        exit
    "
done

