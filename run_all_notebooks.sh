#!/bin/bash

# Replace the following path with the directory containing your Jupyter notebooks
notebook_directory="$1"

# Find all .ipynb files and execute them in sorted order
find "$notebook_directory" -iname "*.ipynb" -type f -print0 | sort -z | while read -d $'\0' notebook
do
    echo "Running $notebook ..."
    jupyter-nbconvert --to notebook --execute --inplace "$notebook"
    if [ $? -eq 0 ]; then
        echo "Successfully ran $notebook"
    else
        echo "Error running $notebook"
        exit 1
    fi
done
