#!/bin/bash

# Replace the following path with the directory containing your Jupyter notebooks
notebook_directory="$1"

# Find all .ipynb files and execute them in sorted order
find "$notebook_directory" -iname "*.ipynb" -type f -print0 | sort -z | while read -d $'\0' notebook
do
    echo "Clearing $notebook ..."
    jupyter-nbconvert --inplace --clear-output "$notebook"
    jq 'del(.metadata.kernelspec)' "$notebook" > notebook_no_kernel.ipynb
    mv notebook_no_kernel.ipynb "$notebook"
    # jupyter kernelspec remove "$notebook"

    if [ $? -eq 0 ]; then
        echo "Successfully cleared $notebook"
    else
        echo "Error running $notebook"
        exit 1
    fi
done
