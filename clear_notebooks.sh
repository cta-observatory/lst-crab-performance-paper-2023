#!/bin/bash

set -eu -o pipefail

# Check if jq is installed
if ! command -v jq &> /dev/null; then
    echo "jq could not be found"
    exit 1
fi

# Check if notebook_directory is supplied
if [ -z "${1:-}" ]; then
    echo "Please provide the directory containing your Jupyter notebooks as an argument."
    exit 1
fi

# Replace the following path with the directory containing your Jupyter notebooks
notebook_directory="$1"

# Check if the notebook_directory exists
if [ ! -d "$notebook_directory" ]; then
    echo "Directory $notebook_directory does not exist."
    exit 1
fi

# Find all .ipynb files and execute them in sorted order
find "$notebook_directory" -iname "*.ipynb" -type f -print0 | sort -z | while read -d $'\0' notebook; do
    echo "Clearing $notebook ..."
    
    jupyter-nbconvert --inplace --clear-output "$notebook"
    
    tmp_notebook=$(mktemp)
    jq 'del(.metadata.kernelspec)' "$notebook" > "$tmp_notebook" && mv "$tmp_notebook" "$notebook"
    
    echo "Successfully cleared $notebook"
done

