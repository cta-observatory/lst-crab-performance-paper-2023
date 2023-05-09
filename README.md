# LST Crab performance paper

This repository includes code and data to reproduce the figures of the paper entitled _Observations of the Crab Nebula and Pulsar with the Large-Sized Telescope prototype of the Cherenkov Telescope Array_ for the CTA LST collaboration.


## How to run this repository

### On mybinder

You may vist https://mybinder.org/v2/gh/cta-observatory/lst-crab-performance-paper-2023/HEAD to run this repository online.

### On your computer

You may run the following on mybinder or on your own laptop/computer/server by installing the provided conda environment:

[Install conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) and create the crablst1 environment using the following command:

```
conda env create -f environment.yml
```

Then activate the crablst1 environment and run the notebooks

```
conda activate crablst1
jupyter-lab
```

## Running the notebooks

The repositoy is organised following the paper sections.    
Each notebook can be run to reproduce the paper's figures.

[![Run Jupyter Notebooks](https://github.com/cta-observatory/lst-crab-performance-paper-2023/actions/workflows/run_notebooks.yml/badge.svg?branch=main)](https://github.com/cta-observatory/lst-crab-performance-paper-2023/actions/workflows/run_notebooks.yml)
