# LST-1 Crab performance paper

This repository includes code and data to reproduce the figures (only the data points, not the entire data analysis) of the paper entitled _Observations of the Crab Nebula and Pulsar with the Large-Sized Telescope prototype of the Cherenkov Telescope Array_ for the CTA LST collaboration.

[![arXiv](https://img.shields.io/badge/arXiv-2306.12960-b31b1b.svg)](https://arxiv.org/abs/2306.12960)


## How to run this repository

### On mybinder

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/cta-observatory/lst-crab-performance-paper-2023/HEAD)

You may visit https://mybinder.org/v2/gh/cta-observatory/lst-crab-performance-paper-2023/HEAD to run this repository online.

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

### Running the notebooks

The repository is organised following the paper sections.    
Each notebook can be run to reproduce the paper's figures.

[![Run Jupyter Notebooks](https://github.com/cta-observatory/lst-crab-performance-paper-2023/actions/workflows/run_notebooks.yml/badge.svg?branch=main)](https://github.com/cta-observatory/lst-crab-performance-paper-2023/actions/workflows/run_notebooks.yml)


## Cite us 

The preferred citation is the published article associated with this repository.

However, if you use the data points or code for this repository for your work, you may in addition cite the exact version of this repository from Zenodo.

[![DOI](https://zenodo.org/badge/617083935.svg)](https://zenodo.org/badge/latestdoi/617083935)
