# Further MNN algorithm development

## Overview

This repository contains code for further development of the mutual nearest neighbours (MNN) batch correction method, as implemented in the `mnnCorrect` and `fastMNN` functions in the [_scran_](https://bioconductor.org/packages/scran) package.
It is based on the code at https://github.com/MarioniLab/MNN2017, which accompanies the paper **Batch effects in single-cell RNA-sequencing data are corrected by matching mutual nearest neighbors** by [Haghverdi _et al. (2018)_](https://doi.org/10.1038/nbt.4091).

## Simulations

To run the simulations, enter the `simulations/` directory and run:

- `cluster_sim.R`, which simulates a variety of scenarios involving orthogonal batch effects.
- `nonorth_sim.R`, which simulates some pathological non-orthogonal batch effects.

## Real data

Three real data analyses are available - `haematopoiesis`, `pancreas` and `droplet`.
Each subdirectory will usually contain:

- `prepareData.R`, to download, pre-process and normalize the data.
- `plotCorrections.R`, to perform the batch correction and visualize the result with t-SNE plots.

Data file downloads are performed using the [_BiocFileCache_](https://bioconductor.org/packages/BiocFileCache) package to save time and bandwidth.
This will cache the files locally after the initial download, and reuse the cached versions when the script is re-run. 
