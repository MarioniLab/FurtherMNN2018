# FurtherMNN2017

## Overview

This repository contains code for further development of the mutual nearest neighbours (MNN) batch correction method, as implemented in the `mnnCorrect` and `fastMNN` functions in the [_scran_](https://bioconductor.org/packages/scran) package.
It is based on the code at https://github.com/MarioniLab/MNN2017, which accompanies the paper **Correcting batch effects in single-cell RNA sequencing data by matching mutual nearest neighbours** by [Haghverdi _et al. (2018)_](https://doi.org/10.1038/nbt.4091).

## Simulations

To run the simulations, enter the `simulations/` directory and run:

- `cluster_sim.R`, which simulates a variety of scenarios involving orthogonal batch effects.
- `nonorth_sim.R`, which simulates some pathological non-orthogonal batch effects.

## Real data

Three real data analyses are available - `haematopoiesis`, `pancreas` and `droplet`.
Each subdirectory will usually contain:

- `download_data.sh`, to download the data locally.
- `prepareData.R`, to pre-process and normalize the data.
- `plotCorrections.R`, to perform the batch correction and visualize the result with t-SNE plots.

