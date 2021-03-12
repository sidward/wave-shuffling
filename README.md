# Wave-Shuffling

Demonstration code for the MRM manuscript, **Wave-encoding and Shuffling Enables Rapid Time Resolved Structural Imaging**.

The data presented in this repository is from a highly accelerated Wave-Shuffling MPRAGE acquisition that was acquired in 1 minute and 21 seconds.

Written by Siddharth Srinivasan. Please feel free to post a GitHub issue if there's any question I can help answer.

Tested with:
  - MATLAB 2020b
  - BART V7 (``https://zenodo.org/record/4570601``)

## Requirements.

This repository requires BART to be installed. In particular, ``script.m`` requires ``system('bart <command>')`` to work from within MATLAB.
The latest version of BART can be grabbed from ``https://github.com/mrirecon/bart/``.

## Organization

- ``script.m``: Main demo script.
- ``lib/``: Folder containing MATLAB helper functions.
- ``data/``: Folder containing data download script. This is also used to store intermediate arrays required by the reconstruction.

## Acknowledgments

- The wave-encoding "phase per centimeter" code was original written by Kawin Setsompop, and was modified as needed.
- The helper functions for reading and writing complex float (cfl) files were copied from the BART repository.

All rights/distribution are the same as for the original code, and should cite the original author and webpage.

## DOI

TODO
