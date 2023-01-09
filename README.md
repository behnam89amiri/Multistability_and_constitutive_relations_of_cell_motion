[![License: CC BY-NC-ND 4.0](https://img.shields.io/badge/License-CC_BY--NC--ND_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-nd/4.0/)

# Multistability of the motile state of MDA MB 231 cells
This code is used to study the multistability and constitutive relations of single cell motion on Fibronectin lanes.
For an automated high-throughput data analysis, the fluorescent signal of the nucleus of single cells is tracked over time, see `track_cells`.
To analyse the motion of cells in more detail, the timelapse movies of the cells are converted to kymographs and the position of the front, the nucleus and the back edge of the cells can be extracted in a semi-automatic fashion, see `kymographs`.
After cells are tracked, the trajectories can be analysed by determining their characteristic parameters, see `trajectory_analysis`.
Cell trajectories can be simulated successfully using our bio-physical model, see `simulations`.


## Structure
### track_cells
This directory contains all scripts necessary to track the nuclear movement. Mainly MATLAB files.

### kymographs
Macros to create kymographs using ImageJ / FIJI and to manually segment the kymographs. Mainly ImageJ macro language.

### trajectory_analysis
Scripts that plot the single cell trajectories, identify change points, classify episodes of the trajectory into a motile state and analyse the motile behaviour at certain events, e.g. direction reversals. All written in MATLAB.

### simulations
Scripts that simulate single cell trajectories based on our bio-physical understanding of the system. The resulting simulated trajectories can be compared to experimental trajectories using scripts in the `trajectory_analysis` directory.
