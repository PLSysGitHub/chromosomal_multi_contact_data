# Physical models for chromosome organization to predict multi-contact statistics
 
## Summary of contents:
1. Data files used for paper, as well as scripts for extracting contact and contact triplet frequencies
2. Functions for creating predictions for contact triplet frequencies, given contact maps (in "ContactTripletPredictions.jl")
3. Scripts for making plots shown in the paper


## Requirements:
- Julia (tested using 1.6)
- Following packages:
  - HDF5 (for help with installation, see https://juliaio.github.io/HDF5.jl/stable/)
  - DelimitedFiles
  - HypothesisTests
  - MultipleTesting
  - Plots
  - LinearAlgebra
  - StatsBase

## Usage:
- All scripts are written presuming that your working directory is this top directory
- Running "./Extract_triplets/Extract_triplets_XXXX.jl" will use the provided raw data files to create both a contact frequency and a contact triplet frequency file. These are saved in "./Contact_files/" and "./Triplet_files/" respectively
	-The large triplet files are not included, but can be generated with the provided code or downloaded from [Zenodo](https://doi.org/10.5281/zenodo.14065130) 

- Code in "./Extract_triplets/Extract_triplets_Brandao_et_al.jl" can be edited to also extract the contact map from non-interacting simulations
- The script "./RandomWalkForP0.jl" can be used to create a contact map for a non-interacting polymer roughly corresponding to the chromosome segment of the Bintu et. al. data set
- "./Example_XXXX.jl" scripts can be used to calculate predictions for each contact triplet data set, and to save plots like shown in the paper in the directory "./Output/XXXX_example/"

## Questions?
- For questions, you can email j.k.harju [at] vu.nl
