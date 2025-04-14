# eitDMD
Collection of dynamic mode decomposition (DMD) algorithms implemented in an object-oriented manner for application in electrical impedance tomography (EIT).

Implemented DMD algorithms:
* Vanilla and Exact DMD
* Batch (sliding window) and online (sample by sample with update Koopman operator) DMD
* Multi-resolution DMD
* Dimensionality reduction of input using compressed sensing
* Continuous and discrete modes
* Advanced mode visualization and selection from [Visualization and selection of Dynamic Mode Decomposition components for unsteady flow](https://www.sciencedirect.com/science/article/pii/S2468502X21000309)

Main DMD class found in ./source/.

# Literature and links for previous work
Literature used for the development found in ./literature/.

Useful links to previous work:
* https://github.com/cwrowley/dmdtools (Streaming DMD v1 in Python and Matlab)
* https://github.com/jaimeliew1/Streaming-DMD (Streaming DMD v2 in Python)
* https://github.com/MSU-dcypherlab/Incremental-DMD-for-EEG-Data (Online DMD v2; code is very messy – leave for last)
* https://github.com/haozhg/odmd (Online DMD v1 in Python)
* https://github.com/haozhg/odmd-matlab (Online DMD v1 in Matlab)
* https://github.com/VArdulov/online_dmd (an indpendent user’s „interpretation“ of the Online DMD code in Python. It might be that he actually managed to improve on something)
* http://www.ece.umn.edu/users/mihailo//software/dmdsp/ (Sparsity DMD)
