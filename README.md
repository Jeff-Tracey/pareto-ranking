# pareto-ranking



Program written by Jeff A. Tracey, PhD (jeff.a.tracey@proton.me)

R and C/C++ code to accompany the article:  Tracey, Jeff A., et al. "Prioritizing conserved areas threatened by wildfire and fragmentation for monitoring and management." PLoS One 13.9 (2018): e0200203.  https://doi.org/10.1371/journal.pone.0200203



**Code for performing Pareto ranking**



Note: This code has been used ONLY on a Linux operating system. You must have the GNU gcc compilers installed along with R and the Rcpp package.

To Use:
1) Build the shared object (*.so) file by running the script `build_fast_pareto.sh`. You must set the permissions on this file so that you can execute it.
2) The file paretoRankFunctions.R provides R wrappers for the C++ functions.
3) In R, run the `rankExample.R` script as an example. This file is a template for running the code on a matrix containing criteria. In the matrix, each row corresponds to a different alternative and each column corresponds to a different criterion. It is assumed that all criteria values are non-negative and that a lower value indicates a better solution.