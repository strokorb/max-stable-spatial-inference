# Max-stable processes for spatial extremes

This is the repository for the code related to the book chapter 

Kirstin Strokorb. Marco Oesting. *Max-stable processes for spatial extremes* \
Chapter 15 of Handbook on Statistics of Extremes (Wiley) - forthcoming 2024+

It contains the following six files:
- *Example.R*
- *Example_Advanced.R*
- *Algorithms2.R*
- *AuxiliaryPlots.R*
- *Chapter15.RData*
- *inlandBRsimu_30.RData*

The analysis constitutes an indicative example for the usage of max-stable processes in spatial extremes for a hypothetical question involving Dutch temperature extremes. 

## Third-Party Code

Overall the analysis follows similar steps as in 

[https://github.com/Oesting/Comparative-tour-through-simulation-algorithms-for-max-stable-processes](https://github.com/Oesting/Comparative-tour-through-simulation-algorithms-for-max-stable-processes)

It uses the same data (see below) and - up to minor modifications - identical auxiliary plot functions and simulation algorithms for max-stable processes.
The two main differences are that

- the fitting procecure in this repository is based on pairwise composite likelihood inference as implemented in [https://spatialextremes.r-forge.r-project.org/](https://spatialextremes.r-forge.r-project.org/) instead of a two-step procedure that separates marginal GEV fitting from M-Estimator based dependence feature fitting,
- the simulation procedure bases Gaussian process simulation now on an incomplete Cholesky decomposition as implemented in [https://cran.r-project.org/package=kernlab](https://cran.r-project.org/package=kernlab) instead of the standard Cholesky decomposition.
  
Some auxiliary functions for the simulation are taken from and parts of simulation algorithms themselves are based on the supplementary material of

C. Dombry, S. Engelke & M. Oesting (2016), Exact simulation of max-stable processes, *Biometrika* 103(2), pp.303-317 

available from [academic.oup.com](https://doi.org/10.1093/biomet/asw008).

## Data

The file *Chapter15.RData* contains daily maximum temperatures from 1990 to 2019 that were measured at 18 inland stations in the Netherlands and are freely available from [knmi.nl](http://projects.knmi.nl/klimatologie/daggegevens/selectie.cgi).



