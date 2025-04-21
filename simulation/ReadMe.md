For an overview of the structure of this directory and a comprehensive demo see `simulation_demo.qmd`.

The main entry point for the simulation is `simulation.R`. Running `source("simulation.R")` will run the current draft version of the simulation, save results, and save plots to a `results/` subdirectory. The simulation study is automatically logged to a file in a `logs/` subdirectory.

Data is generated according to various settings of choice in `generate_data.R`. 

All stan-dependent models are implemented in `BDML_simulation_source.R`, i.e.
- BDML-Basic
- BDML-Hier
- BDML-IW-Hier
all Gibbs-sampling based models are implemented in `BLRs_simulation_source.R`.
- Naive ridge regression
- Hahn (2018)
- Linero (2020)
- FDML-Split and FDML-Full
- BDML-IW

The `test.R` script runs a minimal test of the simulation across all models.

The `bayesm` package dependency is required for Gibbs sampling of the `BDML-IW` model. I found the CRAN version is broken, so until this is fixed, you need to instead install my fork from source: https://github.com/lasse-meixner/bayesm

To do so, clone the repository, navigate to the root directory of the repository, set the working directory to the path of the repository, and run the following command in R:

```R
devtools::install(getwd())
```

This might ask you to also install newer versions of rcpp dependencies like RcppArmadillo from source (this requires a gfortran compiler!). If you already have pre-compiled versions of these packages installed (e.g. from installing the CRAN version of `bayesm`), you can say no to this question, and it should still work.