
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bifurcationEWS

bifurcationEWS is an R package for running an Early Warning Signal
analysis to anticipate deterministic bifurcations in dynamical systems.
The package supports timeseries simulation of many deterministic
bifurcations, computation of some univariate and multivariate EWS, and
an assessment of the performance of EWS.

<!-- badges: start -->
<!-- badges: end -->

## Installation

You can install the development version of bifurcationEWS from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("KCEvers/bifurcationEWS")
#> Downloading GitHub repo KCEvers/bifurcationEWS@HEAD
#> rlang  (1.1.2  -> 1.1.3 ) [CRAN]
#> glue   (1.6.2  -> 1.7.0 ) [CRAN]
#> digest (0.6.33 -> 0.6.34) [CRAN]
#> Rcpp   (1.0.11 -> 1.0.12) [CRAN]
#> scico  (1.3.1  -> 1.5.0 ) [CRAN]
#> plotly (4.10.3 -> 4.10.4) [CRAN]
#> Installing 6 packages: rlang, glue, digest, Rcpp, scico, plotly
#> Installing packages into 'C:/Users/kever/AppData/Local/Temp/RtmpGswm5F/temp_libpath2ea4cc9e89'
#> (as 'lib' is unspecified)
#> package 'rlang' successfully unpacked and MD5 sums checked
#> package 'glue' successfully unpacked and MD5 sums checked
#> package 'digest' successfully unpacked and MD5 sums checked
#> package 'Rcpp' successfully unpacked and MD5 sums checked
#> package 'scico' successfully unpacked and MD5 sums checked
#> package 'plotly' successfully unpacked and MD5 sums checked
#> 
#> The downloaded binary packages are in
#>  C:\Users\kever\AppData\Local\Temp\Rtmpk3DUeB\downloaded_packages
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#>          checking for file 'C:\Users\kever\AppData\Local\Temp\Rtmpk3DUeB\remotes37404e702797\KCEvers-bifurcationEWS-0f7cc4b/DESCRIPTION' ...  ✔  checking for file 'C:\Users\kever\AppData\Local\Temp\Rtmpk3DUeB\remotes37404e702797\KCEvers-bifurcationEWS-0f7cc4b/DESCRIPTION'
#>       ─  preparing 'bifurcationEWS': (1.3s)
#>    checking DESCRIPTION meta-information ...     checking DESCRIPTION meta-information ...   ✔  checking DESCRIPTION meta-information
#>       ─  checking for LF line-endings in source and make files and shell scripts (424ms)
#>       ─  checking for empty or unneeded directories
#>       ─  building 'bifurcationEWS_0.0.0.9000.tar.gz'
#>      
#> 
#> Installing package into 'C:/Users/kever/AppData/Local/Temp/RtmpGswm5F/temp_libpath2ea4cc9e89'
#> (as 'lib' is unspecified)
library(bifurcationEWS)
```

## Generate timeseries

Let’s generate some timeseries for a deterministic Hopf bifurcation:

``` r
# Define control parameter
bifpar_pars = list(
        bifpar_start = .75,
        bifpar_end = .85,
        pre_steps = 0,
        baseline_steps = 0,
        transition_steps = 51,
        post_steps = 0
      )
# Generate Hopf bifurcation timeseries
GLV_Hopf = bifurcation_ts(model = GLV_model,
                         bifpar_pars = bifpar_pars,
                         nr_timesteps = 1000
                        )
#> Specify either a named initial condition X0 (e.g. X0 = (X1 = .1, X2 = .3)) or the names of your variables to set a random initial condition using X_names (e.g. X_names = c('X1', 'X2')).
df = GLV_Hopf$df
```

``` r

# regimes_Hopf = find_regimes(GLV_Hopf_trans, factor_k = 1)
# plot_bifdiag(regimes_Hopf, sel_variables = "X1")
```

<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. -->
<!-- You can also embed plots, for example: -->
<!-- ```{r pressure, echo = FALSE} -->
<!-- plot(pressure) -->
<!-- ``` -->
<!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN. -->
