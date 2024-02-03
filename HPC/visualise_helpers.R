# AUC
main_color= colorspace::sequential_hcl(palette='Greens 3', n =10)[1]
col_facet_labels = scales::viridis_pal(option = "rocket", direction =
                                         -1)(20)[17]

regimes_switch_labels = list(
  "Period-Doubling" = list(
    "PD_2to4" = "2 to 4",
    "PD_4to8" = "4 to 8",
    "PD_8to16" = "8 to 16"
  ),
  "Period-Halving"= list(
    "PH_2to1" = "Hopf (Backward)",
    "PH_4to2"= "4 to 2",
    "PH_8to4"= "8 to 4",
    "PH_16to8"= "16 to 8"
  ),
  "Chaotic" = list(
    "PD_Mixed-Periodic_to_Chaotic1" = "Period-Doubling Cascade", # with a subduction wedged in between, to "Period-6 (x2,x3,x4) AND Period-8 (x1)"
    "PH_Chaotic_to_Mixed-Periodic1" = "Period-Halving Cascade", # with a subduction wedged in between, to "Period-6 (x2,x3,x4) AND Period-8 (x1)"
    "SUBD_Mixed-Periodic_to_Chaotic1" = "Subduction: Periodic to Chaotic", # see no change in advance of bifurcation #"Period-6 (x1,x2,x3,x4)"
    "SUBD_Chaotic_to_Mixed-Periodic1" = "Subduction: Chaotic to Periodic",
    "Boundary-Crisis" = "Boundary Crisis",
    "Interior-Crisis-Merging" = "Expansion Chaos",
    "Interior-Crisis-Separation" = "Reduction Chaos"
  ),
  "Fixed-Point" = list(
    "PD_1to2" = "Hopf",
    "PD_1to1" = "Saddle-Node"
  )
)

regime_switch_to_class = plyr::ldply(1:length(regimes_switch_labels), function(idx){
  return(data.frame(regimes_switch = names(regimes_switch_labels[[idx]])) %>% tibble::add_column(regimes_switch_class = names(regimes_switch_labels[idx])))
}) %>% tibble::deframe()

metric_labels = list(
  "Generic" = list(
    "mean_var1" = latex2exp::TeX("Mean ($x_1$)", output = "character"),
    "mean_var2" = latex2exp::TeX("Mean ($x_2$)", output = "character"),
    "mean_var3" = latex2exp::TeX("Mean ($x_3$)", output = "character"),
    "mean_var4" = latex2exp::TeX("Mean ($x_4$)", output = "character"),
    "variance_var1" = latex2exp::TeX("Variance ($x_1$)", output = "character"),
    "variance_var2" = latex2exp::TeX("Variance ($x_2$)", output = "character"),
    "variance_var3" = latex2exp::TeX("Variance ($x_3$)", output = "character"),
    "variance_var4" = latex2exp::TeX("Variance ($x_4$)", output = "character"),
    "COV_var1" = latex2exp::TeX("COV ($x_1$)", output = "character"),
    "COV_var2" = latex2exp::TeX("COV ($x_2$)", output = "character"),
    "COV_var3" = latex2exp::TeX("COV ($x_3$)", output = "character"),
    "COV_var4" = latex2exp::TeX("COV ($x_4$)", output = "character"),
    "autocorrelation_var1" = latex2exp::TeX("Autocorrelation ($x_1$)", output = "character"),
    "autocorrelation_var2" = latex2exp::TeX("Autocorrelation ($x_2$)", output = "character"),
    "autocorrelation_var3" = latex2exp::TeX("Autocorrelation ($x_3$)", output = "character"),
    "autocorrelation_var4" = latex2exp::TeX("Autocorrelation ($x_4$)", output = "character"),
    "skewness_var1" = latex2exp::TeX("Skewness ($x_1$)", output = "character"),
    "skewness_var2" = latex2exp::TeX("Skewness ($x_2$)", output = "character"),
    "skewness_var3" = latex2exp::TeX("Skewness ($x_3$)", output = "character"),
    "skewness_var4" = latex2exp::TeX("Skewness ($x_4$)", output = "character"),
    "kurtosis_var1" = latex2exp::TeX("Kurtosis ($x_1$)", output = "character"),
    "kurtosis_var2" = latex2exp::TeX("Kurtosis ($x_2$)", output = "character"),
    "kurtosis_var3" = latex2exp::TeX("Kurtosis ($x_3$)", output = "character"),
    "kurtosis_var4" = latex2exp::TeX("Kurtosis ($x_4$)", output = "character")
  ),
  "Multivariate" = list("meanAbsCrossCorr" = latex2exp::TeX("Mean Abs. Cross-Corr.", output = "character"),
                        "largestLambdaCovMatrix" = latex2exp::TeX("First Eigenvalue of $Sigma$", output = "character"),
                        "spatial_variance" = latex2exp::TeX("Spatial Variance", output = "character"),
                        "spatial_skewness" = latex2exp::TeX("Spatial Skewness", output = "character"),
                        "spatial_kurtosis" = latex2exp::TeX("Spatial Kurtosis", output = "character")),
  "Spectral" = list(
    "Smax_X1" = latex2exp::TeX("Max. Spectral Density ($x_1$)", output = "character"),
    "Smax_X2" = latex2exp::TeX("Max. Spectral Density ($x_2$)", output = "character"),
    "Smax_X3" = latex2exp::TeX("Max. Spectral Density ($x_3$)", output = "character"),
    "Smax_X4" = latex2exp::TeX("Max. Spectral Density ($x_4$)", output = "character"),
    "spectral_exp_var1" = latex2exp::TeX("Spectral Exponent ($x_1$)", output = "character"),
    "spectral_exp_var2" = latex2exp::TeX("Spectral Exponent ($x_2$)", output = "character"),
    "spectral_exp_var3" = latex2exp::TeX("Spectral Exponent ($x_3$)", output = "character"),
    "spectral_exp_var4" = latex2exp::TeX("Spectral Exponent ($x_4$)", output = "character"),
    "spectral_ratio_var1" = latex2exp::TeX("Spectral Ratio ($x_1$)", output = "character"),
    "spectral_ratio_var2" = latex2exp::TeX("Spectral Ratio ($x_2$)", output = "character"),
    "spectral_ratio_var3" = latex2exp::TeX("Spectral Ratio ($x_3$)", output = "character"),
    "spectral_ratio_var4" = latex2exp::TeX("Spectral Ratio ($x_4$)", output = "character")
  )
)

metric_to_class = plyr::ldply(1:length(metric_labels), function(idx){
  return(data.frame(metric = names(metric_labels[[idx]])) %>% tibble::add_column(metric_class = names(metric_labels[idx])))
}) %>% tibble::deframe()
