#' Reduce size of dataframe by downsampling
#'
#' @param df Dataframe
#' @param X_names Names of columns in df containing observational timeseries to be downsampled
#' @param win_size Number of data points to be averaged or sampled from in each step
#' @param which_X String of time point in window to select
#' @param seed_nr Seed number for random sampling of data point per window
#'
#' @return Downsampled dataframe
#' @importFrom dplyr .data
#' @export
downsample <- function(df,
                           X_names,
                           win_size = 10,
                           which_X = c("all", "first", "middle", "last", "random")[1],
                           seed_nr = NULL) {
  if (is.numeric(seed_nr)) {
    set.seed(seed_nr)
  }
  df = df %>% as.data.frame()

  if (which_X == "all") {
    X_day_mu = cbind(
      apply(matrix(df$time, ncol=win_size, byrow = T), 1, max),
      lapply(X_names, function(i){
        Matrix::rowMeans(matrix(df[[i]], ncol=win_size, byrow = T))
      }) %>% do.call(cbind, .),
      apply(matrix(df$bifpar_idx, ncol=win_size, byrow = T), 1, max)
    ) %>% magrittr::set_colnames(c("time", X_names, "bifpar_idx"))
  }  else if (which_X %in% c("first", "middle", "last")){
    idx = ifelse(which_X == "first", 1, ifelse(which_X == "last", win_size, ceiling(win_size / 2)))
    X_day_mu = cbind(
      matrix(df$time, ncol=win_size, byrow = T)[,idx],
      lapply(X_names, function(i){
        matrix(df[[i]], ncol=win_size, byrow = T)[,idx]
      }) %>% do.call(cbind, .),
      matrix(df$bifpar_idx, ncol=win_size, byrow = T)[,idx]
    )
  } else if (which_X == "random") {
    # Select random sample in each window
    nrow_mat = (nrow(df)/win_size)
    idxs = cbind(row = 1:nrow_mat, col = sample(1:win_size, size = nrow_mat, replace = T))
    X_day_mu = cbind(
      matrix(df$time, ncol=win_size, byrow = T)[idxs],
      lapply(X_names, function(i){
        matrix(df[[i]], ncol=win_size, byrow = T)[idxs]
      }) %>% do.call(cbind, .),
      matrix(df$bifpar_idx, ncol=win_size, byrow = T)[idxs]
    )
  }
  return(X_day_mu)

}


#' Add observational noise
#'
#' @param df Dataframe
#' @param X_names  Names of columns in df to add observational noise to
#' @param noise_mean Mean of white noise
#' @param noise_sigma Standard deviation of white noise
#' @param noise_constant Constant to add to observed timeseries
#' @param seed_nr Seed number for generating noise
#'
#' @return Noisy dataframe
#' @export
add_obs_noise <-
  function(df,
           X_names,
           noise_mean = 0,
           noise_sigma = .01,
           noise_constant = 0,
           seed_nr = NULL) {
    if (!is.null(seed_nr)) {
      set.seed(seed_nr)
    }
    # Add white noise to each variable
    df[, X_names] = df[, X_names] + pracma::Reshape(
      stats::rnorm(
        n = length(X_names) * nrow(df),
        mean = noise_mean,
        sd = noise_sigma
      ),
      n = nrow(df),
      m = length(X_names)
    ) + noise_constant

    return(df)
  }





#' Modify list
#'
#' @param old_list List to modify
#' @param new_list List to modify old_list with
#'
#' @return Modified list
#' @export
modify_list <- function(old_list, new_list){
  # utils::modifyList() doesn't update nested list as I want them to: if a nested list has a different number of elements, it doesn't seem to be replaced
  utils::modifyList(old_list[setdiff(names(old_list), names(new_list))], new_list)

}



#' Set up parameters
#'
#' @param model_name Chosen dynamical systems model
#' @param pars_add List of parameters to overwrite or add to default parameters
#'
#' @return List of parameters
#' @export
setup_pars <- function(model_name, pars_add = list()) {
  if (model_name == "detGLV") {
    p = 4

    pars_default = list(
      mainDir = getwd(),
      model_name = model_name,
      model = GLV_model,
      p = p,
      X_names = paste0("X", 1:p),
      X0 = c(),
      timestep = .01,
      fs = 10,
      nr_timesteps = 1000,
      model_pars = list(
        s = .85,
        mu = rep(0, p),
        r = c(1, .72, 1.53, 1.27),
        C0 = matrix(
          c(
            1,
            1.09,
            1.52,
            0,
            0,
            1,
            0.44,
            1.36,
            2.33,
            0,
            1,
            0.47,
            1.21,
            0.51,
            0.35,
            1
          ),
          p,
          p,
          byrow = TRUE
        )
      ),
      seed = 123,
      max_iter = 100,
      deSolve_method = "rk4",
      downsample_pars = list(
        type = c("average", "one_sample")[1],
        win_size = 50,
        which_X = c(50, "first", "middle", "last", "random")[1],
        seed = 123
      ),
      stopifregime = get_stopifregime("default"),
      bifpar_pars = list(
        bifpar_start = .96,
        bifpar_end = .97,
        pre_steps = 0,
        baseline_steps = 0,
        transition_steps = 101,
        post_steps = 0
      ),
      X_sigma = 0,
      noise_constant = 0,
      nr_noise_iters = 10,
      downsample_fs = c(10, 1, .1),
      sigma_obs_noise = c(0.0001, .02, .04),
      sigma_crit_step = .25,
      nr_consecutive_warnings = 1
    )

    pars = modify_list(
      modify_list(get_formals(bifurcation_ts), pars_default), pars_add)
    pars$bifpar_list = do.call(get_bifurcation_range, pars$bifpar_pars)
    return(pars)

  } else {
    return("Model not implemented yet!")
  }
}


#' Get formal arguments of function as list
#'
#' @param func Function
#'
#' @return List with formal arguments and default arguments
#' @export
get_formals <- function(func){
  formal_list = purrr::map(as.list(formals(func)),
             function(x){try(eval(x), silent = T)
               })
  return(formal_list[!unlist(purrr::map(formal_list,
                                        function(x){"try-error" %in% class(x)}))])
}


#' Specify parameters for specific regime switches
#'
#' @return List with parameters per regime switch
#' @export
setup_bifpars <- function() {

  pars_template_adjust = list(
    list(
      regime_switch_name = "Saddle-Node",
      data_idxs = 1,
      factor_k = 1,
      nr_timesteps_full = 100,
      nr_noise_iters = 250,
      bifpar_pars = list(
        bifpar_start = 1.15,
        bifpar_end = .9,
        pre_steps = 0,
        baseline_steps = 0,
        transition_steps = 2501,
        post_steps = 0
      ),
      select_regime_switches = c("PD_1to1"),
      GLV_string = c("GLV_SN"),
      stopifregime = get_stopifregime(type = "only_nodes")
    ),
    list(
      regime_switch_name = "Hopf",
      data_idxs = 1,
      factor_k = 1,
      nr_timesteps_full = 200,
      nr_noise_iters = 250,
      bifpar_pars = list(
        bifpar_start = .6,
        bifpar_end = .85,
        pre_steps = 0,
        baseline_steps = 0,
        transition_steps = 251,
        post_steps = 0
      ),
      select_regime_switches = c("PD_1to2"),
      GLV_string = c("GLV_Hopf")
    ),
    list(
      regime_switch_name = "PD_to_Chaos",
      min_length_regime_full = 8,
      bifpar_pars = list(
        bifpar_start = .85,
        bifpar_end = .985,
        pre_steps = 0,
        baseline_steps = 0,
        transition_steps = 1351,
        post_steps = 0
      ),
      select_regime_switches = c(
        "PD_2to4",
        "PD_4to8",
        "PD_8to16",
        "PD_Mixed-Periodic_to_Chaotic1",
        "SUBD_Chaotic_to_Mixed-Periodic1"
      ),
      GLV_string = c(
        "GLV_PD_2to4",
        "GLV_PD_4to8",
        "GLV_PD_8to16",
        "GLV_PD_Chaos",
        "GLV_SUBD_forw")
    ),
    list(
      regime_switch_name = "PH_from_Chaos",
      thresh_coord_spread_trans = .01,
      factor_k = .01,
      min_length_regime_trans = 5,
      bifpar_pars = list(
        bifpar_start = .985,
        bifpar_end = .8,
        pre_steps = 0,
        baseline_steps = 0,
        transition_steps = 1851,
        post_steps = 0
      ),
      select_regime_switches = c(
        "PH_16to8",
        "PH_8to4",
        "PH_4to2",
        "PH_2to1",
        "PH_Chaotic_to_Mixed-Periodic1",
        "SUBD_Mixed-Periodic_to_Chaotic1"
      ),
      GLV_string = c(
        "GLV_PH_16to8",
        "GLV_PH_8to4",
        "GLV_PH_4to2",
        "GLV_PH_2to1",
        "GLV_PH_Chaos",
        "GLV_SUBD_backw"
      ),
      stopifregime = get_stopifregime(type = "no_nodes")
    ),
    list(
      regime_switch_name = "Interior-Crisis-Merging",
      nr_timesteps_full = 2000,
      thresh_full_band_trans = .6,
      nr_smooth_full = 5,
      nr_smooth_trans = 0,
      min_length_regime_full = 5,
      min_length_regime_trans = 1,
      bifpar_pars = list(
        bifpar_start = 1.01,
        bifpar_end = 1.02,
        pre_steps = 0,
        baseline_steps = 0,
        transition_steps = 201,
        post_steps = 0
      ),
      select_regime_switches = c("Interior-Crisis-Merging"),
      GLV_string = c("GLV_IC_forw"),
      stopifregime = get_stopifregime(type = "no_nodes")
    ),
    list(
      regime_switch_name = "Interior-Crisis-Separation",
      nr_timesteps_full = 2000,
      factor_k = 1,
      nr_smooth_full = 10,
      nr_smooth_trans = 5,
      thresh_coord_spread_full = .01,
      thresh_full_band_trans = .6,
      bifpar_pars = list(
        bifpar_start = 1.03,
        bifpar_end = 1.01,
        pre_steps = 0,
        baseline_steps = 0,
        transition_steps = 201,
        post_steps = 0
      ),
      select_regime_switches = c("Interior-Crisis-Separation"),
      GLV_string = c("GLV_IC_backw"),
      stopifregime = get_stopifregime(type = "no_nodes")
    ),
    list(
      regime_switch_name = "Boundary-Crisis",
      factor_k = 1,
      deSolve_method = "euler",
      timestep = .1,
      nr_timesteps_full = 10000,
      keep_nr_timesteps = 1000,
      bifpar_pars = list(
        bifpar_start = 1.1,
        bifpar_end = 1.3,
        pre_steps = 0,
        baseline_steps = 0,
        transition_steps = 201,
        post_steps = 0
      ),
      select_regime_switches = c("Boundary-Crisis"),
      GLV_string = c("GLV_BC"),
      stopifregime = get_stopifregime(type = "no_nodes"),
      min_length_regime_full = 20,
      min_length_regime_trans = 20,
      variable_name = "X4",
      thresh_full_band_trans = 1,
      thresh_coord_spread_full = .1,
      thresh_coord_spread_trans = .1,
      nr_smooth_trans = 0
    )
  )  %>%
    stats::setNames(unlist(purrr::map(., "regime_switch_name"))) %>%
    # Add default arguments of bifurcation_ts() and find_regimes(), add own parameters
    purrr::map(., function(x){modify_list(get_formals(find_regimes), x)}) %>%
    plyr::llply(., function(x){modify_list(x,
                                     list(s_string = sprintf("s%.05f-s%.05f-by%.05f", x$bifpar_pars$bifpar_start, x$bifpar_pars$bifpar_end, c((x$bifpar_pars$bifpar_end-x$bifpar_pars$bifpar_start)/(x$bifpar_pars$transition_steps-1)) ),
                                          bifpar_list = do.call(get_bifurcation_range, x$bifpar_pars)))})

  return(pars_template_adjust)

}

#' Format and recursively create file path
#'
#' @param pars_file List of names of directories and filename components needed to create a directory
#'
#' @return File path
#' @export
format_path <- function(pars_file) {
  pars_file_default = list(
    mainDir = getwd(),
    model_name = "model",
    type_output = "",
    analysis_type = "",
    subfolder1 = "",
    subfolder2 = "",
    subfolder3 = "",
    subfolder4 = "",
    subfolder5 = "",
    filename = "",
    file_ID = "",
    pars_ID = "",
    file_ext = ".RDS"
  )

  pars_file = modify_list(pars_file_default, pars_file)

  # Format parameters correctly
  if (!is.null(pars_file$file_ext)) {
    file_ext <- pars_file$file_ext
    # Add a leading period to file extension
    pars_file$file_ext <-
      if (stringr::str_sub(file_ext, 1, 1) != ".")
        sprintf(".%s", file_ext)
    else
      pars_file$file_ext # Make sure file extension leads with a period
  }


  # Define file paths with information from pars_file
  with(pars_file,
       {
         filepath_base <- file.path(mainDir, model_name, pars_ID)
         filepath_dir <-
           filepath_base %>% file.path(
             type_output,
             analysis_type,
             subfolder1,
             subfolder2,
             subfolder3,
             subfolder4,
             subfolder5
           ) %>%
           normalizePath(mustWork = FALSE)

         # Create directory
         dir.create(filepath_dir,
                    showWarnings = FALSE,
                    recursive = TRUE)

         # Update filename
         filename_ = stringr::str_c(file_ID, filename, sep = "_")

         # Remove leading and trailing separators
         filename = stringr::str_c(
           filename_ %>% stringr::str_replace("^_", "") %>%
             stringr::str_replace("_$", "") %>%
             stringr::str_replace("__", "_"),
           file_ext
         ) # Add file extension

         filepath_final <-
           filepath_dir %>% file.path(filename) %>% normalizePath(mustWork = FALSE)
         return(filepath_final)
       })
}



#' Format parameters
#'
#' @param pars List of parameters
#'
#' @return Updated list of parameters
#' @export
format_pars <- function(pars) {
  pars_add = with(pars, {
    if (!is.null(pars[["data_idx"]])) {
      file_ID = sprintf(
        "nr%d_%s_T%d_ts%.3f_fs%.2f_Xsigma%.5f",
        data_idx,
        deSolve_method,
        nr_timesteps,
        timestep,
        fs,
        X_sigma
      )
    } else {
      file_ID = sprintf(
        "%s_T%d_ts%.3f_fs%.2f_Xsigma%.5f",
        deSolve_method,
        nr_timesteps,
        timestep,
        fs,
        X_sigma
      )

    }
    # pars_ID

    return(list(file_ID = file_ID))
  })

  return(modify_list(pars, pars_add))
}



#' Solve for sampling parameters
#'
#' @param fs Sampling frequency
#' @param sample_interval Sampling interval
#' @param timestep Timestep
#'
#' @return List of sampling parameters
#' @export
solve_sampling_par <-
  function(fs = NULL,
           sample_interval = NULL,
           timestep = NULL) {
    if (sum(is.null(fs) &
            is.null(sample_interval) & is.null(timestep)) > 1) {
      message("Need at least two parameters!")
      return()
    }
    if (is.null(fs)) {
      fs = (1 / sample_interval) * (1 / timestep)
    } else if (is.null(sample_interval)) {
      sample_interval = (1 / timestep) / fs
    } else if (is.null(timestep)) {
      timestep = 1 / (fs * sample_interval)
    }
    return(list(
      fs = fs,
      sample_interval = sample_interval,
      timestep = timestep
    ))
  }

#' Utility function for foreach loop to save memory
#'
#' @param ... All outputs from foreach loop
#'
#' @return NULL
#' @export
cfun <- function(...) {
  NULL
}

#' Generate forloop for foreach loop
#'
#' @param ... Named forloop parameters
#'
#' @return List of forloop parameters
#' @export
get_forloop <- function(...) {
  return(tidyr::expand_grid(...) %>% purrr::transpose() %>% unique())
}


#' Add parameters in list as columns in dataframe
#'
#' @param df Dataframe
#' @param for_par List of parameters
#'
#' @return Updated dataframe
#' @export
add_par_as_cols = function(df, for_par) {

  # Remove any entries that are already in df
  df2 = dplyr::bind_cols(df, for_par[setdiff(names(for_par), names(df))])
  return(df2)
}


#' Download example data sets
#'
#' @return Names of downloaded R objects
#' @export
download_example_data <- function(){

  # Get files from OSF
  anes_project <- osfr::osf_retrieve_node("https://osf.io/85d4h")
  anes_files <- osfr::osf_ls_files(anes_project, n_max = Inf)

  input <- readline(prompt="The complete example data set is about 5.44GB and will take quite a while to download. Would you like to download all files (type 1), some files (type 2), or none (type 3)?")
  if (input == 1){
    osfr::osf_download(anes_files, verbose = T, progress = T, conflicts = "skip")
    print("Download complete! The data can be found in the following R objects:")
    return(stringr::str_replace(anes_files$name, ".RDS", ""))
  } else if (input == 2){
    print(as.data.frame(anes_files["name"]))
    input2 <- readline(prompt="Printed above you'll find the available files. Please type the row numbers of the files you want to download like so: 1,2,3:")
    input2 = as.numeric(stringr::str_split(input2, ",")[[1]])
    if (all(input2 %in% 1:nrow(anes_files))){
      t = osfr::osf_download(anes_files[input2,], verbose = T, progress = T, conflicts = "skip")
      print("Download complete! The data can be found in the following R objects:")
      return(stringr::str_replace(anes_files$name[input2], ".RDS", ""))
    } else {
      print("Your input does not correspond to row indices in the dataframe.")
      return()
    }
  } else {
    return()
  }
}

# Minimal outline building R package
#
# Create .R script
# usethis::use_r("example")
#
# Clear workspace
# devtools::load_all(".")
#
# devtools::build_vignettes()
#
# devtools::document()
# ch = devtools::check()
# ch = devtools::check(vignettes=FALSE)
# devtools::install()
# or:
# devtools::install(build_vignettes = TRUE)

# browseURL("doc/demo.html") # View vignette

# Add dependencies
# pkgs = c("deSolve", "moments", "utils", "tidyr", "ggplot2", "purrr", "casnet", "rlang", "pracma", "dplyr","rgl", "ggh4x", "stats", "stringr", "zoo", "cowplot", 'magrittr', 'scales', "viridis", "grDevices", "tools", "grid", "invctr", "Matrix", "gsignal", "plyr", "plotly", "tibble", "latex2exp", "foreach", "scico", "osfr")
# for (p in pkgs){
# usethis::use_package(p)
# }
# Allow for use of pipes :
# usethis::use_pipe()

# Add data
# https://grasshoppermouse.github.io/posts/2017-10-18-put-your-data-in-an-r-package/
# # # Set up the data-raw directory and data processing script
# # # You can use any name you want for your data
# data_file = "GLV_Hopf_trans"
# usethis::use_data_raw(data_file)
# # # This script in the R directory will contain the documentation.
# # # You can use any name you want.
# file.create(sprintf("R/%s.R", data_file))
# # # Put your data files into the data-raw folder. Your new package directory should look something like this:
# # # Write your data processing code in a data-raw/mydataset.R script. It would look something like this:
# #   # data-raw/mydataset.R
# #   # Data import and processing pipeline
# #   library(readr)
# # library(readxl)
# # mydataset <- read_csv("data-raw/pendulum data.csv")
# # demographics <- read_excel("data-raw/Demographics.xlsx")
# # # Data cleaning code here...
# # # (Do NOT put data analysis code here!)
# # # This should be the last line.
# # # Note that names are unquoted.
# # # I like using overwrite = T so everytime I run the script the
# # # updated objects are saved, but the default is overwrite = F
# # usethis::use_data(GLV_Hopf_trans, overwrite = T)
# sinew::makeOxygen("GLV_Hopf_trans", add_fields = "source")
# sinew::makeOxygen("GLV_Hopf_trans")
# # tools::resaveRdaFiles(sprintf("data/%s.rda", data_file),compress="xz")
# # tools::checkRdaFiles("data/")# To check that the compression has been done successfully
# # GIven a LazyDataCompression warning: Add
# # LazyDataCompression:xz
# # to the description file.

# Citation file
# knitr::write_bib(c(.packages(), "bifurcationEWS"), "packages.bib")

# Readme file
# usethis::use_readme_rmd()
# devtools::build_readme() # when finished

