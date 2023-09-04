# Minimal outline building R package
#
# Create .R script
# usethis::use_r("example")
#
# Clear workspace
# devtools::load_all(".")
#
# devtools::build_vignettes()
# devtools::document()
# devtools::check()
# devtools::check(vignettes=FALSE)
# devtools::install()
# or:
# devtools::install(build_vignettes = TRUE)

# browseURL("doc/demo.html")

# Add dependencies
# pkgs = c("deSolve", "moments", "utils", "tidyr", "ggplot2", "purrr", "casnet", "rlang", "pracma", "dplyr","rgl", "ggh4x", "stats", "stringr", "zoo")
# for (p in pkgs){
# usethis::use_package(p)
# }


#' Reduce size of dataframe by downsampling
#'
#' @param df Dataframe
#' @param X_names Names of columns in df containing observational timeseries to be downsampled
#' @param type Type of downsampling. "average" averages data points in each window; "one_sample" picks one sample in each window.
#' @param win_size Number of data points to be averaged or sampled from in each step
#' @param which_X Time point(s) in window to select
#' @param seed Seed number for random sampling of data point per window
#'
#' @return Downsampled dataframe.
#' @export
#'
#' @examples
downsample <- function(df, X_names, type = c("average", "one_sample")[1],
                       win_size = 10,
                       which_X = c(10, "first", "middle", "last", "random")[1],
                       seed_nr = 123
){
  set.seed(seed_nr)
  df = df %>% as.data.frame()

  if (type == "average"){
    # Compute average every X samples with possibility of random samples drawn in every X
    slice_func = slice_sample

    which_X = as.numeric(which_X)
    if (!is.na(which_X)){
      n = which_X
    } else {
      return("which_X must be a number!")
    }

    if (n > win_size){
      n = win_size
      print("The number of samples per win_size cannot be more than the win_size itself. Setting the number of samples to the win_size.")
    }
  } else if (type == "one_sample"){
    # Select one sample every win_size with the option of choosing which sample
    slice_func = slice

    if (which_X == "first"){
      n = 1
    } else if (which_X == "middle"){
      n = ceiling(win_size/2)
    } else if (which_X == "last"){
      n = win_size
    } else if (which_X == "random"){
      n = 1
      slice_func = slice_sample
    }
  }

  X_day_mu_ = df %>%
    dplyr::mutate(win_nr = rep(1:ceiling(nrow(df)/win_size), each = win_size)[1:nrow(df)]) %>%
    group_by(win_nr) %>%
    slice_func(n = n)
  X_day_mu = merge(X_day_mu_ %>%
                     dplyr::summarise_at(setdiff(names(.), c("win_nr", X_names)), max),
                   X_day_mu_ %>%
                     dplyr::summarise_at(X_names, mean))
  return(X_day_mu)

}


#' Add observational noise
#'
#' @param df Dataframe
#' @param X_names  Names of columns in df to add observational noise to
#' @param noise_mean Mean of white noise
#' @param noise_sigma Standard deviation of white noise
#' @param noise_constant Constant to add to observed timeseries
#'
#' @return Noisy dataframe
#' @export
#'
#' @examples
add_obs_noise <- function(df, X_names, noise_mean = 0, noise_sigma = .01, noise_constant = 0){

  # Add white noise to each variable
  df[,X_names] = df[,X_names] + pracma::Reshape(rnorm(n = length(X_names) * nrow(df), mean = noise_mean, sd = noise_sigma), n = nrow(df), m = length(X_names)) + noise_constant

  return(df)
}


#' Style plot
#'
#' @param pl Original plot
#' @param col_pal Named colours
#' @param fs Font sizes
#'
#' @return Styled plot
#' @export
#'
#' @examples
style_plot <- function(pl,
                       col_pal = c(
                         col_facet_labels = scales::viridis_pal(option = "rocket", direction =
                                                                  -1)(20)[17]),
                       fs = c(
                         "family" = "serif",
                         "strip.text.x" = 12,
                         "strip.text.y" = 12,
                         "plot.title" = 20,
                         "plot.subtitle" = 16,
                         "axis.text" = 8,
                         "axis.title" = 16,
                         "legend.text" = 12,
                         "legend.title" = 14,
                         "legend.spacing.y" = .075
                       )) {
  pl <- pl + theme_bw() +
    theme(
      text = element_text(family = fs["family"]),
      # Change font
      plot.title = element_text(size = fs["plot.title"]),
      plot.subtitle = element_text(size = fs["plot.subtitle"]),
      axis.text = element_text(size = fs["axis.text"]),
      axis.title = element_text(size = fs["axis.title"]),
      legend.text = element_text(size = fs["legend.text"]),
      legend.title = element_text(size = fs["legend.title"]),
      # # Increase font size facet labels
      strip.text.y = element_text(
        size = as.numeric(fs["strip.text.y"]) + 2,
        margin = margin(0.1, 0.1, 0.1, 0.1, "cm")
      ),
      strip.text.x = element_text(
        size = as.numeric(fs["strip.text.x"]) + 2,
        margin = margin(0.1, 0.1, 0.1, 0.1, "cm")
      )
      # panel.spacing.x = unit(0.2, "cm"),
      # panel.spacing.y = unit(0.2, "cm") # Distance between facets
    ) +
    theme(strip.background = element_rect(fill = col_pal["col_facet_labels"], color = col_pal["col_facet_labels"])) + # Change facet rectangle colour
    theme(strip.text = element_text(colour = 'white')) +
    # theme(legend.position = "none") + # Remove group legend
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "pt")) +
    theme(legend.title.align = 0.5,
          # Text label alignment. Number from 0 (left) to 1 (right)
          legend.text.align = 0.5)

  return(pl)
}




#' Set up parameters
#'
#' @param model_name Chosen dynamical systems model
#'
#' @return List of parameters
#' @export
#'
#' @examples
setup_pars <- function(model_name, pars_add = list()){

  if (model_name == "detGLV"){

    p = 4

    pars_default = list(
      mainDir = getwd(),
      model_name = model_name,
      p = p,
      X_names = paste0("X", 1:p),
      timestep = .01,
      nr_timesteps = 1000,
      model_pars = list(s = .85,
                        mu = rep(0, p),
                        r = c(1, .72, 1.53, 1.27),
                        C0 = matrix(c(1,
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
                                      1),
                                    p,
                                    p,
                                    byrow = TRUE)
      ),
      seed = 123,
      downsample_pars = list(type = c("average", "one_sample")[1],
                             win_size = 50,
                             which_X = c(50, "first", "middle", "last", "random")[1],
                             seed = 123),
      times = seq(0,2000,by=.01),
      bifpar_list = purrr::transpose(expand.grid(s = seq(.96,
                                                         .97,
                                                         by = .0001))),
      X_sigma = 0

    )
    pars = utils::modifyList(pars_default, pars_add)
    return(pars)

  } else {
    return("Model not implemented yet!")
  }
}


#' Format and recursively create file path
#'
#' @param pars_file List of names of directories and filename components needed to create a directory
#'
#' @return File path
#' @export
#'
#' @examples
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

  pars_file = utils::modifyList(pars_file_default, pars_file)

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
           filepath_base %>% file.path(type_output, analysis_type,
                                       subfolder1, subfolder2, subfolder3, subfolder4, subfolder5
           ) %>%
           normalizePath(mustWork = FALSE)

         # Create directory
         dir.create(filepath_dir,
                    showWarnings = FALSE,
                    recursive = TRUE)

         # Update filename
         filename_ = paste0(file_ID, filename, sep = "_")

         # Remove leading and trailing separators
         filename = paste0(filename_ %>% stringr::str_replace("^_", "") %>%
                             stringr::str_replace("_$", "") %>%
                             stringr::str_replace("__", ""), file_ext) # Add file extension

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
#'
#' @examples
format_pars <- function(pars){

  pars_add = with(pars, {

    file_ID = sprintf("nr%d_T%d_ts%.3f_Xsigma%.5f",
                      data_idx, nr_timesteps, timestep, X_sigma)

    # pars_ID

    return(list(analysis_type = analysis_type,
                file_ID = file_ID))
  })

  return(utils::modifyList(pars, pars_add))
}



#' Utility function for foreach loop to save memory
#'
#' @param ...
#'
#' @return NULL
#' @export
#'
#' @examples
cfun <- function(...){NULL}

#' Generate forloop for foreach loop
#'
#' @param ...
#'
#' @return List of forloop parameters
#' @export
#'
#' @examples
get_forloop <- function(...){
  return(tidyr::expand_grid(
    ...
  ) %>% purrr::transpose() %>% unique())
}


# get_filepath_GLV <-
#   function(pars) {
#
#     filepath <-
#       utils::modifyList(pars,
#                         list(
#                           # mainDir, model_name,
#                           type_output,
#                           analysis_type,
#                           analysis_subtype,
#                           file_ID =  sprintf("GLV-full_s%.4f-s%.4f-by%.4f",
#                                              pars$ss_[1], dplyr::last(pars$ss_), pars$s_step),
#                           # sprintf("GLV-null_%s", pars$regime_switch),
#                           pars_ID =  sprintf("_%dtransSteps", pars$nr_trans_s_steps),
#                           filename = ,
#                           file_ext = ".RDS"
#                         )) %>% format_path()
#
#     file_ID = sprintf("nr%d_Xsigma%.5f_%ddays_ts%.3f",
#                              pars$data_idx,
#                              pars$X_sigma,
#                              pars$nr_days_per_s, pars$timestep)
#       pars_ID = sprintf("alphaObs%d_sigmaObs%.2f_%s%.3f_%s-%s-%s",
#                               pars_a$alpha_obs_noise,
#                               pars_a$sigma_obs_noise,
#                               pars$fix_emRad_or_RR,
#                               pars$targetValue,
#                               pars$RQA_type, pars$distNorm, pars$rescaleDist
#       )
#       sprintf("data-%s", pars$daily_or_raw)
#       analysis_subtype = pars$separate_or_continued_ts
#
#
#     return(filepath)
#
#   }
#

