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
# ch = devtools::check()
# ch = devtools::check(vignettes=FALSE)
# devtools::install()
# or:
# devtools::install(build_vignettes = TRUE)

# browseURL("doc/demo.html")

# Add dependencies
# pkgs = c("deSolve", "moments", "utils", "tidyr", "ggplot2", "purrr", "casnet", "rlang", "pracma", "dplyr","rgl", "ggh4x", "stats", "stringr", "zoo", "cowplot", 'magrittr', 'scales', "viridis", "grDevices", "tools", "grid", "invctr", "Matrix", "gsignal", "plyr", "plotly", "tibble")
# for (p in pkgs){
# usethis::use_package(p)
# }
# Allow for use of pipes :
# usethis::use_pipe()

# Add data
# https://grasshoppermouse.github.io/posts/2017-10-18-put-your-data-in-an-r-package/
# # Set up the data-raw directory and data processing script
# # You can use any name you want for your data
# usethis::use_data_raw(name = 'GLV_null')
# # This script in the R directory will contain the documentation.
# # You can use any name you want.
# file.create("R/GLV_null.R")
# # Put your data files into the data-raw folder. Your new package directory should look something like this:
# # Write your data processing code in a data-raw/mydataset.R script. It would look something like this:
#   # data-raw/mydataset.R
#   # Data import and processing pipeline
#   library(readr)
# library(readxl)
# mydataset <- read_csv("data-raw/pendulum data.csv")
# demographics <- read_excel("data-raw/Demographics.xlsx")
# # Data cleaning code here...
# # (Do NOT put data analysis code here!)
# # This should be the last line.
# # Note that names are unquoted.
# # I like using overwrite = T so everytime I run the script the
# # updated objects are saved, but the default is overwrite = F
# usethis::use_data(GLV_trans, overwrite = T)
# sinew::makeOxygen(GLV_null, add_fields = "source")
# sinew::makeOxygen(GLV_null)
# tools::resaveRdaFiles("data/GLV_trans.rda",compress="xz")
# tools::checkRdaFiles("data/")# To check that the compression has been done successfully
# GIven a LazyDataCompression warning: Add
# LazyDataCompression:xz
# to the description file.

#' Reduce size of dataframe by downsampling
#'
#' @param df Dataframe
#' @param X_names Names of columns in df containing observational timeseries to be downsampled
#' @param type Type of downsampling. "average" averages data points in each window; "one_sample" picks one sample in each window.
#' @param win_size Number of data points to be averaged or sampled from in each step
#' @param which_X Time point(s) in window to select
#' @param seed_nr Seed number for random sampling of data point per window
#'
#' @return Downsampled dataframe.
#' @importFrom dplyr .data
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
    slice_func = dplyr::slice_sample

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
    slice_func = dplyr::slice

    if (which_X == "first"){
      n = 1
    } else if (which_X == "middle"){
      n = ceiling(win_size/2)
    } else if (which_X == "last"){
      n = win_size
    } else if (which_X == "random"){
      n = 1
      slice_func = dplyr::slice_sample
    }
  }

  X_day_mu_ = df %>%
    dplyr::mutate(win_nr = rep(1:ceiling(nrow(df)/win_size), each = win_size)[1:nrow(df)]) %>%
    dplyr::group_by(.data$win_nr) %>%
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
#' @param seed_nr Seed number for generating noise
#'
#' @return Noisy dataframe
#' @export
#'
#' @examples
add_obs_noise <- function(df, X_names, noise_mean = 0, noise_sigma = .01, noise_constant = 0, seed_nr = NULL){

  if (!is.null(seed_nr)){
    set.seed(seed_nr)
  }
  # Add white noise to each variable
  df[,X_names] = df[,X_names] + pracma::Reshape(stats::rnorm(n = length(X_names) * nrow(df), mean = noise_mean, sd = noise_sigma), n = nrow(df), m = length(X_names)) + noise_constant

  return(df)
}


#' Style plot
#'
#' @param pl Original plot
#' @param col_pal Named colours
#' @param fs Font sizes
#'
#' @return Styled plot
#' @importFrom ggplot2 theme_bw element_text element_rect unit margin theme
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




#' Utility function to save plots
#'
#' @param pl ggplot object
#' @param filepath_image Filepath image needs to be saved to
#' @param w Width
#' @param h Height
#' @param resolution Resolution
#' @param formats File formats
#'
#' @return Success of saving file
#' @export
#'
#' @examples
save_plot <-
  function(pl,
           filepath_image,
           w = 10,
           h = 10,
           resolution = 200,
           formats = c("tiff", "png", "pdf")[1]) {
    if (rlang::is_empty(filepath_image)) {
      print(sprintf("Filepath to image is NULL!"))
      return()
    }
    grDevices::graphics.off()

    # Make sure there are no leading periods in the file formats
    for (f_idx in 1:length(formats)) {
      f = formats[f_idx]
      f_ = if (stringr::str_sub(f, 1, 1) == ".")
        sprintf("%s", stringr::str_sub(f, 2,-1))
      else
        f # Make sure file extension leads with a period
      formats[f_idx] = f_
    }

    # Save plot in different formats
    if ("tiff" %in% formats) {
      filepath_image %>% tools::file_path_sans_ext() %>% paste0(".tiff") %>%
        grDevices::tiff(
          width = w,
          height = h,
          bg = "white",
          # pointsize = 1 / 300,
          units = 'cm',
          res = resolution
        )
      grid::grid.draw(pl) # Make plot
      grDevices::dev.off()
    }
    if ("png" %in% formats) {
      filepath_image %>% tools::file_path_sans_ext() %>% paste0(".png") %>%
        grDevices::png(
        width = w,
        height = h,
        bg = "transparent",
        units = 'cm',
        res = resolution
      )
      grid::grid.draw(pl) # Make plot
      grDevices::dev.off()
    }

    if ("pdf" %in% formats) {
      filepath_image %>% tools::file_path_sans_ext() %>% paste0(".pdf") %>%
        grDevices::pdf(width = w,
            height = h,
            bg = "white")
      grid::grid.draw(pl) # Make plot
      grDevices::dev.off()
    }
    grDevices::graphics.off()

    return(file.exists(filepath_image))
  }




#' Set up parameters
#'
#' @param model_name Chosen dynamical systems model
#' @param pars_add List of parameters to overwrite or add to default parameters
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
      X0 = c(),
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
      max_iter=100,
      deSolve_method = "lsoda",
      downsample_pars = list(type = c("average", "one_sample")[1],
                             win_size = 50,
                             which_X = c(50, "first", "middle", "last", "random")[1],
                             seed = 123),
      nr_timesteps = 1000, timestep = .01,
      stopifregime = function(out){any(apply(out, 2, is.infinite)) | any(apply(out, 2, is.nan))},
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
         filename_ = stringr::str_c(file_ID, filename, sep = "_")

         # Remove leading and trailing separators
         filename = stringr::str_c(filename_ %>% stringr::str_replace("^_", "") %>%
                             stringr::str_replace("_$", "") %>%
                             stringr::str_replace("__", "_"), file_ext) # Add file extension

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

    if (!is.null(pars[["data_idx"]])){
      file_ID = sprintf("nr%d_T%d_ts%.3f_Xsigma%.5f",
                      data_idx, nr_timesteps, timestep, X_sigma)
    } else {
      file_ID = sprintf("T%d_ts%.3f_Xsigma%.5f",
                        nr_timesteps, timestep, X_sigma)

    }
    # pars_ID

    return(list(file_ID = file_ID))
  })

  return(utils::modifyList(pars, pars_add))
}



#' Utility function for foreach loop to save memory
#'
#' @param ... All outputs from foreach loop
#'
#' @return NULL
#' @export
#'
#' @examples
cfun <- function(...){NULL}

#' Generate forloop for foreach loop
#'
#' @param ... Named forloop parameters
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


#' Add parameters in list as columns in dataframe
#'
#' @param df Dataframe
#' @param for_par List of parameters
#'
#' @return Updated dataframe
#' @export
#'
#' @examples
add_par_as_cols = function(df, for_par){
  df2 = cbind(
    df,
    as.data.frame(t(for_par))[rep(1, nrow(df)), ] %>%
      magrittr::set_rownames(NULL))
  df2 <- df2[, !duplicated(colnames(df2))]
  return(df2)
}
