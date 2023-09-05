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
# pkgs = c("deSolve", "moments", "utils", "tidyr", "ggplot2", "purrr", "casnet", "rlang", "pracma", "dplyr","rgl", "ggh4x", "stats", "stringr", "zoo", "cowplot", 'magrittr', 'scales', "viridis", "grDevices", "tools", "grid")
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
#' @param seed_nr Seed number for random sampling of data point per window
#'
#' @return Downsampled dataframe.
#' @importFrom magrittr `%>%`
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
    dplyr::group_by(win_nr) %>%
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
  pl <- pl + ggplot2::theme_bw() +
    ggplot2::theme(
      text = ggplot2::element_text(family = fs["family"]),
      # Change font
      plot.title = ggplot2::element_text(size = fs["plot.title"]),
      plot.subtitle = ggplot2::element_text(size = fs["plot.subtitle"]),
      axis.text = ggplot2::element_text(size = fs["axis.text"]),
      axis.title = ggplot2::element_text(size = fs["axis.title"]),
      legend.text = ggplot2::element_text(size = fs["legend.text"]),
      legend.title = ggplot2::element_text(size = fs["legend.title"]),
      # # Increase font size facet labels
      strip.text.y = ggplot2::element_text(
        size = as.numeric(fs["strip.text.y"]) + 2,
        margin = ggplot2::margin(0.1, 0.1, 0.1, 0.1, "cm")
      ),
      strip.text.x = ggplot2::element_text(
        size = as.numeric(fs["strip.text.x"]) + 2,
        margin = ggplot2::margin(0.1, 0.1, 0.1, 0.1, "cm")
      )
      # panel.spacing.x = unit(0.2, "cm"),
      # panel.spacing.y = unit(0.2, "cm") # Distance between facets
    ) +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill = col_pal["col_facet_labels"], color = col_pal["col_facet_labels"])) + # Change facet rectangle colour
    ggplot2::theme(strip.text = ggplot2::element_text(colour = 'white')) +
    # theme(legend.position = "none") + # Remove group legend
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
          plot.subtitle = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::theme(plot.margin = ggplot2::unit(c(0, 0, 0, 0), "pt")) +
    ggplot2::theme(legend.title.align = 0.5,
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
#' @importFrom magrittr `%>%`
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
#' @importFrom magrittr `%>%`
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
#' @importFrom magrittr `%>%`
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