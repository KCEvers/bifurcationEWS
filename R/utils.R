# Minimal outline building R package
#
# Create .R script
# usethis::use_r("example")
#
# Clear workspace
# devtools::load_all()
#
# devtools::build_vignettes()
# devtools::check()
# devtools::document()
# devtools::install()
# or:
# devtools::install(build_vignettes = TRUE)

# browseURL("doc/demo.html")

# Add dependencies
# pkgs = c("deSolve", "moments", "utils", "tidyr", "ggplot2", "purrr", "casnet", "rlang", "pracma", "dplyr","rgl")
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
