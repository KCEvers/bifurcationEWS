#' Plot 3D landscape
#'
#' @param df Dataframe
#' @param s_idx Index of bifurcation parameter
#' @param s Value of bifurcation parameter
#' @param plot_idx Index of subplot
#' @param colors Vector of colors
#' @param size_marker Size of marker
#' @param size_line Linewidth
#' @param plot_mode Plotting mode
#'
#' @return Plotly object
#' @export
#'
#' @examples
plot_3D_landscape <- function(df,
                      s_idx = 1, s = .6,
                      plot_idx = 1,
                      colors = c("red"), size_marker = 1,
                      size_line =.1,
                      plot_mode = c("markers", "lines+markers")[2]){
    name1 = sprintf("idx: %d, s = %.4f", s_idx, s)
    color = colors[plot_idx]
    scene_name = sprintf('scene%d', plot_idx)
    fig1 <- plotly::plot_ly(df, scene = scene_name, name = name1,
                    type="scatter3d", mode = plot_mode,
                    x = ~X1, y = ~X2, z = ~X3,
                    opacity = 1,
                    marker = list(color = color, size = size_marker),
                    line = list(color = color, width = size_line)
    )

    return(plotly::plotly_build(fig1))
  }

#' Plot transition in 3D
#'
#' @inheritParams plot_3D_landscape
#' @param bifpar_idxs Bifurcation parameter indices to plot
#' @param bifpar_list List of bifurcation parameters to index
#' @param nrows Number of rows
#' @param ncols Number of columns
#' @param x1_range Range of x1-axis
#' @param x2_range Range of x2-axis
#' @param x3_range Range of x3-axis
#'
#' @return Plotly object
#' @importFrom dplyr filter mutate group_map group_by .data
#' @export
#'
#' @examples
plot_3D_landscape_trans <- function(df,
                                   bifpar_idxs,
                                   bifpar_list,
                                   nrows = NULL,
                                   ncols = NULL,
                                   x1_range = c(0,1), x2_range = c(0,1), x3_range = c(0,1),
                                   size_marker = 1, size_line= .1,
                                   plot_mode = c("markers", "lines+markers")[2]){

    # Check if dimensions of plot work
    if (is.null(nrows) & !is.null(ncols)){
      nrows = ceiling(length(bifpar_idxs) / ncols)
    } else if (!is.null(nrows) & is.null(ncols)){
      ncols = ceiling(length(bifpar_idxs) / nrows)
    } else if (is.null(nrows) & is.null(ncols)){
      nrows = 2
      ncols = ceiling(length(bifpar_idxs) / nrows)
    }

    if (length(bifpar_idxs) > (nrows*ncols)){
      message("Too many plots for this number of rows and columns!")
      return()
    }

    scenes = list(xaxis = list(title = 'X1',
                               range = x1_range,
                               aspectmode='cube'),
                  yaxis = list(title = 'X2',
                               range = x2_range,
                               aspectmode='cube'),
                  zaxis = list(title = 'X3',
                               range = x3_range,
                               aspectmode='cube'),
                  camera = list(eye = list(x=2, y=2, z = 2) ))

    colors = viridis::viridis(n = length(bifpar_idxs), end = .85)
    plot_list = df %>%
      filter(bifpar_idx %in% bifpar_idxs) %>%
      group_by(bifpar_idx) %>%
      mutate(plot_idx = which(bifpar_idxs == unique(bifpar_idx))) %>%
      group_map(~ plot_3D_landscape(.x,
                                           s_idx = unique(.y$bifpar_idx),
                                           s = unlist(purrr::map(bifpar_list[unique(.y$bifpar_idx)], "s")),
                                           plot_idx = unique(.x$plot_idx),
                                           colors=colors,
                                 size_marker = size_marker,
                                 size_line =size_line,
                                 plot_mode=plot_mode
      ))

    # Placing plots on a grid
    plot_grid = expand.grid(row_nr = 1:nrows, col_nr = 1:ncols) %>% arrange(.data$row_nr, .data$col_nr) %>% slice(1:length(bifpar_idxs))
    domains_cols_ = seq(0,1,length.out = ncols + 1)
    domains_rows_ = seq(0,1,length.out = nrows + 1)
    domains_cols = plyr::llply(1:ncols, function(i){c(domains_cols_[i], domains_cols_[i+1])})
    domains_rows = plyr::llply(1:nrows, function(i){c(domains_rows_[i], domains_rows_[i+1])}) %>% rev() # reverse to set bottom to top -> top to bottom

    # Merge plots
    fig <- plotly::subplot(plot_list,
                           nrows=nrows,
                           margin = .1, shareX=TRUE, shareY=TRUE) %>%
      plotly::layout(
        legend= list(itemsizing='constant')
      )

    # Update position of subplots
    for (plot_idx in 1:length(bifpar_idxs)){
      xy_domain = list(
        x = domains_cols[[plot_grid[plot_idx, "col_nr"]]], # x = horizontal domain
        y = domains_rows[[plot_grid[plot_idx, "row_nr"]]] # y = vertical domain
      )
      fig$x$layout[[sprintf("scene%s", ifelse(plot_idx != 1, toString(plot_idx), ""))]] = utils::modifyList(scenes, list(domain=xy_domain))
    }

    # Save plotly image: follow instructions on save_image {plotly}
    # htmlwidgets::saveWidget(fig, "p1.html", selfcontained = T, libdir = "lib") #https://plotly-r.com/saving
    # In this case, if you wanted to share "p1.html" and/or "p2.html" with someone else, make sure to include the libdir folder, perhaps via a zip file:
    # zip("p1.zip", c("p1.html", "lib"))


    return(fig)
  }
