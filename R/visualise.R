#' Plot 3D landscape
#'
#' @param df Dataframe
#' @param xvar Name of variable to plot on X-axis
#' @param yvar Name of variable to plot on Y-axis
#' @param zvar Name of variable to plot on Z-axis
#' @param bifpar_idx Index of bifurcation parameter
#' @param bifpar_value Value of bifurcation parameter
#' @param bifpar_name Name of bifurcation parameter
#' @param plot_idx Index of subplot
#' @param color Plotting colour
#' @param size_marker Size of marker
#' @param size_line Linewidth
#' @param plot_mode Plotting mode
#'
#' @return Plotly object
#' @export
#'
#' @examples
plot_3D_landscape <- function(df,
                              xvar = "X1", yvar = "X2", zvar = "X3",
                              bifpar_name = "s",
                      bifpar_idx = 1, bifpar_value = .6,
                      plot_idx = 1,
                      color = c("red"),
                      size_marker = 1,
                      size_line =.1,
                      plot_mode = c("markers", "lines+markers")[2]){

    name1 = sprintf("idx: %d, %s = %.4f", bifpar_idx, bifpar_name, bifpar_value)
    scene_name = sprintf('scene%d', plot_idx)
    fig1 <- plotly::plot_ly(df, scene = scene_name, name = name1,
                    type="scatter3d", mode = plot_mode,
                    x = ~.data[[xvar]], y = ~.data[[yvar]], z = ~.data[[zvar]],
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
#' @param bifpar_values Values of bifurcation parameter
#' @param nrows Number of rows
#' @param ncols Number of columns
#' @param x_range Range of X-axis
#' @param y_range Range of Y-axis
#' @param z_range Range of Z-axis
#' @param colors List or single string of plotting colors
#' @param camera List of camera settings
#' @param alternating_height Logical indicating whether the subplots should be plotted on a grid with alternating height; only works for an uneven number of subplots and nrows = 1
#' @param xaxis_title Title of X-axis
#' @param yaxis_title Title of Y-axis
#' @param zaxis_title Title of Z-axis
#' @param showticklabels Logical indicating whether to show axis tick labels
#' @param showlegend Logical indicating whether to show the legend
#'
#' @return Plotly object
#' @importFrom dplyr filter mutate group_map group_by .data
#' @export
#'
#' @examples
plot_3D_landscape_trans <- function(df,
                                   bifpar_idxs,
                                   bifpar_values,
                                   bifpar_name = "s",
                                   nrows = NULL,
                                   ncols = NULL,
                                   x_range = c(0,1), y_range = c(0,1), z_range = c(0,1),
                                   xaxis_title = NULL,
                                   yaxis_title = NULL,
                                   zaxis_title = NULL,
                                   xvar = "X1", yvar = "X2", zvar = "X3",
                                   showticklabels = TRUE,
                                   showlegend = TRUE,
                                   size_marker = 1, size_line= .1,
                                   plot_mode = c("markers", "lines+markers")[2],
                                   colors = NULL,
                                   camera = list(eye = list(x=2, y=2, z = 2)),
                                   alternating_height = FALSE){

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


  if (is.null(xaxis_title)){
    xaxis_title = "X1"
  }
  if (is.null(yaxis_title)){
    yaxis_title = "X2"
  }
  if (is.null(zaxis_title)){
    zaxis_title = "X3"
  }

  if (is.null(colors)){
    colors = viridis::viridis(n = length(bifpar_idxs), end = .85)
  } else if (length(colors) != length(bifpar_idxs)){
    colors = rep(colors[1], length(bifpar_idxs))
  }

  if (length(size_marker) < length(bifpar_idxs)){
    size_marker = rep(size_marker[1], length(bifpar_idxs))
  }
  if (length(size_line) < length(bifpar_idxs)){
    size_line = rep(size_line[1], length(bifpar_idxs))
  }
   grid_color = "grey100"
    scenes = list(xaxis = list(title=list(text=xaxis_title, standoff = 50,
                                          font=list(size=5, family='Courier')),
                               nticks = 2,
                               zerolinecolor = grid_color,
                               grid_width = .01,
                               zerolinewidth = .75,
                               gridcolor = grid_color,
                               showgrid = F,zeroline = T,
                               showticklabels=showticklabels,
                               aspectmode='cube',
                               range = x_range
                              ),
                  yaxis = list(title=list(text=yaxis_title, standoff = 50,
                                          font=list(size=5, family='Courier')),
                               nticks = 2,
                               zerolinecolor = grid_color,
                               gridcolor = grid_color,
                                grid_width = .01,
                               showgrid = F,zeroline = T,
                               zerolinewidth = .75, showticklabels=showticklabels,
                               aspectmode='cube',
                               range = y_range
                               ),
                  zaxis = list(title=list(text=zaxis_title, standoff = 50,
                                          font=list(size=5, family='Courier')),
                               nticks = 2,
                               zerolinecolor = grid_color,
                               gridcolor =grid_color,
                               grid_width = .01,
                               showgrid = F, zeroline = T,
                               zerolinewidth = .75, showticklabels=showticklabels,
                               aspectmode='cube',
                               range = z_range
                              ),
                  camera = camera)

    plot_list = df %>%
      filter(.data$bifpar_idx %in% bifpar_idxs) %>%
      group_by(.data$bifpar_idx) %>%
      mutate(plot_idx = which(bifpar_idxs == unique(.data$bifpar_idx))) %>%
      group_map(~ plot_3D_landscape(.x,
                                           bifpar_idx = unique(.y$bifpar_idx),
                                    bifpar_value = bifpar_values[unique(.x$plot_idx)],
                                    plot_idx = unique(.x$plot_idx),
                                           color=colors[unique(.x$plot_idx)],
                                 size_marker = size_marker[unique(.x$plot_idx)],
                                 size_line =size_line[unique(.x$plot_idx)],
                                 xvar = xvar, yvar = yvar, zvar = zvar,
                                 plot_mode=plot_mode
      ))

    # Placing plots on a grid
     if (alternating_height & nrows == 1 & length(bifpar_idxs) %% 2 != 0){ # Uneven number
      plot_grid = data.frame(col_nr = 1:ncols) %>% dplyr::mutate(row_nr = rep(1:2, length.out = ncols))

      w_overlap = 1 / (ncols+6) # In order to scale overlap, could be changed
      w = 1 / ncols + (w_overlap*((ncols-1)/ncols))
      domains_cols = plyr::llply(seq.int(ncols), function(i){c(w*(i-1)-(w_overlap*(i-1)), w*i-(w_overlap*(i-1)))})

      h_overlap = .2
      implied_nrows = 2
      h = 1 / implied_nrows + (h_overlap*((implied_nrows-1)/implied_nrows))
      domains_rows = plyr::llply(seq.int(implied_nrows), function(i){c(h*(i-1)-(h_overlap*(i-1)), h*i-(h_overlap*(i-1)))})%>% rev()

     } else {
       plot_grid = expand.grid(row_nr = 1:nrows, col_nr = 1:ncols) %>% arrange(.data$row_nr, .data$col_nr) %>% slice(1:length(bifpar_idxs))
       domains_cols_ = seq(0,1,length.out = ncols + 1)
       domains_rows_ = seq(0,1,length.out = nrows + 1)
       domains_cols = plyr::llply(1:ncols, function(i){c(domains_cols_[i], domains_cols_[i+1])})
       domains_rows = plyr::llply(1:nrows, function(i){c(domains_rows_[i], domains_rows_[i+1])}) %>% rev() # reverse to set bottom to top -> top to bottom

    }


    # Merge plots
    fig <- plotly::subplot(plot_list,
                           nrows=nrows,
                           margin = 0,
                           shareX=TRUE, shareY=TRUE) %>%
      plotly::layout(
     automargin=F,
     # autosize = F,
        margin = list(pad = 5, t = 5, b = 5, l = 5, r = 5),
        legend= list(itemsizing='constant'),
        showlegend = showlegend
      )

    # Update position of subplots
    for (plot_idx in 1:length(bifpar_idxs)){
      xy_domain = list(
        x = domains_cols[[plot_grid[plot_idx, "col_nr"]]], # x = horizontal domain
        y = domains_rows[[plot_grid[plot_idx, "row_nr"]]] # y = vertical domain
      )
      scene_name = sprintf("scene%s", ifelse(plot_idx != 1, toString(plot_idx), ""))
      fig$x$layout[[scene_name]] = utils::modifyList(scenes, list(domain=xy_domain
                                                                  # xaxis = list(showticklabels=showticklabels),
                                                                  # yaxis = list(showticklabels=showticklabels),
                                                                  # zaxis = list(showticklabels=showticklabels)
                                                                  )
                                                     )
    }

    # fig$x$layout$margin = c(t = 5, b = 5, l = 5, r = 5)

    # Save plotly image: follow instructions on save_image {plotly}
    # htmlwidgets::saveWidget(fig, "p1.html", selfcontained = T, libdir = "lib") #https://plotly-r.com/saving
    # In this case, if you wanted to share "p1.html" and/or "p2.html" with someone else, make sure to include the libdir folder, perhaps via a zip file:
    # zip("p1.zip", c("p1.html", "lib"))


    return(fig)
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
                         col_facet_labels = "#3C1A42FF"
                         #scales::viridis_pal(option = "rocket", direction = -1)(20)[17]
                                                                  ),
                       fs = c()) {
  # Font parameters
  fs_default = c(
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
  )
  fs = c(fs_default[setdiff(names(fs_default), names(fs))], fs)

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
    theme(plot.margin = unit(c(10,10,10,10), "pt")) +
    theme(legend.title.align = 0.5,
          # Text label alignment. Number from 0 (left) to 1 (right)
          legend.text.align = 0.5)

  return(pl)
}



#' Plot bifurcation diagram
#'
#' @param regime_list List with detected regimes, output of find_regimes()
#' @param sel_variables Selected variables to plot
#' @param filepath_image Filepath to save image to; if not specified, no plot is saved.
#' @param col_palette_regimes Name of colour palette to colour regimes. Must be one of scico's palettes, see all options with scico::scico_palette_names()
#' @param alpha_regime Transparancy of regime colours
#' @param sz_point_peaks Point size of peaks
#' @param resolution Image resolution
#' @param factor_resolution Factor to multiply image size by to increase resolution
#' @param reverse_x_axis Logical indicating whether the x-axis should be reversed
#' @param nbreaks_x Number of ticks on the x-axis
#'
#' @return Styled bifurcation digram
#' @importFrom ggplot2 coord_cartesian element_blank element_line geom_point geom_rect ggplot scale_x_continuous scale_x_reverse scale_y_continuous label_parsed guide_colorbar labs aes
#' @export
#'
#' @examples
plot_bifdiag = function(regime_list,
                        sel_variables = NULL,
                        filepath_image = NULL,
                        col_palette_regimes = c("lapaz","batlowW")[2],
                        alpha_regime = .99,
                        sz_point_peaks = .0001,
                        resolution = 200,
                        factor_resolution = 3,
                        reverse_x_axis = F, nbreaks_x = 8
){

  if (is.null(sel_variables)){
    sel_variables = unique(regime_list$peaks_df$variable)
  }


  # Relabel variables
  variable_labs <- list("X1" = latex2exp::TeX("$X_1$", output = 'character'),
                        "X2" = latex2exp::TeX("$X_2$", output = 'character'),
                        "X3" = latex2exp::TeX("$X_3$", output = 'character'),
                        "X4" = latex2exp::TeX("$X_4$", output = 'character')
  )

  # Relabel regimes according to the order of periodicity
  regimes_period_nr = regime_list$regimes %>%
    rowwise() %>%
    dplyr::mutate_at(regime_list$X_names, ~ifelse(grepl("Chaotic", .), NA,
                                                   stringr::str_replace_all(., "Period-", "") %>%
                                                     as.character() %>% as.numeric())
    ) %>% ungroup() %>% tidyr::gather(.data$variable, .data$period, -setdiff(colnames(.), regime_list$X_names)) %>%
    dplyr::filter(.data$variable %in% sel_variables) %>%
    mutate(variable_name =
             dplyr::recode(.data$variable, !!!variable_labs))

  # Set up colour bar
  na_colour = dplyr::last(scico::scico(100, palette = col_palette_regimes))
  max_period_nr = max(regimes_period_nr %>% pull(.data$period), na.rm = TRUE)
  if (max_period_nr == 1){
    breaks_colourbar = c(1, 100)
    labels_colourbar = c("Node", "Chaotic")
  } else {
    breaks_colourbar = c(1, 2**seq(1, max(c(ceiling(log2(round(max_period_nr, 0))), 1)), by = 1)) # Logarithmic sequence
    labels_colourbar = c("Node", paste0("Period ", breaks_colourbar[-c(1,length(breaks_colourbar))]), "Chaotic")
  }

  # Set up plot
  pl = regime_list$peaks_df %>%
    dplyr::filter(.data$variable %in% sel_variables) %>%
    mutate(variable_name =
             dplyr::recode(.data$variable, !!!variable_labs)) %>%
    ggplot()

  # Plot regimes
  pl <- pl +
    geom_rect(
      data = regimes_period_nr,
      aes(xmin = .data$start_bifpar_idx, xmax=.data$end_bifpar_idx, ymin = -Inf,ymax=Inf,
          fill = .data$period),
      linewidth = 0 # Remove color
    ) +
    scico::scale_fill_scico(
      palette = col_palette_regimes,
      alpha = alpha_regime,
      begin = .2, end = 1,
      trans = "log2",
      limits = range(breaks_colourbar),
      breaks = breaks_colourbar,
      labels = labels_colourbar,
      na.value = na_colour,
      direction = 1,
      name = "",
      guide = guide_colorbar(
        barheight = 30, #14
        barwidth = 3, #14
        reverse = TRUE,
        title.position = "top",
        title.hjust = 0.5,
        nbin = 20
      ))


  pl <- pl +
    geom_point(aes(x = .data$bifpar_idx, y = .data$X),
               col = "gray10",
               size = sz_point_peaks) +
    labs(y = "",
         x = latex2exp::TeX("$s"),
         title = "") +
    ggh4x::facet_grid2(variable_name ~ .,
                       labeller = label_parsed
    )

  # Set proper axis ranges
  min_x = min(regime_list$peaks_df$bifpar_idx)
  max_x = max(regime_list$peaks_df$bifpar_idx)
  breaks_x = seq(1, length(regime_list$bifpar_list), length.out = nbreaks_x)
  labels_x = unname(unlist(regime_list$bifpar_list[breaks_x]))
  if (reverse_x_axis){
    pl <- pl +
      scale_x_reverse(breaks=rev(breaks_x), labels = labels_x) +
      coord_cartesian(
                      ylim = c(-0.035, 1.035), expand = FALSE, clip = 'on')

  } else {
    pl <- pl +
      scale_x_continuous(breaks=breaks_x, labels = labels_x) +
      coord_cartesian(
                      ylim = c(-0.035, 1.035), expand = FALSE, clip = 'on')
  }

  pl <- pl +
    scale_y_continuous(breaks = c(0, .25, .5, .75, 1))

  factor_fs = 1 #ifelse(length(sel_variables) == 1, 2, 4)
  fs = c(
    "strip.text.x" = 12*factor_fs,
    "strip.text.y" = 12*factor_fs,
    "axis.text" = 8*factor_fs,
    "axis.title" = 16*factor_fs,
    "legend.text" = 12*factor_fs,
    "legend.title" = 14*factor_fs
  )

  # Style plot
  pl = style_plot(pl, fs = fs) +
    # Remove background grid
    theme(axis.line = element_line(color='black'),
          plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank())

  # Save plot
  if (!is.null(filepath_image)){

    save_plot(
      pl,
      filepath_image,
      w = 10*factor_resolution,
      h = 10+2.5*length(sel_variables)*factor_resolution,
      resolution = resolution,
      formats = ".eps")
    print(filepath_image)
  }
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
           formats = c("tiff", "png", "pdf", "eps")[1]) {
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
    if ("svg" %in% formats) {
      filepath_image %>% tools::file_path_sans_ext() %>% paste0(".svg") %>%
        grDevices::svg(width = w,
                       height = h,
                       bg = "white")
      grid::grid.draw(pl) # Make plot
      grDevices::dev.off()
    }
    if ("eps" %in% formats){
      # grDevices::setEPS()
      # grDevices::postscript(filepath_image %>% tools::file_path_sans_ext() %>% paste0(".eps"), fonts = "serif")
      # grid::grid.draw(pl) # Make plot
      # grDevices::dev.off()
      grDevices::cairo_ps(filepath_image %>% tools::file_path_sans_ext() %>% paste0(".eps"),
                          width = w,
                          height = h,
                          bg = "white",pointsize = 2)
      print(pl)
      grDevices::dev.off()
    }
    grDevices::graphics.off()

    return(file.exists(filepath_image))
  }

