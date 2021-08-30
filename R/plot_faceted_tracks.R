# plotting tracks ####
utils::globalVariables(c(":=", "b_rank", "discard_top", "gfp_nb", "time_sec", "id", "parent_id"))

plot_faceted_var_tracks <- function(.df, .var_col=gfp_nb, .time_col=time_sec, 
                                    .cell_col=rlang::expr(id), .parent_col=parent_id, .facet_col=b_rank, 
                                    .log=FALSE, .col=NULL, .alpha_col=discard_top, .show_cellid=FALSE, .show_all=FALSE, 
                                    .facet_labeller=NULL, .gg_theme=ggplot2::theme_get()) {
# programmed following https://dplyr.tidyverse.org/dev/articles/programming.html
# currently with a dirty (?) hack to pass the `id` expression as argument while this is already the name 
# of a dplyr function in the global environment (pass as expr and then "force" with !!)
  # browser()
  if (nrow(.df) == 0) return(ggplot2::ggplot())
  .df_div <- .df %>% dplyr::group_by(!! {{.cell_col}}) %>% # !! forces early evaluation (otherwise `id` is a function name in dplyr)
    dplyr::do((function(.dfgl, .dfc) { # compute_div_between_facets
      # .dfgl: dataframe of the entire GL
      # .dfc: dataframe of the current cell
      
      # browser()
      # case: no parent
      if (unique(dplyr::pull(.dfc, {{.parent_col}})) < 0) return(data.frame())
      .dfp <- dplyr::filter(.dfgl, {{.cell_col}} == unique(dplyr::pull(.dfc, {{.parent_col}})) )
      if (nrow(.dfp) == 0) return(data.frame())
      
      # otherwise a parent exists
      .facet_c <- unique(dplyr::pull(.dfc, {{.facet_col}}))
      .facet_p <- unique(dplyr::pull(.dfp, {{.facet_col}}))
      .var_p_last <- dplyr::filter(.dfp, {{.time_col}} == max({{.time_col}}) ) %>% dplyr::pull({{.var_col}})
      .time_p_last <- dplyr::filter(.dfp, {{.time_col}} == max({{.time_col}}) ) %>% dplyr::pull({{.time_col}})
      .time_c_first <- dplyr::filter(.dfc, {{.time_col}} == min({{.time_col}}) ) %>% dplyr::pull({{.time_col}})
      .dt <- .time_c_first - .time_p_last
      # case: parent cell in the same panel
      if ( .facet_c == .facet_p ) {
        .out <- dplyr::filter(.dfc, {{.time_col}} == min({{.time_col}})) %>% 
          dplyr::select({{.cell_col}}, {{.facet_col}}, {{.time_col}}, {{.var_col}}) %>%
        dplyr::mutate("{{.time_col}}_p" := .time_p_last, "{{.var_col}}_p" := .var_p_last)
      }

      # case: parent cell in a panel under the cell of interest
      if ( .facet_c > .facet_p ) {
        .out <- dplyr::bind_rows(
          dplyr::filter(dplyr::ungroup(.dfc), {{.time_col}} == min({{.time_col}}) ) %>%
            dplyr::select({{.cell_col}}, {{.facet_col}}, {{.time_col}}, {{.var_col}}) %>%
            dplyr::mutate("{{.time_col}}_p" := {{.time_col}} - .dt/2, "{{.var_col}}_p" := ifelse(.log, 0, -Inf) ),
          dplyr::filter(dplyr::ungroup(.dfp), {{.time_col}} == max({{.time_col}}) ) %>%
            dplyr::select({{.cell_col}}, {{.facet_col}}, {{.time_col}}, {{.var_col}}) %>%
            dplyr::mutate("{{.time_col}}_p" := {{.time_col}} + .dt/2, "{{.var_col}}_p" := Inf )
        )
        # add vertical line in intermediate facets if required
        if(diff(dplyr::pull(.out, {{.facet_col}})) < -1) {
          .out <- dplyr::bind_rows(.out, dplyr::tibble(
            "{{.cell_col}}" := unique(dplyr::pull(.dfc, {{.cell_col}})), 
            "{{.facet_col}}" := seq(min(dplyr::pull(.out, {{.facet_col}}))+1, max(dplyr::pull(.out, {{.facet_col}}))-1),
            "{{.time_col}}" := .time_c_first-.dt/2, 
            "{{.time_col}}_p" := .time_c_first-.dt/2,
            "{{.var_col}}" := ifelse(.log, 0, -Inf), 
            "{{.var_col}}_p" := Inf) )
        }
      }

      # case: parent cell in a panel above the cell of interest
      if ( .facet_c < .facet_p ) {
        .out <- dplyr::bind_rows(
          dplyr::filter(dplyr::ungroup(.dfc), {{.time_col}} == min({{.time_col}}) ) %>%
            dplyr::select({{.cell_col}}, {{.facet_col}}, {{.time_col}}, {{.var_col}}) %>%
            dplyr::mutate("{{.time_col}}_p" := {{.time_col}} - .dt/2, "{{.var_col}}_p" := Inf ),

          dplyr::filter(dplyr::ungroup(.dfp), {{.time_col}} == max({{.time_col}}) ) %>%
            dplyr::select({{.cell_col}}, {{.facet_col}}, {{.time_col}}, {{.var_col}}) %>%
            dplyr::mutate("{{.time_col}}_p" := {{.time_col}} + .dt/2, "{{.var_col}}_p" := ifelse(.log, 0, -Inf) )
        )
        # add vertical line in intermediate facets if required
        if(diff(dplyr::pull(.out, {{.facet_col}})) > 1) {
          .out <- dplyr::bind_rows(.out, dplyr::tibble(
            "{{.cell_col}}" := unique(dplyr::pull(.dfc, {{.cell_col}})), 
            "{{.facet_col}}" := seq(min(dplyr::pull(.out, {{.facet_col}}))+1, max(dplyr::pull(.out, {{.facet_col}}))-1),
            "{{.time_col}}" := .time_c_first-.dt/2, 
            "{{.time_col}}_p" := .time_c_first-.dt/2,
            "{{.var_col}}" := ifelse(.log, 0, -Inf), 
            "{{.var_col}}_p" := Inf) )
        }
      }
      return(.out)
    })(.df, .))

  if (nrow(.df_div) == 0)
    warning('unable to compute divisions for this dataset (hint: check that `discard_top` is not pruning the entire tree).')

  # handle plotting arguments after forcing .cell_col
  if (is.null(rlang::enexpr(.col))) { .col_col <- rlang::enexpr(.cell_col) } else { .col_col <- rlang::quo(TRUE) }
  if (is.null(rlang::enexpr(.alpha_col))) .alpha_col <- rlang::quo(FALSE)
  
  # browser()
  .facet_min <- min(dplyr::pull(.df, {{.facet_col}}))
  .facet_max <- max(dplyr::pull(.df, {{.facet_col}}))
  
  .pl <- ggplot2::ggplot() +
    # facets alternated background
    ggplot2::geom_rect(ggplot2::aes(xmin=-Inf, xmax=Inf, ymin=ifelse(.log, 0, -Inf), ymax=Inf), alpha=.05, 
              data=dplyr::tibble("{{.facet_col}}" := seq(.facet_min, .facet_max, 2)) ) +
    # theme(.gg_theme, complete=TRUE) +
    ggplot2::theme(panel.spacing = grid::unit(0, "lines"), panel.border=ggplot2::element_blank()) # ,strip.background = element_blank(), strip.text = element_blank()  
  if (is.null(.facet_labeller)) {
    .pl <- .pl + ggplot2::facet_grid(rows=ggplot2::vars({{.facet_col}}), as.table=FALSE)
  } else {
    .pl <- .pl + ggplot2::facet_grid(rows=ggplot2::vars({{.facet_col}}), as.table=FALSE, labeller=ggplot2::as_labeller(.facet_labeller))
  }
  
  # create the column names for parent's variables
  .time_col_p <- rlang::ensym(.time_col) %>% rlang::as_string() %>% paste0("_p") %>% rlang::parse_expr()
  .var_col_p <- rlang::ensym(.var_col) %>% rlang::as_string() %>% paste0("_p") %>% rlang::parse_expr()
  # show divisions (lty='11' means densely dotted line using hex notation)
  if (nrow(.df_div) > 0)
    .pl <- .pl + 
    ggplot2::geom_segment(ggplot2::aes(
      x={{.time_col}}, xend={{.time_col_p}}, y={{.var_col}}, yend={{.var_col_p}}, col=factor({{.col_col}}), alpha=!{{.alpha_col}} 
    ), alpha=.3, lty='11', data=.df_div)
  
  # show cell traces
  .pl <- .pl + 
    ggplot2::geom_path(ggplot2::aes({{.time_col}}, {{.var_col}}, group=factor({{.cell_col}}), col=factor({{.col_col}}), alpha=!{{.alpha_col}} ), data=.df)
  
  # show cell numbers
  if (.show_cellid)
    .pl <- .pl + ggplot2::geom_text(ggplot2::aes(
      {{.time_col}}, {{.var_col}}, label={{.cell_col}}, col=factor({{.col_col}})
    ), size=2, hjust=0, vjust=1, 
    # data=.df %>% dplyr::group_by({{.cell_col}}) %>% dplyr::filter(row_number()==1) ) #dplyr::row_number()
    data=.df %>% dplyr::group_by({{.cell_col}}) %>% dplyr::slice(1L) )
  
  # show all traces in one panel
  if (.show_all)
    .pl <- .pl + ggplot2::geom_path(ggplot2::aes(
      {{.time_col}}, {{.var_col}}, col=factor({{.col_col}}), alpha=!{{.alpha_col}}
    ), data=.df %>% dplyr::mutate("{{.facet_col}}" := -1) )
  
  # adjust scales
  if (is.null(rlang::enexpr(.col))) {
    .pl <- .pl + ggCustomTJ::scale_colour_periodic_brewer(guide='none')
  } else {
    .pl <- .pl + ggplot2::scale_colour_manual(values=c('TRUE'=.col))
  }
  .pl <- .pl + 
    ggplot2::scale_alpha_manual(values=c('TRUE'=1, 'FALSE'=0.3)) +
    ggplot2::guides(color='none', alpha='none')
  if(.log) .pl <- .pl + ggplot2::scale_y_continuous(trans='log2')
  
  return(.pl)
}

