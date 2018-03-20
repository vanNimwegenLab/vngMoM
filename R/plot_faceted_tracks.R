# plotting tracks ####
utils::globalVariables(c(".", "dt"))

plot_faceted_var_tracks <- function(.df, .var_col='gfp_nb', .time_col='time_sec', 
                                    .cell_col='id', .parent_col='parent_id', .facet_col='b_rank', 
                                    .log=FALSE, .col=NA, .show_cellid=FALSE, .show_all=FALSE, 
                                    .facet_labeller=NULL, .gg_theme=ggplot2::theme_get(), .dt=dt) {
  if (nrow(.df) == 0) return(ggplot2::ggplot())
  
  .df_div <- .df %>% dplyr::group_by_(.cell_col) %>% 
    dplyr::do((function(.dfgl, .dfc) { # compute_div_between_facets
      # .dfgl: dataframe of the entire GL
      # .dfc: dataframe of the current cell
      if (unique(.dfc[[.parent_col]]) < 0) return(data.frame())
      .dfp <- dplyr::filter_(.dfgl, lazyeval::interp(~ cell_id == unique(parent_id), 
                                              cell_id=as.name(.cell_col), parent_id=.dfc[[.parent_col]]))
      if (nrow(.dfp) == 0) return(data.frame())
      
      if (unique(.dfc[[.facet_col]]) == unique(.dfp[[.facet_col]])) {
        .varp <- dplyr::filter_(.dfp, sprintf('%s==max(%s)', .time_col, .time_col))[[.var_col]]
        .out <- dplyr::filter_(.dfc, sprintf('%s==min(%s)', .time_col, .time_col)) %>% 
          dplyr::select_(.cell_col, .facet_col, .time_col, .var_col) %>%
          dplyr::mutate_(.dots=list(lazyeval::interp(~tvar-.dt*60, tvar=as.name(.time_col), dt=.dt), 
                             lazyeval::interp(~.varp, .varp=.varp)) %>% 
                    stats::setNames(paste0(c(.time_col, .var_col), 'p')))
      }
      if (unique(.dfc[[.facet_col]]) > unique(.dfp[[.facet_col]])) {
        .out <- dplyr::bind_rows(
          dplyr::filter_(dplyr::ungroup(.dfc), sprintf('%s==min(%s)', .time_col, .time_col)) %>%
            dplyr::select_(.cell_col, .facet_col, .time_col, .var_col) %>%
            dplyr::mutate_(.dots=list(lazyeval::interp(~tvar-.dt*60/2, tvar=as.name(.time_col), dt=.dt), as.character(ifelse(.log, 0, -Inf))) %>% 
                      stats::setNames(paste0(c(.time_col, .var_col), 'p'))) ,
          dplyr::filter_(dplyr::ungroup(.dfp), sprintf('%s==max(%s)', .time_col, .time_col)) %>%
            dplyr::select_(.cell_col, .facet_col, .time_col, .var_col) %>%
            dplyr::mutate_(.dots=list(lazyeval::interp(~tvar+.dt*60/2, tvar=as.name(.time_col), dt=.dt), "Inf") %>% 
                      stats::setNames(paste0(c(.time_col, .var_col), 'p')))
        )
        if(diff(.out[[.facet_col]]) < -1) {
          .out <- rbind(.out,
                        data.frame(id=unique(.dfc$id), seq(min(.out[[.facet_col]])+1, max(.out[[.facet_col]])-1),
                                   time=min(.dfc[[.time_col]])-.dt*60/2, 
                                   timep=min(.dfc[[.time_col]])-.dt*60/2,
                                   var=ifelse(.log, 0, -Inf), varp=Inf) %>% 
                          stats::setNames(c(.cell_col, .facet_col, .time_col, paste0(.time_col, 'p'), .var_col, paste0(.var_col, 'p'))) )
        }
      }
      
      if (unique(.dfc[[.facet_col]]) < unique(.dfp[[.facet_col]])) {
        .out <- dplyr::bind_rows(
          dplyr::filter_(dplyr::ungroup(.dfc), sprintf('%s==min(%s)', .time_col, .time_col)) %>%
            dplyr::select_(.cell_col, .facet_col, .time_col, .var_col) %>%
            dplyr::mutate_(.dots=list(lazyeval::interp(~tvar-.dt*60/2, tvar=as.name(.time_col), dt=.dt), "Inf") %>% 
                      stats::setNames(paste0(c(.time_col, .var_col), 'p'))) ,
          dplyr::filter_(dplyr::ungroup(.dfp), sprintf('%s==max(%s)', .time_col, .time_col)) %>%
            dplyr::select_(.cell_col, .facet_col, .time_col, .var_col) %>%
            dplyr::mutate_(.dots=list(lazyeval::interp(~tvar+.dt*60/2, tvar=as.name(.time_col), dt=.dt), as.character(ifelse(.log, 0, -Inf))) %>% 
                      stats::setNames(paste0(c(.time_col, .var_col), 'p')))
        )
        # add vertical line in intermediate facets if required
        if(diff(.out[[.facet_col]]) > 1) {
          .out <- rbind(.out,
                        data.frame(id=unique(.dfc$id), seq(min(.out[[.facet_col]])+1, max(.out[[.facet_col]])-1),
                                   time=min(.dfc[[.time_col]])-.dt*60/2, 
                                   timep=min(.dfc[[.time_col]])-.dt*60/2,
                                   var=ifelse(.log, 0, -Inf), varp=Inf) %>% 
                          stats::setNames(c(.cell_col, .facet_col, .time_col, paste0(.time_col, 'p'), .var_col, paste0(.var_col, 'p'))) )
        }
      }
      return(.out)
    })(.df, .))
  
  .facet_min <- min(.df[[.facet_col]])
  .facet_max <- max(.df[[.facet_col]])
  
  .pl <- ggplot2::ggplot() +
    # facets alternated background
    ggplot2::geom_rect(ggplot2::aes(xmin=-Inf, xmax=Inf, ymin=ifelse(.log, 0, -Inf), ymax=Inf), alpha=.05, 
              data=data.frame(seq(.facet_min, .facet_max, 2)) %>% stats::setNames(.facet_col)) +
    # theme(.gg_theme, complete=TRUE) +
    ggplot2::theme(panel.margin = grid::unit(0, "lines"), panel.border=ggplot2::element_blank()) # ,strip.background = element_blank(), strip.text = element_blank()  
  if (is.null(.facet_labeller)) {
    .pl <- .pl + ggplot2::facet_grid(stats::reformulate('.', .facet_col), as.table=FALSE)
  } else {
    .pl <- .pl + ggplot2::facet_grid(stats::reformulate('.', .facet_col), as.table=FALSE, labeller=ggplot2::as_labeller(.facet_labeller))
  }
  
  if (is.na(.col)) {
    #  use cell id to colour cells
    .pl <- .pl + 
      # show divisions (lty='11' means densely dotted line using hex notation)
      ggplot2::geom_segment(ggplot2::aes_(x=as.name(.time_col), xend=as.name(paste0(.time_col, 'p')), 
                        y=as.name(.var_col), yend=as.name(paste0(.var_col, 'p')), 
                        col=lazyeval::interp(~factor(id), id=as.name(.cell_col)) ), 
                   alpha=.3, lty='11', data=.df_div) +
      # show cell traces
      ggplot2::geom_path(ggplot2::aes_(as.name(.time_col), as.name(.var_col), col=lazyeval::interp(~factor(id), id=as.name(.cell_col))), data=.df) +
      ggCustomTJ::scale_colour_periodic_brewer(guide='none')
    
    # show cell numbers
    if (.show_cellid)
      .pl <- .pl + ggplot2::geom_text(ggplot2::aes_(as.name(.time_col), as.name(.var_col), label=as.name(.cell_col), 
                                  col=lazyeval::interp(~factor(id), id=as.name(.cell_col))), 
                             size=2, hjust=0, vjust=1, data=.df %>% dplyr::group_by_(.cell_col) %>% dplyr::filter(row_number()==1) ) #dplyr::row_number()
    # show all traces in one panel
    if (.show_all)
      .pl <- .pl + ggplot2::geom_path(ggplot2::aes_(as.name(.time_col), as.name(.var_col), col=lazyeval::interp(~factor(id), id=as.name(.cell_col))), 
                             data=.df %>% dplyr::mutate_(.dots=list("-1") %>% stats::setNames(.facet_col)) )
  } else {
    #  case of fixed colour
    .pl <- .pl + 
      # show divisions (lty='11' means densely dotted line using hex notation)
      ggplot2::geom_segment(ggplot2::aes_(x=as.name(.time_col), xend=as.name(paste0(.time_col, 'p')), 
                        y=as.name(.var_col), yend=as.name(paste0(.var_col, 'p'))), 
                   col=.col, alpha=.3, lty='11', data=.df_div) +
      # show cell traces
      ggplot2::geom_path(ggplot2::aes_(as.name(.time_col), as.name(.var_col), group=as.name(.cell_col)), col=.col, data=.df)
    
    # show cell numbers
    if (.show_cellid)
      .pl <- .pl + ggplot2::geom_text(ggplot2::aes_(as.name(.time_col), as.name(.var_col), label=as.name(.cell_col)), size=2, hjust=0, vjust=1,
                             data=.df %>% dplyr::group_by_(.cell_col) %>% dplyr::filter(dplyr::row_number()==1) )
    # show all traces in one panel
    if (.show_all)
      .pl <- .pl + ggplot2::geom_path(ggplot2::aes_(as.name(.time_col), as.name(.var_col), group=as.name(.cell_col)), col=.col, 
                             data=.df %>% dplyr::mutate_(.dots=list("-1") %>% stats::setNames(.facet_col)) )
  }
  
  if(.log) .pl <- .pl + ggplot2::scale_y_continuous(trans='log2')
  return(.pl)
}

