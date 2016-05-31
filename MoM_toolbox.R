mylibs <- c("stats", # required to prevent hiding dplyr::filter later on
            "tools", "stringi", "stringr", "tidyr", "dplyr",
            "ggplot2", "scales", "gridExtra",
            "RcppArmadillo", "ccaPP")
invisible( suppressPackageStartupMessages( # don't use %>% before loading dplyr
  lapply(mylibs, library, character.only=TRUE) ))

# ggplot2 housekeeping ####
theme_set(theme_bw())
scale_colour_discrete <- function(...) scale_colour_brewer(..., palette="Set1", na.value='gray50')
scale_fill_discrete <- function(...) scale_fill_brewer(..., palette="Set1", na.value='gray50')
# to revert to the default ggplot2 discrete colour scale, use: + ggplot2::scale_colour_discrete()
brewer_cols <- c(RColorBrewer::brewer.pal(4, 'Set1'), 'gray42')
scale_colour_periodic_brewer <-
  function(...) scale_colour_manual(..., values = rep(brewer_cols, 1000), na.value='gray25')
scale_fill_periodic_brewer <- function(...) 
  scale_fill_manual(..., values = rep(brewer_cols, 1000), na.value='gray25')
scale_shape_periodic <- function(...) 
  scale_shape_manual(..., values = rep(15:18, 5))

formatter_sec_to_h <- function(.x) .x/3600 #%>% format(digits=2)
scale_x_hours <- function(.dh=6, ...) {
  .f_breaks <- function(.lims) 
    seq(.lims[1] %/% (.dh*3600) *.dh*3600, (.lims[2] %/% (.dh*3600) + 1) *.dh*3600, .dh*3600)
  .call <- as.list(match.call()) # avoid pipe
  if("name" %in% names(.call)) {
    scale_x_continuous(..., labels=formatter_sec_to_h, breaks=.f_breaks)
  } else {
    scale_x_continuous(..., name="time (h)", labels=formatter_sec_to_h, breaks=.f_breaks)
  }
}
scale_y_hours <- function(.dh=6, ...) {
  .f_breaks <- function(.lims) 
    seq(.lims[1] %/% (.dh*3600) *.dh*3600, (.lims[2] %/% (.dh*3600) + 1) *.dh*3600, .dh*3600)
  .call <- as.list(match.call()) # avoid pipe
  if("name" %in% names(.call)) {
    scale_y_continuous(..., labels=formatter_sec_to_h, breaks=.f_breaks)
  } else {
    scale_y_continuous(..., name="time (h)", labels=formatter_sec_to_h, breaks=.f_breaks)
  }
}
sec_to_h_trans <- function() trans_new("sec_to_h", function(.x) .x/3600, function(.x) .x*3600)

cluster_assign_obj <- function(.cl, ...) {
  .xs <- list(...)
  .ns <- as.character(eval(substitute(alist(...))))
  for (.i in 1:length(.xs)) {
    stopifnot(exists(.ns[[.i]]))
    .cl <- cluster_assign_value(.cl, .ns[[.i]], .xs[[.i]])
  }
  return(.cl)
}

cluster_assign_func <- function(.cl, ...) {
  .fs <- list(...)
  .ns <- as.character(eval(substitute(alist(...))))
  for (.i in 1:length(.fs)) {
    stopifnot(is.function(.fs[[.i]]))
    .cl <- cluster_assign_value(.cl, .ns[[.i]], .fs[[.i]])
  }
  return(.cl)
}

# general convenience functions ####
between_or <- function(.x, .a, .b) {
  if (length(.a) != length(.b)) error('.a and .b must have the same length.')
  lapply(seq(length(.a)), function(.i) between(.x, .a[.i], .b[.i])) %>%
    do.call(cbind, .) %>%
    apply(1, any)
}

value_occurence_index <- function(.m) {
  # value_occurence_index returns the index of each occurence of a vector level
  .m <- factor(.m)
  .out <- NA
  for (.l in levels(.m)) {
    .idx <- which(.m==.l)
    .out[.idx] <- 1:length(.idx)
  }
  return(.out)
}

find_unique_interval <- function(.xs, .mins, .maxs) {
# find_unique_interval() returns the indices of the interval (defined by [.mins, .maxs]) matching 
# each value of .xs (provided it is unique)
  if (length(.mins) != length(.maxs)) stop(".mins and .maxs have different lengths!")
  .out <-  lapply(seq(length(.mins)), function(.i) between(.xs, .mins[.i], .maxs[.i])) %>% 
    do.call(cbind, .) %>%
    apply(1, which)
  .ls <- sapply(.out, length) 
  if (any(.ls > 1)) warning('some .xs values are matched by several intervals.')
  .out[.ls!=1] <- NA
  return(unlist(.out))
}

as_numeric <- function(.x) as.numeric(as.character(.x))
as_numeric_df <- function(.d) modifyList(.d, lapply(.d, as_numeric))

mycut <- function(.x, ...) {
  .k <- Hmisc::cut2(.x, ...)
  .mids <- sapply(str_match_all(as.character(levels(.k)), "[\\[\\(](.*),(.*)[\\]\\)]"), 
                 function(.m) {
                   mean(c(as.numeric(.m[2]), as.numeric(.m[3])))
                 })
  .mids[.k]
}

# MoMA loading ####
load_timm_data_ini <- function(.path, .scripts_path, 
                           .perl_cmd='perl',
                           .frames_script="get_size_and_fluo_basic.pl", 
                           # .cells_script="estimate_division_growth.pl", 
                           .force=FALSE, .verbose=FALSE,
                           .data2preproc=identity) {
  if (!file.exists(.path)) stop("input path is not valid.")
  .dir <- dirname(.path)
  .name <- basename(.path)
  .frames_script_path <- file.path(.scripts_path, .frames_script)
  # .cells_script_path <- file.path(.scripts_path, .cells_script)
  .outdir <- .dir %>% .data2preproc
  dir.create(.outdir, recursive=TRUE, showWarnings=FALSE)
  .outbase <- basename(.path) %>% sub("ExportedCellStats_", "", .) %>% file_path_sans_ext
  .frames_path <- paste0(.outbase, "_frames.txt") %>% file.path(.outdir, .)
  # .cells_path <- paste0(.outbase, "_cells.txt") %>% file.path(.outdir, .)

  # if (!file.exists(.frames_path))
  #   print(.frames_path)
  # return(list())
  
  # run perl scripts if output doesn't exist
  if (!file.exists(.frames_path) || .force==TRUE) { # || !file.exists(.cells_path) 
    if (.verbose) print(paste("Converting", .path))
    # erik's scripts work only in the working dir
    if (!file.copy(.path, .name, overwrite=TRUE, copy.mode=FALSE)) stop(paste(.path, 'cannot be copied.'))
    .out_frames_path <- system2(.perl_cmd, args=c(.frames_script_path, .name), stdout=TRUE) %>% sub('^>', '', .)
    if (!file.exists(.out_frames_path)) stop("Frames file cannot be found.")
  
#     if (.verbose) print(paste("Converting", .out_frames_path))
#     .out_cells_path <- system2(.perl_cmd, args=c(.cells_script_path, .out_frames_path), stdout=TRUE) %>% sub('^>', '', .)
    
    # housekeeping
    file.rename(.out_frames_path, .frames_path)
    # file.rename(.out_cells_path, .cells_path)
    file.remove(.name, .out_frames_path) #, .out_cells_path
  }

  .frames_out <- parse_frames_stats(.frames_path) %>%
    rename(id=cell_ID, parent_id=parent_ID, end_type=type_of_end) %>%
    mutate(cid=compute_genealogy(id, parent_id, daughter_type) )
    
#   .cells_out <- read.table(.cells_path, comment.char="", header=TRUE) %>%
#     select(-1)
#   names(.cells_out) <- names(.cells_out) %>% 
#     gsub("X\\.", "", .) %>%
#     gsub("\\.", "_", .)
#   .cells_out <- .cells_out %>%
#     rename(id=cellID, parent_id=parID) %>%
#     select(-daughter_type) %>%
#     inner_join(select(.frames_out, id, cid) %>% group_by(id) %>% slice(1), by='id')
#     
#   .frames_out <- .frames_out %>%
#     left_join(select(.cells_out, -parent_id, -end_type, -cid), by='id')
  
  # if (dim(.cells_out)[1] > 0) {
    .m <- str_match(basename(.path), ".*(\\d{8})_.*pos(\\d+).*_GL(\\d+).*")
    return( list(#cells=data.frame(date=as.numeric(.m[1, 2]), pos=as.numeric(.m[1, 3]), gl=as.numeric(.m[1, 4]), .cells_out),
                 frames=data.frame(date=as.numeric(.m[1, 2]), pos=as.numeric(.m[1, 3]), gl=as.numeric(.m[1, 4]), .frames_out)) )
  # }
}


if (!exists("data2preproc_dir", mode="function"))
  data2preproc_dir <- identity
if (!exists("data2preproc_file", mode="function"))
  data2preproc_file <- identity
if (!exists("data2preproc", mode="function"))
  data2preproc <- function(.path)
    file.path(data2preproc_dir(.path), data2preproc_file(.path))

process_moma_data <- function(.x, .data2preproc, .scripts_path, .force=FALSE,
                              .frames_sh_script="get_size_and_fluo_basic.sh", .frames_pl_script="get_size_and_fluo_basic.pl",
                              .qsub_name="MM_pl" # must be shorter than 10 characters
) { # browser()
  # check whether some files need to be preprocessed
  .preprocessed <- .data2preproc(.x) %>% file.exists
  if (all(.preprocessed) && !.force) {
    return(.x) 
  } else {
    # call qsub
    lapply(.x[ which(.force | !.preprocessed) ],
           function(.f) system(sprintf("source /etc/profile.d/sge.sh; qsub -N %s -cwd %s %s %s %s", 
                                       .qsub_name, file.path(.scripts_path, .frames_sh_script), 
                                       file.path(.scripts_path, .frames_pl_script), .f, data2preproc(.f)))) 
    stop("Some files required preprocessing on the cluster. Once they're processed, run the same command again\n",
         sprintf("Hint: use `process_state(\"%s\")` to check if the processing is still running...", .qsub_name),
         call.=FALSE)
  }
}

process_state <- function(.qsub_name="MM_pl") {
  .qs <- system("source /etc/profile.d/sge.sh; qstat", intern=TRUE) 
  
  if (length(.qs) == 0) {
    return("All jobs have been processed.")
  } else {
    paste(.qs[-2], collapse='\n') %>% read_table() %>%
      filter(name==.qsub_name) %>%
      group_by(state) %>% summarise(n_jobs=n()) %>% 
      rowwise() %>% 
      mutate(msg=sprintf('%d job(s) in state %s', n_jobs, state)) %>% 
      .[['msg']] %>% paste(collapse='\n') %>% 
      ifelse(nchar(.)==0, 'All jobs have been processed.', .)
  }
}

load_moma_processed <- function(.path, .verbose=FALSE) {
  if (!file.exists(.path)) stop("input path is not valid.")

  .frames_out <- parse_frames_stats(.path) %>%
    rename(id=cell_ID, parent_id=parent_ID, end_type=type_of_end) %>%
    mutate(cid=compute_genealogy(id, parent_id, daughter_type) )
  .m <- str_match(basename(.path), ".*(\\d{8})_.*pos(\\d+).*_GL(\\d+).*")
  
  data.frame(date=as.numeric(.m[1, 2]), pos=as.numeric(.m[1, 3]), gl=as.numeric(.m[1, 4]), 
             .frames_out)
}


parse_frames_stats <- function(.path) {
  flines <- readLines(.path)
  # browser()

  # find column names for stats
  .fieldnames <- str_split(flines[1], '[ \t]+')[[1]] %>% sub("#", "", .) %>%
    gsub(" ", "", .) %>%
    gsub("-", "_", .) %>%
    gsub("\\(", "_", .) %>%
    gsub("\\)", "", .)
  .colnames <- str_split(flines[2], '[ \t]+')[[1]] %>% sub("#", "", .) %>%
    gsub(" ", "", .) %>%
    gsub("-", "_", .) %>%
    gsub("\\(", "_", .) %>%
    gsub("\\)", "", .) %>%
    gsub("botom", "bottom", .)
  
  # parse all cells
  id_lines <- grep("^>CELL", flines)
  .out <- do.call(rbind, 
                  lapply(seq_along(id_lines)[-1], function(.i) {
                    # browser()
#                     .df <- (id_lines[.i-1] + 1):(id_lines[.i] - 1) %>% # retrieve lines indices
#                       flines[.] %>% paste(collapse='\n') %>%           # concat lines
#                       textConnection %>% read.table(comment.char="", header=FALSE) %>% # read table
#                       setNames(.colnames)
                    .str <- (id_lines[.i-1] + 1):(id_lines[.i] - 1) %>% # retrieve lines indices
                      flines[.] %>% paste(collapse='\n') # concat lines
                    .con <- textConnection(.str)
                    .df <- read.table(.con, comment.char="", header=FALSE) %>% setNames(.colnames)
                    close(.con) # close text connection
                    .info <- flines[id_lines[.i-1]] %>% # extract line
                      strsplit('[ \t]+') %>% .[[1]] %>%     # split to vector
                      .[-1] %>% t %>% data.frame(stringsAsFactors=FALSE) %>%    # convert to df
                      setNames(.fieldnames[-1]) %>%
                      mutate(cell_ID=as.numeric(cell_ID), parent_ID=as.numeric(parent_ID), 
                             first_frame=as.numeric(first_frame), last_frame=as.numeric(last_frame))
                    cbind(.info, .df)
                  }) )
  return(.out)
}


which_touch_exit <- function(.h, .hmin_cutoff) {
  # browser()
  .t <- which(.h<.hmin_cutoff)
  if (length(.t) == 0) {
    return(rep(FALSE, length(.h)))
  } else {
    return(c(rep(FALSE, .t[1]-1),
             rep(TRUE, length(.h)-.t[1]+1)))
  }
}

which_to_progeny <- function(.x, .cid) {
  #   browser()
  .df <- data.frame(x=.x, cid=.cid)
  .out <- .x
  .cs <- unique(.cid[which(.x)])
  for (.c in .cs) {
    .is_daughter <- str_detect(.cid, paste0(.c, ".+"))
    .out[.is_daughter] <- TRUE
  }
  return(.out)  
}


parse_timm_curation <- function(.path) {
  .flines <- readLines(.path)
  .flines <- .flines[6:length(.flines)] %>%
    gsub('\t', '', .) # remove tabs
  # read one table for each type and store all tables in a list
  .id_lines <- grep("^#", .flines)
  .ldata <- lapply(seq(length(.id_lines)), function(.i) {
    .start <- .id_lines[.i] + 1
    .stop <- ifelse(.i==length(.id_lines), length(.flines), .id_lines[.i+1] - 1)
    if (.start <= .stop)
      return ( .flines[.start:.stop] %>% paste(collapse='\n') %>%                             # concat lines
                 textConnection %>% read.table(comment.char="", sep=",", header=FALSE, stringsAsFactors=FALSE) )
  })
  if (all(sapply(.ldata, is.null))) return(data.frame(type='none', frame=-1, action=NA))
  # find the width of the largest table
  .max_cols <- lapply(.ldata, function(.d) dim(.d)[2]) %>% unlist %>% max
  # increase the width of narrower dfs (in order to be able to stack them)
  lapply(.ldata, function(.d) {
    if (is.null(.d) || dim(.d)[2] == .max_cols) return(.d)
    .d[1, .max_cols] <- NA # add NAs columns if narrower dataframe
    .d
  }) %>% 
    # stack all tables
    do.call(rbind, .) %>%
    setNames(c('type', 'frame', 'xx', 'action'))
}


# genealogy utilities ####
compute_genealogy <- function(.id, .pid, .dgtype) {
# browser()
  # keep only one row per cell
  .all_id <- .id
  .idx <- sapply(unique(.id), function(.i) min(which(.id==.i)))
  .id <- .id[.idx]
  .pid <- .pid[.idx]
  .dgtype <- .dgtype[.idx]
  
  # check that all cells are there
  if (sort(unique(.id), decreasing=TRUE)[1] != length(unique(.id)) - 1) 
    stop('compute_genealogy() requires that each cell appear at least once.')
  
  .dgtype <- substr(.dgtype, 1, 1)
  .id <- .id + 1
  .pid <- .pid + 1
  .cid <- character(length(.id))
  .founders <- which(.pid==0)
  
  for (.f in .founders) {
    .cid[.f] <- .id[.f] - 1 # use TIMM index for founder identification
    .daughters <- which(.pid==.id[.f])
    while(length(.daughters) > 0) {
      .cid[.daughters] <- paste0(.cid[.pid[.daughters]], .dgtype[.daughters])
      .daughters <- which(.pid %in% .id[.daughters])
    }
  }
  return(.cid[.all_id+1])
}

get_parent_cid <- function(.cid) {
  stri_sub(.cid, to=-2)
}

get_all_parents_cid <- function(.cid) {
  stri_sub(.cid, to=-(2:nchar(.cid))) 
}

get_daughters_cid <- function(.cid) {
  paste0(.cid, c('B', 'T'))
}

genealogy_ontology <- function(.div_min, .div_max) {
  .rel <- NA
  if (.div_min == 0 && .div_max == 1) {
    .rel <- 'daughter'
  } else if (.div_min == 0 && .div_max == 2) {
    .rel <- 'gd-daughter'
  } else if (.div_min == 0 && .div_max == 3) {
    .rel <- 'gd-gd-daughter'
  } else if (.div_min == 1 && .div_max == 1) {
    .rel <- 'sisters'
  } else if (.div_min == 2 && .div_max == 2) {
    .rel <- 'cousins1'
  } else if (.div_min == 3 && .div_max == 3) {
    .rel <- 'cousins2'
  } else if (.div_min == 1 && .div_max == 2) {
    .rel <- 'niece'
  }
  return(.rel)
}

convert_genealogy_to_poleline <- function(.s, .bottom_char='B', .top_char='T') {
# convert_genealogy_to_poleline() converts a B/T genealogy string (bottom/top) 
# to an O/N one (old pole/new pole cell). It relies on the fact that the old pole cell is 
# on the same side during two consecutive divisions, e.g. for 0BTTBTTT:
#  (0)BTTBTT(T)  B/T genealogy
# (0B)TTBTTT     same, shifted by one letter
#   0BNONNOO     result (O for identical letters, N for different ones)
  .ss <- sprintf('[%s%s]+', .bottom_char, .top_char) %>% str_match(.s, .) %>% strsplit('')
  .p <- lapply(.ss, function(.x) (.x[-length(.x)] == .x[-1]) %>% ifelse('O', 'N') %>% paste(collapse='')) %>% unlist
  
  .pos <- sprintf('[%s%s]+', .bottom_char, .top_char) %>% str_locate(.s, .) %>% .[,'start']
  paste0(str_sub(.s, 0, .pos), .p)
}

get_pole_age <- function(.s, .old_char='O', .new_char='N') {
# get_pole_age() computes the number of consecutive identical letters at the end of the string
# for old pole cells, the age is positive (number of consecutive division inheriting the old pole)
# for new pole cells, the age is negative (number of consecutive division inheriting the youngest old pole)
  .bs <- paste0(.old_char, '+$') %>% str_match(.s, .) %>% str_length
  .ts <- paste0(.new_char, '+$') %>% str_match(.s, .) %>% str_length
  return( ifelse(is.na(.bs), -.ts, .bs) )
}

get_divs_from_mothercell <- function(.s, .old_char='B', .new_char='T') {
# get_divs_from_mothercell() computes the number of divisions since the cell's ancestor was the mother cell of the GL
  .m_length <- sprintf('[^%s%s](%s+)', .old_char, .new_char, .old_char) %>% str_match(.s, .) %>% .[,2] %>% str_length
  .s_length <- sprintf('[%s%s]+', .old_char, .new_char) %>% str_match(.s, .) %>% str_length
  return( .s_length - .m_length )
}

genealogy_relationship <- function(.s, .dist_max=3) {
  # sisters: div_min == div_max == 1
  # 1st cousins: div_min == div_max == 2
  # niece: div_min==1 && div_max==2

  # browser()
  if (length(.s)!=2) stop("the argument must be a vector of 2 strings.")
  .l1 <- str_length(.s[1])
  .l2 <- str_length(.s[2])
  if (abs(.l2 - .l1) > .dist_max) 
    return(NULL)
  
  .min_l <- min(.l1, .l2)
  .ss1 <- str_sub(.s[1], 1, .min_l)
  .ss2 <- str_sub(.s[2], 1, .min_l)
  if (.ss1 == .ss2) {
    .anc_length <- .min_l
  } else {
  .x1 <- .ss1 %>% str_split("") %>% unlist
  .x2 <- .ss2 %>% str_split("") %>% unlist
  .anc_length <- min(which(.x1 != .x2)) - 1 # 0 if no common ancestor
  }
  if (.min_l - .anc_length > .dist_max)
    return(NULL)
  .div_1 <- str_sub(.s[1], .anc_length+1, -1) %>% str_length # -1 to select the last char
  .div_2 <- str_sub(.s[2], .anc_length+1, -1) %>% str_length # -1 to select the last char
  .div_min <- min(.div_1, .div_2)
  .div_max <- max(.div_1, .div_2)
  if (.div_max > .dist_max) 
    return(NULL)
  return(list(rel=genealogy_ontology(.div_min, .div_max), ancestor=str_sub(.s[1], 1, .anc_length), 
              genealogy_1=.s[1], div_1=.div_1, genealogy_2=.s[2], div_2=.div_2, 
              div_min=.div_min, div_max=.div_max, less_div=which.min(c(.div_1, .div_2)) ))
}


# plotting tracks ####

plot_faceted_var_tracks <- function(.df, .var_col='gfp_nb', .time_col='time_sec', 
                                    .cell_col='id', .parent_col='parent_id', .facet_col='b_rank', 
                                    .log=FALSE, .col=NA, .show_cellid=FALSE, .show_all=FALSE, 
                                    .facet_labeller=NULL, .gg_theme=theme_get()) {
  
  .df_div <- .df %>% group_by_(.cell_col) %>% 
    do((function(.dfgl, .dfc) { # compute_div_between_facets
      # .dfgl: dataframe of the entire GL
      # .dfc: dataframe of the current cell
      if (unique(.dfc[[.parent_col]]) < 0) return(data.frame())
      .dfp <- filter_(.dfgl, lazyeval::interp(~ cell_id == unique(parent_id), 
                                      cell_id=as.name(.cell_col), parent_id=.dfc[[.parent_col]]))
      if (dim(.dfp)[1] == 0) return(data.frame())
      
      if (unique(.dfc[[.facet_col]]) == unique(.dfp[[.facet_col]])) {
        .varp <- filter_(.dfp, sprintf('%s==max(%s)', .time_col, .time_col))[[.var_col]]
        .out <- filter_(.dfc, sprintf('%s==min(%s)', .time_col, .time_col)) %>% 
          select_(.cell_col, .facet_col, .time_col, .var_col) %>%
          mutate_(.dots=list(lazyeval::interp(~tvar-dt*60, tvar=as.name(.time_col), dt=dt), 
                             lazyeval::interp(~.varp, .varp=.varp)) %>% 
                    setNames(paste0(c(.time_col, .var_col), 'p')))
      }
      if (unique(.dfc[[.facet_col]]) > unique(.dfp[[.facet_col]])) {
        .out <- bind_rows(
          filter_(ungroup(.dfc), sprintf('%s==min(%s)', .time_col, .time_col)) %>%
            select_(.cell_col, .facet_col, .time_col, .var_col) %>%
            mutate_(.dots=list(lazyeval::interp(~tvar-dt*60/2, tvar=as.name(.time_col), dt=dt), as.character(ifelse(.log, 0, -Inf))) %>% 
                      setNames(paste0(c(.time_col, .var_col), 'p'))) ,
          filter_(ungroup(.dfp), sprintf('%s==max(%s)', .time_col, .time_col)) %>%
            select_(.cell_col, .facet_col, .time_col, .var_col) %>%
            mutate_(.dots=list(lazyeval::interp(~tvar+dt*60/2, tvar=as.name(.time_col), dt=dt), "Inf") %>% 
                      setNames(paste0(c(.time_col, .var_col), 'p')))
        )
        if(diff(.out[[.facet_col]]) < -1) {
          .out <- rbind(.out,
                        data.frame(id=unique(.dfc$id), seq(min(.out[[.facet_col]])+1, max(.out[[.facet_col]])-1),
                                   time=min(.dfc[[.time_col]])-dt*60/2, 
                                   timep=min(.dfc[[.time_col]])-dt*60/2,
                                   var=ifelse(.log, 0, -Inf), varp=Inf) %>% 
                          setNames(c(.cell_col, .facet_col, .time_col, paste0(.time_col, 'p'), .var_col, paste0(.var_col, 'p'))) )
        }
      }
      return(.out)
    })(.df, .))

  .facet_min <- min(.df[[.facet_col]])
  .facet_max <- max(.df[[.facet_col]])
  
  .pl <- ggplot() +
    # facets alternated background
    geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=ifelse(.log, 0, -Inf), ymax=Inf), alpha=.05, 
              data=data.frame(seq(.facet_min, .facet_max, 2)) %>% setNames(.facet_col)) +
    # theme(.gg_theme, complete=TRUE) +
    theme(panel.margin = unit(0, "lines"), panel.border=element_blank()) # ,strip.background = element_blank(), strip.text = element_blank()
  
  if (is.null(.facet_labeller)) {
    .pl <- .pl + facet_grid(reformulate('.', .facet_col), as.table=FALSE)
  } else {
    .pl <- .pl + facet_grid(reformulate('.', .facet_col), as.table=FALSE, labeller=as_labeller(.facet_labeller))
  }
  
  if (is.na(.col)) {
    #  use cell id to colour cells
    .pl <- .pl + 
      # show divisions (lty='11' means densely dotted line using hex notation)
      geom_segment(aes_(x=as.name(.time_col), xend=as.name(paste0(.time_col, 'p')), 
                        y=as.name(.var_col), yend=as.name(paste0(.var_col, 'p')), 
                        col=lazyeval::interp(~factor(id), id=as.name(.cell_col)) ), 
                   alpha=.3, lty='11', data=.df_div) +
      # show cell traces
      geom_path(aes_(as.name(.time_col), as.name(.var_col), col=lazyeval::interp(~factor(id), id=as.name(.cell_col))), data=.df) +
      scale_colour_periodic_brewer(guide='none')
    
    # show cell numbers
    if (.show_cellid)
      .pl <- .pl + geom_text(aes_(as.name(.time_col), as.name(.var_col), label=as.name(.cell_col), 
                                  col=lazyeval::interp(~factor(id), id=as.name(.cell_col))), 
                             size=2, hjust=0, vjust=1, data=.df %>% group_by_(.cell_col) %>% filter(row_number()==1) )
    # show all traces in one panel
    if (.show_all)
      .pl <- .pl + geom_path(aes_(as.name(.time_col), as.name(.var_col), col=lazyeval::interp(~factor(id), id=as.name(.cell_col))), 
                             data=.df %>% mutate_(.dots=list("-1") %>% setNames(.facet_col)) )
  } else {
    #  case of fixed colour
    .pl <- .pl + 
      # show divisions (lty='11' means densely dotted line using hex notation)
      geom_segment(aes_(x=as.name(.time_col), xend=as.name(paste0(.time_col, 'p')), 
                        y=as.name(.var_col), yend=as.name(paste0(.var_col, 'p'))), 
                   col=.col, alpha=.3, lty='11', data=.df_div) +
      # show cell traces
      geom_path(aes_(as.name(.time_col), as.name(.var_col), group=as.name(.cell_col)), col=.col, data=.df)
    
    # show cell numbers
    if (.show_cellid)
      .pl <- .pl + geom_text(aes_(as.name(.time_col), as.name(.var_col), label=as.name(.cell_col)), size=2, hjust=0, vjust=1,
                             data=.df %>% group_by_(.cell_col) %>% filter(row_number()==1) )
    # show all traces in one panel
    if (.show_all)
      .pl <- .pl + geom_path(aes_(as.name(.time_col), as.name(.var_col), group=as.name(.cell_col)), col=.col, 
                             data=.df %>% mutate_(.dots=list("-1") %>% setNames(.facet_col)) )
  }
  
  if(.log) .pl <- .pl + scale_y_log10()
  return(.pl)
}


