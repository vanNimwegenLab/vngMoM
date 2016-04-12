mylibs <- c("stats", # required to prevent hiding dplyr::filter later on
            "tools", "stringi", "stringr", "tidyr", "dplyr",
            "ggplot2", "scales", "RColorBrewer", "gridExtra",
            "RcppArmadillo", "ccaPP")
invisible( suppressPackageStartupMessages( # don't use %>% before loading dplyr
  lapply(mylibs, library, character.only=TRUE) ))

# ggplot2 housekeeping ####
theme_set(theme_bw())
scale_colour_discrete <- function(...) scale_colour_brewer(..., palette="Set1", na.value='gray50')
scale_fill_discrete <- function(...) scale_fill_brewer(..., palette="Set1", na.value='gray50')
# to revert to the default ggplot2 discrete colour scale, use: + ggplot2::scale_colour_discrete()
brewer_cols <- c(brewer.pal(4, 'Set1'), 'gray42')
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

# multidplyr utility ####
# cluster_assign_obj <- function(cl, x) {
#   name <- as.character(substitute(x))
#   stopifnot(exists(name))
#   cluster_assign_value(cl, name, x)
# }
cluster_assign_obj <- function(.cl, ...) {
  .xs <- list(...)
  .ns <- as.character(eval(substitute(alist(...))))
  for (.i in 1:length(.xs)) {
    stopifnot(exists(.ns[[.i]]))
    .cl <- cluster_assign_value(.cl, .ns[[.i]], .xs[[.i]])
  }
  return(.cl)
}

# cluster_assign_func <- function(cl, f) {
#   fname <- as.character(substitute(f))
#   stopifnot(is.function(f))
#   cluster_assign_value(cl, fname, f)
# }
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
load_timm_data <- function(.path, .scripts_path, 
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


