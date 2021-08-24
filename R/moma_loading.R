# MoMA loading ####
utils::globalVariables(c(".", "msg", "name", "n_jobs", "state",
                         "n",
                         "cell_ID", "parent_ID", "first_frame", "last_frame", "type_of_end", "daughter_type",
                         "id", "parent_id", "end_type", "cid", "fluo_background", "fluo_amplitude", "fluo_bg_ch_1", "fluo_ampl_ch_1",
                         "ndgt", "end_type_moma"))

process_moma_data <- function(.x, .data2preproc, .scripts_path=system.file("perl", ".", package="vngMoM"), .force=FALSE, .skip=FALSE,
                              .frames_sh_script="get_size_and_fluo.sh", .frames_pl_script="get_size_and_fluo_basic.pl",
                              .qsub_name="MM_pl" # must be shorter than 10 characters
) { 
  # browser()
  # check whether some files need to be preprocessed
  .out <- .data2preproc(.x)
  .preprocessed <- file.exists(.out)
  if (all(.preprocessed) && !.force) {
    return(.out) 
  } else if (.skip) {
    # return already preprocessed files
    return(ifelse(.preprocessed, .out, NA)) 
  } else {
    # call qsub
    lapply(.x[ which(.force | !.preprocessed) ],
           function(.f) system(sprintf("sbatch --job-name=%s %s %s %s %s", 
                                       .qsub_name, file.path(.scripts_path, .frames_sh_script), 
                                       file.path(.scripts_path, .frames_pl_script), .f, .data2preproc(.f)))) 
    stop("Some files required preprocessing on the cluster. Once they're processed, run the same command again\n",
         sprintf("Hint: use `process_state(\"%s\")` to check if the processing is still running...", .qsub_name),
         call.=FALSE)
  }
}

squeue <- function(.user=Sys.info()[['user']], .args="", .intern=FALSE)
  paste("squeue -u", .user, .args) %>% system(.intern)

sacct <- function(.jobid, .args="")
  paste("sacct -j", .jobid, .args) %>% system

process_state <- function(.qsub_name="MM_pl") {
  .qs <- squeue(.intern=TRUE) 
  
  if (length(.qs) == 1) {
    return("All jobs have been processed.")
  } else {
    .qsn <- .qs[1] %>% stringr::str_to_lower() %>% 
      stringr::str_split("\\s+") %>% 
      (function(.x) .x[[1]][-length(.x[[1]])])
    paste(.qs[-1], collapse='\n') %>% 
      readr::read_fwf(readr::fwf_widths(c(10, 20, 20, 10, 10, 10, 10, 10, 20), col_names=.qsn) ) %>% 
      dplyr::filter(name==.qsub_name) %>%
      dplyr::group_by(state) %>% dplyr::summarise(n_jobs=n()) %>% 
      dplyr::rowwise() %>% 
      dplyr::mutate(msg=sprintf('%d job(s) in state %s', n_jobs, state)) %>% 
      .[['msg']] %>% paste(collapse='\n') %>% 
      ifelse(nchar(.)==0, 'No job currently being processed.', .)
  }
}

parse_frames_stats <- function(.path) {
  if (!file.exists(.path)) stop(sprintf("input path is not valid (%s)", .path))
  flines <- readLines(.path)
  # browser()
  
  # find column names for stats
  .fieldnames <- stringr::str_split(flines[1], '[ \t]+')[[1]] %>% sub("#", "", .) %>%
    gsub(" ", "", .) %>%
    gsub("-", "_", .) %>%
    gsub("\\(", "_", .) %>%
    gsub("\\)", "", .)
  .colnames <- stringr::str_split(flines[2], '[ \t]+')[[1]] %>% sub("#", "", .) %>%
    gsub(" ", "", .) %>%
    gsub("-", "_", .) %>%
    gsub("\\(", "_", .) %>%
    gsub("\\)", "", .) %>%
    gsub("botom", "bottom", .)
  
  # parse all cells
  id_lines <- grep("^>CELL", flines)
  id_lines <- c(id_lines, length(flines)+1) # add the last line index (+1 since readLines skips the last empty line)
  lapply(seq_along(id_lines)[-1], function(.i) {
    # browser()
    #                     .df <- (id_lines[.i-1] + 1):(id_lines[.i] - 1) %>% # retrieve lines indices
    #                       flines[.] %>% paste(collapse='\n') %>%           # concat lines
    #                       textConnection %>% read.table(comment.char="", header=FALSE) %>% # read table
    #                       stats::setNames(.colnames)

    .info <- flines[id_lines[.i-1]] %>% # extract line
      strsplit('[ \t]') %>% .[[1]] %>%     # split to vector
      # strsplit('[ \t]+') %>% .[[1]] %>%     # split to vector
      .[-1] %>% t %>% data.frame(stringsAsFactors=FALSE) %>%    # convert to df
      stats::setNames(.fieldnames[-1]) %>%
      dplyr::mutate(cell_ID=as.numeric(cell_ID), parent_ID=as.numeric(parent_ID), 
                    first_frame=as.numeric(first_frame), last_frame=as.numeric(last_frame))
    
    .str <- (id_lines[.i-1] + 1):(id_lines[.i] - 1) %>% # retrieve lines indices
      flines[.] %>% paste(collapse='\n') # concat lines
    if (stringr::str_detect(.str, "^\\d+$")) 
      return(.info)
    else {
      .con <- textConnection(.str)
      .df <- utils::read.table(.con, comment.char="", header=FALSE, sep="\t", col.names=.colnames)
      close(.con) # close text connection
      cbind(.info, .df)
    }
  }) %>% 
    dplyr::bind_rows() %>% 
    dplyr::rename(id=cell_ID, parent_id=parent_ID, end_type=type_of_end) %>%
    dplyr::mutate(cid=compute_genealogy(id, parent_id, daughter_type),
                  fluo_background=if (exists('fluo_bg_ch_1', where=.)) fluo_bg_ch_1 else fluo_background,
                  fluo_amplitude=if (exists('fluo_ampl_ch_1', where=.)) fluo_ampl_ch_1 else fluo_amplitude ) %>% 
    # fix end_type for pruned cells
    dplyr::mutate(ndgt=compute_daughters_numbers(cid)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(end_type_moma=end_type,
           end_type=ifelse(ndgt==0, "lost", "weird"),
           end_type=ifelse(ndgt==2, "div", end_type))
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
    .is_daughter <- stringr::str_detect(.cid, paste0(.c, ".+"))
    .out[.is_daughter] <- TRUE
  }
  return(.out)  
}

fit_exp_elongation <- function(.t, .l) {
  if (length(.t) != length(.l)) stop('variables of different lengths.')
  if (length(.t) <= 4) return(NA)
  
  stats::lm(log(.l)~.t) %>% # use lm() for predict with se
    stats::predict(se.fit=TRUE) %>% .[['fit']] %>% exp
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
    if (.start <= .stop) {
      .con <- .flines[.start:.stop] %>% paste(collapse='\n') %>% textConnection
      .out <- utils::read.table(.con, comment.char="", sep=",", header=FALSE, stringsAsFactors=FALSE)
      close(.con)
      return (.out)
    }
  })
  if (all(sapply(.ldata, is.null))) return(data.frame(type='none', frame=-1))
  # increase the width of narrower dfs (in order to be able to stack them)
  lapply(.ldata, function(.d) {
    if (is.null(.d)) return(.d)
    return(.d[, 1:2]) # return only `type` and `frame` (avoid handling a variable number of columns)
  }) %>% 
    # stack all tables
    do.call(rbind, .) %>%
    stats::setNames(c('type', 'frame'))
}
