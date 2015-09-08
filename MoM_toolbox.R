suppressPackageStartupMessages( {
  library(Hmisc)
  library(tools)
  library(stringr)
  # library(tidyr)
  library(dplyr)
  
  library(ggplot2)
  library(scales)
  library(RColorBrewer)
} )

theme_set(theme_bw())
scale_colour_discrete <- function(...) scale_colour_brewer(..., palette="Set1")
# to revert to the default ggplot2 discrete colour scale, use: + ggplot2::scale_colour_discrete()
scale_colour_periodic_brewer <-
  function(...) scale_colour_manual(..., values = rep(c(brewer.pal(4, 'Set1'), 'gray42'), 100))
scale_fill_periodic_brewer <-
  function(...) scale_fill_manual(..., values = rep(c(brewer.pal(4, 'Set1'), 'gray42'), 100))
scale_shape_periodic <- 
  function(...) scale_shape_manual(..., values = rep(15:18, 5))


load_timm_data <- function(.path, .scripts_path, 
                           .perl_cmd='perl',
                           .frames_script="get_size_and_fluo.pl", 
                           .cells_script="estimate_division_growth.pl", 
                           .force=FALSE, .verbose=FALSE,
                           .data2preproc=identity) {
  if (!file.exists(.path)) stop("input path is not valid.")
  .dir <- dirname(.path)
  .frames_script_path <- file.path(.scripts_path, .frames_script)
  .cells_script_path <- file.path(.scripts_path, .cells_script)
  .outdir <- .dir %>% .data2preproc
  dir.create(.outdir, recursive=TRUE, showWarnings=FALSE)
  .outbase <- basename(.path) %>% sub("ExportedCellStats_", "", .) %>% file_path_sans_ext
  .out_frames_path <- paste0(.outbase, "_frames.txt") %>% file.path(.outdir, .)
  .out_cells_path <- paste0(.outbase, "_cells.txt") %>% file.path(.outdir, .)
#   browser()
  
  # run perl scripts if output doesn't exist
  if (!file.exists(.out_frames_path) || .force==TRUE) {
    if (.verbose) print(paste("Converting", .path, "to", .out_frames_path))
    system2(.perl_cmd, args=c(.frames_script_path, .path, "0"), stdout=.out_frames_path)
  }
  if (!file.exists(.out_frames_path)) stop("Frames file cannot be found.")

  if (!file.exists(.out_cells_path) || .force==TRUE) {
    if (.verbose) print(paste("Converting", .out_frames_path, "to", .out_cells_path))
    system2(.perl_cmd, args=c(.cells_script_path, .out_frames_path), stdout=.out_cells_path)
  }
  
  .cells_out <- parse_cells_stats(.out_cells_path)$growth
  names(.cells_out) <- names(.cells_out) %>% 
    gsub("X\\.", "", .) %>%
    gsub("\\.", "_", .)
  .cells_out <- .cells_out %>%
    rename(id=cellID, parent_id=parID) %>%
    mutate(cid = compute_genealogy(id, parent_id, daughter_type))
  
  .frames_out <- parse_frames_stats(.out_frames_path) %>%
    select(-daughter_type, -type_of_end) %>%
    rename(id=cell_ID, parent_id=parent_ID) %>%
    left_join(.cells_out)
  
  if (dim(.cells_out)[1] > 0) {
    .m <- str_match(basename(.path), ".*(\\d{8})_.*pos(\\d+).*_GL(\\d+).*")
    return(list(cells=data.frame(date=as.numeric(.m[1, 2]), pos=as.numeric(.m[1, 3]), gl=as.numeric(.m[1, 4]), .cells_out),
                frames=data.frame(date=as.numeric(.m[1, 2]), pos=as.numeric(.m[1, 3]), gl=as.numeric(.m[1, 4]), .frames_out) ))
  }
}

as_numeric <- function(.x) as.numeric(as.character(.x))
as_numeric_df <- function(.d) modifyList(.d, lapply(.d, as_numeric))

parse_cells_stats <- function(.path) {
  .flines <- readLines(.path)
  .id_lines <- grep("^>", .flines)
  .growth_id <- which(.flines[.id_lines] == '>GROWTH')
  .out <- (.id_lines[.growth_id] + 1):(.id_lines[.growth_id+1] - 1) %>% # retrieve lines indices
    .flines[.] %>% paste(collapse='\n') %>%                             # concat lines
    textConnection %>% read.table(comment.char="", header=TRUE)         # read table
  
  return(list(growth=.out))
}

parse_frames_stats <- function(.path) {
  flines <- readLines(.path)

  # find column names for stats
  .fieldnames <- str_split(flines[1], '\t')[[1]] %>% sub("#", "", .) %>%
    gsub("-", "_", .) %>%
    gsub("\\(", "_", .) %>%
    gsub("\\)", "", .)
  .colnames <- str_split(flines[2], '\t')[[1]] %>% sub("#", "", .) %>%
    gsub("-", "_", .) %>%
    gsub("\\(", "_", .) %>%
    gsub("\\)", "", .)
  
  # parse all cells
  id_lines <- grep("^>CELL", flines)
  .out <- do.call(rbind, 
                  lapply(seq_along(id_lines)[-1], function(.i) {
                    .df <- (id_lines[.i-1] + 1):(id_lines[.i] - 1) %>% # retrieve lines indices
                      flines[.] %>% paste(collapse='\n') %>%           # concat lines
                      textConnection %>% read.table(comment.char="", header=FALSE) %>% # read table
                      setNames(.colnames)
                    .info <- flines[id_lines[.i-1]] %>% # extract line
                      strsplit('\t') %>% .[[1]] %>%     # split to vector
                      .[-1] %>% t %>% data.frame(stringsAsFactors=FALSE) %>%    # convert to df
                      setNames(.fieldnames[-1]) %>%
                      mutate(cell_ID=as.numeric(cell_ID), parent_ID=as.numeric(parent_ID), 
                             first_time_sec=as.numeric(first_time_sec), last_time_sec=as.numeric(last_time_sec))
                    cbind(.info, .df)
                  }) )
  return(.out)
}

compute_genealogy <- function(.id, .pid, .dgtype) {
#   browser()
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
  return(.cid)
}


