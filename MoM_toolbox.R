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


load_timm_data <- function(.path, .perl_cmd='perl',
                           .frames_script_path="/import/bc2/home/nimwegen/GROUP/MM_Analysis/scripts/get_height_size_intens.pl", 
                           .cells_script_path="/import/bc2/home/nimwegen/GROUP/MM_Analysis/scripts/estimate_division_growth2.pl", 
                           .force=FALSE, .verbose=FALSE,
                           .data2preproc=identity) {
  if (!file.exists(.path)) stop("input path is not valid.")
  .out_frames_path <- .path %>% .data2preproc %>% file_path_sans_ext %>% paste0("_frames.txt")
  .out_cells_path <- .path %>% .data2preproc %>% file_path_sans_ext %>% paste0("_cells.txt")
  dir.create(dirname(.out_frames_path), recursive=TRUE, showWarnings=FALSE)
  
  # run perl scripts if output doesn't exist
  if (!file.exists(.out_frames_path) || .force==TRUE) {
    if (.verbose) print(paste("Converting", .path, "to", .out_frames_path))
#     browser()
    system2(.perl_cmd, args=c(.frames_script_path, .path), stdout=.out_frames_path)
  }
  if (!file.exists(.out_frames_path)) stop("Frame file cannot be found.")

  if (!file.exists(.out_cells_path) || .force==TRUE) {
    if (.verbose) print(paste("Converting", .out_frames_path, "to", .out_cells_path))
    system2(.perl_cmd, args=c(.cells_script_path, .out_frames_path), stdout=.out_cells_path)
  }
  
  .cells_out <- read.table(.out_cells_path, sep=" ", comment.char="", header=TRUE)
  names(.cells_out) <- gsub("X\\.", "", names(.cells_out))
  .cells_out <- mutate(.cells_out, cid = compute_genealogy(id, parent_id, daughter_type))
  .frames_out <- parse_frames_stats(.out_frames_path) %>%
    inner_join(.cells_out, ., by="id")
  if (dim(.cells_out)[1] > 0) {
    .m <- str_match(basename(.path), ".*(\\d{8})_.*pos(\\d+).*_GL(\\d+).*")
    return(list(cells=data.frame(date=as.numeric(.m[1, 2]), pos=as.numeric(.m[1, 3]), gl=as.numeric(.m[1, 4]), .cells_out),
                frames=data.frame(date=as.numeric(.m[1, 2]), pos=as.numeric(.m[1, 3]), gl=as.numeric(.m[1, 4]), .frames_out) ))
  }
}

as_numeric <- function(.x) as.numeric(as.character(.x))
as_numeric_df <- function(.d) modifyList(.d, lapply(.d, as_numeric))

parse_frames_stats <- function(.path) {
  flines <- readLines(.path)
  # browser()

  # find column names for stats
  .m <- str_match(flines[1], ".* STATS: (.*)$")
  .colnames <- str_split(.m[1,2], '\t')[[1]]
  
  # parse all cells
  id_lines <- grep("^>cell_", flines)
  .out <- do.call(rbind, 
                  lapply(seq_along(id_lines)[-1], function(.i) {
                    .df <- do.call(rbind, 
                                   strsplit(flines[(id_lines[.i-1]+1):(id_lines[.i]-1)], '\t'))
                    colnames(.df) <- .colnames
                    .info_str <- flines[id_lines[.i-1]]
                    .info_spl <- strsplit(.info_str, '_')[[1]]
                    # data.frame(id=.info_spl[2], parent_id=.info_spl[4], birthframe=.info_spl[6], endframe=.info_spl[8], type_of_end=.info_spl[10], .df)
                    data.frame(id=.info_spl[2], .df)
                  }) )
  return(as_numeric_df(.out))
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


