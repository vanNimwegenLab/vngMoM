# MoMA loading ####
utils::globalVariables(c(".", "where", # see https://github.com/r-lib/tidyselect/issues/201
                         "series", "medium", "duration", "interval", "paths", "condition", "step_idx"))

parse_yaml_conditions <- function(.path) {
  yaml::read_yaml(.path) %>% 
    # check `condition` values are unique and of length 1
    # note: this is easier to do before "unnesting" series
    (function(.l) {
      if (!identical(1, unique(purrr::map_dbl(.l, ~length(.$condition))))) stop("All values of `condition` must be scalar.") 
      if (length(unique(purrr::map_chr(.l, ~.$condition))) != length(.l)) stop("Values of `condition` must be unique.")
      return(.l) 
    }) %>% 
    # convert list to a tibble with list-columns and unnest `series` if needed
    purrr::map_dfr(function(.x) {
      .df <- purrr::map(.x, ~{ if (is.list(.)) . else list(.)}) %>% dplyr::as_tibble() # "nest" scalar values to list to allow bind_rows()
      if ('series' %in% names(.df)) .df <- dplyr::mutate(
        .df, series=purrr::map(series, function(.l) {purrr::map(.l, ~list(.)) %>% dplyr::as_tibble()}) ) %>% tidyr::unnest(series)
      return(.df)
    }) %>% 
    # check that `medium`, `duration` and `interval` have same length, and add step index
    dplyr::mutate(step_idx = purrr::pmap(list(medium, duration, interval), ~{ #browser()
      .ls <- unique(c(length(..1), length(..2), length(..3)))
      if (length(setdiff(.ls, 1)) > 1) stop("`medium`, `duration` and `interval` have same length (or be scalar values).")
      return(1:max(.ls))
    })) %>% 
    (function(.df) {
      # check that `paths` exist (warning) and are unique (error)
      .p <- dplyr::pull(.df, paths) %>% unlist()
      if (length(.p) > length(unique(.p))) stop("Values of `paths` must be unique.")
      purrr::walk(.p[!fs::dir_exists(.p)], ~warning("Directory does not exist: ", .))
      # check that other columns have scalar values
      setdiff(names(.df), c("condition", "medium", "duration", "interval", "step_idx", "paths")) %>% 
        purrr::walk(function(.x) if (max( dplyr::pull(.df, .x) %>% purrr::map_dbl(length) ) > 1) stop("Values of `", .x,"` must be scalar") )
      return(.df)
    }) %>% 
    tidyr::unnest(condition) %>% 
    tidyr::unnest(c(medium, duration, interval, step_idx)) %>% 
    tidyr::unnest(where(is.list) & !dplyr::all_of("paths"), keep_empty=TRUE) %>%  # unnest all remaining columns but `paths`
    # select(-paths) %>% 
    # select(condition, treatment, paths) %>% unnest(paths) %>% 
    identity()
}

parse_deepmoma_frames <- function(.path, .reader=readr::read_csv) {
  .reader(.path, skip=2, #lazy=FALSE,
           col_types = readr::cols( # specifying col types to avoid parsing NaN as chr and to speed up the import
             "lane_ID"=readr::col_character(),
             "genealogy"=readr::col_character(),
             "type_of_end"=readr::col_character(),
             .default=readr::col_double(), # if NaN raise problems, install last/dev version of vroom: https://github.com/tidyverse/readr/issues/1277#issuecomment-901297437
           )) %>% 
    tidyr::fill("type_of_end", .direction="up")
}

rename_deepmoma_vars <- function(.df) {
  .df %>% 
    dplyr::rename(
      "id"="cell_ID", #to remove spaces from column names
      "end_type"="type_of_end",
      "parent_id"="parent_ID",
      "vertical_bottom"="bbox_bottom px",
      "vertical_top"="bbox_top px",
      "center_x_px"="center_x px",
      "center_y_px"="center_y px",
      "width_px"="width px",
      "length_px"="length px",
      "tilt_radian"="tilt rad",
      "area_px2"="area px^2",
      "bg_area_px2"="bgmask_area px^2",
    ) %>% 
    dplyr::rename_with(
      function(.name) stringr::str_replace(.name, "(.*)_ch_(\\d+)", "\\1_ch\\2"),
      dplyr::matches("_ch_\\d+") ) %>% 
    dplyr::rename_with(
      function(.name) stringr::str_replace(.name, "fluo_cellmask_(\\d+)", "fluo_cellmask_ch\\1"),
      dplyr::matches("fluo_cellmask_\\d+") ) %>% 
    tidyr::extract("path", c("date"), ".*(20\\d{6})_.*", remove=FALSE, convert=FALSE) %>%
    tidyr::extract("lane_ID", c("pos", "gl"), ".*[Pp]os_(\\d+)_GL_(\\d+)", remove=TRUE, convert=TRUE)
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
