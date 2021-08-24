globalVariables(c("covar"))

.onLoad <- function(libname, pkgname) {
  options("readr.num_columns" = 0) # disable printing when loading data with readr
  
  if (!exists("scale_colour_periodic_brewer", mode="function"))
    scale_colour_periodic_brewer <- function(..., .n=4) 
      ggplot2::scale_colour_manual(..., values = rep(c(RColorBrewer::brewer.pal(9, 'Set1')[1:.n], 'gray42'), 1e4), 
                                   na.value='gray25')

  if (!exists("data2preproc_dir", mode="function"))
    data2preproc_dir <- function(.f) dirname(.f)
  if (!exists("data2preproc_file", mode="function"))
    data2preproc_file <- function(.f)
      basename(.f) %>% tools::file_path_sans_ext %>% paste0("_preproc.txt")
  if (!exists("data2preproc", mode="function"))
    data2preproc <- function(.path)
      file.path(data2preproc_dir(.path), data2preproc_file(.path))
  
  invisible()
}

# general convenience functions ####
# file_path_sans_ext <- tools::file_path_sans_ext

find.files <- function(.path, .name="", .pattern="", .mindepth=NULL, .maxdepth=NULL, .follow_symlinks=FALSE) {
# find.files is a minimalist wrapper around bash's `find` to be used as a faster alternative to list.files
# NB: setting .follow_symlinks=TRUE can hugely slow down the function
  .fargs <- .path
  if (.follow_symlinks) .fargs <- c("-L", .fargs)
  if (.name!="") .fargs <- c(.fargs, paste0("-name \'", .name, "\'"))  # note: escape the pattern string with ' to avoid shell globbing
  if (.pattern!="") .fargs <- c(.fargs, "-regextype posix-basic", paste("-regex", .pattern))
  if (!is.null(.mindepth)) .fargs <- c(.fargs, paste("-mindepth", .mindepth))
  if (!is.null(.maxdepth)) .fargs <- c(.fargs, paste("-maxdepth", .maxdepth))
  system2("find", .fargs, stdout=TRUE)
}

between_or <- function(.x, .a, .b) {
  if (length(.a) != length(.b)) stop('.a and .b must have the same length.')
  lapply(seq(length(.a)), function(.i) dplyr::between(.x, .a[.i], .b[.i])) %>%
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
  .out <-  lapply(seq(length(.mins)), function(.i) dplyr::between(.xs, .mins[.i], .maxs[.i])) %>% 
    do.call(cbind, .) %>%
    apply(1, which)
  .ls <- sapply(.out, length) 
  if (any(.ls > 1)) warning('some .xs values are matched by several intervals.')
  .out[.ls!=1] <- NA
  return(unlist(.out))
}

as_numeric <- function(.x) as.numeric(as.character(.x))
as_numeric_df <- function(.d) utils::modifyList(.d, lapply(.d, as_numeric))

mycut <- function(.x, ...) {
  .k <- Hmisc::cut2(.x, ...)
  .mids <- sapply(stringr::str_match_all(as.character(levels(.k)), "[\\[\\(](.*),(.*)[\\]\\)]"), 
                 function(.m) {
                   mean(c(as.numeric(.m[2]), as.numeric(.m[3])))
                 })
  .mids[.k]
}

gm_mean <- function(x, na.rm=TRUE) {
# computes geometric mean (from http://stackoverflow.com/a/25555105/576684)
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

cor_sym <- function(x, y, ...) {
# compute correlation after duplicating the values to have a symmetric scatter relative to the first diagonale
  stats::cor(c(x, y), c(y, x), ...)
}

cov_pop <- function(x,y=NULL) {
  browser()
  stats::cov(x,y)*(nrow(x)-1)/nrow(x)
}

compute_combn_covar <- function(.df, .groups="str_cond", .vars=
                                  list(l_birth="l_birth", 
                                       l_increase=c("l_div", "dl"),
                                       div_time=c("div_time", "alpha"))) {
  # vars is a named list of named vectors of variable names (provided as strings)
  # for vectors with a single value, the variable is  explicitely 
  stopifnot(!is.null(names(.vars)))
  # check that all variables with more than one value have names
  stopifnot(if (any((names(.vars)==""))) !any( (names(.vars)=="") & (sapply(.vars, length)>1) ) else TRUE) # forbid F/T
  # browser()
  .out <- .df %>% dplyr::group_by_(.groups) %>% 
    dplyr::select(!!unlist(.vars)) 
  
  .names <- names(.vars)
  for (.i in which(.names!=''))
    .out <- tidyr::gather_(.out, paste0(.names[.i], '_var'), .names[.i], .vars[[.i]]) %>% 
    dplyr::group_by_(paste0(.names[.i], '_var'), add=TRUE)
  
  .out %>% 
    dplyr::do(covar=cov_pop(dplyr::select_(., .dots=paste0(.names[which(.names!='')], '_var') %>% c(.groups) %>% paste0('-', .)) ),
       scaled_covar=cov_pop(scale(dplyr::select_(., .dots=paste0(.names[which(.names!='')], '_var') %>% c(.groups) %>% paste0('-', .)) )) ) %>% 
    dplyr::mutate(r2=1-det(covar)/prod(diag(covar)))
}

median_pi <- function(d){
# compute the median and its 95% posterior interval
  d <- sort(d)
  n <- length(d)
  if (n<=1) return(data.frame(y=NA, ymin=NA, ymax=NA))
  
  nmedian <- ceiling(n/2)
  nmin <- round(n/2 - sqrt(n*2))
  nmax <- round(n/2 + sqrt(n*2))
  
  if (nmin<1) nmin <- 1
  if (nmax>n) nmax <- n
  
  data.frame(y=d[nmedian], ymin=d[nmin], ymax=d[nmax])
}
