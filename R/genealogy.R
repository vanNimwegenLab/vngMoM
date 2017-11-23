# genealogy utilities ####
# requires the following cid format: 0BBBTBT (NB: Erik's has an extra colon e.g. 0:BBBTBT)
globalVariables(c("."))

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
  stringi::stri_sub(.cid, to=-2)
}

get_all_parents_cid <- function(.cid) {
  stringi::stri_sub(.cid, to=-(2:nchar(.cid))) 
}

get_daughters_cid <- function(.cid) {
  c(paste0(.cid, 'B'), paste0(.cid, 'T'))
}

has_progeny <- function(.cid, .all_ids) {
  .all_ids <- unique(.all_ids)
  is.element(paste0(.cid, 'B'), .all_ids) | is.element(paste0(.cid, 'T'), .all_ids)
}

is_parent_cid <- function(.cid, .pid) {
  # checks if .pid is a parent of .cid (vectorized over .cid)
  if (length(.pid) != 1) stop(".pid must have length 1.")
  .out <- logical(length(.cid))
  for (.i in seq_along(.cid)) {
    .out[.i] <- .pid %in% get_all_parents_cid(.cid[.i])
  }
  return(.out)
}

which_is_parent_cid <- function(.cid, .pid) {
  # find which value of .pid is a parent of .cid (vectorized over .pid and .cid)
  .out <- logical(length(.cid))
  for (.i in seq_along(.cid)) {
    .idx <- which(.pid %in% get_all_parents_cid(.cid[.i]))
    if (length(.idx)==1) {
      .out[.i] <- .idx
    } else if (length(.idx)==0) {
      .out[.i] <- 0
    } else {
      warning("several parent indices were found...")
      .out[.i] <- NA
    }
  }
  return(.out)
}

compute_daughters_numbers <- function(.cid)
  (paste0(.cid, 'B') %in% .cid) + (paste0(.cid, 'T') %in% .cid)

compute_daughters_numbers_raw <- function(.id, .pid) {
  .idx <- which(!duplicated(.id)) # indices of first occurences in .id
  sapply(.id, function(.i) sum(.pid[.idx]==.i))
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
  .ss <- sprintf('[%s%s]+', .bottom_char, .top_char) %>% stringr::str_match(.s, .) %>% strsplit('')
  .p <- lapply(.ss, function(.x) (.x[-length(.x)] == .x[-1]) %>% ifelse('O', 'N') %>% paste(collapse='')) %>% unlist
  
  .pos <- sprintf('[%s%s]+', .bottom_char, .top_char) %>% stringr::str_locate(.s, .) %>% .[,'start']
  paste0(stringr::str_sub(.s, 0, .pos), .p)
}

get_pole_age <- function(.s, .old_char='O', .new_char='N') {
  # get_pole_age() computes the number of consecutive identical letters at the end of the string
  # for old pole cells, the age is positive (number of consecutive division inheriting the old pole)
  # for new pole cells, the age is negative (number of consecutive division inheriting the youngest old pole)
  .bs <- paste0(.old_char, '+$') %>% stringr::str_match(.s, .) %>% stringr::str_length
  .ts <- paste0(.new_char, '+$') %>% stringr::str_match(.s, .) %>% stringr::str_length
  return( ifelse(is.na(.bs), -.ts, .bs) )
}

get_divs_from_mothercell <- function(.s, .old_char='B', .new_char='T') {
  # get_divs_from_mothercell() computes the number of divisions since the cell's ancestor was the mother cell of the GL
  .m_length <- sprintf('[^%s%s](%s+)', .old_char, .new_char, .old_char) %>% 
    stringr::str_match(.s, .) %>% .[,2] %>% stringr::str_length
  .s_length <- sprintf('[%s%s]+', .old_char, .new_char) %>% stringr::str_match(.s, .) %>% stringr::str_length
  return( .s_length - .m_length )
}

genealogy_relationship <- function(.s, .dist_max=3) {
  # sisters: div_min == div_max == 1
  # 1st cousins: div_min == div_max == 2
  # niece: div_min==1 && div_max==2
  
  # browser()
  if (length(.s)!=2) stop("the argument must be a vector of 2 strings.")
  .l1 <- stringr::str_length(.s[1])
  .l2 <- stringr::str_length(.s[2])
  if (abs(.l2 - .l1) > .dist_max) 
    return(NULL)
  
  .min_l <- min(.l1, .l2)
  .ss1 <- stringr::str_sub(.s[1], 1, .min_l)
  .ss2 <- stringr::str_sub(.s[2], 1, .min_l)
  if (.ss1 == .ss2) {
    .anc_length <- .min_l
  } else {
    .x1 <- .ss1 %>% stringr::str_split("") %>% unlist
    .x2 <- .ss2 %>% stringr::str_split("") %>% unlist
    .anc_length <- min(which(.x1 != .x2)) - 1 # 0 if no common ancestor
  }
  if (.min_l - .anc_length > .dist_max)
    return(NULL)
  .div_1 <- stringr::str_sub(.s[1], .anc_length+1, -1) %>% stringr::str_length # -1 to select the last char
  .div_2 <- stringr::str_sub(.s[2], .anc_length+1, -1) %>% stringr::str_length # -1 to select the last char
  .div_min <- min(.div_1, .div_2)
  .div_max <- max(.div_1, .div_2)
  if (.div_max > .dist_max) 
    return(NULL)
  return(list(rel=genealogy_ontology(.div_min, .div_max), ancestor=stringr::str_sub(.s[1], 1, .anc_length), 
              genealogy_1=.s[1], div_1=.div_1, genealogy_2=.s[2], div_2=.div_2, 
              div_min=.div_min, div_max=.div_max, less_div=which.min(c(.div_1, .div_2)) ))
}
