utils::globalVariables(c(".", "str", "line", "id", "pid", "frame", "top", "bottom"))
                         
parse_moma_tracks <- function(.path) {
  if (!file.exists(.path)) stop(sprintf("input path is not valid (%s)", .path))
  flines <- readLines(.path, warn=FALSE)
  # browser()
  
  offset <- stringr::str_which(flines, "^trackRegionInterval") %>% 
    flines[[.]] %>% 
    stringr::str_match("^trackRegionInterval = \\[(\\d+),\\d+]") %>% 
    .[2] %>% as.numeric()
 
  dplyr::bind_rows(
    stringr::str_which(flines, "^id=") %>% 
      dplyr::tibble(line=., str=flines[line]) %>% 
      tidyr::extract(str, c("id", "pid"), "^\\h*id=(\\d+); pid=(-?\\d+)", convert=TRUE),
    stringr::str_which(flines, "^\\h*frame=") %>% 
      dplyr::tibble(line=., str=flines[line]) %>% 
      tidyr::extract(str, c("frame", "top", "bottom"), "^\\h*frame=(\\d+);.*pixel_limits=\\[(\\d+),(\\d+)\\]", convert=TRUE),
    stringr::str_which(flines, "^\\h*(DIVISION|EXIT|ENDOFDATA)$") %>% 
      dplyr::tibble(line=., str=flines[line]) %>% 
      tidyr::extract(str, c("end_type"), "^\\h*(DIVISION|EXIT|ENDOFDATA)$"),
  ) %>% 
    dplyr::arrange(line) %>% 
    tidyr::fill(id, pid) %>% 
    tidyr::fill(end_type, .direction="up") %>% 
    dplyr::filter(!is.na(frame)) %>%
    dplyr::mutate(top=top+offset, bottom=bottom+offset, center=top+(bottom-top)/2, line=NULL, 
                  end_type=factor(end_type) %>% forcats::fct_recode('div'='DIVISION', 'exit'='EXIT', 'eod'='ENDOFDATA'))
}
