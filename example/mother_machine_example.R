# SET ENVIRONMENT ####
library(tidyverse)
library(vngMoM)

# # set a parallel environment to run multidplyr (ALL packages explicitely loaded before will be loaded too)
# library(multidplyr)
# mycluster <- min(30, parallel::detectCores()-1) %>%  # do not use more than 30 cores
#   new_cluster(.options = callr::r_session_options(user_profile = FALSE)) %>% # disable user profile to avoid endless startup with `renv`
#   cluster_library( # load currently loaded packages on each core (but multidplyr)
#     sessionInfo()$otherPkgs %>% names() %>% setdiff("multidplyr")
#   )

# Packages loaded from now on won't be attached to the multidplyr cluster
library(svglite)

dl <- 0.065 # µm/pixel (pixel size)

# LOAD MoMA DATA ####
# find raw data files
myfiles <- 
  parse_yaml_conditions(here::here("mother_machine_example.yml")) %>% 
  group_by(condition, antibio) %>% 
  slice(1L) %>% 
  unnest(paths) %>% 
  select(condition, antibio, data_path=paths) %>% 
  # for each path, find all files matched by the pattern .*\\d+\\.csv (e.g. *20151023.csv)
  mutate(path=map(data_path, ~vngMoM::find.files(., "ExportedCellStats_*.csv.zip"))) %>% 
  unnest(path) %>% 
  # TODO: check that file paths are unique
  identity()

# myconditions describes acquisition times and temporal change of each condition 
myconditions <-
  parse_yaml_conditions(here::here("mother_machine_example.yml")) %>% 
  group_by(condition, step_idx) %>%
  slice(1L) %>% 
  select(-antibio) %>% 
  mutate(m_start=cumsum(c(0, duration[-(length(duration))])) * 60,
         m_end=cumsum(duration) * 60 - 1e-5, 
         interval = interval*60, duration=NULL,
         m_cycle=vngMoM::value_occurence_index(medium), 
  ) %>% 
  mutate(times=pmap(list(m_start, m_end, interval), 
                   ~tibble(time_sec = seq(..1, ..2, ..3)) %>% 
                     mutate(frame=row_number()-1) )) %>% 
  unnest(times)

myframes <-
  myfiles %>% 
  mutate(data=map(path, parse_deepmoma_frames)) %>%
  # alternative to read non-standard csv files
  # mutate(data=map(path, ~parse_deepmoma_frames(
  #   ., .reader=(function(.p, ...) readr::read_delim(.p, delim=";", ...)) ))) %>% 
  unnest(data) %>% 
  mutate(regrow=NULL) %>%  # get rid of a column added manually in some input files
  rename_deepmoma_vars() %>% 
  # propagate time and medium info
  left_join(myconditions %>% ungroup() %>% unnest(paths) %>% 
              rename(data_path=paths) %>% select(-condition, -interval),
            by=c("data_path", "frame")) %>%
  select(-data_path) %>% 
  # append useful variables #(per cell)
  mutate(
    cell=paste(date, pos, gl, id, sep='.'),
    ugen=paste(date, pos, gl, genealogy, sep='.'),
    ugl=paste(date, pos, gl, sep='.'),
    vertical_center=(vertical_bottom + vertical_top)/2,
    length_um=length_px*dl,
    width_um=width_px*dl,
    volume_um=compute_pill_vol(length_um, width_um),
  ) %>%
  group_by(date, pos, gl, id) %>%
  mutate(
    time_start=first(time_sec), time_end=last(time_sec),
    avg_rank=round(mean(cell_rank)),
  ) %>% 
  identity()

# CONVERT FLUO UNITS ####
autofluo_per_um <- 422.8 # for MG1655 with 2 sec illumination at 17% intensity with ND4 filter
fp_per_dn <- 1 # dummy placeholder

myframes <- myframes %>%
  mutate(fluogfp_amplitude = fluo_ampl_ch1 - autofluo_per_um*length_um,
         gfp_nb = fluogfp_amplitude * fp_per_dn)



# RENDER ANALYSIS FILES ####
# calling `render()` or `render_site()` from the command line allows to execute the function 
# in the global envt (hence inheriting existing variables and keeping the newly created ones)...

rmarkdown::render('./mother_machine_example_analysis.Rmd')


# SAVE ENVIRONMENT ####
# save the current environment without the huge plots variables
ls(all.names = TRUE) %>% 
  setdiff(c("pls", "myplots")) %>% 
  save(list=., file=".RData", envir=.GlobalEnv)
