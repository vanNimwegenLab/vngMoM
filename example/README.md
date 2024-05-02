
After image acquisition, Mother Machine data are processed using deepMoMA. Here we describe the analysis of deepMoMA's output using R, in particular:
- how to load MoMA data
- how to concatenate all data analysed in a project in the same dataframe
- simple examples of data selection and plotting.

Most of this R code rely heavily on the `tidyverse` libraries. Getting familiar with those first, e.g. reading [R for data science](https://r4ds.hadley.nz), will definitely help you follow this analysis.

`.` is the project directory. All paths to data files are defined locally.


# Loading data in R

## R environment

First, load the `mother_machine_example.Rproj` file of the template project (in `example`). 
The R environment used for this project is managed using `renv` (running `renv::init()` should restore all necessary package, but you need to make sure that you use an appropriate R version — at best the same as described in the file `renv.lock`). 
Learn more about collaborating with `renv` at https://rstudio.github.io/renv/articles/collaborating.html#collaborating-with-renv.

Run `src/mother_machine_example.R` to load the data and render the analysis files to html.
Note that calling `render()` or `render_site()` from the command line allows to execute the function in the global env() (hence inheriting existing variables and keeping newly created ones).


## Defining experimental conditions and paths to deepMoMA export files

Since experimental conditions are not recorded in deepMoMA output, they are documented, together with paths to deepMoMA export files, in a YAML file. Each condition is defined as an associative array which is expected to have the following properties:

- `condition` must be a unique string
- the combination of`medium`, `duration` and `interval` define "steps"; those variables can be lists (of the same length for a given condition) or scalars (value is then recycled)
- `paths` are lists of strings
- other custom fields can be added (must be scalar)
- for a given condition, several `series` can be declared as a list containing paths and other custom fields (see example below). 
NB: since series are *not* named, (at least one) custom field must be used to distinguish them.

For the rest, any valid YAML syntax should be accepted (mind tabs vs spaces for indentation).

This example defines two "conditions": the first one is specified using the series syntax (with a custom field `antibio`) and has a custom field `strain`; the second one only has a list of (one) replicate specified using `paths` and a custom field `antibio` (for consistency with the first condition).

```
- condition: AB_before
  duration: 360
  interval: 3
  medium: M9+glucose
  strain: MG1655_hi2
  series:
  - antibio: CIP
    paths:
    - ./data/20200922_curated/20200922_CIP_cell_stats
    - ./data/20201119_curated/20201119_CIP_cell_stats
  - antibio: CEF
    paths:
    - ./data/
    - ./data/20200922_curated/20200922_CEF_cell_stats
    - ./data/20201119_curated/20201119_CEF_cell_stats

- condition: AB_2h
  duration: [360, 6, 960, 120, 960]
  interval: 3
  medium: [M9, antibio, M9, antibio, M9]
  antibio: TMP
  paths:
  - ./data/20200922_curated/20200922_TMP_cell_stats
```

This is parsed to the following dataframe by `parse_yaml_conditions()`, with one row per condition and step:

```
# A tibble: 7 × 8
  condition duration interval medium     strain     antibio paths     step_idx
  <chr>        <int>    <int> <chr>      <chr>      <chr>   <list>       <int>
1 AB_before      360        3 M9+glucose MG1655_hi2 CIP     <chr [2]>        1
2 AB_before      360        3 M9+glucose MG1655_hi2 CEF     <chr [3]>        1
3 AB_2h          360        3 M9         NA         TMP     <chr [1]>        1
4 AB_2h            6        3 antibio    NA         TMP     <chr [1]>        2
5 AB_2h          960        3 M9         NA         TMP     <chr [1]>        3
6 AB_2h          120        3 antibio    NA         TMP     <chr [1]>        4
7 AB_2h          960        3 M9         NA         TMP     <chr [1]>        5
```

Note that it is possible to change the order of the fields in the YAML file, e.g. to document experiment in a more replicate-centered way.


## How is the information organised in dataframes

The data loading step produces three dataframes:

- myfiles has one row per data file, and keeps track of the corresponding condition (and series if relevant)
- myconditions has one row per condition and frame - it features all info related to time and environmental conditions, and is left-joined to myframes (by path and frame) in order to propagate those.
- myframes gathers all observations measured by deepMoMA with one row per cell and per frame.
NB: series don't need to be kept track of in myfiles and myconditions

`myframes` has the following variables (and more!) by default

- `id`, `gl`: the ID and growth lane of the given cell.
- `path`: the path of the CSV file containing the information of the given cell.
- `time\_sec`: time in seconds when the measurement has been done.
- `length\_um` length of the cell in μm at the time of the measurement.
- `end\_type`: end of a cell and it can be division (*div*), end of the measurement (*eod`), "pruned" or exit.
- `start\_time` and `end\_time`: time in seconds when the cell was born and died.
- `cell_rank`: cell rank counted from the channel closed end


`myframes` is designed to be enriched with new columns during the analysis; in particular dummy variables are helpful to filter the subset of the data that matches a given criteria. 

- `discard\_start`: when acquiring fluorescence, it is likely that cells experience some phototoxicity. In the case of the glucose/lactose switch, we identified that the doubling time of cells increase during the first 2 hours. Hence `discard\_start==TRUE` for the first 2 hours (aka "preexperiment step").

`filter(myframes, !discard_start)` should be used in any case to select valid observations.


## Estimating the GFP fluorescence

- `fluo\_amplitude` and `fluo\_background`: the fluorescence of the cell and the background as fitted from perl postprocessing scripts. The intensity of a cell is fitted with a Cauchy distribution mixed with a uniform distribution. The fluorescence amplitude is given by the Cauchy, while the background is given by the uniform.
- `fluo\_gfp\_amplitude`: the fluorescence of the GFP contained in the cell and it is defined as *fluo\_amplitude* minus the estimated autofluorescence (inferred from the length, using calibration data from a GFP-less strain).
- `gfp\_nb`: number of GFP molecules i.e. the GFP fluorescence rescaled with the parameter `fp\_per\_dn`, discussed later.


# Rudimentary data exploration

Once data are stored in the dataframe `myframes`, this tibble supports all the manipulations allowed by the dplyr library. Here are a few basic examples:

1. We can select cells satisfying specific requirements with the function `filter`. For example, if we want to isolate usable cells (not in the pre-experiment step, not excluded at the top of the channel), from the growth channel #14, and from a given condition (e.g. `test`), we can write
```
empty <- filter(myframes, !discard_start, condition=='test', gl==14)
```

2. If we want to stratify the cell by medium and date we can write
```
stratified <- group_by(myframes, medium, date)
```

3. We can compute statistics with the function *summarize*. If the data are stratified, the statistics is computed for
every group
```
stratified <- group_by(myframes, medium, date)
summarize(stratified, mu=mean(gfp_nb))
```

4. Finally we can combine different commands with the pipe operator %>%. Let's suppose we want to find the mean and variance of the number of gfp molecules discarding cells at the beginning of the experiment or that exited and stratified by medium and date. We can write
```
myframes %>% 
  filter(!discard_start) %>% 
  group_by(medium, date) %>% 
  summarize(mu=mean(gfp_nb), sd=sd(gfp_nb))
```


# Sharing the notebooks output online

## Preparing the rendered files

Once a project is completed, it is convenient to produce a static output of all its notebooks, which can later be browsed without initialising an R session. To facilitate this, the project template is designed as a Rmarkdown "site" (which requires rmarkdown ≥ 1.0).  
An example of this type of notebook rendering can be browsed at https://julou.github.io/MoM_lacInduction.

In _site.yml, shared properties are set, in particular the structure of the navigation bar.  
`exclude: ["*"]` is required to prevent all subdirectories to be copied (all the more so as symlinks are followed!)

The output is rendered in `/docs`, which is supported by Github Pages as website root; a `.nojekyll` empty file ensures that files are served as they are.  
NB: the `index.Rmd` file is required for site_render() to execute.

Here is an example of the minimal YAML header to put in each Rmarkdown file.
NB: date syntax from http://stackoverflow.com/questions/23449319

```
---
title: "My relevant title"
author: Thomas Julou
date: "`r format(Sys.time(), '%d %B, %Y')`"
---
```


## Deploying the rendered files using Github Pages

If your repo is synchronised with github at the url `https://github.com/myAccount/myProject", the corresponding Github Pages website wil be available at https://myAccount.github.io/myProject.

The corresponding files must be pushed to the github repository, and the deployment must be configured in the repository's settings (under "Code and automation > Pages").
There are primarily two ways to achieve this, and the settings must be adjusted accordingly.

### First option (not recommended): adding the rendered files to the project source code
Arguably the simplest option is to leave the rendered files in the `docs` subdirectory where they have been output by Rmarkdown. The obvious drawback is that the code is now cluttered by this additional subdirectory.

To do the versioning this way, simply delete the `.gitignore` file in `docs`, then commit and push the files as usual.

### Second option: moving the rendered files to another branch
In order to separate the code and the rendered files, it is possible to copy the content of the `docs` directory (of the code branch `master`) to the root of a dedicated branch (typically named `gh-branch`).

```
# create branch if it does not exist
git checkout -b gh-pages 

git checkout master # make sure you are on the "code" branch 
git checkout gh-pages "docs/*"
```

After that, don't forget to commit and push the content of both branches.

