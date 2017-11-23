---
title: "vngMoM"
---


After image acquisition, Mother Machine data are processed using MoMA. Here we describe the analysis of MoMA's output using R, in particular:
- how to load MoMA data and apply its postprocessing scripts in perl
- how to concatenate all data analysed in a project in the same dataframe
- simple examples of data selection and plotting

Most of this R code rely heavily on Hadley Wickham's libraries, including stringr, tidyr, dplyr, multidplyr and ggplot2. Getting started with dplyr and ggplot2 first will definitely help you follow this analysis.

`.` is the project directory. All paths to data files (raw or preprocessed) must be defined locally.


# Getting started
## Installation

It should be as simple as running the following:

```
# install.packages('devtools') # only if devtools isn't already installed
devtols::install_github('')
```

## Contribution
In case you're granted read-only access to https://github.com/vanNimwegenLab/vngWetLabR repository, the recommend use is the following:

- fork this repository in your own account (it will be kept private even if you use a free github plan; by the way, you can get 5 private repositories for free at https://education.github.com/discount_requests/new).
- clone your repository (e.g. https://github.com/myAccount/vngWetLabR) to your local computer:  
`git clone git@github.com:myAccount/vngWetLabR.git`
- add https://github.com/vanNimwegenLab/vngWetLabR as a new remote to your local repository, calling it `vng`:  
`git remote add vng git@github.com:vanNimwegenLab/vngWetLabR.git`
- you can now pull from the vanNimwegenLab/vngWetLabR repository and push to yours if you want to update it / share new code:  
`git pull vng master`  
`git pull vng dev` (if you want to pull the `dev` branch)  
`git push origin master` (or `git push` if you're lazy)
- after pushing your changes to your local repository (using an up-to-date version of vanNimwegenLab/vngWetLabR), you can create a pull request in the github web interface.
- (optionnal, and discouraged!) if you want to pull from `vng` and push to `origin` by default, you can change the git config locally:  
`git config --local branch.master.remote vng` (set `vng` as default remote)  
`git config --local branch.master.pushremote origin` (override the default push remote)




# Converting MoMA output using postprocessing perl scripts

In a first step, postprocessing perl scripts (producing more precise length and fluorescence estimates) are applied to the following raw output files from MoMA: 

```
./data/date1/ (cond1/cond2/) pos0/pos0GL01.csv
                                  pos0GL02.csv
                                  ...
                                  pos1GL01.csv
                                  ...
```

the resulting "preprocessed" data files are written to the directory of your choice (e.g. `preproc`; the "preprocessing" nature of this step is seen from this R analysis perspective, not from the MoMA perspective where it corresponds to "postprocessing"):

```
./preproc/date1/ (cond1/cond2/) pos0/pos0GL01_frames.txt
                                     pos0GL02_frames.txt
                                     ...
                                     pos1GL01_frames.txt
                                     ...
```

> NB: The directory structure of the raw data directory can be replicated in the preproc directory. This is controlled by the function `data2preproc`: it converts raw data path to preproc path. If you want cache files to be stored along with raw data, simply use `identity`. For the sake of flexibility and readability, `data2preproc` is split into two subfunctions `data2preproc_dir` and `data2preproc_file` which handle the name conversion of the directory and file respectively.

`process_moma_data()` handles this preprocessing step (using the computing power of the BC2 cluster) and returns a vector containing the paths to all preprocessed files. For any raw data file whose corresponding preprocessed file is not found, the preprocessing is started using the BC2 queue and the function will return an error message; the user must run this function once more after all jobs are executed.
Alternatively, if all preprocessed files are found, the function returns a vector containing the paths to all preprocessed files. This approach uses the preprocessed files as cache when an analysis is rerun later by the user; the analysis can be rerun from scratch by setting `force=TRUE`.

`parse_frames_stats()` parses one preprocessed file to an R dataframe. It is used over the output of `process_moma_data()` to load them all and concatenate the dataframes. The result is a large dataframe `myframes` gathering all data considered in the analysis (with one row per cell per frame and one column per variable).

`myframes` is designed to be enriched with new columns during the analysis; in particular dummy variables are helpful to distinguish the setup of the data that matches a given driteria. `myframes` has the following variables (and more!) by default

- `id`, `gl`: the ID and growth lane of the given cell.
- `path`: the path of the CSV file containing the information of the given cell.
- `time\_sec`: time in seconds when the measurement has been done.
- `length\_um` length of the cell in μm at the time of the measurement.
- `end\_type`: end of a cell and it can be division (*div*), end of the measurement (*eod`) or exit.
- `start\_time` and `end\_time`: time in seconds when the cell was born and died.
- `b\_rank`: cell rank counted from the channel closed end ("b" stands for "bottom")


# Computing additional biologically relevant variables

## Appending experimental conditions

By default, the naming of MoMA files don't encode details of the conditions used in each experiment. For this purpose, we define a dataframe (`myconditions`) that describes the details of a given condition (that might have been used over several replicate experiments) and where corresponding data are located (so as to map each experiment to a given condition).

This information is appended to `myframes` and parsed to identify the following variables for each row: in which medium the cells currently are (`medium`), how long has elapsed since the last medium switch (`m_start`), and how many time they have seen this medium already (`m_cycle`).

With this syntax, all data (several conditions and controls) can be loaded in the same step and then used at will for downstream analysis.


## Appending dummy variables

- `discard\_start`: when acquiring fluorescence, it is likely that cells experience some phototoxicity. In the case of the glucose/lactose switch, we identified that the doubling time of cells increase during the first 2 hours. Hence `discard\_start==TRUE` for the first 2 hours (aka "preexperiment step").

- `discard\_top`: because phase contrast ilumination is very uneven next to the growth channel exit, it is important to exclude cells once they have reached a minimum distnace to the exit. All later frames of the same cell as well as its progeny is excluded. \
NB: during the preexperiment, measurements are not used and bottom cells are likely to touch the exit threshold because they progressively enter the channel. Hence we don't label any cell as `discard\_top==TRUE` since it would be likely to exclude all the cells in the growth channel (as being the progeny of the bottom cell which touched the threshold once).

`filter(myframes, !discard_start, !discard_top)` should be used in any case to select valid observations.


## Estimating the GFP fluorescence

- `fluo\_amplitude` and `fluo\_background`: the fluorescence of the cell and the background as fitted from perl postprocessing scripts. The intensity of a cell is fitted with a Cauchy distribution mixed with a uniform distribution. The fluorescence amplitude is given by the Cauchy, while the background is given by the uniform.
- `fluo\_gfp\_amplitude`: the fluorescence of the GFP contained in the cell and it is defined as *fluo\_amplitude* minus the estimated autofluorescence (inferred from the length, using calibration data from a GFP-less strain).
- `gfp\_nb`: number of GFP molecules i.e. the GFP fluorescence rescaled with the parameter `fp\_per\_dn`, discussed later.


# Data exploration

Once data are stored in the dataframe `myframes`, this tibble supports all the manipulations allowed by the dplyr library. Here are a few basic examples:

1. We can select cells satisfying specific requirements with the function `filter`. For example, if we want to isolate usable cells (not in the pre-experiment step, not excluded at the top of the channel), from the growth channel #14, and from a given condition (e.g. `test`), we can write
```
empty <- filter(myframes, !discard_top, !discard_start, condition=='test', gl==14)
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
  filter(!discard_top, !discard_start) %>% 
  group_by(medium, date) %>% 
  summarize(mu=mean(gfp_nb), sd=sd(gfp_nb))
```


