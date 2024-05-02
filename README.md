

This repository is the source code for the R package `vngMoM` which provides helper functions for analysing mother machine data obtaine with the image analysis software (deep)MoMA, as well as a template project for analysing such data (in the `example` subdirectory).
Most of this R code rely heavily on the `tidyverse` libraries. Getting familiar with those first, e.g. reading [R for data science](https://r4ds.hadley.nz), will definitely help you follow this analysis.


# Installation

It should be as simple as running the following: 

```
# install.packages('remotes') # only if devtools isn't already installed
remotes::install_github('nimwegenLab/vngMoM')
```

Note that to analyse data obtained from Legacy MoMA, one should use the `v0.1` branch of this package:

```
remotes::install_github('nimwegenLab/vngMoM', "v0.1")
```

# Usage

This repository contains an example project, which can be used as a template for new projects, and also illustrate useful features of vngMoM.

Detailed information about the usage of `vngMoM` is provided in the [readme file of `example`](example/README.md).


# Contribution
In case you're granted read-only access to https://github.com/nimwegenLab/vngMoM repository, the recommend use is the following:

- fork this repository in your own account (it will be kept private even if you use a free github plan; by the way, you can get 5 private repositories for free at https://education.github.com/discount_requests/new).
- clone your repository (e.g. https://github.com/myAccount/vngMoM) to your local computer:  
`git clone https://github.com/myAccount/vngMoM.git`
- add https://github.com/nimwegenLab/vngMoM as a new remote to your local repository, calling it `vng`:  
`git remote add vng https://github.com/nimwegenLab/vngMoM.git`
- you can now pull from the nimwegenLab/vngMoM repository and push to yours if you want to update it / share new code:  
`git pull vng master`  
`git pull vng dev` (if you want to pull the `dev` branch)  
`git push origin master` (or `git push` if you're lazy)
- after pushing your changes to your local repository (using an up-to-date version of nimwegenLab/vngMoM), you can create a pull request in the github web interface.
- (optionnal, and discouraged!) if you want to pull from `vng` and push to `origin` by default, you can change the git config locally:  
`git config --local branch.master.remote vng` (set `vng` as default remote)  
`git config --local branch.master.pushremote origin` (override the default push remote)

