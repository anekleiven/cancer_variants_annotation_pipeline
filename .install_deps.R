# Script to install all R dependencies for the pipeline
if (!require("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if (!require("remotes", quietly = TRUE)) install.packages("remotes")

# Installing the specific bioinformatics package 
if (!require("geneOncoX", quietly = TRUE)) remotes::install_github("sigven/geneOncoX")

