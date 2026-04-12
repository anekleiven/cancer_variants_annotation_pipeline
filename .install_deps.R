# Script to install all R dependencies for the pipeline
if (!require("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if (!require("remotes", quietly = TRUE)) install.packages("remotes")
if (!require("data.table", quietly = TRUE)) install.packages("data.table")

# Installing the specific bioinformatics package 
if (!require("geneOncoX", quietly = TRUE)) remotes::install_github("sigven/geneOncoX")

