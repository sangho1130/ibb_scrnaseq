### R version 4.2.2 or version > 4.x.x recommended 
### execute following commands line by line

# Seurat - scRNA-seq analysis pipeline in R
# https://satijalab.org/seurat/articles/install.html
install.packages("Seurat")

# Data integration tool
# https://github.com/immunogenomics/harmony
install.packages("harmony")

# Cran packages required for workshop
# https://cran.r-project.org/web/packages/plyr/index.html
install.packages("plyr")
# https://cran.r-project.org/web/packages/dplyr/
install.packages("dplyr")
# https://cran.r-project.org/web/packages/ggplot2/
install.packages("ggplot2")

# check installed packages
# if packages are not properly loaded, you will see error messages
# re-install and see stdout messages during installation
library(Seurat)
library(plyr)
library(dplyr)
library(ggplot2)
library(harmony)
