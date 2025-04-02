#Instal and load the required packages
install.packages(Signac)
library(Signac)
install.packages(Seurat)
library(Seurat)
install.packages(GenomeInfoDb)
library(GenomeInfoDb)
install.packages(EnsDb.Hsapiens.v75)
library(EnsDb.Hsapiens.v75)
install.packages(ggplot2)
library(ggplot2)
install.packages(patchwork)
library(patchwork)
install.packages(Rsamtools)
library(Rsamtools)

#Setting seed to make results reproducible
set.seed(1234)