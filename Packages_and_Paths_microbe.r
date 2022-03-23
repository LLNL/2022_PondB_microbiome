#  Packages and Paths

#### Load Packages ####
packages <- c("dada2", "ggplot2", "phyloseq", "vegan", "tidyr", "viridis", "reshape2", "stringr", "ggthemes", "pander", "plyr", "ranacapa", "ade4", "FactoMineR", "factoextra", "ggrepel", "dbstats", "Rcpp", "ape", "dplyr", "forcats", "colorspace", "ggsci", "microbiome", "data.table", "metagMisc", "ggpubr", "decontam", "Hmisc", "RColorBrewer", "cowplot", "ALDEx2", "CoDaSeq", "zCompositions", "qiime2R", "ampvis", "biomformat", "picante", "hablar", "sigminer", "microbiomeMarker")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  BiocManager::install(packages[!installed_packages])
}

# Install packages not yet installed
remotes::install_github("ShixiangWang/sigminer")
remotes::install_github("yiluheihei/microbiomeMarker")
BiocManager::install("philr")

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

#### Set Paths ####
# Change the sequence batch number
batch <- c("SeqBatch4")
#SeqBatch3
#SeqBatch2
#SeqBatch1
#Sediment

path <- <path to main folder for project>
quality <- file.path(paste0(path, "/01.quality/", batch))
filter <- file.path(paste0(path, "/02.filter_trim/", batch))
error <- file.path(paste0(path, "/03.error/", batch))
dereplication <- file.path(paste0(path, "/04.dereplication/", batch))
merging <- file.path(paste0(path, "/05.merging/", batch))
qiime_otu <- file.path(paste0(path, "/06.qiime_otu97/input"))
path_phy <- file.path(paste0(path, "/07.phyloseq/"))
input <- file.path(paste0(path, "/input/", batch))

setwd(path)

dir.create(quality, showWarnings = FALSE, recursive = TRUE)
dir.create(filter, showWarnings = FALSE, recursive = TRUE)
dir.create(error, showWarnings = FALSE, recursive = TRUE)
dir.create(dereplication, showWarnings = FALSE, recursive = TRUE)
dir.create(merging, showWarnings = FALSE, recursive = TRUE)
dir.create(qiime_otu, showWarnings = FALSE, recursive = TRUE)
dir.create(path_phy, showWarnings = FALSE)

# set save file
save.image(paste0(path,"/PondB_Microbiome_2022.RData"))

# load file
load(paste0(path,"/PondB_Microbiome_2022.RData"))
