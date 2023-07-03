################################  HEADER  ####################################
setwd('~/work/tri-channel-plots/')
source("R/plot_channels.R")

# Load or automatically install libraries
libs <- c("data.table", "ggplot2", "scales", "GenomicRanges", "readr",
          "tidyverse", "cowplot", "ggformula", "reshape2", "GenomicAlignments",
          "Rsamtools")
if (!require(libs, quietly = TRUE)) {
  BiocManager::install(libs)
}

# After the installation process completes, we load all packages.
lapply(libs, require, character.only = TRUE)

#Set working dir to tri-channel-plots folder


##############################################################################

# Load count and hap data for all cells if available
d <- fread("data/raw/AGLCD.txt.gz") # scTRIP count output
d.hap <- fread("data/proc/P1530_singleCell_haplotagData.txt") # generate below

# Run only once if not done on files yet
# d.hap <- haplotaggeR(haplotag.bams.path = "haplotag/bam/",
#                      chromosomes = paste0("chr", c(1:22, "X", "Y")),
#                      output = "P1530_singleCell_haplotagData.txt") 

plot_channels(cell_ID = "i567",
              chromosome = "chr2",
              channels = 4,
              plot_range = NULL,
              roi = NULL)
