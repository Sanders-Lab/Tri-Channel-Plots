################################  HEADER  ####################################
setwd('~/work/tri-channel-plots/')
source("R/plot_channels.R")

# Load or automatically install libraries
libs <- c("data.table", "ggplot2", "scales", "GenomicRanges", "readr",
          "tidyverse", "cowplot", "ggformula", "reshape2", "GenomicAlignments",
          "Rsamtools", "utils")
if (!require(libs, quietly = TRUE)) {
  BiocManager::install(libs)
}

# After the installation process completes, we load all packages.
lapply(libs, require, character.only = TRUE)

#Set working dir to tri-channel-plots folder


##############################################################################

# Load count and hap data for all cells if available
d <- fread("/fast/groups/ag_sanders/scratch/kiwi_tmp/Martina/UC17_D/UC17_D/counts/UC17_D.txt.gz") # scTRIP count output
d.hap <- fread("data/proc/P1530_singleCell_haplotagData.txt") # generate below

# Run only once if not done on files yet
# d.hap <- haplotagger(haplotag.bams.path = "/fast/groups/ag_sanders/scratch/kiwi_tmp/Martina/UC17_D/UC17_D/haplotag/bam/",
#                      output = "UC17_D_singleCell_haplotagData.txt")

plot_channels(cell_ID = "P1530_i414",
              chromosome = "chr16",
              channels = 4,
              plot_range = c(4000000,9000000),
              roi = c(6019024,7713340))
