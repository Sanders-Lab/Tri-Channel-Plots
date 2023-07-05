################################  HEADER  ####################################
setwd('~/work/tri-channel-plots/')
source("R/plot_channels.R")

# Load or automatically install libraries
libs <- c("data.table", "ggplot2", "scales", "GenomicRanges", "readr",
          "tidyverse", "cowplot", "ggformula", "reshape2", "GenomicAlignments",
          "Rsamtools", "R.utils")
if (!require(libs, quietly = TRUE)) {
  BiocManager::install(libs)
}

# After the installation process completes, we load all packages.
lapply(libs, require, character.only = TRUE)


#################################  SETUP  ####################################
#Set working dir to tri-channel-plots folder
# setwd()
# Run only once if not done on files yet
# d.hap <- haplotagger(haplotag.bams.path = "/fast/groups/ag_sanders/scratch/kiwi_tmp/Martina/UC17_D/UC17_D/haplotag/bam/",
#                      output = "UC17_D_singleCell_haplotagData.txt")

# Load count and hap data (if previously genereated) for all cells if available
d <- fread("/fast/groups/ag_sanders/scratch/kiwi_tmp/Martina/UC17_D/UC17_D/counts/UC17_D.txt.gz") # scTRIP count output
d.hap <- fread("data/proc/P1530_singleCell_haplotagData.txt") # generate below


################################# FUNCTION ###################################


plot_channels(cell_ID = "P1530_i402",
              chromosome = "chr1",
              plot_range = c(4019024,9713340),
              channels = 4,
              count = F,
              roi = c(5900000,8700000,8800000))
