################################  HEADER  ####################################
# Source functions
source(file = "/data/gpfs-1/groups/ag_sanders/work/projects/kiwi/scripts/dev/haplotaggeR_dev.R")

# Check if BiocManager is installed
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Load or automatically install libraries
libs <- c("data.table", "ggplot2", "scales", "GenomicRanges", "readr",
          "tidyverse", "cowplot", "ggformula", "reshape2", "GenomicAlignments",
          "Rsamtools")
if (!require(libs, quietly = TRUE)) {
  BiocManager::install(libs)
}

# After the installation process completes, we load all packages.
lapply(libs, require, character.only = TRUE)

#Set working dir to scTRIP output folder
setwd('/fast/groups/ag_sanders/scratch/bendy_tmp/20220705_scTRIP_Martina_P1530/') # <-- hierarchy example

##############################################################################


#################################  DATA  #####################################
#Generate and load data frames from scTRIP output files

# Load count data for all cells
d <- as.data.frame(read_table(file = '/fast/groups/ag_sanders/scratch/bendy_tmp/20220705_scTRIP_Martina_P1530/counts/P1530/100000_fixed.txt.gz', col_names = T))

d.hap <- haplotaggeR(haplotag.bams.path = "/fast/groups/ag_sanders/scratch/bendy_tmp/20220705_scTRIP_Martina_P1530/haplotag/bam/P1530/100000_fixed_norm.selected_j0.1_s0.5_scedist20",
                     #MUST be the bam files in haplotag folder !!!
                     chromosomes = paste0("chr", c(1:22, "X")),
                     #here I run all chromosomes, default is all chr so you can skip this if wanted
                     output = "/fast/groups/ag_sanders/scratch/bendy_tmp/20220705_scTRIP_Martina_P1530/haplotag/P1530_singleCell_haplotagData.txt") 
#provide path to filename for output of the dataframe so you can load it next time, default is NULL

##############################################################################


###############################  FUNCTION  ###################################
#For now only works for whole chromosomes, so dont use plot_range please!
plot_channels(cell_ID = "i563", #keep string format!!!!
                  roi = c(74000000:120000000), #highlight region with dotted lines, if not supplied vanishes
                  plot_range = , # DO NOT USE FOR NOW!!
                  chromosome =  "chr7",
                  channels = 3) #choose chromosome to plot
##############################################################################


###############################   SAVING   ###################################
#If you wanna save your tri channel plot you could store it in an obj and use ggsave:
#
# plot <- plot_3channel_chr(cell_ID = "P1530_i505_",
#                           roi = c(5239802:7711458),
#                           chromosome =  "chr7")
#
#ggsave(filename = "my_plot.pdf", plot = obj, dpi = 300) #carefull saves in working dir
##############################################################################