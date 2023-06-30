################################  HEADER  ####################################
# Source functions
source("./haplotagger_dev.R")

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

library(tidyverse)
library(data.table)


##############################################################################


#################################  DATA  #####################################
#Generate and load data frames from scTRIP output files

# Load count data for all cells
d <- as.data.frame(read_table(file = '/fast/groups/ag_sanders/scratch/kiwi_tmp/AGLCD_hg38_hgsvc/AGLCD/counts/AGLCD.txt.gz',
                              col_names = T))
write.table(d,
            '../../data/raw/AGLCD.txt.gz')
d <- fread('../../data/raw/AGLCD.txt.gz')
d_dt <- as.data.table(d)
d_dt


d_dt[c != 0]


d_dt[, .SD, by = .(start, 
                   end, 
                   cell,
                   chrom)]




d.hap <- haplotaggeR(haplotag.bams.path = "/fast/groups/ag_sanders/scratch/kiwi_tmp/AGLCD_hg38_hgsvc/AGLCD/haplotag/bam/",
                     #MUST be the bam files in haplotag folder !!!
                     chromosomes = paste0("chr", c(1:22, "X")),
                     #here I run all chromosomes, default is all chr so you can skip this if wanted
                     output = "./P1530_singleCell_haplotagData.txt") 
write.table(d.hap,
            "../../data/proc/P1530_singleCell_haplotagData.txt"
)
#provide path to filename for output of the dataframe so you can load it next time, default is NULL

##############################################################################
d.hap <- fread("../../data/proc/P1530_singleCell_haplotagData.txt")

d.hap %>% head()
d  %>% head()
in.d$cell  %>% unique()
in.d$



###############################  FUNCTION  ###################################
#For now only works for whole chromosomes, so dont use plot_range please!
tri_plot <- plot_trichannel(cell_ID = "i563", #keep string format!!!!
                  roi = NULL, #highlight region with dotted lines, if not supplied vanishes
                  plot_range = c(74000000:120000000), # DO NOT USE FOR NOW!!
                  chromosome =  "chr7") #choose chromosome to plot
##############################################################################

pdf("../../plots/tmp_plot.pdf")
plot(tri_plot)
dev.off()


###############################   SAVING   ###################################
#If you wanna save your tri channel plot you could store it in an obj and use ggsave:
#
# plot <- plot_3channel_chr(cell_ID = "P1530_i505_",
#                           roi = c(5239802:7711458),
#                           chromosome =  "chr7")
#
#ggsave(filename = "my_plot.pdf", plot = obj, dpi = 300) #carefull saves in working dir
##############################################################################
