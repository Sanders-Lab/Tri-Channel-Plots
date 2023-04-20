library(data.table)
library(ggplot2)
library(scales)
library(GenomicRanges)
library(readr)
library(tidyverse)
library(cowplot)
setwd('~/Documents/Data/Sequencing/AGLCD/P1530_rerun/') #scTRIP output folder

## assign cell ID and region of interest/genomic range:
cell_ID="P1530_i563_" # cell of interest
genomic_range = c(140000000:145000000) #plot ylim, like a padding around roi
roi = c(142299011:142813287) #highlighted dotted line area in genomic_range
chromosome = "chr7" #specify chromosome for whole chr filtering here

## load in template strand data for all cells from count folder in working dir
File <- ('~/Documents/Data/Sequencing/AGLCD/P1530_rerun/counts/P1530/100000_fixed.txt.gz')
inputFile <- paste0("zcat <", File)
d = fread(cmd = inputFile) # reads in the file
d <- merge(d, d[,.(total = sum(w+c)), by = cell], by = "cell")
d <- merge(d, d[,.(mean = median(w+c)), by = cell], by = "cell")
#d <- d[total > 1e5,] # filters based on minimum read number

#filter out roi data here depending on what youre analysing (remove/add arguments as needed)
ind <- d %>% 
  filter(cell == cell_ID & start %in% genomic_range & chrom == chromosome)
in.d <- melt(ind, measure.vars = c("c", "w"))

#generate hap_data from haplotagData_plottingTable.R here
haplotag.bams.path <- '~/Documents/Data/Sequencing/AGLCD/P1530_rerun/haplotag/bam/P1530/100000_fixed_norm.selected_j0.1_s0.5_scedist20/'
haplotag.bams <- list.files(path = haplotag.bams.path, pattern = "\\.bam$", full.names = T)
region<- paste0("chr", c(1:22, "X"))

all.counts <- data.frame()

for (i in 1:length(haplotag.bams)) {
  bam <- haplotag.bams[i]
  message("Processing bamfile ", bam, " ...")
  filename <- basename(bam)
  cell.id <- unlist(strsplit(filename, "\\."))[1]
  fragments <- bamregion2GRanges(bamfile = bam, chromosomes = region, pairedEndReads = T, min.mapq = 10, filterAltAlign = TRUE)
  # source function below for bamregion2Granges
  cell.hp<- fragments[!is.na(fragments$HP)]
  
  c<- which(strand(cell.hp)=="+")
  w <- which(strand(cell.hp)=="-") 
  
  out.df<- as.data.frame(cbind(cell=cell.id, chrom=as.character(seqnames(cell.hp)), start=start(cell.hp), end=end(cell.hp) ))
  out.df$w<-0
  out.df$c<-0
  out.df$c[c] <- 1
  out.df$w[w] <- 1
  out.df$hp <- cell.hp$HP
  
  all.counts<- rbind(all.counts, out.df)
}

write.table(all.counts, file = "~/Documents/Data/Sequencing/AGLCD/P1530_rerun/haplotag/P1530_singleCell_haplotagData.txt", quote = FALSE, row.names = FALSE)

# read in the hap_data generated from haplotagData_plottingTable.R
# d.hap<- as.data.table(read.table("/g/korbel2/StrandSeq/20180628_HaplotaggedData/C7/C7_singleCell_haplotagData.txt", header=T))
d.hap <- read.table("~/Documents/Data/Sequencing/AGLCD/P1530_rerun/haplotag/P1530_singleCell_haplotagData.txt", sep = ' ',header = TRUE)
d.hap[nrow(d.hap)+1,] <- c(cell_ID, chromosome, min(genomic_range), min(genomic_range), 0, 0, 0)
d.hap[nrow(d.hap)+1,] <- c(cell_ID, chromosome, max(genomic_range), max(genomic_range), 0, 0, 0)
d.hap$chrom <- as.factor(d.hap$chrom)
d.hap$start <- as.numeric(d.hap$start)
d.hap$end <- as.numeric(d.hap$end)
d.hap$w <- as.numeric(d.hap$w)
d.hap$c <- as.numeric(d.hap$c)
d.hap$hp <- as.numeric(d.hap$hp)

in.d.hap <- d.hap %>% 
  filter(cell == cell_ID & start %in% (genomic_range+100000) & chrom == chromosome)


#set plotting params up
format_Mb <- function(x) {
  paste(comma(x/1e6), "Mb")
}
theme_update(plot.title = element_text(hjust = 0.5))
bar_width = median(in.d$end - in.d$start)

############################PLOT W:C ratios for genomic_range ###############################
(plt1 <- ggplot(in.d) +
  geom_col(aes(x = ((start+end)/2), y = value, fill = variable), width=bar_width, position = "fill")+
  scale_fill_manual(values=c("paleturquoise4", "sandybrown"), name='Strand')+
  geom_vline(xintercept = c(min(roi),max(roi)), linetype="dotted") +
  #coord_flip(expand = F, ylim=c(-scale_fac*y_lim, scale_fac*y_lim)) +
  coord_flip(expand=F) +
  # formatting
  ggtitle(cell_ID) +
  xlab("Genomic position")+
  ylab("W:C ratio")+
  #scale_x_continuous(limits=c(min(d$x), max(d$x)), breaks = pretty_breaks(10), labels = format_Mb) +
  xlim(c(min(d$x), max(d$x)))+ scale_x_reverse(breaks = pretty_breaks(10), labels = format_Mb) +
  scale_y_continuous(breaks = pretty_breaks(3)) +
  theme_bw() +
  theme(panel.spacing = unit(0.5, "lines"),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        strip.placement = 'outside',
        strip.background = element_rect(fill = NA, colour=NA),
        plot.title = element_text(hjust = 0.5)) + 
  guides(fill = FALSE))
###############################################################################



###############################PLOT read depth#################################
y_lim = max(ind[,quantile(w+c, 0.9),])
bar_width = median(ind$end - ind$start)
(plt2 <- ggplot(ind) +
  aes(x = ((start+end)/2)) +
  geom_hline(aes(yintercept = median(w+c)), col = "gold4", linetype = "dotdash")+ # adds horizontal line for median reads
  geom_bar(aes(y = w+c), width = bar_width, stat='identity', position = 'identity', fill='wheat3') +
  geom_vline(xintercept = c(min(roi),max(roi)), linetype="dotted") +
  #geom_rect(aes(ymin=-(max(w+c))*0.1, ymax=0, xmin=min(start), xmax=max(end)), fill='grey48', alpha=0.75)  + # ideogram bar
  
  # formatting
  #coord_flip(expand = F, ylim=c(-scale_fac*y_lim/2, scale_fac*1.1*y_lim)) + # moves data toward middle of pl
  coord_flip(expand = F, ylim=c(-2, 1.4*1.1*y_lim)) +

  xlab("Genomic Position")+ylab("Depth")+
  xlim(c(min(d$x), max(d$x)))+ scale_x_reverse(breaks = pretty_breaks(10), labels = format_Mb) +
  #scale_x_continuous(limits=c(min(d$x), max(d$x)), breaks = pretty_breaks(10), labels = format_Mb) +
  scale_y_continuous(breaks = pretty_breaks(5)) + 
  theme_bw() +
  ggtitle(cell_ID) +
  theme(panel.spacing = unit(0.5, "lines"),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        strip.placement = 'outside',
        strip.background = element_rect(fill = NA, colour=NA),
        plot.title = element_text(hjust = 0.5)) + 
  guides(fill = FALSE))
###############################################################################


#############################PLOT haplotag lollis##############################
(plt3 <- ggplot(in.d.hap) +
  aes(x = ((start+end)/2)) +
  ### REPLACED WITH LINE+BALL 
  geom_linerange(data=in.d.hap[in.d.hap$w!=0,], aes(ymin=0, ymax=-w, color=hp), size=1.5, alpha=0.3) + # SNPs on W reads (right)
  geom_linerange(data=in.d.hap[in.d.hap$c!=0,], aes(ymin=0, ymax=c, color=hp), size=1.5, alpha=0.3) + # SNPs on C reads (left)
  geom_point(data=in.d.hap[in.d.hap$w!=0,], aes(y=-w, color=hp), size=2) +
  geom_point(data=in.d.hap[in.d.hap$c!=0,],aes(y=c, color=hp), size=2) +
  geom_vline(xintercept = c(min(roi),max(roi)), linetype="dotted") +
  coord_flip(expand = F)  +
  xlim(c(min(in.d.hap$w), max(in.d.hap$w))) +
  ylab("Phase") +
  scale_x_reverse(breaks = pretty_breaks(10), labels = format_Mb) +
  scale_y_continuous(breaks = pretty_breaks(2), limits = c(-1.2,1.2)) +
  theme_bw() +
  ggtitle(cell_ID) +
  theme(panel.spacing = unit(0.5, "lines"),
        axis.title.y = element_blank(),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        strip.placement = 'outside',
        strip.background = element_rect(fill = NA, colour=NA),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
   guides(fill = "none"))


############################triple plot#########################################
dev.off()
plot_grid(plt1,plt2,plt3, nrow = 1, align = "v")
 