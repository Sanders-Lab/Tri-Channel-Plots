bamregion2GRanges <- function(bamfile, bamindex=bamfile, chromosomes=NULL, pairedEndReads=FALSE, min.mapq=10, filterAltAlign=TRUE) {
  
  ## Check if bamindex exists
  bamindex.raw <- sub('\\.bai$', '', bamindex)
  bamindex <- paste0(bamindex.raw,'.bai')
  if (!file.exists(bamindex)) {
    bamindex.own <- Rsamtools::indexBam(bamfile)
    warning("Couldn't find BAM index-file ",bamindex,". Creating our own file ",bamindex.own," instead.")
    bamindex <- bamindex.own
  }
  
  file.header <- Rsamtools::scanBamHeader(bamfile)[[1]]
  chrom.lengths <- file.header$targets
  chroms.in.data <- names(chrom.lengths)
  if (is.null(chromosomes)) {
    chromosomes <- chroms.in.data
  }
  chroms2use <- intersect(chromosomes, chroms.in.data)
  if (length(chroms2use)==0) {
    chrstring <- paste0(chromosomes, collapse=', ')
    stop('The specified chromosomes ', chrstring, ' do not exist in the data. Please try ', paste(paste0('chr', chromosomes), collapse=', '), ' instead.')
  }
  ## Issue warning for non-existent chromosomes
  diff <- setdiff(chromosomes, chroms.in.data)
  if (length(diff)>0) {
    diffs <- paste0(diff, collapse=', ')
    warning(paste0('Not using chromosomes ', diffs, ' because they are not in the data.'))
  }
  ## Import the file into GRanges
  gr <- GenomicRanges::GRanges(seqnames=chroms2use, ranges=IRanges(start=rep(1, length(chroms2use)), end=chrom.lengths[chroms2use]))
  ## Check if bam is truly paired ended in case pairedEndReads set to TRUE
  is.Paired <- Rsamtools::testPairedEndBam(file = bamfile, index = bamindex)
  if (pairedEndReads) {
    if (!is.Paired) {
      warning("You are trying to process single-ended BAM as paired-ended, Please set proper BAM directioanlity!!!")
    } 
  } else {
    if (is.Paired) {
      warning("You are trying to process paired-ended BAM as single-ended, Please set proper BAM directioanlity!!!")
    }  
  }
  
  ## read in reads data
  if (pairedEndReads) {
    #suppressWarnings( data.raw <- GenomicAlignments::readGAlignmentPairs(bamfile, index=bamindex, param=Rsamtools::ScanBamParam(tag="XA", which=range(region), what=c('seq', 'qual','mapq','cigar'), flag=scanBamFlag(isDuplicate=F))) )
    suppressWarnings( data.raw <- GenomicAlignments::readGAlignmentPairs(bamfile, index=bamindex, param=Rsamtools::ScanBamParam(tag=c("XA","HP"), which=range(gr), what=c('mapq','flag'), )) )	
  } else {
    suppressWarnings( data.raw <- GenomicAlignments::readGAlignments(bamfile, index=bamindex, param=Rsamtools::ScanBamParam(tag=c("XA","HP"), which=range(gr), what=c('mapq'), flag=scanBamFlag(isDuplicate=F))) )
  } 
  
  ## Second mate of the pair will inherit directionality from the first mate of the pair
  if (pairedEndReads) {
    data.first <- as(GenomicAlignments::first(data.raw), 'GRanges')
    data.last <- as(GenomicAlignments::last(data.raw), 'GRanges')
    strand(data.last) <- strand(data.first)
    data <- GenomicRanges::sort(c(data.first, data.last), ignore.strand=TRUE)
  } else {
    data <- as(data.raw, 'GRanges')
  }
  
  ## Filter duplicates for pairedEndReads
  if (pairedEndReads) {
    bit.flag <- bitwAnd(1024, data$flag)
    mask <- bit.flag == 0 	
    data <- data[mask]
  }  
  
  ## Filter by mapping quality
  if (!is.null(min.mapq)) {
    if (any(is.na(mcols(data)$mapq))) {
      warning(paste0(file,": Reads with mapping quality NA (=255 in BAM file) found and removed. Set 'min.mapq=NULL' to keep all reads."))
      mcols(data)$mapq[is.na(mcols(data)$mapq)] <- -1
    }
    data <- data[mcols(data)$mapq >= min.mapq]
  }
  
  ## filter XA tag
  if (filterAltAlign) {
    data <- data[is.na(mcols(data)$XA)]
  }    
  
  #data <- data[seqnames(data) %in% seqlevels(region)]
  #seqlevels(data) <- seqlevels(region)
  data <- GenomeInfoDb::keepSeqlevels(data, seqlevels(gr), pruning.mode="coarse")
  #data <- GenomeInfoDb::keepSeqlevels(data, seqlevels(region))	
  
  return(data)
}


################################################################################

#generate hap_data from haplotagData_plottingTable.R here
haplotaggeR <- function(haplotag.bams.path, chromosomes = paste0("chr", c(1:22, "X")), output = NULL) {
  haplotag.bams <- list.files(path = haplotag.bams.path, pattern = "\\.bam$", full.names = T) #read in bam paths
  #create empty starting data frame for loop
  d.hap <- data.frame()
  pb = txtProgressBar(min = 0, max = length(haplotag.bams), initial = 0, style = 3, title = "Running HaplotaggeR", label = "Progress: ")
  for (i in 1:length(haplotag.bams)) {
    setTxtProgressBar(pb,i, title = "Running HaplotaggeR", label = "Progress: ")
    bam <- haplotag.bams[i]
    filename <- basename(bam)
    cell.id <- unlist(strsplit(filename, "\\."))[1]
    
    fragments <- bamregion2GRanges(bamfile = bam, chromosomes = chromosomes, pairedEndReads = T, min.mapq = 10, filterAltAlign = TRUE)
    # source function below for bamregion2Granges
    cell.hp <- fragments[!is.na(fragments$HP)]
    
    c <- which(strand(cell.hp) == "+")
    w <- which(strand(cell.hp) == "-") 
    
    out.df <- as.data.frame(cbind(cell = cell.id,
                                  chrom = as.character(seqnames(cell.hp)),
                                  start = start(cell.hp),
                                  end = end(cell.hp)))
    out.df$w <- 0
    out.df$c <- 0
    out.df$c[c] <- 1
    out.df$w[w] <- 1
    out.df$hp <- cell.hp$HP
    
    d.hap <- rbind(d.hap, out.df)
  }
  #close progress bar after last iteration to print new lines from here
  close(pb)
  #write table all.counts to outdir after wrangling data types
  d.hap$cell <- as.character(d.hap$cell)
  d.hap$chrom <- as.character(d.hap$chrom)
  d.hap$w <- as.numeric(d.hap$w)
  d.hap$c <- as.numeric(d.hap$c)
  d.hap$start <- as.numeric(d.hap$start)
  d.hap$end <- as.numeric(d.hap$end)
  d.hap$hp <- as.numeric(d.hap$hp)
  return(d.hap)
  if (!is.null(output)) {
    write.table(d.hap, file = output, quote = FALSE, row.names = FALSE)
  }
}
###############################################################################

format_Mb <- function(x) {
  paste(comma(x/1e6), "Mb")
}

###############################################################################
#dev of tri channel plotting function
plot_trichannel <- function(cell_ID, plot_range = NULL, channels = c(4,3,2,1), chromosome = NULL, roi = NULL) {
  
  #check whether count data needs to be wrangled for plotting (safes time in function later i assume)
  message("Setting filters ...")
  if(!"total" %in% colnames(d)) {
    d <- d %>%
      group_by(cell) %>% 
      mutate(total = sum(w+c), mean = median(w+c))
  }
  #Check on channels input
  if (channels > 4) {
    stop("Abort mission! Channels parameter larger than 4 detected!")
  } else {
    message("Channels: ", channels)}
  ########################## Filter for params #################################
  if (is.null(plot_range)) {
    # split haplotypes into plotting df
    message("Fetching haplotype data for cell ", cell_ID, " on ", chromosome, " ...")
    in.d.hap1 <- d.hap %>% 
      filter(cell == cell_ID & chrom == chromosome & hp == 1)
    in.d.hap2 <- d.hap %>% 
      filter(cell == cell_ID & chrom == chromosome & hp == 2)
    
    #filter out roi data here depending on what youre analysing (remove/add arguments as needed)
    message("Fetching count data for ", cell_ID, " on ", chromosome, " ...")
    ind <- d %>% 
      filter(cell == cell_ID & chrom == chromosome)
    in.d <- melt(ind, measure.vars = c("c", "w"))
  } else {
    # split haplotypes into plotting df
    message("Fetching haplotype data for cell ", cell_ID, " in range ", min(plot_range), " to ", max(plot_range), " ...")
    in.d.hap1 <- d.hap %>% 
      filter(cell == cell_ID & start %in% plot_range & end %in% plot_range & hp == 1)
    in.d.hap2 <- d.hap %>% 
      filter(cell == cell_ID & start %in% plot_range & end %in% plot_range & hp == 2)
    
    #filter out roi data here depending on what youre analysing (remove/add arguments as needed)
    message("Fetching count data for ", cell_ID, " in range ", min(plot_range), " to ", max(plot_range), " ...")
    ind <- d %>% 
      filter(cell == cell_ID & start %in% plot_range & end %in% plot_range)
    in.d <- melt(ind, measure.vars = c("c", "w"))
  }
  message("Data fetch complete!")
  ############################PLOT W:C ratios for genomic_range##################
  message("Plotting strand states ...")
  bar_width = median(in.d$end - in.d$start)
  plt1 <- ggplot(in.d) +
    geom_col(aes(x = (start+end)/2, y = value, fill = variable), width=bar_width, position = "fill")+
    scale_fill_manual(values=c("paleturquoise4", "sandybrown"), name='Strand') +
    coord_flip(expand=F) +
    # formatting
    ggtitle(cell_ID) +
    xlab("Genomic position")+
    ylab("W:C ratio")+
    scale_x_continuous(breaks = pretty_breaks(15), labels = format_Mb) +
    scale_y_continuous(breaks = c(0,0.5,1.0)) +
    theme_bw() +
    theme(panel.spacing = unit(0.4, "lines"),
          strip.placement = 'outside',
          strip.background = element_rect(fill = NA, colour=NA),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 11)) +
    guides(fill = FALSE)
  
  if (!is.null(roi)) {
    plt1 <- plt1 +
      geom_vline(xintercept = c(min(roi),max(roi)), linetype="dotted", size = 1.1)
  }
  
  ###############################PLOT read depth#################################
  message("Plotting read depth ...")
  y_lim = quantile(ind$c + ind$w, seq(0,1,0.1))[10] * 1.4
  bar_width = median(ind$end - ind$start)
  plt2 <- ggplot(ind) +
    aes(x = ((start+end)/2)) +
    geom_hline(aes(yintercept = median(w+c)), col = "gold4", linetype = "dotdash")+ # adds horizontal line for median reads
    geom_bar(aes(y = w+c), width = bar_width, stat='identity', position = 'identity', fill='wheat3') +
    # formatting
    coord_flip(expand = F) +
    xlab("Genomic Position")+ylab("Depth") +
    scale_x_continuous(breaks = pretty_breaks(15), labels = format_Mb) +
    scale_y_continuous(breaks = pretty_breaks(5), limits = c(0, y_lim)) + 
    theme_bw() +
    theme(panel.spacing = unit(0.4, "lines"),
          axis.title.y = element_blank(),
          strip.placement = 'outside',
          strip.background = element_rect(fill = NA, colour=NA),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 11)) +
    ggtitle(cell_ID) +
    guides(fill = FALSE)
  
  if (!is.null(roi)) {
    plt2 <- plt2 +
      geom_vline(xintercept = c(min(roi),max(roi)), linetype="dotted", size = 1.1)
  }
  
  #############################PLOT haplotag lollis H1 ##########################
  message("Plotting H1 ...")
  plt3 <- ggplot(in.d.hap1) +
    aes(x = ((start+end)/2)) +
    ### REPLACED WITH LINE+BALL 
    geom_linerange(data=in.d.hap1[in.d.hap1$w!=0,], aes(ymin=0, ymax=-w), size=.5, color="red", alpha = 0.5) + # SNPs on W reads (right)
    geom_linerange(data=in.d.hap1[in.d.hap1$c!=0,], aes(ymin=0, ymax=c),  size=.5, color="red", alpha = 0.5) + # SNPs on C reads (left)
    geom_point(data=in.d.hap1[in.d.hap1$w!=0,], aes(y=-w), size=2, color="red") +
    geom_point(data=in.d.hap1[in.d.hap1$c!=0,], aes(y=c),  size=2, color="red") +
    coord_flip(expand = F)  +
    #formatting
    ylab("H1") +
    scale_x_continuous(breaks = pretty_breaks(15), labels = format_Mb) +
    scale_y_continuous(breaks = pretty_breaks(2), limits = c(-1.2,1.2)) +
    theme_bw() +
    theme(panel.spacing = unit(0.4, "lines"),
          axis.title.y = element_blank(),
          strip.placement = 'outside',
          strip.background = element_rect(fill = NA, colour=NA),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 11)) +
    ggtitle(cell_ID) +
    guides(fill = "none")
  
  if (!is.null(roi)) {
    plt3 <- plt3 +
      geom_vline(xintercept = c(min(roi),max(roi)), linetype="dotted", size = 1.1)
  }
  
  #############################PLOT haplotag lollis H2 ##########################
  message("Plotting H2 ...")
  plt4 <- ggplot(in.d.hap2) +
    aes(x = ((start+end)/2)) +
    ### REPLACED WITH LINE+BALL 
    geom_linerange(data=in.d.hap2[in.d.hap2$w!=0,], aes(ymin=0, ymax=-w), size=.5, color="blue", alpha = 0.5) + # SNPs on W reads (right)
    geom_linerange(data=in.d.hap2[in.d.hap2$c!=0,], aes(ymin=0, ymax=c), size=.5, color="blue", alpha = 0.5) + # SNPs on C reads (left)
    geom_point(data=in.d.hap2[in.d.hap2$w!=0,], aes(y=-w), size=2, color="blue") +
    geom_point(data=in.d.hap2[in.d.hap2$c!=0,],aes(y=c), size=2, color="blue") +
    coord_flip(expand = F)  +
    #formatting
    ylab("H2") +
    scale_x_continuous(breaks = pretty_breaks(15), labels = format_Mb) +
    scale_y_continuous(breaks = pretty_breaks(2), limits = c(-1.2,1.2)) +
    theme_bw() +
    theme(panel.spacing = unit(0.4, "lines"),
          axis.title.y = element_blank(),
          strip.placement = 'outside',
          strip.background = element_rect(fill = NA, colour=NA),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 11)) +
    ggtitle(cell_ID) +
    guides(fill = "none")
  
  if (!is.null(roi)) {
    plt4 <- plt4 +
      geom_vline(xintercept = c(min(roi),max(roi)), linetype="dotted", size = 1.1)
  }
  
  ############################triple plot########################################
  if (channels == 1) {
    tri_plot <- plot_grid(plt1, nrow = 1, align = "v")
    return(tri_plot)
  } else {
    if (channels == 2) {
      tri_plot <- plot_grid(plt1,plt2, nrow = 1, align = "v")
      return(tri_plot)
    } else {
      if (channels == 3) {
        tri_plot <- plot_grid(plt1,plt2,plt3, nrow = 1, align = "v")
        return(tri_plot)
      } else {
        if (channels == 4) {
          tri_plot <- plot_grid(plt1,plt2,plt3,plt4, nrow = 1, align = "v")
          return(tri_plot)
        }
      }
    }
  }
}
