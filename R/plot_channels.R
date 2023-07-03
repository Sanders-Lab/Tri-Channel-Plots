###############################################################################
source("R/dev/bamregion2ranges.R")
source("R/dev/haplotagger.R")
source("R/dev/misc.R")

#dev of tri channel plotting function
plot_channels <- function(cell_ID, plot_range = NULL, channels = c(4,3,2,1), chromosome = NULL, roi = NULL) {
  
  #check whether count data needs to be wrangled for plotting (safes time in function later i assume)
  message("Setting filters:")
  if(!"total" %in% colnames(d)) {
    d <- d %>%
      group_by(cell) %>% 
      mutate(total = sum(w+c), mean = median(w+c))
  }
  #Check on channels input
  if (channels > 4) {
    stop("Abort mission! Channels parameter larger than 4 detected!")
  } else {
    message("Cell ID      [", cell_ID, "]")
    message("Channels     [", channels, "]")
    message("ROI          [", round(min(roi/1000000), digits = 2),":",
            round(max(roi/1000000), digits = 2)," Mb]")
    message("Chromosome   [", chromosome, "]")}
  ########################## Filter for params #################################
  if (is.null(plot_range)) {
    # split haplotypes into plotting df
    message("Fetching haplotype data for cell ", cell_ID, " on ", chromosome, " ...")
    in.d.hap1 <- d.hap %>% 
      filter(cell == cell_ID & chrom == chromosome & hp == 1)
    in.d.hap2 <- d.hap %>% 
      filter(cell == cell_ID & chrom == chromosome & hp == 2)
    in.d.hap <- d.hap %>% 
      filter(cell == cell_ID & chrom == chromosome)
    #filter out roi data here depending on what youre analysing (remove/add arguments as needed)
    message("Fetching count data for ", cell_ID, " on ", chromosome, " ...")
    ind <- d %>% 
      filter(cell == cell_ID & chrom == chromosome)
    in.d <- melt(ind, measure.vars = c("c", "w"))
  } else {
    # split haplotypes into plotting df
    message("Fetching haplotype data for cell ", cell_ID, " in range ", min(plot_range), " to ", max(plot_range), " ...")
    in.d.hap1 <- d.hap %>% 
      filter(cell == cell_ID &
             chrom == chromosome &
             start > min(plot_range) &
             end < max(plot_range) &
             hp == 1)
    in.d.hap2 <- d.hap %>% 
      filter(cell == cell_ID &
               chrom == chromosome &
               start > min(plot_range) &
               end < max(plot_range) &
               hp == 2)
    in.d.hap <- d.hap %>% 
      filter(cell == cell_ID &
             chrom == chromosome &
             start > min(plot_range) &
             end < max(plot_range))
    #filter out roi data here depending on what youre analysing (remove/add arguments as needed)
    message("Fetching count data for ", cell_ID, " in range ", min(plot_range), " to ", max(plot_range), " ...")
    ind <- d %>% 
      filter(cell == cell_ID &
             chrom == chromosome &
             start > min(plot_range) &
             end < max(plot_range))
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

#############################PLOT haplotag lollis merged ######################
  message("Plotting strand phase ...")
  suppressMessages(
    plt3 <- ggplot(in.d.hap, aes(x = ((start+end)/2))) +
    geom_linerange(data=in.d.hap[in.d.hap$w!=0,],
                   aes(ymin=0, ymax=-w, color=hp),
                   size=1.5,
                   alpha=0.3) +
    geom_linerange(data=in.d.hap[in.d.hap$c!=0,],
                   aes(ymin=0, ymax=c, color=hp),
                   size=1.5,
                   alpha=0.3) +
    geom_point(data=in.d.hap[in.d.hap$w!=0,],
               aes(y=-w, color=hp),
               size=3) +
    geom_point(data=in.d.hap[in.d.hap$c!=0,],
               aes(y=c, color=hp),
               size=3) +
    geom_vline(xintercept = c(min(roi),max(roi)), linetype="dotted") +
    coord_flip(expand = F)  +
    xlim(c(min(in.d.hap$w), max(in.d.hap$w))) +
    ylab("Phase") +
    scale_x_continuous(breaks = pretty_breaks(10), labels = format_Mb) +
    scale_y_continuous(breaks = pretty_breaks(2), limits = c(-1.2,1.2)) +
    scale_color_gradient(low = "red", high =  "blue") +
    theme_bw() +
    ggtitle(cell_ID) +
    theme(panel.spacing = unit(0.5, "lines"),
          axis.title.y = element_blank(),
          strip.placement = 'outside',
          strip.background = element_rect(fill = NA, colour=NA),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5)) +
    guides(fill = "none")
  )
#############################PLOT haplotag lollis H1 ##########################
  message("Plotting H1 ...")
  plt4 <- ggplot(in.d.hap1) +
    aes(x = ((start+end)/2)) +
    ### REPLACED WITH LINE+BALL 
    geom_linerange(data=in.d.hap1[in.d.hap1$w!=0,], aes(ymin=0, ymax=-w), size=.4, color="red", alpha = 0.2) + # SNPs on W reads (right)
    geom_linerange(data=in.d.hap1[in.d.hap1$c!=0,], aes(ymin=0, ymax=c),  size=.4, color="red", alpha = 0.2) + # SNPs on C reads (left)
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
    plt4 <- plt4 +
      geom_vline(xintercept = c(min(roi),max(roi)), linetype="dotted", size = 1.1)
  }
  
#############################PLOT haplotag lollis H2 ##########################
  message("Plotting H2 ...")
  plt5 <- ggplot(in.d.hap2) +
    aes(x = ((start+end)/2)) +
    ### REPLACED WITH LINE+BALL 
    geom_linerange(data=in.d.hap2[in.d.hap2$w!=0,], aes(ymin=0, ymax=-w), size=.4, color="blue", alpha = 0.2) + # SNPs on W reads (right)
    geom_linerange(data=in.d.hap2[in.d.hap2$c!=0,], aes(ymin=0, ymax=c), size=.4, color="blue", alpha = 0.2) + # SNPs on C reads (left)
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
    plt5 <- plt5 +
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
          tri_plot <- plot_grid(plt1,plt2,plt4,plt5, nrow = 1, align = "v")
          return(tri_plot)
        }
      }
    }
  }
}
