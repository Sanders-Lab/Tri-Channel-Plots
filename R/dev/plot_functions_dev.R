############################PLOT W:C ratios for genomic_range##################
plot_state <- function(in.d, cell_ID) {
  # Print user message
  message("Plotting strand states ...")
  # Grab call to match cell_ID name
  arg <- match.call()
  # Set bar width according to data
  bar_width = median(in.d$end - in.d$start)

  p1 <- ggplot(in.d, aes(x = (start+end)/2, y = value, fill = variable)) + 
    geom_col(width = bar_width, position = "fill") +
    scale_fill_manual(values= c("paleturquoise4", "sandybrown"), 
                      name='Strand') +
    coord_flip(expand=F) +
    ggtitle(arg$cell_ID) +
    xlab("Genomic position") +
    ylab("W:C ratio") +
    scale_x_continuous(breaks = pretty_breaks(15), labels = format_Mb) +
    scale_y_continuous(breaks = c(0,0.5,1.0)) +
    theme_bw() +
    theme(panel.spacing = unit(0.4, "lines"),
          strip.placement = 'outside',
          strip.background = element_rect(fill = NA, colour=NA),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 11)) +
    guides(fill = FALSE)

  return(p1)
}

###############################PLOT read depth#################################
plot_depth <- function(ind, cell_ID) {
  # Print user message
  message("Plotting read depth ...")
  # Set ylim to quantiles of data
  y_lim = quantile(ind$c + ind$w, seq(0,1,0.1))[10] * 1.4
  # Set bar width to data
  bar_width = median(ind$end - ind$start)
  
  p2 <- ggplot(ind, aes(x = ((start+end)/2))) +
    geom_hline(aes(yintercept = median(w+c)),
               col = "gold4",
               linetype = "dotdash") +
    geom_bar(aes(y = w+c),
             width = bar_width,
             stat='identity',
             position = 'identity',
             fill='wheat3') +
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
    p2 <- p2 +
      geom_vline(xintercept = c(min(roi),max(roi)),
                 linetype="dotted",
                 size = 1.1)}
  
  return(p2)
}

######################## PLOT haplotag lollis merged ##########################
plot_phase <- function(in.d.hap, cell_ID) {
  # Print user message
  message("Plotting H merged ...")
  
  p3 <- ggplot(in.d.hap, aes(x = ((start+end)/2))) +
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
    scale_x_reverse(breaks = pretty_breaks(10), labels = format_Mb) +
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
  
  return(p3)
}

#############################PLOT haplotag lollis H1 ##########################
plot_H1 <- function(in.d.hap1, cell_ID) {
  # Print user message
  message("Plotting H1 ...")
  
  p4 <- ggplot(in.d.hap1, aes(x = ((start+end)/2))) +
    geom_linerange(data=in.d.hap1[in.d.hap1$w!=0,],
                   aes(ymin=0, ymax=-w),
                   size=.5,
                   color="red",
                   alpha = 0.3) +
    geom_linerange(data=in.d.hap1[in.d.hap1$c!=0,],
                   aes(ymin=0, ymax=c),
                   size=.5,
                   color="red",
                   alpha = 0.3) +
    geom_point(data=in.d.hap1[in.d.hap1$w!=0,],
               aes(y=-w),
               size=2,
               color="red") +
    geom_point(data=in.d.hap1[in.d.hap1$c!=0,],
               aes(y=c),
               size=2,
               color="red") +
    coord_flip(expand = F)  +
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
  p4 <- p4 +
    geom_vline(xintercept = c(min(roi),max(roi)),
               linetype="dotted",
               size = 1.1)}
  
  return(p4)
}

#############################PLOT haplotag lollis H2 ##########################
plot_H2 <- function(in.d.hap2, cell_ID) {
  # Print user message
  message("Plotting H2 ...")
  
  p5 <- ggplot(in.d.hap2, aes(x = ((start+end)/2))) +
    geom_linerange(data=in.d.hap2[in.d.hap2$w!=0,],
                   aes(ymin=0, ymax=-w),
                   size=.5,
                   color="blue",
                   alpha = 0.3) +
    geom_linerange(data=in.d.hap2[in.d.hap2$c!=0,],
                   aes(ymin=0, ymax=c),
                   size=.5,
                   color="blue",
                   alpha = 0.3) +
    geom_point(data=in.d.hap2[in.d.hap2$w!=0,],
               aes(y=-w),
               size=2,
               color="blue") +
    geom_point(data=in.d.hap2[in.d.hap2$c!=0,],
               aes(y=c),
               size=2,
               color="blue") +
    coord_flip(expand = F)  +
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
  p5 <- p5 +
    geom_vline(xintercept = c(min(roi),max(roi)),
               linetype="dotted",
               size = 1.1)}
  
  return(plt5)
}

