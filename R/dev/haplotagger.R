#generate hap_data from haplotagData_plottingTable.R here
haplotagger <- function(haplotag.bams.path, chromosomes = paste0("chr", c(1:22, "X")), output = NULL) {
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
