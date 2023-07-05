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