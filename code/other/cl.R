
readAlignmentFromBed <- function(filename, paired){
  if (is.installed('readr')){
    tmp <- readr::read_tsv(file = filename, col_names = FALSE)
  } else{
    tmp <- read.delim(file = filename, col.names = FALSE, sep = "\t", stringsAsFactors = FALSE)
  }
  strand_col <- which(apply(tmp[1:min(100, nrow(tmp)),], 2, function(x) all_true(x %in% c("+","-","*"))))
  if (length(strand_col) == 1){
    tmp <- tmp[,c(1:3,strand_col)]
    colnames(tmp) <- c("chr", "start", "end", "strand")
    tmp[,"start"] <- tmp[,"start"] + 1
    tmp <- with(tmp, GRanges(tmp$chr, ranges = IRanges(tmp$start, tmp$end), strand = strand))
  } else{
    tmp <- tmp[,1:6]
    colnames(tmp) <- c("chr", "start", "end")
    tmp[,"start"] <- tmp[,"start"] + 1
    tmp <- with(tmp, GRanges(tmp$chr, ranges = IRanges(tmp$start, tmp$end)))
  }
  if (paired){
    left <- resize(tmp, width = 1, fix = "start", ignore.strand = TRUE)
    right <- resize(tmp, width = 1, fix = "end", ignore.strand = TRUE)
    out <- left_right_to_grglist(left,right)
  } else{
    out <- resize(tmp, width = 1, ignore.strand = FALSE)
  }
  return(out)
}

left_right_to_grglist <- function(left, right){
  stopifnot(length(left) == length(right))
  if (length(left) == 0){
    return(GenomicRangesList())
  }
  x = c(left,right)[as.vector(matrix(seq_len(2L * length(left)), nrow=2L, byrow=TRUE))]
  p = IRanges::PartitioningByEnd(cumsum(rep(2,length(x)/2)))
  out = BiocGenerics::relist(x,p)
  return(out)
}

#' @export
bamToFragmentsByRG <- function(bamfile, paired){
  
  if (paired){
    scanned <- Rsamtools::scanBam(bamfile,
                                  param = Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isMinusStrand=FALSE,
                                                                                                isProperPair = TRUE),
                                                                  what = c("rname","pos","isize"),
                                                                  tag = "RG"))[[1]]
    RG_tags <- mxsort(unique(scanned$tag$RG))
    
    out <- BiocParallel::bplapply(RG_tags, function(RG){
      match_RG = which(scanned$tag$RG == RG)
      scanned_left <- GRanges(seqnames = scanned$rname[match_RG],
                              IRanges::IRanges(start = scanned$pos[match_RG],
                                               width = 1),
                              strand = "+")
      scanned_right <- GRanges(seqnames = scanned$rname[match_RG],
                               IRanges::IRanges(start = scanned$pos[match_RG] +
                                                  abs(scanned$isize[match_RG]) - 1,
                                                width = 1),
                               strand = "-")
      return(left_right_to_grglist(scanned_left, scanned_right))
    })
  } else{
    scanned <- Rsamtools::scanBam(bamfile,
                                  param = Rsamtools::ScanBamParam(what = c("rname","pos","strand","qwidth"),
                                                                  tag = "RG"))[[1]]
    RG_tags <- mxsort(unique(scanned$tag$RG))
    
    out <- BiocParallel::bplapply(RG_tags, function(RG){
      match_RG = which(scanned$tag$RG == RG)
      return(GRanges(seqnames = scanned$rname[match_RG],
                     IRanges::IRanges(start = ifelse(scanned$strand[match_RG] == "-", scanned$pos[match_RG] + scanned$qwidth[match_RG]-1, scanned$pos[match_RG]),
                                      width = 1)))
    })
  }
  
  names(out) = RG_tags
  
  return(out)
}

#' @export
bamToFragments <- function(bamfile, paired){
  if (paired){
    scanned <- Rsamtools::scanBam(bamfile, param = Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isMinusStrand=FALSE, isProperPair = TRUE), what = c("rname","pos","isize")))[[1]]
    scanned_left <- GRanges(seqnames = scanned$rname, IRanges::IRanges(start = scanned$pos, width = 1), strand = "+")
    scanned_right <- GRanges(seqnames = scanned$rname, IRanges::IRanges(start = scanned$pos + abs(scanned$isize) - 1, width = 1), strand = "-")
    out <- left_right_to_grglist(scanned_left, scanned_right)
  } else{
    scanned <- Rsamtools::scanBam(bamfile,
                                  param = Rsamtools::ScanBamParam(what = c("rname","pos","strand","qwidth")))[[1]]
    out <- GRanges(seqnames = scanned$rname,
                   IRanges::IRanges(start = ifelse(scanned$strand == "-", scanned$pos + scanned$qwidth-1, scanned$pos),
                                    width = 1))
  }
  return(out)
  
}

getFragmentCountsByRG <- function(bam, peaks, paired){
  message(paste("Reading in file: ",bam, sep="",collapse=""))
  rg_fragments <- bamToFragmentsByRG(bam, paired)
  
  tmpfun <- function(frags){
    overlaps = as.data.frame(GenomicRanges::findOverlaps(peaks, frags, type="any", ignore.strand=TRUE))
    return(overlaps)
  }
  
  all_overlaps <-  BiocParallel::bplapply(rg_fragments,tmpfun)
  counts_mat <- sparseMatrix(i = do.call(rbind,all_overlaps)$queryHits,
                             j = unlist(lapply(seq_along(all_overlaps),function(y) rep(y,nrow(all_overlaps[[y]]))),use.names=FALSE),
                             x = 1, dims = c(length(peaks),length(rg_fragments)), dimnames = list(NULL,names(rg_fragments)))
  
  return(list(counts = counts_mat, depths = sapply(rg_fragments, length)))
}

# Read in depths from bam ------------------------------------------------------

#' get_sample_depths
#'
#' makes vector of read depths in bam files or RG groups within bam files
#' @param alignment_files filenames for bam or bed file(s) with aligned reads
#' @param by_rg use RG tags to separate groups?
#' @param paired paired end data?
#' @param format bam or bed format? default is bam
#' @return numeric vector
#' @seealso \code{\link{get_counts}},  \code{\link{get_inputs}}, \code{\link{filter_samples}}
#' @export
get_sample_depths <- function(alignment_files, paired = TRUE, by_rg = FALSE, format = c("bam","bed")){
  format = match.arg(format)
  if (format == "bam"){
    return(get_sample_depths_from_bams(alignment_files, paired, by_rg))
  } else{
    return(get_sample_depths_from_beds(alignment_files))
  }
}
get_sample_depths_from_bams <- function(bams, paired = TRUE, by_rg = FALSE){
  if (by_rg){
    out <- do.call(c, lapply(bams, getSampleDepthsByRG, paired = paired))
  } else{
    if (paired){
      out <- sapply(bams,
                    Rsamtools::countBam,
                    param = Rsamtools::ScanBamParam(Rsamtools::scanBamFlag(isMinusStrand=FALSE,
                                                                           isProperPair = TRUE)))
    } else{
      out <- sapply(bams, Rsamtools::countBam)
    }
    names(out) <- sapply(bams, basename)
  }
  return(out)
}

get_sample_depths_from_beds <- function(beds){
  if (is.installed('readr')){
    out = do.call(c, BiocParallel::bplapply(beds, function(filename) nrow(readr::read_tsv(file = filename, col_names = FALSE))))
  } else{
    out = do.call(c, BiocParallel::bplapply(beds, function(filename) nrow(read.delim(file = filename, header = FALSE, sep = "\t", stringsAsFactors = FALSE))))
  }
  names(out) <- sapply(beds, basename)
  return(out)
}

getSampleDepthsByRG <- function(bamfile, paired = TRUE){
  if (paired){
    tags <- Rsamtools::scanBam(bamfile,
                               param = Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isMinusStrand=FALSE,
                                                                                             isProperPair = TRUE),
                                                               tag = "RG"))[[1]]$tag$RG
  } else{
    tags <- Rsamtools::scanBam(bamfile,
                               param = Rsamtools::ScanBamParam(tag = "RG"))[[1]]$tag$RG
  }
  
  RG_tags <- mxsort(unique(tags))
  out <- tabulate(factor(tags, levels = RG_tags, ordered = TRUE))
  names(out) <- RG_tags
  return(out)
}

get_counts_from_bams_cl <- function(bams, peaks, paired = TRUE, by_rg = FALSE, sample_annotation = NULL){
  
  if (by_rg){
    tmp <- lapply(bams, getFragmentCountsByRG, peaks = peaks, paired = paired)
    if (!is.null(colData) && nrow(sample_annotation) == length(bams)){
      l <- sapply(tmp, function(x) length(x$depths))
      sample_annotation <- do.call(rbind, lapply(seq_along(l), function(x) rep(sample_annotation[x,,drop=FALSE],l[x])))
    }
    counts_mat <- do.call(cBind, lapply(tmp, function(x) x$counts))
    depths <- do.call(c, lapply(tmp, function(x) x$depths))
  } else{
    mat = matrix(nrow = length(peaks), ncol = length(bams))
    depths = vector("numeric", length(bams))
    
    for (i in seq_along(bams)){
      message(paste("Reading in file: ",bams[i], sep="",collapse=""))
      fragments = bamToFragments(bams[i], paired = paired)
      depths[i] = length(fragments)
      mat[,i] = GenomicRanges::countOverlaps(peaks, fragments, type="any", ignore.strand=T)
    }
    colnames(mat) = basename(bams)
    return(mat)
  }
}


get_peaks <- function(filename, extra_cols = c(), sort_peaks = FALSE){
    bed <- as.data.frame(suppressMessages(readr::read_tsv(file = filename, col_names = FALSE)[, c(1:3, extra_cols)]))
  colnames(bed) <- c("chr", "start", "end", names(extra_cols))
  bed[,"start"] <- bed[,"start"] + 1
  bed <- makeGRangesFromDataFrame(bed, keep.extra.columns = TRUE)
  if (!isDisjoint(bed)){
    warning("Peaks are overlapping!
            After getting peak counts, peaks can be reduced to non-overlapping set
              using filter_peaks function")
  }
  if (sum(width(bed) == width(bed[1])) != length(bed)){
    warning('Peaks are not equal width!
            Use resize(peaks, width = x, fix = "center") to make peaks equal in size,
            where x is the desired size of the peaks)')
  }
  bed <- GenomeInfoDb::sortSeqlevels(bed)
  sorted_bed = sort(bed, ignore.strand = TRUE)
  if (sort_peaks){
    if (!isTRUE(all.equal(sorted_bed, bed))){
      message("Peaks sorted")
    }
    return(sorted_bed)
  } else{
    if (!isTRUE(all.equal(sorted_bed, bed))){
      warning("Peaks not sorted")
    }
    return(bed)
  }
}