############################################################################
#                                                                          #
#                               Define functions                           #
#                                                                          #
############################################################################


# Suppress the messages
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 


# Convert GRanges into data frame
gr2dt <- function(gr, basic = FALSE) {
  if (any(class(gr) == 'data.table'))
    return(gr)
  out <- with(gr, data.table(seqnames = as.character(seqnames(gr)),
                             start=start(gr), end=end(gr), strand = as.character(strand(gr))))
  if (!basic && ncol(mcols(gr)))
    out <- cbind(out, as.data.frame(mcols(gr)))
  
  return(out)
}


# identify TF-peak relations
# TBD: convert strings into GRange objects
add_TF <- function(peaks, hg.mus.jaspar.list, org = "hg38") {
  
  # peaks : the peaks
  # hg.mus.jaspar.list : TF binding sites in JASPAR database for human or mouse
  # org : the organism, "hg38" by default
  
  if (grepl("mm", org)) {
    jaspar.sites <- hg.mus.jaspar.list[["Mouse"]]
  } else {
    jaspar.sites <- hg.mus.jaspar.list[["Human"]]
  }
  
  overlap <- findOverlaps(query = StringToGRanges(peaks), subject = jaspar.sites$peak)
  # overlap <- findOverlaps(query = strToGR(peaks), subject = strToGR(jaspar.sites$peak))
  
  if (length(overlap) < 1) {
    stop ("--------> No TF is found to bind the peaks.\n")
  }
  
  TF.peak.df <- do.call(rbind.data.frame, pbmclapply(seq_along(overlap), function(x) {
    return(data.frame(TF = jaspar.sites$TF[overlap@to[x]], peak = peaks[overlap@from[x]]))
  }, mc.cores = detectCores())) # build the TF-peak data frame
  message ("--------> Overlapping the TF binding sites with JASPAR sites ...")
  
  peak.dfs <- split(x = TF.peak.df$TF, f = TF.peak.df$peak)
  # nest list, where the names are peaks and the elements are TF lists
  
  TF.dfs <- split(x = TF.peak.df$peak, f = TF.peak.df$TF)
  # nest list, where the names are TFs and the elements are peak lists
  
  peak.TFs <- I(mclapply(peaks, function(x) {
    ifelse (x %in% names(peak.dfs), return(unique(peak.dfs[[x]])), return('NA'))
  }, mc.cores = detectCores())) # annotate peaks using TFs
  names(peak.TFs) <- peaks # add the peak names
  
  TF.peaks <- I(mclapply(names(TF.dfs), function(x) {
    return(unique(TF.dfs[[x]]))
  }, mc.cores = detectCores())) # annotate peaks using TFs
  TF.peaks[[length(TF.peaks) + 1]] <- names(peak.TFs[which(peak.TFs == "NA")])
  names(TF.peaks) <- c(names(TF.dfs), "NA") # add the TF names
  
  return(list(peak.hash = peak.TFs, TF.hash = TF.peaks))
}
