# Generate peak-gene linkage pairs
generate_link_pairs <- function(pbmc, distance = 500000, peak.assay = "ATAC") {
  
  # pbmc : Seurat Object
  # distance : the distance to build peak-gene linkages, 500000 by default
  # peak.assay : the assay of chromatin accessibility, ATAC by default
  
  # calculate the nearby genes
  gene.coords <- CollapseToLongestTranscript(ranges = Annotation(object = pbmc[[peak.assay]]))
  # get gene coordinates
  
  peaks <- StringToGRanges(rownames(pbmc[[peak.assay]])) # get peak coordinates
  distance.df <- summary(DistanceToTSS(peaks = peaks, genes = gene.coords, 
                                       distance = distance)) # distance matrix
  
  # convert link ids into peak-gene relations
  peak.names <- rownames(pbmc[[peak.assay]]) # get peak names
  gene.names <- gene.coords$gene_name # get gene names
  
  # links.df <- rbindlist(pblapply(1:nrow(distance.df), function(i) {
  # message (i, "\n")
  return(rbindlist(pbmclapply(1:nrow(distance.df), function(i) {
    return(list(peak = peak.names[distance.df[i, 1]], gene = gene.names[distance.df[i, 2]]))
  }, mc.cores = min(detectCores(), nrow(distance.df))), fill = T) %>% 
    dplyr::filter(gene %in% rownames(pbmc[["RNA"]])))
}


# Calculate Pearson correlation coefficients between 


# Get the most coherent peak-gene pairs
get_coherent_peak_gene_pairs <- function(peak_distance_matrix, HBC.rna, HBC.atac) {
  
  # peak_distance_matrix : the peak-gene distances
  # HBC.rna : the RNA matrix
  # HBC.atac : the ATAC matrix
  
  row.col <- which(peak_distance_matrix > 0, arr.ind = T) # get the nonzero elements
  
  peak.gene <- rbindlist(apply(row.col, 1, function(rr) {
    return(list(rownames(peak_distance_matrix)[rr[1]], colnames(peak_distance_matrix)[rr[2]]))
  })) # convert row and column ids into peaks and rows
  colnames(peak.gene) <- c("peak", "gene")
  
  peak.gene.weight <- cbind(peak.gene, weight = apply(peak.gene, 1, function(rr) {
    # return (length(HBC.rna[rr[2], ]))
    return(length(which(HBC.atac[rr[1], ] > 0 & HBC.rna[rr[2], ] > 0)))
  })) # add weights
  
  peak.gene.weight <- peak.gene.weight[with(peak.gene.weight, order(gene, -weight))] # sort the data frame
  
  return(tapply(peak.gene.weight$peak, peak.gene.weight$gene, function(i) {
    head(i, n = 1) # select the pairs of peaks and genes
  }))
}


LinksToGRanges <- function(linkmat, gene.coords, sep = c("-", "-"))
{
  # linkmat     : the correlation matrix between peaks and genes
  # gene.coords : the gene coordinates
  # sep         : the separators to segment peaks
  
  tss <- resize(gene.coords, width = 1, fix = 'start')# the TSS
  gene.idx <- sapply(
    X = rownames(x = linkmat), # the gene names
    FUN = function(x) {
      which(x = x == tss$gene_name)[[1]]
    }
  )
  tss <- tss[gene.idx] # select the TSS for each gene
  
  # get midpoint of each peak
  peak.ranges <- StringToGRanges(
    regions = colnames(x = linkmat),
    sep = sep
  )
  midpoints <- start(x = peak.ranges) + (width(x = peak.ranges) / 2) # the midpoint of the peaks
  
  # convert to triplet form
  dgtm <- as(object = linkmat, Class = "dgTMatrix")
  
  # create data frame
  df <- data.frame(
    chromosome = as.character(x = seqnames(x = peak.ranges)[dgtm@j + 1]),
    tss = start(x = tss)[dgtm@i + 1],
    pk = midpoints[dgtm@j + 1],
    score = dgtm@x,
    gene = rownames(x = linkmat)[dgtm@i + 1],
    peak = colnames(x = linkmat)[dgtm@j + 1]
  )
  
  # work out start and end coords (I did not understand these lines)
  df$start <- ifelse(test = df$tss < df$pk, yes = df$tss, no = df$pk)
  df$end <- ifelse(test = df$tss < df$pk, yes = df$pk, no = df$tss)
  df$tss <- NULL
  df$pk <- NULL
  
  # convert to granges
  gr.use <- makeGRangesFromDataFrame(df = df, keep.extra.columns = TRUE)
  return(sort(x = gr.use))
}


# I will search for the usage
CollapseToLongestTranscript <- function(ranges) {

  # ranges : the annotation of peaks in a Seurat object

  range.df <- as.data.table(x = ranges) # transform a GRanges object into a data frame
  range.df$strand <- as.character(x = range.df$strand) # transform GRanges into character strings
  range.df$strand <- ifelse( # check whether the strand information is available
    test = range.df$strand == "*", # ambiguous
    yes = "+", # treat all sequences by only using their positive strands
    no = range.df$strand # use the provided strands
  )
  collapsed <- range.df[
    , .(unique(seqnames),
        min(start),
        max(end),
        strand[[1]],
        gene_biotype[[1]]),
    "gene_name"
    ] # merge exons into the longest transcripts as a representative of the gene
  colnames(x = collapsed) <- c(
    "gene_name", "seqnames", "start", "end", "strand", "gene_biotype"
  )
  gene.ranges <- makeGRangesFromDataFrame(
    df = collapsed,
    keep.extra.columns = TRUE # the information not used to build the data frame
    # will be retained in meta data
  )
  return(gene.ranges)
}


# get the distance to TSS
DistanceToTSS <- function(peaks, genes, distance = 200000, sep = c("-", "-")) {

  # peaks    : the peak list
  # genes    : the gene coordinates
  # distance : the distance cutoff to retain peaks near the gene TSSs
  # sep      : the separators to parse the peak coordinates

  tss <- resize(x = genes, width = 1, fix = 'start') # find the TSS location
  genes.extended <- suppressWarnings(
    expr = Extend(
      x = tss, upstream = distance, downstream = distance
    )
  ) # extand the genomic range from the TSS till downstream/upstream 200000 bp
  overlaps <- findOverlaps(
    query = peaks,
    subject = genes.extended,
    type = 'any',
    select = 'all'
  ) # find the peaks overlaped with the extended genomic ranges of genes
  hit_matrix <- sparseMatrix(
    i = queryHits(x = overlaps),
    j = subjectHits(x = overlaps),
    x = 1,
    dims = c(length(x = peaks), length(x = genes.extended))
  ) # build a sparse matrix to record the overlaps between peaks and extended genomic ranges of genes
  rownames(x = hit_matrix) <- GRangesToString(grange = peaks, sep = sep) # use peak names as the row names
  colnames(x = hit_matrix) <- genes.extended$gene_name # use gene names as the column names
  return(hit_matrix)
}


# extract the peak-peak linkages from Cicero
cicero_peak_to_peak <- function(coaccess, peak.coords, cutoff = 0) {
  
  # coaccess : the peak-peak co-accessibility linkages saved in a three-column data frame
  # peak.coords : the list of peak coordinates
  # cutoff : the cutoff for co-accessibility scores, 0 by default
  
  coaccess.use <- coaccess[!is.na(coaccess$coaccess) & # admit no N/A values
                             coaccess$Peak1 %in% peak.coords & 
                             coaccess$Peak2 %in% peak.coords, ]
  # select linkages using coaccessibility score and whether the two peaks are
  # from atac, respectively
  
  link.set <- rbindlist(apply(coaccess.use, 1, function(i) {
    ifelse (i[1] > i[2], return(list(Peak1 = i[1], Peak2 = i[2], coaccess = i[3])), 
            return(list(Peak1 = i[2], Peak2 = i[1], coaccess = i[3])))
    
    # ifelse (i$Peak1 > i$Peak2, return(list(Peak1 = i$Peak1, Peak2 = i$Peak2, coaccess = i$coaccess)), 
    #         return(list(Peak1 = i$Peak2, Peak2 = i$Peak1, coaccess = i$coaccess)))
  }), fill = T) # rearrange the order of Peak1 and Peak2
  
  # check the character strings of peak names
  link.set$Peak1 <- gsub(pattern = '_', replacement = '-', x = link.set$Peak1) # replace '_' with '-
  link.set$Peak2 <- gsub(pattern = '_', replacement = '-', x = link.set$Peak2) # replace '_' with '-
  link.set$coaccess <- as.numeric(link.set$coaccess) # transform characters into numeric values
  
  return(link.set[!duplicated(link.set[ , 1:2]) & 
                    link.set$coaccess >= cutoff, 1:2]) # retain the effective rows
}


# # extract the peak-gene linkages from Cicero
# cicero_peak_to_gene <- function(coaccess, gene.coords, peak.coords, cutoff = 0, coord.hash) {
#   
#   # coaccess : the peak-peak co-accessibility linkages saved in a three-column data frame
#   # gene.coords : the list of gene coordinates
#   # peak.coords : the list of peak coordinates
#   # cutoff : the cutoff for co-accessibility scores
#   # coord.hash : a hash table to map gene symbols to coordinates
#   
#   coaccess.use <- coaccess[!is.na(coaccess$coaccess) & # admit no N/A values
#                              (coaccess$Peak1 %in% gene.coords & 
#                                 coaccess$Peak2 %in% peak.coords | 
#                                 coaccess$Peak2 %in% gene.coords & 
#                                 coaccess$Peak1 %in% peak.coords), ]
#   # select linkages using coaccessibility score and whether the two peaks are
#   # from rna and atac, respectively
#   
#   # if (nrow(coaccess.use) < 1) {
#   #   cat('No peak-gene linkage is found via Cicero.\n')
#   #   
#   #   return(NULL)
#   # }
#   
#   # coaccess.use$Peak2 <- as.character(coaccess.use$Peak2) # format reansformation
#   symbol.hash <- invert(coord.hash)
#   # create a new hash table by inverting keys and values
#   peak <- sapply(1:nrow(coaccess.use), function(i) {
#     return(ifelse(coaccess.use[i, 1] %in% peak.coords, coaccess.use[i, 1], 
#                   coaccess.use[i, 2]))
#   }) # the vector to hold peaks
#   gene <- sapply(1:nrow(coaccess.use), function(i) {
#     return(ifelse(coaccess.use[i, 1] %in% peak.coords, symbol.hash[[coaccess.use[i, 2]]], 
#                   symbol.hash[[coaccess.use[i, 1]]]))
#   }) # the vector to hold genes
#   
#   link.set <- data.frame(peak = peak, gene = gene, coaccess = coaccess.use$coaccess)
#   # transform the gene coordinates into gene symbols
#   
#   link.set$peak <- gsub(pattern = '_', replacement = '-', x = link.set$peak) # replace '_' with '-
#   
#   return(link.set[!duplicated(link.set[ , 1:2]) & 
#                     link.set$coaccess >= cutoff, 1:2]) # retain the effective rows
#   # return(link.set[abs(link.set$coaccess) > cutoff])
# }


# link peaks to genes using heuristics
LinkPeaks2 <- function(object, peak.assay = "ATAC", expression.assay = "RNA", expression.slot = "data", 
                       gene.coords = NULL, distance = 5e+05, min.distance = NULL, 
                       min.cells = 5, method = "pearson", genes.use = NULL, n_sample = 200, 
                       pvalue_cutoff = 0.05, score_cutoff = 0.05, verbose = TRUE) {
  
  # object      : the discretized matrices
  # distance     : the size of the window to build linkages
  # signac.score : the cutoff to select linkages, 0 by default.
  # method : method to calculate correlation, 'Pearson' by default
  # n.sample : the number of samples to select, 200 by default
  # signac.pval : the cutoff of p-values, 0.5 by default
  # org.gs     : the genome annotations, e.g., BSgenome.Hsapiens.UCSC.hg38
  # min.cells : the minimum number of cells in which genes or peaks are expressed 
  #             or accessible, 10 by default
  # cell.weight : use the number of cells to calculate weight, T by default
  # peak.assay : the assay to save the ATAC data, 'ATAC' by default; 'peaks' can
  #              also be used
  
  if (!inherits(x = object[[peak.assay]], what = "ChromatinAssay")) {
    stop("--------> The requested assay is not a ChromatinAssay")
  }
  
  if (!is.null(x = min.distance)) {
    if (!is.numeric(x = min.distance)) {
      stop("--------> min.distance should be a numeric value")
    }
    
    if (min.distance < 0) {
      warning("--------> Requested a negative min.distance value, setting min.distance to zero")
      min.distance <- NULL
    } else if (min.distance == 0) {
      min.distance <- NULL
    }
  }
  
  if (is.null(x = gene.coords)) {
    gene.coords <- CollapseToLongestTranscript(ranges = Annotation(object = object[[peak.assay]]))
  }
  
  meta.features <- GetAssayData(object = object, assay = peak.assay, 
                                slot = "meta.features")
  features.match <- c("GC.percent", "count")
  
  if (!("GC.percent" %in% colnames(x = meta.features))) {
    stop("--------> GC content per peak has not been computed.\n", 
         "Run RegionStats before calling this function.")
  }
  
  peak.data <- GetAssayData(object = object, assay = peak.assay, 
                            slot = "counts")
  
  if (!("count" %in% colnames(x = meta.features))) {
    hvf.info <- FindTopFeatures(object = peak.data)
    hvf.info <- hvf.info[rownames(x = meta.features), ]
    meta.features <- cbind(meta.features, hvf.info)
  }
  
  expression.data <- GetAssayData(object = object, assay = expression.assay, 
                                  slot = expression.slot)
  peakcounts <- meta.features[rownames(x = peak.data), "count"]
  genecounts <- rowSums(x = expression.data > 0)
  peaks.keep <- peakcounts > min.cells # top-ranked peaks, e.g., 3000
  genes.keep <- genecounts > min.cells
  peak.data <- peak.data[peaks.keep, ]
  
  if (is.null(x = genes.use)) {
    expression.data <- expression.data[genes.keep, ]
  } else {
    genes.keep <- intersect(x = names(x = genes.keep[genes.keep]), 
                            y = genes.use)
    expression.data <- expression.data[genes.keep, , drop = FALSE]
  }
  
  if (verbose) {
    message("--------> Testing ", nrow(x = expression.data), " genes and ", 
            sum(peaks.keep), " peaks")
  }
  
  genes <- rownames(x = expression.data)
  gene.coords.use <- gene.coords[gene.coords$gene_name %in% 
                                   genes, ]
  # peaks <- granges(x = object[[peak.assay]])
  peaks <- StringToGRanges(rownames(object[[peak.assay]]))
  peaks <- peaks[peaks.keep]
  peak_distance_matrix <- DistanceToTSS(peaks = peaks, genes = gene.coords.use, 
                                        distance = distance)
  
  if (!is.null(x = min.distance)) {
    peak_distance_matrix_min <- DistanceToTSS(peaks = peaks, 
                                              genes = gene.coords.use, distance = min.distance)
    peak_distance_matrix <- peak_distance_matrix - peak_distance_matrix_min
  }
  
  if (sum(peak_distance_matrix) == 0) {
    stop("No peaks fall within distance threshold\n", "Have you set the proper genome and seqlevelsStyle for ", 
         peak.assay, " assay?")
  }
  
  genes.use <- colnames(x = peak_distance_matrix)
  all.peaks <- rownames(x = peak.data)
  peak.data <- t(x = peak.data)
  coef.vec <- c()
  gene.vec <- c()
  zscore.vec <- c()
  
  if (nbrOfWorkers() > 1) {
    mylapply <- future_lapply
  } else {
    mylapply <- ifelse(test = verbose, yes = pblapply, no = lapply)
  }
  
  res <- mylapply(X = seq_along(along.with = genes.use), FUN = function(i) {
    # message (i, "\n")
    peak.use <- as.logical(x = peak_distance_matrix[, genes.use[[i]]])
    gene.expression <- t(x = expression.data[genes.use[[i]], 
                                             , drop = FALSE])
    gene.chrom <- as.character(x = seqnames(x = gene.coords.use[i]))
    
    if (sum(peak.use) < 2) {
      return(list(gene = NULL, coef = NULL, zscore = NULL))
    } else {
      peak.access <- peak.data[, peak.use, drop = FALSE]
      coef.result <- corSparse(X = peak.access, Y = gene.expression)
      rownames(x = coef.result) <- colnames(x = peak.access)
      coef.result <- coef.result[abs(x = coef.result) > 
                                   score_cutoff, , drop = FALSE]
      
      if (nrow(x = coef.result) == 0) {
        return(list(gene = NULL, coef = NULL, zscore = NULL))
      } else {
        peaks.test <- rownames(x = coef.result)
        trans.peaks <- all.peaks[!grepl(pattern = paste0("^", 
                                                         gene.chrom), x = all.peaks)]
        meta.use <- meta.features[trans.peaks, ]
        pk.use <- meta.features[peaks.test, ]
        bg.peaks <- lapply(X = seq_len(length.out = nrow(x = pk.use)), 
                           FUN = function(x) {
                             MatchRegionStats(meta.feature = meta.use, 
                                              query.feature = pk.use[x, , drop = FALSE], 
                                              features.match = c("GC.percent", "count", 
                                                                 "sequence.length"), 
                                              n = n_sample, verbose = FALSE)
                           })
        n_sample <- min(n_sample, min(sapply(bg.peaks, length)))
        bg.access <- peak.data[, unlist(x = bg.peaks), 
                               drop = FALSE]
        bg.coef <- corSparse(X = bg.access, Y = gene.expression)
        rownames(bg.coef) <- colnames(bg.access)
        zscores <- vector(mode = "numeric", length = length(x = peaks.test))
        
        for (j in seq_along(along.with = peaks.test)) {
          coef.use <- bg.coef[(((j - 1) * n_sample) + 
                                 1):(j * n_sample), ]
          z <- (coef.result[j] - mean(x = coef.use))/sd(x = coef.use)
          zscores[[j]] <- z
        }
        
        names(x = coef.result) <- peaks.test
        names(x = zscores) <- peaks.test
        zscore.vec <- c(zscore.vec, zscores)
        gene.vec <- c(gene.vec, rep(i, length(x = coef.result)))
        coef.vec <- c(coef.vec, coef.result)
      }
      
      gc(verbose = FALSE)
      pval.vec <- pnorm(q = -abs(x = zscore.vec))
      links.keep <- pval.vec < pvalue_cutoff
      
      if (sum(x = links.keep) == 0) {
        return(list(gene = NULL, coef = NULL, zscore = NULL))
      } else {
        gene.vec <- gene.vec[links.keep]
        coef.vec <- coef.vec[links.keep]
        zscore.vec <- zscore.vec[links.keep]
        
        return(list(gene = gene.vec, coef = coef.vec, 
                    zscore = zscore.vec))
      }
    }
  })
  
  gene.vec <- do.call(what = c, args = lapply(X = res, FUN = `[[`, 
                                              1))
  coef.vec <- do.call(what = c, args = lapply(X = res, FUN = `[[`, 
                                              2))
  zscore.vec <- do.call(what = c, args = lapply(X = res, FUN = `[[`, 
                                                3))
  if (length(x = coef.vec) == 0) {
    if (verbose) {
      message("No significant links found")
    }
    
    return(object)
  }
  
  peak.key <- seq_along(along.with = unique(x = names(x = coef.vec)))
  names(x = peak.key) <- unique(x = names(x = coef.vec))
  coef.matrix <- sparseMatrix(i = gene.vec, j = peak.key[names(x = coef.vec)], 
                              x = coef.vec, dims = c(length(x = genes.use), max(peak.key)))
  rownames(x = coef.matrix) <- genes.use
  colnames(x = coef.matrix) <- names(x = peak.key)
  links <- LinksToGRanges(linkmat = coef.matrix, gene.coords = gene.coords.use)
  z.matrix <- sparseMatrix(i = gene.vec, j = peak.key[names(x = zscore.vec)], 
                           x = zscore.vec, dims = c(length(x = genes.use), max(peak.key)))
  rownames(x = z.matrix) <- genes.use
  colnames(x = z.matrix) <- names(x = peak.key)
  z.lnk <- LinksToGRanges(linkmat = z.matrix, gene.coords = gene.coords.use)
  links$zscore <- z.lnk$score
  links$pvalue <- pnorm(q = -abs(x = links$zscore))
  links <- links[links$pvalue < pvalue_cutoff]
  Links(object = object[[peak.assay]]) <- links
  
  return(object)
}


# # Match DNA sequence characteristics
# MatchRegionStats2 <- function(meta.feature, query.feature, features.match = c("GC.percent"), 
#           n = 10000, verbose = TRUE, ...) {
#   if (!inherits(x = meta.feature, what = "data.frame")) {
#     stop("meta.feature should be a data.frame")
#   }
#   if (!inherits(x = query.feature, what = "data.frame")) {
#     stop("query.feature should be a data.frame")
#   }
#   if (length(x = features.match) == 0) {
#     stop("Must supply at least one sequence characteristic to match")
#   }
#   if (nrow(x = meta.feature) < n) {
#     n <- nrow(x = meta.feature)
#     warning("Requested more features than present in supplied data.\n            Returning ", 
#             n, " features")
#   }
#   for (i in seq_along(along.with = features.match)) {
#     featmatch <- features.match[[i]]
#     if (!(featmatch %in% colnames(x = query.feature))) {
#       if (featmatch == "GC.percent") {
#         stop("GC.percent not present in meta.features.", 
#              " Run RegionStats to compute GC.percent for each feature.")
#       }
#       else {
#         stop(featmatch, " not present in meta.features")
#       }
#     }
#     if (verbose) {
#       message("Matching ", featmatch, " distribution")
#     }
#     density.estimate <- density(x = query.feature[[featmatch]], 
#                                 kernel = "gaussian", bw = 1)
#     weights <- approx(x = density.estimate$x, y = density.estimate$y, 
#                       xout = meta.feature[[featmatch]], yright = 1e-04, 
#                       yleft = 1e-04)$y %>% na.approx
#     if (i > 1) {
#       feature.weights <- feature.weights * weights
#     }
#     else {
#       feature.weights <- weights
#     }
#   }
#   feature.select <- sample.int(n = nrow(x = meta.feature), 
#                                size = n, prob = feature.weights)
#   feature.select <- rownames(x = meta.feature)[feature.select]
#   return(feature.select)
# }


# link peaks to genes using Signac
link_signac <- function(x, distance = 500000, 
                        signac.score = 0, 
                        method = 'pearson', n.sample = 200,  
                        signac.pval = 1, 
                        min.cells = 10, cell.weight = T, 
                        peak.assay = 'ATAC') {
  
  # x          : the discretized matrices
  # distance     : the size of the window to build linkages
  # signac.score : the cutoff to select linkages, 0 by default.
  # method : method to calculate correlation, 'Pearson' by default
  # n.sample : the number of samples to select, 200 by default
  # signac.pval : the cutoff of p-values, 0.5 by default
  # org.gs     : the genome annotations, e.g., BSgenome.Hsapiens.UCSC.hg38
  # min.cells : the minimum number of cells in which genes or peaks are expressed 
  #             or accessible, 10 by default
  # cell.weight : use the number of cells to calculate weight, T by default
  # peak.assay : the assay to save the ATAC data, 'ATAC' by default; 'peaks' can
  #              also be used
  
  DefaultAssay(x) <- peak.assay # set 'ATAC' as he default assay
  xx <- x # back up
  x <- tryCatch(LinkPeaks2(object = x, distance = distance, 
                                        min.cells = min.cells,
                 peak.assay = peak.assay, expression.assay = 'RNA', 
                 method = method, n_sample = n.sample, pvalue_cutoff = signac.pval, 
                 score_cutoff = signac.score, verbose = T), 
                 error = function(e) {
                   0 }) # build linkages
  
  # if (is.numeric(x)) {
  #   x <- tryCatch(LinkPeaks2(object = xx, distance = distance, 
  #                            min.cells = min.cells,
  #                            peak.assay = peak.assay, expression.assay = 'SCT', 
  #                            method = method, n_sample = n.sample, pvalue_cutoff = signac.pval, 
  #                            score_cutoff = signac.score, verbose = T), 
  #                 error = function(e) {
  #                   0 }) # modified function of LinkPeaks
  # } # no peaks are linked to the genes
  
  if (!is.numeric(x)) {
    signac.links <- data.frame(node1 = x[[peak.assay]]@links$peak, 
                               node2 = x[[peak.assay]]@links$gene, 
                               weight = x[[peak.assay]]@links$score)
    # arrange the data frame to hold peak-gene linkages
  } else {
    x <- xx
    signac.links <- generate_link_pairs(pbmc = x) # link peaks to genes using heuristics

    signac.links <- cbind(signac.links, mclapply(1:nrow(signac.links), function(i) {
     vx <- as.vector(x[["RNA"]][signac.links$gene[i]])
     vy <- as.vector(x[[atac.assay]][signac.links$peak[i]])
     
     if (sum(vx > 0) <= 0 | sum(vy > 0) <= 0) {
       return(0)
     } else if (sd(vx) == 0 | sd(vy) == 0) {
       return(1)
     } else {
       return(cor(x = vx, y = vy, method = method))
     }
    }, mc.cores = min(detectCores(), nrow(signac.links))) %>% unlist)
    colnames(signac.links) <- c("node1", "node2", "weight")
    signac.links <- signac.links[signac.links$weight > signac.score,] # filter linkages
    # return(NULL)
  } # no peaks are linked to the genes
  
  
  if (nrow(signac.links) < 1) {
    return(NULL)
  }
  
  # signac.links <- signac.links[signac.links$weight >= signac.score, ]
  # # retain the linkages with weight larger than the cutoff
  
  # calculate weights based on the number of cells
  if (cell.weight) {
    signac.links <- rbindlist(lapply(1:nrow(signac.links), function(i) {
      p <- signac.links$node1[i] # the peak
      g <- signac.links$node2[i] # the gene
      p.vv <- x[[peak.assay]][p, ][1, ] # the accessibility vector of the peak
      g.vv <- x[['RNA']][g, ][1, ] # the expression vector of the gene
      
      return(list(node1 = p, node2 = g, weight = length(names(p.vv)[p.vv > 0 & g.vv > 0])))
    }), fill = T)
  }
  
  max.weight <- max(signac.links$weight)
  min.weight <- min(signac.links$weight)
  diff <- max.weight - min.weight
  
  ifelse(diff > 0, signac.links$weight <- (max.weight - signac.links$weight) /
                                  diff, signac.links$weight <- 0) # normalize the weights
  
  message ('--------> Finished generating ', nrow(signac.links), ' peak-gene linkages.\n')
  
  return(signac.links) # remove negative correlations
}


# # retian the top-ranked or all peak-gene linkages for hybrid biclustering
# retain_cis_links <- function(hybrid.links, var.genes = NULL) {
#   
#   # hybrid.links : the list of peak-gene linkages predicted by Cicero and Signac
#   # var.genes : the highly variable genes, NULL by default
#   
#   if (!is.null(var.genes)) {
#     top.links <- dplyr::intersect(hybrid.links$cicero.cis, hybrid.links$signac) # intersection
#     
#     return(top.links[top.links$gene %in% var.genes, ])
#     # only retain the peak-gene linkages connected to variable genes
#   } else {
#     return(dplyr::union(hybrid.links$cicero.cis, hybrid.links$signac)) # HBC linkages
#   }
# }


# select the linkages with no TF bound
na_links <- function(HBC.links, TF.links) {
  
  # HBC.links : all the peak-gene linkages (two-column data.frame format)
  # TF.links : the peak-gene linkages bound by at least one TF (three-column data.frame)
  
  return(cbind(TF = 'NA', dplyr::setdiff(HBC.links, unique(TF.links[, 2 : 3])))) # return the 
  # peak-gene linkages bound by no TF
}

# get the seeds, i.e., genes, peaks, and cells
get_GPCT <- function(rna.var, atac.var, TF.top.links, i, j, seed.ratio = 0, genes) {
  
  # rna.var : the expression matrix subsetted based on variable genes
  # atac.var : the accessibility matrix composed of peaks linked to the variable genes
  # TF.top.links : the TF-peak-gene relations of variable genes
  # i : the id of the first gene
  # j : the id of the second gene
  # seed.ratio : the consistency cutoff between genes and peaks
  # genes : the set of genes having at least one linked peak in TF.top.links
  
  # avoid redundant comparison
  if (i <= j) {
    return(NA)
    # return(list(gene1 = NA, gene2 = NA, 
    #             cells = NA, TFs = NA, 
    #             peak1 = NA, peak2 = NA))
  }
  
  TFs <- intersect(unique(TF.top.links$TF[TF.top.links$gene == genes[i]]), 
                   unique(TF.top.links$TF[TF.top.links$gene == genes[j]]))
  # the TFs regulating the two genes simultaneously
  
  # if there is no intersected TFs regulating the two genes simultaneously
  if (length(TFs) < 1) {
    return(NA)
  }
  
  i.v <- rna.var[genes[i], ] # the expression profile of the first gene
  j.v <- rna.var[genes[j], ] # the expression profile of the second gene
  # i.v <- rna.var[genes[i], , drop = F] # the expression profile of the first gene
  # j.v <- rna.var[genes[j], , drop = F] # the expression profile of the second gene
  same.cells <- colnames(rna.var)[i.v > 0 & j.v > 0] # the cells where two genes 
  # are coexpressed
  
  # two few cells where the two genes are coexpressed
  if (length(same.cells) < 3) {
    return(NA)
  }
  
  # get the peaks for genes 1 and 2
  peaks1 <- unique(TF.top.links$peak[TF.top.links$gene == genes[i] & 
                                       TF.top.links$TF %in% TFs])
  peaks2 <- unique(TF.top.links$peak[TF.top.links$gene == genes[j] & 
                                       TF.top.links$TF %in% TFs])
  
  # skip the situations that no peaks are linked to any one of the genes
  if (length(peaks1) * length(peaks2) == 0) {
    return(NA)
  }
  
  # initialize an empty data frame
  GPCT.df <- list()
  id <- 0
  
  # GPCT.df <- data.frame(matrix(ncol = 6, nrow = 0))
  # colnames(GPCT) <- c('gene1', 'gene2', 'cells', 'TFs', 
  #                     'peak1', 'peak2')
  
  # enumerate the peak pairs
  for (k in 1 : length(peaks1)) {
    # k.v <- atac.var[peaks1[k], same.cells, drop = F] # the accessibility profile of the first peak
    k.v <- atac.var[peaks1[k], same.cells] # the accessibility profile of the first peak
    
    if (length(k.v[k.v > 0]) / length(same.cells) <= seed.ratio) {
      next
    }
    
    TFs1 <- TF.top.links$TF[TF.top.links$gene == genes[i] & 
                              TF.top.links$peak == peaks1[k]]
    
    for (l in 1 : length(peaks2)) {
      # l.v <- atac.var[peaks2[l], same.cells, drop = F] # the accessibility profile of the first peak
      l.v <- atac.var[peaks2[l], same.cells] # the accessibility profile of the first peak
      
      # the consistency level between peaks and genes are too lower
      if (length(colnames(atac.var[ , same.cells])[k.v > 0 & l.v > 0]) / 
          length(same.cells) <= seed.ratio) {
        next
      }
      
      TFs2 <- TF.top.links$TF[TF.top.links$gene == genes[j] & 
                                TF.top.links$peak == peaks2[l]]
      
      same.TFs <- intersect(TFs1, TFs2) # the intersected TFs
      
      # no common TFs bind the two peaks
      if (length(same.TFs) < 1) {
        next
      }
      
      id <- id + 1
      GPCT.df[[id]] <- list(gene1 = genes[i], gene2 = genes[j], 
                            cells = length(same.cells), TFs = same.TFs[[1]],
                            peak1 = peaks1[k], peak2 = peaks2[l])
      # GPCT.df <- rbind(GPCT.df, c(gene1 = genes[i], gene2 = genes[j], 
      #             cells = length(same.cells), TFs = same.TFs, 
      #             peak1 = peaks1[k], peak2 = peaks2[l])) # all an row
    }
  }
  
  ifelse(length(GPCT.df) < 1, return(NA), return(GPCT.df)) # a n * 6 data frame
  # ifelse(nrow(GPCT.df) < 1, return(NA), return(GPCT.df)) # a n * 6 data frame
}


# seeding generates seed set in the format of nested list
# enumerate gene pairs first and then TF
# return a nested list, each element of which is gene1, gene2, cells, TFs, peak1, and peak2
seeding <- function(rna.var, atac.var, TF.top.links, seed.ratio = 0, parallel = T) {
  
  # rna.var : the RNA matrix containing variable genes
  # atac.var : the ATAC matrix composed of peaks linked to the variable genes
  # TF.top.links : the top-ranked TF-CRE-gene relations associated to variable genes
  # seed.ratio : the consistency ratio between the gene and the peak connected by the same linkage
  # parallel : perform parallel computing, T by default
  
  genes <- unique(TF.top.links$gene) # the genes having at least one linked peak via
  # TF.top.links
  
  if (parallel) {
    # load the libraries for parallel computing
    suppressMessages(library(foreach))
    suppressMessages(library(doParallel))
    cl <- makeCluster(detectCores() - 1)
    registerDoParallel(cl)
    
    GPCT <- foreach (i = 1 : length(genes), .combine = 'c', .packages = 'Matrix', .export = 'get_GPCT') %:%
      foreach (j = 1 : length(genes), .combine = 'c') %dopar% {
        # source(paste0("/fs/ess/PCON0022/liyang/tune_stream/11_no_regulonDB_KL/codes/cisAnalysis.R"))
        get_GPCT(rna.var = rna.var, atac.var = atac.var,
                 TF.top.links = TF.top.links,
                 i = i, j = j, seed.ratio = seed.ratio,
                 genes = genes) # get the TF list based on
      }
    
    stopCluster(cl)
    sort_nested_list(list = GPCT[!is.na(GPCT)], 
                     str = 'cells', decreasing = T)
    # sort nested list according to elements in the third id
  } else {
    GPCT <- list()
    id <- 0
    
    for (i in seq_along(genes)[-1]) {
      for (j in 1:(i - 1)) {
        GPCT.ll <- get_GPCT(rna.var = rna.var, atac.var = atac.var,
                            TF.top.links = TF.top.links,
                            i = i, j = j, seed.ratio = seed.ratio,
                            genes = genes) %>% '[[' (1)
      }
      
      if (is.na(GPCT.ll)) {
        next
      }
      
      id <- id + 1
      GPCT[[id]] <- GPCT.ll
    }
    
    sort_nested_list(list = GPCT, 
                     str = 'cells', decreasing = T)
  }
}