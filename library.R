############################################################################
#                                                                          #
#                           Load the R libraries                           #
#                                                                          #
############################################################################

suppressMessages(library(Seurat)) # single-cell RNA-Seq and multimodal analysis tools
suppressMessages(library(Signac)) # single-cell ATAC-Seq or multimodal analysis tools
suppressMessages(library(GenomeInfoDb)) # the genome information annotation
suppressMessages(library(ggplot2)) # for plotting
suppressMessages(library(patchwork)) # parallel computing
suppressMessages(library(hash)) # hash table
suppressMessages(library(dplyr)) # pipe functions
suppressMessages(library(cicero)) # co-accessibility analysis
suppressMessages(library(pbmcapply)) # parallel computing
suppressMessages(library(qs)) # fast loading and saving
suppressMessages(library(Repitools)) # epigenetic tools
suppressMessages(library(GenomicRanges)) # genomic ranges
suppressMessages(library(data.table)) # conveert nested list into data frame
# suppressMessages(library(SeuratWrappers)) # a collection of community-provided methods and extensions for Seurat

# add specific libraries
specLib <- function(org = 'hg38', org.anno = 'EnsDb.Hsapiens.v86', 
                    org.gs = 'BSgenome.Hsapiens.UCSC.hg38', 
                    org.db = 'org.Hs.eg.db') {
  suppressMessages(library(org.anno, character.only = T))
  suppressMessages(library(org.gs, character.only = T))
  suppressMessages(library(org.db, character.only = T))
}