############################################################################
#                                                                          #
#                       Define the global variables                        #
#                                                                          #
############################################################################


# build global hash tables to save global variables
getGlobalVar <- function() {
  orgId.hash <<- hash(
    'hg19' = 9606,
    'hg38' = 9606,
    'mm10' = 10090,
    'mm9' = 10090
  ) # the correspondence between organisms and ids
  
  orgPathway.hash <<- hash(
    'hg19' = "hsa",
    'hg38' = "hsa",
    'mm10' = "mmu", 
    'mm9' = "mmu"
  ) # the correspondence between organisms and ids
  
  orgAnno.hash <<- hash(
    'hg19' = 'EnsDb.Hsapiens.v86',
    'hg38' = 'EnsDb.Hsapiens.v86',
    'mm10' = 'EnsDb.Mmusculus.v75',
    'mm9' = 'EnsDb.Mmusculus.v75'
  ) # the correspondence between organisms and annotations
  
  orgGS.hash <<- hash(
    'hg19' = 'BSgenome.Hsapiens.UCSC.hg19',
    'hg38' = 'BSgenome.Hsapiens.UCSC.hg38',
    'mm10'= 'BSgenome.Mmusculus.UCSC.mm10', 
    'mm9'= 'BSgenome.Mmusculus.UCSC.mm9'
  ) # the correspondence between organisms and genome sequences
  orgDB.hash <<- hash(
    'hg19' = 'org.Hs.eg.db',
    'hg38' = 'org.Hs.eg.db',
    'mm10' = 'org.Mm.eg.db', 
    'mm9' = 'org.Mm.eg.db'
  ) # the name database for an organism
} # '<<-' means generate a global variable


# get the name list of all chromosomes
getChrlist <- function(org) {
  
  # org : the name of organism
  
  chrlist <- vector() # an empty vector
  if (grepl(pattern = 'mm', org)) # if it's a mouse genome
  {
    chrlist <- c(as.character(1:19), 'X', 'Y') # mouse chromosomes
  }
  else
  {
    chrlist <- c(as.character(1:23), 'X', 'Y') # human chromosomes
  }
  chrlist <- paste0('chr', chrlist) # add the common prefix 'chr'
  return(chrlist)
}
