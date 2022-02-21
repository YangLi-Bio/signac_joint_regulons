############################################################################
#                                                                          #
#                           Load the R codes                               #
#                                                                          #
############################################################################


start_time <- Sys.time() # get the start time
message ("----> Loading libraries ...\n")
source("/fs/ess/PCON0022/liyang/STREAM/benchmarking/Signac_joint_regulons/library.R") # libraries
source("/fs/ess/PCON0022/liyang/STREAM/benchmarking/Signac_joint_regulons/globalVar.R") # global variables
source("/fs/ess/PCON0022/liyang/STREAM/benchmarking/Signac_joint_regulons/functions.R") # define functions
source("/fs/ess/PCON0022/liyang/STREAM/benchmarking/Signac_joint_regulons/cisAnalysis.R") # link peaks to genes
getGlobalVar() # build global variables
set.seed(1234) # set the seed for random computing


############################################################################
#                                                                          #
#                           Load the dataset                               #
#                                                                          #
############################################################################

# 
# Run : Rscript Signac_joint.R /fs/ess/PCON0022/liyang/STREAM/benchmarking/Signac_joint_regulons/example.RDS hg38 ./
# 
args <- commandArgs(T) # get parameters from Shell command
rds.path <- args[1] # the path to RDS file, e.g, "/fs/ess/PCON0022/liyang/STREAM/benchmarking/10x_human_brain_3k.RDS"
org <- args[2] # the organism, e.g., "hg38"
out.dir <- args[3] # the ourput directory


out.dir.array <- strsplit(x = rds.path, split = "[/|\\.]") %>% `[[` (1) # split the input file
out.dir <- paste0(out.dir, out.dir.array[length(out.dir.array) - 1], "_dir/")
dir.create(path = out.dir)
# create the output directory


org.pathway <- orgPathway.hash[[org]] # organism ID in KEGG
org.gs <- orgGS.hash[[org]] # genome sequences
org.anno <- orgAnno.hash[[org]] # the annotation
org.db <- orgDB.hash[[org]] # database


# convert characters into objects
specLib(org, org.anno, org.gs, org.db) # load specific libraries
org.gs <- get(org.gs)
org.anno <- get(org.anno)
org.db <- get(org.db)


pbmc <- readRDS(file = rds.path) # read the RDS file
message ("----> Finished loading the dataset composed of ", ncol(pbmc), " cells.\n")
DefaultAssay(pbmc) <- "ATAC"


############################################################################
#                                                                          #
#                          Link peaks to genes                             #
#                                                                          #
############################################################################

message ("----> Began running Signac to build peak-to-gene linkages ...\n")
pbmc <- RegionStats(pbmc, genome = org.gs) # first compute the GC content for each peak

# link peaks to genes
pbmc <- LinkPeaks2(
  object = pbmc,
  peak.assay = "ATAC",
  expression.assay = "RNA",
)

qs::qsave(Links(pbmc), paste0(out.dir, "/pbmc_links.qsave"))


# Annotate peaks using TF binding sites in JASPAR
message ("----> Annotating peaks using JASPAR TF binding sites ...\n")
hg.mus.jaspar.list <- qs::qread('/fs/ess/PCON0022/liyang/STREAM/databases/hg_mus_JASPAR_TF_binding_sites.qsave')
# load the JASPAR TF binding sites

double.hash <- add_TF(peaks = da_peaks$feature, hg.mus.jaspar.list = hg.mus.jaspar.list, 
                      org = org) # add TF binding sites
peak.TFs <- double.hash$peak.hash # TFs binding each peak
TF.peaks <- double.hash$TF.hash # peaks bound by each TF
rm(hg.mus.jaspar.list)


# Identify regulons
message ("----> Began identifying regulons ...\n")
links <- Links(pbmc) # get peak-gene linkages

regulons <- pbmclapply(names(TF.peaks), function(tf) {
  peaks <- TF.peaks[[tf]] # peaks bound by the TF
  genes <- unique(links$gene[links$peak %in% peaks]) # genes linked to the peaks
  
  return(list(TF = tf, genes = genes, peaks = intersect(peaks, unique(links$peaks)))) # return the regulon
}, mc.cores = min(ddetectCores(), length(TF.peaks)))

message ("----> Finished identifying ", length(regulons), " regulons.\n")


end_time <- Sys.time()
run_time <- end_time - start_time
message ("----> Running time: ", run_time, " min.\n")