# This script prepares data for the haematopoiesis analysis.

library(BiocFileCache)
bfc <- BiocFileCache(ask=FALSE)    

##########################################
##########################################

# Read the counts, metadata of Nestorowa et al. 2016
library(scater)
fname.F <- bfcrpath(bfc, "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81682/suppl/GSE81682_HTSeq_counts.txt.gz")
dataF <- readSparseCounts(fname.F)

library(SingleCellExperiment)
sceF <- SingleCellExperiment(list(counts=dataF))
sceF <- sceF[!grepl("^_", rownames(sceF)),]
sceF <- sceF[!grepl("^ERCC-", rownames(sceF)),] # removing spike-ins completely.
dim(sceF)

# Loading in the metadata.
meta.fname <- bfcrpath(bfc, "http://blood.stemcells.cam.ac.uk/data/all_cell_types.txt")
metaF <- read.table(meta.fname, stringsAsFactors = FALSE, header=TRUE, check.names=FALSE)
metainds <- match(colnames(dataF), rownames(metaF))
metaF <- metaF[metainds,] # This will contain NA's... which is okay, at this point, to preserve length.

# Defining the cell type based on the metadata.
metatypeF <- rep("other", nrow(metaF))
for (col in rev(colnames(metaF))) { # reverse, so earlier columns end up overwriting later ones.
    chosen <- metaF[,col]==1
    metatypeF[chosen] <- sub("[0-9]?_.*", "", col)
}
metatypeF[metatypeF=="ESLAM"] <- "HSPC"

# Filling in metadata from the cell sorting label, if metadata was missing.
missing.meta <- is.na(metainds)
metatypeF[missing.meta] <- sub("_.*", "", colnames(dataF)[missing.meta])
metatypeF[metatypeF=="LT-HSC"] <- "LTHSC"
metatypeF[metatypeF=="Prog"] <- "other"
sceF$CellType <- metatypeF

# Perform size factor normalization within this data set.
library(scran)
clustF <- quickCluster(sceF, method="igraph", min.mean=1)
sceF <- computeSumFactors(sceF, cluster=clustF, min.mean=1)

# Cleaning up memory.
gc() 

##########################################
##########################################

# Download and read the counts and meta data of Paul et al. 2015
fname.A <- bfcrpath(bfc, "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE72nnn/GSE72857/suppl/GSE72857_umitab.txt.gz")
dataA <- readSparseCounts(fname.A, quote='"')
sceA <- SingleCellExperiment(list(counts=dataA))
dim(sceA)

# Only selecting cells that are in the metadata.
map.fname <- bfcrpath(bfc, "https://raw.githubusercontent.com/MarioniLab/MNN2017/master/Haematopoiesis/MAP.csv")
metaA <- read.csv2(map.fname, sep=",", stringsAsFactors = FALSE, head=TRUE, row.names=1)
metainds <- match(rownames(metaA), colnames(sceA))
sceA <- sceA[,metainds]

# Organizing cell type labels.
metatypeA <- character(nrow(metaA))
metatypeA[metaA[,1]<7] <- "ERY"
metatypeA[metaA[,1]>6 & metaA[,1]<12] <- "CMP"
metatypeA[metaA[,1]>11] <- "GMP"
sceA$CellType <- metatypeA

# Perform size factor normalization within this data set.
clustA <- quickCluster(sceA, method="igraph", min.mean=0.1)
sceA <- computeSumFactors(sceA, cluster=clustA, min.mean=0.1)

# Cleaning up memory.
gc() 

##########################################
##########################################

# Download list of highly variable genes identified by Nestrowa et al. 2016
hvg.fname <- bfcrpath(bfc, "http://blood.stemcells.cam.ac.uk/data/coordinates_gene_counts_flow_cytometry.txt.gz")
TFs <- read.table(hvg.fname, nrows=1, stringsAsFactors=FALSE)
features <- as.character(unlist(TFs))
features <- features[grep("ENSMUS", features)]

# Pull down IDs from BioMaRt.
library(biomaRt)
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host="www.ensembl.org" )
out <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"), values = features, mart = mart,filters = "ensembl_gene_id")

# Select features that are HVGs _and_ present in both data sets.
mF <- match(out$ensembl_gene_id, rownames(sceF))
mA <- pmatch(out$mgi_symbol, rownames(sceA)) # partial, due to use of concatenated gene symbols.
keep <- !is.na(mF) & !is.na(mA)

sceA <- sceA[mA[keep],]
sceF <- sceF[mF[keep],]
rownames(sceA) <- rownames(sceF)

# Save results to file.
saveRDS(sceA, file="haem_data_A.rds")
saveRDS(sceF, file="haem_data_F.rds")

###########
# END
