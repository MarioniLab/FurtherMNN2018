library(scran)
library(scater)
library(DropletUtils)

library(BiocFileCache)
bfc <- BiocFileCache(ask=FALSE)    

##########################################
##########################################

# Pre-processing the 68K PBMC dataset.
path.68 <- bfcrpath(bfc, "http://cf.10xgenomics.com/samples/cell-exp/1.1.0/fresh_68k_pbmc_donor_a/fresh_68k_pbmc_donor_a_filtered_gene_bc_matrices.tar.gz")
tmp.68 <- tempfile()
untar(path.68, exdir=tmp.68)
sce.68 <- read10xCounts(file.path(tmp.68, "filtered_matrices_mex/hg19/")) 

library(org.Hs.eg.db)
symb <- mapIds(org.Hs.eg.db, keys=rownames(sce.68), keytype="ENSEMBL", column="SYMBOL")
rowData(sce.68)$Symbol <- symb

# Adding locational annotation (using a slightly off-version ensembl, but chromosome assignment shouldn't change).
library(EnsDb.Hsapiens.v86)
loc <- mapIds(EnsDb.Hsapiens.v86, keys=rownames(sce.68), keytype="GENEID", column="SEQNAME")
rowData(sce.68)$Chr <- loc

# Brief quality control.
sce.68 <- calculateQCMetrics(sce.68, compact=TRUE, feature_controls=list(Mt=which(loc=="MT")))
lowlib <- isOutlier(sce.68$scater_qc$all$log10_total_counts, type="lower", nmads=3)
lowfeat <- isOutlier(sce.68$scater_qc$all$log10_total_features_by_counts, type="lower", nmads=3)
highmito <- isOutlier(sce.68$scater_qc$feature_control_Mt$pct_counts, type="higher", nmads=3)
discard <- lowlib | lowfeat | highmito
##summary(discard)
sce.68 <- sce.68[,!discard]

# Performing normalization, breaking the problem up into smaller blocks and subclustering within them.
blocks <- rep(seq_len(10), length.out=ncol(sce.68))
clusters <- quickCluster(sce.68, min.mean=0.1, block=blocks, method="igraph", block.BPPARAM=MulticoreParam(2))
##table(clusters)
##table(clusters, blocks)
sce.68 <- computeSumFactors(sce.68, clusters=clusters, min.mean=0.1, BPPARAM=MulticoreParam(2))
##plot(sce.68$scater_qc$all$total_counts, sizeFactors(sce.68), log="xy")

saveRDS(file="sce.pbmc68k.rds", sce.68)

# Modelling the mean-variance trend.
tmp <- normalize(sce.68) 
fit.68 <- trendVar(tmp, use.spikes=FALSE, loess.args=list(span=0.1))
##plot(fit.68$mean, fit.68$vars)
##curve(fit.68$trend(x), add=TRUE, col="red")
dec.68 <- decomposeVar(fit=fit.68)
saveRDS(file="dec.pbmc68k.rds", dec.68)

# Cleaning out the memory.
rm(list=ls())
gc()

##########################################
##########################################

# Pre-processing the 4K T-cell dataset.
path.4 <- bfcrpath(bfc, "http://cf.10xgenomics.com/samples/cell-exp/2.1.0/t_4k/t_4k_filtered_gene_bc_matrices.tar.gz")
tmp.4 <- tempfile()
untar(path.4, exdir=tmp.4)
sce.4 <- read10xCounts(file.path(tmp.4, "filtered_gene_bc_matrices/GRCh38/"))

symb <- mapIds(org.Hs.eg.db, keys=rownames(sce.4), keytype="ENSEMBL", column="SYMBOL")
rowData(sce.4)$Symbol <- symb

# Adding locational annotation.
library(EnsDb.Hsapiens.v86)
loc <- mapIds(EnsDb.Hsapiens.v86, keys=rownames(sce.4), keytype="GENEID", column="SEQNAME")
rowData(sce.4)$Chr <- loc

# Brief quality control.
sce.4 <- calculateQCMetrics(sce.4, compact=TRUE, feature_controls=list(Mt=which(loc=="MT")))
lowlib <- isOutlier(sce.4$scater_qc$all$log10_total_counts, type="lower", nmads=3)
lowfeat <- isOutlier(sce.4$scater_qc$all$log10_total_features_by_counts, type="lower", nmads=3)
highmito <- isOutlier(sce.4$scater_qc$feature_control_Mt$pct_counts, type="higher", nmads=3)
discard <- lowlib | lowfeat | highmito
##summary(discard)

sce.4 <- sce.4[,!discard]

# Performing normalization.
clusters <- quickCluster(sce.4, min.mean=0.1, method="igraph")
##table(clusters)
sce.4 <- computeSumFactors(sce.4, clusters=clusters, min.mean=0.1, BPPARAM=MulticoreParam(2))
##plot(sce.4$scater_qc$all$total_counts, sizeFactors(sce.4), log="xy")

saveRDS(file="sce.t4k.rds", sce.4)

# Modelling the mean-variance trend.
tmp <- normalize(sce.4)
fit.4 <- trendVar(tmp, use.spikes=FALSE, loess.args=list(span=0.1, control=loess.control(iterations=10)))
##plot(fit.4$mean, fit.4$vars)
##curve(fit.4$trend(x), add=TRUE, col="red")
dec.4 <- decomposeVar(fit=fit.4)
saveRDS(file="dec.pbmc4k.rds", dec.4)

