library(scran)
library(Rtsne)
dir.create("results", showWarnings=FALSE)

# Normalizing both batches to the same scale.

sce.68 <- readRDS("sce.pbmc68k.rds")
sce.4 <- readRDS("sce.t4k.rds")
universe <- intersect(rownames(sce.68), rownames(sce.4))

normed <- multiBatchNorm(sce.68[universe,], sce.4[universe,])
sce.68 <- normed[[1]]
sce.4 <- normed[[2]]

# Taking the top 50000 genes with the largest biological components.
dec.68 <- readRDS("dec.pbmc68k.rds")
dec.4 <- readRDS("dec.t4k.rds")
combined <- combineVar(dec.68[universe,], dec.4[universe,])
to.use <- rownames(combined)[order(combined$bio, decreasing=TRUE)[1:5000]]

sce.68 <- sce.68[to.use,]
sce.4 <- sce.4[to.use,]

# Organizing colors.
batch.id <- rep(1:2, c(ncol(sce.68), ncol(sce.4)))
batchcolor <- c("lavender", "lightcoral")[batch.id]

cd3e <- c(logcounts(sce.68)["ENSG00000198851",], logcounts(sce.4)["ENSG00000198851",])
cd3e <- pmin(cd3e, 2) # to preserve dynamic range.
allcolors <- viridis::viridis(20)[cut(cd3e, 20)]

# Making a plotting function.

plotFUN <- function(fname, Y, subset=NULL, ..., xlab="tSNE 1",ylab="tSNE 2",main="") {
    if (is.null(subset)) {
        subset <- seq_len(nrow(Y))
    }
    png(fname,width=900,height=700)
    par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
    plot(Y[,1], Y[,2], cex=2, pch=21, bg=allcolors[subset], xlab=xlab, ylab=ylab, main=main) 
    dev.off()
}

plotFUNb <- function(fname, Y, subset=NULL, ...) {
    if (is.null(subset)) {
        subset <- seq_len(nrow(Y))
    }
    png(fname,width=900,height=700)
    par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
    plot(Y[,1], Y[,2], cex=2, pch=21, bg=batchcolor[subset], ...)
    dev.off()
}

######################## 
# Running it without any correction (but with a PCA, to make things manageable w.r.t. time).

set.seed(1000)
out <- scran:::.multi_pca(list(logcounts(sce.68), logcounts(sce.4)), d=50, approximate=TRUE, use.crossprod=TRUE)
t.unc <- do.call(rbind, out)

# Generating a t-SNE plot.
set.seed(0)
tsne.unc <- Rtsne(t.unc, perplexity = 30)
plotFUN("results/tsne_unc_type.png", tsne.unc$Y, main="Uncorrected", xlab="tSNE 1",ylab="tSNE 2")
plotFUNb("results/tsne_unc_batch.png", tsne.unc$Y, main="Uncorrected", xlab="tSNE 1",ylab="tSNE 2")

rm(t.unc)
gc()

######################## 
# Performing the correction with faster MNN.

set.seed(1000)
mnn.out2 <- fastMNN(logcounts(sce.68), logcounts(sce.4), k=20, approximate=TRUE, cos.norm=TRUE)# BPPARAM=MulticoreParam(1))
t.mnn <- mnn.out2$corrected

# Generating a t-SNE plot.
set.seed(0)
tsne.mnn <- Rtsne(t.mnn, perplexity = 30)
plotFUN("results/tsne_mnn2_type.png", tsne.mnn$Y, main="Fast MNN", xlab="tSNE 1",ylab="tSNE 2")
plotFUNb("results/tsne_mnn2_batch.png", tsne.mnn$Y, main="Fast MNN", xlab="tSNE 1",ylab="tSNE 2")

rm(t.mnn)
gc()

