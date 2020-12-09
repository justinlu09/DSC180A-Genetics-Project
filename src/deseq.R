# install packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

# ---------------- IMPORTS ---------------- #

# import packages
library("DESeq2")

install.packages("corrplot")
library("corrplot")

install.packages("VennDiagram")
library("VennDiagram")

############## ANCG CONTROL VS BPD #######################
# ---------------- DATA ---------------- #

# get data
cts <- as.matrix(read.csv('~/R/dsc180a_genetics/ancg_control_bpd.csv', row.names = 1, header = TRUE))
head(cts)

coldata <- read.csv('~/R/dsc180a_genetics/ancg_control_bpd_coldata.csv', row.names=1)
coldata$disorder <- factor(coldata$disorder)
head(coldata)

norm_vec <- function(x) sqrt(sum(x^2))

vsd <- vst(round(cts))
vsd_df <- data.frame(vsd)
vsd_df$sqerr <- apply(vsd_df, 1, norm_vec)
sorted_vsd_df = vsd_df[order(-vsd_df$sqerr), ]
final_vsd_cts = head(sorted_vsd_df, 8000)
final_vsd_cts$sqerr = NULL

head(final_vsd_cts,5)

cts = cts[rownames(final_vsd_cts), ]
all(rownames(coldata) == colnames(cts))

#create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = coldata,
                              design = ~ age + PMI + pH + disorder)

# factors in R!
dds$disorder <- factor(dds$disorder, levels = c("AnCgControl","AnCgBipolarDisorder"))

featureData <- data.frame(gene = rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)

dds = DESeq(dds, test = "LRT", reduced = ~ age + PMI + pH)


res_ancg_control_bpd = results(dds, alpha = 0.05)
res_ancg_control_bpd = res_ancg_control_bpd[res_ancg_control_bpd$baseMean >= 10,]
res_ancg_control_bpd
hist(res_ancg_control_bpd$pvalue, ylim=c(0, 1500), col = 'red', main = 'AnCg Control vs BPD', xlab = 'P-value', breaks = 20)




############### ANCG CONTROL VS SZ ################################
# ---------------- DATA ---------------- #

# get data
cts <- as.matrix(read.csv('~/R/dsc180a_genetics/ancg_control_sz.csv', row.names = 1, header = TRUE))
head(cts)

coldata <- read.csv('~/R/dsc180a_genetics/ancg_control_sz_coldata.csv', row.names=1)
coldata$disorder <- factor(coldata$disorder)
head(coldata)

vsd <- vst(round(cts))
vsd_df <- data.frame(vsd)
vsd_df$sqerr <- apply(vsd_df, 1, norm_vec)
sorted_vsd_df = vsd_df[order(-vsd_df$sqerr), ]
final_vsd_cts = head(sorted_vsd_df, 8000)
final_vsd_cts$sqerr = NULL
cts = cts[rownames(final_vsd_cts), ]

all(rownames(coldata) == colnames(cts))


#create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = coldata,
                              design = ~ age + PMI + pH + disorder)
dds

featureData <- data.frame(gene = rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)


# factors in R!
dds$disorder <- factor(dds$disorder, levels = c("AnCgControl", "AnCgSchizophrenia"))

dds <- DESeq(dds, test = "LRT", reduced = ~age + PMI + pH)

res_ancg_control_sz = results(dds, alpha = 0.05)
res_ancg_control_sz = res_ancg_control_sz[res_ancg_control_sz$baseMean >= 10,]
res_ancg_control_sz
hist(res_ancg_control_sz$pvalue, ylim = c(0,1500), col = 'red', main = 'AnCg Control vs SZ', xlab = 'P-value', breaks = 20)


############## ANCG CONTROL VS MDD ################

#----------- DATA -------------------
cts <- as.matrix(read.csv('~/R/dsc180a_genetics/ancg_control_mdd.csv', row.names = 1, header = TRUE))
head(cts)

coldata <- read.csv('~/R/dsc180a_genetics/ancg_control_mdd_coldata.csv', row.names=1)
coldata$disorder <- factor(coldata$disorder)
head(coldata)

vsd <- vst(round(cts))
vsd_df <- data.frame(vsd)
vsd_df$sqerr <- apply(vsd_df, 1, norm_vec)
sorted_vsd_df = vsd_df[order(-vsd_df$sqerr), ]
final_vsd_cts = head(sorted_vsd_df, 8000)
final_vsd_cts$sqerr = NULL
cts = cts[rownames(final_vsd_cts), ]

all(rownames(coldata) == colnames(cts))

#create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = coldata,
                              design = ~ age + PMI + pH + disorder)


featureData <- data.frame(gene = rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)

# factors in R!
dds$disorder <- factor(dds$disorder, levels = c("AnCgControl","AnCgMajorDepression"))

dds <- DESeq(dds, test = "LRT", reduced = ~ age + PMI + pH)

res_ancg_control_mdd <- results(dds, alpha = 0.05)
res_ancg_control_mdd = res_ancg_control_mdd[res_ancg_control_mdd$baseMean >= 10,]
hist(res_ancg_control_mdd$pvalue, ylim = c(0,1500), col = 'red', main = 'AnCg Control vs MDD', xlab = 'P-value', breaks = 20)



############## DLPFC CONTROL VS BPD #################

#----------- DATA -------------------
cts <- as.matrix(read.csv('~/R/dsc180a_genetics/dlpfc_control_bpd.csv', row.names = 1, header = TRUE))
head(cts)

coldata <- read.csv('~/R/dsc180a_genetics/dlpfc_control_bpd_coldata.csv', row.names=1)
coldata$disorder <- factor(coldata$disorder)
head(coldata)

vsd <- vst(round(cts))
vsd_df <- data.frame(vsd)
vsd_df$sqerr <- apply(vsd_df, 1, norm_vec)
sorted_vsd_df = vsd_df[order(-vsd_df$sqerr), ]
final_vsd_cts = head(sorted_vsd_df, 8000)
final_vsd_cts$sqerr = NULL
cts = cts[rownames(final_vsd_cts), ]

all(rownames(coldata) == colnames(cts))

#create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = coldata,
                              design = ~ age + PMI + pH + disorder)


featureData <- data.frame(gene = rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)

# factors in R!
dds$disorder <- factor(dds$disorder, levels = c("DLPFCControl","DLPFCBipolarDisorder"))

dds <- DESeq(dds, test = "LRT", reduced = ~ age + PMI + pH)

res_dlpfc_control_bpd <- results(dds, alpha = 0.05)
res_dlpfc_control_bpd = res_dlpfc_control_bpd[res_dlpfc_control_bpd$baseMean >= 10,]
hist(res_dlpfc_control_bpd$pvalue, ylim = c(0,1500), col = 'blue', main = 'DLPFC Control vs BPD', xlab = 'P-value', breaks = 20)


############## DLPFC CONTROL VS SZ #################

#----------- DATA -------------------
cts <- as.matrix(read.csv('~/R/dsc180a_genetics/dlpfc_control_sz.csv', row.names = 1, header = TRUE))
head(cts)

coldata <- read.csv('~/R/dsc180a_genetics/dlpfc_control_sz_coldata.csv', row.names=1)
coldata$disorder <- factor(coldata$disorder)
head(coldata)

vsd <- vst(round(cts))
vsd_df <- data.frame(vsd)
vsd_df$sqerr <- apply(vsd_df, 1, norm_vec)
sorted_vsd_df = vsd_df[order(-vsd_df$sqerr), ]
final_vsd_cts = head(sorted_vsd_df, 8000)
final_vsd_cts$sqerr = NULL
cts = cts[rownames(final_vsd_cts), ]

all(rownames(coldata) == colnames(cts))

#create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = coldata,
                              design = ~ age + PMI + pH + disorder)

featureData <- data.frame(gene = rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)

# factors in R!
dds$disorder <- factor(dds$disorder, levels = c("DLPFCControl","DLPFCSchizophrenia"))

dds <- DESeq(dds, test = "LRT", reduced = ~ age + PMI + pH)

res_dlpfc_control_sz <- results(dds, alpha = 0.05)
res_dlpfc_control_sz = res_dlpfc_control_sz[res_dlpfc_control_sz$baseMean >= 10,]
hist(res_dlpfc_control_sz$pvalue, ylim = c(0,1500), col = 'blue', main = 'DLPFC Control vs SZ', xlab = 'P-value', breaks = 20)


############## DLPFC CONTROL VS MDD #################

#----------- DATA -------------------
cts <- as.matrix(read.csv('~/R/dsc180a_genetics/dlpfc_control_mdd.csv', row.names = 1, header = TRUE))
head(cts)

coldata <- read.csv('~/R/dsc180a_genetics/dlpfc_control_mdd_coldata.csv', row.names=1)
coldata$disorder <- factor(coldata$disorder)
head(coldata)

vsd <- vst(round(cts))
vsd_df <- data.frame(vsd)
vsd_df$sqerr <- apply(vsd_df, 1, norm_vec)
sorted_vsd_df = vsd_df[order(-vsd_df$sqerr), ]
final_vsd_cts = head(sorted_vsd_df, 8000)
final_vsd_cts$sqerr = NULL
cts = cts[rownames(final_vsd_cts), ]

all(rownames(coldata) == colnames(cts))

#create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = coldata,
                              design = ~ age + PMI + pH + disorder)

featureData <- data.frame(gene = rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)

# factors in R!
dds$disorder <- factor(dds$disorder, levels = c("DLPFCControl","DLPFCMajorDepression"))

dds <- DESeq(dds, test = "LRT", reduced = ~ age + PMI + pH)

res_dlpfc_control_mdd <- results(dds, alpha = 0.05)
res_dlpfc_control_mdd = res_dlpfc_control_mdd[res_dlpfc_control_mdd$baseMean >= 10,]
hist(res_dlpfc_control_mdd$pvalue, ylim = c(0,1500), col = 'blue', main = 'DLPFC Control vs MDD', xlab = 'P-value', breaks = 20)


############## NACC CONTROL VS BPD #################

#----------- DATA -------------------
cts <- as.matrix(read.csv('~/R/dsc180a_genetics/nacc_control_bpd.csv', row.names = 1, header = TRUE))
head(cts)

coldata <- read.csv('~/R/dsc180a_genetics/nacc_control_bpd_coldata.csv', row.names=1)
coldata$disorder <- factor(coldata$disorder)
head(coldata)

vsd <- vst(round(cts))
vsd_df <- data.frame(vsd)
vsd_df$sqerr <- apply(vsd_df, 1, norm_vec)
sorted_vsd_df = vsd_df[order(-vsd_df$sqerr), ]
final_vsd_cts = head(sorted_vsd_df, 8000)
final_vsd_cts$sqerr = NULL
cts = cts[rownames(final_vsd_cts), ]

all(rownames(coldata) == colnames(cts))

#create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = coldata,
                              design = ~ age + PMI + pH + disorder)


featureData <- data.frame(gene = rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)

# factors in R!
dds$disorder <- factor(dds$disorder, levels = c("nAccControl","nAccBipolarDisorder"))

dds <- DESeq(dds, test = "LRT", reduced = ~ age + PMI + pH)

res_nacc_control_bpd <- results(dds, alpha = 0.05)
res_nacc_control_bpd = res_nacc_control_bpd[res_nacc_control_bpd$baseMean >= 10,]
hist(res_nacc_control_bpd$pvalue, ylim = c(0,1500), col = 'green', main = 'nAcc Control vs BPD', xlab = 'P-value', breaks = 20)


############## NACC CONTROL VS SZ #################

#----------- DATA -------------------
cts <- as.matrix(read.csv('~/R/dsc180a_genetics/nacc_control_sz.csv', row.names = 1, header = TRUE))
head(cts)

coldata <- read.csv('~/R/dsc180a_genetics/nacc_control_sz_coldata.csv', row.names=1)
coldata$disorder <- factor(coldata$disorder)
head(coldata)

vsd <- vst(round(cts))
vsd_df <- data.frame(vsd)
vsd_df$sqerr <- apply(vsd_df, 1, norm_vec)
sorted_vsd_df = vsd_df[order(-vsd_df$sqerr), ]
final_vsd_cts = head(sorted_vsd_df, 8000)
final_vsd_cts$sqerr = NULL
cts = cts[rownames(final_vsd_cts), ]

all(rownames(coldata) == colnames(cts))

#create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = coldata,
                              design = ~ age + PMI + pH + disorder)

featureData <- data.frame(gene = rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)

# factors in R!
dds$disorder <- factor(dds$disorder, levels = c("nAccControl","nAccSchizophrenia"))

dds <- DESeq(dds, test = "LRT", reduced = ~ age + PMI + pH)

res_nacc_control_sz <- results(dds, alpha = 0.05)
res_nacc_control_sz = res_nacc_control_sz[res_nacc_control_sz$baseMean >= 10,]
hist(res_nacc_control_sz$pvalue, ylim = c(0,1500), col = 'green', main = 'nAcc Control vs SZ', xlab = 'P-value', breaks = 20)


############## NACC CONTROL VS MDD #################

#----------- DATA -------------------
cts <- as.matrix(read.csv('~/R/dsc180a_genetics/nacc_control_mdd.csv', row.names = 1, header = TRUE))
head(cts)

coldata <- read.csv('~/R/dsc180a_genetics/nacc_control_mdd_coldata.csv', row.names=1)
coldata$disorder <- factor(coldata$disorder)
head(coldata)

vsd <- vst(round(cts))
vsd_df <- data.frame(vsd)
vsd_df$sqerr <- apply(vsd_df, 1, norm_vec)
sorted_vsd_df = vsd_df[order(-vsd_df$sqerr), ]
final_vsd_cts = head(sorted_vsd_df, 8000)
final_vsd_cts$sqerr = NULL
cts = cts[rownames(final_vsd_cts), ]


all(rownames(coldata) == colnames(cts))

#create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = coldata,
                              design = ~ age + PMI + pH + disorder)

featureData <- data.frame(gene = rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)

# factors in R!
dds$disorder <- factor(dds$disorder, levels = c("nAccControl","nAccMajorDepression"))

dds <- DESeq(dds, test = "LRT", reduced = ~ age + PMI + pH)

res_nacc_control_mdd <- results(dds, alpha = 0.05)
res_nacc_control_mdd = res_nacc_control_mdd[res_nacc_control_mdd$baseMean >= 10,]
hist(res_nacc_control_mdd$pvalue, ylim = c(0,1500), col = 'green', main = 'nAcc Control vs MDD', xlab = 'P-value', breaks = 20)


################# nAcc SPEARMAN CORRELATION PLOTTING ################

res_nacc_control_bpd_LFC = res_nacc_control_bpd$log2FoldChange
res_nacc_control_sz_LFC = res_nacc_control_sz$log2FoldChange
res_nacc_control_mdd_LFC = res_nacc_control_mdd$log2FoldChange
nacc_df = data.frame(res_nacc_control_sz_LFC, res_nacc_control_bpd_LFC, res_nacc_control_mdd_LFC)


cor_nacc = cor(nacc_df, method = "spearman", use = "complete.obs")
#cor_nacc = abs(cor_nacc)
corrplot(cor_nacc, method = "circle", diag = FALSE)




############# 3x3 for pvalues ###################
jpeg("~/R/dsc180a_genetics/histplot.jpg", width = 1080, height = 720)
par(mfrow = c(3,3), mai = c(0.7, 0.7, 0.2, 0.1))
hist(res_ancg_control_sz$pvalue, ylim = c(0,1500), col = 'red', main = 'AnCg', xlab = NULL, ylab = expression(bold("SZ")), breaks = 20, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
hist(res_dlpfc_control_sz$pvalue, ylim = c(0,1500), col = 'blue', main = 'DLPFC', xlab = NULL, ylab = "Frequency", breaks = 20, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
hist(res_nacc_control_sz$pvalue, ylim = c(0,1500), col = 'green', main = 'nAcc', xlab = NULL, ylab = "Frequency", breaks = 20, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
hist(res_ancg_control_bpd$pvalue, ylim=c(0, 1500), col = 'red', main = NULL, xlab = NULL, ylab = expression(bold("BPD")), breaks = 20, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
hist(res_dlpfc_control_bpd$pvalue, ylim = c(0,1500), col = 'blue', main = NULL, xlab = NULL, ylab = "Frequency", breaks = 20, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
hist(res_nacc_control_bpd$pvalue, ylim = c(0,1500), col = 'green', main = NULL, xlab = NULL, ylab = "Frequency", breaks = 20, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
hist(res_ancg_control_mdd$pvalue, ylim = c(0,1500), col = 'red', main = NULL, xlab = expression(bold('P-value')), ylab = expression(bold("MDD")), breaks = 20, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
hist(res_dlpfc_control_mdd$pvalue, ylim = c(0,1500), col = 'blue', main = NULL, xlab = expression(bold('P-value')), ylab = "Frequency", breaks = 20, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
hist(res_nacc_control_mdd$pvalue, ylim = c(0,1500), col = 'green', main = NULL, xlab = expression(bold('P-value')), ylab = "Frequency", breaks = 20, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
dev.off()



############# 3x3 for spearman #####################
jpeg("~/R/dsc180a_genetics/corrplot.jpg", width = 1080, height = 720)
all = data.frame(ancg_df, dlpfc_df, nacc_df)
all_cor = cor(all, method = "spearman", use = "complete.obs")
colnames(all_cor) <- c("SZ", "BPD", "MDD", "SZ", "BPD", "MDD", "SZ", "BPD", "MDD")
rownames(all_cor) <- c("SZ", "BPD", "MDD", "SZ", "BPD", "MDD", "SZ", "BPD", "MDD")
corrplot(all_cor, method = "circle", diag = FALSE, tl.col = "black", mar = c(5,5,5,5))
mtext(text = c(expression(bold("nAcc")), expression(bold("DLPFC")), expression(bold("AnCg"))), side = 2, line = -12, at = c(2,5,8), las = 0)
mtext(text = c(expression(bold("AnCg")), expression(bold("DLPFC")), expression(bold("nAcc"))), side = 1, line = -42, at = c(2,5,8), las = 0)
dev.off()


################ Venn diagram #######################
ancg_bpd_diffs = rownames(res_ancg_control_bpd[res_ancg_control_bpd$pvalue < 0.05,])
ancg_sz_diffs = rownames(res_ancg_control_sz[res_ancg_control_sz$pvalue < 0.05,])
ancg_mdd_diffs = rownames(res_ancg_control_mdd[res_ancg_control_mdd$pvalue < 0.05,])

jpeg("~/r/dsc180a_genetics/vennDiagram.jpg", width = 1080, height = 1080)
VD = venn.diagram(x = list(SZ = ancg_sz_diffs, BPD = ancg_bpd_diffs, MDD = ancg_mdd_diffs), filename = NULL, scaled = FALSE,
                  fill = c("red", "blue", "green"), cat.fontface = "bold", output = TRUE)
grid.draw(VD)
dev.off()
