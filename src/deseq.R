myArgs <- commandArgs(trailingOnly = TRUE)


# install packages
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#BiocManager::install("DESeq2")
#BiocManager::install("apeglm")
#install.packages("corrplot")
#install.packages("VennDiagram")
#install.packages("pheatmap")

# ---------------- IMPORTS ---------------- #
library("DESeq2")

#library("corrplot")

#library("VennDiagram")

#library("pheatmap")

# GETTING VST, RANKING BY L2 NORM
norm_vec <- function(x) sqrt(sum(x^2))

gene_matrix <- as.matrix(read.csv(myArgs[1], row.names = 1, header = TRUE))
vsd <- vst(round(gene_matrix))
vsd_df <- data.frame(vsd)
vsd_df$sqerr <- apply(vsd_df, 1, norm_vec)
sorted_vsd_df = vsd_df[order(-vsd_df$sqerr), ]
final_vsd_cts = head(sorted_vsd_df, 8000)
final_vsd_cts$sqerr = NULL



############## ANCG CONTROL VS BPD #######################
# ---------------- DATA ---------------- #
print("Performing DESeq2 on AnCg control group vs BPD")
# get data
cts <- as.matrix(read.csv(myArgs[2], row.names = 1, header = TRUE))

coldata <- read.csv(myArgs[3], row.names=1)
coldata$disorder <- factor(coldata$disorder)

cts = cts[rownames(final_vsd_cts), ]

#create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = coldata,
                              design = ~ age + PMI + pH + disorder)

# factors in R!
dds$disorder <- factor(dds$disorder, levels = c("AnCgControl","AnCgBipolarDisorder"))

featureData <- data.frame(gene = rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)

dds = DESeq(dds, test = "LRT", reduced = ~ age + PMI + pH)

res_ancg_control_bpd = results(dds, alpha = 0.05)
res_ancg_control_bpd = res_ancg_control_bpd[res_ancg_control_bpd$baseMean >= 10,]
#hist(res_ancg_control_bpd$pvalue, ylim=c(0, 1500), col = 'red', main = 'AnCg Control vs BPD', xlab = 'P-value', breaks = 20)



############### ANCG CONTROL VS SZ ################################
# ---------------- DATA ---------------- #
print("Performing DESeq2 on AnCg control group vs SZ")
# get data
cts <- as.matrix(read.csv(myArgs[4], row.names = 1, header = TRUE))

coldata <- read.csv(myArgs[5], row.names=1)
coldata$disorder <- factor(coldata$disorder)

cts = cts[rownames(final_vsd_cts), ]


#create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = coldata,
                              design = ~ age + PMI + pH + disorder)

featureData <- data.frame(gene = rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)


# factors in R!
dds$disorder <- factor(dds$disorder, levels = c("AnCgControl", "AnCgSchizophrenia"))

dds <- DESeq(dds, test = "LRT", reduced = ~age + PMI + pH)

res_ancg_control_sz = results(dds, alpha = 0.05)
res_ancg_control_sz = res_ancg_control_sz[res_ancg_control_sz$baseMean >= 10,]
#hist(res_ancg_control_sz$pvalue, ylim = c(0,1500), col = 'red', main = 'AnCg Control vs SZ', xlab = 'P-value', breaks = 20)


############## ANCG CONTROL VS MDD ################
print("Performing DESeq2 on AnCg control group vs MDD")
#----------- DATA -------------------
cts <- as.matrix(read.csv(myArgs[6], row.names = 1, header = TRUE))

coldata <- read.csv(myArgs[7], row.names=1)
coldata$disorder <- factor(coldata$disorder)

cts = cts[rownames(final_vsd_cts), ]

#create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = coldata,
                              design = ~ age + PMI + pH + disorder)


featureData <- data.frame(gene = rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)

# factors in R!
dds$disorder <- factor(dds$disorder, levels = c("AnCgControl","AnCgMajorDepression"))

dds <- DESeq(dds, test = "LRT", reduced = ~ age + PMI + pH)

res_ancg_control_mdd = results(dds, alpha = 0.05)
res_ancg_control_mdd = res_ancg_control_mdd[res_ancg_control_mdd$baseMean >= 10,]
#hist(res_ancg_control_mdd$pvalue, ylim = c(0,1500), col = 'red', main = 'AnCg Control vs MDD', xlab = 'P-value', breaks = 20)


############## DLPFC CONTROL VS BPD #################
print("Performing DESeq2 on DLPFC control group vs BPD")
#----------- DATA -------------------
cts <- as.matrix(read.csv(myArgs[8], row.names = 1, header = TRUE))

coldata <- read.csv(myArgs[9], row.names=1)
coldata$disorder <- factor(coldata$disorder)

cts = cts[rownames(final_vsd_cts), ]

#create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = coldata,
                              design = ~ age + PMI + pH + disorder)


featureData <- data.frame(gene = rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)

# factors in R!
dds$disorder <- factor(dds$disorder, levels = c("DLPFCControl","DLPFCBipolarDisorder"))

dds <- DESeq(dds, test = "LRT", reduced = ~ age + PMI + pH)

res_dlpfc_control_bpd = results(dds, alpha = 0.05)
res_dlpfc_control_bpd = res_dlpfc_control_bpd[res_dlpfc_control_bpd$baseMean >= 10,]
#hist(res_dlpfc_control_bpd$pvalue, ylim = c(0,1500), col = 'blue', main = 'DLPFC Control vs BPD', xlab = 'P-value', breaks = 20)


############## DLPFC CONTROL VS SZ #################
print("Performing DESeq2 on DLPFC control group vs SZ")
#----------- DATA -------------------
cts <- as.matrix(read.csv(myArgs[10], row.names = 1, header = TRUE))

coldata <- read.csv(myArgs[11], row.names=1)
coldata$disorder <- factor(coldata$disorder)

cts = cts[rownames(final_vsd_cts), ]

#create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = coldata,
                              design = ~ age + PMI + pH + disorder)

featureData <- data.frame(gene = rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)

# factors in R!
dds$disorder <- factor(dds$disorder, levels = c("DLPFCControl","DLPFCSchizophrenia"))

dds <- DESeq(dds, test = "LRT", reduced = ~ age + PMI + pH)

res_dlpfc_control_sz = results(dds, alpha = 0.05)
res_dlpfc_control_sz = res_dlpfc_control_sz[res_dlpfc_control_sz$baseMean >= 10,]
#hist(res_dlpfc_control_sz$pvalue, ylim = c(0,1500), col = 'blue', main = 'DLPFC Control vs SZ', xlab = 'P-value', breaks = 20)


############## DLPFC CONTROL VS MDD #################
print("Performing DESeq2 on DLPFC control group vs MDD")
#----------- DATA -------------------
cts <- as.matrix(read.csv(myArgs[12], row.names = 1, header = TRUE))

coldata <- read.csv(myArgs[13], row.names=1)
coldata$disorder <- factor(coldata$disorder)

cts = cts[rownames(final_vsd_cts), ]

#create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = coldata,
                              design = ~ age + PMI + pH + disorder)

featureData <- data.frame(gene = rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)

# factors in R!
dds$disorder <- factor(dds$disorder, levels = c("DLPFCControl","DLPFCMajorDepression"))

dds <- DESeq(dds, test = "LRT", reduced = ~ age + PMI + pH)

res_dlpfc_control_mdd = results(dds, alpha = 0.05)
res_dlpfc_control_mdd = res_dlpfc_control_mdd[res_dlpfc_control_mdd$baseMean >= 10,]
#hist(res_dlpfc_control_mdd$pvalue, ylim = c(0,1500), col = 'blue', main = 'DLPFC Control vs MDD', xlab = 'P-value', breaks = 20)


############## NACC CONTROL VS BPD #################
print("Performing DESeq2 on nAcc control group vs BPD")
#----------- DATA -------------------
cts <- as.matrix(read.csv(myArgs[14], row.names = 1, header = TRUE))

coldata <- read.csv(myArgs[15], row.names=1)
coldata$disorder <- factor(coldata$disorder)

cts = cts[rownames(final_vsd_cts), ]

#create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = coldata,
                              design = ~ age + PMI + pH + disorder)


featureData <- data.frame(gene = rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)

# factors in R!
dds$disorder <- factor(dds$disorder, levels = c("nAccControl","nAccBipolarDisorder"))

dds <- DESeq(dds, test = "LRT", reduced = ~ age + PMI + pH)

res_nacc_control_bpd = results(dds, alpha = 0.05)
res_nacc_control_bpd = res_nacc_control_bpd[res_nacc_control_bpd$baseMean >= 10,]
#hist(res_nacc_control_bpd$pvalue, ylim = c(0,1500), col = 'green', main = 'nAcc Control vs BPD', xlab = 'P-value', breaks = 20)


############## NACC CONTROL VS SZ #################
print("Performing DESeq2 on nAcc control group vs SZ")
#----------- DATA -------------------
cts <- as.matrix(read.csv(myArgs[16], row.names = 1, header = TRUE))

coldata <- read.csv(myArgs[17], row.names=1)
coldata$disorder <- factor(coldata$disorder)

cts = cts[rownames(final_vsd_cts), ]


#create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = coldata,
                              design = ~ age + PMI + pH + disorder)

featureData <- data.frame(gene = rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)

# factors in R!
dds$disorder <- factor(dds$disorder, levels = c("nAccControl","nAccSchizophrenia"))

dds <- DESeq(dds, test = "LRT", reduced = ~ age + PMI + pH)

res_nacc_control_sz = results(dds, alpha = 0.05)
res_nacc_control_sz = res_nacc_control_sz[res_nacc_control_sz$baseMean >= 10,]
#hist(res_nacc_control_sz$pvalue, ylim = c(0,1500), col = 'green', main = 'nAcc Control vs SZ', xlab = 'P-value', breaks = 20)


############## NACC CONTROL VS MDD #################
print("Performing DESeq2 on nAcc control group vs MDD")
#----------- DATA -------------------
cts <- as.matrix(read.csv(myArgs[18], row.names = 1, header = TRUE))

coldata <- read.csv(myArgs[19], row.names=1)
coldata$disorder <- factor(coldata$disorder)

cts = cts[rownames(final_vsd_cts), ]


#create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = coldata,
                              design = ~ age + PMI + pH + disorder)

featureData <- data.frame(gene = rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)

# factors in R!
dds$disorder <- factor(dds$disorder, levels = c("nAccControl","nAccMajorDepression"))

dds <- DESeq(dds, test = "LRT", reduced = ~ age + PMI + pH)

res_nacc_control_mdd = results(dds, alpha = 0.05)
res_nacc_control_mdd = res_nacc_control_mdd[res_nacc_control_mdd$baseMean >= 10,]
#hist(res_nacc_control_mdd$pvalue, ylim = c(0,1500), col = 'green', main = 'nAcc Control vs MDD', xlab = 'P-value', breaks = 20)



############# 3x3 for pvalues ###################
jpeg(paste(myArgs[20], "histplot.jpg", sep = "/"), width = 1080, height = 720)
par(mfrow = c(3,3), mai = c(0.7, 0.7, 0.2, 0.1))
hist(res_ancg_control_sz$pvalue, ylim = c(0,1500), col = 'red', main = 'AnCg', xlab = NULL, ylab = expression(bold("SZ")), breaks = 20, cex.main = 2.5, cex.lab = 2, cex.axis = 1.5)
hist(res_dlpfc_control_sz$pvalue, ylim = c(0,1500), col = 'blue', main = 'DLPFC', xlab = NULL, ylab = "Frequency", breaks = 20, cex.main = 2.5, cex.lab = 2, cex.axis = 1.5)
hist(res_nacc_control_sz$pvalue, ylim = c(0,1500), col = 'green', main = 'nAcc', xlab = NULL, ylab = "Frequency", breaks = 20, cex.main = 2.5, cex.lab = 2, cex.axis = 1.5)
hist(res_ancg_control_bpd$pvalue, ylim=c(0, 1500), col = 'red', main = NULL, xlab = NULL, ylab = expression(bold("BPD")), breaks = 20, cex.main = 2.5, cex.lab = 2, cex.axis = 1.5)
hist(res_dlpfc_control_bpd$pvalue, ylim = c(0,1500), col = 'blue', main = NULL, xlab = NULL, ylab = "Frequency", breaks = 20, cex.main = 2.5, cex.lab = 2, cex.axis = 1.5)
hist(res_nacc_control_bpd$pvalue, ylim = c(0,1500), col = 'green', main = NULL, xlab = NULL, ylab = "Frequency", breaks = 20, cex.main = 2.5, cex.lab = 2, cex.axis = 1.5)
hist(res_ancg_control_mdd$pvalue, ylim = c(0,1500), col = 'red', main = NULL, xlab = expression(bold('P-value')), ylab = expression(bold("MDD")), breaks = 20, cex.main = 2.5, cex.lab = 2, cex.axis = 1.5)
hist(res_dlpfc_control_mdd$pvalue, ylim = c(0,1500), col = 'blue', main = NULL, xlab = expression(bold('P-value')), ylab = "Frequency", breaks = 20, cex.main = 2.5, cex.lab = 2, cex.axis = 1.5)
hist(res_nacc_control_mdd$pvalue, ylim = c(0,1500), col = 'green', main = NULL, xlab = expression(bold('P-value')), ylab = "Frequency", breaks = 20, cex.main = 2.5, cex.lab = 2, cex.axis = 1.5)
dev.off()



# ############# 3x3 for spearman #####################
# res_ancg_control_bpd_LFC = res_ancg_control_bpd$log2FoldChange
# res_ancg_control_sz_LFC = res_ancg_control_sz$log2FoldChange
# res_ancg_control_mdd_LFC = res_ancg_control_mdd$log2FoldChange
# ancg_df = data.frame(res_ancg_control_sz_LFC, res_ancg_control_bpd_LFC, res_ancg_control_mdd_LFC)

# res_dlpfc_control_bpd_LFC = res_dlpfc_control_bpd$log2FoldChange
# res_dlpfc_control_sz_LFC = res_dlpfc_control_sz$log2FoldChange
# res_dlpfc_control_mdd_LFC = res_dlpfc_control_mdd$log2FoldChange
# dlpfc_df = data.frame(res_dlpfc_control_sz_LFC, res_dlpfc_control_bpd_LFC, res_dlpfc_control_mdd_LFC)

# res_nacc_control_bpd_LFC = res_nacc_control_bpd$log2FoldChange
# res_nacc_control_sz_LFC = res_nacc_control_sz$log2FoldChange
# res_nacc_control_mdd_LFC = res_nacc_control_mdd$log2FoldChange
# nacc_df = data.frame(res_nacc_control_sz_LFC, res_nacc_control_bpd_LFC, res_nacc_control_mdd_LFC)

# all = data.frame(ancg_df, dlpfc_df, nacc_df)
# all_cor = cor(all, method = "spearman", use = "complete.obs")
# colnames(all_cor) <- c("SZ", "BPD", "MDD", "SZ", "BPD", "MDD", "SZ", "BPD", "MDD")
# rownames(all_cor) <- c("SZ", "BPD", "MDD", "SZ", "BPD", "MDD", "SZ", "BPD", "MDD")

# jpeg(paste(myArgs[20], "corrplot.jpg", sep = "/"), width = 1080, height = 720)
# corrplot(all_cor, method = "circle", diag = FALSE, tl.col = "black", tl.cex = 2, mar = c(4,4,4,4))
# mtext(text = c(expression(bold("nAcc")), expression(bold("DLPFC")), expression(bold("AnCg"))), side = 2, line = -11, at = c(2,5,8), las = 0, cex = 3)
# mtext(text = c(expression(bold("AnCg")), expression(bold("DLPFC")), expression(bold("nAcc"))), side = 1, line = -43, at = c(2,5,8), las = 0, cex = 3)
# dev.off()


# ################ Venn diagram #######################
# ancg_bpd_diffs = rownames(res_ancg_control_bpd[res_ancg_control_bpd$pvalue < 0.05,])
# ancg_sz_diffs = rownames(res_ancg_control_sz[res_ancg_control_sz$pvalue < 0.05,])
# ancg_mdd_diffs = rownames(res_ancg_control_mdd[res_ancg_control_mdd$pvalue < 0.05,])

# jpeg(paste(myArgs[20], "vennDiagram.jpg", sep = "/"), width = 1080, height = 1080)
# VD = venn.diagram(x = list(SZ = ancg_sz_diffs, BPD = ancg_bpd_diffs, MDD = ancg_mdd_diffs), filename = NULL, scaled = FALSE,
#                   fill = c("red", "blue", "green"), cex = 6, cat.cex = 5, cat.fontface = "bold", output = TRUE)
# dev.off()


# ################### HEAT MAP ####################
# # get data
# gene_matrix <- as.matrix(read.csv(myArgs[1], row.names = 1, header = TRUE))
# vsd <- vst(round(gene_matrix))
# vsd_df <- data.frame(vsd)
# vsd_df$sqerr <- apply(vsd_df, 1, norm_vec)
# sorted_vsd_df = vsd_df[order(-vsd_df$sqerr), ]
# final_vsd_cts = head(sorted_vsd_df, 8000)
# final_vsd_cts$sqerr = NULL

# cts <- as.matrix(read.csv(myArgs[4], row.names = 1, header = TRUE))

# coldata <- read.csv(myArgs[5], row.names=1)
# coldata$disorder <- factor(coldata$disorder)

# cts = cts[rownames(final_vsd_cts), ]

# #create DESeqDataSet object
# dds <- DESeqDataSetFromMatrix(countData = round(cts),
#                               colData = coldata,
#                               design = ~ age + PMI + pH + disorder)

# featureData <- data.frame(gene = rownames(cts))
# mcols(dds) <- DataFrame(mcols(dds), featureData)



# # factors in R!
# dds$disorder <- factor(dds$disorder, levels = c("AnCgControl", "AnCgSchizophrenia"))

# dds <- DESeq(dds, test = "LRT", reduced = ~age + PMI + pH)
# lfc <- lfcShrink(dds, coef = "disorder_AnCgSchizophrenia_vs_AnCgControl", type = "apeglm")

# res_ancg_control_sz = lfc
# res_ancg_control_sz = res_ancg_control_sz[res_ancg_control_sz$baseMean >= 10,]
# #hist(res_ancg_control_sz$pvalue, ylim = c(0,1500), col = 'red', main = 'AnCg Control vs SZ', xlab = 'P-value', breaks = 20)

# col_df <- as.data.frame((colData(dds))[,c("disorder", "age")])
# col_df$age <- NULL
# ann_colors = list(disorder = c("AnCgSchizophrenia" = "red", "AnCgControl" = "black"))
# res_hm = head(res_ancg_control_sz[res_ancg_control_sz$pvalue < 0.05,], 1003)
# filter_rn = rownames(res_hm)
# final_cts_hm = cts[filter_rn,]

# dev.off()
# jpeg(paste(myArgs[20], "heatMap.jpg", sep = "/"), width = 1080, height = 720)
# pheatmap(final_cts_hm * res_hm$log2FoldChange, scale = "row", 
#          cluster_rows = FALSE, cluster_cols=TRUE, show_rownames=FALSE, show_colnames = FALSE, 
#          annotation_col = col_df, annotation_colors = ann_colors, fontsize = 15)
# dev.off()
