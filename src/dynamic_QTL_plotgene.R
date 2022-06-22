#Genes (MR positive): NAPSA, RALGDS, RAB2A, ADAM15, HIP1, IFNAR2, JD275616 (not found), ABO,
#    NAPSB (not found), OAS3, ENSG00000236263 (not found)

#get data
metainfo <- as.data.frame(fread("COVID19Donors_scRNA/lung.scp.metadata.txt"))
metainfo <- metainfo[-1,]
rownames(metainfo) <- metainfo$NAME

umap <- as.data.frame(fread("COVID19Donors_scRNA/upload.scp.X_umap.coords.txt"))
umap <- umap[-1,]
rownames(umap) <- umap$NAME
barcode <- as.data.frame(fread("COVID19Donors_scRNA/upload.scp.barcodes.tsv",header = F))
gene <- as.data.frame(fread("COVID19Donors_scRNA/upload.scp.features2.tsv",header = F))
exp <- as.data.frame(fread("COVID19Donors_scRNA/upload.scp.matrix.mtx.gz"))

#start gene specific work
gene.name = 'OAS3'
selectexp <- exp[exp$`%`==which(gene$V1==gene.name),]
rownames(selectexp) <- barcode$V1[selectexp$V2]

pltdat <- metainfo
pltdat$gene_count <- 0
pltdat[rownames(selectexp),]$gene_count <- selectexp$V3
pltdat <- pltdat[rownames(umap),]
pltdat$umap_x <- as.numeric(umap$X)
pltdat$umap_y <- as.numeric(umap$Y)

#plot showing umap for all tissues, not gene specific
ggplot(pltdat, aes(x = umap_x, y = umap_y, color = Cluster)) +
  geom_point() +
  theme_classic() +
  xlab("UMAP dimension 1") +
  ylab("UMAP dimension 2") +
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 14)) +
  labs(color = "") +
  scale_color_manual(values = c(rainbow(12)))


pltdat$predicted_celltype[pltdat$predicted_celltype==""] <- "unknown"
pltdat$predicted_celltype <- factor(pltdat$predicted_celltype, levels = sort(unique(pltdat$predicted_celltype))[c(1:27,29,28)])
pltdat$Viral <- ifelse(pltdat$`Viral+`,"Viral infection +","Viral infection -")
ggplot(pltdat, aes(x = predicted_celltype, y = gene_count, color = predicted_celltype)) +
  geom_boxplot() +
  facet_grid(Viral~.) +
  theme_classic() +
  stat_summary(fun=mean, geom="point",show.legend = FALSE, pch = '*', size = 9) +
  theme(legend.position = "none",
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 14),
        axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.5)) +
  xlab("") +
  labs(color = "Predicted cell type") +
  ylab("log(1 + TP10K)") +
  scale_color_manual(values = c(rainbow(28),"grey"))

################## statistical tests for enrichment
# Fisher-exact test was performed for viral+ cells, due to the relatively small total sample size
# Chi-square test was performed for viral- cells - this is very likely to yield a significant p-value due to the very large sample size, thus the fold of enrichment (odds ratio) is also important
# Odds ratio: (# selected cells expressing NPNT / # selected cells not expressing NPNT) / (# other cells expressing NPNT / # other cells not expressing NPNT)
# Wilcoxon signed-rank test (a non-parametric version of t-test) was performed to test whether the expression level in selected cells was higher than that in other cells
fishertest <- function(cell, label, expression) {
  exp_cell = sum(label == cell & expression != 0)
  noexp_cell = sum(label == cell & expression == 0)
  exp_other = sum(label != cell & expression != 0)
  noexp_other = sum(label != cell & expression == 0)
  results = fisher.test(matrix(c(exp_cell, noexp_cell, exp_other, noexp_other),2,2))
  cat("p-value",results$p.value,collapse="\n")
  cat("% expressing gene in selected cells",exp_cell / (exp_cell + noexp_cell) * 100,collapse="\n")
  cat("% expressing gene in other cells",exp_other / (exp_other + noexp_other) * 100,collapse="\n")
  cat("odds ratio",exp_cell / noexp_cell / (exp_other / noexp_other),collapse="\n")
}
chisqtest <- function(cell, label, expression) {
  exp_cell = sum(label == cell & expression != 0)
  noexp_cell = sum(label == cell & expression == 0)
  exp_other = sum(label != cell & expression != 0)
  noexp_other = sum(label != cell & expression == 0)
  results = chisq.test(matrix(c(exp_cell, noexp_cell, exp_other, noexp_other),2,2))
  cat("p-value",results$p.value,collapse="\n")
  cat("% expressing gene in selected cells",exp_cell / (exp_cell + noexp_cell) * 100,collapse="\n")
  cat("% expressing gene in other cells",exp_other / (exp_other + noexp_other) * 100,collapse="\n")
  cat("odds ratio",exp_cell / noexp_cell / (exp_other / noexp_other),collapse="\n")
}
wilcoxtest <- function(cell, label, expression) {
  exp_cell = expression[label == cell]
  exp_other = expression[label != cell]
  results = wilcox.test(exp_cell, exp_other, alternative = "greater") # note the alternative hypothesis is one-sided since we do not really care which cell types are depleted
  cat("p-value",results$p.value,collapse="\n")
  cat("quantiles (0,0.25,0.5,0.75,1) in selected cells",quantile(exp_cell),collapse="\n")
  cat("quantiles (0,0.25,0.5,0.75,1) in in other cells",quantile(exp_other),collapse="\n")
}
​
fishertest("AT1", pltdat$predicted_celltype[pltdat$Viral=="Viral infection +"], pltdat$gene_count[pltdat$Viral=="Viral infection +"])
chisqtest("AT1", pltdat$predicted_celltype[pltdat$Viral=="Viral infection -"], pltdat$gene_count[pltdat$Viral=="Viral infection -"])
wilcoxtest("AT1", pltdat$predicted_celltype[pltdat$Viral=="Viral infection +"], pltdat$gene_count[pltdat$Viral=="Viral infection +"])
wilcoxtest("AT1", pltdat$predicted_celltype[pltdat$Viral=="Viral infection -"], pltdat$gene_count[pltdat$Viral=="Viral infection -"])
​
fishertest("AT2", pltdat$predicted_celltype[pltdat$Viral=="Viral infection +"], pltdat$gene_count[pltdat$Viral=="Viral infection +"])
chisqtest("AT2", pltdat$predicted_celltype[pltdat$Viral=="Viral infection -"], pltdat$gene_count[pltdat$Viral=="Viral infection -"])
wilcoxtest("AT2", pltdat$predicted_celltype[pltdat$Viral=="Viral infection +"], pltdat$gene_count[pltdat$Viral=="Viral infection +"])
wilcoxtest("AT2", pltdat$predicted_celltype[pltdat$Viral=="Viral infection -"], pltdat$gene_count[pltdat$Viral=="Viral infection -"])

fishertest("fibroblast", pltdat$predicted_celltype[pltdat$Viral=="Viral infection +"], pltdat$gene_count[pltdat$Viral=="Viral infection +"])
chisqtest("fibroblast", pltdat$predicted_celltype[pltdat$Viral=="Viral infection -"], pltdat$gene_count[pltdat$Viral=="Viral infection -"])
wilcoxtest("fibroblast", pltdat$predicted_celltype[pltdat$Viral=="Viral infection +"], pltdat$gene_count[pltdat$Viral=="Viral infection +"])
wilcoxtest("fibroblast", pltdat$predicted_celltype[pltdat$Viral=="Viral infection -"], pltdat$gene_count[pltdat$Viral=="Viral infection -"])

fishertest("unknown", pltdat$predicted_celltype[pltdat$Viral=="Viral infection +"], pltdat$gene_count[pltdat$Viral=="Viral infection +"])
chisqtest("unknown", pltdat$predicted_celltype[pltdat$Viral=="Viral infection -"], pltdat$gene_count[pltdat$Viral=="Viral infection -"])
wilcoxtest("unknown", pltdat$predicted_celltype[pltdat$Viral=="Viral infection +"], pltdat$gene_count[pltdat$Viral=="Viral infection +"])
wilcoxtest("unknown", pltdat$predicted_celltype[pltdat$Viral=="Viral infection -"], pltdat$gene_count[pltdat$Viral=="Viral infection -"])
##################

pltdat <- pltdat[order(pltdat$gene_count),]
ggplot(pltdat, aes(x = umap_x, y = umap_y, color = gene_count)) +
  geom_point() +
  theme_classic() +
  scale_color_gradient(low = "lightgrey", high = "darkred") +
  xlab("UMAP dimension 1") +
  ylab("UMAP dimension 2") +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 14)) +
  labs(color = "log(1 + TP10K)")


pltdat <- pltdat[order(pltdat$`Viral+`),]
pltdat$`Viral+` <- ifelse(pltdat$`Viral+`,"Positive","Negative")
ggplot(pltdat, aes(x = umap_x, y = umap_y, color = `Viral+`)) +
  geom_point() +
  theme_classic() +
  xlab("UMAP dimension 1") +
  ylab("UMAP dimension 2") +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 14)) +
  scale_color_manual(values = c("darkblue","pink")) +
  labs(color = "Viral infection")
