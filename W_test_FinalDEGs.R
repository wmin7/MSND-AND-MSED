#nonparametric Wilcoxon test on the obtained tumor-related genes by thresholds

rm(list = ls())
setwd("E:/gse81861/ALL");
AllDEG <- read.csv("gene_HVG.csv", header = T, check.names = F)
FPKM <- read.csv("dupSymbol_FPKM.csv", header = T, check.names = F)

DEGData <- sapply(AllDEG$x, function(x){
  FPKM[which(FPKM[,1] == as.vector(x)), ]
})
DEGData <- as.data.frame(t(DEGData))
DEGData <- DEGData[,-1]

NormalSymbol <- read.csv("celltypeNM.csv", header = T, check.names = F)
TumorSymbol <- read.csv("AllTumorCellName.csv", header = T, check.names = F)

NormalData <- sapply(NormalSymbol$sample, function(x){
  DEGData[, which(as.vector(colnames(DEGData)) == as.vector(x))]
})
NormalData <- as.data.frame(NormalData)

TumorData <- sapply(TumorSymbol$x, function(x){
  DEGData[, which(as.vector(colnames(DEGData)) == as.vector(x))]
})
TumorData <- as.data.frame(TumorData)



##wilcoxon_ test
data <- cbind(NormalData, TumorData)
mydata <- data.matrix(data) 
mydata <- log2(mydata +1)


p_values <- apply(mydata, 1, function(x){
  wilcox.test(x[1:ncol(NormalData)], x[(ncol(NormalData)+1):ncol(mydata)])$p.value
})
# adjust p values
adj_p <- p.adjust(p_values,"BH")

pValues <- cbind(AllDEG, cbind(p_values, adj_p))
colnames(pValues) <- c("symbol", "p.value", "adj.p")
pValues <- pValues[order(pValues$adj.p),]
FinalDEG <- pValues[pValues$adj.p < 0.05, ]

Marker <- read.table("Marker.txt", header = T, check.names = F)
commonGene <- intersect(FinalDEG$symbol, Marker$Marker)
View(commonGene)

#output final tumor-related genes
write.csv(FinalDEG, "FinalDEG.csv", row.names = F)

