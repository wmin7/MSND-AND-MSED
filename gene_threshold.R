#determine threshold combinations

rm(list = ls())
path <- "E:/gse81861/ALL/gene"  
fileNames <- dir(path)
filePath <- sapply(fileNames,function(x){
  paste(path, x, sep='/')
})  

threshold <- matrix(nrow = length(fileNames), ncol = 4)
colnames(threshold) <- c("group", "geneNumber", "insectNumber", "gene")

marker <- read.table("E:/gse81861/ALL/Marker.txt", header = F, check.names = F)
marker <- t(marker)[-1]

for(j in 1:length(fileNames)){
  genepath <- dir(filePath[j])
  genefile <- paste(filePath[j], genepath, sep = "/")
  data <- lapply(genefile,function(x){
    read.table(x, header = F, check.names = F)
  })  
  
  if(length(data) == 11){
    xve <- sapply(data, function(x){
      as.vector(as.matrix(x))
    })
    
    Epithelial <- intersect(xve[[1]], intersect(xve[[4]], intersect(xve[[7]], xve[[9]])))
    Fibroblast <- intersect(xve[[10]], intersect(xve[[2]], xve[[5]]))
    Immune <- intersect(xve[[11]], intersect(xve[[8]], intersect(xve[[6]], xve[[3]])))
    
    gene <- intersect(Epithelial, intersect(Fibroblast, Immune))
    
    commGene <- intersect(gene, marker)
    
    threshold[j, 1] <- fileNames[j]
    threshold[j, 2] <- length(gene)
    threshold[j, 3] <- length(commGene)
    
    if(length(commGene) == 1){
      threshold[j, 4] <- commGene
    }else{
      insectGENE <- c()
      for(m in 1:length(commGene)){
        insectGENE <- paste(insectGENE, commGene[m], sep = "/")
        threshold[j, 4] <- insectGENE
      }
    }
  }
  
}

write.csv(threshold, "E:/gse81861/ALL/threshold_select.csv", row.names = F)

