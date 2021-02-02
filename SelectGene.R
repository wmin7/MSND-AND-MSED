#select tumor-related genes for 25 combination thresholds

#install.packages("stabledist")
library(stabledist)
#install.packages("StableEstim")
library(StableEstim)

rm(list = ls())
path <- "E:/gse81861/ALL/screen"  
fileNames <- dir(path)
filePath <- sapply(fileNames,function(x){
  paste(path, x, sep='/')
})  
data <- lapply(filePath,function(x){
  read.table(x, header = T, check.names = F)
})  

NM <- read.table("E:/gse81861/ALL/FPKM_NM.txt", header = T, check.names = F)
para <- read.csv("E:/gse81861/ALL/allpar.csv", header = T, row.names = 1, check.names = F)
qat <- c(0.1, 0.15, 0.2, 0.25, 0.3)

for(m in 1:5){
  l <- m
  for(n in 1:5){
    r <- n
    for (i in 1:nrow(para)) {
      sh1 <- qstable(qat[l], para[i,1], para[i,2], para[i,3], para[i,4])
      sh2 <- qexp(qat[r], para[i,6])
      
      num <- apply(data[[i]], 1, function(x){
        l1 <- length(which(x < sh1))
        l2 <- length(which(x > sh2))
        l1 + l2
      })
      num <- as.data.frame(num)
      row.names(num) <- row.names(NM)
      
      HVG <- apply(num, 2, function(x){
        which(x > 0.8*ncol(data[[i]]))
      })
      HVG <- row.names(HVG)
      
      if (length(HVG) > 1) {
        path <- paste("E:/gse81861/ALL/gene/", "l", qat[l], "r", qat[r], sep = "")
        outfile <- paste(paste(path, row.names(para)[i], sep = "/"), "gene.txt", sep = "_")
        write.table(HVG, outfile, row.names = F, sep = "\t", col.names = F)
      }
      
    }
  }
}
