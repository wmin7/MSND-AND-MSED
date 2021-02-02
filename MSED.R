#estimate parameter of MSED

#install.packages("stabledist")
library(stabledist)
#install.packages("StableEstim")
library(StableEstim)

rm(list = ls())
#set workspace
setwd("E:/gse81861/ALL/screen");

#import data
path <- getwd() 
fileNames <- dir(path)  
filePath <- sapply(fileNames, function(x){
  paste(path, x, sep='/')
}) 
data <- lapply(filePath, function(x){
  read.table(x, header = T, check.names = F)
})  

data <- lapply(data,function(x){
  as.vector(as.matrix(x))
})

#create null data frame to save parameter
allpar <- as.data.frame(matrix(nrow = length(data), ncol = 6))

for(i in 1:length(data)){
  subname <- substr(names(data)[i], 1, nchar(names(data)[i])-4)
  row.names(allpar)[i] <- subname
}
colnames(allpar) <- c("alpha", "beta", "gamma", "delta", "lambda", "eta")

#set t_points and u_points to estimate alpha, beta,
# gamma and delta
#It's key to obtain appropriate parameter.
t1 <- 0.1; t2 <- 2.5; tout = 20;
u1 <- 0.2; u2 <- 1.85; uout = 12;
Tpoints <- paste(paste(t1, t2, sep = ","), tout, sep = ",")
Upoints <- paste(paste(u1, u2, sep = ","), uout, sep = ",")
UTpoints <- paste(Tpoints, Upoints, sep = "&")

for(i in 1:length(data)){
  mpar <- KoutParametersEstim(x = data[[i]], pm = 1, spacing = "free",
                              t_points = seq(t1, t2, length.out = tout),
                              u_points = seq(u1, u2, length.out = uout),
                              PrintTime = TRUE)
  
  allpar[i, 1:4] <- round(mpar[["Estim"]][["par"]], 4)
}
#View(allpar)

## estimate lambda
for(m in 1:length(data)){
  p1 <- allpar[m, 1]
  p2 <- allpar[m, 2]
  p3 <- allpar[m, 3]
  p4 <- allpar[m, 4]
  if(p4 > 0){
    s0 <- p4
  }else{
    s0 <- 0
  }
  xxve <- data[[m]][data[[m]] > s0]
  
  ep <- fitdistr(xxve, densfun = "exponential") 
  allpar[m, 5] <- round(ep[["estimate"]][["rate"]], 4)
}

## tuning eta
#j is length(data), 
j <- 1 #j<-2(3,4,...,length(data))
p1 <- allpar[j, 1]
p2 <- allpar[j, 2]
p3 <- allpar[j, 3]
p4 <- allpar[j, 4]
p5 <- allpar[j, 5]

eta <- 0.41

opar <- par(no.readonly = T)
par(cex.axis = 1.5,cex.lab = 1.5,mfrow = c(1, 2),mar = c(5,6,3.5,2)-0.7)
plot(density(xve[[j]]), lty = 1, bty = "l", xlab = "", main = "", ylab = "", ylim = c(0,0.2),  col = "green")

ss <- seq(min(xve[[j]]), max(xve[[j]]), 0.01)
s1 <- ss[ss <= s0]
s2 <- ss[ss > s0]

den <- c(dstable(s1, p1, p2, p3, p4), eta*dexp(s2, p5))
lines(ss, den, lty = 1, col = "red")

plot(density(xve[[m]]), lty = 1, bty = "l", xlab = "", 
     main = "", ylab = "", col = "green")
lines(ss, den, lty = 2, col = "red")
par(opar)

allpar[j, 6] <- eta


#output parameter
outname <- paste(paste("E:/gse81861/ALL/parameter/p", UTpoints, sep = ""), ".csv", sep = "")
write.csv(allpar, outname, row.names = T)