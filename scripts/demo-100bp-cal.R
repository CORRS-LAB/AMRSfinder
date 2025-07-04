library(dmrff)
library(AMRs.finder)
## define covariate matrix
## Load Demo data
intput_dat <- readRDS("data/bulk.sub.txt.20.Rds")
## generate y
y <- data.frame(y = rpois(20, 50))

source('demo-100bp-generation.R')


cortest<-function(intput_dat,y,cov.mod,a,b){
  x<-as.numeric(colMeans(intput_dat[a:b,-c(1:3)]))
  y <-as.numeric(y[,1])
  if (!is.null(cov.mod)){
       lm.dat <- data.frame(y,x,cov.mod)
  }else{
       lm.dat <- data.frame(y,x)
  }         
  fit <- summary(lm(y~.,data = lm.dat))
  p_value<-fit$coef[2,4]
  coef_lm<-fit$coef[2,1]
  cor_est<-cor(y,x) 
  ret_ks<-c(p_value,coef_lm, cor_est)
  return(ret_ks)
}


ewas.regions <- function(intput_dat, y, cov.mod){
  out <- NULL
  for(i in 1:nrow(intput_dat)){
    out <- rbind(out, cortest(intput_dat,y, cov.mod, i,i))
  }
  colnames(out) <- c("p_value", "coef_lm", "cor_est")
  res <- data.frame(intput_dat[,1:3], out)
  res$FDR <- p.adjust(res$p_value, method = "BH")
  return(res)
}

####
intput_dat <-read.table("bulk.sub.100bpWindows.20.txt", header = T)

nfo <- ewas.regions(intput_dat, y, cov.mod)
