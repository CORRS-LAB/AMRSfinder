cortest<-function(intput_dat,y,cov.mod,a,b){
  x<-as.numeric(colMeans(intput_dat[a:b,-c(1,2)]))
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


ewas <- function(intput_dat, y, cov.mod){
  out <- NULL
  for(i in 1:nrow(intput_dat)){
    out <- rbind(out, cortest(intput_dat,y, cov.mod, i,i))
  }
  colnames(out) <- c("p_value", "coef_lm", "cor_est")
  res <- data.frame(intput_dat[,1:2], out)
  res$FDR <- p.adjust(res$p_value, method = "BH")
  return(res)
}




args <- commandArgs(trailingOnly = TRUE)
reps <- as.numeric(args[1])
num <- args[2]

cov.mod <- NULL
dat <- read.table(paste("bulk.sub.txt.", num, sep = ""), header = T)
row_sd <- apply(dat[,-c(1:2)], 1, sd)
intput_dat <- dat[row_sd != 0, ]

  for(hh in c(1:50)){
  
       DifAMRs <- (reps-1)*50 + hh
       y <- read.table(paste(paste("y.0R.50.",DifAMRs,sep=""),num, sep = "."), header = T)
       colnames(y) <- "y"       
       nfo <- ewas(intput_dat, y, cov.mod)
       write.table(nfo,paste(paste("nfo.EWAS.out.",DifAMRs,sep=""),num, sep= "."), row.names = F, sep = "\t", quote = F)
 } 


