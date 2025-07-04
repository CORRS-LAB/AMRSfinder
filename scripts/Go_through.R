## Go through


## Level 1: Naive Case
intput_dat = readRDS("data/bulk.sub.txt.20.Rds")
y <- data.frame(y=rpois(20, 50))
cov.mod <- NULL 

library(AMRs.finder)
nfo <- AMR.finder(intput_dat, y, cov.mod)
head(nfo)


