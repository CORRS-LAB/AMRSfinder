library(AMRfinder)

intput_dat <- readRDS("data/bulk.sub.txt.20.Rds")
y <- data.frame(y=rpois(20, 50))
## define covariate matrix
cov.mod <- NULL
AMRnfo <- AMRfinder(intput_dat, y, cov.mod)
AMRnfo$FDR <- p.adjust(AMRnfo$p_value, method = "BH")

head(AMRnfo)

# Expected results (leading six rows)
#     chr   start     end #CpGs    cor_est    coef_lm    p_value     methX methY       FDR
# 1 chr21 9699132 9699607     9 -0.1910182 -33.462156 0.41982590 0.4832250 48.45 0.5359946
# 2 chr21 9704369 9704508     9  0.5501849  57.124946 0.01195525 0.6864139 48.45 0.1432859
# 3 chr21 9708981 9709116    11  0.2116939  34.762865 0.37026316 0.6740095 48.45 0.4987533
# 4 chr21 9825538 9825579     7  0.1874799  18.927252 0.42864898 0.2264622 48.45 0.5420658
# 5 chr21 9825581 9825646     8  0.0998579   8.445006 0.67530903 0.3084667 48.45 0.7377364
# 6 chr21 9825651 9825690     8  0.1837555  10.722467 0.43804174 0.3885482 48.45 0.5517968