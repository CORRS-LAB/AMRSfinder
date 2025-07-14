library(AMRfinder)
## define covariate matrix
## Load Demo data
intput_dat <- readRDS("data/bulk.sub.txt.20.Rds")
## generate y
y <- data.frame(y = rpois(20, 50))
cov.mod <- NULL
cortest <- function(intput_dat, y, cov.mod, a, b) {
  x <- as.numeric(colMeans(intput_dat[a:b, -c(1, 2)]))
  y <- as.numeric(y[, 1])
  if (!is.null(cov.mod)) {
    lm.dat <- data.frame(y, x, cov.mod)
  } else {
    lm.dat <- data.frame(y, x)
  }
  fit <- summary(lm(y ~ ., data = lm.dat))
  p_value <- fit$coef[2, 4]
  coef_lm <- fit$coef[2, 1]
  se_lm <- fit$coef[2, 2]
  cor_est <- cor(y, x)
  ret_ks <- c(p_value, coef_lm, se_lm, cor_est)
  return(ret_ks)
}


ewas <- function(intput_dat, y, cov.mod) {
  out <- NULL
  for (i in 1:nrow(intput_dat)) {
    out <- rbind(out, cortest(intput_dat, y, cov.mod, i, i))
  }
  colnames(out) <- c("p_value", "coef_lm", "se_lm", "cor_est")
  res <- data.frame(intput_dat[, 1:2], out)
  res$FDR <- p.adjust(res$p_value, method = "BH")
  return(res)
}

nfo <- ewas(intput_dat, y, cov.mod)

# Expected results (leading six rows)
#     chr     pos   p_value    coef_lm     se_lm     cor_est       FDR
# 1 chr21 9437432 0.4236649  -7.161063  8.746876 -0.18947384 0.9893652
# 2 chr21 9437444 0.4116557 -11.764707 13.997314 -0.19433045 0.9893652
# 3 chr21 9437458 0.7138942  -3.599939  9.665112 -0.08745503 0.9893652
# 4 chr21 9437461 0.9052367  -1.058834  8.769744 -0.02844652 0.9893652
# 5 chr21 9437464 0.7004885  -3.118398  7.978302 -0.09173806 0.9893652
# 6 chr21 9437470 0.2603281 -15.754025 13.555210 -0.26420227 0.9893652