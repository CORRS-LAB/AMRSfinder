library(dmrff)
library(AMRfinder)
## define covariate matrix
## Load Demo data
intput_dat <- readRDS("data/bulk.sub.txt.20.Rds")
## generate y
y <- data.frame(y = rpois(20, 50))

dmrff_region_origin <- dmrff(
    estimate = nfo$coef_lm,
    se = nfo$se_lm,
    p.value = nfo$p_value,
    methylation = as.matrix(intput_dat[, -c(1, 2)]),
    chr = nfo$chr,
    pos = nfo$pos,
    maxgap = 300,
    verbose = T
)
dmrff_region_origin <- dmrff_region_origin[, c("chr", "start", "end", "n", "estimate", "se", "p.value")]


dmrff_region_origin$start.idx <- match(dmrff_region_origin$start, intput_dat$pos)
dmrff_region_origin$end.idx <- match(dmrff_region_origin$end, intput_dat$pos)

result <- t(apply(dmrff_region_origin, 1, function(row) {
    a <- as.integer(row["start.idx"])
    b <- as.integer(row["end.idx"])
    cortest(intput_dat, y, cov.mod, a, b)
}))

dmrff_region_origin[, c("p_value_lm", "coef_lm", "cor_est")] <- result[, c(1, 2, 4)]


dmrff_region <- dmrff_region_origin[, c(
    "chr", "start", "end", "n",
    "estimate", "p.value", "p_value_lm",
    "coef_lm", "cor_est"
)]

names(dmrff_region)[c(4:6)] <- c("N.CpGs", "estimate_dmrff", "p_value_dmrff")

# Expected results (leading six rows)
#     chr    start      end N.CpGs estimate_dmrff p_value_dmrff p_value_lm    coef_lm    cor_est
# 1 chr21 47581998 47581998     1     -656.99191  1.758439e-05  0.1772850 -22.719491 -0.3142063
# 2 chr21 47581500 47581744     2      -77.69146  3.437775e-04  0.2475903   8.093722  0.2711178
# 3 chr21 47580627 47580627     1     -146.12122  1.234963e-03  0.5533747   5.199087  0.1409451
# 4 chr21 47580575 47580575     1      -39.14796  3.634867e-02  0.2749670   5.719088  0.2565139
# 5 chr21 47580592 47580592     1       40.88840  5.975636e-01  0.2247231  12.227191  0.2841291
# 6 chr21 47581497 47581497     1      -40.80384  3.113511e-03  0.1086950   7.413681  0.3696475
