library(dmrff)
library(AMRs.finder)
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

names(dmrff_region)[c(4:6)] <- c("#CpGs", "estimate_dmrff", "p_value_dmrff")
