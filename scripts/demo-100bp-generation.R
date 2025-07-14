library(GenomicRanges)
library(IRanges)


dat <- readRDS("data/bulk.sub.txt.20.Rds")
gr <- GRanges(data.frame(
         chr = "chr21",
         start = dat[, 2],
         end = dat[, 2],
         dat[, -c(1, 2)]))

cutpoint <- seq(dat[1, 2], dat[nrow(dat), 2], by = 100)
start <- cutpoint
end <- c(cutpoint[-length(cutpoint)] + 99, dat[nrow(dat), 2])

wds <- GRanges(data.frame(chr = "chr21", start, end))

hits.df <- as.data.frame(findOverlaps(wds, gr))
revmap <- IntegerList(tapply(hits.df[, 2], hits.df[, 1], unique))

wds <- wds[as.numeric(revmap@partitioning@NAMES), ]
wds$revmap <- revmap

ids <- colnames(dat)[-c(1:2)]
sum.count <- NULL
for (i in seq_along(ids)) {
   tmp <- mean(extractList(dat[, i + 2], wds$revmap))
   sum.count <- cbind(sum.count, tmp)
}

colnames(sum.count) <- colnames(dat)[-c(1:2)]

wds.out <- data.frame(wds)
out <- data.frame(chr = wds.out[, 1], wds.out[, c(2, 3)], sum.count)
write.table(out, "bulk.sub.100bpWindows.20.txt", row.names = F, sep = "\t", quote = F)
