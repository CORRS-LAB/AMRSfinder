## Load Demo data

```R
s = readRDS("data/bulk.sub.txt.20.Rds")
```

## generate y

```R
for(i in c(1:500)){
  y <- data.frame(y=rpois(20, 50))
  write.table(y, paste("y.0R.50.",i,sep=""), row.names = F, sep = "\t", quote = F)
}
```

##### run AMR.finder.R, ewas.R, ewas.100bpWindows.R separately

```sh
R CMD BATCH "--args $n 20" ewas.100bpWindows.R  > ewas.100bpWindows.R.out &
R CMD BATCH "--args $n 20" ewas.R  > ewas.R.out &
R CMD BATCH "--args $n 20" AMR.finder.R  > AMR.finder.R.out &
```
