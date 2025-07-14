# AMRSfinder

AMRfinder (Associated Methylation Region finder) is a novel algorithm specifically designed to identify phenotype-associated methylation regions (especially for RRBS data). AMRfinder integrates circular binary segmentation with a correlation-based scoring function to dynamically detect methylated blocks most strongly associated with phenotype.


## Install

```sh
git pull https://github.com/CORRS-LAB/AMRfinder
cd AMRfinder
## Use Pure R version
R CMD INSTALL package-r-naive
## Use Rcpp version
R CMD INSTALL package
```

For advanced users, we recommend use the Rcpp version for better performance.

## Demo

With `AMRfinder::AMRfinder`, one could perform the following calculation:

```R
library(AMRfinder)
## Load Demo data
intput_dat = readRDS("data/bulk.sub.txt.20.Rds")
## generate y
y <- data.frame(y=rpois(20, 50))
## define covariate matrix
cov.mod <- NULL
nfo <- AMRfinder(intput_dat, y, cov.mod)
head(nfo)
```

```R

The Result is as following:

|chr   |   start|     end| #CpGs|    cor_est|    coef_lm|   p_value|     methX| methY|       FDR|
|:-----|-------:|-------:|-----:|----------:|----------:|---------:|---------:|-----:|---------:|
|chr21 | 9437431| 9437538|    14|  0.5832074|  38.328544| 0.0069528| 0.6017233|    49| 0.2759106|
|chr21 | 9825466| 9825568|    12|  0.1596132|  11.463980| 0.5014693| 0.2783314|    49| 0.6066556|
|chr21 | 9825571| 9825600|     5|  0.1243991|   9.382826| 0.6012908| 0.2705365|    49| 0.6540356|
|chr21 | 9825839| 9825864|     5| -0.2768194| -16.275038| 0.2373935| 0.3259030|    49| 0.4405461|
|chr21 | 9825870| 9825894|     6| -0.2385190| -13.892837| 0.3111941| 0.4334588|    49| 0.4909643|
|chr21 | 9826120| 9826207|    20| -0.0306179|  -2.452290| 0.8980380| 0.2583344|    49| 0.8980380|

The detailed description of the above output is in the following:

| Output | Description | 
|:-----|:--------|
| chr | Chromosome | 
| start | The start position of a methylation region. | 
| end | The start position of a methylation region. | 
| #CpGs | The number of CpG sites.|
| cor_est | The estimated correlation coefficient. | 
| coef_lm | The methylation coefficient derived from linear regression analysis. | 
| p_value |  The methylation P value derived from linear regression analysis. | 
| methX | Mean DNA methylation level within each region across all samples.|
| methY | Mean phenotype value within each region across all samples.|
| FDR | False discovery rate using BH approach.|


## Further examples to test ewas and others

| Demo | Link | Comments |
|:-----|:-----|:---------|
| AMRfinder | [AMRfinder](./scripts/AMR.finder.R) | Perform AMRfinder, and calculate the BH type FDR.|
| ewas | [ewas](./scripts/demo-ewas.R) |  Perform correlation and linear regression analyses to each CpG site.|
| dmrff | [dmrff](./scripts/demo-dmrff.R) | Dynamically identify methylation regions associated with phenotype using `dmrff`. Install [dmrff](https://github.com/perishky/dmrff) (referred to https://doi.org/10.1101/508556) first. Run `demo-ewas.R` first to initialize the objetcs.|
| 100bp-gen | [Data-generation](./scripts/demo-100bp-generation.R) |  Segment CpG sites into non-overlapping 100bp genomic regions.|
| 100bp-cal | [Calculation](./scripts/demo-100bp-cal.R) | Perform correlation and linear regression analyses to 100bp genomic regions. Run `emo-100bp-generation.R` first to generate data.|

```mermaid
graph TD;
  A["CpGs & Phenotype"] --> B[EWAS]
  A --> C[dmrff]
  A --> D["100bp (Data Generation)"]
  D --> E["100bp (Calculation)"]
  A --> F["AMRfinder"]
```
