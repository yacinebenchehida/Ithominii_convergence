# Mapping, SNPs calling and phasing: 

This folder contains all the scripts used to perform all the step from the initial mapping of the raw reads to the phasing of the final VCF.

## 1) Mapping

``` bash
library("devtools")
devtools::install_github("yacinebenchehida/SAMPLE")
```

## Dependencies

-   R (\>= 4.3.0)

SAMPLES requires: `ggplot2`, `dplyr`, `Rmisc`,`RColorBrewer`, `magrittr`.

## Example usage
### Run the full pipeline
``` bash
library(SAMPLE)
data("coral_symbionts")
set.seed(812)
SAMPLE(input = coral_symbionts, output_N = "Example", replicates = 50, stability_thresh = 2, sucess_points = 10, diff = 1)
```
