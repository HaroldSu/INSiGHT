# INSiGHT
INSiGHT (INtegrating multi-Samples for robust Spatial Gene Heterogeneity deTection) is a non-parametric statistical framework designed to robustly analyze spatial variability across multiple tissue slices.
![The overview of INSiGHT](figures/INSiGHT.jpg)
For each gene, INSiGHT ranks expression levels within samples, combines the rank scores, and constructs a gene expression similarity matrix, which is compared with a spatial similarity matrix derived from distance-based kernels. A rank-based $U$-statistic is then used to test the association between expression and spatial similarity, with inference based on its asymptotic distribution. Because it relies only on relative spatial distances, INSiGHT is rotation-invariant and ensures consistent results

## Installation
Please run the following codes in R to install STANCE package from GitHub.
```
if (!require("devtools", quietly = TRUE)){
  install.packages("devtools")
}
devtools::install_github("HaroldSu/INSiGHT")
```

## Normalization before ranking

`processINSiGHT()` defaults to `normalization = "library-size"`: each spot is
divided by its total input count and multiplied by 10,000 before ranking.
Library totals are stored before intersecting or filtering genes. Quality
control still uses the original input, and stored expression values are not
overwritten by normalization.

```r
object <- createINSiGHTobject(data.input = counts_list, spatial.locs = positions_list)
object <- processINSiGHT(object, normalization = "library-size")
# Other options: "log-library-size" applies log1p after normalization;
# "none" preserves the previous behavior and supports already-normalized input.
```

Optionally supply `size.factors` as a list of numeric vectors in sample order.
Named vectors are aligned by spot barcode; unnamed vectors must match the
current spot order. For matrices containing only a subset of genes, supply
factors computed from the full counts. Retained spots require positive factors.
The log transformation preserves normalized expression ordering, so the two
normalization modes give equivalent ranks apart from numerical precision.
