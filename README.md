<!-- README.md is generated from README.Rmd. Please edit that file -->
R package rankdist
==================

Implements distance based probability models for ranking data. The supported models include Mallow's Phi-model (Kendall distance), Weighted Kendall distance model and Phi-Component model. Mixture models are also supported.

Examples
========

Fit a single-cluster Mallow's Phi-model to APA election data.

``` r
library(rankdist)
testctrl <- new("RankControlKendall")
testinit <- new("RankInit",param.init=list(1),modal_ranking.init=list(c(2,3,4,1,5)),clu=1L)
ken_c1 <- RankDistanceModel(apa_obj,testinit,testctrl)
```

Fit a three-cluster Weighted Kendall distance model to APA election data.

``` r
testinit <- new("RankInit",param.init=list(rep(0.1,4),rep(0.1,4),rep(0.1,4)),modal_ranking.init=list(c(3,4,5,1,2),c(2,3,1,5,4),c(4,2,5,3,1)),clu=3L,p.init=rep(1,3)/3)
testctrl <- new("RankControlWeightedKendall")
# wken_c3 <- RankDistanceModel(apa_obj,testinit,testctrl)
```
