## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(tidy=FALSE, collapse = TRUE, comment = "#>", message=FALSE, error=FALSE, warning=TRUE)

## ----Gene expression data------------------------------------------------
library(DIGREsyn)
head(geneExp.demo)

## ----Dose response data--------------------------------------------------
head(doseRes.demo)

## ----Gene connectivity data----------------------------------------------
head(geneNetLymph)

## ----DIGRE profile gene expression---------------------------------------
geneExpDiff <- profileGeneExp(geneExp = geneExp.demo)

## ----DIGRE gene network--------------------------------------------------
geneNetLymph.mat <- constGeneNet(geneNet = geneNetLymph)

## ----DIGRE prediction----------------------------------------------------
pred.res <- DIGREscore(geneExpDiff = geneExpDiff, doseRes = doseRes.demo, pathway = "GeneNet", geneNet = geneNetLymph.mat, fold = 0.6)

## ----Visualization heatmap, fig.align='center'---------------------------
vis.heat <- DIGREvis(pred.pair = pred.res$scoreRank, type = "heat")
plot(vis.heat)

## ----Visualization barplot, fig.align='center'---------------------------
vis.bar <- DIGREvis(pred.pair = pred.res$scoreRank, type = "bar")
plot(vis.bar)

## ----sessionInfo---------------------------------------------------------
sessionInfo()

