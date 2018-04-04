![DIGREsyn_logo](QBRC.jpg)

# DIGREsyn
Drug Induced Genomic Residual Effect (DIGRE) model for predicting compound pair synergy

## Introduction
The Drug-Induced Genomic Residual Effect (DIGRE) algorithm is developed to predict compound pair synergistic effect. The algorithm ranked first in the NCI DREAM challenge competition of predicting 91 compound pairs from the most synergistic to the most antagonistic. 

Cite our papers:  
["A community computational challenge to predict the activity of pairs of compounds."](http://www.nature.com/nbt/journal/v32/n12/full/nbt.3052.html) Nature biotechnology 32.12 (2014): 1213-1222.  
["DIGRE: Drug-Induced Genomic Residual Effect Model for Successful Prediction of Multidrug Effects."](http://onlinelibrary.wiley.com/doi/10.1002/psp4.1/abstract;jsessionid=2874AA1DB0B0E048BA1041B8D46BA07D.f03t04) CPT: pharmacometrics & systems pharmacology 4.2 (2015): 91-97.

## Before started
To use the DIGRE model to predict your drug pair synergism, we suggest you read the original DIGRE model paper and get famillar the experimental settings about how to generate corresponding experimental data that DIGRE can accept. The Vignette document will also help for preparing input data. If you just want an idea of how DIGRE works, you can run DIGRE model on the demo data we provided.


## Get started
### Install ###
DIGREsyn is still under development. You can install it from `github` using the `devtool`. `DIGREsyn` depends on other R packages. In some machines, it may take additional effort to get all the dependencies ready.

```{r}
library(devtools)
install_github("DIGREsyn", "Minzhe")
library(DIGREsyn)
```

### Prepare input data ###
The DIGRE model takes three forms of input data to predict the compound synergistic effect. Users will first need to read their drug treated gene expression data, drug dose response data and gene interaction data (optional) into R. A specific format is preferred. Load the example data file to see the format.

**1. Read and profile gene expression data**
Read the gene expression data of cells treated with each compound as well as the negative controls into R, and parse it with the `profileGeneExp` function. <u>*(Notice: Do not need to specify column names and row names by drug names and gene names for gene expression data as duplicated drug names and gene names are common.)*</u>
```{r}
data(geneExp.demo)
geneExpDiff <- profileGeneExp(geneExp = geneExp.demo)
```

**2. Read dose response data**  
Parse and plot the dose response curve of each compound by yourself (the functions are not provided in this package), and extract two important values from the curve for each drug, which DIGREsyn will need to do the prediction. The two values of each compound are: 1) cell viability reduction (*percentage of cells killed in a certain dose*) under a single dose (*compound concentration used for genomic profiling*) and 2) cell viability reduction under double dose. In the demo file, we use IC20 for the drug-treated gene expression data, so we look for viability reduction under IC20 and 2 \* IC20 concentrations. <u>*(Notice: Please check the drug name in the dose response data. It should exactly match the drug name in the gene expression data. Be careful when reading those data into R, R will check column names and will convert space to dot. You should disable it.)*</u>
```{r}
data(doseRes.demo)
```

**3. Read and parse gene interaction data (optional)**
Read the gene connectivity data to construct the gene network that DIGRE uses to compare compound effects on upstream and downstream genes. DIGRE uses the constructed KEGG pathway information by default, but *de novo* construction of a gene network from the user's own data is also supported. We provide here gene interaction data refined from lymphoma patients as demo files. <u>*(Notice: Be careful to set stringAsFactor = FALSE when reading data to prevent undefined behavior.)*
```{r}
data(geneNetLymph)
geneNetLymph.mat <- constGeneNet(geneNet = geneNetLymph)
```

### Predicting drug pair synergy ###
In the `DIGREsyn` package, `DIGREscore` is the core function that calculates all the possible compound pair synergistic scores and their ranks. <u>*(Notice: for the following, we are using a default cut off of 0.6 for gene expression difference, but you can also set your own preferred value.)*</u>
```{r}
res.KEGG <- DIGREscore(geneExpDiff = geneExpDiff, doseRes = doseRes.demo, pathway = "KEGG", fold = 0.6)
res.geneNet <- DIGREscore(geneExpDiff = geneExpDiff, doseRes = doseRes.demo, pathway = "GeneNet", geneNet = geneNetLymph.mat)
```

## Result visualization
The `DIGREvis` function is designed to visualize prediction results, so specify the parameter `type` as either `heat` or `bar`.
```{r}
vis.heat <- DIGREvis(pred.pair = res.KEGG$scoreRank, type = "heat")
vis.bar <- DIGREvis(pred.pair = res.KEGG$scoreRank, type = "bar")
plot(vis.heat)
plot(vis.bar)
```

## Version update
0.1.0: First release. (07-14-2016)  
0.2.0: Support de novo construction of gene network. (10-7-2016)  
0.2.1: Improve documentation. (04-04-2018)
