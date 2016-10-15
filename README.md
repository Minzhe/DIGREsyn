![DIGREsyn_logo](QBRC.jpg)

# DIGREsyn
Drug Induced Genomic Residual Effect (DIGRE) model for predicting compound pair synergy

## Introduction
The Drug-Induced Genomic Residual Effect (DIGRE) algorithm is developed to predict the compound pair synergistic effect. It ranked first in the NCI DREAM challenge competition of predicting 91 compound pairs from the most synergistic to the most antagonist. 

Cite our paper:  
["A community computational challenge to predict the activity of pairs of compounds."](http://www.nature.com/nbt/journal/v32/n12/full/nbt.3052.html) Nature biotechnology 32.12 (2014): 1213-1222.  
["DIGRE: Drug-Induced Genomic Residual Effect Model for Successful Prediction of Multidrug Effects."](http://onlinelibrary.wiley.com/doi/10.1002/psp4.1/abstract;jsessionid=2874AA1DB0B0E048BA1041B8D46BA07D.f03t04) CPT: pharmacometrics & systems pharmacology 4.2 (2015): 91-97.

## Before started
If user want to use DIGRE model to predict your own compounds data, we suggest you read our DIGRE model paper first to get famillar with all experiment settings and generate corresponding experimental data that DIGRE could recognize. If you just want to have an idea about how DIGRE works, you could look at this document and try to run demo files.


## Get started
### Install ###
DIGREsyn is still under development. You can install it from `github` using the `devtool` package. `DIGREsyn` depend on other R package. In some machine, it may take some effort to get all the dependency ready.

```{r}
library(devtools)
install_github("DIGREsyn", "Minzhe")
library(DIGREsyn)
```

### Prepare input data ###
DIGRE model take three input data to predict compound synergistic effect. User will first need to read their drug treated gene expression data, drug dose response data and gene interaction data (optional) into R. Specific format is desired, load example data file to see the format.

**1. Read and profile gene expression data**
Read gene expression data of cells cheated by each compound as well as negative control into R, and parse it with `profileGeneExp` function. <u>*(Notice: Do not specify column names and row names by drug names and gene names for gene expression data as duplicated drug names and gene names are usual.)*</u>
```{r}
data(geneExp.demo)
geneExpDiff <- profileGeneExp(geneExp = geneExp.demo)
```

**2. Read dose response data**  
Parse and plot your dose response curve of each compound by yourself (we are not providing functions here), and extract two important value from the curve for each drug which DIGREsyn need to do prediction. The two value of each compound are 1) cell viability reduction(*percentage of cell killed in a certain dose*) under single dose(*compound concentration used for genomic profiling*) 2) cell viability reduction under double dose. In the demo file, we are using IC20 for the drug-treated gene expression data, so we look for viability reduction under IC20 and 2 \* IC20 concentration. <u>*(Notice: Please check the drug name in dose response data, it should exactly match the drug name in gene expression data. Be careful when read those data into R, R will check column names and will convert space to dot.)*</u>
```{r}
data(doseRes.demo)
```

**3. Read and parse gene interaction data (optional)**
Read gene connnectivity data to construct gene network that DIGRE use to compare compound effect on upstream and downstream genes. DIGRE use constructed KEGG pathway information by default, *de novo* construction of gene network from user's own data is also supported. We provided gene interaction data refined from lymphoma patients as demo files. <u>*(Notice: Be careful to set stringAsFactor = FALSE when reading data to prevent strange things happening)*
```{r}
data(geneNetLymph)
geneNetLymph.mat <- constGeneNet(geneNet = geneNetLymph)
```

### Scoring pair synergy ###
`DIGREscore` is the core function in `DIGREsyn` package to calculate all the possible compound pairs synergistic score and their ranks. <u>*(Notice: following we are using default cut off 0.6 for gene expression difference, we could also set your own preferable value.)*</u>
```{r}
res.KEGG <- DIGREscore(geneExpDiff = geneExpDiff, doseRes = doseRes.demo, pathway = "KEGG", fold = 0.6)
res.geneNet <- DIGREscore(geneExpDiff = geneExpDiff, doseRes = doseRes.demo, pathway = "GeneNet", geneNet = geneNetLymph.mat)
```

## Result visualization
`DIGREvis` function is to visualize prediction result, specify parameter `type` to be either `heat` or `bar`.
```{r}
vis.heat <- DIGREvis(pred.pair = res.KEGG$scoreRank, type = "heat")
vis.bar <- DIGREvis(pred.pair = res.KEGG$scoreRank, type = "bar")
print(vis.heat)
print(vis.bar)
```

## Version update
0.1.0: First release. (7-14-2016)  
0.2.0: Support de novo construction of gene network. (10-7-2016)