#' Drug Induced Genomic Residual Effect (DIGRE) model
#'
#' This function estimate the compound pair synergism/antagonism.
#' It takes drug treated gene expression data and drug dose response curve data as input,
#' and caculate the synergistic score using Drug Induced Genomic Residual Effect (DIGRE) model.
#'
#' @param geneExpDiff a matrix of drug treated gene expression data. Each column represent one drug,
#' each row represent one gene. The value represent the fold change of the expression of a particular
#' gene after drug treated compared to negative control.
#' @param doseRes a matrix contains drug dose response data. Each column represent one drug, two rows
#' of two different drug dose response curve. See `doseRes.demo` for example.
#' @param pathway pathway information used in DIGRE model. User would specify either "KEGG" to use the KEGG
#' pathway information or "GeneNet" to use the gene network information.
#' @param geneNet optional parameter. If pathway parameter is "GeneNet", then specify this parameter to use 
#' your own gene network data. If pathway parameter is "KEGG", then do not set this parameter.
#' @param fold a value between 0 and 1. Gene expression fold change above this value would be considered
#' as upregulated, below the opposite number would be considered as downregulated, otherwise would be
#' considered as no effect. The default value is 0.6.
#' @return a list contains two matrices. One gives the drug pair synergistic score and their rank, the other
#' contains the raw data to calculate the score.
#'
#' @author Jichen Yang, Sangin Lee, Minzhe Zhang(\email{zenroute.mzhang@gmail.com})
#' @details
#' This function takes drug treated gene expression data, dose-response curve data and pathway information as inputs,
#' and calculate pair synergistic score of all the possible combination of the compound you provided, and their rank
#' from the most synergistic to the most antagonist. Larger score indicates high possibility of the pair to have
#' synergistic effect, and vice versa. (Notice that this algorithm focus more on predicting the relative rank of your
#' compound pairs not the exact synergistic strength. If you want to do that, maybe you should involve positive control
#' in your experiment. And also the score calculated by two pathway information is not comparable.)
#' @seealso
#' Bansal M, Yang J, Karan C, et al. A community computational challenge to predict the activity of pairs of compounds[J].
#' Nature biotechnology, 2014, 32(12): 1213-1222.
#' @seealso
#' Yang J, Tang H, Li Y, et al. DIGRE: Drug Induced Genomic Residual Effect Model for Successful Prediction of Multidrug
#' Effects[J]. CPT: pharmacometrics & systems pharmacology, 2015, 4(2): 91-97.
#' 
#' @examples 
#' ### profile gene expression data
#' geneExpDiff <- profileGeneExp(geneExp = geneExp.demo)
#' ### DIGRE prediction
#' res.KEGG <- DIGREscore(geneExpDiff = geneExpDiff, doseRes = doseRes.demo, pathway = "KEGG", fold = 0.6) # KEGG pathway 
#' res.geneNet <- DIGREscore(geneExpDiff = geneExpDiff, doseRes = doseRes.demo, pathway = "GeneNet", geneNet = geneNetLymph.mat) # Gene network
#'
#' @export
#' @import org.Hs.eg.db KEGGgraph AnnotationDbi




DIGREscore <- function(geneExpDiff, doseRes, pathway = "KEGG", geneNet, fold = 0.60) {

      cat("Start scoring compound pairs by DIGRE model ...\n")
      cat("------------\n")
      
      ### load pathway matrix
      data("CGP.mat", "KEGGnet.mat")
      if (pathway == "KEGG") {      # using KEGG pathway
            if (missing(geneNet)) {
                  cat("Using KEGG pathway(default) information\n")
                  GP.mat <- KEGGnet.mat
            } else {
                  stop("Error: Conflict of gene network. Do not specify geneNet parameter if using KEGG pathway, or set pathway to be GeneNet to continue use your own gene network.")
            }
      } else if (pathway == "GeneNet") {  # using own gene network
            GP.mat <- geneNet
            cat("Using self-constructed gene network information\n")
      } else {
            stop('Pathway must be either "KEGG" or "GeneNet"\n')
      }
      
      cat("Fold change cut off used:", fold, "(default:0.6)\n")

      ### Translate gene name to KEGG id
      KEGG.id <- gene2KEGGid(geneExpDiff$Gene)
      idx <- gsub("hsa:", "", KEGG.id, fixed = TRUE) != "NA"
      geneExpDiff <- geneExpDiff[idx,-1]
      row.names(geneExpDiff) <- KEGG.id[idx]


      ### Estimate similarities between two drugs
      ## generate drug pairs
      drugName <- colnames(geneExpDiff)
      drugPair <- combn(drugName, 2)

      ## construct similarity score table
      sim.score.mat <- matrix(NA, nrow = ncol(drugPair), ncol = 10, dimnames = list(1:ncol(drugPair), c("drugA", "drugB", "Similarity", "mPositive", "mNegative", "mFalse", "iPositive", "iNegative", "iFalse", "Score")))
      sim.score.mat <- data.frame(sim.score.mat)

      ## loop through all drug pairs
      for (i in 1:ncol(drugPair)) {
            drugA <- drugPair[1,i]
            drugB <- drugPair[2,i]
            sim.score.mat[i, "drugA"] <- drugA
            sim.score.mat[i, "drugB"] <- drugB

            ## calculate (drugA, drugB) and (drugB, drugA) respectively
            sim.score.pair <- matrix(NA, nrow = 2, ncol = 10, dimnames = list(1:2, c("drugA", "drugB", "Similarity", "mPositive", "mNegative", "mFalse", "iPositive", "iNegative", "iFalse", "Score")))
            sim.score.pair <- data.frame(sim.score.pair)
            for (j in 1:2) {
                  if (j == 1) tempPair <- c(drugA, drugB)
                  if (j == 2) tempPair <- c(drugB, drugA)
                  sim.score.pair[j,"drugA"] <- tempPair[1]
                  sim.score.pair[j,"drugB"] <- tempPair[2]
                  geneExpDiff.pair <- geneExpDiff[,tempPair]

                  # compare to fold change cut off
                  geneExpDiff.pair[geneExpDiff.pair > fold] <- 1
                  geneExpDiff.pair[geneExpDiff.pair < -fold] <- -1
                  geneExpDiff.pair[abs(geneExpDiff.pair) <= fold] <- 0

                  # marginal relationship: effect on same gene in CGP
                  mSim.pair <- mSim.score(geneExpDiff.pair = geneExpDiff.pair, CGP.mat = CGP.mat)
                  sim.score.pair[j,"mPositive"] <- mSim.pair$mPos
                  sim.score.pair[j,"mNegative"] <- mSim.pair$mNeg
                  sim.score.pair[j,"mFalse"] <- mSim.pair$mNon

                  # interaction relationship: effect of upstream gene in GP
                  iSim.pair <- iSim.score(geneExpDiff.pair = geneExpDiff.pair, CGP.mat = CGP.mat, GP.mat = GP.mat)
                  sim.score.pair[j,"iPositive"] <- iSim.pair$iPos
                  sim.score.pair[j,"iNegative"] <- iSim.pair$iNeg
                  sim.score.pair[j,"iFalse"] <- iSim.pair$iNon

                  # calculate the ratio
                  denom <- iSim.pair$countB.updown
                  if (denom > 0) {
                        ratio <- (mSim.pair$mSim + iSim.pair$iSim) / denom
                  } else {ratio <- 0}
                  if (ratio == 0) {
                        ratio <- 0.01
                  }
                  # adjust ratio to avoid same ratio for different drug pair (use geneExp data before fold change curation)
                  ratio <- ratio * max(geneExpDiff[row.names(geneExpDiff)%in%row.names(CGP.mat),tempPair])

                  # calculate residual effect
                  fA <- doseRes[1,tempPair[1]]
                  fB <- doseRes[1,tempPair[2]]
                  f2B <- doseRes[2,tempPair[2]]
                  effect <- 1-(1-fA)*(1-ratio*f2B)*(1-(1-ratio)*fB)
                  effect <- effect-(1-(1-fA)*(1-fB))

                  sim.score.pair[j,"Similarity"] <- round(ratio, 5)
                  sim.score.pair[j,"Score"] <- round(effect, 5)
            }
            sim.score.mat[i,-(1:2)] <- colMeans(sim.score.pair[,-(1:2)])
      }
      pair.rank <- data.frame(sim.score.mat[,c("drugA", "drugB", "Score")], Rank = rank(-sim.score.mat$Score))
      cat("Done scoring.")
      return(list(scoreRank = pair.rank, rawTable = sim.score.mat))

} # function


### translate gene name to KEGG id
gene2KEGGid <- function(geneName) {

      entre.id <- sapply(mget(geneName, org.Hs.egSYMBOL2EG, ifnotfound = NA), "[[" , 1)
      kegg.id <- translateGeneID2KEGGID(entre.id, organism = "hsa")

      return(kegg.id)
}


### Marginal similarity: effect on same gene in CGP
mSim.score <- function(geneExpDiff.pair, CGP.mat) {

      # filter genes in CGP
      geneExpDiff.pair <- geneExpDiff.pair[row.names(geneExpDiff.pair) %in% row.names(CGP.mat),]

      # marginal relationship
      mPos <- sum(rowSums(geneExpDiff.pair) == 2)
      mNeg <- sum(rowSums(geneExpDiff.pair) == -2)
      mNon <- sum(apply(geneExpDiff.pair, 1, prod) == -1)
      mSim <- mPos + mNeg - mNon
      mPara <- list(mPos = mPos, mNeg = mNeg, mNon = mNon, mSim = mSim)

      return(mPara)
}


### Interaction relationship: effect of upstream gene in GP
iSim.score <- function(geneExpDiff.pair, CGP.mat, GP.mat) {

      # filter genes
      idx <- (row.names(geneExpDiff.pair) %in% row.names(CGP.mat)) & (row.names(geneExpDiff.pair) %in% row.names(GP.mat))
      geneExpDiff.pair <- geneExpDiff.pair[idx,]
      idx <- row.names(GP.mat) %in% row.names(geneExpDiff.pair)
      sGP.mat <- GP.mat[idx,idx]

      # sort genes before comparsion
      geneExpDiff.pair <- geneExpDiff.pair[order(row.names(geneExpDiff.pair)),]
      sGP.mat <- sGP.mat[order(row.names(sGP.mat)), order(colnames(sGP.mat))]

      # calculate interation relationship
      geneExp.drugA <- geneExpDiff.pair[,1]
      geneExp.drugB <- geneExpDiff.pair[,2]

      int.mat <- outer(geneExp.drugA, geneExp.drugB, FUN = "*")
      diag(sGP.mat) <- diag(int.mat) <- 0

      iPos <- sum(sGP.mat + int.mat == 2)
      iNeg <- sum(sGP.mat + int.mat == -2)
      iNon <- sum(sGP.mat * int.mat == -1)
      iSim <- iPos + iNeg

      # count all up- and down-regulated genes induced by drug B (used as denominator for normalization in the next step)
      countB.updown <- sum(abs(geneExp.drugB) > 0)
      iPara <- list(iPos = iPos, iNeg = iNeg, iNon = iNon, iSim = iSim, countB.updown = countB.updown)

      return(iPara)
}






#' Profiling Drug Treated Gene Expression Data
#'
#' This function profile the drug treated gene expression data to prepare the input for DIGREscore function.
#'
#' @param geneExp a data frame contains the drug treated gene expression data with each column representing 
#' one drug, and each row representing one gene. See `geneExp.demo` for example.
#' @return a matrix of processed gene expression data.
#'
#' @author Jichen Yang, Sangin Lee, Minzhe Zhang(\email{zenroute.mzhang@gmail.com})
#' @details
#' Gene expression is measured after cell is treated by a single compound (or negative control).
#' Raw data (micro array or RNA-Seq) should already be log transformed to have proper scale. 
#' This function will average duplicated data with same drug name, and collapse multiple probes to gene level.
#'
#' @examples 
#' geneExpDiff <- profileGeneExp(geneExp = geneExp.demo)
#'
#' @export
#' @importFrom preprocessCore normalize.quantiles


profileGeneExp <- function(geneExp) {
      
      cat("Start parse drug treated gene expression data ...\n")
      cat("------------\n")
      
      ### Average duplicated drug data
      geneExp.mat <- geneExp.aveDrug(geneExp)

      ### Normalize data
      geneExp.mat <- normGeneExp(geneExp.mat)

      ### Subtract nagtive control effect
      geneExp.mat <- geneExp.diffCtl(geneExp.mat)

      ### Convert probe level to gene level
      geneExp.profile <- prob2gene(geneExp.mat)

      cat("Done parsing drug treated gene expression data.\n")
      return(geneExp.profile)

}


### Average duplicated drug data in gene expression file
geneExp.aveDrug <- function(geneExp) {

      cat("Checking drug name ...\n")
      ### Check duplicated drugs
      if (any(duplicated(as.character(geneExp[1,-1])))) {
            cat("Duplicated drug name found, average duplicated data.\n")
      }

      ### Print drug list
      drugName <- sort(unique(as.character(geneExp[1,-1])))
      if (!("Neg_control" %in% drugName)) {
            stop("Error: Negative control not found. Plase involve Neg_control column in the file.\n")
      }
      cat("Drug list you provided:\n")
      idx <- 1
      for (i in 1:length(drugName)) {
            if (drugName[i] != "Neg_control") {
                  cat(idx, ". ", drugName[i], "\n", sep = "")
                  idx <- idx + 1
            }
      }

      ### Average dupliacte drug data (if any)
      geneExp.mat <- data.frame(matrix(NA, nrow = nrow(geneExp)-1, ncol = length(drugName)+1))
      colnames(geneExp.mat) <- c("Gene", drugName)
      geneExp.mat$Gene <- geneExp[-1,1]

      for (i in 1:length(drugName)) {
            tmpDrug.data <- geneExp[-1, geneExp[1,] == drugName[i]]
            tmpDrug.data <- apply(as.data.frame(tmpDrug.data), 2, as.numeric)
            tmpDrug.mean <- apply(tmpDrug.data, 1, mean)
            geneExp.mat[drugName[i]] <- tmpDrug.mean
      }

      return(geneExp.mat)
}


### Quantile normalize gene expression data
normGeneExp <- function(geneExp) {

      cat("Normalizing data ...\n")
      drugName <- colnames(geneExp)[-1]; geneName <- geneExp[,1]
      geneExp.mat <- normalize.quantiles(as.matrix(geneExp[,-1]))
      colnames(geneExp.mat) <- drugName
      geneExp.mat <- data.frame(Gene = geneName, geneExp.mat, stringsAsFactors = FALSE, check.names = FALSE)

      return(geneExp.mat)
}


### Subtract nagtive control effect
geneExp.diffCtl <- function(geneExp) {

      cat("Measuring gene expression difference ...\n")
      neg_control <- geneExp[,"Neg_control"]
      for (i in 2:ncol(geneExp)) {
            geneExp[,i] <- geneExp[,i] - neg_control
      }

      ### delete Neg_control column
      geneExp$Neg_control <- NULL

      return(geneExp)
}


### Conver probe level to gene level
prob2gene <- function(geneExp) {

      cat("Collapse multiple probes to genes ...\n")

      geneName <- sapply(strsplit(geneExp$Gene, " ", fixed = TRUE), "[[", 1)
      geneExp$Gene <- geneName

      geneExp.mean <- aggregate(geneExp[,-1], by = list(Gene = geneExp$Gene), FUN = mean)
      geneExp.max <- aggregate(geneExp[,-1], by = list(Gene = geneExp$Gene), FUN = max)
      geneExp.min <- aggregate(geneExp[,-1], by = list(Gene = geneExp$Gene), FUN = min)

      geneExp.profile <- list(geneExp.mean = geneExp.mean, geneExp.max = geneExp.max, geneExp.min = geneExp.min)
      return(geneExp.profile$geneExp.max)
}


#' Construct Gene Network Matrix
#'
#' This function is to convert a gene-gene interaction table to gene network matrix, which is an input of DIGREscore function.
#'
#' @param geneNet a data frame contains gene-gene interaction. See `geneNetLymph` for example.
#' @return a matrix of gene connectivity matrix
#'
#' @author Jichen Yang, Sangin Lee, Minzhe Zhang(\email{zenroute.mzhang@gmail.com})
#' @details
#' Input gene-gene interaction table should have two columns with gene SYMBOL names. 
#' Each raw represents two connected genes. The interaction is regarded as undirected.
#'
#' @examples 
#' geneNetLymph.mat <- constGeneNet(geneNet = geneNetLymph)
#'
#' @export
#' @import org.Hs.eg.db KEGGgraph

constGeneNet <- function(geneNet) {
      cat("Parsing gene network data ...\n")
      cat("------------\n")
      
      # convert gene SYMBOL name to KEGG id
      source.gene <- sapply(mget(geneNet[,1], org.Hs.egSYMBOL2EG, ifnotfound = NA), "[[", 1)
      target.gene <- sapply(mget(geneNet[,2], org.Hs.egSYMBOL2EG, ifnotfound = NA), "[[", 1)
      geneNet.id <- data.frame(source = source.gene, target = target.gene, stringsAsFactors = FALSE)
      
      # delete NA
      geneNet.id <- geneNet.id[(!(grepl("NA", source.gene))) & (!(grepl("NA", target.gene))),]
      
      # paste "hsa"
      geneNet.id$source <- paste("hsa:", geneNet.id$source, sep = "")
      geneNet.id$target <- paste("hsa:", geneNet.id$target, sep = "")
      
      # construct gene network matrix
      gene.list <- union(geneNet.id$source, geneNet.id$target)
      geneNet.mat <- matrix(0, length(gene.list), length(gene.list), dimnames = list(gene.list, gene.list))
      
      cat("Total", length(gene.list), "nodes, total", nrow(geneNet.id), "connnection.\n")
      
      # loop through each raw in geneNet.id to impute connected genes in matrix
      for (i in 1:nrow(geneNet.id)) {
            # symmetric matrix
            geneNet.mat[geneNet.id[i,1],geneNet.id[i,2]] <- 1
            geneNet.mat[geneNet.id[i,2],geneNet.id[i,1]] <- 1
      }
      
      # delete self connectivity (if any)
      diag(geneNet.mat) <- 0
      
      cat("Done constructing gene network matrix.\n")
      return(geneNet.mat)
}



#' Visualization of compound pairs synergistic score
#'
#' This function plot heatmap and barplot of compound pairs synergistic score.
#'
#' @param pred.pair data frame of score rank table generated by DIGREscore function.
#' @param type specify "heat" or "bar" to plot whether heatmap or barplot
#' @return ggplot object of the heatmap or barplot of compound pairs synergistic score.
#'
#' @author Jichen Yang, Sangin Lee, Minzhe Zhang(\email{zenroute.mzhang@gmail.com})
#' @details
#' Please directly apply this funtion to the result got from DIGREscore funtion.
#'
#' @examples 
#' vis.heat <- DIGREvis(pred.pair = res.KEGG$scoreRank, type = "heat")
#' vis.bar <- DIGREvis(pred.pair = res.KEGG$scoreRank, type = "bar")
#' plot(vis.heat)
#' plot(vis.bar)
#'
#' @export
#' @import ggplot2

DIGREvis <- function(pred.pair, type) {
      
      if (type == "heat") {
            
            # plot heatmap
            
            p <- ggplot(pred.pair, aes(drugA, drugB)) + geom_tile(aes(fill = Score), colour = "white") + scale_fill_gradient(low = "white", high = "steelblue")
            p.title <- ggtitle("Heatmap of compound pair synergistic scores\n")
            p.theme <- theme(plot.title = element_text(size = 18, lineheight = 0.8, face = "bold", hjust = 0.5), 
                       axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"), 
                       axis.text.y = element_text(colour = "black"), 
                       axis.title.x = element_text(face = "bold"), 
                       axis.title.y = element_text(face = "bold"))
            p.axis <- scale_x_discrete(labels = abbreviate)
      
            return(p + p.title + p.theme + p.axis)
            
      } else if (type == "bar") {
            
            # plot barplot
            
            
            rankScore <- pred.pair[with(pred.pair, order(Rank)),]
            
            if (nrow(rankScore) <= 20) {
                  num_plot <- nrow(rankScore)
            } else {num_plot <- 20}
            
            rankScore <- data.frame(DrugPair = mapply(function(x,y) 
                  paste(x, y, sep = " & "), abbreviate(rankScore$drugA), abbreviate(rankScore$drugB)), 
                  Score = rankScore$Score, stringsAsFactors = FALSE)[1:num_plot,]
            
            rankScore$DrugPair <- factor(rankScore$DrugPair, levels = rankScore$DrugPair)
            
            p <- ggplot(rankScore, aes(DrugPair, Score)) + geom_bar(stat = "identity", fill = "steelblue")
            p.title <- ggtitle("Top compound pair synergistic scores\n")
            p.theme <- theme(plot.title = element_text(size = 18, lineheight = 0.8, face = "bold", hjust = 0.5), 
                             axis.text.x = element_text(size = 12, angle = 45, hjust = 1, colour = "black"), 
                             axis.text.y = element_text(size = 12, colour = "black"), 
                             axis.title.x = element_text(size = 12, face = "bold"), 
                             axis.title.y = element_text(size = 12, face = "bold"), 
                             plot.margin = unit(c(0.5,0.5,0.2,0.2), "cm"))
            
            return(p + p.title + p.theme)
      }
}



#' Drug treated gene expression data
#'
#' Examplary gene expression data of DLBCL cell after treated with 14 different drugs.
#'
#' @details
#' This data is from NCI-DREAM challenge competition for predicting drug pairs synergy. 
#' OCI-LY3 human diffuse large B-cell lymphoma (DLBCL) cell line was treated by 14 different 
#' drugs in its dose of IC20. 24 hours after Perturbation, gene expression level was measured.
#' @seealso
#' \url{http://www.nature.com/nbt/journal/v32/n12/full/nbt.3052.html}
#'
#' @docType data
#' @usage data(geneExp.demo)
#' @keywords datasets

"geneExp.demo"


#' Dose response data
#'
#' Examplary drug dose response data of 14 drugs in NCI-DREAM challenge.
#'
#' @details
#' This demo data contains dose response data of 14 drugs from NCI-DREAM challenge. 
#' It contain the cell viability reduction when treated with drug in two different dose. 
#' One dose is IC20 of the drug, therefore the cell viability reduction is always 0.2 for all drugs. 
#' The other dose is double dose of IC20, this value is infered from the dose response curve of each drug.
#' @seealso
#' \url{http://www.nature.com/nbt/journal/v32/n12/full/nbt.3052.html}
#'
#' @docType data
#' @usage data(doseRes.demo)
#' @keywords datasets

"doseRes.demo"


#' Lymphoma Specific Gene Network
#'
#' Gene-gene interaction information refined from lymphoma patients gene expression data.
#'
#' @details
#' We use gene expression data from lymphoma patients (GSE10846) to construct Lymphoma-specific 
#' gene network. The statistical algorithm sparse partial correlation estimation (SPACE) is used 
#' to infer the network structure from the expression data.
#' @seealso
#' \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10846}
#' \url{https://cran.r-project.org/web/packages/space/index.html}
#'
#' @docType data
#' @usage data(geneNetLymph)
#' @keywords datasets

"geneNetLymph"


#' Universal KEGG pathway Gene Network
#'
#' A matrix contains the global pathway (GP) gene-gene interaction network merged from KEGG pathway database.
#'
#' @details
#' We selected KEGG pathways belonging to genetic information processing, environmental 
#' information processing, cellular processing, and cancer disease. From this set of selected pathways we removed 
#' any pathway with fewer than 10 edges. Finally we merged the remaining 32 KEGG pathways into a global pathway (GP) 
#' which included 11642 interactions among 2322 genes.
#'
#' @docType data
#' @usage data(KEGGnet.mat)
#' @keywords datasets

"KEGGnet.mat"


#' Cancer Growth Pathway (CGP)
#'
#' A matrix contains the Cancer Growth Pathway (CGP) gene-gene interaction network information
#'
#' @details
#' The CGP was built using pathways empirically selected from KEGG pathway database based on our knowledge
#' firstly, and then refined by a small drug combination training dataset. In particular, we selected 12
#' pathways which were highly related to cell growth. Then we remove one of the 12 pathways each time and
#' merged the remaining to build a CGP with which we applied our approach to the external training dataset.
#' Based on the results, we screened out 8 pathways that contributed mostly to the performance and merged
#' them to build final CGP. The 8 KEGG pathways are: aminoacyl-tRNA biosynthesis, MAPK signaling pathway,
#' NF-kappa B signaling pathway, Cell Cycle, p53 signaling pathway, Apoptosis, TGF-beta signaling pathway,
#' Cancer pathway.
#'
#' @docType data
#' @usage data(CGP.mat)
#' @keywords datasets

"CGP.mat"