### Preparation ================================================================
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#
#BiocManager::install("DRviaSPCN")

###Load DRviaSPCN package
library(DRviaSPCN)
###Clarify the target cancer
job <- "UCEC"

## Download DRviaSPCN Data package ----------------------------------------------

# DRviaSPCNData package: contains the essential data DrugSPESCMatrix and 
# DrugSPPvalueMatrix which are subpathways weighted-ES induced by all drugs and 
# statistic significance (p-value) of centrality score of subpathways regulated 
# by all drugs.

###Download DRviaSPCNData package from GitHub
library(devtools)
install_github("hanjunwei-lab/DRviaSPCNData",force = TRUE)
library(DRviaSPCNData)
###Get weighted-ES of subpathways
DrugSPESCMatrix<-GetData("DrugSPESCMatrix")
###Get p-value of subpathways centrality score
DrugSPPvalueMatrix<-GetData("DrugSPPvalueMatrix")


### Data =======================================================================
## Example 1 : Calculating the centrality scores of subpathways ----------------

# The function “CalCentralityScore” was used to calculate the centrality scores 
# of SPs to reflect the crosstalk influence, which were used as weights in the 
# calculation of drug-disease reverse association score.

###Load depend package
library(data.table) # help read files
library(tidyverse) # help adjust formats
###Read csv file of gene expression profile and adjust the format
GEP <- fread("UCEC_HiSeqV2.csv")
GEP <- GEP %>%
  remove_rownames() %>%
  column_to_rownames(var = 'gene_name')
GEP <- as.matrix(GEP)
###Read csv file of the corresponding sample class label and adjust the format
Slabel <- fread("UCEC_experiment.csv")
Slabel <- Slabel %>%
  remove_rownames() %>%
  column_to_rownames(var = 'id')
Slabel <- as.vector(Slabel$experiment)
Slabel <- ifelse(Slabel == "Normal", 0, 1)

###Load depend package
library(igraph)
###Obtain input data
#print(class(GEP))   # "data.frame"
#rint(class(Slabel))   # "numeric"
###Run the function
CentralityScoreResult <- CalCentralityScore(ExpData=GEP,Label=Slabel,nperm = 1000)
###View all subpathways result
CentralityScoreResult
###Store subpathways results
write.csv(CentralityScoreResult,paste0(job,"_","CentralityScoreResult.csv"))

## Example 2 : Calculating the drug-disease reverse association score and ------ 
## corresponding p-value of drugs 

# The function Optimaldrugs is used to calculate the DES and statistic 
# significance of drugs. The detailed algorithm can be seen in the introduction 
# part. Users could screen out the optimal therapeutic drugs according to a 
# specific threshold. Here we provide weighted and unweighted methods to 
# calculate the score, which can be selected by parameters weight = ’’ . 
# The screening criteria of the up- and down-regulated subpathways can be 
# changed through the parameters pcut = ’’ and topcut = ’’. 

###Run the function
Opdrugresult<-Optimaldrugs(ExpData=GEP,Label=Slabel,DrugSPESC=DrugSPESCMatrix,
                           CentralityScore=CentralityScoreResult,nperm=1000,topcut=10,
                           pcut=0.01,weight=FALSE)
###view all drugs result
Opdrugresult
###Store optimal drug results
write.csv(Opdrugresult,paste0(job,"_","Opdrugresult.csv"))

### Plot a heatmap of the subpathways that are regulated by disease ============

# The function Disease2SPheatmap plots a heat map of the subpathways that are 
# regulated by disease. The input is the result of function CalCentralityScore, 
# disease gene expression profile and sample class in the expression profile. 
# We map subpathways to the disease gene expression through ssgsea to get a 
# subpathway abundance matrix. Then we visualize the matrix by heatmap. Users 
# could change the threshold that is used to screen significant subpathways 
# through the param pcut.

###Load depend package
library(GSVA)
library(pheatmap)

###Run the function
#png(paste0(job,"_regulated_","pheatmap.png"))  # resolution is a problem
Disease2SPheatmap(CentralityScore=CentralityScoreResult,ExpData=GEP,Label=Slabel,pcut=0.05,
                  bk=c(-2,2),cluster.rows=FALSE,cluster.cols=FALSE,
                  show.rownames=TRUE,show.colnames=FALSE,
                  col=c("navy","firebrick3"),cell.width=NA,
                  cell.height=NA,scale="row",fontsize=7,
                  fontsize.row=7,fontsize.col=10)
dev.off()

### References =================================================================
# https://cran.r-project.org/web/packages/DRviaSPCN/vignettes/DRviaSPCN.html
# https://rdrr.io/cran/DRviaSPCN/f/vignettes/DRviaSPCN.Rmd
