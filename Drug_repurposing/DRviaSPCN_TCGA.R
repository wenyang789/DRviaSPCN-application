# Preparation ----
getwd()
setwd("/Users/wendy_1/DRviaSPCN-application/Drug_repurposing")

# If not DRviaSPCN, install the package first.
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#
#BiocManager::install("DRviaSPCN")

# Setting the parameter
job <- "UCEC"  # Set "UCEC"as an example.


# Load required libraries ----
library(DRviaSPCN)
library(devtools)
library(data.table)
library(tidyverse)
library(igraph)
library(pheatmap)


# Downloading DRviaSPCN Data package ----
install_github("hanjunwei-lab/DRviaSPCNData",force = TRUE)
library(DRviaSPCNData)

DrugSPESCMatrix<-GetData("DrugSPESCMatrix")  # Get weighted-ES of subpathways
DrugSPPvalueMatrix<-GetData("DrugSPPvalueMatrix")  # Get p-value of subpathways centrality score


# Calculating the centrality scores of subpathways ----

# Read csv file of gene expression profile and adjust the format
GEP <- fread("UCEC_HiSeqV2.csv")
colnames(GEP)[1] <- "gene_name"
GEP <- GEP %>%
  remove_rownames() %>%
  column_to_rownames(var = 'gene_name')
GEP <- as.matrix(GEP)

# Read csv file of the corresponding sample class label and adjust the format
Slabel <- fread("UCEC_experiment.csv")
colnames(Slabel)[1] <- "id"
Slabel <- Slabel %>%
  remove_rownames() %>%
  column_to_rownames(var = 'id')
Slabel <- as.vector(Slabel$experiment)
Slabel <- ifelse(Slabel == "Normal", 0, 1)

# Calculating the centrality scores
CentralityScoreResult <- CalCentralityScore(ExpData=GEP,Label=Slabel,nperm = 1000)
head(CentralityScoreResult)
write.csv(CentralityScoreResult,paste0(job,"_","CentralityScoreResult.csv"))


# Plotting heatmaps of subpathways that are regulated by disease ----
Disease2SPheatmap(CentralityScore=CentralityScoreResult,ExpData=GEP,Label=Slabel,pcut=0.05,
                  bk=c(-2,2),cluster.rows=FALSE,cluster.cols=FALSE,
                  show.rownames=TRUE,show.colnames=FALSE,
                  col=c("navy","firebrick3"),cell.width=NA,
                  cell.height=NA,scale="row",fontsize=8,
                  fontsize.row=3,fontsize.col=5)
dev.off()


# References ----
# https://cran.r-project.org/web/packages/DRviaSPCN/vignettes/DRviaSPCN.html
# https://rdrr.io/cran/DRviaSPCN/f/vignettes/DRviaSPCN.Rmd
