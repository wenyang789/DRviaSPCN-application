# Preparation ----
getwd()
setwd("/Users/wendy_1/DRviaSPCN-application/DEG-analysis")

# Set parameters
options(stringsAsFactors = F)
rm(list=ls())  # Clear the variables
job <- "UCEC"  # Set "UCEC" as an example

# Load libraries
library(limma)
library(ggplot2)
library(ggVolcano) 

# Below list other cancer data available on TCGA database.
# ValidRawData <- list('BLCA_HiSeqV2.csv', 'BRCA_HiSeqV2.csv', 'CESC_HiSeqV2.csv', 'CHOL_HiSeqV2.csv', 
#                      'COAD_HiSeqV2.csv', 'COADREAD_HiSeqV2.csv', 'ESCA_HiSeqV2.csv', 'HNSC_HiSeqV2.csv', 
#                      'KICH_HiSeqV2.csv', 'KIRC_HiSeqV2.csv', 'KIRP_HiSeqV2.csv', 'LIHC_HiSeqV2.csv', 
#                      'LUAD_HiSeqV2.csv', 'LUNG_HiSeqV2.csv', 'LUSC_HiSeqV2.csv', 'PADD_HiSeqV2.csv', 
#                      'PCPG_HiSeqV2.csv', 'PRAD_HiSeqV2.csv', 'READ_HiSeqV2.csv', 'SARC_HiSeqV2.csv', 
#                      'SKCM_HiSeqV2.csv', 'STAD_HiSeqV2.csv', 'THCA_HiSeqV2.csv', 'THYM_HiSeqV2.csv', 
#                      'UCEC_HiSeqV2.csv')
# Experiments <- list('BLCA_experiment.csv', 'BRCA_experiment.csv', 'CESC_experiment.csv', 'CHOL_experiment.csv', 
#                     'COAD_experiment.csv', 'COADREAD_experiment.csv', 'ESCA_experiment.csv', 'HNSC_experiment.csv', 
#                     'KICH_experiment.csv', 'KIRC_experiment.csv', 'KIRP_experiment.csv', 'LIHC_experiment.csv', 
#                     'LUAD_experiment.csv', 'LUNG_experiment.csv', 'LUSC_experiment.csv', 'PADD_experiment.csv', 
#                     'PCPG_experiment.csv', 'PRAD_experiment.csv', 'READ_experiment.csv', 'SARC_experiment.csv', 
#                     'SKCM_experiment.csv', 'STAD_experiment.csv', 'THCA_experiment.csv', 'THYM_experiment.csv', 
#                     'UCEC_experiment.csv')
# job_names <- list('BLCA', 'BRCA', 'CESC', 'CHOL', 
#             'COAD', 'COADREAD', 'ESCA', 'HNSC', 
#             'KICH', 'KIRC', 'KIRP', 'LIHC', 
#             'LUAD', 'LUNG', 'LUSC', 'PADD', 
#             'PCPG', 'PRAD', 'READ', 'SARC', 
#             'SKCM', 'STAD', 'THCA', 'THYM', 
#             'UCEC')

  
# Reading data ----

# Read the TPM raw expression data with row × col = gene × samples
expr_data <- read.csv("UCEC_HiSeqV2.csv", 
                    header = TRUE, sep = "", 
                    stringsAsFactors = FALSE, 
                    check.names = FALSE)

# Delete genes that have no expression
expr_data <- expr_data[which(rowSums(expr_data)!=0),]

# Perform log2 transformation and replace negative infinity values with 0
expr_data = log2(expr_data)
expr_data[expr_data == -Inf] = 0
head(expr_data)

# Read group labels of samples
group<-read.csv("UCEC_experiment.csv",
                header = T, sep = "")
head(group)


# Constructing grouping matrix (design matrix) ----
design <- model.matrix(~0+factor(group$experiment))  # The design matrix is composed of 0 & 1 and is a diagonal matrix.
head(design)
colnames(design) <- levels(factor(group$experiment))
head(design)
rownames(design) <- colnames(expr_data)
head(design)


# Constructing contrast matrix ----
# The benefit of this step is to focus more on the differences between the normal and cancer groups, rather than comparing all groups pairwise.
# Comparisons within the cancer/normal group are unnecessary.

# Set the sample comparison method
contrast.matrix <- makeContrasts(Cancer - Normal, levels = design)  # The contrast matrix is a matrix composed of 1 & -1.
contrast.matrix


# Linear hybrid simulation (KEY STEP of limma) ----

# Perform nonlinear least squares analysis
fit <- lmFit(expr_data,design)
fit2 <- contrasts.fit(fit, contrast.matrix)

# Use empirical Bayes to adjust the variance part in the t-test to obtain the differential expression results
fit2 <- eBayes(fit2)


# DEG results ----
DEG <- topTable(fit2, coef = 1,n = Inf,sort.by="logFC")
DEG <- na.omit(DEG)
write.csv(DEG, paste0(job,"_","allGene.csv"))
colnames(DEG)  # Check the column names of the result table.
