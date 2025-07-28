## 准备环节 --------------------------------------------------------------------
# 载入R包，设置参数，其中job变量用于项目输出文件前缀标识，可以自定义修改。
# getwd()
# setwd("/Users/wendy/Desktop/final_edition/limma_allGene")
options(stringsAsFactors = F)
rm(list=ls()) # 清空变量
job <- "UCEC"
library(limma)
library(ggplot2) # 用于绘制火山图
library(ggVolcano) 

## 整理需要循环的所有癌症和变量
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

# for (n in 1:25) {
#   job <- job_names[[n]]
  
## 数据导入 --------------------------------------------------------------------
# 导入样本信息和表达量数据，然后进行删除表达量之和为0基因、log2化、替换异常值等步骤，得到原始数据矩阵。
expr_data <- read.csv("UCEC_HiSeqV2.csv", 
                    header = TRUE, sep = "", 
                    stringsAsFactors = FALSE, 
                    check.names = FALSE) #输入文件TPM原始值，行名是基因，列名是样本
# expr_data <- as.numeric(expr_data)
expr_data <- expr_data[which(rowSums(expr_data)!=0),] #删除表达量为0的基因
expr_data = log2(expr_data) #log化处理
expr_data[expr_data == -Inf] = 0 #将log化后的负无穷值替换为0
head(expr_data)
group<-read.csv("UCEC_experiment.csv",
                header = T, sep = "") #输入文件，样本信息表，包含分组信息
head(group)
# tutorial: extract the summarized experiment object.
# table(colData(data)$sample_type)

## 构建分组矩阵--design --------------------------------------------------------
# 根据样本的分组信息，构建分组矩阵，最终得到的design矩阵由0和1构成，为斜对角矩阵。
design <- model.matrix(~0+factor(group$experiment))
head(design)

colnames(design) <- levels(factor(group$experiment))
head(design)
rownames(design) <- colnames(expr_data)
head(design)

## 构建比较矩阵--contrast ------------------------------------------------------
# 这一步的好处是更加focus在normal和cancer组之间的差异，而不是所有的组之间的两两比较（e.g.cancer/normal组内部的比较是不必要的）
# 设置样本的比较方式，这里为CK对照比HT处理，该步骤生成的文件为1和-1构成的矩阵。
contrast.matrix <- makeContrasts(Cancer - Normal, levels = design) #根据实际的样本分组修改，这里对照组CK，处理组HT
contrast.matrix

## 线性混合模拟-----------------------------------------------------------------
# 该步骤是limma包的核心步骤，首先使用lmFit函数进行非线性最小二乘法分析，然后用经验贝叶斯调整t-test中方差部分，得到差异表达结果。
fit <- lmFit(expr_data,design) #非线性最小二乘法
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)#用经验贝叶斯调整t-test中方差的部分

## allGene =====================================================================
DEG <- topTable(fit2, coef = 1,n = Inf,sort.by="logFC")
DEG <- na.omit(DEG)
write.csv(DEG, paste0(job,"_","allGene.csv"))
# 最终生成的DEG文件包含以下几列信息：
# > colnames(DEG)
# [1] "logFC"     "AveExpr"   "t"         "P.Value"   "adj.P.Val" "B"        
# [7] "regulate"  "Genes" 
# DEG
