library(SummarizedExperiment)
library(TCGAbiolinks)
library(EDASeq)
library(sesameData)
library(minfi)
library(limma)
library(MultiAssayExperiment)
library(dplyr)

# Read Gene expression data
query <- GDCquery(
  project = "TCGA-READ",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor")

GDCdownload(query = query)
dataPrep1_4 <- GDCprepare(query = query, save = F)
dataPrep_4 <- TCGAanalyze_Preprocessing(object = dataPrep1_4,
                                        cor.cut = 0.6)
dataNorm_4 <- TCGAanalyze_Normalization(tabDF = dataPrep_4,
                                        geneInfo = geneInfoHT,
                                        method = "geneLength")
dataFilt_4 <- TCGAanalyze_Filtering(tabDF = dataNorm_4,
                                    method = "quantile",
                                    qnt.cut =  0.25)
sampleTypes_4 <- c(rep('READ', ncol(dataFilt_4)))

# READ DNA methylation
query <- GDCquery(
  project = "TCGA-READ",
  data.category = "DNA Methylation",
  platform = "Illumina Human Methylation 450",
  data.type = "Methylation Beta Value",
  sample.type = "Primary Tumor")
GDCdownload(query = query)
dataPrep1_9 <- GDCprepare(query = query, save = F)
betaMatrix <- assays(dataPrep1_9)[[1]]
betaMatrixClean <- na.omit(betaMatrix)
cpgs <- rowSds(as.matrix(betaMatrixClean))
quantileCutoff <- quantile(cpgs, probs = 0.75)
informativeCpGs <- names(cpgs)[cpgs > quantileCutoff]
featureMatrix <- betaMatrixClean[informativeCpGs, ]
dim(featureMatrix)
normFeatureMatrix <- normalizeBetweenArrays(as.matrix(featureMatrix), method = "quantile")
DNA_READ = t(normFeatureMatrix)


# merge gene expression and DNA methylation for READ
colnames_dataFilt_4 <- colnames(dataFilt_4)
colnames_normFeatureMatrix <- colnames(normFeatureMatrix)
participant_ids_dataFilt_4 <- substr(colnames_dataFilt_4, 9, 12)
participant_ids_normFeatureMatrix <- substr(colnames_normFeatureMatrix, 9, 12)
common_participant_ids <- intersect(participant_ids_dataFilt_4, participant_ids_normFeatureMatrix)
common_columns_dataFilt_4 <- colnames_dataFilt_4[participant_ids_dataFilt_4 %in% common_participant_ids]
common_columns_normFeatureMatrix <- colnames_normFeatureMatrix[participant_ids_normFeatureMatrix %in% common_participant_ids]
# dataFilt_4_common <- dataFilt_4[, common_columns_dataFilt_4]
# normFeatureMatrix_common <- normFeatureMatrix[, common_columns_normFeatureMatrix]
common_columns_dataFilt_4_sorted <- sort(common_columns_dataFilt_4)
common_columns_normFeatureMatrix_sorted <- sort(common_columns_normFeatureMatrix)
dataFilt_4_common_sorted <- dataFilt_4_common[, common_columns_dataFilt_4_sorted]
normFeatureMatrix_common_sorted <- normFeatureMatrix_common[, common_columns_normFeatureMatrix_sorted]
merged_data <- rbind(dataFilt_4_common_sorted, normFeatureMatrix_common_sorted)
merged_READ <- t(merged_data)
dim(merged_READ)

# Read clinical data
query <- GDCquery(
  project       = "TCGA-READ",
  data.category = "Clinical",
  data.type     = "Clinical Supplement",
  data.format   = "BCR Biotab"    
)

GDCdownload(query = query)
clinical_data <- GDCprepare(query = query)

fu <- clinical_data$clinical_follow_up_v1.0_read

time <- with(fu, ifelse(is.na(death_days_to), death_days_to, 
                        last_contact_days_to))[- c(1,2)]
status <- fu$vital_status[- c(1,2)]
samples <- fu$bcr_patient_barcode[- c(1,2)]

time <- as.numeric(time)
status <- ifelse(status == "Dead", 1, 0) 

tmp <- data.frame(samples, time, status)
tmp <- tmp[!is.na(tmp$time), ]

df <- tmp %>%
  group_by(samples) %>%
  # 选 OS.time 最大的那条记录
  slice_max(time, n = 1, with_ties = FALSE) %>%
  ungroup()

ge <- t(dataFilt_4)

# match samples with ge
# 提取ge的每一个行名的前12位
rownames(ge) <- substr(rownames(ge), 1, 12)

matched_samples <- intersect(df$samples, rownames(ge))
matched_df <- df[df$samples %in% matched_samples, ]
matched_ge <- ge[matched_samples, ]

# 检查是否对齐
all(rownames(matched_ge) == matched_df$samples)




# ===== LAML =====
query <- GDCquery(
  project = "TCGA-BLCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor")


GDCdownload(query = query)
dataPrep1_4 <- GDCprepare(query = query, save = F)
dataPrep_4 <- TCGAanalyze_Preprocessing(object = dataPrep1_4,
                                        cor.cut = 0.6)
dataNorm_4 <- TCGAanalyze_Normalization(tabDF = dataPrep_4,
                                        geneInfo = geneInfoHT,
                                        method = "geneLength")
dataFilt_4 <- TCGAanalyze_Filtering(tabDF = dataNorm_4,
                                    method = "quantile",
                                    qnt.cut =  0.25)

query <- GDCquery(
  project       = "TCGA-BLCA",
  data.category = "Clinical",
  data.type     = "Clinical Supplement",
  data.format   = "BCR Biotab"    
)

GDCdownload(query = query)

clinical_data <- GDCprepare(query = query)

fu <- clinical_data$clinical_follow_up_v4.0_nte_blca






X <- matched_ge
time <- matched_df$time
status <- matched_df$status

X_center <- colMeans(X_final_noNA)
X <- sweep(X_final_noNA, 2, X_center)

READ_drpls <- val_DRPLS_cv(X, time, status, k = 5,
                           ncomp_candidates = 1:10, auc_time_grid = NULL,
                           use_preselection = TRUE, p_thresh = 0.05,
                           auc_method = "NNE")

READ_drpcapls <- val_DRPCAPLS_cv(X, time, status, k = 5,
                                 ncomp_candidates = 1:10, auc_time_grid = NULL,
                                 alpha_candidates = seq(0, 1, by = 0.1),
                                 use_preselection = TRUE, p_thresh = 0.05,
                                 auc_method = "NNE")

READ_pcox <- val_partial_cox_cv(X, time, status, k =5,
                                ncomp_candidates = 2:2, auc_time_grid = NULL,
                                use_preselection = TRUE, p_thresh = 0.05,
                                auc_method = "NNE")

READ_pcapcox <- val_partial_cox_pca_cv(X, time, status, k = 5,
                                 ncomp_candidates = 2:10, auc_time_grid = NULL,
                                 alpha_candidates = seq(0, 1, by = 0.1),
                                 use_preselection = TRUE, p_thresh = 0.05,
                                 auc_method = "NNE")

READ_wpca_pcox <- val_partial_cox_pca_weight_cv(X, time, status, k = 5,
                                         ncomp_candidates = 2:10, auc_time_grid = NULL,
                                         use_preselection = TRUE, p_thresh = 0.05,
                                         auc_method = "NNE")


