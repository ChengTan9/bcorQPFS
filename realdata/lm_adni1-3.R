# Preprocess and combine ADNI1, ADNI2, and ADNI3 data.

rm(list = ls())

find_script_directory <- function() {
  local_candidates <- c(".", file.path("code_V1", "realdata"))
  local_candidates <- local_candidates[file.exists(
    file.path(local_candidates, "lm_adni1-3.R")
  )]
  if (length(local_candidates)) return(local_candidates[[1L]])

  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    path <- rstudioapi::getSourceEditorContext()$path
    if (nzchar(path)) return(dirname(path))
  }
  for (frame in rev(sys.frames())) {
    if (!is.null(frame$ofile)) return(dirname(frame$ofile))
  }
  file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_argument)) {
    return(dirname(sub("^--file=", "", file_argument[[1L]])))
  }
  stop("Cannot locate the script directory.", call. = FALSE)
}

# Define input and output directories.
SCRIPT_DIR <- find_script_directory()
RAW_DATA_DIR <- file.path(SCRIPT_DIR, "data", "raw")
ADNI1_SAMPLE_DIR <- file.path(RAW_DATA_DIR, "Hippocampus_all", "sample information")
ADNI1_HIPP_DIR <- file.path(RAW_DATA_DIR, "Hippocampus_all", "hipp_combine_information")
ADNI23_DIR <- file.path(RAW_DATA_DIR, "ADNI2-3")
ADNI23_HIPP_DIR <- file.path(ADNI23_DIR, "hipp_information_radial")
ADNI2_SNP_DIR <- file.path(ADNI23_DIR, "ADNI2_SNP")
ADNI3_SNP_DIR <- file.path(ADNI23_DIR, "ADNI3_SNP")
ADNI13_AUX_DIR <- file.path(RAW_DATA_DIR, "ADNI1-3")
PROCESSED_DATA_DIR <- file.path(SCRIPT_DIR, "..", "..", "code", "data")
ID_OUTPUT_DIR <- file.path(PROCESSED_DATA_DIR, "hipp_id")
dir.create(PROCESSED_DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(ID_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Load and align ADNI1 genotype, phenotype, and clinical data.
hipp_inf <- read.csv(file.path(ADNI1_SAMPLE_DIR, "ADNI1_hip_clinical_information.csv"),sep = ",")
snp_id <- read.csv(file.path(ADNI1_SAMPLE_DIR, "ADNI1_SNP_ID.csv"),sep = ",")
hip_id <- read.csv(file.path(ADNI1_SAMPLE_DIR, "ADNI1_hip_ID.csv"),sep = ",")

apoe_id <- read.csv(file.path(ADNI23_DIR, "ADNIMERGE_25Jan2024.csv"))

match_id <- intersect(as.vector(snp_id$Var1),as.vector(hipp_inf[hipp_inf[,"RACE"] == "White","PTID"]))
match_id <- match_id[-match(c("012_S_0634","032_S_0187","002_S_0295","002_S_1268","031_S_0321","131_S_0409","137_S_0631","141_S_1137"),match_id)]

load(file = file.path(ADNI1_HIPP_DIR, "ADNI1_MATCH_id23.Rdata"))
adni1_snp_index <- match(match_id,ADNI1_MATCH_id23[,1])
ADNI1_MATCH_id23 <- ADNI1_MATCH_id23[adni1_snp_index,]

adni1id <- match_id
save(adni1id,file = file.path(ID_OUTPUT_DIR, "adni1id.Rdata"))

top_pc <- read.table(file = file.path(ADNI1_SAMPLE_DIR, "ADNI1Top20Pcs.txt"))
top_pc <- top_pc[match(match_id,top_pc[,"V1"]),]

cov_inf <- read.csv(file = file.path(ADNI1_SAMPLE_DIR, "ADNI1_hip_clinical_information.csv"),sep = ",",header = TRUE)
cov_inf <- cov_inf[match(match_id,cov_inf[,"PTID"]),]
cov_inf[,"GENDER"] <- ifelse(cov_inf[,"GENDER"] == "Male",1,0)

left_hipp <- read.csv(file = file.path(ADNI1_HIPP_DIR, "left_hipp.csv"))
right_hipp <- read.csv(file = file.path(ADNI1_HIPP_DIR, "right_hipp.csv"))

left_hipp1 <- left_hipp[na.omit(match(match_id,left_hipp[,1])),]
right_hipp1 <- right_hipp[na.omit(match(match_id,right_hipp[,1])),]

save(right_hipp1,file = file.path(PROCESSED_DATA_DIR, "right_hipp.RData"))
save(left_hipp1,file = file.path(PROCESSED_DATA_DIR, "left_hipp.RData"))

row.names(right_hipp1) <- right_hipp1[,1]
row.names(left_hipp1) <- left_hipp1[,1]
right_hipp1 <- right_hipp1[,-1]
left_hipp1 <- left_hipp1[,-1]

top_pc1 <- apply(as.matrix(top_pc[,-c(1,2)]), 2, function(x) as.numeric(x))
cov_inf1.1 <- apply(as.matrix(as.matrix(cov_inf[,c("GENDER","AGE")])), 2, function(x) as.numeric(x))
apoe_index <- match(match_id,apoe_id[,"PTID"])
apoe_inf <- apoe_id[apoe_index,"APOE4"]

cov_inf1 <- scale(cbind(cov_inf1.1,top_pc1[,seq(5)],apoe_inf))

# Load ADNI2 and ADNI3 hippocampal measurements and genotype data.
load(file.path(ADNI23_HIPP_DIR, "left_hipp_adni2.Rdata"))
load(file.path(ADNI23_HIPP_DIR, "right_hipp_adni2.Rdata"))
adni2_left_hipp <- left_hipp_matrix
adni2_right_hipp <- right_hipp_matrix
load(file.path(ADNI23_HIPP_DIR, "left_hipp_adni3.Rdata"))
load(file.path(ADNI23_HIPP_DIR, "right_hipp_adni3.Rdata"))
adni3_left_hipp <- left_hipp_matrix
adni3_right_hipp <- right_hipp_matrix

load(file = file.path(ADNI2_SNP_DIR, "ADNI2_MATCH_id13.Rdata"))
load(file = file.path(ADNI3_SNP_DIR, "ADNI3_MATCH_id12.Rdata"))

# Load and standardize ADNI2 and ADNI3 covariates.
hip_id2 <- read.csv(file.path(ADNI2_SNP_DIR, "ADNI2_hipMADNI1.csv"),sep = ",")
hip_id3 <- read.csv(file.path(ADNI3_SNP_DIR, "ADNI3_SamplehippMADNI1.csv"),sep = ",")
hip_id2 <- na.omit(hip_id2)
hip_id3 <- na.omit(hip_id3)

hip_id2[,"PTGENDER"] <- ifelse(hip_id2[,"PTGENDER"] == "Male",1,0)
hip_id3[,"PTGENDER"] <- ifelse(hip_id3[,"PTGENDER"] == "Male",1,0)

cov_inf2 <- hip_id2[,c("PTGENDER","AGE","PCs1","PCs2","PCs3","PCs4","PCs5")]
cov_inf3 <- hip_id3[,c("PTGENDER","AGE","PCs1","PCs2","PCs3","PCs4","PCs5")]

cov_inf2 <- scale(cov_inf2)
cov_inf3 <- scale(cov_inf3)

# Align subjects across genotype, phenotype, and covariate data.
adni2_inter_index <- intersect(ADNI2_MATCH_id13[,1],intersect(row.names(adni2_left_hipp),hip_id2$PTID))
adni3_inter_index <- intersect(ADNI3_MATCH_id12[,1],intersect(row.names(adni3_left_hipp),hip_id3$PTID))

adni2_inter_index <- adni2_inter_index[-na.omit(match(c("018_S_4399","019_S_4835","022_S_4266","023_S_4034","041_S_4138","941_S_4377","941_S_5124"),adni2_inter_index))]

adni2id <- adni2_inter_index
save(adni2id,file = file.path(ID_OUTPUT_DIR, "adni2id.Rdata"))
adni3id <- adni3_inter_index
save(adni3id,file = file.path(ID_OUTPUT_DIR, "adni3id.Rdata"))

adni2_hipp_index <- match(adni2_inter_index,row.names(adni2_left_hipp))
adni2_snp_index <- match(adni2_inter_index,ADNI2_MATCH_id13[,1])
adni2_cov_index <- match(adni2_inter_index,hip_id2$PTID)
adni2_apoe_index <- match(adni2_inter_index,apoe_id[,"PTID"])

adni3_hipp_index <- match(adni3_inter_index,row.names(adni3_left_hipp))
adni3_snp_index <- match(adni3_inter_index,ADNI3_MATCH_id12[,1])
adni3_cov_index <- match(adni3_inter_index,hip_id3$PTID)
adni3_apoe_index <- match(adni3_inter_index,apoe_id[,"PTID"])

adni2_left_hipp <- adni2_left_hipp[adni2_hipp_index,]
adni2_right_hipp <- adni2_right_hipp[adni2_hipp_index,]
ADNI2_MATCH_id13_index <- ADNI2_MATCH_id13[adni2_snp_index,]

cov_inf2 <- cov_inf2[adni2_cov_index,]
apoe_inf2 <- scale(apoe_id[adni2_apoe_index,"APOE4"])
cov_inf2 <- cbind(cov_inf2,apoe_inf2)

adni3_left_hipp <- adni3_left_hipp[adni3_hipp_index,]
adni3_right_hipp <- adni3_right_hipp[adni3_hipp_index,]
ADNI3_MATCH_id12_index <- ADNI3_MATCH_id12[adni3_snp_index,]

cov_inf3 <- cov_inf3[adni3_cov_index,]
apoe_inf3 <- scale(apoe_id[adni3_apoe_index,"APOE4"])
cov_inf3 <- cbind(cov_inf3,apoe_inf3)

# Combine covariates, hippocampal measurements, and genotypes across phases.
cov_inf1_3 <- rbind(cov_inf1,cov_inf2,cov_inf3)

colnames(right_hipp1) <- colnames(adni2_right_hipp)
colnames(left_hipp1) <- colnames(adni2_left_hipp)
right_hipp_inf1_3 <- rbind(right_hipp1,adni2_right_hipp,adni3_right_hipp)
left_hipp_inf1_3 <- rbind(left_hipp1,adni2_left_hipp,adni3_left_hipp)

colnames(ADNI1_MATCH_id23) <-  colnames(ADNI2_MATCH_id13_index)
colnames(ADNI3_MATCH_id12_index) <-  colnames(ADNI2_MATCH_id13_index)
ADNI_COMBINE <- rbind(ADNI1_MATCH_id23,ADNI2_MATCH_id13_index,ADNI3_MATCH_id12_index)

save(ADNI_COMBINE,file = file.path(PROCESSED_DATA_DIR, "ADNI_1_3_COMBINE1215.Rdata"))

# Residualize hippocampal measurements with the combined covariate matrix.
cov_inf1_3 <- cbind(1,cov_inf1_3)
write.csv(cov_inf1_3,file = file.path(PROCESSED_DATA_DIR, "cov_inf_1215_1_3.csv"))

XX <- cov_inf1_3%*%solve(t(cov_inf1_3)%*%cov_inf1_3)%*%t(cov_inf1_3)
right_res_ori <- as.matrix(right_hipp_inf1_3) - as.matrix(XX)%*%as.matrix(right_hipp_inf1_3)
left_res_ori <- as.matrix(left_hipp_inf1_3) - as.matrix(XX)%*%as.matrix(left_hipp_inf1_3)

save(right_res_ori,file = file.path(PROCESSED_DATA_DIR, "right_res_ori.Rdata"))
save(left_res_ori,file = file.path(PROCESSED_DATA_DIR, "left_res_ori.Rdata"))

# Export combined and residualized data as CSV files.
write.table(ADNI_COMBINE, file = file.path(PROCESSED_DATA_DIR, "ADNI_1_3_COMBINE1215.csv"),sep = ",",row.names = FALSE,col.names = FALSE)
load(file.path(ADNI13_AUX_DIR, "left_G_21.Rdata"))
write.table(y_G, file = file.path(PROCESSED_DATA_DIR, "left_G_21.csv"),sep = ",",row.names = FALSE,col.names = FALSE)
load(file.path(ADNI13_AUX_DIR, "right_G_21.Rdata"))
write.table(y_G, file = file.path(PROCESSED_DATA_DIR, "right_G_21.csv"),sep = ",",row.names = FALSE,col.names = FALSE)

write.table(right_res_ori, file = file.path(PROCESSED_DATA_DIR, "right_res_ori.csv"),sep = ",",row.names = FALSE,col.names = FALSE)
write.table(left_res_ori, file = file.path(PROCESSED_DATA_DIR, "left_res_ori.csv"),sep = ",",row.names = FALSE,col.names = FALSE)

