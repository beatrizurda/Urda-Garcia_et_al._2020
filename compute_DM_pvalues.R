#########################################################################################
# ANALYZE THE VARIANCE DISTRIBUTION OF GENES AT THE DISEASE LEVEL
# 
# FINAL VERSION STORES IN ~/Desktop/ANALYSIS
#
# Beatriz Urda García 2020
#########################################################################################

# To run in the Mare
args <- commandArgs(trailingOnly = TRUE)
dis <- args[1]

# library(nortest)
library(scran) # INSTALL IN THE MARE
library(ggplot2)
library(fitdistrplus)
library(data.table)
library(forcats) # FALTA
library(grid)
library(ggpubr) # FALTA
library(fda) # INSTALL IN THE MARE

# To test in my computer
# dis = 'BreastCancer'      # dis ='ChronicLymphocyticLeukemia'  # dis ='Cardiomyopathy'
# print(dis)
# setwd("~/Desktop/ANALYSIS/")


meta <- readRDS("new_meta_jon_grein_good_quality.rds")

dismeta <- read.csv("dl_general_info_umls_ncolors_fixed.csv",header=T, sep="\t",stringsAsFactors = F, row.names = NULL)
ch1 <- which(dismeta$tis_cat_colors == 'fuchsia')
dismeta$tis_cat_colors[ch1] <- 'magenta3'
ch1 <- which(dismeta$dis_cat_colors == 'grey')
dismeta$dis_cat_colors[ch1] <- 'black'

combat_files <- list.files(path = "DL_RESULTS_GREIN/after_combat_counts/")
combat_diseases <- gsub("_combat_counts.rds","",combat_files)
normalized_files <- list.files(path = "DL_RESULTS_GREIN/normalized_counts/", pattern = "_se_normalized.rds$")


set.seed(5)

cmeta <- meta[meta$disease == dis,]
ccontrols <- cmeta[cmeta$type == 'control',]
ccontrols <- as.character(ccontrols$sample)

# Get counts
if(dis %in% combat_diseases){
  counts <- readRDS(paste("DL_RESULTS_GREIN/after_combat_counts/",dis,"_combat_counts.rds",sep=""))
  used_genes='combat_counts'
}else{
  cobject <- readRDS(paste("DL_RESULTS_GREIN/normalized_counts/",dis,"_se_normalized.rds",sep=""))
  used_genes='normalized_counts'
  counts <- assays(cobject)$logCPM.norm
}

# Since voom nomalization transforms counts to log2
func = function(x){
  return(2^x)
}
crownames <- rownames(counts)
counts <- sapply(as.data.frame(counts),func)
rownames(counts) <- crownames

# Dividing into patients and controls
dim(counts)
N_samples <- dim(counts)[2]
to_remove <- which(colnames(counts) %in% ccontrols)
pcounts <- counts[,-to_remove]
ccounts <- counts[,to_remove]
dim(pcounts)
dim(ccounts)
N_genes <- dim(counts)[1]

# Open empirical dm_dif values
dm_dif_dir <- "Disease_gene_variability/Gene_deltaDM_tables/"
dm_dif <- read.csv(paste(dm_dif_dir,dis,"_deltadm_table.csv",sep=""),header=T, sep="\t",stringsAsFactors = F, row.names = NULL)


# RANDOMIZATION TESTS TO OBTAIN THE P-VALUE -- THE GOOD ONE ----
N_rand <- 100000 #10000 # 100000
N_patients <- dim(pcounts)[2]
N_controls <- dim(ccounts)[2]
# boot <- rep(0,N_boot)
boot_meansp <- rep(0,N_rand)
boot_cv2p <- rep(0,N_rand)
boot_meansc <- rep(0,N_rand)
boot_cv2c <- rep(0,N_rand)
DDMM <- matrix(nrow = N_genes, ncol=N_rand); # dim(DDMM)
CV2 <- function(var,mean){
  return(var / (mean^2))
}
random_pvalues_drcha <- rep(NA,N_genes)
random_pvalues_izqda <- rep(NA,N_genes)
for(i in 1:N_rand){ # For each bootstrap replication
  to_sample <- sample(c(1:N_samples), replace=T)
  p <- counts[, to_sample[1:N_patients]];
  c <- counts[, to_sample[(N_patients+1):N_samples]];
  
  # Computing delta DM , HACERLO PARA CADA GEN! CON UN APPLY
  b_varp <- apply(p, 1, var); b_varc <- apply(c, 1, var)
  b_meanp <- apply(p, 1, mean); b_meanc <- apply(c, 1, mean)
  b_cv2p <- b_varp / (b_meanp^2) ; b_cv2c <- b_varc / (b_meanc^2)
  
  if(length(unique(colnames(p))) == 1){
    DMp <- rep(0,N_genes)
  }else{
    b_varp <- apply(p, 1, var);
    b_meanp <- apply(p, 1, mean)
    b_cv2p <- b_varp / (b_meanp^2)
    DMp <- DM(b_varp,b_cv2p)
  }
  
  if(length(unique(colnames(c))) == 1){
    DMc <- rep(0,N_genes)
  }else{
    b_varc <- apply(c, 1, var)
    b_meanc <- apply(c, 1, mean)
    b_cv2c <- b_varc / (b_meanc^2)
    DMc <- DM(b_varc,b_cv2c)
  }
  
  # Computing DM and DDM
  DDMM[,i] <- DMp - DMc
  # if(length(which(is.na(DDMM[,i]))) > 0){
  #   break
  # }
  # print(i)
}
# Saving the randomization table
write.table(DDMM,file=paste("Disease_gene_variability/Randomized_matrix/",dis,"randomized_matrix.txt",sep=""),sep='\t',row.names=FALSE, quote = FALSE)
# Computing p-values
for (i in 1:N_genes){
  random_pvalues_drcha[i] <- sum(DDMM[i,] > dm_dif[i,2])/N_rand
  random_pvalues_izqda[i] <- sum(DDMM[i,] < dm_dif[i,2])/N_rand
}
result <- cbind(dm_dif,random_pvalues_izqda)
result <- cbind(result,random_pvalues_drcha)
result$random_pvalues <- random_pvalues_drcha
result[result$dm_dif < 0, ]$random_pvalues <- result[result$dm_dif < 0, ]$random_pvalues_izqda

FDR_corrected <- p.adjust(result$random_pvalues, method="fdr")
# head(FDR_corrected, 50)
length(which(FDR_corrected < 0.05))
# FDR_corrected[which(FDR_corrected < 0.05)]
result <- cbind(result,FDR_corrected)
write.table(result,file=paste("Disease_gene_variability/FDR_corrected/",dis,"FDR_corrected_ddm.txt",sep=""),sep='\t',row.names=FALSE, quote = FALSE)
# </>-- THE GOOD ONE ----c

pdf(paste("Disease_gene_variability/DM_histograms/",dis,"dm_histograms.pdf", sep=""))
# plot(density(FDR_corrected))
hist(random_pvalues)
hist(FDR_corrected)
dev.off()






