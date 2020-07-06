#########################################################################################
# ANALYZE THE VARIANCE DISTRIBUTION OF GENES AT THE DISEASE LEVEL
# 
# FINAL VERSION STORES IN ~/Desktop/ANALYSIS
#
# Beatriz Urda García 2020
#########################################################################################

setwd("Analysis/")

library(nortest)
library(scran)
library(ggplot2)
library(fitdistrplus)
library(data.table)
library(forcats)
library(grid)
library(ggpubr)
library(fda)

meta <- readRDS("new_meta_jon_grein_good_quality.rds")

dismeta <- read.csv("dl_general_info_umls_ncolors.csv",header=T, sep="\t",stringsAsFactors = F, row.names = NULL)
ch1 <- which(dismeta$tis_cat_colors == 'fuchsia')
dismeta$tis_cat_colors[ch1] <- 'magenta3'
ch1 <- which(dismeta$dis_cat_colors == 'grey')
dismeta$dis_cat_colors[ch1] <- 'black'

set.seed(5)

combat_files <- list.files(path = "DL_RESULTS_GREIN/after_combat_counts/")
combat_diseases <- gsub("_combat_counts.rds","",combat_files)

normalized_files <- list.files(path = "DL_RESULTS_GREIN/normalized_counts/", pattern = "_se_normalized.rds$")
disease_list <- gsub("_se_normalized.rds","",normalized_files)
disease_list <- disease_list

dis_summary <- data.frame(matrix(ncol=16, nrow=0))
colnames(dis_summary) <- c('disease',
                           'var_ks_ntest', 'var_ad_ntest', 'var_ks_ntestlog', 'var_ad_ntestlog',
                           'cv_ks_ntest', 'cv_ad_ntest', 'cv_ks_ntestlog', 'cv_ad_ntestlog',
                           'se_ks_ntest', 'se_ad_ntest', 'se_ks_ntestlog', 'se_ad_ntestlog',
                           'n_patients','n_controls','n_genes')

# Preparing dfs with all the genes (union) as rows and their variability in the diseases as columns (3 df: VAR, CV and SE) 
# 2020 genes and 52 diseases
alldf <- read.csv("Network_building/all_genes_all_diseases_logFC_table.csv",header=T, sep=",",stringsAsFactors = F, row.names = 1)
dim(alldf) # 45 19992
gene_union <- as.character(colnames(alldf))

dis_var <- data.frame(matrix(ncol=0, nrow=19992))
dis_cv <- data.frame(matrix(ncol=0, nrow=19992))
dis_se <- data.frame(matrix(ncol=0, nrow=19992))
dis_dm <- data.frame(matrix(ncol=0, nrow=19992))
dis_var$genes <- gene_union
dis_cv$genes <- gene_union
dis_se$genes <- gene_union
dis_dm$genes <- gene_union

for(dis in disease_list){
  
  pdf(paste("Disease_gene_variability/Plots/",dis,"_gene_variability_figs.pdf",sep=""))
  
  dis = 'BreastCancer'
  # dis = 'AdenomatousPolyps'
  # dis = disease_list[2]
  print(dis)
  
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
  cnsamples <- dim(counts)[2]
  to_remove <- which(colnames(counts) %in% ccontrols)
  pcounts <- counts[,-to_remove]
  ccounts <- counts[,to_remove]
  dim(pcounts)
  dim(ccounts)
  cnpatients <- dim(pcounts)[2]
  cncontrols <- dim(ccounts)[2]
  
  # Computing the variance for each gene in the patients and the controls
  # print(dim(counts))
  # counts[1:5,1:5]
  n_genes <- dim(counts)[1]
  
  # Means
  meanp <- rep(NA,n_genes)
  meanc <- rep(NA,n_genes)
  
  # Variance
  varp <- rep(NA,n_genes)
  varc <- rep(NA,n_genes)
  varrat <- rep(NA,n_genes)
  
  # Coefficient of variation
  coefvarp <- rep(NA,n_genes)
  coefvarc <- rep(NA,n_genes)
  coefvar_rat <- rep(NA,n_genes)
  
  # Standard Error
  sep <- rep(NA,n_genes)
  sec <- rep(NA,n_genes)
  se_rat <- rep(NA,n_genes)
  
  for(k in 1:n_genes){
    meanp[k] <- mean(pcounts[k,])
    meanc[k] <- mean(ccounts[k,])
    
    varp[k] <- var(pcounts[k,])
    varc[k] <- var(ccounts[k,])
    
    coefvarp[k] <- sd(pcounts[k,]) / mean(pcounts[k,])
    coefvarc[k] <- sd(ccounts[k,]) / mean(ccounts[k,])
    
    sep[k] <- sd(pcounts[k,]) / sqrt(cnpatients)
    sec[k] <- sd(ccounts[k,]) / sqrt(cncontrols)
  }
  
  # Computing DM patients and DM controls
  cv2p <- varp / (meanp^2) # squared CV
  cv2c <- varc / (meanc^2)
  dmp <- DM(meanp, cv2p)
  dmc <- DM(meanc, cv2c)
  
  varrat <- varp / varc
  coefvar_rat <- coefvarp / coefvarc
  se_rat <- sep / sec
  dm_dif <- dmp - dmc
  
  tmpdf <-  data.frame(genes = rownames(counts), varrat, stringsAsFactors = FALSE) ; colnames(tmpdf)[2] <- dis
  dis_var = merge(dis_var, tmpdf, by='genes', all.x=TRUE)
  tmpdf <-  data.frame(genes = rownames(counts), coefvar_rat, stringsAsFactors = FALSE) ; colnames(tmpdf)[2] <- dis
  dis_cv = merge(dis_cv, tmpdf, by='genes', all.x=TRUE)
  tmpdf <-  data.frame(genes = rownames(counts), se_rat, stringsAsFactors = FALSE) ; colnames(tmpdf)[2] <- dis
  dis_se = merge(dis_se, tmpdf, by='genes', all.x=TRUE)
  tmpdf <-  data.frame(genes = rownames(counts), dm_dif, stringsAsFactors = FALSE) ; colnames(tmpdf)[2] <- dis
  dis_dm = merge(dis_dm, tmpdf, by='genes', all.x=TRUE)
  
  var_result <- counts
  dim(var_result)
  var_result <- cbind(var_result,meanp)
  var_result <- cbind(var_result,meanc)
  var_result <- cbind(var_result,varp)
  var_result <- cbind(var_result,varc)
  var_result <- cbind(var_result,coefvarp)
  var_result <- cbind(var_result,coefvarc)
  var_result <- cbind(var_result,sep)
  var_result <- cbind(var_result,sec)
  var_result <- cbind(var_result,dmp)
  var_result <- cbind(var_result,dmc)
  var_result <- cbind(var_result,varrat)
  var_result <- cbind(var_result,coefvar_rat)
  var_result <- cbind(var_result,se_rat)
  var_result <- cbind(var_result,dm_dif)
  dim(var_result) # After adding 14 columns
  
  lastcol = dim(var_result)[2]
  just_var_results <- var_result[,c((lastcol-13):lastcol)]
  # just_var_results[,3] <- just_var_results[,1] / just_var_results[,2]
  
  vardf <- as.data.frame(just_var_results)
  vardf$log_varrat <- log(vardf$varrat)
  vardf$log_coefvar_rat <- log(vardf$coefvar_rat)
  vardf$log_se_rat <- log(vardf$se_rat)
  setDT(vardf, keep.rownames = TRUE)[] # Put row names as first column 
  
  # Write the table
  write.table(vardf,file=paste('Disease_gene_variability/Gene_variability_tables/',dis,"_gene_variability_table.csv", sep=""),sep='\t',row.names=FALSE, quote = FALSE)
  write.table(vardf[,c('rn','dm_dif')],file=paste('Disease_gene_variability/Gene_deltaDM_tables/',dis,"_deltadm_table.csv", sep=""),sep='\t',row.names=FALSE, quote = FALSE)
  # Plot the distribution of most and less variable genes diseases vs control
  # par(mfrow=c(3,1))
  require(gridExtra)
  # VAR
  p1 <- ggplot(vardf, aes(varrat)) + 
    geom_density(fill='blue',alpha=0.25) + geom_vline(xintercept = 1, col='blue') +
    xlab("Variance of the patients / Variance of the controls") +
    labs(title = paste("Variance ratio distribution in ",dis,sep=""))
  # print(p1)
  
  # COEFF VARIATION
  p2 <- ggplot(vardf, aes(coefvar_rat)) + 
    geom_density(fill='blue',alpha=0.25) + geom_vline(xintercept = 1, col='blue') +
    xlab("CV of the patients / CV of the controls") +
    labs(title = paste("Coefficient of variation ratio distribution in ",dis,sep=""))
  # print(p2)
  
  # DM 
  p3 <- ggplot(vardf, aes(dm_dif)) + 
    geom_density(fill='blue',alpha=0.25) + geom_vline(xintercept = 0, col='blue') +
    xlab("DM of the patients / DM of the controls") +
    labs(title = paste("Delta DM in ",dis,sep=""))
  # print(p3)
  
  # # STANDARD ERROR
  # p3 <- ggplot(vardf, aes(cv2_rat)) + 
  #   geom_density(fill='blue',alpha=0.25) + geom_vline(xintercept = 1, col='blue') +
  #   xlab("SD of the patients / SD of the controls") +
  #   labs(title = paste("SE ratio distribution in ",dis,sep=""))
  # # print(p3)
  
  grid.arrange(p1, p2, p3, ncol=1)
  
  ### Alternatively we can select the genes with a variability higher than mean + 2standard deviation (top 5% of the curve)
  # Check normality
  # shapiro.test(vardf$varrat) # can't be run -- too many samples 
  par(mfrow=c(3,1))
  print(qqnorm(vardf$varrat, main="Normal Q-Q Plot - Variance ratio")); print(qqline(vardf$varrat, col = 2)) # We have maybe gamma distribution
  print(qqnorm(vardf$coefvar_rat, main="Normal Q-Q Plot - CV ratio")); print(qqline(vardf$coefvar_rat, col = 2))
  print(qqnorm(vardf$dm_dif, main="Normal Q-Q Plot - Delta DM")); print(qqline(vardf$dm_dif, col = 2))
  # print(qqnorm(vardf$se_rat, main="Normal Q-Q Plot - SE ratio")); print(qqline(vardf$se_rat, col = 2))
  
  print(descdist(vardf$varrat, boot=1000)) # More kurtosis and less skewness than gamma distribution and lognormal
  print(descdist(vardf$coefvar_rat, boot=1000))
  print(descdist(as.numeric(vardf$dm_dif), boot=1000))
  # print(descdist(vardf$se_rat, boot=1000))
  
  var_ks_ntest <- try(ks.test(x=vardf$varrat,"pnorm",mean=mean(vardf$varrat), sd=sd(vardf$varrat)), TRUE) # p-value < 2.2e-16 < 0.05 The data is not normally distributed
  cv_ks_ntest <- try(ks.test(x=vardf$coefvar_rat,"pnorm",mean=mean(vardf$coefvar_rat), sd=sd(vardf$coefvar_rat)), TRUE)
  se_ks_ntest <- try(ks.test(x=vardf$se_rat,"pnorm",mean=mean(vardf$se_rat), sd=sd(vardf$se_rat)), TRUE)
  dm_ks_ntest <- try(ks.test(x=vardf$dm_dif,"pnorm",mean=mean(vardf$dm_dif), sd=sd(vardf$se_rat)), TRUE)
  if(class(var_ks_ntest) == "try-error"){
    var_ks_ntest <- NA
  }
  if(class(cv_ks_ntest) == "try-error"){
    cv_ks_ntest <- NA
  }
  if(class(se_ks_ntest) == "try-error"){
    se_ks_ntest <- NA
  }
  if(class(dm_ks_ntest) == "try-error"){
    dm_ks_ntest <- NA
  }
  var_ad_ntest <- ad.test(vardf$varrat)
  cv_ad_ntest <- ad.test(vardf$coefvar_rat)
  se_ad_ntest <- ad.test(vardf$se_rat)
  dm_ad_ntest <- ad.test(vardf$se_rat)
  
  par(mfrow=c(3,1))
  print(qqnorm(vardf$log_varrat, main="Normal Q-Q Plot - log of the variance ratio")); print(qqline(vardf$log_varrat, col = 2))
  print(qqnorm(vardf$log_coefvar_rat, main="Normal Q-Q Plot - log of the CV ratio")); print(qqline(vardf$log_coefvar_rat, col = 2))
  # print(qqnorm(vardf$log_se_rat, main="Normal Q-Q Plot - log of the SE ratio")); print(qqline(vardf$log_se_rat, col = 2))
  
  # Plot the LOG distribution of most and least variable genes diseases vs control
  # VAR
  p1 <- ggplot(vardf, aes(log_varrat)) + 
    geom_density(fill='blue',alpha=0.25) + geom_vline(xintercept = mean(vardf$log_varrat), col='blue') +
    xlab("Variance of the patients / Variance of the controls") +
    labs(title = paste("log variance ratio distribution in ",dis,sep=""))
  # print(p1)
  
  # COEFF VARIATION
  p2 <- ggplot(vardf, aes(log_coefvar_rat)) + 
    geom_density(fill='blue',alpha=0.25) + geom_vline(xintercept = mean(vardf$log_coefvar_rat), col='blue') +
    xlab("CV of the patients / CV of the controls") +
    labs(title = paste("log CV ratio distribution in ",dis,sep=""))
  # print(p2)
  
  # STANDARD ERROR
  p3 <- ggplot(vardf, aes(log_se_rat)) + 
    geom_density(fill='blue',alpha=0.25) + geom_vline(xintercept = mean(vardf$log_se_rat), col='blue') +
    xlab("SE of the patients / SE of the controls") +
    labs(title = paste("log SE ratio distribution in ",dis,sep=""))
  # print(p3)
  
  grid.arrange(p1, p2, p3, ncol=1)
  
  var_ks_ntestlog <-try(ks.test(x=vardf$log_varrat,"pnorm",mean=mean(vardf$log_varrat), sd=sd(vardf$log_varrat)), TRUE)
  cv_ks_ntestlog <-try(ks.test(x=vardf$log_coefvar_rat,"pnorm",mean=mean(vardf$log_coefvar_rat), sd=sd(vardf$log_coefvar_rat)), TRUE)
  se_ks_ntestlog <-try(ks.test(x=vardf$log_se_rat,"pnorm",mean=mean(vardf$log_se_rat), sd=sd(vardf$log_se_rat)), TRUE)
  if(class(var_ks_ntestlog) == "try-error"){
    var_ks_ntestlog <- NA
  }
  if(class(cv_ks_ntestlog) == "try-error"){
    cv_ks_ntestlog <- NA
  }
  if(class(se_ks_ntestlog) == "try-error"){
    se_ks_ntestlog <- NA
  }
  var_ad_ntestlog <- ad.test(vardf$log_varrat)
  cv_ad_ntestlog <- ad.test(vardf$log_coefvar_rat)
  se_ad_ntestlog <- ad.test(vardf$log_se_rat)
  
  to_append = data.frame(disease=dis, 
                         var_ks_ntest=var_ks_ntest$p.value, var_ad_ntest=var_ad_ntest$p.value, var_ks_ntestlog=var_ks_ntestlog$p.value, var_ad_ntestlog=var_ad_ntestlog$p.value,
                         cv_ks_ntest=cv_ks_ntest$p.value, cv_ad_ntest=var_ad_ntest$p.value, cv_ks_ntestlog=cv_ks_ntestlog$p.value, cv_ad_ntestlog=cv_ad_ntestlog$p.value,
                         se_ks_ntest=se_ks_ntest$p.value, se_ad_ntest=var_ad_ntest$p.value, se_ks_ntestlog=se_ks_ntestlog$p.value, se_ad_ntestlog=se_ad_ntestlog$p.value,
                         dm_ks_ntest = dm_ks_ntest$p.value,dm_ad_ntest = dm_ad_ntest$p.value,
                         n_patients=dim(pcounts)[2],n_controls=length(ccontrols),n_genes=n_genes)
  
  # to_append = data.frame(disease=dis,ntest=ntest$p.value,ntest_log=ntestlog$p.value,n_patients=dim(pcounts)[2],n_controls=length(ccontrols),n_genes=n_genes)
  dis_summary <- rbind(dis_summary,  to_append)
  # dis_summary <- rbind(dis_summary, c(dis,ntest$p.value,ntestlog$p.value,dim(pcounts)[2],length(ccontrols),n_genes))
  dev.off()
  
}

write.table(dis_summary,file="Disease_gene_variability/summary_disease_variability.csv",sep='\t',row.names=FALSE, quote = FALSE)

write.table(dis_var,file="Disease_gene_variability/all_diseases_var_ratio.csv",sep='\t',row.names=FALSE, quote = FALSE)
write.table(dis_cv,file="Disease_gene_variability/all_diseases_cv_ratio.csv",sep='\t',row.names=FALSE, quote = FALSE)
write.table(dis_se,file="Disease_gene_variability/all_diseases_se_ratio.csv",sep='\t',row.names=FALSE, quote = FALSE)
write.table(dis_dm,file="Disease_gene_variability/all_diseases_dm_ratio.csv",sep='\t',row.names=FALSE, quote = FALSE)

#### VISUALIZING GENE VARIABILITY AT THE DISEASE LEVEL ####

# Selecting the melted data to plot
#c(1,c(10:20))
melted <- melt(dis_cv[,c(1,c(40:45))],id.var='genes')
melted <- melt(dis_cv,id.var='genes')
melted_se <- melt(dis_se,id.var='genes')

# Subset of data to plot
to_see <- dis_var[,c("genes","CrohnsDisease","BreastCancer","Schizophrenia","ThyroidCancer_Papillary")]
melted_se <- melt(to_see,id.var='genes')
head(to_see)

# Final data to plot -- VAR RATIO
dis_var <- read.csv("Disease_gene_variability/all_diseases_var_ratio.csv",header=T, sep="\t",stringsAsFactors = F)
melted_se <- melt(dis_var,id.var='genes')
melted <- melted_se

# Final data to plot -- MD DIFFERENCE
dis_dm <- read.csv("Disease_gene_variability/all_diseases_dm_ratio.csv",header=T, sep="\t",stringsAsFactors = F)
melted_se <- melt(dis_dm,id.var='genes')
melted <- melted_se
dis_var <- dis_dm

# Adding the dis and tissue category and colors to the melted object
to_color <- dismeta[,c("disease_name","tissue_cat","disease_cat","tis_cat_colors","dis_cat_colors")]
to_merge <- data.frame(disease_name=unique(melted$variable))
to_color <- merge(to_merge, to_color, by = "disease_name", all.x=TRUE)
melted <- merge(melted, to_color, by.x = "variable", by.y ="disease_name",all.x=TRUE )

# Adding the dis and tissue category and colors to the bidimensional df
to_color <- dismeta[,c("disease_name","tissue_cat","disease_cat","tis_cat_colors","dis_cat_colors")]
dis_cols <- bididf$dis_cat_colors
tis_cols <- bididf$tis_cat_colors

### Al density distributions together log
tgp <- ggplot(melted,aes(log(value), fill=variable)) + geom_density(alpha=0.25) +
  xlab("log of the variance ratio") + ylab("density") +
  # theme(legend.position = 'none') + 
  labs(title = 'Density distribution of gene variability', color="Distance and gene set")
tgp

### Al density distributions together
tgp <- ggplot(melted,aes(value, fill=variable)) + geom_density(alpha=0.25) +
  xlab("log of the variance ratio") + ylab("density") +
  # theme(legend.position = 'none') + 
  labs(title = 'Density distribution of gene variability', color="Distance and gene set")
tgp

# Violin plots
tgp <- ggplot(melted,aes(variable,log(value), fill=variable)) + geom_violin(scale="area") +
  xlab("Number of clusters") + ylab("Average Silhouette") +
  # theme(legend.position = 'none') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = 'none') + 
  labs(title = 'Silhouette distribution with different distances and gene sets', color="Distance and gene set")
tgp

# Bloxplot
tgp <- ggplot(melted,aes(x = reorder(variable,value, FUN=mean),y=log(value), fill=variable)) + geom_boxplot() +
  xlab("Number of clusters") + ylab("Average Silhouette") +
  # theme(legend.position = 'none') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = 'none') + 
  labs(title = 'Silhouette distribution with different distances and gene sets', color="Distance and gene set")
tgp

# means <- melt(summarize_all(dis_var,mean))
# variance <- melt(summarize_all(dis_var,var))
# bidi <- ggplot(data = melt(dis_var)) + geom_point(data=means)
# bidi

### BIDIMENSIONAL PLOTS --- NOT TRANSPOSING
# For each disease, the distribution of gene variability in log scale. Here we see the VAR / MEAN
# FOR VAR
median <- as.data.frame(apply(dis_var, 2, function(x)  median(log(as.numeric(as.character(x))),na.rm=TRUE)))
means <- as.data.frame(apply(dis_var, 2, function(x)  mean(log(as.numeric(as.character(x))),na.rm=TRUE)))
vars <- as.data.frame(apply(dis_var, 2, function(x) var(log(as.numeric(as.character(x))),na.rm=TRUE)))
# FOR DM DIFF
median <- as.data.frame(apply(dis_var, 2, function(x)  median(as.numeric(x),na.rm=TRUE)))
means <- as.data.frame(apply(dis_var, 2, function(x)  mean(as.numeric(x),na.rm=TRUE)))
vars <- as.data.frame(apply(dis_var, 2, function(x) var(as.numeric(x),na.rm=TRUE)))

# median(dis_var$AdenomatousPolyps,na.rm=TRUE) # It is OK
head(median)
head(means)
head(vars)
colnames(median) <- "median" ; colnames(means) <- "mean"; colnames(vars) <- "var";
setDT(median, keep.rownames = TRUE)[] ; setDT(vars, keep.rownames = TRUE)[]; setDT(means, keep.rownames = TRUE)[];
median <- median[-c(1),] ; vars <- vars[-1,] ; means <- means[-1,];
bididf <- cbind(median,means$mean)
bididf <- cbind(bididf,vars$var)
colnames(bididf) <- c('diseases', 'median', 'mean', 'var')

require("ggrepel")

to_color <- dismeta[,c("disease_name","tissue_cat","disease_cat","tis_cat_colors","dis_cat_colors")]
bididf <- merge(bididf, to_color, by.x = "diseases", by.y="disease_name", all.x=TRUE)

# 1.1. BIDIMENSIONAL PLOT WITH THE MEAN | COLOURED BY DISEASE CATEGORY
lev_order <- unique(bididf$disease_cat)
dis_cols <- unique(bididf$dis_cat_colors)
names(dis_cols) <- lev_order
bidiplot <- ggplot(bididf,aes(bididf$mean,sqrt(bididf$var), color=disease_cat, label=bididf$diseases)) + 
  geom_point() + 
  geom_label_repel(size=3.5) +
  scale_color_manual(breaks=lev_order,values=dis_cols, name="Disease Category") +
  # geom_text(size=3) +
  guides(color=guide_legend(nrow=5,byrow=TRUE)) +
  xlab("mean of the log density distribution") + ylab("sd of the log density distribution") +
  theme(legend.position = 'bottom') + 
  labs(title = 'Standard deviation versus Mean of the gene variability density distributions') +
  theme(plot.title = element_text(hjust = 0.5))
bidiplot

# 1.2. BIDIMENSIONAL PLOT WITH THE MEAN | COLOURED BY TISSUE CATEGORY
lev_order <- unique(bididf$tissue_cat)
dis_cols <- unique(bididf$tis_cat_colors)
names(dis_cols) <- lev_order
bidiplot <- ggplot(bididf,aes(bididf$mean,sqrt(bididf$var), color=tissue_cat, label=bididf$diseases)) + 
  geom_point() + 
  geom_label_repel(size=3.5) +
  scale_color_manual(breaks=lev_order,values=dis_cols, name="Tissue Category") +
  # geom_text(size=3) +
  guides(color=guide_legend(nrow=5,byrow=TRUE)) +
  xlab("mean of the log density distribution") + ylab("sd of the log density distribution") +
  theme(legend.position = 'bottom') + 
  labs(title = 'Standard deviation versus Mean of the gene variability density distributions') +
  theme(plot.title = element_text(hjust = 0.5))

bidiplot

# 2.1. BIDIMENSIONAL PLOT WITH THE MEDIAN | COLOURED BY DISEASE CATEGORY
lev_order <- unique(bididf$disease_cat)
dis_cols <- unique(bididf$dis_cat_colors)
names(dis_cols) <- lev_order
bidiplot <- ggplot(bididf,aes(bididf$median,sqrt(bididf$var), color=disease_cat, label=bididf$diseases)) + 
  geom_point() + 
  geom_label_repel(size=3.5) +
  scale_color_manual(breaks=lev_order,values=dis_cols, name="Disease Category") +
  # geom_text(size=3) +
  guides(color=guide_legend(nrow=5,byrow=TRUE)) +
  xlab("mean of the density distribution") + ylab("sd of the density distribution") +
  theme(legend.position = 'bottom') + 
  labs(title = 'Standard deviation vs median of the DM diference distribution') +
  theme(plot.title = element_text(hjust = 0.5))
bidiplot

# 2.2. BIDIMENSIONAL PLOT WITH THE MEDIAN | COLOURED BY TISSUE CATEGORY
lev_order <- unique(bididf$tissue_cat)
dis_cols <- unique(bididf$tis_cat_colors)
names(dis_cols) <- lev_order
bidiplot <- ggplot(bididf,aes(bididf$median,sqrt(bididf$var), color=tissue_cat, label=bididf$diseases)) + 
  geom_point() + 
  geom_label_repel(size=3.5) +
  scale_color_manual(breaks=lev_order,values=dis_cols, name="Tissue Category") +
  # geom_text(size=3) +
  guides(color=guide_legend(nrow=5,byrow=TRUE)) +
  xlab("median of the density distribution") + ylab("sd of the density distribution") +
  theme(legend.position = 'bottom') + 
  labs(title = 'Standard deviation vs median of the DM diference distribution') +
  theme(plot.title = element_text(hjust = 0.5))

bidiplot

# BOXPLOTS FOR DM

# 1. Boxplot sorted by median
############ <parece el bueno>
melted$variable <- with(melted, reorder(x = variable, X = value, FUN = function(x) median(x, na.rm = TRUE)))

# 1.1. Coloured by disease
lev_order <- unique(melted$disease_cat)
dis_cols <- unique(melted$dis_cat_colors)
names(dis_cols) <- lev_order

tgp <- ggplot(melted,aes(variable,value, color=disease_cat)) + geom_boxplot() +
  xlab("Diseases") + ylab("density distribution") +
  # theme(legend.position = 'none') + 
  scale_color_manual(breaks=lev_order,values=dis_cols, name="Disease Category") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(color=guide_legend(nrow=5,byrow=TRUE)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1,2,"cm")) +
  theme(legend.position = 'bottom') +
  labs(title = "DM difference distribution | sorted by median")
tgp

# 1.2. Coloured by tissue
lev_order <- unique(melted$tissue_cat)
tis_cols <- unique(melted$tis_cat_colors)
names(tis_cols) <- lev_order

tgp <- ggplot(melted,aes(variable,value, color=tissue_cat)) + geom_boxplot() +
  xlab("Diseases") + ylab("density distribution") +
  # theme(legend.position = 'none') + 
  scale_color_manual(breaks=lev_order,values=tis_cols, name="Disease Category") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(color=guide_legend(nrow=5,byrow=TRUE)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1,2,"cm")) +
  theme(legend.position = 'bottom') +
  labs(title = "DM difference distribution | sorted by median")
tgp
############ </> parece el bueno

# 2. Boxplot sorted by VAR
############ <parece el bueno>
melted$variable <- with(melted, reorder(x = variable, X = value, FUN = function(x) var(x, na.rm = TRUE)))

# 1.1. Coloured by dis
lev_order <- unique(melted$disease_cat)
dis_cols <- unique(melted$dis_cat_colors)
names(dis_cols) <- lev_order

tgp <- ggplot(melted,aes(variable,value, color=disease_cat)) + geom_boxplot() +
  xlab("Diseases") + ylab("density distribution") +
  # theme(legend.position = 'none') +
  scale_color_manual(breaks=lev_order,values=dis_cols, name="Tissue Category") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(color=guide_legend(nrow=5,byrow=TRUE)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1,2,"cm")) + # top, right, down, left
  theme(legend.position = 'bottom') +
  labs(title = 'DM difference distribution | sorted by sd')
tgp

# 1.2. Coloured by tissue
lev_order <- unique(melted$tissue_cat)
tis_cols <- unique(melted$tis_cat_colors)
names(tis_cols) <- lev_order

tgp <- ggplot(melted,aes(variable,value, color=tissue_cat)) + geom_boxplot() +
  xlab("Diseases") + ylab("density distribution") +
  # theme(legend.position = 'none') +
  scale_color_manual(breaks=lev_order,values=tis_cols, name="Tissue Category") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(color=guide_legend(nrow=5,byrow=TRUE)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1,2,"cm")) + # top, right, down, left
  theme(legend.position = 'bottom') +
  labs(title = 'DM difference distribution | sorted by sd')
tgp

#### BOXPLOTS FOR VAR
# 1. Boxplot sorted by median
############ <parece el bueno>
melted$variable <- with(melted, reorder(x = variable, X = value, FUN = function(x) median(x, na.rm = TRUE)))

# 1.1. Coloured by disease
lev_order <- unique(melted$disease_cat)
dis_cols <- unique(melted$dis_cat_colors)
names(dis_cols) <- lev_order

tgp <- ggplot(melted,aes(variable,log(value), color=disease_cat)) + geom_boxplot() +
  xlab("Diseases") + ylab("log of the density distribution") +
  # theme(legend.position = 'none') + 
  scale_color_manual(breaks=lev_order,values=dis_cols, name="Disease Category") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(color=guide_legend(nrow=5,byrow=TRUE)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1,2,"cm")) +
  theme(legend.position = 'bottom') +
  labs(title = "Summary of the density distribution of gene variability of diseases | sorted by median")
tgp

# 1.2. Coloured by tissue
lev_order <- unique(melted$tissue_cat)
tis_cols <- unique(melted$tis_cat_colors)
names(tis_cols) <- lev_order

tgp <- ggplot(melted,aes(variable,log(value), color=tissue_cat)) + geom_boxplot() +
  xlab("Diseases") + ylab("log of the density distribution") +
  # theme(legend.position = 'none') + 
  scale_color_manual(breaks=lev_order,values=tis_cols, name="Disease Category") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(color=guide_legend(nrow=5,byrow=TRUE)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1,2,"cm")) +
  theme(legend.position = 'bottom') +
  labs(title = "Summary of the density distribution of gene variability of diseases | sorted by median")
tgp
############ </> parece el bueno

# 2. Boxplot sorted by VAR
############ <parece el bueno>
melted$variable <- with(melted, reorder(x = variable, X = value, FUN = function(x) var(x, na.rm = TRUE)))

# 1.1. Coloured by dis
lev_order <- unique(melted$disease_cat)
dis_cols <- unique(melted$dis_cat_colors)
names(dis_cols) <- lev_order

tgp <- ggplot(melted,aes(variable,log(value), color=disease_cat)) + geom_boxplot() +
  xlab("Diseases") + ylab("log of the density distribution") +
  # theme(legend.position = 'none') +
  scale_color_manual(breaks=lev_order,values=dis_cols, name="Tissue Category") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(color=guide_legend(nrow=5,byrow=TRUE)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1,2,"cm")) + # top, right, down, left
  theme(legend.position = 'bottom') +
  labs(title = 'Summary of the density distribution of gene variability of diseases | sorted by variance')
tgp

# 1.2. Coloured by tissue
lev_order <- unique(melted$tissue_cat)
tis_cols <- unique(melted$tis_cat_colors)
names(tis_cols) <- lev_order

tgp <- ggplot(melted,aes(variable,log(value), color=tissue_cat)) + geom_boxplot() +
  xlab("Diseases") + ylab("log of the density distribution") +
  # theme(legend.position = 'none') +
  scale_color_manual(breaks=lev_order,values=tis_cols, name="Tissue Category") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(color=guide_legend(nrow=5,byrow=TRUE)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1,2,"cm")) + # top, right, down, left
  theme(legend.position = 'bottom') +
  labs(title = 'Summary of the density distribution of gene variability of diseases | sorted by variance')
tgp
############ </> parece el bueno

### EXPLORING THE CORRELATIONS BETWEEN VARIANCE AND THE N_SAMPLES OR THE MEAN EXPRESSION OF THE GENES
summar <- read.csv("Disease_gene_variability/summary_disease_variability.csv",header=T, sep="\t",stringsAsFactors = F)
dis_nsamples <- summar[,c("disease","n_patients","n_controls","n_genes")]

var_rats <- read.csv("Disease_gene_variability/all_diseases_var_ratio.csv",header=T, sep="\t",stringsAsFactors = F)

cdis <- "AdenomatousPolyps"
cdis <- "CrohnsDisease"
k <- 1
pdf("Disease_gene_variability/var_mean_correlation_controls.pdf")
for(page in c(1:5)){
  par(mfrow=c(3,4))
  while(k %% 12 != 0){
    cdis <- disease_list[k]
    dis_info <- read.csv(paste("Disease_gene_variability/Gene_variability_tables/",cdis,"_gene_variability_table.csv",sep=""),header=T, sep="\t",stringsAsFactors = F)
    
    p1 <- ggplot(dis_info,aes(meanc, varc)) + geom_point() +
      geom_smooth() +
      xlab("Mean expression of the gene") + ylab("Variance") +
      labs(title = paste('Variance vs mean expression of the controls | ', cdis, sep="")) +
      theme(plot.title = element_text(hjust = 0.5))
    print(p1)
    k <- k+1
    
  }
}
dev.off()

cdis <- "Schizophrenia"
dis_info <- read.csv(paste("Disease_gene_variability/Gene_variability_tables/",cdis,"_gene_variability_table.csv",sep=""),header=T, sep="\t",stringsAsFactors = F)

### CORRELATION BETWEEN MEAN EXPRESSION OF GENES AND DM FOR CONTROLS AND PATIENTS

require(gridExtra)
p1 <- ggplot(dis_info,aes(dmp, meanp)) + geom_point() +
  geom_smooth(method="lm") +
  xlab("DM") + ylab("Mean expression of the gene") +
  labs(title = 'Mean gene expression vs DM - patients') +
  theme(plot.title = element_text(hjust = 0.5))
p1

# library(ggpubr)
# ggscatter(dis_info, x = "dmp", y = "meanp", add = "reg.line") +
#   stat_cor(label.x = 3, label.y = 34) +
#   stat_regline_equation(label.x = 2, label.y = 1)


p2 <- ggplot(dis_info,aes(dmc, meanc)) + geom_point() +
  geom_smooth(method=lm) +
  xlab("DM") + ylab("Mean expression of the gene") +
  labs(title = 'Mean gene expression vs DM - controls') +
  theme(plot.title = element_text(hjust = 0.5))
p2

### CORRELATION BETWEEN MEAN EXPRESSION OF GENES AND VAR FOR CONTROLS AND PATIENTS
p3 <- ggplot(dis_info,aes(varp, meanp)) + geom_point() +
  geom_smooth(method=lm) +
  xlab("VAR") + ylab("Mean expression of the gene") +
  labs(title = 'Mean gene expression vs var - patients') +
  theme(plot.title = element_text(hjust = 0.5))


p4 <- ggplot(dis_info,aes(varc, meanc)) + geom_point() +
  geom_smooth(method=lm) +
  xlab("VAR") + ylab("Mean expression of the gene") +
  labs(title = 'Mean gene expression vs var - controls') +
  theme(plot.title = element_text(hjust = 0.5))


grid.arrange(p3, p4, p1, p2, ncol=2, top="Effect of gene expression in gene variability metrics")

### Same with ggscatter
p1 <- ggscatter(dis_info, x = "dmp", y = "meanp",
          color = "black", shape = 21, size = 3, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.sep = "\n"),
          xlab = "DM", ylab = "Mean expression",
          title = 'Mean gene expression vs DM | patients'
            
) 

p2 <- ggscatter(dis_info, x = "dmc", y = "meanc",
                color = "black", shape = 21, size = 3, # Points color, shape and size
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE, # Add confidence interval
                cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                cor.coeff.args = list(method = "pearson", label.sep = "\n"),
                xlab = "DM", ylab = "Mean expression",
                title = 'Mean gene expression vs DM | controls'
)

dim(dis_info)
wo_outliers <- dis_info[dis_info$varp < 9e+07, ]
dim(wo_outliers)

p3 <- ggscatter(wo_outliers, x = "varp", y = "meanp",
                color = "black", shape = 21, size = 3, # Points color, shape and size
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE, # Add confidence interval
                cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                cor.coeff.args = list(method = "pearson", label.sep = "\n"),
                xlab = "Var", ylab = "Mean expression",
                title = 'Mean gene expression vs Var | patients'
)

p4 <- ggscatter(wo_outliers, x = "varc", y = "meanc",
                color = "black", shape = 21, size = 3, # Points color, shape and size
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE, # Add confidence interval
                cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                cor.coeff.args = list(method = "pearson", label.sep = "\n"),
                xlab = "Var", ylab = "Mean expression",
                title = 'Mean gene expression vs Var | controls'
)

grid.arrange(p3, p4, p1, p2, ncol=2, top="Effect of gene expression in gene variability metrics")


dm_df <- c("Disease"=c(), "dmp"=c() , "dmc"=c())
for (cdis in disease_list){
  dis_info <- read.csv(paste("Disease_gene_variability/Gene_variability_tables/",cdis,"_gene_variability_table.csv",sep=""),header=T, sep="\t",stringsAsFactors = F)
  to_add <- cbind(rep(cdis,dim(dis_info)[1]), dis_info[,c("dmp","dmc")])
  colnames(to_add) <- c("Disease", "dmp", "dmc")
  dm_df <- rbind(dm_df,to_add)
}

# lev_order <- unique(melted$disease_cat)
# dis_cols <- unique(melted$dis_cat_colors)
# names(dis_cols) <- lev_order

tgp <- ggplot(dm_df,aes(Disease,dmp)) + geom_boxplot() +
  xlab("Diseases") + ylab("DM distribution") +
  # theme(legend.position = 'none') +
  # scale_color_manual(breaks=lev_order,values=tis_cols, name="Tissue Category") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(color=guide_legend(nrow=5,byrow=TRUE)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1,2,"cm")) + # top, right, down, left
  theme(legend.position = 'bottom') +
  labs(title = 'Boxplot of the DM values for the patients among diseases')
tgp

# FILENAME PLOT: boxplot_patients_DM_values.png

tgp <- ggplot(dm_df,aes(Disease,dmc)) + geom_boxplot() +
  xlab("Diseases") + ylab("DM distribution") +
  # theme(legend.position = 'none') +
  # scale_color_manual(breaks=lev_order,values=tis_cols, name="Tissue Category") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(color=guide_legend(nrow=5,byrow=TRUE)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1,2,"cm")) + # top, right, down, left
  theme(legend.position = 'bottom') +
  labs(title = 'Boxplot of the DM values for the controls among diseases')
tgp

# FILENAME PLOT: boxplot_controls_DM_values.png

# GENERATE .rnks FOR ALL DISEASES USING DELTADM (DIFFERENCE)
for (cdis in disease_list){
  dis_info <- read.csv(paste("Disease_gene_variability/Gene_variability_tables/",cdis,"_gene_variability_table.csv",sep=""),header=T, sep="\t",stringsAsFactors = F)
  dm_way <- dis_info[c("rn","dm_dif")]
  dm_way <- dm_way[order(-dm_way$dm_dif),]
  colnames(dm_way) <- ""
  write.table(dm_way,file=paste("Disease_gene_variability/ranked_diseases_delta_dm/",cdis,"_delta_dm_values.rnk",sep=""),sep='\t',row.names=FALSE, col.names=FALSE,quote = FALSE)
}


# PREPARE .rnk for Jon for SCHIZOPRENIA
cdis <- "Schizophrenia"
dis_info <- read.csv(paste("Disease_gene_variability/Gene_variability_tables/",cdis,"_gene_variability_table.csv",sep=""),header=T, sep="\t",stringsAsFactors = F)
dm_way <- dis_info[c("rn","dm_dif")]
var_way <- dis_info[c("rn","varrat")]
dm_way <- dm_way[order(-dm_way$dm_dif),]
var_way <- var_way[order(-var_way$varrat),]
colnames(dm_way) <- ""
colnames(var_way) <- ""
write.table(dm_way,file="Disease_gene_variability/Schizophrenia/dm_difference_values.rnk",sep='\t',row.names=FALSE, quote = FALSE)
write.table(var_way,file="Disease_gene_variability/Schizophrenia/var_ratio_values.rnk",sep='\t',row.names=FALSE, quote = FALSE)



