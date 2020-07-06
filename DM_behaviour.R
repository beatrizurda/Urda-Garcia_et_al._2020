#########################################################################################
# ANALYZING DM BEHAVIOUR
# 
# FINAL VERSION STORES IN ~/Desktop/ANALYSIS
#
# Beatriz Urda García 2020
#########################################################################################

library(nortest)
library(scran)
library(ggplot2)
library(fitdistrplus)
library(data.table)
library(forcats)
library(grid)
library(ggpubr)

setwd("Analysis/")

dismeta <- read.csv("dl_general_info_umls_ncolors.csv",header=T, sep="\t",stringsAsFactors = F, row.names = NULL)
ch1 <- which(dismeta$tis_cat_colors == 'fuchsia')
dismeta$tis_cat_colors[ch1] <- 'magenta3'
ch1 <- which(dismeta$dis_cat_colors == 'grey')
dismeta$dis_cat_colors[ch1] <- 'black'
dislist <- dismeta$disease_name

summar <- read.csv("Disease_gene_variability/summary_disease_variability.csv",header=T, sep="\t",stringsAsFactors = F)
dis_nsamples <- summar[,c("disease","n_patients","n_controls","n_genes")]
alldm <- read.csv("Disease_gene_variability/all_diseases_dm_ratio.csv",header=T, sep="\t",stringsAsFactors = F)
dim(alldm) # dim: 19992    46

patdm <- data.frame("genes" = as.character(alldm$genes)) # dim: 20020     1
condm <- data.frame("genes" = as.character(alldm$genes)) # dim: 20020     1
for(cdis in dislist){
  # cdis <- "Schizophrenia"
  dis_info <- read.csv(paste("Disease_gene_variability/Gene_variability_tables/",cdis,"_gene_variability_table.csv",sep=""),header=T, sep="\t",stringsAsFactors = F)
  # dis_info <- dis_info[,c("rn","dmp","dmc","dm_dif")]
  to_pat <- dis_info[,c("rn","dmp")]
  to_con <- dis_info[,c("rn","dmc")]
  colnames(to_pat)[2] <- cdis
  colnames(to_con)[2] <- cdis
  patdm <- merge(patdm,to_pat, by.x="genes", by.y="rn", all.x=TRUE, all.y=FALSE)
  condm <- merge(condm,to_con, by.x="genes", by.y="rn", all.x=TRUE, all.y=FALSE)
}


melted <- melt(as.data.frame(patdm), id.vars= "genes")
meltedplus <- merge(melted,dis_nsamples, by.x = "variable", by.y = "disease")
# meltedplus <- meltedplus[1:100000,]
shist <- ggplot(meltedplus, aes(genes,value, color=n_patients)) + geom_point()+
  # geom_bar(stat="identity", position="dodge") +
  # geom_bar(stat="identity") +
  xlab("") + ylab("")
shist

require(gridExtra)
p1 <- ggplot(meltedplus,aes(n_patients, value)) + geom_point() +
  geom_smooth(method="lm", se=TRUE, formula= y~x) +
  # stat_cor(label.x.npc = "center", label.y.npc = "bottom") +
  stat_regline_equation(aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
                        formula=y~x, 
                        label.x.npc = "center", label.y.npc = "bottom")+
  xlab("Number of patients") + ylab("DMp") +
  labs(title = 'DMp vs Number of patients') + 
  # theme(legend.position="bottom") +
  theme(plot.title = element_text(hjust = 0.5))
p1

# PLOT FILENAME: Correlation_DMp_Npatients.png

melted <- melt(as.data.frame(condm), id.vars= "genes")
meltedplus <- merge(melted,dis_nsamples, by.x = "variable", by.y = "disease")
p2 <- ggplot(meltedplus,aes(n_controls, value)) + geom_point() +
  geom_smooth(method="lm", se=TRUE, formula= y~x) +
  # stat_cor(label.x.npc = "center", label.y.npc = "bottom") +
  stat_regline_equation(aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
                        formula=y~x, 
                        label.x.npc = "center", label.y.npc = "bottom")+
  xlab("Number of controls") + ylab("DMc") +
  labs(title = 'DMc vs Number of controls') + 
  # theme(legend.position="bottom") +
  theme(plot.title = element_text(hjust = 0.5))
p2

# PLOT FILENAME: Correlation_DMc_Ncontrols.png

