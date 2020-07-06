#####################################################################################
# COMPUTING NETWORK OVERLAP AT THE DISEASE LEVEL
#
#   
#  
#  
# 
# Beatriz Urda García 2020
######################################################################################


library(igraph)
library(ggplot2)

# COMMAND-LINE ARGUMENTS
args = commandArgs(trailingOnly=TRUE)

# Network filename
filename <- args[1] 

# Keep all positive interactions OR consistent
if(args[2] == 'keep_all_pos'){kepp_all_pos_interactions <- TRUE}else{kepp_all_pos_interactions <- FALSE}
if (kepp_all_pos_interactions == TRUE){
  handle_dup_interactions <- '_allpos'
}else{
  handle_dup_interactions <- '_consistent'
}

# Universe option
universe_option <- args[3]
randomizing_labels <-  args[4]


setwd("Analysis/")

# Read the networks
# MANUAL ARGUMENTS
# filename <- 'pairwise_spearman_sDEGs_pos_network.txt'        # PEARWISE SPEARMAN sDEGs
# filename <- 'pairwise_spearman_distance_sDEGs_pos_network.txt'        # PEARWISE SPEARMAN sDEGs
# # filename <- 'union_spearman_sDEGs_pos_network.txt'           # SPEARMAN UNION of sDEGs for ALL DISEASES
# # filename <- 'pairwise_union_spearman_sDEGs_pos_network.txt'    # SPEARMAN UNION of sDEGs for ALL DISEASES
# kepp_all_pos_interactions <- TRUE
# if (kepp_all_pos_interactions == TRUE){
#   handle_dup_interactions <- '_allpos'
# }else{
#   handle_dup_interactions <- '_consistent'
# }
# universe_option <- "B"
# randomizing_labels <-  TRUE

# SCRIPT BEGINS
filedir <- 'Network_building/Defined_networks/'
filepath <- paste(filedir,filename,sep="")
if(randomizing_labels == TRUE){
  outpath <- paste('Network_building/Overlapping_results/Shuffling_labels/',gsub('_network.txt',"",filename),handle_dup_interactions,'_',universe_option,'_',randomizing_labels,'_overlap.txt',sep="")
}else{
  outpath <- paste('Network_building/Overlapping_results/Sampling_randomly/',gsub('_network.txt',"",filename),handle_dup_interactions,'_',universe_option,'_',randomizing_labels,'_overlap.txt',sep="")
}

network <- read.csv(filepath,header=T, sep="\t",stringsAsFactors = F)
filecon <- file(outpath)
writeLines(paste('OVERLAP ANALYSIS OF THE NETWORK:\t',filename,sep=""),filecon)
close(filecon)
write(paste('\n',args[2],'\t','Universe option: ',universe_option,sep=""),file=outpath, append=TRUE)

# ALL - JON
# network <- read.csv("Network_building/Defined_networks/union_spearman_sDEGs_network.txt",header=T, sep="\t",stringsAsFactors = F)
# network <- read.csv("Network_building/Defined_networks/pairwise_spearman_sDEGs_network.txt",header=T, sep="\t",stringsAsFactors = F)
# dim(network) # 241   5



dismeta <- read.csv("dl_general_info_umls_ncolors.csv",header=T, sep="\t",stringsAsFactors = F, row.names = NULL)
dismeta$icd9 <- as.character(dismeta$icd9)

# Add the icds that are missing
dismeta[dismeta$disease_name == 'ThyroidCancer_Papillary', ]$icd9 <- '193' 
dismeta[dismeta$disease_name == 'ThyroidCancer_Papillary', ]$icd10 <- 'C73'

# Put three letters where there are only 2 before the point.
dismeta[dismeta$disease_name == 'HIV', ]$icd9 <- '042' 
dismeta[dismeta$disease_name == 'BorreliaBurdorferiInfection', ]$icd9 <- '088' 

dismeta$icd9 <- substr(dismeta$icd9, start = 1, stop = 3)
dismeta$icd10 <- substr(dismeta$icd10, start = 1, stop = 3)

# Save new dismeta
# write.table(dismeta,file="dl_general_info_umls_ncolors_fixed.csv",sep="\t",row.names=F)

dis_wo_sDEGs <- dismeta[dismeta$Number_sDEGs == 0,]
dis_wo_sDEGs <- dis_wo_sDEGs[, c('disease_name','icd10','icd9')]

# Create a df with diseases and their icds, keeping only first 3 letters for the icd codes.
icds <- dismeta[,c('disease_name','icd10','icd9')]


# Obtain the networks in icds
merged <- merge(network, icds, by.x = "Dis1", by.y="disease_name", all.x = TRUE, all.y=FALSE)
dim(merged)
colnames(merged)[6:7] <- c("Dis1_icd10","Dis1_icd9")

merged <- merge(merged, icds, by.x = "Dis2", by.y="disease_name", all.x = TRUE, all.y=FALSE)
dim(merged)
colnames(merged)[8:9] <- c("Dis2_icd10","Dis2_icd9")
head(merged)


write(paste('\nN interactions in the network:\t',dim(merged)[1],sep=""),file=outpath, append=TRUE)

##### Select network based on icd10
# icd10_network <- merged[,c('Dis1_icd10','Dis2_icd10','Distance')]
# for_jon <- icd10_network
# for_jon$Distance[which(for_jon$Distance > 0)] <- 1
# for_jon$Distance[which(for_jon$Distance < 0)] <- -1 
# for (k in 1:length(for_jon$Distance)){
#   if (for_jon$Dis1_icd10[k] > for_jon$Dis2_icd10[k]){
#     dis1 <- for_jon$Dis1_icd10[k]
#     for_jon$Dis1_icd10[k] <- for_jon$Dis2_icd10[k]
#     for_jon$Dis2_icd10[k] <- dis1
#   }
#     
# }
# # Union sDEGs
# write.table(for_jon,file=paste("Network_building/Defined_networks/","icd10_binarized_union_sDEGs.txt",sep=""),sep="\t",row.names=F)
# # Pairwise sDEGs
# write.table(for_jon,file=paste("Network_building/Defined_networks/","icd10_binarized_pairwise_sDEGs.txt",sep=""),sep="\t",row.names=F)

##### Select network based on icd9
sort_icd_interactions <- function(df){
  # Sort interactions such that icd of Dis1 < icd of Dis2
  
  for (k in 1:length(df$Dis1)){
    if (df$Dis1[k] > df$Dis2[k]){
      dis1 <- df$Dis1[k]
      df$Dis1[k] <- df$Dis2[k]
      df$Dis2[k] <- dis1
    }
  }
  return(df)
}

# Read and filter barabasi network
barabasi <- read.csv("Network_building/PDN_3_digits.net",header=F, sep="\t",stringsAsFactors = F)
colnames(barabasi) <- c('Dis1','Dis2','Prev1','Prev2','Co-ocurrence','RR','RR_99perc_left_bound','RR_99perc_right_bound','phi_corr','p_value')
barabasi <- barabasi[,c('Dis1','Dis2','RR_99perc_left_bound')]
barabasi <- barabasi[which(barabasi$RR_99perc_left_bound > 1 ), ]
dim(barabasi)
write(paste('Barabasi network initial interactions:\t',dim(barabasi)[1],sep=""),file=outpath, append=TRUE)
# Sort barabasi network
barabasi <- sort_icd_interactions(barabasi)

# Obtain icd9-based network with sorted interactions. That is: icd9 Dis1 < icd9 Dis2
cnetwork <- merged[,c('Dis1_icd9','Dis2_icd9','Distance')]
names(cnetwork) <- c('Dis1','Dis2','Distance')
cnetwork <- sort_icd_interactions(cnetwork)
length(cnetwork$Distance) # 241

# Remove interactions where icd1 = icd2 ¡¡NEW!!
cnetwork <- cnetwork[!(cnetwork$Dis1 == cnetwork$Dis2), ]
dim(cnetwork)

# Remove duplicated interactions
if(kepp_all_pos_interactions == TRUE){
  concat <- paste(cnetwork$Dis1,cnetwork$Dis2,sep="_")
  length(concat)             # 241
  length(unique(concat))     # 167
  cnetwork$concat <- concat
  cnetwork <- cnetwork[!duplicated(cnetwork$concat),]
  dim(cnetwork)
}else{
  c_icds <- union(unique(cnetwork$Dis1), unique(cnetwork$Dis2)) # los que tengo en la red: con interaccion signif con alguna enfermedad
  length(c_icds) # 26
  n_times_icd9 <- table(dismeta$icd9)
  n_times_icd9 <- n_times_icd9[n_times_icd9 > 1]
  dups_icd9 <- names(n_times_icd9)
  dups_icd9 <- intersect(dups_icd9,c_icds)
  
  concat <- paste(cnetwork$Dis1,cnetwork$Dis2,sep="_")
  length(concat)             # 241
  length(unique(concat))     # 167
  cnetwork$concat <- concat
  
  dup_cnetwork <- cnetwork[(cnetwork$Dis1 %in% dups_icd9) | (cnetwork$Dis2 %in% dups_icd9),]
  dim(dup_cnetwork) # 155
  
  k <- 1
  to_remove <- c()
  for (dup in dups_icd9){
    # dup <- dups_icd9[1]
    cnet <- dup_cnetwork[(dup_cnetwork$Dis1 == dup) | (dup_cnetwork$Dis2 == dup),]
    table <- table(cnet$concat)
    n_times_icd9[[dup]]
    ftable <- table[table < n_times_icd9[[dup]]]
    to_remove <- append(to_remove, names(ftable))
  }
  to_remove <- unique(to_remove)
  length(to_remove)
  
  cnetwork <- cnetwork[!cnetwork$concat %in% to_remove,]
  cnetwork <- cnetwork[!duplicated(cnetwork$concat),]
}

cnetwork <- cnetwork[,c("Dis1","Dis2","Distance")]
dim(cnetwork)

# Compute number of interactions
length(cnetwork$Dis1) # 167
My_interactions <- length(cnetwork$Dis1)
write(paste('Unique icd9 interactions:\t',My_interactions,sep=""),file=outpath, append=TRUE)

# Obtain overlap with Barabasi
# my_icds <- union()
barabasi_icds <- unique(union(barabasi$Dis1,barabasi$Dis2))
length(barabasi_icds)

# Quitar las enfermedades que no existen en la otra red -- Remove from Barabasi the interactions that contain a disease that I don't have

if(universe_option == "A"){
  # OPTION A. If I consider all the diseases for which I have run the analysis
  all_my_icds <- unique(dismeta$icd9)
  length(all_my_icds) # 44 
  write(paste('icd9 codes in our data:\t',length(all_my_icds),sep=""),file=outpath, append=TRUE)
  if(filename =='pairwise_spearman_distance_sDEGs_pos_network.txt'){
    all_my_icds <- setdiff(all_my_icds, dis_wo_sDEGs$icd9) # Remove the diseases with 0 sDEGs ¿¿¿¿¿DO I HAVE TO DO THIS??? 
    cnetwork <- cnetwork[(cnetwork$Dis1 %in% all_my_icds) & (cnetwork$Dis2 %in% all_my_icds), ] # Filtering my network ¡¡ NEW !!
    write(paste('icd9 codes in our data for which we have sDEGs:\t',length(all_my_icds),sep=""),file=outpath, append=TRUE)
  }
  length(all_my_icds)
}else if(universe_option == "B"){
  # OPTION B. Considering the diseases in my network currently
  c_icds <- unique(union(cnetwork$Dis1, cnetwork$Dis2))
  all_my_icds <- c_icds
  write(paste('icd9 codes in our network (option B):\t',length(all_my_icds),sep=""),file=outpath, append=TRUE)
}else{
  # OPTION C. Considering the diseases in my network with at least 3sDEGs
  min_sDEGs <- 1
  power_dis <- dismeta[dismeta$Number_sDEGs >= min_sDEGs, ] ; power_dis <- unique(power_dis$icd9)
  length(power_dis) # 36 always --> 28 sin duplicates
  all_my_icds <-power_dis
  num_all_my_icds <- as.numeric(all_my_icds) 
  cnetwork <- cnetwork[(cnetwork$Dis1 %in% num_all_my_icds) & (cnetwork$Dis2 %in% num_all_my_icds), ] # Filtering my network
  My_interactions <- length(cnetwork$Dis1)
  num_all_my_icds <- intersect(num_all_my_icds, as.numeric(unique(union(cnetwork$Dis1, cnetwork$Dis2))))
}

dim(barabasi)
num_all_my_icds <- as.numeric(all_my_icds); num_all_my_icds <- unique(num_all_my_icds)
num_all_my_icds <- intersect(num_all_my_icds,unique(union(barabasi$Dis1,barabasi$Dis2))) # Remove icd codes not in Barabasi network
barabasi <- barabasi[(barabasi$Dis1 %in% num_all_my_icds) & (barabasi$Dis2 %in% num_all_my_icds), ] # Filtering Barabasi's network
dim(barabasi) # 160
write(paste('Barabasi interaction within our icd9 codes:\t',dim(barabasi)[1],sep=""),file=outpath, append=TRUE)

cnetwork$Dis1 <- as.numeric(cnetwork$Dis1)
cnetwork$Dis2 <- as.numeric(cnetwork$Dis2)

# Saving the networks
finalnet_filename <- gsub('_overlap.txt',"_final_network.txt",outpath)
finalbarabasi_filename <- gsub('_overlap.txt',"_final_barabasi.txt",outpath)
write.table(cnetwork,file=finalnet_filename,sep="\t",row.names=F, quote = FALSE)
write.table(barabasi,file=finalbarabasi_filename,sep="\t",row.names=F, quote = FALSE)


mygraph <- graph_from_data_frame(cnetwork, directed=FALSE)
barabasigraph <- graph_from_data_frame(barabasi, directed=FALSE)
intersect <- intersection(mygraph,barabasigraph)
# write_graph(intersect ,file="last_computed_overlap.txt", format ="ncol")
gsize(intersect) # Number of edges
write(paste('OVERLAP - Interactions in both networks:\t',gsize(intersect),sep=""),file=outpath, append=TRUE)
perc_overlap <- (gsize(intersect)/dim(cnetwork)[1])*100
write(paste('Percentage of overlap:\t',perc_overlap,sep=""),file=outpath, append=TRUE)
# 64 interactions in the olverlap

# Generate a network with all my possible interactions (52 diseases)
N <- length(all_my_icds) # 31
# n_possible_interactions <- sum(range(N-1)) # 60

connected_net <- data.frame(Dis1=c(),Dis2=c(),Distance=c())
set.seed(5)
for (k1 in 1:(N-1)){
  for (k2 in (k1+1):N){
    if(is.na(num_all_my_icds[k2])== T){
      print(k1)
      print(k2)
      break
    }
    connected_net <- rbind(connected_net,data.frame(num_all_my_icds[k1],num_all_my_icds[k2],1))
  }
}
# Sort the connected network
colnames(connected_net) <- c('Dis1','Dis2','Distance')
connected_net <- sort_icd_interactions(connected_net)
dim(connected_net) # ((31*31) - 31)/2
# 465 possible interactions
write(paste('Dimension of the connected network:\t',dim(connected_net)[1],sep=""),file=outpath, append=TRUE)

# Generate 10.000 iterations of: select random network with my nº of interactions and compute overlap with Barabasi
n_iters <- 10000 
random_overlaps <- c()
set.seed(5) ### ESTABLECER SEMILLA

if(randomizing_labels == FALSE){
  for(k in 1:n_iters){
    crandom_indexes <- sample(c(1:length(connected_net$Dis1)),My_interactions)
    crandom_net <- connected_net[crandom_indexes,]
    cgraph <- graph_from_data_frame(crandom_net, directed=FALSE)
    cintersect <- intersection(cgraph,barabasigraph)
    gsize(cintersect) # Number of edges
    random_overlaps <- append(random_overlaps,gsize(cintersect))
  }
}

# Randomizing labels
if(randomizing_labels == TRUE){
  for(k in 1:n_iters){
    cgraph <- rewire(mygraph, with=keeping_degseq(loops=FALSE, niter = ecount(mygraph)*10))
    # str(cgraph)
    cintersect <- intersection(cgraph,barabasigraph)
    # cintersect
    gsize(cintersect)
    random_overlaps <- append(random_overlaps,gsize(cintersect))
    # print_all(rewire(mygraph, with = keeping_degseq(loops=FALSE, niter = vcount(mygraph) * 10)))
  }
}


pval <- length(which(random_overlaps>=gsize(intersect)))/length(random_overlaps)
pval
write(paste('P-value:\t',pval,sep=""),file=outpath, append=TRUE)
# 0.1076
             
hist(random_overlaps, main = "Random overlaps - pairwise spearman sDEGs")
abline(v = gsize(intersect), col="red")

# t.test(random_overlaps,gsize(intersect))
# p <- ggplot(as.data.frame(random_overlaps), aes()) + 




