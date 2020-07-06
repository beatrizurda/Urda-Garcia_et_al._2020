#########################################################################################
# ANALYZE NETWORK'S TOPOLOGY
# 
# 
#
# Beatriz Urda García 2020
#########################################################################################

library(igraph)
library(ggplot2)
library(gplots)

setwd("Analysis/")

# RUN ANALYSIS FOR SSN
network_file <- "metapatients_and_disease/metap_dis_pairwise_union_spearman_distance_sDEGs_network.txt"
output_filename <- "SSN"
obtain_cliques <- FALSE; cex_axis <- 4

# RUN ANALYSIS FOR DSN
network_file <- "pairwise_union_spearman_distance_sDEGs_network.txt"
output_filename <- "DSN"
obtain_cliques <- TRUE; cex_axis <- 6

# RUN ANALYSIS FOR JON'S DATA
barabasi <- read.csv("Network_building/PDN_3_digits.net",header=TRUE, sep="\t",stringsAsFactors = F)
colnames(barabasi) <- c('Dis1','Dis2','Prev1','Prev2','Co-ocurrence','RR','RR_99perc_left_bound','RR_99perc_right_bound','phi_corr','p_value')
barabasi <- barabasi[,c('Dis1','Dis2','RR_99perc_left_bound')]
barabasi <- barabasi[which(barabasi$RR_99perc_left_bound > 1 ), ]
colnames(barabasi) <- c('Dis1','Dis2','Distance')
barabasi$Dis1 <- as.character(barabasi$Dis1); barabasi$Dis2 <- as.character(barabasi$Dis2)
df <- barabasi
output_filename <- "Barabasi"
obtain_cliques <- FALSE; cex_axis <- 4

# RUN ANALYSIS FOR JON'S DATA
# network_file <- "../ICD9_RMS_Sanchez_et_al_2020.txt"
# output_filename <- "Sanchez_et_al"
# obtain_cliques <- FALSE; cex_axis <- 4

net_filepath <- paste("Network_building/Defined_networks/",network_file, sep="")
df <- read.csv(net_filepath,header=T, sep="\t",stringsAsFactors = F); dim(df)
dfpos <- df[df$Distance >= 0, ]
dfneg <- df[df$Distance < 0, ]
graph_list <- list(df,dfpos,dfneg)


k <- 1
for(df in graph_list){
  # head(df)
  # class(df)
  pdf(paste("Network_building/Network_analysis/Topology/",output_filename,k,"_topology_plots.pdf",sep=""))
  graph <- graph_from_data_frame(df, directed=FALSE)
  E(graph)$weight <- df$Distance
  
  connected_components <- clusters(graph)
  n_connected_components <- connected_components$no
  size_connected_components <- connected_components$csize
  if(length(size_connected_components) > 1){
    size_connected_components <- paste(size_connected_components,collapse="_")
  }
  
  deg <- degree(graph, mode="all")
  mean_deg <- mean(deg) 
  
  #Density of a graph: The proportion of present edges from all possible edges in the network
  edge_density <- edge_density(graph, loops=F)
  edge_density
  
  # Reciprocity: it only makes sense for directed graphs
  reciprocity(graph)
  
  # dyad_census(graph) # Mutual, asymmetric, and nyll node pairs
  # 2*dyad_census(graph)$mut/ecount(graph) # Calculating reciprocity
  
  # Transitivity: Transitivity measures the probability that the adjacent vertices of a vertex are connected
  transitivity <- transitivity(graph, type="global")
  local_trans <- transitivity(graph, vids=V(graph), type="local")
  transitivity_table <- data.frame(nodes=names(V(graph)), transitivity=local_trans)
  print(ggplot(transitivity_table, aes(reorder(nodes,transitivity), transitivity)) + geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis),plot.margin = margin(0.7, 0.7, 0.7, 0.7, "cm"))+
    # ggtitle("Transitivity") + 
    # xlab("Diseases") + 
    ylab("Transitivity"))
  
  ### Diameter: longest geodesic distance (length of the shortest path between two nodes) in the network
  # With weights
  diameter <- diameter(graph, directed=F, weights=abs(E(graph)$weight))
  get_diameter(graph, directed=F, weights=abs(E(graph)$weight))
  # Without weights
  diameter_wo_weights <- diameter(graph, directed=F, weights=NA)
  get_diameter(graph, directed=F, weights=NA)
  
  ###Degree
  deg <- degree(graph, mode="all")
  plot(graph, vertex.size=deg*0.4,layout=layout.circle)
  hist(deg, breaks=1:vcount(graph)-1, main="Histogram of node degree")
  
  # Degree distribution
  deg.dist <- degree_distribution(graph, cumulative=T, mode="all")
  
  plot( x=0:max(deg), y=1-deg.dist, pch=19, cex=1.2, col="orange", 
        xlab="Degree", ylab="Cumulative Frequency")
  
  ### CENTRALITY AND CENTRALIZATION
  # Degree (number of ties) ### hacer un plot de esto con ggplot quizá!
  degree_table <- data.frame('nodes'=names(degree(graph, mode="in")), 'feature'=degree(graph, mode="in"))
  centr_degree(graph, mode="in", normalized=T)
  print(ggplot(degree_table, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.7, 0.7, 0.7, 0.7, "cm"))+
    # ggtitle("Transitivity") + 
    # xlab("Diseases") + 
    ylab("Degree"))
  
  # Closeness (centrality based on distance to others in the graph)
  #Inverse of the node’s average geodesic distance to others in the network.
  closeness(graph, mode="all", weights=NA) 
  centr_clo(graph, mode="all", normalized=T)
  degree_table <- data.frame('nodes'=names(closeness(graph, mode="all", weights=NA)), 'feature'=closeness(graph, mode="all", weights=NA))
  print(ggplot(degree_table, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.7, 0.7, 0.7, 0.7, "cm"))+
    # ggtitle("Transitivity") + 
    # xlab("Diseases") + 
    ylab("Closeness"))
  
  # Eigenvector (centrality proportional to the sum of connection centralities)
  # Values of the first eigenvector of the graph matrix.
  
  eigen_centrality(graph, directed=T, weights=NA)
  centr_eigen(graph, directed=T, normalized=T) 
  
  # Betweenness (centrality based on a broker position connecting others)
  # Number of geodesics that pass through the node or the edge.
  betweenness <- betweenness(graph, directed=FALSE, weights=abs(E(graph)$weight))
  betweenness_wo_weights <- betweenness(graph, directed=FALSE, weights=NA)
  edge_betweenness(graph, directed=T, weights=NA)
  centr_betw(graph, directed=T, normalized=T)
  betweenness <- data.frame('nodes'=names(betweenness), 'feature'=betweenness)
  print(ggplot(betweenness, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.7, 0.7, 0.7, 0.7, "cm"))+
    ggtitle("Betweenness with weights") +
    # xlab("Diseases") + 
    ylab("Betweenness"))
  
  betweenness <- data.frame('nodes'=names(betweenness_wo_weights), 'feature'=betweenness_wo_weights)
  print(ggplot(betweenness, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.7, 0.7, 0.7, 0.7, "cm"))+
    ggtitle("Betweenness wo weights") +
    # xlab("Diseases") + 
    ylab("Betweenness"))
  
  ### HUBS AND AUTHORITIES
  # HUBS with weights
  hs <- hub_score(graph, weights=abs(E(graph)$weight))$vector
  hubs_table <- data.frame('nodes'=names(hs), 'feature'=hs)
  print(ggplot(hubs_table, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.7, 0.7, 0.7, 0.7, "cm"))+
    ggtitle("Hubs with weights") +
    # xlab("Diseases") + 
    ylab("Hubs"))
  # wo weights
  hs <- hub_score(graph, weights=NA)$vector
  hubs_table <- data.frame('nodes'=names(hs), 'feature'=hs)
  print(ggplot(hubs_table, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.7, 0.7, 0.7, 0.7, "cm"))+
    ggtitle("Hubs wo weights") +
    # xlab("Diseases") + 
    ylab("Hubs"))
  
  # AUTHORITIES with weights
  as <- authority_score(graph, weights=abs(E(graph)$weight))$vector
  auth_table <- data.frame('nodes'=names(as), 'feature'=as)
  print(ggplot(auth_table, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.7, 0.7, 0.7, 0.7, "cm"))+
    ggtitle("Authorities with weights") +
    # xlab("Diseases") + 
    ylab("Authority"))
  
  # wo weights
  as <- authority_score(graph, weights=NA)$vector
  auth_table <- data.frame('nodes'=names(as), 'feature'=as)
  print(ggplot(auth_table, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.7, 0.7, 0.7, 0.7, "cm"))+
    ggtitle("Authorities wo weights") +
    # xlab("Diseases") + 
    ylab("Authority"))
  
  par(mfrow=c(1,2))
  plot(graph, vertex.size=hs*20, main="Hubs")
  plot(graph, vertex.size=as*10, main="Authorities")
  
  ### DISTANCES AND PATHS
  # Average path length: the mean of the shortest distance between 
  #each pair of nodes in the network (in both directions for directed graphs).
  mean_distance <- mean_distance(graph, directed=F)
  # distances(graph) # with edge weights # YIELDS AND ERROR
  distances(graph, weights=NA) # ignore weights
  
  ### FINDING COMMUNITIES
  # cliques: complete subgraphs of an undirected graph
  if(obtain_cliques == TRUE){
    graph.sym <- as.undirected(graph, mode= "collapse",
                               edge.attr.comb=list(weight="sum", "ignore"))
    cliques(graph) # list of cliques       
    sapply(cliques(graph), length) # clique sizes
    saveRDS(largest_cliques(graph), file=paste("Network_building/Network_analysis/Topology/Cliques/",output_filename,k,'.rds', sep=""))
  }
  
  
  # # Community detection algorithms
  # if (input$comm_algorithm == 'greedy'){
  #   # Using greedy optimization of modularity
  #   fc <- fastgreedy.community(graph)
  #   V(graph)$community <- fc$membership
  # }else if(input$comm_algorithm == 'rand_walks'){
  #   # Using random walks
  #   fc <- cluster_walktrap(graph)
  #   V(graph)$community <- fc$membership #membership(fc)
  # }
  
  # K-core decomposition
  # The k-core is the maximal subgraph in which every node has degree of 
  # at least k. The result here gives the coreness of each vertex in the 
  # network. A node has coreness D if it belongs to a D-core but not 
  # to (D+1)-core.
  kc <- coreness(graph, mode="all")
  # plot(graph, vertex.size=kc*6, vertex.label=kc, vertex.color=colrs[kc]) # FIXX
  kc_table <- data.frame('nodes'=names(kc), 'feature'=kc)
  print(ggplot(kc_table, aes(reorder(nodes,feature), feature)) + geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 45,hjust=1, size=cex_axis), plot.margin = margin(0.7, 0.7, 0.7, 0.7, "cm"))+
    # ggtitle("Authorities with weights") +
    # xlab("Diseases") + 
    ylab("K-core decomposition"))
  
  ### ASSORTATIVITY AND HOMOPHILY
  # assortativity_nominal(graph, V(graph)$media.type, directed=F)
  assortativity_degree <- assortativity_degree(graph, directed=F)
  
  cfeat <- c('N conn comp','Size conn comp','Mean deg',
             'edge_density','transitivity','diameter',
             'diameter_wo_weights','mean_distance','assortativity_degree') 
  cvals <- c(n_connected_components,size_connected_components,mean_deg,
             edge_density,transitivity,diameter,
             diameter_wo_weights,mean_distance,assortativity_degree)
  cdf <- data.frame(Feature=cfeat,Value = cvals)
  
  dev.off()
  write.table(cdf,file=paste("Network_building/Network_analysis/Topology/",output_filename,k,"_topology.txt",sep=""),sep="\t",row.names=F, quote=FALSE)
  k <- k+1
}

# ### ANALYZE CLIQUES
# cliq <- readRDS(file="Network_building/Network_analysis/Topology/Cliques/DSN2.rds")
# cliq

