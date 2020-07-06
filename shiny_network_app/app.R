#########################################################################################
# SHINY APP THAT PERFORMS COMMUNITY DETECTION ALGORITHMS TO FIND META-PATIENT'S COMMUNITIES
# 
# 
#
# Beatriz Urda García 2020
#########################################################################################

# shiny::runApp(system.file("shiny", package = "visNetwork"))

library(shiny) # runExample("01_hello")
library(igraph)
library(magrittr)
library(visNetwork)
library(data.table)
library(DT)

# To share the library
library(rsconnect)
# rsconnect::deployApp('/home/beatriz/Desktop/ANALYSIS/shiny_network_app')

# net_filename <- 'pairwise_union_spearman_distance_sDEGs_network.txt'
# df <- read.csv(net_filename,header=T, sep="\t",stringsAsFactors = F); dim(df)
set.seed(5)

ui <- fluidPage(
  titlePanel("Community detection in disease similarity networks"),
  sidebarLayout(
    sidebarPanel(
      radioButtons("net_choice", "Select network:",
                   choices=c("Disease similarity network (DSN)" = "pairwise_union_spearman_distance_sDEGs_network.txt",
                             "Stratified similarity network (SSN)" = "metap_dis_pairwise_union_spearman_distance_sDEGs_network.txt")),
      radioButtons("pos_neg_interactions", "Select interactions:",
                   choices=c("Positive" = "pos",
                             "Negative" = "neg")),
      sliderInput("corr_slider","Edge's weight (|Spearman's correlation|):",min=0, max=1, value=c(0,1)),
      radioButtons("comm_algorithm", "Select community detection algorithm:",
                   choices=c("Greedy modularity optimization algorithm" = "greedy",
                             "Random walks (walktrap)" = "rand_walks")),
      textOutput("net_choicePlot"),
      width = 3
    ),
    mainPanel(
      # textOutput("txt"),
      htmlOutput("txt"),
      visNetworkOutput("net_plot", height = "700px"),
      # plotOutput("net_plot"),
      br(), br(), br(),
      # tableOutput("net_table"),
      DTOutput("net_dt_table",width='90%'),
      align="center"
    )
  )
)

server <- function(input, output) {
  
  # output$txt <- renderText({
  #   paste(input$current_node_id,"-",input$net_plot_selected)
  # })
  
  output$txt <- renderUI({
    HTML(paste('<b>','Select node:','</b>'))
    # paste(input$current_node_id,"-",input$net_plot_selected)
  })
  
  output$net_plot <- renderVisNetwork({
    df <- read.csv(input$net_choice,header=T, sep="\t",stringsAsFactors = F)
    if(input$pos_neg_interactions == 'pos'){
      # Selecting only positive interactions
      df <- df[df$Distance >= 0,]
    }else{
      # Selecting only negative interactions
      df <- df[df$Distance < 0,]
      df$Distance <- abs(df$Distance)
    }
    df <- df[((df$Distance >= input$corr_slider[1]) & (df$Distance <= input$corr_slider[2])), ]
    graph <- graph_from_data_frame(df, directed=FALSE)
    # is_weighted(graph)
    E(graph)$weight <- df$Distance
    # is_weighted(graph)
    #### Detect the communities and add the community as a vertex attribute
    if (input$comm_algorithm == 'greedy'){
      # Using greedy optimization of modularity
      fc <- fastgreedy.community(graph)
      V(graph)$community <- fc$membership
    }else if(input$comm_algorithm == 'rand_walks'){
      # Using random walks
      fc <- cluster_walktrap(graph)
      V(graph)$community <- fc$membership #membership(fc)
    }
    
    # Visualize the communities
    nodes <- data.frame(id = V(graph)$name, title = V(graph)$name, group = V(graph)$community)
    nodes <- nodes[order(nodes$id, decreasing = F),]
    edges <- get.data.frame(graph, what="edges")[1:2]
    edges$color <- rep("grey",length(edges$from))
    edges$value <- df$Distance
    
    visNetwork(nodes, edges) %>%
      visExport() %>%
      visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%   # list(enabled=TRUE, selected="BreastCancer")
      visIgraphLayout() %>%
      visInteraction(multiselect = T) %>%
      visEvents(select = "function(nodes) {
            Shiny.onInputChange('current_node_id', nodes.nodes);
            ;}")
  })
  observeEvent(input$current_node_id, {
    visNetworkProxy("net_plot") %>%
      visGetSelectedNodes()
      # visGetNodes()
  })
  
  # PRO TABLE
  output$net_dt_table <- renderDT({
    df <- read.csv(input$net_choice,header=T, sep="\t",stringsAsFactors = F)
    if(input$pos_neg_interactions == 'pos'){
      # Selecting only positive interactions
      df <- df[df$Distance >= 0,]
    }else{
      # Selecting only negative interactions
      df <- df[df$Distance < 0,]
    }
    df <- df[((abs(df$Distance) >= input$corr_slider[1]) & (abs(df$Distance) <= input$corr_slider[2])), ]
    if(is.null(input$net_plot_selected) | input$net_plot_selected == ''){
      colnames(df) <- c("Disease 1", "Disease 2","Spearman's correlation", "p-value", "adj.p-value")
      df
    }else{
      filtered <- df[(df$Dis1 == input$net_plot_selected | df$Dis2 == input$net_plot_selected),]
      # filtered <- df[(df$Dis1 == input$networkid_selected | df$Dis2 == input$networkid_selected),]
      # df %>%
      # filter((Dis1 == input$dis_choice | Dis2 == input$dis_choice))
      colnames(filtered) <- c("Disease 1", "Disease 2","Spearman's correlation", "p-value", "adj.p-value")
      filtered
    }
  })
  
  # BASIC TABLE
  output$net_table <- renderTable({
     
    df <- read.csv(input$net_choice,header=T, sep="\t",stringsAsFactors = F)
    if(input$pos_neg_interactions == 'pos'){
      # Selecting only positive interactions
      df <- df[df$Distance >= 0,]
    }else{
      # Selecting only negative interactions
      df <- df[df$Distance < 0,]
    }
    filtered <- df[(df$Dis1 == input$net_plot_selected | df$Dis2 == input$net_plot_selected),]
    # filtered <- df[(df$Dis1 == input$networkid_selected | df$Dis2 == input$networkid_selected),]
      # df %>%
      # filter((Dis1 == input$dis_choice | Dis2 == input$dis_choice))
    filtered
  })
}
shinyApp(ui = ui, server = server)


