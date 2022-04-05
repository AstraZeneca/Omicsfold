#Build feature interaction network for visualization in cytoscape based on correlations between features selected in final OmicsFold model

rm(list = ls())

detach("package:OmicsFold", unload = TRUE)
library(ggplot2)

setwd("~/Documents/COCKTAIL/")

###read in diablo model generated in step 3 of OmicsFold
diablo.model<-readRDS("3_trained_model_object.rds")

selected.factors <- list()
mylist<-list()

#set number of blocks in model
block.count<-4


#######################
###edit find.feature.associations function to only include top 20 features per block

find.feature.associations<-function (diablo.tuned, block.count) 
{
  pdf(file = NULL)
  circos <- mixOmics::circosPlot(diablo.tuned, cutoff = 0.7, 
                                 line = TRUE, size.labels = 1.5)
  dev.off()
  selected.factors <- list()
  for (i in 1:block.count) {
    dev.new(width = 3000, height = 3000, unit = "px")
    loadings <- mixOmics::plotLoadings(diablo.model, block = i, 
                                       comp = 1, contrib = "max", method = "median", ndisplay = 20)
    dev.off()
    top.features<-length(rownames(loadings))
    list.top.features<-list(top.features)
    names(list.top.features)<-names(diablo.model$keepX)[i]
    mylist <- c(mylist,list.top.features)
    selected.factors <- c(selected.factors, rownames(loadings))
  }
  circos.selected <- circos[rownames(circos) %in% selected.factors, 
                            colnames(circos) %in% selected.factors]
  diag(circos.selected) <- 1
  circos.selected.long <- reshape2::melt(circos.selected)
  circos.dendro <- stats::as.dendrogram(hclust(d = dist(x = circos.selected)))
  dendro.plot <- ggdendro::ggdendrogram(data = circos.dendro, 
                                        rotate = TRUE)
  dendro.plot <- dendro.plot + ggplot2::theme(axis.text.y = element_blank())
  circos.order <- stats::order.dendrogram(circos.dendro)
  circos.selected.long$Var1 <- factor(x = circos.selected.long$Var1, 
                                      levels = rownames(circos.selected)[circos.order], ordered = TRUE)
  circos.selected.long$Var2 <- factor(x = circos.selected.long$Var2, 
                                      levels = rownames(circos.selected)[circos.order], ordered = TRUE)
  heatmap.plot <- ggplot2::ggplot(data = circos.selected.long, 
                                  ggplot2::aes(x = Var1, y = Var2)) + ggplot2::geom_tile(ggplot2::aes(fill = value)) + 
    viridis::scale_fill_viridis(option = "plasma") + ggplot2::theme(legend.position = "top")
  grid::grid.newpage()
  print(heatmap.plot, vp = grid::viewport(x = 0.4, y = 0.5, 
                                          width = 0.8, height = 1))
  print(dendro.plot, vp = grid::viewport(x = 0.9, y = 0.465, 
                                         width = 0.2, height = 0.99))
  return(circos.selected)
}

###generate associations file
associations <- find.feature.associations(diablo.model, block.count = 4)

###generate mylist file which is a list of lists (top 20 features per block)
dev.off()
for (i in 1:block.count) {
  dev.new(width = 3000, height = 3000)
  loadings <- mixOmics::plotLoadings(diablo.model, block = i, 
                                     comp = 1, contrib = "max", method = "median", ndisplay = 20)
  dev.off()
  top.features<-length(rownames(loadings))
  list.top.features<-list(top.features)
  names(list.top.features)<-names(diablo.model$keepX)[i]
  mylist <- c(mylist,list.top.features)
  selected.factors <- c(selected.factors, rownames(loadings))
}

###generate block names for top 20 features
block.association <- character(0)
for (block.name in names(mylist)) {
  print(block.name)
  num.features.in.block <- mylist[[block.name]][1]
  block.association <- append(block.association,
                              rep(block.name, num.features.in.block))
}


#generate block names for all features in model
feature_block_names <- data.frame()
for (i in 2:length(diablo.model$loadings)-1){
active.block <- diablo.model$loadings[[i]]
active.block.names <- data.frame(rownames(active.block), names(diablo.model$loadings[i]))
feature_block_names <- rbind(feature_block_names, active.block.names)
}


library(OmicsFold)
###export network file with 0.5 correlation cut off cutoff
network <- export.matrix.as.network(associations, filename = "network_0.5.csv", cutoff= 0.5,
                         block.association = block.association)

##examining interactions for selected features per block


#create list with features to set as a source node in interaction network

list<-c("Alistipes_putredinis", "Streptococcus_parasanguinis")

library(readr)
library(data.table)
library(dplyr)

#subset features of interest with correct blocks for target nodes
network_subset <- subset(network, feature.1 == list | feature.2 == list)

#convert data into format for visualization of network in cytoscape 

#conditionally switch values where the target node is in the first column
network_subset_switch <-network_subset
names(network_subset_switch)[1]<-"target_node"
names(network_subset_switch)[2]<-"source_node"
names(network_subset_switch)[4]<-"source_node_block"
#identify which rows need to be switched and do not have the correct source node
i <- which(!network_subset_switch$source_node %in% list)
#create vectors with new row values
new_source <- network_subset_switch$target_node[i]
new_target <- network_subset_switch$source_node[i]
#change row values for the rows that need to be switched
network_subset_switch$source_node[i] <- new_source
network_subset_switch$target_node[i] <- new_target
#merge in block names for target nodes
network_subset_merged <- merge(network_subset_switch, feature_block_names, by.x="target_node", by.y="rownames.active.block.", all.x = TRUE)
names(network_subset_merged)[5] <- "target_node_block"

#remove intrablock connections
'%!in%' <- function(x,y)!('%in%'(x,y))
#create primary key for each id
network_subset_merged$id <- rownames(network_subset_merged)
#identify intrablock connections to remove
network_to_remove <- network_subset_merged[network_subset_merged$source_node_block == network_subset_merged$target_node_block 
              & network_subset_merged$target_node %!in% list,]
#remove the unwanted intrablock connections
network_filtered <- network_subset_merged[network_subset_merged$id %!in% network_to_remove$id,]
#remove primary key
network_filtered$id <- NULL


write.csv(network_filtered,"network.final.csv")

