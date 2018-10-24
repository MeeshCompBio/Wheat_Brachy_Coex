# Should have made a funciton for all of these

#Set working directory
setwd("C:/Users/Jean-Michel/Dropbox/Work/Manuscripts/Collaborations/WheatBrachyCoex/")
#Load in the differentially expressed data with FDR less than 0.05
ATBT2 <- read.delim("ATBT2.txt", stringsAsFactors = FALSE)
ATBT4 <- read.delim("ATBT4.txt", stringsAsFactors = FALSE)
ATBT6 <- read.delim("ATBT6.txt", stringsAsFactors = FALSE)


#Load in the co-expression clusters for both of the networks
W2691 <- read.csv("W2691_clusters.csv", stringsAsFactors = FALSE)
Sr9b <- read.csv("Sr9b_clusters.csv", stringsAsFactors = FALSE)

#calculate the total number of genes per cluster
WClust_totals <- table(W2691$cluster)
SClust_totals <- table(Sr9b$cluster)


#subset clusters based off DE gene list
W2691_DE_Overlap_ATBT2 <- W2691[W2691$X %in% toupper(ATBT2$gene_id), ]
WClust_totals_DE_ATBT2 <- table(W2691_DE_Overlap_ATBT2$cluster)
#Pull out anything great than two genes in a cluster
WClust_totals_DE_ATBT2_GT2 <- W2691_DE_Overlap_ATBT2[W2691_DE_Overlap_ATBT2$cluster %in% names(which(WClust_totals_DE_ATBT2 >= 2)),]


DEinfo <- c()
#add some statistics to each gene
for (i in 1:length(WClust_totals_DE_ATBT2_GT2[,1])){
  #get total number of genes in cluster
  ClusterTotal <- WClust_totals[which(names(WClust_totals)==WClust_totals_DE_ATBT2_GT2$cluster[i])]
  ClusterDETotal <- WClust_totals_DE_ATBT2[which(names(WClust_totals_DE_ATBT2)==WClust_totals_DE_ATBT2_GT2$cluster[i])]
  TotNumGenes <- length(W2691[,1])
  TotNumDEGenes <- length(ATBT2[toupper(ATBT2$gene_id) %in% W2691$X, ][,1])
  GOinfo <- ATBT2$GO.IDs..Description..via.Interpro[which(WClust_totals_DE_ATBT2_GT2$X[i]==toupper(ATBT2$gene_id))]
  HYPGEO <- as.numeric(phyper((ClusterDETotal-1),TotNumDEGenes, TotNumGenes,ClusterTotal, lower.tail = FALSE))
  Network <- c("W2691")
  DEGeneList <- c("ATBT2")
  DEinfo <- rbind(DEinfo, c(ClusterTotal,ClusterDETotal, TotNumGenes, TotNumDEGenes, GOinfo,HYPGEO, Network, DEGeneList))
  
}

Cluster_DE_ATBT2_W2691_Summary <- cbind(WClust_totals_DE_ATBT2_GT2, DEinfo)
colnames(Cluster_DE_ATBT2_W2691_Summary) <- c("Gene", "ClusterID",
                                        "ClusterTotal","ClusterDETotal", "NumNetworkGenes","TotalDEGenes",
                                        "GOinfo","HYPGEO", "Network", "DEGeneList")

#Lets do this on one cluster per line format
Reduce_Cluster_DE_ATBT2_W2691_Summary <- c()
for (i in 1:length(unique(Cluster_DE_ATBT2_W2691_Summary$ClusterID))){
  temp <- which(Cluster_DE_ATBT2_W2691_Summary$ClusterID == unique(Cluster_DE_ATBT2_W2691_Summary$ClusterID)[i])
  genes <- Cluster_DE_ATBT2_W2691_Summary$Gene[temp]
  ClusterTotal <- WClust_totals[which(names(WClust_totals)==unique(Cluster_DE_ATBT2_W2691_Summary$ClusterID)[i])]
  ClusterDETotal <- WClust_totals_DE_ATBT2[which(names(WClust_totals_DE_ATBT2)==unique(Cluster_DE_ATBT2_W2691_Summary$ClusterID)[i])]
  TotNumGenes <- length(W2691[,1])
  TotNumDEGenes <- length(ATBT2[toupper(ATBT2$gene_id) %in% W2691$X, ][,1])
  GOinfo <- ATBT2$GO.IDs..Description..via.Interpro[toupper(ATBT2$gene_id) %in% genes]
  HYPGEO <- as.numeric(phyper((ClusterDETotal-1),TotNumDEGenes, TotNumGenes,ClusterTotal, lower.tail = FALSE))
  Network <- c("W2691")
  DEGeneList <- c("ATBT2")
  Reduce_Cluster_DE_ATBT2_W2691_Summary <- rbind(Reduce_Cluster_DE_ATBT2_W2691_Summary
                                                        , c(unique(Cluster_DE_ATBT2_W2691_Summary$ClusterID)[i],ClusterTotal,ClusterDETotal,
                                                            TotNumGenes, TotNumDEGenes, paste(GOinfo, sep=",", collapse = ""), HYPGEO,
                                                            Network,DEGeneList, paste(genes, sep=",", collapse = ",")))
}

colnames(Reduce_Cluster_DE_ATBT2_W2691_Summary) <- c("ClusterID",
                                                     "ClusterTotal","ClusterDETotal", "NumNetworkGenes","TotalDEGenes",
                                                     "GOinfo","HYPGEO", "Network", "DEGeneList","Genes")

#Repeat the same for Sr9b
#subset clusters based off DE gene list
Sr9b_DE_Overlap_ATBT2 <- Sr9b[Sr9b$X %in% toupper(ATBT2$gene_id), ]
SClust_totals_DE_ATBT2 <- table(Sr9b_DE_Overlap_ATBT2$cluster)
#Pull out anything great than two genes in a cluster
SClust_totals_DE_ATBT2_GT2 <- Sr9b_DE_Overlap_ATBT2[Sr9b_DE_Overlap_ATBT2$cluster %in% names(which(SClust_totals_DE_ATBT2 >= 2)),]


DEinfo <- c()
#add some statistics to each gene
for (i in 1:length(SClust_totals_DE_ATBT2_GT2[,1])){
  #get total number of genes in cluster
  ClusterTotal <- SClust_totals[which(names(SClust_totals)==SClust_totals_DE_ATBT2_GT2$cluster[i])]
  ClusterDETotal <- SClust_totals_DE_ATBT2[which(names(SClust_totals_DE_ATBT2)==SClust_totals_DE_ATBT2_GT2$cluster[i])]
  TotNumGenes <- length(Sr9b[,1])
  TotNumDEGenes <- length(ATBT2[toupper(ATBT2$gene_id) %in% Sr9b$X, ][,1])
  GOinfo <- ATBT2$GO.IDs..Description..via.Interpro[which(SClust_totals_DE_ATBT2_GT2$X[i]==toupper(ATBT2$gene_id))]
  HYPGEO <- phyper(ClusterDETotal,TotNumDEGenes, TotNumGenes,ClusterTotal, lower.tail = FALSE)
  Network <- c("Sr9b")
  DEGeneList <- c("ATBT2")
  DEinfo <- rbind(DEinfo, c(ClusterTotal,ClusterDETotal, TotNumGenes, TotNumDEGenes, GOinfo,HYPGEO, Network, DEGeneList))
  
}

Cluster_DE_ATBT2_Sr9b_Summary <- cbind(SClust_totals_DE_ATBT2_GT2, DEinfo)
colnames(Cluster_DE_ATBT2_Sr9b_Summary) <- c("Gene", "ClusterID",
                                              "ClusterTotal","ClusterDETotal", "NumNetworkGenes","TotalDEGenes",
                                              "GOinfo","HYPGEO", "Network", "DEGeneList")

#Lets do this on one cluster per line format
Reduce_Cluster_DE_ATBT2_Sr9b_Summary <- c()
for (i in 1:length(unique(Cluster_DE_ATBT2_Sr9b_Summary$ClusterID))){
  temp <- which(Cluster_DE_ATBT2_Sr9b_Summary$ClusterID == unique(Cluster_DE_ATBT2_Sr9b_Summary$ClusterID)[i])
  genes <- Cluster_DE_ATBT2_Sr9b_Summary$Gene[temp]
  ClusterTotal <- SClust_totals[which(names(SClust_totals)==unique(Cluster_DE_ATBT2_Sr9b_Summary$ClusterID)[i])]
  ClusterDETotal <- SClust_totals_DE_ATBT2[which(names(SClust_totals_DE_ATBT2)==unique(Cluster_DE_ATBT2_Sr9b_Summary$ClusterID)[i])]
  TotNumGenes <- length(Sr9b[,1])
  TotNumDEGenes <- length(ATBT2[toupper(ATBT2$gene_id) %in% Sr9b$X, ][,1])
  GOinfo <- ATBT2$GO.IDs..Description..via.Interpro[toupper(ATBT2$gene_id) %in% genes]
  HYPGEO <- as.numeric(phyper((ClusterDETotal-1),TotNumDEGenes, TotNumGenes,ClusterTotal, lower.tail = FALSE))
  Network <- c("Sr9b")
  DEGeneList <- c("ATBT2")
  Reduce_Cluster_DE_ATBT2_Sr9b_Summary <- rbind(Reduce_Cluster_DE_ATBT2_Sr9b_Summary
                                                 , c(unique(Cluster_DE_ATBT2_Sr9b_Summary$ClusterID)[i],ClusterTotal,ClusterDETotal,
                                                     TotNumGenes, TotNumDEGenes, paste(GOinfo, sep=",", collapse = ""), HYPGEO,
                                                     Network,DEGeneList, paste(genes, sep=",", collapse = ",")))
}

colnames(Reduce_Cluster_DE_ATBT2_Sr9b_Summary) <- c("ClusterID",
                                                     "ClusterTotal","ClusterDETotal", "NumNetworkGenes","TotalDEGenes",
                                                     "GOinfo","HYPGEO", "Network", "DEGeneList","Genes")


############################################
#ATBT4
############################################
#subset clusters based off DE gene list
W2691_DE_Overlap_ATBT4 <- W2691[W2691$X %in% toupper(ATBT4$gene_id), ]
WClust_totals_DE_ATBT4 <- table(W2691_DE_Overlap_ATBT4$cluster)
#Pull out anything great than two genes in a cluster
WClust_totals_DE_ATBT4_GT2 <- W2691_DE_Overlap_ATBT4[W2691_DE_Overlap_ATBT4$cluster %in% names(which(WClust_totals_DE_ATBT4 >= 2)),]


DEinfo <- c()
#add some statistics to each gene
for (i in 1:length(WClust_totals_DE_ATBT4_GT2[,1])){
  #get total number of genes in cluster
  ClusterTotal <- WClust_totals[which(names(WClust_totals)==WClust_totals_DE_ATBT4_GT2$cluster[i])]
  ClusterDETotal <- WClust_totals_DE_ATBT4[which(names(WClust_totals_DE_ATBT4)==WClust_totals_DE_ATBT4_GT2$cluster[i])]
  TotNumGenes <- length(W2691[,1])
  TotNumDEGenes <- length(ATBT4[toupper(ATBT4$gene_id) %in% W2691$X, ][,1])
  GOinfo <- ATBT4$GO.IDs..Description..via.Interpro[which(WClust_totals_DE_ATBT4_GT2$X[i]==toupper(ATBT4$gene_id))]
  HYPGEO <- phyper(ClusterDETotal,TotNumDEGenes, TotNumGenes,ClusterTotal, lower.tail = FALSE)
  Network <- c("W2691")
  DEGeneList <- c("ATBT4")
  DEinfo <- rbind(DEinfo, c(ClusterTotal,ClusterDETotal, TotNumGenes, TotNumDEGenes, GOinfo,HYPGEO, Network, DEGeneList))
  
}

Cluster_DE_ATBT4_W2691_Summary <- cbind(WClust_totals_DE_ATBT4_GT2, DEinfo)
colnames(Cluster_DE_ATBT4_W2691_Summary) <- c("Gene", "ClusterID",
                                              "ClusterTotal","ClusterDETotal", "NumNetworkGenes","TotalDEGenes",
                                              "GOinfo","HYPGEO", "Network", "DEGeneList")

#Lets do this on one cluster per line format
Reduce_Cluster_DE_ATBT4_W2691_Summary <- c()
for (i in 1:length(unique(Cluster_DE_ATBT4_W2691_Summary$ClusterID))){
  temp <- which(Cluster_DE_ATBT4_W2691_Summary$ClusterID == unique(Cluster_DE_ATBT4_W2691_Summary$ClusterID)[i])
  genes <- Cluster_DE_ATBT4_W2691_Summary$Gene[temp]
  ClusterTotal <- WClust_totals[which(names(WClust_totals)==unique(Cluster_DE_ATBT4_W2691_Summary$ClusterID)[i])]
  ClusterDETotal <- WClust_totals_DE_ATBT4[which(names(WClust_totals_DE_ATBT4)==unique(Cluster_DE_ATBT4_W2691_Summary$ClusterID)[i])]
  TotNumGenes <- length(W2691[,1])
  TotNumDEGenes <- length(ATBT4[toupper(ATBT4$gene_id) %in% W2691$X, ][,1])
  GOinfo <- ATBT4$GO.IDs..Description..via.Interpro[toupper(ATBT4$gene_id) %in% genes]
  HYPGEO <- as.numeric(phyper((ClusterDETotal-1),TotNumDEGenes, TotNumGenes,ClusterTotal, lower.tail = FALSE))
  Network <- c("W2691")
  DEGeneList <- c("ATBT4")
  Reduce_Cluster_DE_ATBT4_W2691_Summary <- rbind(Reduce_Cluster_DE_ATBT4_W2691_Summary
                                                 , c(unique(Cluster_DE_ATBT4_W2691_Summary$ClusterID)[i],ClusterTotal,ClusterDETotal,
                                                     TotNumGenes, TotNumDEGenes, paste(GOinfo, sep=",", collapse = ""), HYPGEO,
                                                     Network,DEGeneList, paste(genes, sep=",", collapse = ",")))
}

colnames(Reduce_Cluster_DE_ATBT4_W2691_Summary) <- c("ClusterID",
                                                     "ClusterTotal","ClusterDETotal", "NumNetworkGenes","TotalDEGenes",
                                                     "GOinfo","HYPGEO", "Network", "DEGeneList","Genes")





#Repeat the same for Sr9b
#subset clusters based off DE gene list
Sr9b_DE_Overlap_ATBT4 <- Sr9b[Sr9b$X %in% toupper(ATBT4$gene_id), ]
SClust_totals_DE_ATBT4 <- table(Sr9b_DE_Overlap_ATBT4$cluster)
#Pull out anything great than two genes in a cluster
SClust_totals_DE_ATBT4_GT2 <- Sr9b_DE_Overlap_ATBT4[Sr9b_DE_Overlap_ATBT4$cluster %in% names(which(SClust_totals_DE_ATBT4 >= 2)),]


DEinfo <- c()
#add some statistics to each gene
for (i in 1:length(SClust_totals_DE_ATBT4_GT2[,1])){
  #get total number of genes in cluster
  ClusterTotal <- SClust_totals[which(names(SClust_totals)==SClust_totals_DE_ATBT4_GT2$cluster[i])]
  ClusterDETotal <- SClust_totals_DE_ATBT4[which(names(SClust_totals_DE_ATBT4)==SClust_totals_DE_ATBT4_GT2$cluster[i])]
  TotNumGenes <- length(Sr9b[,1])
  TotNumDEGenes <- length(ATBT4[toupper(ATBT4$gene_id) %in% Sr9b$X, ][,1])
  GOinfo <- ATBT4$GO.IDs..Description..via.Interpro[which(SClust_totals_DE_ATBT4_GT2$X[i]==toupper(ATBT4$gene_id))]
  HYPGEO <- phyper(ClusterDETotal,TotNumDEGenes, TotNumGenes,ClusterTotal, lower.tail = FALSE)
  Network <- c("Sr9b")
  DEGeneList <- c("ATBT4")
  DEinfo <- rbind(DEinfo, c(ClusterTotal,ClusterDETotal, TotNumGenes, TotNumDEGenes, GOinfo,HYPGEO, Network, DEGeneList))
  
}

Cluster_DE_ATBT4_Sr9b_Summary <- cbind(SClust_totals_DE_ATBT4_GT2, DEinfo)
colnames(Cluster_DE_ATBT4_Sr9b_Summary) <- c("Gene", "ClusterID",
                                              "ClusterTotal","ClusterDETotal", "NumNetworkGenes","TotalDEGenes",
                                              "GOinfo","HYPGEO", "Network", "DEGeneList")
#Lets do this on one cluster per line format
Reduce_Cluster_DE_ATBT4_Sr9b_Summary <- c()
for (i in 1:length(unique(Cluster_DE_ATBT4_Sr9b_Summary$ClusterID))){
  temp <- which(Cluster_DE_ATBT4_Sr9b_Summary$ClusterID == unique(Cluster_DE_ATBT4_Sr9b_Summary$ClusterID)[i])
  genes <- Cluster_DE_ATBT4_Sr9b_Summary$Gene[temp]
  ClusterTotal <- SClust_totals[which(names(SClust_totals)==unique(Cluster_DE_ATBT4_Sr9b_Summary$ClusterID)[i])]
  ClusterDETotal <- SClust_totals_DE_ATBT4[which(names(SClust_totals_DE_ATBT4)==unique(Cluster_DE_ATBT4_Sr9b_Summary$ClusterID)[i])]
  TotNumGenes <- length(Sr9b[,1])
  TotNumDEGenes <- length(ATBT4[toupper(ATBT4$gene_id) %in% Sr9b$X, ][,1])
  GOinfo <- ATBT4$GO.IDs..Description..via.Interpro[toupper(ATBT4$gene_id) %in% genes]
  HYPGEO <- as.numeric(phyper((ClusterDETotal-1),TotNumDEGenes, TotNumGenes,ClusterTotal, lower.tail = FALSE))
  Network <- c("Sr9b")
  DEGeneList <- c("ATBT4")
  Reduce_Cluster_DE_ATBT4_Sr9b_Summary <- rbind(Reduce_Cluster_DE_ATBT4_Sr9b_Summary
                                                , c(unique(Cluster_DE_ATBT4_Sr9b_Summary$ClusterID)[i],ClusterTotal,ClusterDETotal,
                                                    TotNumGenes, TotNumDEGenes, paste(GOinfo, sep=",", collapse = ""), HYPGEO,
                                                    Network,DEGeneList, paste(genes, sep=",", collapse = ",")))
}

colnames(Cluster_DE_ATBT4_Sr9b_Summary) <- c("Gene", "ClusterID",
                                             "ClusterTotal","ClusterDETotal", "NumNetworkGenes","TotalDEGenes",
                                             "GOinfo","HYPGEO", "Network", "DEGeneList")
                                                    

############################################
#ATBT6
############################################
#subset clusters based off DE gene list
W2691_DE_Overlap_ATBT6 <- W2691[W2691$X %in% toupper(ATBT6$gene_id), ]
WClust_totals_DE_ATBT6 <- table(W2691_DE_Overlap_ATBT6$cluster)
#Pull out anything great than two genes in a cluster
WClust_totals_DE_ATBT6_GT2 <- W2691_DE_Overlap_ATBT6[W2691_DE_Overlap_ATBT6$cluster %in% names(which(WClust_totals_DE_ATBT6 >= 2)),]


DEinfo <- c()
#add some statistics to each gene
for (i in 1:length(WClust_totals_DE_ATBT6_GT2[,1])){
  #get total number of genes in cluster
  ClusterTotal <- WClust_totals[which(names(WClust_totals)==WClust_totals_DE_ATBT6_GT2$cluster[i])]
  ClusterDETotal <- WClust_totals_DE_ATBT6[which(names(WClust_totals_DE_ATBT6)==WClust_totals_DE_ATBT6_GT2$cluster[i])]
  TotNumGenes <- length(W2691[,1])
  TotNumDEGenes <- length(ATBT6[toupper(ATBT6$gene_id) %in% W2691$X, ][,1])
  GOinfo <- ATBT6$GO.IDs..Description..via.Interpro[which(WClust_totals_DE_ATBT6_GT2$X[i]==toupper(ATBT6$gene_id))]
  HYPGEO <- phyper(ClusterDETotal,TotNumDEGenes, TotNumGenes,ClusterTotal, lower.tail = FALSE)
  Network <- c("W2691")
  DEGeneList <- c("ATBT6")
  DEinfo <- rbind(DEinfo, c(ClusterTotal,ClusterDETotal, TotNumGenes, TotNumDEGenes, GOinfo,HYPGEO, Network, DEGeneList))
  
}

Cluster_DE_ATBT6_W2691_Summary <- cbind(WClust_totals_DE_ATBT6_GT2, DEinfo)
colnames(Cluster_DE_ATBT6_W2691_Summary) <- c("Gene", "ClusterID",
                                              "ClusterTotal","ClusterDETotal", "NumNetworkGenes","TotalDEGenes",
                                              "GOinfo","HYPGEO", "Network", "DEGeneList")


#Lets do this on one cluster per line format
Reduce_Cluster_DE_ATBT6_W2691_Summary <- c()
for (i in 1:length(unique(Cluster_DE_ATBT6_W2691_Summary$ClusterID))){
  temp <- which(Cluster_DE_ATBT6_W2691_Summary$ClusterID == unique(Cluster_DE_ATBT6_W2691_Summary$ClusterID)[i])
  genes <- Cluster_DE_ATBT6_W2691_Summary$Gene[temp]
  ClusterTotal <- WClust_totals[which(names(WClust_totals)==unique(Cluster_DE_ATBT6_W2691_Summary$ClusterID)[i])]
  ClusterDETotal <- WClust_totals_DE_ATBT6[which(names(WClust_totals_DE_ATBT6)==unique(Cluster_DE_ATBT6_W2691_Summary$ClusterID)[i])]
  TotNumGenes <- length(W2691[,1])
  TotNumDEGenes <- length(ATBT6[toupper(ATBT6$gene_id) %in% W2691$X, ][,1])
  GOinfo <- ATBT6$GO.IDs..Description..via.Interpro[toupper(ATBT6$gene_id) %in% genes]
  HYPGEO <- as.numeric(phyper((ClusterDETotal-1),TotNumDEGenes, TotNumGenes,ClusterTotal, lower.tail = FALSE))
  Network <- c("W2691")
  DEGeneList <- c("ATBT6")
  Reduce_Cluster_DE_ATBT6_W2691_Summary <- rbind(Reduce_Cluster_DE_ATBT6_W2691_Summary
                                                 , c(unique(Cluster_DE_ATBT6_W2691_Summary$ClusterID)[i],ClusterTotal,ClusterDETotal,
                                                     TotNumGenes, TotNumDEGenes, paste(GOinfo, sep=",", collapse = ""), HYPGEO,
                                                     Network,DEGeneList, paste(genes, sep=",", collapse = ",")))
}

colnames(Reduce_Cluster_DE_ATBT6_W2691_Summary) <- c("ClusterID",
                                                     "ClusterTotal","ClusterDETotal", "NumNetworkGenes","TotalDEGenes",
                                                     "GOinfo","HYPGEO", "Network", "DEGeneList","Genes")







#Repeat the same for Sr9b
#subset clusters based off DE gene list
Sr9b_DE_Overlap_ATBT6 <- Sr9b[Sr9b$X %in% toupper(ATBT6$gene_id), ]
SClust_totals_DE_ATBT6 <- table(Sr9b_DE_Overlap_ATBT6$cluster)
#Pull out anything great than two genes in a cluster
SClust_totals_DE_ATBT6_GT2 <- Sr9b_DE_Overlap_ATBT6[Sr9b_DE_Overlap_ATBT6$cluster %in% names(which(SClust_totals_DE_ATBT6 >= 2)),]


DEinfo <- c()
#add some statistics to each gene
for (i in 1:length(SClust_totals_DE_ATBT6_GT2[,1])){
  #get total number of genes in cluster
  ClusterTotal <- SClust_totals[which(names(SClust_totals)==SClust_totals_DE_ATBT6_GT2$cluster[i])]
  ClusterDETotal <- SClust_totals_DE_ATBT6[which(names(SClust_totals_DE_ATBT6)==SClust_totals_DE_ATBT6_GT2$cluster[i])]
  TotNumGenes <- length(Sr9b[,1])
  TotNumDEGenes <- length(ATBT6[toupper(ATBT6$gene_id) %in% Sr9b$X, ][,1])
  GOinfo <- ATBT6$GO.IDs..Description..via.Interpro[which(SClust_totals_DE_ATBT6_GT2$X[i]==toupper(ATBT6$gene_id))]
  HYPGEO <- phyper(ClusterDETotal,TotNumDEGenes, TotNumGenes,ClusterTotal, lower.tail = FALSE)
  Network <- c("Sr9b")
  DEGeneList <- c("ATBT6")
  DEinfo <- rbind(DEinfo, c(ClusterTotal,ClusterDETotal, TotNumGenes, TotNumDEGenes, GOinfo,HYPGEO, Network, DEGeneList))
  
}

Cluster_DE_ATBT6_Sr9b_Summary <- cbind(SClust_totals_DE_ATBT6_GT2, DEinfo)
colnames(Cluster_DE_ATBT6_Sr9b_Summary) <- c("Gene", "ClusterID",
                                              "ClusterTotal","ClusterDETotal", "NumNetworkGenes","TotalDEGenes",
                                              "GOinfo","HYPGEO", "Network", "DEGeneList")

#Lets do this on one cluster per line format
Reduce_Cluster_DE_ATBT6_Sr9b_Summary <- c()
for (i in 1:length(unique(Cluster_DE_ATBT6_Sr9b_Summary$ClusterID))){
  temp <- which(Cluster_DE_ATBT6_Sr9b_Summary$ClusterID == unique(Cluster_DE_ATBT6_Sr9b_Summary$ClusterID)[i])
  genes <- Cluster_DE_ATBT6_Sr9b_Summary$Gene[temp]
  ClusterTotal <- SClust_totals[which(names(SClust_totals)==unique(Cluster_DE_ATBT6_Sr9b_Summary$ClusterID)[i])]
  ClusterDETotal <- SClust_totals_DE_ATBT6[which(names(SClust_totals_DE_ATBT6)==unique(Cluster_DE_ATBT6_Sr9b_Summary$ClusterID)[i])]
  TotNumGenes <- length(Sr9b[,1])
  TotNumDEGenes <- length(ATBT6[toupper(ATBT6$gene_id) %in% Sr9b$X, ][,1])
  GOinfo <- ATBT6$GO.IDs..Description..via.Interpro[toupper(ATBT6$gene_id) %in% genes]
  HYPGEO <- as.numeric(phyper((ClusterDETotal-1),TotNumDEGenes, TotNumGenes,ClusterTotal, lower.tail = FALSE))
  Network <- c("Sr9b")
  DEGeneList <- c("ATBT6")
  Reduce_Cluster_DE_ATBT6_Sr9b_Summary <- rbind(Reduce_Cluster_DE_ATBT6_Sr9b_Summary
                                                , c(unique(Cluster_DE_ATBT6_Sr9b_Summary$ClusterID)[i],ClusterTotal,ClusterDETotal,
                                                    TotNumGenes, TotNumDEGenes, paste(GOinfo, sep=",", collapse = ""), HYPGEO,
                                                    Network,DEGeneList, paste(genes, sep=",", collapse = ",")))
}

colnames(Cluster_DE_ATBT6_Sr9b_Summary) <- c("Gene", "ClusterID",
                                             "ClusterTotal","ClusterDETotal", "NumNetworkGenes","TotalDEGenes",
                                             "GOinfo","HYPGEO", "Network", "DEGeneList")





FinalData <- rbind(Cluster_DE_ATBT2_Sr9b_Summary,Cluster_DE_ATBT2_W2691_Summary, Cluster_DE_ATBT4_Sr9b_Summary,
                   Cluster_DE_ATBT4_W2691_Summary, Cluster_DE_ATBT6_Sr9b_Summary, Cluster_DE_ATBT6_W2691_Summary)
write.table(FinalData, "DE_Coex_Wheat_Cluster_Enrichment_Results.txt", sep = "\t", quote = FALSE, row.names = FALSE)



FinalReduce <- rbind(Reduce_Cluster_DE_ATBT2_Sr9b_Summary, Reduce_Cluster_DE_ATBT2_W2691_Summary, Reduce_Cluster_DE_ATBT4_Sr9b_Summary,
                     Reduce_Cluster_DE_ATBT4_W2691_Summary, Reduce_Cluster_DE_ATBT6_Sr9b_Summary, Reduce_Cluster_DE_ATBT6_W2691_Summary)

FDR1 <- p.adjust(Reduce_Cluster_DE_ATBT2_Sr9b_Summary[,7]  , method = "fdr")
FDR2 <- p.adjust(Reduce_Cluster_DE_ATBT2_W2691_Summary[,7]  , method = "fdr")
FDR3 <- p.adjust(Reduce_Cluster_DE_ATBT4_Sr9b_Summary[,7]  , method = "fdr")
FDR4 <- p.adjust(Reduce_Cluster_DE_ATBT4_W2691_Summary[,7]  , method = "fdr")
FDR5 <- p.adjust(Reduce_Cluster_DE_ATBT6_Sr9b_Summary[,7]  , method = "fdr")
FDR6 <- p.adjust(Reduce_Cluster_DE_ATBT6_W2691_Summary[,7]  , method = "fdr")

FDR <- c(FDR1, FDR2, FDR3, FDR4, FDR5, FDR6)

FinalReduce <- cbind(FinalReduce,FDR)

FinalReduce <- cbind(FinalReduce)
write.table(FinalReduce, "DE_Coex_Wheat_Cluster_Enrichment_Results_Corrected_Reduced_FDR_V2.txt", sep = "\t", quote = FALSE, row.names = FALSE)












#Lets take all genes from all the DE lists into one analysis
UniqueDELIST <- rbind(ATBT2,ATBT4,ATBT6)
UniqueDELIST <- UniqueDELIST[!duplicated(UniqueDELIST$gene_id),]



############################################
#UniqueDELIST
############################################
#subset clusters based off DE gene list
W2691_DE_Overlap_UniqueDELIST <- W2691[W2691$X %in% toupper(UniqueDELIST$gene_id), ]
WClust_totals_DE_UniqueDELIST <- table(W2691_DE_Overlap_UniqueDELIST$cluster)
#Pull out anything great than two genes in a cluster
WClust_totals_DE_UniqueDELIST_GT2 <- W2691_DE_Overlap_UniqueDELIST[W2691_DE_Overlap_UniqueDELIST$cluster %in% names(which(WClust_totals_DE_UniqueDELIST >= 2)),]


DEinfo <- c()
#add some statistics to each gene
for (i in 1:length(WClust_totals_DE_UniqueDELIST_GT2[,1])){
  #get total number of genes in cluster
  ClusterTotal <- WClust_totals[which(names(WClust_totals)==WClust_totals_DE_UniqueDELIST_GT2$cluster[i])]
  ClusterDETotal <- WClust_totals_DE_UniqueDELIST[which(names(WClust_totals_DE_UniqueDELIST)==WClust_totals_DE_UniqueDELIST_GT2$cluster[i])]
  TotNumGenes <- length(W2691[,1])
  TotNumDEGenes <- length(UniqueDELIST[,1])
  GOinfo <- UniqueDELIST$GO.IDs..Description..via.Interpro[which(WClust_totals_DE_UniqueDELIST_GT2$X[i]==toupper(UniqueDELIST$gene_id))]
  HYPGEO <- phyper(ClusterDETotal,TotNumDEGenes, TotNumGenes,ClusterTotal, lower.tail = FALSE)
  Network <- c("W2691")
  DEGeneList <- c("UniqueDELIST")
  DEinfo <- rbind(DEinfo, c(ClusterTotal,ClusterDETotal, TotNumGenes, TotNumDEGenes, GOinfo,HYPGEO, Network, DEGeneList))
  
}

Cluster_DE_UniqueDELIST_W2691_Summary <- cbind(WClust_totals_DE_UniqueDELIST_GT2, DEinfo)
colnames(Cluster_DE_UniqueDELIST_W2691_Summary) <- c("Gene", "ClusterID",
                                              "ClusterTotal","ClusterDETotal", "NumNetworkGenes","TotalDEGenes",
                                              "GOinfo","HYPGEO", "Network", "DEGeneList")


#Repeat the same for Sr9b
#subset clusters based off DE gene list
Sr9b_DE_Overlap_UniqueDELIST <- Sr9b[Sr9b$X %in% toupper(UniqueDELIST$gene_id), ]
SClust_totals_DE_UniqueDELIST <- table(Sr9b_DE_Overlap_UniqueDELIST$cluster)
#Pull out anything great than two genes in a cluster
SClust_totals_DE_UniqueDELIST_GT2 <- Sr9b_DE_Overlap_UniqueDELIST[Sr9b_DE_Overlap_UniqueDELIST$cluster %in% names(which(SClust_totals_DE_UniqueDELIST >= 2)),]


DEinfo <- c()
#add some statistics to each gene
for (i in 1:length(SClust_totals_DE_UniqueDELIST_GT2[,1])){
  #get total number of genes in cluster
  ClusterTotal <- SClust_totals[which(names(SClust_totals)==SClust_totals_DE_UniqueDELIST_GT2$cluster[i])]
  ClusterDETotal <- SClust_totals_DE_UniqueDELIST[which(names(SClust_totals_DE_UniqueDELIST)==SClust_totals_DE_UniqueDELIST_GT2$cluster[i])]
  TotNumGenes <- length(Sr9b[,1])
  TotNumDEGenes <- length(UniqueDELIST[,1])
  GOinfo <- UniqueDELIST$GO.IDs..Description..via.Interpro[which(SClust_totals_DE_UniqueDELIST_GT2$X[i]==toupper(UniqueDELIST$gene_id))]
  HYPGEO <- phyper(ClusterDETotal,TotNumDEGenes, TotNumGenes,ClusterTotal, lower.tail = FALSE)
  Network <- c("Sr9b")
  DEGeneList <- c("UniqueDELIST")
  DEinfo <- rbind(DEinfo, c(ClusterTotal,ClusterDETotal, TotNumGenes, TotNumDEGenes, GOinfo,HYPGEO, Network, DEGeneList))
  
}

Cluster_DE_UniqueDELIST_Sr9b_Summary <- cbind(SClust_totals_DE_UniqueDELIST_GT2, DEinfo)
colnames(Cluster_DE_UniqueDELIST_Sr9b_Summary) <- c("Gene", "ClusterID",
                                             "ClusterTotal","ClusterDETotal", "NumNetworkGenes","TotalDEGenes",
                                             "GOinfo","HYPGEO", "Network", "DEGeneList")


#####Brachy Data

setwd("C:/Users/Jean-Michel/Dropbox/Work/Manuscripts/Collaborations/WheatBrachyCoex/")
CMCT2 <- read.csv("CMCT2_DE.csv", stringsAsFactors = FALSE)
CMCT4 <- read.csv("CMCT4_DE.csv", stringsAsFactors = FALSE)
CMCT6 <- read.csv("CMCT6_DE.csv", stringsAsFactors = FALSE)

Bd21 <- read.csv("Bd21_NGFF_Clusters.csv", stringsAsFactors = FALSE)

WClust_totals <- table(Bd21$cluster)
#subset clusters based off DE gene list
Bd21_DE_Overlap_CMCT2 <- Bd21[Bd21$X %in% toupper(CMCT2$gene), ]
WClust_totals_DE_CMCT2 <- table(Bd21_DE_Overlap_CMCT2$cluster)
#Pull out anything great than two genes in a cluster
WClust_totals_DE_CMCT2_GT2 <- Bd21_DE_Overlap_CMCT2[Bd21_DE_Overlap_CMCT2$cluster %in% names(which(WClust_totals_DE_CMCT2 >= 2)),]


DEinfo <- c()
#add some statistics to each gene
for (i in 1:length(WClust_totals_DE_CMCT2_GT2[,1])){
  #get total number of genes in cluster
  ClusterTotal <- WClust_totals[which(names(WClust_totals)==WClust_totals_DE_CMCT2_GT2$cluster[i])]
  ClusterDETotal <- WClust_totals_DE_CMCT2[which(names(WClust_totals_DE_CMCT2)==WClust_totals_DE_CMCT2_GT2$cluster[i])]
  TotNumGenes <- length(Bd21[,1])
  TotNumDEGenes <- length(CMCT2[toupper(CMCT2$gene) %in% Bd21$X, ][,1])
  HYPGEO <- as.numeric(phyper((ClusterDETotal-1),TotNumDEGenes, TotNumGenes,ClusterTotal, lower.tail = FALSE))
  Network <- c("Bd21")
  DEGeneList <- c("CMCT2")
  DEinfo <- rbind(DEinfo, c(ClusterTotal,ClusterDETotal, TotNumGenes, TotNumDEGenes,HYPGEO, Network, DEGeneList))
  
}

Cluster_DE_CMCT2_Bd21_Summary <- cbind(WClust_totals_DE_CMCT2_GT2, DEinfo)
colnames(Cluster_DE_CMCT2_Bd21_Summary) <- c("Gene", "ClusterID",
                                              "ClusterTotal","ClusterDETotal", "NumNetworkGenes","TotalDEGenes","HYPGEO", "Network", "DEGeneList")





#Lets do this on one cluster per line format
Reduce_Cluster_DE_CMCT2_Bd21_Summary <- c()
for (i in 1:length(unique(Cluster_DE_CMCT2_Bd21_Summary$ClusterID))){
  temp <- which(Cluster_DE_CMCT2_Bd21_Summary$ClusterID == unique(Cluster_DE_CMCT2_Bd21_Summary$ClusterID)[i])
  genes <- Cluster_DE_CMCT2_Bd21_Summary$Gene[temp]
  ClusterTotal <- WClust_totals[which(names(WClust_totals)==unique(Cluster_DE_CMCT2_Bd21_Summary$ClusterID)[i])]
  ClusterDETotal <- WClust_totals_DE_CMCT2[which(names(WClust_totals_DE_CMCT2)==unique(Cluster_DE_CMCT2_Bd21_Summary$ClusterID)[i])]
  TotNumGenes <- length(Bd21[,1])
  TotNumDEGenes <- length(CMCT2[toupper(CMCT2$gene) %in% Bd21$X, ][,1])
  HYPGEO <- as.numeric(phyper((ClusterDETotal-1),TotNumDEGenes, TotNumGenes,ClusterTotal, lower.tail = FALSE))
  Network <- c("Bd21")
  DEGeneList <- c("CMCT2")
  Reduce_Cluster_DE_CMCT2_Bd21_Summary <- rbind(Reduce_Cluster_DE_CMCT2_Bd21_Summary
                                                 , c(unique(Cluster_DE_CMCT2_Bd21_Summary$ClusterID)[i],ClusterTotal,ClusterDETotal,
                                                     TotNumGenes, TotNumDEGenes, HYPGEO,
                                                     Network,DEGeneList, paste(genes, sep=",", collapse = ",")))
}

colnames(Reduce_Cluster_DE_CMCT2_Bd21_Summary) <- c("ClusterID",
                                                     "ClusterTotal","ClusterDETotal", "NumNetworkGenes","TotalDEGenes","HYPGEO", "Network", "DEGeneList","Genes")


####CMCT4
#subset clusters based off DE gene list
Bd21_DE_Overlap_CMCT4 <- Bd21[Bd21$X %in% toupper(CMCT4$gene), ]
WClust_totals_DE_CMCT4 <- table(Bd21_DE_Overlap_CMCT4$cluster)
#Pull out anything great than two genes in a cluster
WClust_totals_DE_CMCT4_GT2 <- Bd21_DE_Overlap_CMCT4[Bd21_DE_Overlap_CMCT4$cluster %in% names(which(WClust_totals_DE_CMCT4 >= 2)),]


DEinfo <- c()
#add some statistics to each gene
for (i in 1:length(WClust_totals_DE_CMCT4_GT2[,1])){
  #get total number of genes in cluster
  ClusterTotal <- WClust_totals[which(names(WClust_totals)==WClust_totals_DE_CMCT4_GT2$cluster[i])]
  ClusterDETotal <- WClust_totals_DE_CMCT4[which(names(WClust_totals_DE_CMCT4)==WClust_totals_DE_CMCT4_GT2$cluster[i])]
  TotNumGenes <- length(Bd21[,1])
  TotNumDEGenes <- length(CMCT4[toupper(CMCT4$gene) %in% Bd21$X, ][,1])
  HYPGEO <- as.numeric(phyper((ClusterDETotal-1),TotNumDEGenes, TotNumGenes,ClusterTotal, lower.tail = FALSE))
  Network <- c("Bd21")
  DEGeneList <- c("CMCT4")
  DEinfo <- rbind(DEinfo, c(ClusterTotal,ClusterDETotal, TotNumGenes, TotNumDEGenes,HYPGEO, Network, DEGeneList))
  
}

Cluster_DE_CMCT4_Bd21_Summary <- cbind(WClust_totals_DE_CMCT4_GT2, DEinfo)
colnames(Cluster_DE_CMCT4_Bd21_Summary) <- c("Gene", "ClusterID",
                                             "ClusterTotal","ClusterDETotal", "NumNetworkGenes","TotalDEGenes","HYPGEO", "Network", "DEGeneList")





#Lets do this on one cluster per line format
Reduce_Cluster_DE_CMCT4_Bd21_Summary <- c()
for (i in 1:length(unique(Cluster_DE_CMCT4_Bd21_Summary$ClusterID))){
  temp <- which(Cluster_DE_CMCT4_Bd21_Summary$ClusterID == unique(Cluster_DE_CMCT4_Bd21_Summary$ClusterID)[i])
  genes <- Cluster_DE_CMCT4_Bd21_Summary$Gene[temp]
  ClusterTotal <- WClust_totals[which(names(WClust_totals)==unique(Cluster_DE_CMCT4_Bd21_Summary$ClusterID)[i])]
  ClusterDETotal <- WClust_totals_DE_CMCT4[which(names(WClust_totals_DE_CMCT4)==unique(Cluster_DE_CMCT4_Bd21_Summary$ClusterID)[i])]
  TotNumGenes <- length(Bd21[,1])
  TotNumDEGenes <- length(CMCT4[toupper(CMCT4$gene) %in% Bd21$X, ][,1])
  HYPGEO <- as.numeric(phyper((ClusterDETotal-1),TotNumDEGenes, TotNumGenes,ClusterTotal, lower.tail = FALSE))
  Network <- c("Bd21")
  DEGeneList <- c("CMCT4")
  Reduce_Cluster_DE_CMCT4_Bd21_Summary <- rbind(Reduce_Cluster_DE_CMCT4_Bd21_Summary
                                                , c(unique(Cluster_DE_CMCT4_Bd21_Summary$ClusterID)[i],ClusterTotal,ClusterDETotal,
                                                    TotNumGenes, TotNumDEGenes, HYPGEO,
                                                    Network,DEGeneList, paste(genes, sep=",", collapse = ",")))
}

colnames(Reduce_Cluster_DE_CMCT4_Bd21_Summary) <- c("ClusterID",
                                                    "ClusterTotal","ClusterDETotal", "NumNetworkGenes","TotalDEGenes","HYPGEO", "Network", "DEGeneList","Genes")


######CMCT6
#subset clusters based off DE gene list
Bd21_DE_Overlap_CMCT6 <- Bd21[Bd21$X %in% toupper(CMCT6$gene), ]
WClust_totals_DE_CMCT6 <- table(Bd21_DE_Overlap_CMCT6$cluster)
#Pull out anything great than two genes in a cluster
WClust_totals_DE_CMCT6_GT2 <- Bd21_DE_Overlap_CMCT6[Bd21_DE_Overlap_CMCT6$cluster %in% names(which(WClust_totals_DE_CMCT6 >= 2)),]


DEinfo <- c()
#add some statistics to each gene
for (i in 1:length(WClust_totals_DE_CMCT6_GT2[,1])){
  #get total number of genes in cluster
  ClusterTotal <- WClust_totals[which(names(WClust_totals)==WClust_totals_DE_CMCT6_GT2$cluster[i])]
  ClusterDETotal <- WClust_totals_DE_CMCT6[which(names(WClust_totals_DE_CMCT6)==WClust_totals_DE_CMCT6_GT2$cluster[i])]
  TotNumGenes <- length(Bd21[,1])
  TotNumDEGenes <- length(CMCT6[toupper(CMCT6$gene) %in% Bd21$X, ][,1])
  HYPGEO <- as.numeric(phyper((ClusterDETotal-1),TotNumDEGenes, TotNumGenes,ClusterTotal, lower.tail = FALSE))
  Network <- c("Bd21")
  DEGeneList <- c("CMCT6")
  DEinfo <- rbind(DEinfo, c(ClusterTotal,ClusterDETotal, TotNumGenes, TotNumDEGenes,HYPGEO, Network, DEGeneList))
  
}

Cluster_DE_CMCT6_Bd21_Summary <- cbind(WClust_totals_DE_CMCT6_GT2, DEinfo)
colnames(Cluster_DE_CMCT6_Bd21_Summary) <- c("Gene", "ClusterID",
                                             "ClusterTotal","ClusterDETotal", "NumNetworkGenes","TotalDEGenes","HYPGEO", "Network", "DEGeneList")





#Lets do this on one cluster per line format
Reduce_Cluster_DE_CMCT6_Bd21_Summary <- c()
for (i in 1:length(unique(Cluster_DE_CMCT6_Bd21_Summary$ClusterID))){
  temp <- which(Cluster_DE_CMCT6_Bd21_Summary$ClusterID == unique(Cluster_DE_CMCT6_Bd21_Summary$ClusterID)[i])
  genes <- Cluster_DE_CMCT6_Bd21_Summary$Gene[temp]
  ClusterTotal <- WClust_totals[which(names(WClust_totals)==unique(Cluster_DE_CMCT6_Bd21_Summary$ClusterID)[i])]
  ClusterDETotal <- WClust_totals_DE_CMCT6[which(names(WClust_totals_DE_CMCT6)==unique(Cluster_DE_CMCT6_Bd21_Summary$ClusterID)[i])]
  TotNumGenes <- length(Bd21[,1])
  TotNumDEGenes <- length(CMCT6[toupper(CMCT6$gene) %in% Bd21$X, ][,1])
  HYPGEO <- as.numeric(phyper((ClusterDETotal-1),TotNumDEGenes, TotNumGenes,ClusterTotal, lower.tail = FALSE))
  Network <- c("Bd21")
  DEGeneList <- c("CMCT6")
  Reduce_Cluster_DE_CMCT6_Bd21_Summary <- rbind(Reduce_Cluster_DE_CMCT6_Bd21_Summary
                                                , c(unique(Cluster_DE_CMCT6_Bd21_Summary$ClusterID)[i],ClusterTotal,ClusterDETotal,
                                                    TotNumGenes, TotNumDEGenes, HYPGEO,
                                                    Network,DEGeneList, paste(genes, sep=",", collapse = ",")))
}

colnames(Reduce_Cluster_DE_CMCT6_Bd21_Summary) <- c("ClusterID",
                                                    "ClusterTotal","ClusterDETotal", "NumNetworkGenes","TotalDEGenes","HYPGEO", "Network", "DEGeneList","Genes")



FinalReduce <- rbind(Reduce_Cluster_DE_CMCT2_Bd21_Summary,Reduce_Cluster_DE_CMCT4_Bd21_Summary,Reduce_Cluster_DE_CMCT6_Bd21_Summary)

FDR1 <- p.adjust(Reduce_Cluster_DE_CMCT2_Bd21_Summary[,6]  , method = "fdr")
FDR2 <- p.adjust(Reduce_Cluster_DE_CMCT4_Bd21_Summary[,6]  , method = "fdr")
FDR3 <- p.adjust(Reduce_Cluster_DE_CMCT6_Bd21_Summary[,6]  , method = "fdr")


FDR <- c(FDR1, FDR2, FDR3)

FinalReduce <- cbind(FinalReduce,FDR)


write.table(FinalReduce, "DE_Coex_Brachy_Cluster_Enrichment_Results_Corrected_Reduced_FDR.txt", sep = "\t", quote = FALSE, row.names = FALSE)



