library(RcisTarget)
library(GENIE3)

## calculated the weighted co-expression network
exprMatr <- all.data

weightMat <- GENIE3(exprMatr)
# Use Extra-Trees (ET) method
# 7 randomly chosen candidate regulators at each node of a tree

weightMat <- GENIE3(exprMatr, treeMethod="ET", K=7, nTrees=50)

linkList <- getLinkList(weightMat, threshold=0.1)
linkList <- getLinkList(weightMat, reportMax=5)

### prune the regulons by enriched motif
# Select motif database to use (i.e. organism and distance around TSS)
library(RcisTarget.hg19.motifDatabases.20k)
data(hg19_10kbpAroundTss_motifRanking)
data(hg19_direct_motifAnnotation)
motifRankings <- hg19_10kbpAroundTss_motifRanking

# Motif enrichment analysis
motifEnrichmentTable_wGenes <- cisTarget(geneLists, motifRankings,
                                         motifAnnot_direct=hg19_direct_motifAnnotation)

# Calculate AUC
motifs_AUC <- calcAUC(geneLists, motifRankings)

# Select significant motifs, add TF annotation & format as table
motifEnrichmentTable <- addMotifAnnotation(motifs_AUC, motifAnnot_direct=hg19_direct_motifAnnotation)

# Identify significant genes for each motif 
motifEnrichmentTable_wGenes <- addSignificantGenes(motifEnrichmentTable, 
                                                   geneSets=geneLists,
                                                   rankings=motifRankings, 
                                                   nCores=1,
                                                   method="iCisTarget")


motifs_AUC <- calcAUC(geneLists, motifRankings, nCores=1)



### rebuild the regulatory network

signifMotifNames <- motifEnrichmentTable$motif[1:3]

incidenceMatrix <- getSignificantGenes(geneLists, 
                                       motifRankings,
                                       signifRankingNames=signifMotifNames,
                                       plotCurve=TRUE, maxRank=5000-20, 
                                       genesFormat="incidMatrix",
                                       method="aprox")$incidMatrix

library(reshape2)
edges <- melt(incidenceMatrix)
edges <- edges[which(edges[,3]==1),1:2]
colnames(edges) <- c("from","to")

library(visNetwork)
motifs <- unique(as.character(edges[,1]))
genes <- unique(as.character(edges[,2]))
nodes <- data.frame(id=c(motifs, genes),   
                    label=c(motifs, genes),    
                    title=c(motifs, genes), # tooltip 
                    shape=c(rep("diamond", length(motifs)), rep("elypse", length(genes))),
                    color=c(rep("purple", length(motifs)), rep("skyblue", length(genes))))
