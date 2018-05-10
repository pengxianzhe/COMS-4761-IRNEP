library(RcisTarget)
library(GENIE3)

# calculated the weighted co-expression network
exprMatr <- all.data

weightMat <- GENIE3(exprMatr)
# Use Extra-Trees (ET) method
# 7 randomly chosen candidate regulators at each node of a tree
# 5 trees per ensemble
weightMat <- GENIE3(exprMatr, treeMethod="ET", K=7, nTrees=50)

linkList <- getLinkList(weightMat, threshold=0.1)
linkList <- getLinkList(weightMat, reportMax=5)


# Select motif database to use (i.e. organism and distance around TSS)
library(RcisTarget.hg19.motifDatabases.20k)
data(hg19_10kbpAroundTss_motifRanking)
data(hg19_direct_motifAnnotation)
motifRankings <- hg19_10kbpAroundTss_motifRanking

# Motif enrichment analysis:
motifEnrichmentTable_wGenes <- cisTarget(geneLists, motifRankings,
                                         motifAnnot_direct=hg19_direct_motifAnnotation)


# Calculate AUC
motifs_AUC <- calcAUC(geneLists, motifRankings)

# Select significant motifs, add TF annotation & format as table
motifEnrichmentTable <- addMotifAnnotation(motifs_AUC, motifAnnot_direct=hg19_direct_motifAnnotation)

# Identify significant genes for each motif (i.e. genes from the gene set in the top of the ranking)
# Note: Method 'iCisTarget' is more accurate, but slower
motifEnrichmentTable_wGenes <- addSignificantGenes(motifEnrichmentTable, 
                                                   geneSets=geneLists,
                                                   rankings=motifRankings, 
                                                   nCores=1,
                                                   method="iCisTarget")


motifs_AUC <- calcAUC(geneLists, motifRankings, nCores=1)
