library(limma)
library(plyr)
library(edgeR)

master.data <- read.csv("HFN1a_counts.csv", sep=",", row.names=1, header = T, stringsAsFactors = F)
umi_count <- apply(master.data,2,sum)
gene_count <- apply(master.data, 2, function(c)sum(c!=0))
master.data <- rbind(master.data,umi_count,gene_count)
master.data <- master.data[,which(master.data[33695,]>5000 & master.data[33695,]<40000 & master.data[33696,]>2000 )]
all.data <- master.data[-c(33695,33696),]



### Filter Genes
# Import data and mapping info
#master.data <- read.csv("HFN1a_counts.csv", sep=",", row.names=1, header = T, stringsAsFactors = F)
gene.symbol <- read.csv("HFN1a_counts_downsample_addsymbol.csv", sep=",", row.names=1, header = T, stringsAsFactors = F)
master.data$symbol <- gene.symbol$symbol

all.data <- master.data
all.data <- all.data[-which(duplicated(all.data$symbol)),]
row.names(all.data) <- all.data$symbol
all.data <- all.data[,-ncol(all.data)]

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

cell_count <- table(sapply(colnames(all.data),function(x) substrRight(x,2)))

# Label each cell
lbl.data <- data.frame(t(c(rep("R200Q.1",2793),rep("KO2.1",3212),rep("WT.1",2115),
                           rep("WT.2",2146),rep("WT.3",2588),rep("WT.4",2489),
                           rep("KO2.2",1965),rep("KO2.3",2118),rep("R200Q.2",1983),
                           rep("R200Q.3",1948))), stringsAsFactors = F) # Batch
lbl.data <- rbind(lbl.data,data.frame(t(c(rep("R200Q",2793),rep("KO2",3212),
                                          rep("WT",2115),rep("WT",2146),rep("WT",2588),
                                          rep("WT",2489),rep("KO2",1965),rep("KO2",2118),
                                          rep("R200Q",1983),rep("R200Q",1948))), stringsAsFactors = F)) # Genotype
colnames(lbl.data) <- colnames(all.data) # to corresponding cell
row.names(lbl.data) <- c("Batch","Genotype")

 

# Remove if
#a. < 1,000 genes
#b. genes in less than <4 cells
#c. > 0.1 mitochondria

## index for each criteria

#b. genes in <10 cells
gene.count <- apply(all.data, 1, function(x) sum(x>0)) # number of cells with non-zero expression for each gene

pdf("gene_count_in_cells.pdf")
hist(log2(gene.count), breaks = 50, main="", freq = T, col = "grey", ylab="frequency", xlab="gene count")
dev.off()
keep.gene <- which(gene.count >= 33) # keep genes if expressed in >= 10 cells
all.data <- all.data[keep.gene,] # all cells expressing more than >1000 genes with non-zero expression, for all genes expressed in more than 10 cells.


#c. > 0.1 mitochondria
mito.idx <- integer()
via.idx <- numeric()
cell.via <- numeric()
mito.idx <- grep("MT-",rownames(all.data),value=F)
for (i in 1:ncol(all.data)) {
        via.idx <- (sum(all.data[mito.idx,i])/sum(all.data[,i]))
        cell.via <- c(cell.via, via.idx)
}
# plot viability histogram and filter data based on viability
pdf("percentage.pdf")
hist(cell.via, breaks = seq(0,.6,0.005), main="", freq = F, col = "grey", ylab="", xlab="percentage")
abline(v = 0.1, lty = 3) # <98% cuz the majority have <0.3 viability score
dev.off()
all.data <- all.data[-mito.idx,which(cell.via <= 0.1)]
lbl.data <- data.frame(lbl.data[,which(cell.via <= 0.1)], stringsAsFactors = F)

# Organized data
save(all.data,file="all_cells_clean.rda")
save(lbl.data,file="label_clean.rda")
