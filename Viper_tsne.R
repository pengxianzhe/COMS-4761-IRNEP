########## viper analysis ##########
source("https://bioconductor.org/biocLite.R")
#biocLite("viper")
library(viper)

#### load the TF-target relation matrix
chea <- read.csv("gene_chea_matrix.txt",sep="\t",stringsAsFactors = F)
chea <- chea[-c(1,2),-c(2,3)]
rownames(chea) <- chea[,1]
chea <- chea[,-1]

## pancreas data
load("all_cells_clean.rda")
chea_sub <- chea[intersect(rownames(all.data),rownames(chea)),]
all.data_sub <- all.data[intersect(rownames(all.data),rownames(chea)),]
#write.table(chea,"gene_chea_regulation.txt",quote = F,sep="\t",col.names = NA)

## brain data
chea_sub <- chea[intersect(rownames(ESM),rownames(chea)),]
ESM_sub <- ESM[intersect(rownames(ESM),rownames(chea)),]

### bulid the regulons
temp <- sapply(chea_sub, as.numeric)
rownames(temp)<- rownames(chea_sub)

reg <- apply(temp, 2, function(x) {
  tfmode <- x[x>0]
  list(tfmode=tfmode/abs(tfmode), likelihood=rep(1, length(tfmode)))
})

vpres <- viper(ESM_sub, reg, verbose = T)



#### tsne plot
library(ggplot2)
size = 1.2
width = 6
height = 4
theme_update(plot.title = element_text(hjust = 0.5))

pca_plot <- function(tpm=vpres,title=title,meta=meta){
  pr<-prcomp(t(tpm), scale = F, center = T)
  #eigs <- pr$sdev^2
  IDc <- c("KO2"="red", "WT" = "blue", "R200Q" = "green") 
  p = ggplot(as.data.frame(pr$x)) +
    geom_point(aes(x=PC1, y=PC2, color=meta$Genotype), size=0.2) +
    scale_color_manual(values=IDc, name = "Genotype") +
    
    labs(x="PC1", 
         y="PC2") +
    theme(legend.position = "right", legend.title = element_text(size=10), 
          axis.title = element_text(size=14),title = element_text(size=14),
          legend.text = element_text(size=10),legend.key.size = unit(.5, "cm"))
  ggsave(title, plot=p, width=width, height=height)
}
pca_plot(t=vpres,title="PCA for nes.pdf",meta=meta)


IDc <- c("KO2"="red", "WT" = "blue", "R200Q" = "green")

tsne_sep = list()
set.seed(42)
tsne_sep[['p50']] = Rtsne(t(vpres), theta=0,perplexity=30,max_iter=5000)



###############################

library(Seurat)
seu.data <- CreateSeuratObject(ESM, project = "single cell analysis")
seu.data <- AddMetaData(seu.data, metadata = ESM_meta)
seu.data <- NormalizeData(seu.data)
seu.data <- ScaleData(seu.data, model.use = "negbinom")


seu.data <- FindVariableGenes(seu.data)
seu.data <- RunPCA(seu.data, pcs.compute = 20, weight.by.var = F , do.print = F) 




###########################
seu.data@scale.data <- vpres
seu.data <- RunPCA(seu.data,pc.genes = rownames(vpres), pcs.compute = 20, weight.by.var = F , do.print = F) 

#elbow plot

PCAPlot(seu.data,1,2,pt.size=0.8,group.by = "Cell_type")

#tSNE
seu.data <- RunTSNE(seu.data, dims.use = 1:10, do.fast = F,perplexity=30)

pdf("tsne.pdf")
TSNEPlot(seu.data, group.by = "Cell_type", pt.size = 1.0)
dev.off()

pdf("elbowplot.pdf")
PCElbowPlot(seu.data)
dev.off()

####################################

#### additional analysis on multiple networks
library(org.Hs.eg.db)
library(aracne.networks)
hs <- org.Hs.eg.db

#temp <- sapply(strsplit(row.names(ESM), split='_', fixed=T), function(x) (x[1]))
id <- select(hs, keys = rownames(seu.data.sub.0),columns = c("ENTREZID", "SYMBOL"),keytype = "SYMBOL")
id <- id[-which(duplicated(id$SYMBOL)),]
id <- id[complete.cases(id),]

tmp <- seu.data.sub@scale.data
tmp <- tmp[id$SYMBOL,]
rownames(tmp)<- id$ENTREZID

vpres <- viper(tmp, regulon = regulonpaad, verbose = T)
meta_vpres <- viper(tmp, regulon = c(regulongbm,regulonbrca), verbose = T)




