####prework####
#library and input dataset
library(Seurat) # Tools for Single cell Genomics; version 4.01
library(patchwork)# composer of plots; version 1.1.1
library(Matrix) # sparse and dense matrix classes and methods; version 1.2-18
library(tidyverse) # a collection of dataset manipulation tools; version 1.3.1.9000
library(clustree) #visualise clustering at different resolutions; verison 0.4.3 
library(clusterProfiler) #analyze and visualize functional profiles (GO and KEGG) of gene and gene clusters; version 3.18.1
library(org.Mm.eg.db) # genome wide annotation for mouse; version 3.12.0
library(presto) #Fast functions for differential expression using wicox and auc version 1.0.0
library(msigdbr) # MSigDB gene sets; version 7.4.1
library(fgsea) # fast gene set enrichment analysis; version 1.16.0
library(ggnewscale) #multiple fill and color scales in "ggplot2"; version 0.4.5
setwd("D:/gp/sc/")
meta <- read.csv("pancreas_FACS_metadata.csv")
anno <- read.csv("pancreas_annotations_FACS.csv",row.names = 1)
dt <- read.csv("Pancreas-counts.csv",row.names = 1)
# use the infomation in sample infomation and annotations to make meta data for next part
cellname <- colnames(dt)
sex <- c()
subtissue <- c()
# catch the sex and subtissue infomation
for (i in 1:length(cellname)){
  if(unlist(strsplit(cellname[i], "[.]"))[2]=="MAA000574"){
    sex <- c(sex,"M")
    subtissue <- c(subtissue,"Exocrine")
  } else if (unlist(strsplit(cellname[i], "[.]"))[2]== "MAA000577"){
    sex <- c(sex,"M")
    subtissue <- c(subtissue,"Endocrine")
  } else if (unlist(strsplit(cellname[i], "[.]"))[2]== "MAA000884"){
    sex <- c(sex,"M")
    subtissue <- c(subtissue,"Endocrine")
  } else if (unlist(strsplit(cellname[i], "[.]"))[2]== "MAA000910"){
    sex <- c(sex,"M")
    subtissue <- c(subtissue,"Exocrine")
  } else if (unlist(strsplit(cellname[i], "[.]"))[2]== "MAA001857"){
    sex <- c(sex,"F")
    subtissue <- c(subtissue,"Endocrine")
  } else if (unlist(strsplit(cellname[i], "[.]"))[2]== "MAA001861"){
    sex <- c(sex,"F")
    subtissue <- c(subtissue,"Exocrine")
  } else if (unlist(strsplit(cellname[i], "[.]"))[2]== "MAA001862"){
    sex <- c(sex,"F")
    subtissue <- c(subtissue,"Endocrine")
  } else if (unlist(strsplit(cellname[i], "[.]"))[2]== "MAA001862"){
    sex <- c(sex,"F")
    subtissue <- c(subtissue,"Exocrine")
  } else {
    sex <- c(sex,NA)
    subtissue <-c(subtissue,NA)
  }
}
mt<-data.frame(row.names = cellname,sex,subtissue)
mt$key = colnames(dt)
anno$key = rownames(anno)
# join the two given tables into one meta data
nmt<-full_join(mt, anno, by = "key")
rownames(nmt)<-nmt$key 
nmt<-select(nmt, -tissue,-cell_ontology_term_iri,-key)
p<-CreateSeuratObject(counts = dt,meta.data = nmt,
                      project = "pancreas", min.cells = 3, min.features = 200)
#view the object created
head(p@meta.data)
summary(p@meta.data)#1872 cells

# Someone suggested calculating the percentage of red blood cells, but I did not do it
#HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
#HB_m <- match(HB.genes, rownames(p@assays$RNA)) 
#HB.genes <- rownames(p@assays$RNA)[HB_m] 
#HB.genes <- HB.genes[!is.na(HB.genes)] 
#scRNA[["percent.HB"]]<-PercentageFeatureSet(scRNA, features=HB.genes)

# quality control by looking at the percentage of mitochonria genes
p[["percent.mt"]] = PercentageFeatureSet(p, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(p, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
# All cells have the same value of percent.mt.
FeatureScatter(p, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# Some cells are highly expressed?
# we found that nFeature_RNA should be within [200,6500] 
summary(p[["nFeature_RNA"]])
summary(p[["nCount_RNA"]])
p = subset(p, subset = nFeature_RNA > 200 & nFeature_RNA < 6500 ) # filter the outlier
# save the work here (not neccessary)
# saveRDS(p, file = "rds/qc")
#p = read_rds("rds/qc")
# normalise
hist(colSums(p$RNA@data),
     breaks = 100,
     main = "Total expression before normalisation",
     xlab = "Sum of expression")
# 3 ways to normalise, normalization results from RC (relative counts) looks better, 
# but cluster result from "LogNormalize" was better
p = NormalizeData(p, normalization.method = "LogNormalize",scale.factor = 10000)
#p = NormalizeData(p, normalization.method = "CLR", 
#                  scale.factor = 10000)
#p = NormalizeData(p, normalization.method = "RC",scale.factor = 10000)
hist(colSums(p[["RNA"]]@data),
     breaks = 100,
     main = "Total expression after normalisation",
     xlab = "Sum of expression")

#scale
scale.genes <-rownames(p)
p <- ScaleData(p, features = scale.genes)



# Identification of highly variable features
p <- FindVariableFeatures(p, selection.method = "vst", nfeatures = 2000)
#the 10 most highly variable genes
top10 = head(VariableFeatures(p), 10)
# plot variable features with labels
LabelPoints(plot = VariableFeaturePlot(p), points = top10, repel = TRUE)

#####dimensional reduction####

# before cluster
###Some cell cycle related genes can express differently at different stages, 
###so we first see whether they can affect clustering.
head(c(cc.genes$s.genes,cc.genes$g2m.genes))
# [1] "MCM5" "PCNA" "TYMS" "FEN1" "MCM2" "MCM4"
#match with the highly virable features
CaseMatch(c(cc.genes$s.genes,cc.genes$g2m.genes),VariableFeatures(p))
#genes associated with S phase
g2m_genes = cc.genes$g2m.genes
g2m_genes = CaseMatch(search = g2m_genes, match = rownames(p))
#genes associated with G2M phase
s_genes = cc.genes$s.genes
s_genes = CaseMatch(search = s_genes, match = rownames(p))
# add S.Score, G2M.Score, and Phase to the object meta data
p <- CellCycleScoring(object=p,  g2m.features=g2m_genes,  s.features=s_genes)
head(p@meta.data)
p <- RunPCA(p, features = c(s_genes, g2m_genes))
DimPlot(p, reduction = "pca", group.by = "Phase")
# we found that cell cycle did not affect clustering.
# If you want to remove, you can run:
# scRNAb <- ScaleData(p, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(p))

####cluster####
p <- RunPCA(p, features = VariableFeatures(p)) 
# view results
DimPlot(p, reduction = "pca") 
VizDimLoadings(p, dims = 1:2, reduction = "pca")
print(p[["pca"]], dims = 1:5, nfeatures = 5)
# heat map allows for nice visualization of sources of heterogeneity
DimHeatmap(p, dims = 1, cells = 500, balanced = TRUE)
# The positive and negative components are not so distinctive
DimHeatmap(p, dims = 1:15, cells = 500, balanced = TRUE)
# ElbowPlot: a ranking of principle components based on the percentage 
# of variance explained by each one 
# I used this to decide the number of principle components used later
#Instead of JackStrawPlot, ElbowPlot can be used to reduce computation time
ElbowPlot(p, ndims=25, reduction="pca") 
pc.num=1:20 # I chose the number 20 

# Cluster cells
# Computing SNN
p <- FindNeighbors(p, dims = pc.num) 
# Identify clusters of cells by a shared nearest neighbor (SNN) modularity optimization 
# based clustering algorithm.
# We should choose resolution and we run iteration to decide which one is good.
# resolution parameter: effect the number of clusters, 
# with increased values leading to a greater number of clusters.
pc <- FindClusters(
  object = p,
  resolution = c(seq(.4,1,.1))
)
clustree(pc@meta.data, prefix = "RNA_snn_res.")
# I think resolution of 0.5 was enough for clustering
p <- FindClusters(p, resolution = 0.5)
table(p@meta.data$seurat_clusters) #there were 12 cluters, 0-12
#save the work here, it is not necessary for you
saveRDS(p, file="rds/cluster2")
p <- read_rds("rds/cluster2")
####non-linear dimensional reduction
#use T-SNE  to do the dimension reduction
p = RunTSNE(p, dims = pc.num)
# first cooridinates of tsne
embed_tsne <- Embeddings(p, 'tsne')
plot1<-DimPlot(p, reduction = "tsne")
#UMAP
p <- RunUMAP(p, dims = pc.num)
embed_umap <- Embeddings(p, 'umap')
plot2<-DimPlot(p, reduction = "umap")
CombinePlots(plots = list(plot1, plot2),legend="bottom")

#In meta data, we have infomation about cell ontology, we can compare it with cluster we made
unique(p@meta.data$cell_ontology_class)
p1<-DimPlot(p, group.by="cell_ontology_class", label=T, label.size=5, reduction='umap')
p2<-DimPlot(p, label=T,reduction = "umap")
ptsen <- DimPlot(p, label=T,reduction = "tsne")
p1+p2
p1tsen<-p1<-DimPlot(p, group.by="cell_ontology_class", label=T, label.size=5, reduction='tsne')
p1tsen+ptsen
DimPlot(p,group.by="cell_ontology_class", label=T, label.size=5,reduction = "pca")

# Ins1, Ins2
FeaturePlot(p, features = c("Ins1"))
VlnPlot(p, features = c("Ins1"), slot = "counts", log = TRUE) #cluster1,3
#Ins1 is highly expressed in cluster 1 and 3
cluster1.markers = FindMarkers(p, ident.1 = 0, min.pct = 0.25)
cluster3.markers = FindMarkers(p, ident.1 = 0, min.pct = 0.25)
markergene1<-rownames(top_n(cluster1.markers ,n = 6, wt = avg_log2FC))
p1<-subset(p, idents="1")

# find all markers
#Finds markers (differentially expressed genes) for each of the identity classes
#only.pos:Only return positive markers 
#min.pct:only test genes that show a minimum difference 
#in the fraction of detection between the two groups.
all.markers <- FindAllMarkers(p, only.pos = TRUE, 
                              min.pct = 0.25, logfc.threshold = 0.25)
head(all.markers)
all.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
#we are interested in in ins1 and ins2 in cluster 1
select_genes <- c('Ins1','Ins2','Pyy','Ppy')
VlnPlot(p, features = select_genes, pt.size=0, ncol=2)
FeaturePlot(p, features = select_genes)

#save the work here, it is not necessary for you
#saveRDS(p, file="rds/markers")


# we find that cell type in cluster 2 are unknow.


##### compare female to male by enrichment analysis###
#p <- read_rds("rds/markers")
ps <- subset(p, cell_ontology_class == 'type B pancreatic cell')
dge.sex <- FindMarkers(ps, ident.1 = "F", ident.2 = "M", group.by = 'sex')
sig.dge.sex <- subset(dge.sex, p_val_adj<0.05)
sig.dge.sex %>% arrange(desc(avg_log2FC)) %>% head(10)


# library(psych)
# library(pheatmap)
# AverageExp=AverageExpression(ps)
# coorda<-corr.test(AverageExp$RNA,AverageExp$RNA,method="spearman")
# pheatmap(coorda$r)

# GO
ego_CC1 <- enrichGO(gene          = row.names(sig.dge.sex),
                    #universe     = row.names(dge.celltype),
                    OrgDb         = 'org.Mm.eg.db',
                    keyType       = 'SYMBOL',
                    ont           = "CC",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
ego_MF1 <- enrichGO(gene          = row.names(sig.dge.sex),
                    #universe     = row.names(dge.celltype),
                    OrgDb         = 'org.Mm.eg.db',
                    keyType       = 'SYMBOL',
                    ont           = "MF",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
ego_BP1 <- enrichGO(gene          = row.names(sig.dge.sex),
                    #universe     = row.names(dge.celltype),
                    OrgDb         = 'org.Mm.eg.db',
                    keyType       = 'SYMBOL',
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
BPdt <- data.frame(ego_BP1)
MFdt <- data.frame(ego_MF1)
CCdt <- data.frame(ego_CC1)
#substring the description for plotting
ego_CC1@result$Description <- substring(ego_CC1@result$Description,1,70)
ego_MF1@result$Description <- substring(ego_MF1@result$Description,1,70)
ego_BP1@result$Description <- substring(ego_BP1@result$Description,1,70)
p_BP1 <- barplot(ego_BP1,showCategory = 10) + ggtitle("barplot for Biological process")
p_CC1 <- barplot(ego_CC1,showCategory = 10) + ggtitle("barplot for Cellular component")
p_MF1 <- barplot(ego_MF1,showCategory = 10) + ggtitle("barplot for Molecular function")
plotc1 <- p_BP1/p_CC1/p_MF1
plotc1
# GO DAG graph
goplot(ego_BP1)
# plot Gene-Concept Network
cnetplot(ego_BP1, foldChange=genelist)
cnetplot(ego_CC1, foldChange=genelist)
cnetplot(ego_MF1, foldChange=genelist)
cnetplot(ego_BP1, foldChange=genelist, circular = TRUE, colorEdge = TRUE)

# KEGG
# change the format
genelist <- bitr(row.names(sig.dge.sex), fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Mm.eg.db') 
#7.06% of input gene IDs are fail to map...
genelist <- pull(genelist,ENTREZID)               
ekegg <- enrichKEGG(gene = genelist, organism = 'mouse')
ekegg.df=as.data.frame(ekegg)
barplot(ekegg, showCategory=20)
dotplot(ekegg, showCategory=20)
head(ekegg.df[,1:2],6)
# "mmu04917" "mmu04022" "mmu04935" "mmu04728" "mmu04910" "mmu04911"
# see more infomation on websites
browseKEGG(ekegg,'mmu04917') #Prolactin signaling pathway
browseKEGG(ekegg,"mmu04911") #Insulin secretion
browseKEGG(ekegg,"mmu04910") #Insulin signaling pathway


#fgsea
# get the gene set from MSigDB
msigdbr_species()#Mus musculus
mdb_c2 <- msigdbr(species = "Mus musculus", category = "C2")
fgsea_sets<- mdb_c2 %>% split(x = .$gene_symbol, f = .$gs_name)
# get ranks arranged by fold change
fsig.dge.sex = sig.dge.sex
fsig.dge.sex$genes = rownames(fsig.dge.sex)
fgenes = fsig.dge.sex %>% arrange(desc(avg_log2FC)) %>% dplyr::select(genes,avg_log2FC)
franks<- deframe(fgenes)
fgseaRes<- fgsea(fgsea_sets, stats = franks, nperm = 1000)
# see the top pathways
fgseaRes %>% as_tibble() %>% arrange(desc(NES)) %>% filter(pval < 0.05) %>% 
  head(n= 20)
# plot some pathways
plotEnrichment(fgsea_sets[["LU_AGING_BRAIN_UP"]],
               ranks) + labs(title="LU_AGING_BRAIN_UP_PATHWAY")
plotEnrichment(fgsea_sets[["KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY"]],
               ranks) + labs(title="KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY")
plotEnrichment(fgsea_sets[["KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY"]],
               ranks) + labs(title="KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY")
# show pathways with pval < 0.05 and top absolute NES value
ggplot(fgseaRes %>% as_tibble() %>% arrange(desc(NES)) %>% filter(pval < 0.05) %>% 
         head(n= 20), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES)) +
  coord_flip() +
  labs(x="Pathways", y="Normalized Enrichment Score",title="Gene sets NES from GSEA")

# poor result
fgseaRes[ES > 0][padj<0.05] # empty
fgseaRes[ES < 0][padj<0.05] # empty
#no plot can be drawn from here
topPathwaysUp <- fgseaRes[ES > 0][paj<0.05][head(order(NES,decreasing = T), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][padj<0.05][head(order(NES,decreasing = T), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(fgsea_sets[topPathways], franks, fgseaRes,
              gseaParam = 0.5)

