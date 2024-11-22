library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(ggsci)
options(future.globals.maxSize = 60000 * 1024^2)

# 1.2 load 10X genomics datasets
samplename = c(
  "/human_Liver/2.process_data/N001/outs/filtered_feature_bc_matrix",
  "/human_Liver/2.process_data/T001/outs/filtered_feature_bc_matrix",
  "/human_Liver/2.process_data/N0111/outs/filtered_feature_bc_matrix",
  "/human_Liver/2.process_data/T0111/outs/filtered_feature_bc_matrix",
  "/human_Liver/2.process_data/N0118/outs/filtered_feature_bc_matrix",
  "/human_Liver/2.process_data/T0118/outs/filtered_feature_bc_matrix",
  "/human_Liver/2.process_data/N1118/outs/filtered_feature_bc_matrix",
  "/human_Liver/2.process_data/T1118/outs/filtered_feature_bc_matrix",
  "/human_Liver/2.process_data/Ct0422/outs/filtered_feature_bc_matrix",
  "/human_Liver/2.process_data/Tu0422/outs/filtered_feature_bc_matrix",
  "/human_Liver/2.process_data/GX40T1/outs/filtered_feature_bc_matrix",
  "/human_Liver/2.process_data/GX40T2/outs/filtered_feature_bc_matrix",
  "/human_Liver/2.process_data/P0525/outs/filtered_feature_bc_matrix",
  "/human_Liver/2.process_data/T0525/outs/filtered_feature_bc_matrix",
  "/human_Liver/2.process_data/NC095/outs/filtered_feature_bc_matrix",
  "/human_Liver/2.process_data/NC104/outs/filtered_feature_bc_matrix",
  "/human_Liver/2.process_data/NC106/outs/filtered_feature_bc_matrix",
  "/human_Liver/2.process_data/NC114/outs/filtered_feature_bc_matrix",
  "/human_Liver/2.process_data/NC119/outs/filtered_feature_bc_matrix",
  "/human_Liver/2.process_data/NC713/outs/filtered_feature_bc_matrix",
  "/human_Liver/2.process_data/NC725/outs/filtered_feature_bc_matrix",
  "/human_Liver/2.process_data/NC740/outs/filtered_feature_bc_matrix")

names(samplename) = c('N001','T001','N0111','T0111','N0118','T0118','N1118','T1118','Ct0422','Tu0422',
                      'GX40T1','GX40T2','P0525','T0525',
                      'T095','T104','T106','T114','T119','T713','T725','T740')

print(samplename)

# 1.3 read matrix
proj <- list()
for(i in 1:length(samplename)){
  print(names(samplename[i]))
  counts <- Read10X(data.dir = samplename[i])
  newproj <- CreateSeuratObject(counts = counts, min.cells = 3,min.features = 200, project = names(samplename[i]))
  #print(names(samplename[i]))
  #print(newproj) # print datasets
  newproj$sample <- names(samplename[i])
  newproj[["percent.mt"]] <- PercentageFeatureSet(object = newproj, pattern = "^MT-")
  newproj[["percent.hba"]] <- PercentageFeatureSet(object = newproj, pattern = "HBA-")
  newproj[["percent.hbb"]] <- PercentageFeatureSet(object = newproj, pattern = "HBB-")
  #saveRDS(newproj,file=paste0(names(samplename[i]),".rds"))
  proj[[names(samplename[i])]] <- newproj
}

# filter of each sample
pdf("0-feature.pdf")
VlnPlot(proj[[1]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.hba"), ncol = 3,pt.size=0)
VlnPlot(proj[[2]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.hba"), ncol = 3,pt.size=0)
VlnPlot(proj[[3]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.hba"), ncol = 3,pt.size=0)
VlnPlot(proj[[4]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.hba"), ncol = 3,pt.size=0)
VlnPlot(proj[[5]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.hba"), ncol = 3,pt.size=0)
VlnPlot(proj[[6]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.hba"), ncol = 3,pt.size=0)
VlnPlot(proj[[7]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.hba"), ncol = 3,pt.size=0)
VlnPlot(proj[[8]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.hba"), ncol = 3,pt.size=0)
VlnPlot(proj[[9]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.hba"), ncol = 3,pt.size=0)
VlnPlot(proj[[10]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.hba"), ncol = 3,pt.size=0)
VlnPlot(proj[[11]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.hba"), ncol = 3,pt.size=0)
VlnPlot(proj[[12]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.hba"), ncol = 3,pt.size=0)
VlnPlot(proj[[13]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.hba"), ncol = 3,pt.size=0)
VlnPlot(proj[[14]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.hba"), ncol = 3,pt.size=0)
VlnPlot(proj[[15]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.hba"), ncol = 3,pt.size=0)
VlnPlot(proj[[16]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.hba"), ncol = 3,pt.size=0)
VlnPlot(proj[[17]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.hba"), ncol = 3,pt.size=0)
VlnPlot(proj[[18]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.hba"), ncol = 3,pt.size=0)
VlnPlot(proj[[19]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.hba"), ncol = 3,pt.size=0)
VlnPlot(proj[[20]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.hba"), ncol = 3,pt.size=0)
VlnPlot(proj[[21]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.hba"), ncol = 3,pt.size=0)
VlnPlot(proj[[22]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.hba"), ncol = 3,pt.size=0)
dev.off()

# 1.4 merge
T0040 <- merge(proj[[11]],proj[[12]])
T0040$sample <- 'T0040'
T0040@meta.data$orig.ident <- 'T0040'
scc_integrated<-merge(proj[[1]],y=c(proj[[2]],proj[[3]],proj[[4]],proj[[5]],proj[[6]],proj[[7]],proj[[8]],proj[[9]],proj[[10]],
                                    T0040,proj[[13]],proj[[14]],proj[[15]],proj[[16]],proj[[17]],proj[[18]],proj[[19]],
                                    proj[[20]],proj[[21]],proj[[22]]))
Idents(scc_integrated) <- 'orig.ident'
print(dim(scc_integrated))

#1.5 plot data quality before
pdf("multi_sample_scRNA_before_filter_new.pdf",height = 15,width = 15)
VlnPlot(scc_integrated, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.hba","percent.hbb"), ncol = 1,pt.size = 0)
dev.off()

#1.6 plot data quality after
scc_integrated <- subset(x = scc_integrated, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & nCount_RNA > 500 & nCount_RNA < 150000 & percent.mt < 15)

pdf("multi_sample_scRNA_after_filter_new.pdf",height=10,width=10)
VlnPlot(scc_integrated, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 1,pt.size = 0)
dev.off()

print(dim(scc_integrated))

#1.7 sct-pc-umap

scc_integrated<-SCTransform(scc_integrated,return.only.var.genes = FALSE,assay = "RNA",verbose = FALSE)
scc_integrated <- RunPCA(object = scc_integrated,verbose = FALSE)
scc_integrated <- FindNeighbors(scc_integrated,dim=1:50)
scc_integrated <- FindClusters(scc_integrated,resolution = 0.8)
scc_integrated <- RunUMAP(scc_integrated,reduction="pca", dims = 1:50)
scc_integrated <- RunTSNE(scc_integrated,dims = 1:50)
print(scc_integrated)

#1.8

orgin <- scc_integrated$sample
orgin <- as.character(orgin)
orgin[scc_integrated$sample %in% c('N001','T001','N0111','T0111','N0118','T0118','N1118','T1118','Ct0422','Tu0422')] <- 'HBV+Fluke+'
orgin[scc_integrated$sample %in% c('T0040','P0525','T0525')] <- 'HBV-Fluke+'
orgin[scc_integrated$sample %in% c('T095','T104','T106','T114','T119','T713','T725','T740')] <- 'HBV+Fluke-'
table(orgin)

orgin <- factor(orgin,levels = c('HBV+Fluke+','HBV-Fluke+','HBV+Fluke-'))
scc_integrated$orgin <- orgin

type <- scc_integrated$sample
type <- as.character(type)
type[scc_integrated$sample %in% c('N001','T001','N0111','T0111','N0118','T0118','N1118','T1118','Ct0422','Tu0422')] <- 'Control'
type[scc_integrated$sample %in% c('T0040','P0525','T0525')] <- 'Fluke'
type[scc_integrated$sample %in% c('T095','T104','T106','T114','T119','T713','T725','T740')] <- 'HBV'
table(type)

type <- factor(type,levels = c("Control","Fluke","HBV"))
scc_integrated$type <- type

#  1.9 saveRDS and plot the Umap
saveRDS(scc_integrated,"multi_sample_scRNA_integrated_Umap_new.rds")

