###########################################################
##########MULTIOME ANALYSIS
##########################################################

#Load libraries 
library(JASPAR2022)
library(motifmatchr)
library(Seurat)
library(Signac)
library(monocle3)
library(SeuratWrappers)
library(Gviz)
library(GenomicRanges)
library(rtracklayer)
library(devtools)
library(ggplot2)
library(patchwork)
library(cicero)
library(dplyr)
library(tidyverse)
library(WriteXLS)
library(scCATCH)
library(SeuratDisk)
library(SeuratObject)
library(ggthemes)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)
library(GenomeInfoDb)
library(chromVAR)
library(BiocParallel)
library(Rsamtools)
library(GenomicRanges)
library(DropletQC)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(DESeq2)
library(ggrepel)
library(Rsamtools)
library(DropletQC)
library(ggpubr)
library(tibble)
library(ggridges)
library(spdep)

#Read in cellbender output files
sample1.data <- Read10X_h5(filename="cellbender/KO1/KO1_output_filtered_seurat.h5")
sample1 <- CreateSeuratObject(counts = sample1.data$'Gene Expression', assay = "RNA")
sample1$id <- "KO1"
sample1$treatment <- "KnockOut"

sample2.data <- Read10X_h5(filename="cellbender/KO2/KO2_output_filtered_seurat.h5")
sample2 <- CreateSeuratObject(counts = sample2.data$'Gene Expression', assay = "RNA")
sample2$id <- "KO2"
sample2$treatment <- "KnockOut"

sample3.data <- Read10X_h5(filename="cellbender/KO3/KO3_output_filtered_seurat.h5")
sample3 <- CreateSeuratObject(counts = sample3.data$'Gene Expression', assay = "RNA")
sample3$id <- "KO3"
sample3$treatment <- "KnockOut"

sample4.data <- Read10X_h5(filename="cellbender/KO4/KO4_output_filtered_seurat.h5")
sample4 <- CreateSeuratObject(counts = sample4.data$'Gene Expression', assay = "RNA")
sample4$id <- "KO4"
sample4$treatment <- "KnockOut"

sample5.data <- Read10X_h5(filename="cellbender/WT1/WT1_output_filtered_seurat.h5")
sample5 <- CreateSeuratObject(counts = sample5.data$'Gene Expression', assay = "RNA")
sample5$id <- "WT1"
sample5$treatment <- "WildType"

sample6.data <- Read10X_h5(filename="cellbender/WT2/WT2_output_filtered_seurat.h5")
sample6 <- CreateSeuratObject(counts = sample6.data$'Gene Expression', assay = "RNA")
sample6$id <- "WT2"
sample6$treatment <- "WildType"

sample7.data <- Read10X_h5(filename="cellbender/WT3/WT3_output_filtered_seurat.h5")
sample7 <- CreateSeuratObject(counts = sample7.data$'Gene Expression', assay = "RNA")
sample7$id <- "WT3"
sample7$treatment <- "WildType"

sample8.data <- Read10X_h5(filename="cellbender/WT4/WT4_output_filtered_seurat.h5")
sample8 <- CreateSeuratObject(counts = sample8.data$'Gene Expression', assay = "RNA")
sample8$id <- "WT4"
sample8$treatment <- "WildType"

#Pre-processing for the ATAC-Seq data
annotations_v79 <- GetGRangesFromEnsDb(ensdb=EnsDb.Mmusculus.v79)
genome(annotations_v79) <- "mm10"

#Change chromosome names to "chr#" or TSSEnrichment won't work.
chrStyle <- mapSeqlevels(seqlevels(annotations_v79), "UCSC")
annotations <- renameSeqlevels(annotations_v79, chrStyle)


#Create ATAC assays for each sample
fragpath <- "data/KO1/atac_fragments.tsv.gz"
sample1[["ATAC"]] <- CreateChromatinAssay(
  counts = sample1.data$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotations
)
fragpath <- "data/KO2/atac_fragments.tsv.gz"
sample2[["ATAC"]] <- CreateChromatinAssay(
  counts = sample2.data$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotations
)
fragpath <- "data/KO3/atac_fragments.tsv.gz"
sample3[["ATAC"]] <- CreateChromatinAssay(
  counts = sample3.data$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotations
)
fragpath <- "data/KO4/atac_fragments.tsv.gz"
sample4[["ATAC"]] <- CreateChromatinAssay(
  counts = sample4.data$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotations
)
fragpath <- "data/WT1/atac_fragments.tsv.gz"
sample5[["ATAC"]] <- CreateChromatinAssay(
  counts = sample5.data$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotations
)
fragpath <- "data/WT2/atac_fragments.tsv.gz"
sample6[["ATAC"]] <- CreateChromatinAssay(
  counts = sample6.data$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotations
)
fragpath <- "data/WT3/atac_fragments.tsv.gz"
sample7[["ATAC"]] <- CreateChromatinAssay(
  counts = sample7.data$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotations
)
fragpath <- "data/WT4/atac_fragments.tsv.gz"
sample8[["ATAC"]] <- CreateChromatinAssay(
  counts = sample8.data$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotations
)

#Merge Seurat Objects
scObj <- merge(sample1, y = c(sample2, sample3, sample4, sample5, sample6, sample7, sample8),
               add.cell.ids = c("sample1", "sample2", "sample3", "sample4", "sample5", "sample6", "sample7", "sample8"),
               project = "sc-Seurat")

#Add percent Mito
scObj[["percent.mt"]] <- PercentageFeatureSet(scObj, pattern = "^mt-")
scObj$mitoRatio <- PercentageFeatureSet(scObj, pattern = "^mt-")
scObj$mitoRatio <- scObj@meta.data$mitoRatio / 100

#Save raw data object
saveRDS(scObj, file = "data/RData/project-raw-cellbender.rds")

###############################################################
#######NUCLEAR FRACTION
###################################################

scObj <- readRDS("data/RData/project-raw-cellbender.rds")
metadata<-scObj@meta.data
metadata$cells <- rownames(metadata)
scObj@meta.data<-metadata

#Split objects into samples do this manually for each sample. below is an example for ko4
Idents(scObj)<-"id"

scObj_KO4 <-subset(scObj, id="KO4")
metadata<-scObj_KO4@meta.data
barcodes<-strsplit(scObj_KO4$cells,split="_")

k <- 2
x <- c()
for (item in barcodes) {
  x <- append(x, item[k])
}

nf <- nuclear_fraction_tags(
  bam = "data/KO4/gex_possorted_bam.bam",
  barcodes = x,
  tiles = 1,
  cores = 36,
  verbose = TRUE)

metadata$nuclear_fraction<-nf$nuclear_fraction
scObj_KO4@meta.data<-metadata

#once it is done for each sample, merge them
scObj_merged <- merge(scObj_KO1, y = c(scObj_KO2, scObj_KO3, scObj_KO4, scObj_WT1, scObj_WT2, scObj_WT3, scObj_WT4),project = "sc-Seurat")
metadata_merged<-scObj_merged@meta.data
write.table(metadata_merged,file="results/metadata_nf_cellbender.txt",sep="\t",row.names=TRUE,col.names=TRUE,quote=FALSE)

#Retain cellbender counts for filtered cells (from core). This will be used for furtehr analysis
scObj <- readRDS("data/RData/project-raw-cellbender.rds")
scObj_right <- readRDS("data/RData/project-raw.rds")

metadata1<-scObj@meta.data
metadata1$cells<-rownames(metadata1)
metadata2<-scObj_right@meta.data
metadata2$cells<-rownames(metadata2)

metadata<-merge(x=metadata1,y=metadata2,by="cells")
scObj_intersection<-subset(scObj,cells=metadata$cells)

saveRDS(scObj_intersection, file = "data/RData/project-raw-intersection.rds")

################################################################################
##################RNA ANALYSIS : FILTERING
##############################################################################

scObj <- readRDS("data/RData/project-raw-intersection.rds")
DefaultAssay(scObj) <- "RNA"
scObj@meta.data
# Add number of genes per UMI for each cell
scObj$RNA_log10GenesPerUMI <- log10(scObj$nFeature_RNA) / log10(scObj$nCount_RNA)

# look at the total cells, KO, and WT:
metadata <- scObj@meta.data
KO <- metadata %>% dplyr::filter(treatment=="KnockOut") %>% tally()
WT <- metadata %>% dplyr::filter(treatment=="WildType") %>% tally()
KO
WT

# See thresholds and make plots
p <- ggplot(metadata, aes(x=orig.ident, y=percent.mt,fill=orig.ident)) + 
  geom_violin() + 
  scale_fill_manual(values=c("steelblue"))+ 
  theme_classic()+ 
  xlab("All cells") +
  ylab("Mitochondrial percentage")+
  ggtitle("Filtering : Mitochondrial percentage")+
  theme(legend.position="none")+ geom_hline(yintercept = 2,linetype=2,color="darkred") + geom_hline(yintercept = 0,linetype=2,color="darkred")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.x= element_text(face="bold", size=12, colour="black"))+
  theme(title= element_text(face="bold", size=12, colour="black"))+
  theme(aspect.ratio=1)
p
png(file="filtering_mitochondrial_percentage.png",width=800,height=800,res=100)
print(p)
dev.off()

p <- ggplot(metadata, aes(x=orig.ident, y=nCount_RNA,fill=orig.ident)) + 
  geom_violin() + 
  scale_fill_manual(values=c("steelblue"))+ 
  theme_classic()+ 
  xlab("All cells") +
  ylab("Number of UMI")+
  ggtitle("Filtering : UMI threshold")+
  theme(legend.position="none")+ geom_hline(yintercept = 1000,linetype=2,color="darkred") + geom_hline(yintercept = 10000,linetype=2,color="darkred")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.x= element_text(face="bold", size=12, colour="black"))+
  theme(title= element_text(face="bold", size=12, colour="black"))+
  theme(aspect.ratio=1)+
  ylim(0,100000)
p
png(file="filtering_umi_threshold.png",width=800,height=800,res=100)
print(p)
dev.off()

p <- ggplot(metadata, aes(x=orig.ident, y=nFeature_RNA,fill=orig.ident)) + 
  geom_violin() + 
  scale_fill_manual(values=c("steelblue"))+ 
  theme_classic()+ 
  xlab("All cells") +
  ylab("Number of Genes")+
  ggtitle("Filtering : Gene threshold")+
  theme(legend.position="none")+ geom_hline(yintercept = 300,linetype=2,color="darkred") + geom_hline(yintercept = 2500,linetype=2,color="darkred")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.x= element_text(face="bold", size=12, colour="black"))+
  theme(title= element_text(face="bold", size=12, colour="black"))+
  theme(aspect.ratio=1)
p
png(file="filtering_gene_threshold.png",width=800,height=800,res=100)
print(p)
dev.off()

# Filter
scObj <- subset(scObj, subset = nCount_RNA >= 1000 &  nCount_RNA <10000 & nFeature_RNA >= 300 & nFeature_RNA <= 2500 & percent.mt <= 2)
metadata <- scObj@meta.data
KO <- metadata %>% dplyr::filter(treatment=="KnockOut") %>% tally()
WT <- metadata %>% dplyr::filter(treatment=="WildType") %>% tally()
KO
WT
stats <- metadata %>% group_by(id) %>% tally()

saveRDS(scObj, file = "data/RData/postcb_rna_filtered.rds")

scObj <- readRDS("data/RData/postcb_rna_filtered.rds")
metaa <- scObj@meta.data
metaa <- metaa %>%
  rename(seq_folder = orig.ident,
         nUMI = nCount_RNA,
         nGene = nFeature_RNA,
         Treatment=treatment)

p<- metaa %>%
  ggplot(aes(x=Treatment, fill=Treatment)) + 
  geom_bar() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("Number of Cells (After)") +
  theme(axis.text.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.text.x= element_text(face="bold", size=12, colour="black"))+
  theme(legend.text= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.x= element_text(face="bold", size=12, colour="black"))+
  theme(title= element_text(face="bold", size=12, colour="black"))+
  ylab("Count")+
  theme(aspect.ratio=1)+
  ylim(0,80000)
p
png(file="cells_treatment.png",width=800,height=800,res=100)
print(p)
dev.off()

p<- metaa %>%
  ggplot(aes(x=id, fill=id)) + 
  geom_bar() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("Number of Cells (After)") +
  theme(axis.text.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.text.x= element_text(face="bold", size=12, colour="black"))+
  theme(legend.text= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.x= element_text(face="bold", size=12, colour="black"))+
  theme(title= element_text(face="bold", size=12, colour="black"))+
  ylab("Count")+
  xlab("")+
  ylim(0,20000)+
  theme(aspect.ratio=1)
p
png(file="cells_id.png",width=800,height=800,res=100)
print(p)
dev.off()

Idents(scObj) <- "id"
p <- VlnPlot(scObj, features = "nCount_RNA", pt.size = 0)
p <- p + xlab("") + ylab("Number of UMIs") + ggtitle("Number of UMIs per sample (After)")+
  theme(axis.text.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.text.x= element_text(face="bold", size=12, colour="black"))+
  theme(legend.text= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.x= element_text(face="bold", size=12, colour="black"))+
  theme(title= element_text(face="bold", size=12, colour="black"))+
  theme(aspect.ratio=1)
p
png(file="umi_id.png",width=800,height=800,res=100)
print(p)
dev.off()

p <- VlnPlot(scObj, features = "nFeature_RNA", pt.size = 0)
p <- p + xlab("") + ylab("Number of genes") + ggtitle("Number of genes per sample (After)")+
  theme(axis.text.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.text.x= element_text(face="bold", size=12, colour="black"))+
  theme(legend.text= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.x= element_text(face="bold", size=12, colour="black"))+
  theme(title= element_text(face="bold", size=12, colour="black"))+
  theme(aspect.ratio=1)
p
png(file="genes_id.png",width=800,height=800,res=100)
print(p)
dev.off()

left<-scObj@meta.data
left$cells<-rownames(left)
right<-read.table(file="results/metadata_nf_cellbender.txt",sep="\t",header=TRUE)
nf<-right[,c(10,11)]
metadata<-merge(left,nf,by=c("cells"))
metadata<-metadata %>% mutate(UMI_bin = cut(nCount_RNA, breaks=c(0, 999,2000,3000,4000,5000,6000,7000,8000,9000,10000)))
p<-ggplot(metadata,aes(x=factor(UMI_bin),y=nuclear_fraction,fill=factor(UMI_bin)))+
  geom_violin()+ 
  theme_classic()+
  ylab("Intronic reads ratio") + 
  xlab("UMI bins")+
  theme(axis.text.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.text.x= element_text(face="bold", size=12, colour="black"))+
  theme(legend.text= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.x= element_text(face="bold", size=12, colour="black"))+
  theme(title= element_text(face="bold", size=12, colour="black"))  +
  theme(legend.title=element_blank())+
  ylim(0,1)+
  scale_x_discrete(labels=c('1000-2000','2000-3000','3000-4000','4000-5000','5000-6000','6000-7000','7000-8000','8000-9000','9000-10000' ))+
  scale_fill_discrete(labels=c('1000-2000','2000-3000','3000-4000','4000-5000','5000-6000','6000-7000','7000-8000','8000-9000','9000-10000'))
p
png(file="irr.png",width=800,height=800,res=100)
print(p)
dev.off()

##################################################################################################
##### RNA:  NORMALIZE THE DATA and INTEGRATE.  THIS CREATES THE SCT and INTEGRATED ASSAYS
##################################################################################################
scObj <- readRDS("data/RData/postcb_rna_filtered.rds")

DefaultAssay(scObj) <- "RNA"
scObj <- NormalizeData(object = scObj, assay = 'RNA', normalization.method = 'LogNormalize')
# split the object by id and normalize each.
scObj.objects <- SplitObject(scObj, split.by = "id")

for (i in 1:length(scObj.objects)) {
  scObj.objects[[i]] <- SCTransform(scObj.objects[[i]], method = "glmGamPoi", verbose = TRUE)
}

# Run SelectIntegrationFeatures to select the genes of interest.  If specific genes are
# required in the output, add them here.
features <- SelectIntegrationFeatures(object.list = scObj.objects, nfeatures = 3500)

# send the features to integration.
scObj.list <- PrepSCTIntegration(object.list = scObj.objects, anchor.features = features, verbose = TRUE)

# Save these variables so you don't have to recreate later.
save(list = c("scObj.list", "features"), file = "data/RData/postcb_rna_vars.Rdata")
load (file="data/RData/postcb_rna_vars.Rdata")

#Save data object
saveRDS(scObj, file = "data/RData/postcb_before_integration.rds")
scObj <- readRDS("data/RData/postcb_before_integration.rds")


reference_dataset <- which(names(scObj.list) %in% c("KO1", "WT1"))
scObj.anchors2 <- FindIntegrationAnchors(object.list = scObj.list, normalization.method = "SCT", 
                                         anchor.features = features, reference = reference_dataset)
save(scObj.anchors2, file="data/RData/postcb_rna_anchors.Rdata")
load (file="data/RData/postcb_rna_anchors.Rdata")

scObj <- IntegrateData(new.assay.name="integrated", anchorset = scObj.anchors2, normalization.method = "SCT", verbose=TRUE)
saveRDS(scObj, file = "data/RData/postcb_integrated.rds") 
scObj <- readRDS("data/RData/postcb_integrated.rds")

DefaultAssay(scObj) <- "integrated"
scObj <- RunPCA(scObj, verbose = FALSE, assay = "integrated",npcs=100)
scObj <- RunUMAP(scObj, dims = 1:30, verbose = FALSE, assay = "integrated")
scObj <- FindNeighbors(object = scObj, dims = 1:30, assay = "integrated")
scObj <- FindClusters(object = scObj, resolution = 0.3)
saveRDS(scObj, file = "data/RData/postcb_afterclustering.rds") 
scObj <- readRDS("data/RData/postcb_afterclustering.rds")

##Visualize UMAP for resolution of your choice
DefaultAssay(scObj)<-"RNA"
Idents(scObj) <- "integrated_snn_res.0.3" 
p<- DimPlot(scObj, reduction = "umap", raster=FALSE,label=TRUE,order=sort(as.integer(levels(scObj)),decreasing=TRUE))+theme(aspect.ratio=1)+coord_fixed(ratio=1)
p
png(file="UMAP_allcells.png",width=800,height=800,res=100)
print(p)
dev.off()
# split by treatment
p<- DimPlot(scObj, reduction="umap", label=TRUE, raster=FALSE,order=sort(as.integer(levels(scObj)),decreasing=TRUE), split.by="treatment")+theme(aspect.ratio=1)+coord_fixed(ratio=1)
p
png(file="UMAP_split.png",width=1600,height=800,res=100)
print(p)
dev.off()

f_cortex<-c("Slc17a7",	"Nrgn",	"Syt17",	"Cdh13",	"Adamts2",	"Stard8",	"Rrad",	"Cux2",	"Calb1",	"Rorb",	"Inhba",	"Arf5",	"Itga7",	"Whrn",	"Col27a1",	"Pld5",	"Lypd1",	"Gpr88",	"Tmem163",	"Arhgap25",	"Cpa6",	"Pcdh19",	"Satb2",	"Tcerg1I",	"Foxo1",	"Deptor",	"Pde1c",	"Nnat",	"Col6a1",	"Foxp2",	"Oprk1",	"Penk",	"Col27a1",	"Col23a1",	"Car3",	"Sulf1",	"Sla",	"Mgp",	"Car12",	"Rgs12",	"Gad1",	"Gad2",	"Pvalb",	"Sst",	"Vip",	"Reln",	"Cck",	"Cnr1",	"Npy",	"Ndnf",	"Nos1",	"Hepacam",	"Pdgfra",	"Vcan",	"Cspg4",	"Sox10",	"C1ql1",	"Mbp",	"Mobp",	"Mog",	"Gja1",	"F3",	"Aqp4",	"Gfap",	"Csf1r",	"Inpp5d",	"Cx3cr1",	"Apbb1ip",	"Itgam",	"Ctss","Flt1",	"Cldn5",	"Slc2a1",	"Fli1",	"Vtn",	"Nes",	"Prom1")
f_cortex<-unique(f_cortex)
p<-DotPlot(scObj, features=f_cortex, assay="RNA", dot.min=0.1, cols=c("gray95", "purple")) + theme(axis.text.x = element_text(angle = 90, hjust=1))
p
pdf(file="dotplot.pdf",height=10,width=15)
print(p)
dev.off()

######################################################
############SPLITTING CLUSTERS AND ANNOTATING
######################################################
scObj <- readRDS("data/RData/postcb_afterclustering.rds")
DefaultAssay(scObj) <- "integrated"
Idents(scObj) <- "integrated_snn_res.0.3"

DimPlot(scObj, reduction = "umap", label=TRUE, raster=FALSE)

scObj$split_clusters <- scObj$integrated_snn_res.0.3
cluster<-18
clus <- subset(scObj, idents=cluster)
clus <- FindNeighbors(object = clus, dims = 1:30, assay = "integrated")
res<-0.05
clus <- FindClusters(object = clus, resolution = res, assay = "integrated")
Idents(clus) <- paste("integrated_snn_res.",res,sep="")
DimPlot(clus, reduction = "umap", label=TRUE, raster=FALSE)

new_clus <- c(paste(cluster,"a",sep=""), paste(cluster,"b",sep=""))
levels(clus$integrated_snn_res.0.05) <- new_clus
clus_meta <- clus@meta.data
clusa <- clus_meta %>% dplyr::filter(integrated_snn_res.0.05==paste(cluster,"a",sep="")) 
clusb <- clus_meta %>% dplyr::filter(integrated_snn_res.0.05==paste(cluster,"b",sep="")) 

scObj_meta <- scObj@meta.data
scObj_meta <- scObj_meta %>% mutate(split_clusters = case_when(
  split_clusters=="18" & rownames(scObj_meta) %in% rownames(clusa) ~ "18a",
  split_clusters=="18" & rownames(scObj_meta) %in% rownames(clusb) ~ "18b",
  TRUE ~ as.character(split_clusters)
)) 
scObj_meta$split_clusters <- as.factor(scObj_meta$split_clusters)
scObj@meta.data <- scObj_meta

metadata<-scObj@meta.data
metadata$split_clusters<-as.character(metadata$split_clusters)
metadata$annotation[metadata$split_clusters=="0"]<-"layer 2/3"
metadata$annotation[metadata$split_clusters=="1"]<-"layer 4/5"
metadata$annotation[metadata$split_clusters=="2"]<-"layer 6"
metadata$annotation[metadata$split_clusters=="3"]<-"layer 5/6"
metadata$annotation[metadata$split_clusters=="4"]<-"layer 2/3/4 and olig"
metadata$annotation[metadata$split_clusters=="5"]<-"olig and astro"
metadata$annotation[metadata$split_clusters=="6"]<-"inh pvalb"
metadata$annotation[metadata$split_clusters=="7"]<-"layer 5a"
metadata$annotation[metadata$split_clusters=="8"]<-"layer 2"
metadata$annotation[metadata$split_clusters=="9"]<-"olig"
metadata$annotation[metadata$split_clusters=="10"]<-"inh vip reln cnr1"
metadata$annotation[metadata$split_clusters=="11"]<-"inh sst reln"
metadata$annotation[metadata$split_clusters=="12"]<-"layer 5"
metadata$annotation[metadata$split_clusters=="13"]<-"OPC"
metadata$annotation[metadata$split_clusters=="14"]<-"layer 5a/6"
metadata$annotation[metadata$split_clusters=="15"]<-"layer 6"
metadata$annotation[metadata$split_clusters=="16"]<-"layer 5/6"
metadata$annotation[metadata$split_clusters=="17"]<-"layer 5"
metadata$annotation[metadata$split_clusters=="18a"]<-"layer 6 "
metadata$annotation[metadata$split_clusters=="18b"]<-"micro"
metadata$annotation[metadata$split_clusters=="19"]<-"endo"
metadata$annotation[metadata$split_clusters=="20"]<-"layer 2/3 and 5 and 6"
metadata$annotation[metadata$split_clusters=="21"]<-"micro"

metadata$split_clusters[metadata$split_clusters=="18a"]<-"18"
metadata$split_clusters[metadata$split_clusters=="18b"]<-"22"

metadata$split_clusters <- as.factor(metadata$split_clusters)
scObj@meta.data <- metadata

Idents(scObj) <- "split_clusters"
levels(scObj$split_clusters)
levels(scObj)<-c("8","0","20","4","1","12","17","7","3","14","16","2","15","18","11","6","10","13","9","5","22","21","19")
DimPlot(scObj, reduction = "umap", label=TRUE, raster=FALSE)

f_cortex<-c("Slc17a7",	"Nrgn",	"Syt17",	"Cdh13",	"Adamts2",	"Stard8",	"Rrad",	"Cux2",	"Calb1",	"Rorb",	"Inhba",	"Arf5",	"Itga7",	"Whrn",	"Col27a1",	"Pld5",	"Lypd1",	"Gpr88",	"Tmem163",	"Arhgap25",	"Cpa6",	"Pcdh19",	"Satb2",	"Tcerg1",	"Foxo1",	"Deptor",	"Pde1c",	"Nnat",	"Col6a1",	"Foxp2",	"Oprk1",	"Penk",	"Col27a1",	"Col23a1",	"Car3",	"Sulf1",	"Sla",	"Mgp",	"Car12",	"Rgs12",	"Gad1",	"Gad2",	"Pvalb",	"Sst",	"Vip",	"Reln",	"Cck",	"Cnr1",	"Npy",	"Ndnf",	"Nos1",	"Hepacam",	"Pdgfra",	"Vcan",	"Cspg4",	"Sox10",	"C1ql1",	"Mbp",	"Mobp",	"Mog",	"Gja1",	"F3",	"Aqp4",	"Gfap",	"Csf1r",	"Inpp5d",	"Cx3cr1",	"Apbb1ip",	"Itgam",	"Ctss","Flt1",	"Cldn5",	"Slc2a1",	"Fli1",	"Vtn",	"Nes",	"Prom1")
f_cortex<-unique(f_cortex)
p<-DotPlot(scObj, features=f_cortex, assay="RNA", dot.min=0.1, cols=c("gray95", "purple")) + theme(axis.text.x = element_text(angle = 90, hjust=1))
p
pdf(file="split_dotplot.pdf",height=10,width=15)
print(p)
dev.off()

saveRDS(scObj, file = "data/RData/postcb_split.rds")
scObj <- readRDS("data/RData/postcb_split.rds")

DefaultAssay(scObj)<-"RNA"
Idents(scObj) <- "split_clusters"
levels(scObj)<-c("8","0","20","4","1","12","17","7","3","14","16","2","15","18","11","6","10","13","9","5","22","21","19")
our_colors<-c("#00e5d3","#15ee76","#b2eeb1", "#09cd65", "#66cdaa",  "#118941",   "#3f6cee", "#58acee","#0c9acd", "skyblue","steelblue","royalblue4","turquoise1","lightblue","#ee7942","darkred", "#FF6600","yellow2","gold","#ffbbff","#ff69b5","#ff34b3","#b3ee3a" )
p<- DimPlot(scObj, reduction = "umap", raster=FALSE,label=TRUE,cols=our_colors)+theme(aspect.ratio=1)+coord_fixed(ratio=1)+xlim(-16,15)+ylim(-16,15)
p
png(file="UMAP_allcells_label.png",width=800,height=800,res=100)
print(p)
dev.off()

r<-DimPlot(scObj, reduction = "umap", raster=FALSE,cols=our_colors)+theme(aspect.ratio=1)+xlim(-16,15)+ylim(-16,15)
r
png(file="UMAP_allcells_wolabel.png",width=800,height=800,res=100)
print(r)
dev.off()

Idents(scObj)<-"treatment"
scObj_wt <- subset(scObj, idents="WildType")
scObj_ko <- subset(scObj, idents="KnockOut")
Idents(scObj_wt)<-"split_clusters"
Idents(scObj_ko)<-"split_clusters"
levels(scObj_wt)<-c("8","0","20","4","1","12","17","7","3","14","16","2","15","18","11","6","10","13","9","5","22","21","19")
levels(scObj_ko)<-c("8","0","20","4","1","12","17","7","3","14","16","2","15","18","11","6","10","13","9","5","22","21","19")


r<-DimPlot(scObj_wt, reduction = "umap", raster=FALSE,cols=our_colors)+theme(aspect.ratio=1)+xlim(-16,15)+ylim(-16,15)
r
png(file="UMAP_wtcells_wolabel.png",width=800,height=800,res=100)
print(r)
dev.off() 

r<-DimPlot(scObj_wt, reduction = "umap",label=TRUE, raster=FALSE,cols=our_colors)+theme(aspect.ratio=1)+xlim(-16,15)+ylim(-16,15)
r
png(file="UMAP_wtcells_label.png",width=800,height=800,res=100)
print(r)
dev.off() 

r<-DimPlot(scObj_ko, reduction = "umap", raster=FALSE,cols=our_colors)+theme(aspect.ratio=1)+xlim(-16,15)+ylim(-16,15)
r
png(file="UMAP_kocells_wolabel.png",width=800,height=800,res=100)
print(r)
dev.off() 

r<-DimPlot(scObj_ko, reduction = "umap",label=TRUE, raster=FALSE,cols=our_colors)+theme(aspect.ratio=1)+xlim(-16,15)+ylim(-16,15)
r
png(file="UMAP_kocells_label.png",width=800,height=800,res=100)
print(r)
dev.off() 

#####################################################################
########DOWNSAMPLING
######################################################################

#we downsample ko to be comaprable to wt 
scObj_full <- readRDS("data/RData/postcb_split.rds")
full_meta <- scObj_full@meta.data
Idents(scObj_full) <- "treatment"
scObj_KO <- subset(scObj_full, idents="KnockOut")
KO_meta <- scObj_KO@meta.data

full_cells <- as.data.frame(cbind(rownames(scObj_KO@meta.data), as.character(scObj_KO@meta.data$split_clusters)))
colnames(full_cells) <- c("barcode", "celltype")
full_prop <- full_cells %>% dplyr::group_by(celltype) %>% tally() 
full_prop$prop <- (full_prop$n / sum(full_prop$n))*100
write.table(full_prop, file = "data/full_prop.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
# The above code determines the proper distribution of cells across clusters.
  
# props is the cluster % calculated in the full_prop.txt.  
# Note:  due to rounding errors, it was more beneficial to round these stats before sampling.
props<-c(15.5,	12.5,	3.49,	3.4,	2.1,	1.92,	0.4,	3.09,	1.78,	0.78,	1.22,	0.88,	9.9,	1.06,	0.79,	0.24,	8.52,	8.8,	5.69,	4.98,	4.48,	4.63,	3.84)
props <- props/100
types <- c("0",	"1",	"10",	"11",	"12",	"13",	"14",	"15",	"16",	"17",	"18",	"19",	"2",	"20",	"21",	"22",	"3",	"4",	"5",	"6",	"7",	"8",	"9")
# type_props contains our target proportions.  40658 is our target number of cells.
type_props <- as.data.frame(cbind(types, props))

FirstPass <- TRUE
for (t in 1:23) {
  set.seed=100
  print(t)
  ko_type <- KO_meta %>% dplyr::filter(split_clusters==type_props[t,1])
  ko_sample <- sample(rownames(ko_type), as.numeric(type_props[t,2])*40658)
  if (FirstPass) {
    ko_barcodes <- ko_sample
    FirstPass <- FALSE
  }
  else {ko_barcodes <- c(ko_barcodes, ko_sample)} 
}

# pull in all wildtype cells.
wt_cell_names <- rownames(full_meta %>% dplyr::filter(treatment=="WildType"))
sampled_barcodes <- c(ko_barcodes, wt_cell_names)
scObj <- subset(scObj_full, cells=sampled_barcodes)
saveRDS(scObj, file = "data/RData/postcb_downsampled.rds")

############################################################################
##############PLOTS
##########################################################################
scObj<-readRDS("data/RData/postcb_downsampled.rds")
DefaultAssay(scObj)<-"RNA"

metaa <- scObj@meta.data
metaa <- metaa %>%
  rename(seq_folder = orig.ident,
         nUMI = nCount_RNA,
         nGene = nFeature_RNA,
         Treatment=treatment)

p<- metaa %>%
  ggplot(aes(x=Treatment, fill=Treatment)) + 
  geom_bar() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("Number of Cells (After)") +
  theme(axis.text.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.text.x= element_text(face="bold", size=12, colour="black"))+
  theme(legend.text= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.x= element_text(face="bold", size=12, colour="black"))+
  theme(title= element_text(face="bold", size=12, colour="black"))+
  ylab("Count")+
  theme(aspect.ratio=1)
p
png(file="cells_treatment_ds.png",width=800,height=800,res=100)
print(p)
dev.off()

p<- metaa %>%
  ggplot(aes(x=id, fill=id)) + 
  geom_bar() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("Number of Cells (After)") +
  theme(axis.text.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.text.x= element_text(face="bold", size=12, colour="black"))+
  theme(legend.text= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.x= element_text(face="bold", size=12, colour="black"))+
  theme(title= element_text(face="bold", size=12, colour="black"))+
  ylab("Count")+
  xlab("")+
  theme(aspect.ratio=1)
p
png(file="cells_id_ds.png",width=800,height=800,res=100)
print(p)
dev.off()

Idents(scObj) <- "id"
p <- VlnPlot(scObj, features = "nCount_RNA", pt.size = 0)
p <- p + xlab("") + ylab("Number of UMIs") + ggtitle("Number of UMIs per sample (After)")+
  theme(axis.text.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.text.x= element_text(face="bold", size=12, colour="black"))+
  theme(legend.text= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.x= element_text(face="bold", size=12, colour="black"))+
  theme(title= element_text(face="bold", size=12, colour="black"))+
  theme(aspect.ratio=1)
p
pdf(file="umi_id_ds.pdf")
print(p)
dev.off()

p <- VlnPlot(scObj, features = "nFeature_RNA", pt.size = 0)
p <- p + xlab("") + ylab("Number of genes") + ggtitle("Number of genes per sample (After)")+
  theme(axis.text.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.text.x= element_text(face="bold", size=12, colour="black"))+
  theme(legend.text= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.x= element_text(face="bold", size=12, colour="black"))+
  theme(title= element_text(face="bold", size=12, colour="black"))+
  theme(aspect.ratio=1)
p
pdf(file="genes_id_ds.pdf")
print(p)
dev.off()

Idents(scObj) <- "split_clusters"
levels(scObj)<-c("8","0","20","4","1","12","17","7","3","14","16","2","15","18","11","6","10","13","9","5","22","21","19")
p <- VlnPlot(scObj, features = "nCount_RNA", pt.size = 0)
p <- p + xlab("") + ylab("Number of UMIs") + ggtitle("Number of UMIs per sample (After)")+
  theme(axis.text.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.text.x= element_text(face="bold", size=12, colour="black"))+
  theme(legend.text= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.x= element_text(face="bold", size=12, colour="black"))+
  theme(title= element_text(face="bold", size=12, colour="black"))+
  scale_x_discrete(labels=c("L2","L2/3","L2/3/5/6","L2/3/4","L4/5","L5(1)","L5(2)","L5(3)","L5/6(1)","L5/6(2)","L5/6(3)","L6(1)","L6(2)","L6(3)","Inh1","Inh2","Inh3","OPC","Olig","OLig Astro","Micro1","Micro2","Endo"))
p
pdf(file="umi_clusters_ds.pdf",height=5,width=10)
print(p)
dev.off()

p <- VlnPlot(scObj, features = "nFeature_RNA", pt.size = 0)
p <- p + xlab("") + ylab("Number of genes") + ggtitle("Number of genes per sample (After)")+
  theme(axis.text.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.text.x= element_text(face="bold", size=12, colour="black"))+
  theme(legend.text= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.x= element_text(face="bold", size=12, colour="black"))+
  theme(title= element_text(face="bold", size=12, colour="black"))+
  scale_x_discrete(labels=c("L2","L2/3","L2/3/5/6","L2/3/4","L4/5","L5(1)","L5(2)","L5(3)","L5/6(1)","L5/6(2)","L5/6(3)","L6(1)","L6(2)","L6(3)","Inh1","Inh2","Inh3","OPC","Olig","OLig Astro","Micro1","Micro2","Endo"))
p
pdf(file="genes_clusters_ds.pdf",height=5,width=10)
print(p)
dev.off()

left<-scObj@meta.data
left$cells<-rownames(left)
right<-read.table(file="results/metadata_nf_cellbender.txt",sep="\t",header=TRUE)
nf<-right[,c(10,11)]
metadata<-merge(left,nf,by=c("cells"))
metadata<-metadata %>% mutate(UMI_bin = cut(nCount_RNA, breaks=c(0, 999,2000,3000,4000,5000,6000,7000,8000,9000,10000)))
p<-ggplot(metadata,aes(x=factor(UMI_bin),y=nuclear_fraction,fill=factor(UMI_bin)))+
  geom_violin()+ 
  theme_classic()+
  ylab("Intronic reads ratio") + 
  xlab("UMI bins")+
  theme(axis.text.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.text.x= element_text(face="bold", size=12, colour="black"))+
  theme(legend.text= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.x= element_text(face="bold", size=12, colour="black"))+
  theme(title= element_text(face="bold", size=12, colour="black"))  +
  theme(legend.title=element_blank())+
  ylim(0,1)+
  scale_x_discrete(labels=c('1000-2000','2000-3000','3000-4000','4000-5000','5000-6000','6000-7000','7000-8000','8000-9000','9000-10000' ))+
  scale_fill_discrete(labels=c('1000-2000','2000-3000','3000-4000','4000-5000','5000-6000','6000-7000','7000-8000','8000-9000','9000-10000'))
p
pdf(file="irr_ds.pdf",height=6,width=12)
print(p)
dev.off()

level_order<-c("8","0","20","4","1","12","17","7","3","14","16","2","15","18","11","6","10","13","9","5","22","21","19")
p<-ggplot(metadata,aes(x=factor(split_clusters),y=nuclear_fraction,fill=factor(split_clusters)))+
  geom_violin()+ 
  theme_classic()+
  ylab("Intronic reads ratio") + 
  xlab("Clusters")+
  theme(axis.text.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.text.x= element_text(face="bold", size=12, colour="black",angle=45,hjust=1))+
  theme(legend.text= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.x= element_text(face="bold", size=12, colour="black"))+
  theme(title= element_text(face="bold", size=12, colour="black"))  +
  theme(legend.title=element_blank())+
  ylim(0,1)+
  scale_x_discrete(labels=c("L2","L2/3","L2/3/5/6","L2/3/4","L4/5","L5(1)","L5(2)","L5(3)","L5/6(1)","L5/6(2)","L5/6(3)","L6(1)","L6(2)","L6(3)","Inh1","Inh2","Inh3","OPC","Olig","OLig Astro","Micro1","Micro2","Endo"))
p$data$split_clusters <- factor(x=p$data$split_clusters, levels=level_order)
p
pdf(file="irr_clusters_ds.pdf",width=10,height=5)
print(p)
dev.off()

# Plot per treatment, plot KO on top of WT and vice-versa to easily see differences.
p<-DimPlot(scObj, reduction="umap", label=FALSE, raster=FALSE, group.by="treatment", order=c("KnockOut", "WildType"), cols=c("purple", "yellow"))+ ggtitle("")+theme(aspect.ratio=1)+xlim(-16,15)+ylim(-16,15)
p
png(file="umap_ko_over_wt.png",width=800,height=800,res=100)
print(p)
dev.off()

q<-DimPlot(scObj, reduction="umap", label=FALSE, raster=FALSE, group.by="treatment", order=c("WildType", "KnockOut"), cols=c("yellow", "purple"))+ ggtitle("")+theme(aspect.ratio=1)+xlim(-16,15)+ylim(-16,15)
q
png(file="umap_wt_over_ko.png",width=800,height=800,res=100)
print(q)
dev.off()

Idents(scObj) <- "split_clusters"
levels(scObj)<-c("8","0","20","4","1","12","17","7","3","14","16","2","15","18","11","6","10","13","9","5","22","21","19")
our_colors<-c("#00e5d3","#15ee76","#b2eeb1", "#09cd65", "#66cdaa",  "#118941",   "#3f6cee", "#58acee","#0c9acd", "skyblue","steelblue","royalblue4","turquoise1","lightblue","#ee7942","darkred", "#FF6600","yellow2","gold","#ffbbff","#ff69b5","#ff34b3","#b3ee3a" )

p<- DimPlot(scObj, reduction = "umap", raster=FALSE,label=TRUE,cols=our_colors)+theme(aspect.ratio=1)+coord_fixed(ratio=1)+xlim(-16,15)+ylim(-16,15)
p
png(file="UMAP_allcells_label_ds.png",width=800,height=800,res=100)
print(p)
dev.off()

r<-DimPlot(scObj, reduction = "umap", raster=FALSE,cols=our_colors)+theme(aspect.ratio=1)+xlim(-16,15)+ylim(-16,15)
r
png(file="UMAP_allcells_wolabel_ds.png",width=800,height=800,res=100)
print(r)
dev.off()

Idents(scObj)<-"treatment"
scObj_wt <- subset(scObj, idents="WildType")
scObj_ko <- subset(scObj, idents="KnockOut")
Idents(scObj_wt)<-"split_clusters"
Idents(scObj_ko)<-"split_clusters"
levels(scObj_wt)<-c("8","0","20","4","1","12","17","7","3","14","16","2","15","18","11","6","10","13","9","5","22","21","19")
levels(scObj_ko)<-c("8","0","20","4","1","12","17","7","3","14","16","2","15","18","11","6","10","13","9","5","22","21","19")

DefaultAssay(scObj_wt)<-"RNA"
p<- FeaturePlot(scObj_wt, features = c("Kdm5a"),raster=FALSE,cols=c("grey","darkred")) & WhiteBackground() & theme(aspect.ratio=1) &xlim(-16,15)&ylim(-16,15)
p<- FeaturePlot(scObj_wt, features = c("Kdm5a"),raster=FALSE,cols=c("lightgrey","#006600")) & WhiteBackground() & theme(aspect.ratio=1) &xlim(-16,15)&ylim(-16,15)
p
png(file="kdm5a_exp_wt_green.png",width=800,height=800,res=100)
print(p)
dev.off() 

r<-DimPlot(scObj_wt, reduction = "umap", raster=FALSE,cols=our_colors)+theme(aspect.ratio=1)+xlim(-16,15)+ylim(-16,15)
r
png(file="UMAP_wtcells_wolabel_ds.png",width=800,height=800,res=100)
print(r)
dev.off() 

r<-DimPlot(scObj_wt, reduction = "umap",label=TRUE, raster=FALSE,cols=our_colors)+theme(aspect.ratio=1)+xlim(-16,15)+ylim(-16,15)
r
png(file="UMAP_wtcells_label_ds.png",width=800,height=800,res=100)
print(r)
dev.off() 

r<-DimPlot(scObj_ko, reduction = "umap", raster=FALSE,cols=our_colors)+theme(aspect.ratio=1)+xlim(-16,15)+ylim(-16,15)
r
png(file="UMAP_kocells_wolabel_ds.png",width=800,height=800,res=100)
print(r)
dev.off() 

r<-DimPlot(scObj_ko, reduction = "umap",label=TRUE, raster=FALSE,cols=our_colors)+theme(aspect.ratio=1)+xlim(-16,15)+ylim(-16,15)
r
png(file="UMAP_kocells_label_ds.png",width=800,height=800,res=100)
print(r)
dev.off() 

Idents(scObj)<-"split_clusters"
levels(scObj)<-c("8","0","20","4","1","12","17","7","3","14","16","2","15","18","11","6","10","13","9","5","22","21","19")
f_cortex<-c("Slc17a7",	"Nrgn",	"Syt17",	"Cdh13",	"Adamts2",	"Stard8",	"Rrad",	"Cux2",	"Calb1",	"Rorb",	"Inhba",	"Arf5",	"Itga7",	"Whrn",	"Col27a1",	"Pld5",	"Lypd1",	"Gpr88",	"Tmem163",	"Arhgap25",	"Cpa6",	"Pcdh19",	"Satb2",	"Tcerg1",	"Foxo1",	"Deptor",	"Pde1c",	"Nnat",	"Col6a1",	"Foxp2",	"Oprk1",	"Penk",	"Col27a1",	"Col23a1",	"Car3",	"Sulf1",	"Sla",	"Mgp",	"Car12",	"Rgs12",	"Gad1",	"Gad2",	"Pvalb",	"Sst",	"Vip",	"Reln",	"Cck",	"Cnr1",	"Npy",	"Ndnf",	"Nos1",	"Hepacam",	"Pdgfra",	"Vcan",	"Cspg4",	"Sox10",	"C1ql1",	"Mbp",	"Mobp",	"Mog",	"Gja1",	"F3",	"Aqp4",	"Gfap",	"Csf1r",	"Inpp5d",	"Cx3cr1",	"Apbb1ip",	"Itgam","Ctss",	"Flt1",	"Cldn5",	"Slc2a1",	"Fli1",	"Vtn",	"Nes",	"Prom1")
f_cortex<-unique(f_cortex)
p<-DotPlot(scObj, features=f_cortex, assay="RNA", dot.min=0.1, cols=c("gray95", "purple")) + theme(axis.text.x = element_text(angle = 90, hjust=1))
p
pdf(file="dotplot_ds.pdf",height=10,width=15)
print(p)
dev.off()

scObj_meta <- scObj@meta.data
scObj_df <- scObj_meta %>%
  dplyr::group_by(split_clusters, treatment) %>% dplyr::tally()
colnames(scObj_df) <- c("Cluster", "Treatment", "Number")

level_order<-c("8","0","20","4","1","12","17","7","3","14","16","2","15","18","11","6","10","13","9","5","22","21","19")
p<-ggplot(scObj_df, aes(fill=Treatment, x=factor(Cluster, level=rev(level_order)), y=Number)) +
  geom_bar(position="fill", stat="identity") + theme_tufte() +
  coord_flip() + geom_abline(intercept=0.5, slope=0, color="black", size=.75) +
  geom_abline(intercept=0.4, slope=0, color="gray", linetype="dashed", size=.75) +
  geom_abline(intercept=0.6, slope=0, color="gray", linetype="dashed", size=.75) + 
  ggtitle("Cell Proportion") +  scale_fill_manual(values=c("#ddb2e0","#ffef8b")) +
  xlab("Cluster") +
  guides(fill=guide_legend(reverse=TRUE))
p
pdf(file=paste(region,"_chisq_barchart.pdf",sep=""))
print(p)
dev.off()

scObj_df_1 <- scObj_meta %>%
  dplyr::group_by(split_clusters, id,treatment) %>% dplyr::tally()
colnames(scObj_df_1) <- c("Cluster", "Replicate", "Treatment","Number")

p<-ggplot(scObj_df_1, aes(fill=Replicate, x=factor(Cluster, level=rev(level_order)), y=Number)) +
  geom_bar(position="fill", stat="identity") + theme_tufte() +
  coord_flip() + geom_abline(intercept=0.5, slope=0, color="black", size=.75) +
  geom_abline(intercept=0.4, slope=0, color="gray", linetype="dashed", size=.75) +
  geom_abline(intercept=0.6, slope=0, color="gray", linetype="dashed", size=.75) + 
  ggtitle("Cell Proportion") +  scale_fill_manual(values=c("lightpink1","lightpink2","lightpink3","lightpink4","khaki1","khaki2","khaki3","khaki4")) +
  xlab("Cluster") +
  guides(fill=guide_legend(reverse=TRUE))
p
pdf(file="cell_prop_replicate_barchart.pdf")
print(p)
dev.off()

write.table(scObj_df, file = "clusterbarnums_treat.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
write.table(scObj_df_1, file = "clusterbarnums_replicate.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

#for barplots that represent each replicate as a % of cells in each cluster
df<- read.table("results/downsampled/percentage_replicate_celltype.txt", header = TRUE, sep = "\t")
df$Cell.Type<-factor(df$Cell.Type,levels=c("Excitatory","Inhibitory","Glia","Endothelial"))
level_order=c("WT1","WT2","WT3","WT4","KO1","KO2","KO3","KO4")
our_colors <- c("#009E73","#D55E00","#F0E442","black")
r<-ggplot(df,aes(fill=factor(Cell.Type),x=factor(Replicate,level=level_order),y=Percentage))+
  geom_bar(stat="identity")+ggtitle("Cell Proportions") +
  xlab("Replicate") +
  guides(fill=guide_legend(title="Cell type"))+
  theme_classic()+
  scale_fill_manual(values=our_colors) +
  theme(axis.text.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.text.x= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.x= element_text(face="bold", size=12, colour="black"))+
  theme(legend.title= element_text(face="bold", size=12, colour="black"))+
  theme(legend.text= element_text(face="bold", size=12, colour="black"))+
  theme(title= element_text(face="bold", size=12, colour="black"))
r
pdf(file="percentage_replicate_celltype.pdf")
print(r)
dev.off()  

df<- read.table("results/downsampled/allcells_cluster_comp.txt", header = TRUE, sep = "\t")
level_order<-c("8","0","20","4","1","12","17","7","3","14","16","2","15","18","11","6","10","13","9","5","22","21","19")
our_colors<-c("#00e5d3","#15ee76","#b2eeb1", "#09cd65", "#66cdaa",  "#118941",   "#3f6cee", "#58acee","#0c9acd", "skyblue","steelblue","royalblue4","turquoise1","lightblue","#ee7942","darkred", "#FF6600","yellow2","gold","#ffbbff","#ff69b5","#ff34b3","#b3ee3a" )
r<-ggplot(df,aes(fill=factor(cluster),x="",y=Percentage))+
  geom_bar(stat="identity",width=1)+
  coord_polar("y",start=0)+
  ggtitle("Cell Composition") +
  guides(fill=guide_legend(title="Cluster"))+
  theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank())+
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_blank()) +
  scale_fill_manual(values=our_colors) +
  theme(legend.title= element_text(face="bold", size=12, colour="black"))+
  theme(legend.text= element_text(face="bold", size=12, colour="black"))+
  theme(title= element_text(face="bold", size=12, colour="black"))
r$data$cluster <- factor(x=r$data$cluster, levels=level_order)
r
pdf(file="allcells_cluster_comp_piechart.pdf")
print(r)
dev.off()

df<- read.table("results/downsampled/allcells_celltypegrp_comp.txt", header = TRUE, sep = "\t")
our_colors <- c("#009E73","#D55E00","#F0E442","black")
level_order=c("Excitatory","Inhibitory","Glia","Endothelial")
r<- ggplot(df,aes(fill=factor(Cell.Type),x="",y=Percentage))+
  geom_bar(stat="identity",width=1)+
  coord_polar("y",start=0)+
  ggtitle("Cell Composition") +
  guides(fill=guide_legend(title="Cell type"))+
  theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank())+
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_blank()) +
  scale_fill_manual(values=our_colors) +
  theme(legend.title= element_text(face="bold", size=12, colour="black"))+
  theme(legend.text= element_text(face="bold", size=12, colour="black"))+
  theme(title= element_text(face="bold", size=12, colour="black"))
r$data$Cell.Type <- factor(x=r$data$Cell.Type, levels=level_order)
r
pdf(file="allcells_celltypegrp_comp_piechart.pdf")
print(r)
dev.off()  


levels(scObj_wt)<-c("8","0","20","4","1","12","17","7","3","14","16","2","15","18","11","6","10","13","9","5","22","21","19")
levels(scObj_ko)<-c("8","0","20","4","1","12","17","7","3","14","16","2","15","18","11","6","10","13","9","5","22","21","19")
our_colors<-c("grey","grey","grey","#09cd65","grey", "grey", "#3f6cee","grey",   "grey", "skyblue","grey", "grey","turquoise1","grey","grey","grey","grey","grey","grey","grey", "grey","grey","grey")

r<-DimPlot(scObj_wt, reduction = "umap", raster=FALSE,cols=our_colors)+theme(aspect.ratio=1)+xlim(-16,15)+ylim(-16,15)
r
png(file="affected_clusters_wtcells_wolabel.png",width=800,height=800,res=100)
print(r)
dev.off() 

r<-DimPlot(scObj_wt, reduction = "umap",label=TRUE, raster=FALSE,cols=our_colors)+theme(aspect.ratio=1)+xlim(-16,15)+ylim(-16,15)
r
png(file="affected_clusters_wtcells_label.png",width=800,height=800,res=100)
print(r)
dev.off() 

r<-DimPlot(scObj_ko, reduction = "umap", raster=FALSE,cols=our_colors)+theme(aspect.ratio=1)+xlim(-16,15)+ylim(-16,15)
r
png(file="affected_clusters_kocells_wolabel.png",width=800,height=800,res=100)
print(r)
dev.off() 

r<-DimPlot(scObj_ko, reduction = "umap",label=TRUE, raster=FALSE,cols=our_colors)+theme(aspect.ratio=1)+xlim(-16,15)+ylim(-16,15)
r
png(file="affected_clusters_kocells_label.png",width=800,height=800,res=100)
print(r)
dev.off() 

#####################################################################################
##### FIND GENE MARKERS RNA
#####################################################################################

scObj<-readRDS("data/RData/postcb_downsampled.rds")
DefaultAssay(scObj)<-"RNA"

Idents(scObj) <- "split_clusters"   

clusters<-levels(scObj)
for (clust in clusters) {
  scObj.markers <- FindMarkers(scObj, assay = "RNA", ident.1=clust, test.use = "MAST")
  scObj.markers$cluster<-clust
  write.table(scObj.markers, file = paste(clust,"_marker_sc.txt",sep=""), quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")
}

#After manually merging each file produced above into one file (final_markers.txt), you read it back in. 
data <- read.table(file="results/markers/final_markers_sc_res.txt", header=TRUE, sep="\t")
genes <- data %>% dplyr::filter(data$avg_log2FC > 0) %>%
  dplyr::group_by(gene) %>% 
  dplyr::filter(n()==1)
write.table(genes, file="markers_singles_sc.txt", row.names=FALSE, col.names=TRUE, sep="\t",quote=FALSE)

########################################################################
###########DIFFERENTIALLY EXPRESSED GENES
#######################################################################
scObj<-readRDS("data/RData/postcb_downsampled.rds")
DefaultAssay(scObj)<-"RNA"
Idents(scObj) <- "treatment"   

##Pseudobulk using DESeq2
#for ko v/s wt in all cells
cts<- AggregateExpression(scObj,group.by = c("id"),assays="RNA",slot="counts",return.seurat=FALSE)
cts<-cts$RNA
colData<-data.frame(samples=colnames(cts))
rownames(colData)<-NULL
colData<-colData %>%
  mutate(condition=ifelse(grepl("KO",samples),"KnockOut","WildType")) %>%
  column_to_rownames(var="samples")
dds<- DESeqDataSetFromMatrix(countData = cts,colData = colData,design= ~ condition)
dds<-DESeq(dds)
res<-results(dds,contrast=c("condition","KnockOut","WildType"))
write.table(res, file = "pseudobulk_allcells_kovswt.txt", quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")

result <- read.csv (file="results/DEGs/pseudobulk_allcells_kovswt.txt",sep="\t")
##Make a volcano plot
result<-as.data.frame(res)
result$Differential_Expression<-"No"
result$Differential_Expression[result$log2FoldChange>0.3 & result$padj<0.05 & result$pvalue<0.05]<-"Upregulated"
result$Differential_Expression[result$log2FoldChange< -0.3 & result$padj<0.05 & result$pvalue<0.05]<-"Downregulated"
up<-result[result$Differential_Expression=="Upregulated",]
up<-up[order(up$log2FoldChange,decreasing = TRUE),]
down<-result[result$Differential_Expression=="Downregulated",]
down<-down[order(down$log2FoldChange),]
top_genes<-append((rownames(up[1:10,])),(rownames(down[1:10,])))
top_genes <- top_genes[!(endsWith(top_genes,"Rik"))]
result$gene<-rownames(result)
result$delabel<-NA
result$delabel[result$gene %in% top_genes]<-result$gene[result$gene %in% top_genes]
mycolors <- c("cornflowerblue", "brown1", "black")
names(mycolors) <- c("Downregulated", "Upregulated", "No")

p <- ggplot(data=result, aes(x=log2FoldChange, y=-log10(padj),col=Differential_Expression,label=delabel)) + 
  geom_point() +
  geom_text_repel()+
  scale_colour_manual(values = mycolors,breaks = c("Upregulated","Downregulated"))+
  geom_vline(xintercept=c(-0.3, 0.3), col="darkgreen",linetype=2) +
  geom_hline(yintercept=-log10(0.05), col="darkgreen",linetype=2) +
  xlab("Log2 fold change")+
  ylab("-log10(adjusted p-value)")+
  ggtitle("Differentially expressed genes KO v/s WT") +
  theme_minimal()+
  theme(axis.text.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.text.x= element_text(face="bold", size=12, colour="black"))+
  theme(legend.text= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.x= element_text(face="bold", size=12, colour="black"))+
  theme(title= element_text(face="bold", size=12, colour="black"))+
  theme(legend.title=element_blank())+
  theme(aspect.ratio=1)
p
png(file="pseudobulk_volcano_allcells.png",width=800,height=800,res=100)
print(p)
dev.off()

##for ko v/s wt in each cluster
cts<- AggregateExpression(scObj,group.by = c("split_clusters","id"),assays="RNA",slot="counts",return.seurat=FALSE)
cts<-cts$RNA
cts.t<-t(cts)
cts.t<-as.data.frame(cts.t)
splitRows<-gsub('_.*','',rownames(cts.t))
cts.split<-split.data.frame(cts.t,f=factor(splitRows))
cts.split.modified<-lapply(cts.split,function(x){
  rownames(x)<-gsub('.*_(.*)','\\1',rownames(x))
  t(x)
})
#run for each cluster
cts<-cts.split.modified$'22'
colData<-data.frame(samples=colnames(cts))
rownames(colData)<-NULL
colData<-colData %>%
  mutate(condition=ifelse(grepl("KO",samples),"KnockOut","WildType")) %>%
  column_to_rownames(var="samples")
dds<- DESeqDataSetFromMatrix(countData = cts,colData = colData,design= ~ condition)
dds<-DESeq(dds)
res<-results(dds,contrast=c("condition","KnockOut","WildType"))
write.table(res, file = "22_pseudobulk_kovswt.txt", quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")
summary(res)

###########################################################################
#################ATAC FILTERING
###########################################################################
scObj <- readRDS("data/RData/project-raw-intersection.rds")
DefaultAssay(scObj) <- "ATAC"

# Add number of genes per UMI for each cell
scObj$ATAC_log10GenesPerUMI <- log10(scObj$nFeature_ATAC) / log10(scObj$nCount_ATAC)

#Run additional QC on ATACSeq:
scObj <- NucleosomeSignal(scObj)
scObj <- TSSEnrichment(scObj)

metadata<-scObj@meta.data
write.table(metadata,file="results/metadata_atac.txt",sep="\t",row.names=TRUE,col.names=TRUE,quote=FALSE)

# See thresholds and make plots
p <- ggplot(metadata, aes(x=orig.ident, y=percent.mt,fill=orig.ident)) + 
  geom_violin() + 
  scale_fill_manual(values=c("steelblue"))+ 
  theme_classic()+ 
  xlab("All cells") +
  ylab("Mitochondrial percentage")+
  ggtitle("Filtering : Mitochondrial percentage")+
  theme(legend.position="none")+ geom_hline(yintercept = 2,linetype=2,color="darkred") + geom_hline(yintercept = 0,linetype=2,color="darkred")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.x= element_text(face="bold", size=12, colour="black"))+
  theme(title= element_text(face="bold", size=12, colour="black"))+
  theme(aspect.ratio=1)
p
png(file="atacfiltering_mitochondrial_percentage.png",width=800,height=800,res=100)
print(p)
dev.off()

p <- ggplot(metadata, aes(x=orig.ident, y=nCount_ATAC,fill=orig.ident)) + 
  geom_violin() + 
  scale_fill_manual(values=c("steelblue"))+ 
  theme_classic()+ 
  xlab("All cells") +
  ylab("Number of UMI")+
  ggtitle("Filtering : UMI threshold")+
  theme(legend.position="none")+ geom_hline(yintercept = 1000,linetype=2,color="darkred") + geom_hline(yintercept = 10000,linetype=2,color="darkred")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.x= element_text(face="bold", size=12, colour="black"))+
  theme(title= element_text(face="bold", size=12, colour="black"))+
  theme(aspect.ratio=1)
p
png(file="atacfiltering_umi_threshold.png",width=800,height=800,res=100)
print(p)
dev.off()

p <- ggplot(metadata, aes(x=orig.ident, y=nFeature_ATAC,fill=orig.ident)) + 
  geom_violin() + 
  scale_fill_manual(values=c("steelblue"))+ 
  theme_classic()+ 
  xlab("All cells") +
  ylab("Number of Genes")+
  ggtitle("Filtering : Gene threshold")+
  theme(legend.position="none")+ geom_hline(yintercept = 300,linetype=2,color="darkred") + geom_hline(yintercept = 10000,linetype=2,color="darkred")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.x= element_text(face="bold", size=12, colour="black"))+
  theme(title= element_text(face="bold", size=12, colour="black"))+
  theme(aspect.ratio=1)
p
png(file="atacfiltering_gene_threshold.png",width=800,height=800,res=100)
print(p)
dev.off()

# Nucleosomal data:
p <- ggplot(metadata, aes(x=orig.ident, y=nucleosome_signal,fill=orig.ident)) + 
  geom_violin() + 
  scale_fill_manual(values=c("steelblue"))+ 
  theme_classic()+ 
  xlab("All cells") +
  ylab("Nucleosome Signal")+
  ggtitle("Filtering : Nucleosome SIgnal")+
  theme(legend.position="none")+  geom_hline(yintercept = 0.8,linetype=2,color="darkred")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.x= element_text(face="bold", size=12, colour="black"))+
  theme(title= element_text(face="bold", size=12, colour="black"))+
  theme(aspect.ratio=1)
p
png(file="atacfiltering_nucleosome_signal.png",width=800,height=800,res=100)
print(p)
dev.off()

# TSS Enrichment values:
p <- ggplot(metadata, aes(x=orig.ident, y=TSS.enrichment,fill=orig.ident)) + 
  geom_violin() + 
  scale_fill_manual(values=c("steelblue"))+ 
  theme_classic()+ 
  xlab("All cells") +
  ylab("TSS enrichment")+
  ggtitle("Filtering : TSS enrichment")+
  theme(legend.position="none")+  geom_hline(yintercept = 2,linetype=2,color="darkred")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.x= element_text(face="bold", size=12, colour="black"))+
  theme(title= element_text(face="bold", size=12, colour="black"))+
  theme(aspect.ratio=1)
p
png(file="atacfiltering_TSS_enrichment.png",width=800,height=800,res=100)
print(p)
dev.off()

# Use above data to determin final quality filtering and subset below:
scObj <- subset(scObj, subset = nFeature_ATAC >= 1000 & nFeature_ATAC < 10000 & nCount_ATAC >= 300  & nCount_ATAC < 10000 & percent.mt <= 2 & TSS.enrichment > 2 & nucleosome_signal < 0.8)
metadata <- scObj@meta.data
KO <- metadata %>% dplyr::filter(treatment=="KnockOut") %>% tally()
WT <- metadata %>% dplyr::filter(treatment=="WildType") %>% tally()
KO
WT
stats <- metadata %>% group_by(id) %>% tally()

saveRDS(scObj, file = "data/RData/postcb_atac_filtered.rds")
scObj <- readRDS("data/RData/postcb_atac_filtered.rds")

DefaultAssay(scObj)<-"ATAC"

Idents(scObj) <- "id"
p <- VlnPlot(scObj, features = "nCount_ATAC", pt.size = 0)
p <- p + xlab("") + ylab("Number of UMIs") + ggtitle("ATAC- Number of UMIs per sample (After)")+
  theme(axis.text.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.text.x= element_text(face="bold", size=12, colour="black"))+
  theme(legend.text= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.x= element_text(face="bold", size=12, colour="black"))+
  theme(title= element_text(face="bold", size=12, colour="black"))+
  theme(aspect.ratio=1)
p
pdf(file="umi_id_ATAC.pdf")
print(p)
dev.off()

p <- VlnPlot(scObj, features = "nFeature_ATAC", pt.size = 0)
p <- p + xlab("") + ylab("Number of genes") + ggtitle("ATAC- Number of genes per sample (After)")+
  theme(axis.text.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.text.x= element_text(face="bold", size=12, colour="black"))+
  theme(legend.text= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.x= element_text(face="bold", size=12, colour="black"))+
  theme(title= element_text(face="bold", size=12, colour="black"))+
  theme(aspect.ratio=1)
p
pdf(file="genes_id_ATAC.pdf")
print(p)
dev.off()

scObj<-TSSEnrichment(scObj)
p<-TSSPlot(scObj)

#####################################################################################
##### Peak-calling with Atac-Seq data
#####################################################################################
# Before running this section, open a terminal and type "module load macs/2.1.2".

scObj <- readRDS("data/RData/postcb_atac_filtered.rds")
DefaultAssay(scObj) <- "ATAC"
npeaks <- CallPeaks(scObj, assay="ATAC", macs2.path="/cm/shared/apps/macs/2.1.2/bin/macs2")

# clean the peaks - remove peaks on nonstandard chromosomes and in genomic blacklist regions.
npeaks <- keepStandardChromosomes(npeaks, pruning.mode="coarse")
npeaks <- subsetByOverlaps(x=npeaks, ranges = blacklist_mm10, invert=TRUE)

save(list = c("npeaks"), file = "data/RData/npeaks.Rdata")
load (file="data/RData/npeaks.Rdata")

# Quantify the counts in each peak.
DefaultAssay(scObj) <- "ATAC"
ncounts <- FeatureMatrix(fragments = Fragments(scObj), features=npeaks, cells=colnames(scObj))
save(list = c("ncounts"), file = "data/RData/ncounts.Rdata")
load (file="data/RData/ncounts.Rdata")

frags <- GetAssayData(scObj, assay="ATAC", slot="fragments")

annotations_v79 <- GetGRangesFromEnsDb(ensdb=EnsDb.Mmusculus.v79)
genome(annotations_v79) <- "mm10"
chrStyle <- mapSeqlevels(seqlevels(annotations_v79), "UCSC")
annotations <- renameSeqlevels(annotations_v79, chrStyle)

# Create new assays in the object. 
scObj[["npeaks"]] <- CreateChromatinAssay(
  counts=ncounts,
  annotation=annotations,
  fragments=frags
)

# perform latent semantic indexing (LSI) on narrow peaks.
DefaultAssay(scObj) <- "npeaks"
# FindTopFeatures sets VariableFeatures metadata
scObj <- FindTopFeatures(scObj, min.cutoff="q0")
scObj <- RunTFIDF(scObj)
scObj <- RunSVD(scObj)

saveRDS(scObj, file = "data/RData/postcb_atac_peaks.rds")
scObj <- readRDS("data/RData/postcb_atac_peaks.rds")

#nonlinear dimenion reduction and clustering.
DefaultAssay(scObj) <- "npeaks"
scObj <- RunUMAP(object=scObj, reduction='lsi', dims=2:30, reduction.name="nplsi")
scObj <- FindNeighbors(object=scObj, reduction='lsi', graph.name="npeaks_g", dims=2:30)
scObj <- FindClusters(object=scObj, algorithm=3, graph.name="npeaks_g", resolution=c(0.05,0.1,0.2,0.3,0.4,0.5), verbose=TRUE)

DefaultAssay(scObj) <- "npeaks"
Idents(scObj) <- "npeaks_g_res.0.4"
p<- DimPlot(scObj, reduction = "nplsi", raster=FALSE,label=TRUE,order=sort(as.integer(levels(scObj)),decreasing=TRUE))+theme(aspect.ratio=1)+coord_fixed(ratio=1)
p
png(file="UMAP_atac_0.4.png",width=800,height=800,res=100)
print(p)
dev.off()

write.table(scObj@assays$npeaks@annotation, file = "npeaks_annotations.txt", quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")

saveRDS(scObj, file = "data/RData/postcb_atac_clustered.rds")   

############################################################################
#############TRANSFER CLUSTER INFORMATION
#############################################################
scObj <- readRDS("data/RData/postcb_atac_clustered.rds")

# cells are now filtered according to ATAC data.  Go back to get the split_clusters from the RNA analysis and assign to these cells.  
scObj_full <- readRDS("data/RData/postcb_split.rds")

atac_meta <- scObj@meta.data
bar_codes <- rownames(atac_meta)

scObj_full_atac <- subset(scObj_full, cells=bar_codes)
full_meta <- scObj_full_atac@meta.data

# Because RNA was filtered independently, the total full cells (which we have cell types for) 
# may not include every ATAC cell. For any unassigned cells, use cell type="other".
atac_meta_cols <- as.data.frame(rownames(atac_meta))
colnames(atac_meta_cols) <- c("barcode")
cell_types_meta <- data.frame(full_meta$split_clusters)
full_meta_cols <- data.frame(cbind(rownames(full_meta), cell_types_meta), stringsAsFactors = TRUE)
colnames(full_meta_cols) <- c("barcode", "orig_cell_type")
missing_cells <- atac_meta_cols %>% dplyr::filter(!barcode %in% full_meta_cols$barcode) 
missing_cells$orig_cell_type <- "other"
cells_meta <- merge(full_meta_cols, atac_meta_cols)
cells_meta_all <- rbind(cells_meta, missing_cells)
atac_meta$barcode <- rownames(atac_meta)
atac_meta_merged <- merge(atac_meta, cells_meta_all)
rownames(atac_meta_merged) <- atac_meta_merged$barcode
scObj@meta.data <- atac_meta_merged

# distribution numbers of all ATAC cells:
stats <- cells_meta_all %>% group_by(orig_cell_type) %>% tally()
write.table(stats, file = "ATAC_stats_clusters.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

Idents(scObj) <- "orig_cell_type"
levels(scObj)<-c("8","0","20","4","1","12","17","7","3","14","16","2","15","18","11","6","10","13","9","5","22","21","19","other")
our_colors<-c("#00e5d3","#15ee76","#b2eeb1", "#09cd65", "#66cdaa",  "#118941",   "#3f6cee", "#58acee","#0c9acd", "skyblue","steelblue","royalblue4","turquoise1","lightblue","#ee7942","darkred", "#FF6600","yellow2","gold","#ffbbff","#ff69b5","#ff34b3","#b3ee3a","white" )
p<- DimPlot(scObj, reduction = "nplsi", raster=FALSE,label=TRUE,cols=our_colors)+theme(aspect.ratio=1)+coord_fixed(ratio=1)
p

# save ATAC with cell types:
saveRDS(scObj, file = "data/RData/postcb_atac_clustered.rds")   

############################################################################
############DIFFERENTIALLY ACCESSIBLE PEAKS
#############################################################
scObj <- readRDS("data/RData/postcb_atac_clustered.rds")

DefaultAssay(scObj) <- 'npeaks'

#comparison of knockout vs wildtype
Idents(scObj) <- "treatment"
da_peaks <- FindMarkers(scObj, ident.1="KnockOut", test.use = "LR",logfc.threshold = 0.1,min.pct=0.05)
da_peaks$cluster<-"KnockOut"
da_peaks$gene<-rownames(da_peaks)

positions<-rownames(da_peaks)
closest<- ClosestFeature(scObj,positions)
colnames(closest)[7]<-"gene"
df_merged<-merge(x=da_peaks,y=closest, by="gene")
df_merged<-distinct(df_merged)

write.table(df_merged, file = "dap_kovswt.txt", quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")

#########################################################################
##########MOTIF ANALYSIS
##########################################################################
scObj <- readRDS("data/RData/postcb_atac_clustered.rds")
DefaultAssay(scObj) <- 'npeaks'

#over represented motifs
pfm <- getMatrixSet(x = JASPAR2022,opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))

# add motif information
scObj <- AddMotifs(scObj,genome = BSgenome.Mmusculus.UCSC.mm10,pfm = pfm)
motif.matrix<-CreateMotifMatrix(features=granges(scObj),pwm=pfm,genome=BSgenome.Mmusculus.UCSC.mm10)
mtx= as.data.frame(as.matrix(motif.matrix))
write.table(mtx,file="motif_matrix.txt",quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")
#same matrix can be accessed via scObj@assays$npeaks@motifs@data

saveRDS(scObj, file = "data/RData/postcb_atac_motif.rds")      
scObj <- readRDS("data/RData/postcb_atac_motif.rds")

#comparison of knockout vs wildtype - open peaks
Idents(scObj) <- "treatment"
da_peaks <- FindMarkers(scObj,ident.1="KnockOut",ident.2="WildType",logfc.threshold = 0.1,min.pct=0.05,test.use = 'LR')

da_peaks_open <- da_peaks[da_peaks$p_val_adj < 0.05 & da_peaks$avg_log2FC > 0,]
top.da.peak <- rownames(da_peaks_open)

#Finding background peaks
open.peaks <- AccessiblePeaks(scObj, idents = c("KnockOut", "WildType"))
meta.feature <- GetAssayData(scObj, assay = "npeaks", slot = "meta.features")
peaks.matched <- MatchRegionStats(meta.feature = meta.feature[open.peaks, ],query.feature = meta.feature[top.da.peak, ],n=50000)

enriched.motifs <- FindMotifs(scObj,features = top.da.peak,background=peaks.matched)
write.table(enriched.motifs,file="enriched_motif_positive.txt",quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

#comparison of knockout vs wildtype - closed peaks
da_peaks_closed <- da_peaks[da_peaks$p_val_adj < 0.05 & da_peaks$avg_log2FC < 0,]
top.da.peak <- rownames(da_peaks_closed)
open.peaks <- AccessiblePeaks(scObj, idents = c("KnockOut", "WildType"))
meta.feature <- GetAssayData(scObj, assay = "npeaks", slot = "meta.features")
peaks.matched <- MatchRegionStats(meta.feature = meta.feature[open.peaks, ],query.feature = meta.feature[top.da.peak, ],n=50000)
enriched.motifs <- FindMotifs(scObj,features = top.da.peak,background=peaks.matched)
write.table(enriched.motifs,file="enriched_motif_negative.txt",quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

MotifPlot(scObj,motifs = head(rownames(enriched.motifs)))

#Motif activity
library(BiocParallel)
register(SerialParam())
register(MulticoreParam(8))

scObj <- RunChromVAR(scObj,genome = BSgenome.Mmusculus.UCSC.mm10)

saveRDS(scObj, file = "data/RData/postcb_atac_chromvar.rds")      
scObj <- readRDS("data/RData/postcb_atac_chromvar.rds")

DefaultAssay(scObj) <- 'chromvar'
Idents(scObj) <- "treatment"
differential.activity <- FindMarkers(scObj,ident.1 = 'KnockOut',ident.2 = 'WildType',mean.fxn = rowMeans,fc.name = "avg_diff")
differential.activity$id<-rownames(differential.activity)
motif_name<-scObj@assays$npeaks@motifs@motif.names
motif_name<-as.matrix(motif_name)
motif_name<-as.data.frame(motif_name)
motif_name$id<-rownames(motif_name)
annotated_differential_activity<-merge(x=differential.activity,y=motif_name,by="id",all.x=TRUE)
annotated_differential_activity<-as.matrix(annotated_differential_activity)
write.table(annotated_differential_activity,file="motif_differential_activity.txt",quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")

df<-differential.activity[order(differential.activity$avg_diff,decreasing = TRUE),]

p<-MotifPlot(scObj,motifs = c("MA0489.2"),assay = 'npeaks')
p
pdf(file="Jun.pdf")
print(p)
dev.off()

##########################################################################
################GENE ACTIVITY AND ANNOTATE "OTHER"
#########################################################################
scObj <- readRDS("data/RData/postcb_atac_chromvar.rds")
DefaultAssay(scObj)<-"npeaks"
# compute gene activities
gene.activities <- GeneActivity(scObj)

# add the gene activity matrix to the Seurat object as a new assay
scObj[['GA']] <- CreateAssayObject(counts = gene.activities)
scObj <- NormalizeData(scObj,assay = 'GA',normalization.method = 'LogNormalize',scale.factor = median(scObj$nCount_GA))

saveRDS(scObj, file = "data/RData/postcb_atac_geneactivity.rds")   
scObj <- readRDS("data/RData/postcb_atac_geneactivity.rds")

ref <- readRDS("data/RData/postcb_downsampled.rds")

DefaultAssay(scObj)<-"GA"
DefaultAssay(ref)<-"RNA"
ref <- FindVariableFeatures(object = ref,nfeatures = 5000)

# Identify anchors
transfer.anchors <- FindTransferAnchors(reference = ref, query = scObj,reduction = "cca")
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = ref$split_clusters,
                                     weight.reduction = scObj[["lsi"]], dims = 2:30)

scObj <- AddMetaData(scObj, metadata = celltype.predictions)
meta <- scObj@meta.data
meta$orig_cell_type <- as.character(meta$orig_cell_type)
meta$final_clusters <- ifelse(meta$orig_cell_type != "other", meta$orig_cell_type, meta$predicted.id)
scObj@meta.data <- meta

Idents(scObj) <- "predicted.id"
levels(scObj)<-c("8","0","20","4","1","12","17","7","3","14","16","2","15","18","11","6","10","13","9","5","22","21","19")
our_colors<-c("#00e5d3","#15ee76","#b2eeb1", "#09cd65", "#66cdaa",  "#118941",   "#3f6cee", "#58acee","#0c9acd", "skyblue","steelblue","royalblue4","turquoise1","lightblue","#ee7942","darkred", "#FF6600","yellow2","gold","#ffbbff","#ff69b5","#ff34b3","#b3ee3a")
p<- DimPlot(scObj, reduction = "nplsi", raster=FALSE,label=TRUE,cols=our_colors)+theme(aspect.ratio=1)+coord_fixed(ratio=1)
p
pdf(file="ATAC_UMAP_predicted.id.pdf")
print(p)
dev.off()

Idents(scObj) <- "final_clusters"
levels(scObj)<-c("8","0","20","4","1","12","17","7","3","14","16","2","15","18","11","6","10","13","9","5","22","21","19")
our_colors<-c("#00e5d3","#15ee76","#b2eeb1", "#09cd65", "#66cdaa",  "#118941",   "#3f6cee", "#58acee","#0c9acd", "skyblue","steelblue","royalblue4","turquoise1","lightblue","#ee7942","darkred", "#FF6600","yellow2","gold","#ffbbff","#ff69b5","#ff34b3","#b3ee3a")
p<-DimPlot(scObj, reduction = "nplsi", raster=FALSE,label=FALSE,cols=our_colors)+theme(aspect.ratio=1)+xlim(-12,15)+ylim(-12,15)
pdf(file="ATAC_UMAP_finalclusters.pdf")
print(p)
dev.off()

p<-DimPlot(scObj, reduction = "nplsi", raster=FALSE,label=FALSE,cols=our_colors,split.by = "treatment")+theme(aspect.ratio=1)+xlim(-12,15)+ylim(-12,15)
pdf(file="ATAC_UMAP_finalclusters_split.pdf",height=8,width=16)
print(p)
dev.off()

saveRDS(scObj, file = "data/RData/postcb_atac_geneactivity.rds")  