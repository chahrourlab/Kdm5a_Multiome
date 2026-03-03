######################################################################
############hdWGCNA
######################################################################

# Load libraries
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(corrplot)
library(igraph)
library(WGCNA)
library(hdWGCNA)
library(enrichR)
library(GeneOverlap)
library(stringr)
library(Signac)
theme_set(theme_cowplot())
set.seed(12345)
enableWGCNAThreads(nThreads = 8)

#Read input file. This is the downsampled Kdm5a cortex file without the ATAC assay
scObj<-readRDS("RData/hdwgcna_input.rds")
DefaultAssay(scObj)<-"RNA"

##################################################
#######ALL CELLS
####################################################
scObj <- SetupForWGCNA(scObj,gene_select = "fraction", fraction = 0.05, wgcna_name = "metacells_allcells" )

#construct metacells  in each group
scObj <- MetacellsByGroups(seurat_obj = scObj,group.by = c("orig.ident", "id"), reduction = 'umap', k = 50, max_shared = 20, ident.group = 'orig.ident' )

#normalize metacell expression matrix:
scObj <- NormalizeMetacells(scObj)
scObj <- SetDatExpr(scObj,group_name = "SeuratProject", group.by='orig.ident', assay = 'RNA', slot = 'data')

# Test different soft powers:
scObj <- TestSoftPowers(scObj,networkType = 'signed')

# plot the results:
plot_list <- PlotSoftPowers(scObj)
p<- wrap_plots(plot_list, ncol=2)
pdf(file="metacells_allcells_soft_powers.pdf")
print(p)
dev.off()

# construct co-expression network:
scObj <- ConstructNetwork(scObj, soft_power=14,setDatExpr=FALSE,tom_name = 'SeuratProject')

p<- PlotDendrogram(scObj, main='All cells hdWGCNA Dendrogram')

scObj <- ScaleData(scObj, features=VariableFeatures(scObj))
scObj <- ModuleEigengenes(scObj,group.by.vars="id")
# harmonized module eigengenes:
hMEs <- GetMEs(scObj)
# module eigengenes:
MEs <- GetMEs(scObj, harmonized=FALSE)
# compute eigengene-based connectivity (kME):
scObj <- ModuleConnectivity(scObj,group.by = 'orig.ident', group_name = 'SeuratProject')

# rename the modules
scObj <- ResetModuleNames(scObj,new_name = "AllCells-M")
# plot genes ranked by kME for each module
p <- PlotKMEs(scObj, ncol=3)
pdf(file="metacells_allcells_kme.pdf")
print(p)
dev.off()

# get the module assignment table:
modules <- GetModules(scObj)
head(modules)
write.table(modules,file="metacells_allcells_module_assignment.txt",sep="\t",row.names=FALSE,quote=FALSE)

# get hub genes
hub_df <- GetHubGenes(scObj, n_hubs = 10)
write.table(hub_df,file="metacells_allcells_hubgenes.txt",sep="\t",row.names=FALSE,quote=FALSE)

saveRDS(scObj, file='RData/metacells_allcells_hdWGCNA_object.rds')

# make a featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(scObj,features='hMEs', order=TRUE )

# stitch together with patchwork
p<- wrap_plots(plot_list, ncol=3)
pdf(file="metacells_allcells_featureplots.pdf")
print(p)
dev.off()

##FEATURE PLOTS
right<- scObj@misc$metacells_treatment_allcells$hMEs
right$cells<-rownames(right)
left<- scObj@meta.data
left$cells<-rownames(left)
joined<- merge(x=left,y=right,by=c("cells"))
rownames(joined)<-joined$cells
scObj@meta.data<-joined

colours<-c("turquoise","magenta","green","red","blue","brown","yellow","pink","black")
plot_list<- list()
for (i in 1:9) {
  p <- FeaturePlot(scObj, paste("AllCells-M",i,sep=""),raster=FALSE,min.cutoff=0,col=c("lightgrey",colours[i]))
  p <- p & theme(legend.position = "right")
  pdf(file=paste("metacells_allcells_FeatureM",i,".pdf",sep=""),height=8,width=8)
  print(p)
  dev.off()
  plot_list[[i]] <- p
  p <- FeaturePlot(scObj, paste("AllCells-M",i,sep=""),raster=FALSE,min.cutoff=0,split.by="treatment",col=c("lightgrey",colours[i]))
  p <- p & theme(legend.position = "right")
  pdf(file=paste("metacells_allcells_SplitFeatureM",i,".pdf",sep=""),height=8,width=16)
  print(p)
  dev.off()
}
p<- wrap_plots(plot_list, ncol=3)
pdf(file="metacells_allcells_featureplots_ag.pdf",height=15,width=15)
print(p)
dev.off()

ModuleCorrelogram(scObj)

ModuleNetworkPlot(scObj)

scObj <- RunModuleUMAP(scObj,n_hubs = 10, n_neighbors=15, min_dist=0.1 )
# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(scObj)

# plot with ggplot
p<- ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) + geom_point(color=umap_df$color,size=umap_df$kME*2 ) + umap_theme()
pdf(file="metacells_allcells_moduleUMAP.pdf")
print(p)
dev.off()

sfari_df <- read.csv('SFARI-Gene_genes_01-16-2024release_02-23-2024export.csv')
sfari_df$gene.symbol<- str_to_title(sfari_df$gene.symbol)
sfari_df <- subset(sfari_df, gene.symbol %in% rownames(scObj))
hub_genes <- GetHubGenes(scObj, 10)
label_genes <- intersect(hub_genes$gene_name, sfari_df$gene.symbol)
set.seed(12345)
pdf(file="metacells_allcells_moduleUMAP_genes.pdf", width=10, height=10)
ModuleUMAPPlot(scObj,edge.alpha=1,sample_edges=TRUE,edge_prop=0.1, label_hubs=0 ,keep_grey_edges=FALSE,label_genes=label_genes)
dev.off()
set.seed(12345)
pdf(file="metacells_allcells_moduleUMAP_nolabel.pdf", width=10, height=10)
ModuleUMAPPlot(scObj,edge.alpha=1,sample_edges=TRUE,edge_prop=0.1, label_hubs=0 ,keep_grey_edges=FALSE)
dev.off()

##DME
group1 <- scObj@meta.data %>% subset(treatment == "KnockOut") %>% rownames
group2 <- scObj@meta.data %>% subset(treatment == "WildType") %>% rownames
DMEs <- FindDMEs(scObj,barcodes1 = group1,barcodes2 = group2,test.use='wilcox')
write.table(DMEs,file="metacells_allcells_DME_KOvsWT.txt",sep="\t",row.names=FALSE,quote=FALSE)

p<- PlotDMEsLollipop(scObj, DMEs, pvalue = "p_val_adj")
p <- ggplot(df, aes(x = avg_log2FC, y = module)) +
  geom_segment(aes(xend = 0, yend = module),col=rev(c("blue","red","yellow","green","cyan","pink",'magenta',"brown","black"))) +
  geom_vline(xintercept = 0)+
  geom_point(aes(size = n_genes),col=rev(c("blue","red","yellow","green","cyan","pink",'magenta',"brown","black"))) +
  scale_size(name = "Number of Genes", range = c(1, 5),limits=c(50,250)) +
  theme_classic()+ ylab("") + xlab("Average log2 Fold Change")
p <- p + theme(axis.text.x= element_text(face="bold", size=12, colour="black"))+ 
  theme(axis.text.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.x= element_text(face="bold", size=12, colour="black"))+
  theme(title= element_text(face="bold", size=12, colour="black"))
p
pdf(file="metacells_allcells_lollipop_AG.pdf",height=5,width=7.5)
print(p)
dev.off()

p<- PlotDMEsVolcano(scObj,DMEs,wgcna_name = 'metacells_allcells')
pdf(file="metacells_allcells_volcano.pdf")
print(p)
dev.off()

# list of clusters to loop through
clusters <- c("8","0","20","4","1","12","17","7","3","14","16","2","15","18","11","6","10","13","9","5","22","21","19")

# set up an empty dataframe for the DMEs
DMEs <- data.frame()

# loop through the clusters
for(cur_cluster in clusters){
  
  # identify barcodes for group1 and group2 in eadh cluster
  group1 <- scObj@meta.data %>% subset(split_clusters == cur_cluster & treatment == "KnockOut") %>% rownames
  group2 <- scObj@meta.data %>% subset(split_clusters == cur_cluster & treatment == "WildType") %>% rownames
  
  # run the DME test
  cur_DMEs <- FindDMEs(scObj,barcodes1 = group1,barcodes2 = group2,test.use='wilcox',pseudocount.use=0.01, wgcna_name = 'metacells_treatment_allcells')
  
  # add the cluster info to the table
  cur_DMEs$cluster <- cur_cluster
  
  # append the table
  DMEs <- rbind(DMEs, cur_DMEs)
}
head(DMEs)
write.table(DMEs,file="metacells_allcells_DME_KOvsWT_percluster.txt",sep="\t",row.names=FALSE,quote=FALSE)

# get the modules table:
modules <- GetModules(scObj)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

# make a copy of the DME table for plotting
plot_df <- DMEs
plot_df <- read.csv(file="metacells_allcells_DME_KOvsWT_percluster.txt",sep="\t")
plot_df$annotation[plot_df$cluster=="0"]<-"L2/3"
plot_df$annotation[plot_df$cluster=="1"]<-"L4/5"
plot_df$annotation[plot_df$cluster=="2"]<-"L6 (1)"
plot_df$annotation[plot_df$cluster=="3"]<-"L5/6 (1)"
plot_df$annotation[plot_df$cluster=="4"]<-"L2/3/4"
plot_df$annotation[plot_df$cluster=="5"]<-"Olig Astro"
plot_df$annotation[plot_df$cluster=="6"]<-"Inh2"
plot_df$annotation[plot_df$cluster=="7"]<-"L5 (3)"
plot_df$annotation[plot_df$cluster=="8"]<-"L2"
plot_df$annotation[plot_df$cluster=="9"]<-"Olig"
plot_df$annotation[plot_df$cluster=="10"]<-"Inh3"
plot_df$annotation[plot_df$cluster=="11"]<-"Inh1"
plot_df$annotation[plot_df$cluster=="12"]<-"L5 (1)"
plot_df$annotation[plot_df$cluster=="13"]<-"OPC"
plot_df$annotation[plot_df$cluster=="14"]<-"L5/6 (2)"
plot_df$annotation[plot_df$cluster=="15"]<-"L6 (2)"
plot_df$annotation[plot_df$cluster=="16"]<-"L5/6 (3)"
plot_df$annotation[plot_df$cluster=="17"]<-"L5 (2)"
plot_df$annotation[plot_df$cluster=="18"]<-"L6 (3)"
plot_df$annotation[plot_df$cluster=="22"]<-"Micro1"
plot_df$annotation[plot_df$cluster=="19"]<-"Endo"
plot_df$annotation[plot_df$cluster=="20"]<-"L2/3/5/6"
plot_df$annotation[plot_df$cluster=="21"]<-"Micro2"
# set the factor level for the modules so they plot in the right order:
plot_df$module <- factor(as.character(plot_df$module), levels=mods)
plot_df$cluster <- factor(as.character(plot_df$cluster), levels=clusters)
annotations <- c("L2","L2/3","L2/3/5/6","L2/3/4","L4/5","L5 (1)","L5 (2)","L5 (3)","L5/6 (1)","L5/6 (2)","L5/6 (3)","L6 (1)","L6 (2)","L6 (3)","Inh1","Inh2","Inh3","OPC","Olig","Olig Astro","Micro1","Micro2","Endo")
plot_df$annotation <- factor(as.character(plot_df$annotation), levels=rev(annotations))
# set a min/max threshold for plotting
maxval <- 2; minval <- -2
plot_df$avg_log2FC <- ifelse(plot_df$avg_log2FC > maxval, maxval, plot_df$avg_log2FC)
plot_df$avg_log2FC <- ifelse(plot_df$avg_log2FC < minval, minval, plot_df$avg_log2FC)

# add significance levels
plot_df$Significance <- gtools::stars.pval(plot_df$p_val_adj)
plot_df$Significance <- ifelse(plot_df$avg_log2FC > 0.3 | plot_df$avg_log2FC < -0.3 , plot_df$Significance,'')

# change the text color to make it easier to see 
plot_df$textcolor <- ifelse(plot_df$avg_log2FC > 0.3 | plot_df$avg_log2FC < -0.3 ,  'white','black')

# make the heatmap with geom_tile
p <- plot_df %>% ggplot(aes(y=cluster, x=module, fill=avg_log2FC)) +geom_tile() 
p <- plot_df %>% ggplot(aes(y=annotation, x=module, fill=avg_log2FC)) +geom_tile() 
# add the significance levels
p <- p + geom_text(label=plot_df$Significance, color=plot_df$textcolor) 

# customize the color and theme of the plot
p <- p + 
  scale_fill_gradient2(low='darkblue', mid='white', high='darkred') +
  RotatedAxis() +
  theme(
    panel.border = element_rect(fill=NA, color='black', size=1),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    plot.margin=margin(0,0,0,0)
  ) + xlab('') + ylab('') +
  coord_equal()

pdf(file="metacells_allcells_heatmap_celltypes.pdf")
print(p)
dev.off()

# enrichr databases to test
dbs <- c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023')

# perform enrichment tests
scObj <- RunEnrichr(scObj,dbs=dbs, max_genes = Inf )

# retrieve the output table
enrich_df <- GetEnrichrTable(scObj)
write.table(enrich_df,file="metacells_allcells_EnrichR.txt",sep="\t",row.names=FALSE,quote=FALSE)

# make GO term plots:
EnrichrBarPlot(scObj,outdir = "enrichr_plots", n_terms = 10, plot_size = c(5,7), logscale=TRUE)

# enrichr dotplot
p<-EnrichrDotPlot(scObj,mods = "all", database = "GO_Biological_Process_2023", n_terms=1)
pdf(file="metacells_allcells_EnrichR_1.pdf")
print(p)
dev.off()

p<-EnrichrDotPlot(scObj,mods = "all", database = "GO_Cellular_Component_2023", n_terms=1)
pdf(file="metacells_allcells_EnrichR_2.pdf")
print(p)
dev.off()

p<-EnrichrDotPlot(scObj,mods = "all", database = "GO_Molecular_Function_2023", n_terms=1)
pdf(file="metacells_allcells_EnrichR_3.pdf")
print(p)
dev.off()

saveRDS(scObj, file='RData/metacells_allcells_hdWGCNA_object2.rds')
scObj<-readRDS("RData/metacells_allcells_hdWGCNA_object2.rds")

######Module projection and preservation
seurat_ref<-readRDS("RData/metacells_allcells_hdWGCNA_object2.rds")
seurat_atac<-readRDS("RData/hdwgcna_input_geneactivity.rds")
DefaultAssay(seurat_atac) <- 'GA'
seurat_atac <- ScaleData(seurat_atac, features=VariableFeatures(seurat_atac))

# project modules
seurat_atac <- ProjectModules(seurat_obj = seurat_atac,seurat_ref = seurat_ref,group.by.vars = 'id',wgcna_name = "metacells_treatment_allcells",wgcna_name_proj="RNA_projected")
# compute module hub scores for projected modules:
seurat_atac <- ModuleExprScore(seurat_atac,n_genes = 25,method='Seurat')

right<- seurat_atac@misc$RNA_projected$module_scores
right$cells<-rownames(right)
left<- seurat_atac@meta.data
left$cells<-rownames(left)
joined<- merge(x=left,y=right,by=c("cells"))
rownames(joined)<-joined$cells
seurat_atac@meta.data<-joined
colours<-c("turquoise","magenta","green","red","blue","brown","yellow","pink","black")
plot_list<- list()
for (i in 1:9) {
  p <- FeaturePlot(seurat_atac, paste("AllCells-M",i,sep=""),reduction="nplsi",raster=FALSE,min.cutoff=0,col=c("lightgrey",colours[i]))
  p <- p & theme(legend.position = "right")
  pdf(file=paste("metacells_allcells_projected_FeatureM",i,".pdf",sep=""),height=8,width=8)
  print(p)
  dev.off()
  plot_list[[i]] <- p
  p <- FeaturePlot(seurat_atac, paste("AllCells-M",i,sep=""),reduction="nplsi",raster=FALSE,min.cutoff=0,split.by="treatment",col=c("lightgrey",colours[i]))
  p <- p & theme(legend.position = "right")
  pdf(file=paste("metacells_allcells_projected_SplitFeatureM",i,".pdf",sep=""),height=8,width=16)
  print(p)
  dev.off()
}
p<- wrap_plots(plot_list, ncol=3)
pdf(file="metacells_allcells_projected_featureplots_ag.pdf",height=15,width=15)
print(p)
dev.off()

# set dat expr for single-cell dataset:
# compute metacells
seurat_atac <- MetacellsByGroups(seurat_obj = seurat_atac,group.by = c("orig.ident", "id"), reduction = 'nplsi',k = 50,max_shared = 20, ident.group = 'orig.ident')
seurat_atac <- NormalizeMetacells(seurat_atac)

# set expression matrix for query dataset:
seurat_atac <- SetDatExpr(seurat_atac,group_name = "SeuratProject",group.by = "orig.ident")
#seurat_atac <- SetDatExpr(seurat_atac,group_name = "SeuratProject",group.by = "orig.ident",use_metacells = FALSE)

# run module preservation function
seurat_atac <- ModulePreservation(seurat_atac,seurat_ref = seurat_ref,name="allcells_preservation")
# plot the summary stats
plot_list <- PlotModulePreservation(seurat_atac,name="allcells_preservation",statistics = "summary")
p<- wrap_plots(plot_list, ncol=2)
pdf(file="metacells_allcells_preservation.pdf",height=8,width=16)
print(p)
dev.off()

plot_list$Zsummary.qual$data$module <- str_split(plot_list$Zsummary.qual$data$module, "-", simplify = TRUE)[, 2]
df_qual <- plot_list$Zsummary.qual$data
p <- ggplot(df_qual, aes(x = size, y = value, color = color, label = module)) +
  geom_point(size = 4,col=c("cyan",'magenta',"green","red","blue","brown","yellow","pink","black")) + # Adjust point size if needed
  scale_x_continuous(limits = c(50, 200)) + # Set minimum x-axis value to 50
  scale_y_continuous(limits = c(0, 60),breaks=c(0,10,20,30,40,50,60)) +
  theme_minimal() + # Adjust theme as needed
  geom_text_repel(col="black")+ theme_classic() +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 2, fill = "#999999", alpha = 0.2,linetype=0) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 10, fill = "#CCCCCC", alpha = 0.05,linetype=0)+
  xlab("Module size")+
  ylab("Z summary statistics - Quality")+ ggtitle("") +
  theme(axis.text.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.text.x= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.x= element_text(face="bold", size=12, colour="black"))+
  theme(legend.position="none")
pdf(file="zstats_qual.pdf",height=8,width=8)
print(p)
dev.off()

plot_list$Zsummary.pres$data$module <- str_split(plot_list$Zsummary.pres$data$module, "-", simplify = TRUE)[, 2]
df_pres <- plot_list$Zsummary.pres$data
p <- ggplot(df_pres, aes(x = size, y = value, color = color, label = module)) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 2, fill = "#999999", alpha = 0.2,linetype=0) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 10, fill = "#CCCCCC", alpha = 0.05,linetype=0)+
  geom_point(size = 4,col=c("cyan",'magenta',"green","red","blue","brown","yellow","pink","black")) + # Adjust point size if needed
  scale_x_continuous(limits = c(50, 200)) + # Set minimum x-axis value to 50
  scale_y_continuous(limits = c(0, 20),breaks=c(0,5,10,15,20)) +
  theme_minimal() + # Adjust theme as needed
  geom_text_repel(col="black")+ theme_classic() +
  xlab("Module size")+
  ylab("Z summary statistics - Preservation")+ ggtitle("") +
  theme(axis.text.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.text.x= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.y= element_text(face="bold", size=12, colour="black"))+
  theme(axis.title.x= element_text(face="bold", size=12, colour="black"))+
  theme(legend.position="none")
pdf(file="zstats_pres.pdf",height=8,width=8)
print(p)
dev.off()

# save the results:
saveRDS(seurat_atac, file='RData/metacells_allcells_projected.rds')
seurat_atac<-readRDS("RData/metacells_allcells_projected.rds")