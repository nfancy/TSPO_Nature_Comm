library(SingleCellExperiment)
library(Seurat)

sce<-qs::qread('sce.qs')

sce<-as.Seurat(sce,data=NULL)

sce<-NormalizeData(sce)

table(sce$clusters)

# Visualization

DimPlot(sce, reduction = "UMAP_Liger", group.by = "clusters", label=T)
DimPlot(sce, reduction = "UMAP_Liger", label = TRUE)
DimPlot(sce, reduction = "UMAP_Liger", group.by = "manifest")

dev.off()


FeaturePlot(object = sce, features = c("TSPO",'APOE','C1QA','C1QB','P2RY12','TMEM119','HLA-DRB1','CD68'), pt.size = 0.05,  reduction='UMAP_Liger')


ggsave(filename = '.tiff',height=20,width=20)

VlnPlot(object = sce, features = c("TSPO"), pt.size = 0.05, group.by = 'clusters')
VlnPlot(object = sce, features = c('APOE'), pt.size = 0.05, group.by = 'clusters')
VlnPlot(object = sce, features = c('C1QA'), pt.size = 0.05, group.by = 'clusters')
VlnPlot(object = sce, features = c('C1QB'), pt.size = 0.05, group.by = 'clusters')
VlnPlot(object = sce, features = c('P2RY12'), pt.size = 0.05, group.by = 'clusters')
VlnPlot(object = sce, features = c('TMEM119'), pt.size = 0.05, group.by = 'clusters')
VlnPlot(object = sce, features = c('HLA-DRB1'), pt.size = 0.05, group.by = 'clusters')
VlnPlot(object = sce, features = c('CD68'), pt.size = 0.05, group.by = 'clusters')
dev.off()


sce<-SetIdent(sce,value = 'clusters')
cluster_markers<-FindAllMarkers(object = sce,only.pos = T,test.use = 'wilcox',logfc.threshold = 0.01)
saveRDS(cluster_markers,'cluster_markers_wilcoxon.RDS')


gc()


###### clusters 1, 2,5 seem to be the most activated
sce<-qs::qread('sce.qs')

a<-as.data.frame(as.character(sce$clusters))
colnames(a)<-'clusters'

a$new_clusters<-'other'


  a$new_clusters[a$clusters %in% c('5','6')]<-'Activ'
  a$new_clusters[a$clusters %in% c('3','4')]<-'Non-activ'
  
  sce$new_clusters<-as.character(a$new_clusters)
  sce$new_clusters<-factor(sce$new_clusters)

  sce<-as.Seurat(sce,data=NULL)
  sce<-NormalizeData(sce)
  qs::qsave(sce,'Seurat.qs')
  
  #####
  sce<-SetIdent(sce,value = 'new_clusters')
  
 activ_vs_non_activ<-FindMarkers(object = sce,only.pos = F,test.use = 'wilcox',slot = 'counts',logfc.threshold = 0.1,ident.1 = 'Activ',ident.2 = 'Non-activ')
  saveRDS(activ_vs_non_activ,'Activ_vs_Non_Activ.RDS')
  
  gc()
  
  
  
  ######## AUCell with known gene sets
  library(AUCell)
  library(Seurat)
  
  object<-qs::qread('Seurat.qs')
  
  
  cells_rankings<-AUCell_buildRankings(as.matrix(object@assays$originalexp@counts))
  qs::qsave(cells_rankings,'cells_rankings.qs')
  
  Important_GeneSets <- readRDS("Important_GeneSets.RDS")
  
  list<-list()
  list[['DAM']]<-Important_GeneSets$DAM_logFC0.25_padj0.05
  list[['Homeostatic']]<-Important_GeneSets$Rangaraju_micro$Homeostatic
  list[['Homeostatic']]<-Important_GeneSets$Rangaraju_micro$Homeostatic
  list[['Proinflammatory']]<-Important_GeneSets$Rangaraju_micro$Pro_inflammatory
  list[['Core_microglia']]<-Important_GeneSets$Butovsky_micro
  
  
  
  
  cells_AUC<-AUCell_calcAUC(list,cells_rankings)
  
  saveRDS(cells_AUC,'AUCell_gene_sets.RDS')
  
  object@meta.data<-cbind(object@meta.data,t(getAUC(cells_AUC)))

  VlnPlot(object = object, features = c("DAM"), pt.size = 0.05, group.by = 'clusters')
  VlnPlot(object = object, features = c("Homeostatic"), pt.size = 0.05, group.by = 'clusters')
  VlnPlot(object = object, features = c("Proinflammatory"), pt.size = 0.05, group.by = 'clusters')
  VlnPlot(object = object, features = c("Core_microglia"), pt.size = 0.05, group.by = 'clusters')
  dev.off()
  
  
  ##### limma comparisons
  library(limma)
  object<-qs::qread('Seurat.qs')



 cells_AUC <- readRDS("AUCell_gene_sets.RDS")





  # ### quantification using limma
  #
  #
  #

 limma_cluster<-'6'
 a<-as.data.frame(as.character(object$clusters))
 colnames(a)<-'clusters'
 
 a$limma_clusters<-NA
 
 
 a$limma_clusters[a$clusters == limma_cluster]<-'test'
 a$limma_clusters[a$clusters !=limma_cluster]<-'CTR'
 
 object$limma_clusters<-as.character(a$limma_clusters)

  nFeature<-scale(as.matrix(object$total_features_by_counts),scale=T,center = T)
  pc.mito<-scale(as.matrix(object$pc_mito),scale = T,center = T)
individual<-as.factor(object$manifest)
  
  cluster<-as.factor(object$limma_clusters)
  cluster<-relevel(cluster, ref = 'CTR')


  mm <- model.matrix(~0+cluster+nFeature+pc.mito+individual)

  fit <- lmFit(log2(cells_AUC@assays@data@listData$AUC[,colnames(object)]), mm)

  contr <- makeContrasts(contrasts = 'clustertest-clusterCTR', levels = colnames(fit$coefficients))
  tmp <- contrasts.fit(fit,contrasts = contr) #coefficients = 3)

  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, sort.by = 'none', n = Inf)


   limma_results_clusters<-list()
  limma_results_clusters[[cluster1_vs_others]]<-top.table1
   limma_results_clusters[['cluster1_vs_others']]<-top.table1
   limma_results_clusters[['cluster2_vs_others']]<-top.table2
   limma_results_clusters[['cluster3_vs_others']]<-top.table3
   limma_results_clusters[['cluster4_vs_others']]<-top.table4
   limma_results_clusters[['cluster5_vs_others']]<-top.table5
   limma_results_clusters[['cluster6_vs_others']]<-top.table6
  
  saveRDS(limma_results_clusters,'limma_results_clusters.RDS')
  
  
  
  #### comparison to homeostatic subclusters
  
  limma_cluster<-c('5')
  a<-as.data.frame(as.character(object$clusters))
  colnames(a)<-'clusters'
  
  a$limma_clusters<-a$clusters
  
  
  a$limma_clusters[a$clusters%in%limma_cluster]<-'Activated'
  a$limma_clusters[a$clusters %in% c('3','4')]<-'Homeostatic'
  
  object$limma_clusters<-as.character(a$limma_clusters)
  
  nFeature<-scale(as.matrix(object$total_features_by_counts),scale=T,center = T)
  pc.mito<-scale(as.matrix(object$pc_mito),scale = T,center = T)
  individual<-as.factor(object$manifest)
  
  cluster<-as.factor(object$limma_clusters)
  cluster<-relevel(cluster, ref = 'Homeostatic')
  
  
  mm <- model.matrix(~0+cluster+nFeature+pc.mito)
  
  dc <- duplicateCorrelation(log2(cells_AUC@assays@data@listData$AUC[,colnames(object)]), design=mm, block=individual)
  fit <- lmFit(log2(cells_AUC@assays@data@listData$AUC[,colnames(object)]), design=mm, block=individual, correlation=dc$consensus)
  
  contr <- makeContrasts(contrasts = 'clusterActivated-clusterHomeostatic', levels = colnames(fit$coefficients))
  tmp <- contrasts.fit(fit,contrasts = contr) #coefficients = 3)
  
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, sort.by = 'none', n = Inf)
  
  limma_clusters_vs_homeost<-list()
  limma_clusters_vs_homeost[['cluster5_vs_homeost']]<-top.table
  saveRDS(limma_clusters_vs_homeost,'limma_results_clusters_vs_Homeost.RDS')
  
  
  
  
  ######check individual genes
  
  
  
  
  
  
  #### boxplots
  library(ggpubr)
  
  results<-object@meta.data
  results$TSPO<-object@assays$originalexp@data['TSPO',]
  results$APOE<-object@assays$originalexp@data['APOE',]
  results$C1QB<-object@assays$originalexp@data['C1QB',]
  results$HLA_DRA<-object@assays$originalexp@data['HLA-DRA',]
  results$P2RY12<-object@assays$originalexp@data['P2RY12',]
  results$TMEM119<-object@assays$originalexp@data['TMEM119',]
  
  results$new_clusters<-factor(results$new_clusters,levels = c('Non-activ','Activ'))

  ggboxplot(results[results$clusters%in%c('5','3','4'),], x = "new_clusters",
            y = 'APOE',
            combine = TRUE,
            color = "new_clusters", palette = c(  '#FF0000','#00A08A'),
            ylab = 'AUCell enrichment value',
            xlab='Microglia subset',
            add = c("jitter"),    # Add jittered points
            add.params = list(size = 0.1, jitter = 0.2) )+theme(legend.position = 'right')+theme(plot.title = element_text(hjust = 0.5,face = 'bold'))+
  scale_x_discrete(labels=c('Homeostatic','Activated'))+NoLegend()

  
  ggsave('APOE_cluster5_vs_homeost.tiff',height=4,width=4)
  
  
  
  ## vlnplots
  object<-subset(object, subset=clusters %in%c('3','4','5'))
  object$new_clusters<-factor(object$new_clusters,levels = c('Non-activ','Activ'))
  object@meta.data<-cbind(object@meta.data,t(getAUC(cells_AUC)[,colnames(object)]))
  
  VlnPlot(object = object, features = c("P2RY12"), pt.size = 0.01, group.by = 'new_clusters') +  scale_x_discrete(labels=c('Homeostatic','Activated'))+NoLegend()+
labs(x='' )#,y='Gene set enrichment')
    

  
  
  
  
  ggsave('P2RY12_vln_cluster5_vs_homeost.tiff',height=5,width=3)
  
  
 