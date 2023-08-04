


#######
object <- qs::qread("sce.qs")


object<-object[,colData(object)$cluster_celltype=='Micro']
gc()


#######
library(Seurat)
object <- readRDS('object.RDS')
object<-SetIdent(object)
object<-subset(object, subset=new_celltype=='Microglia')
object<-as.SingleCellExperiment(object)


source('~/Basic_R_scRNAseq/make_pseudobulk.R')

object<-make_pseudobulk(object, pseudobulk_ID = 'Case.Region', pb_columns = c('Region','Braak','Disease','Sex','p_tau','Amyloid_beta','HLA'))
gc()
qs::qsave(object,'micro_pseudobulk.qs')



gc()

########
         
           












#########

data<-object$sumDat
metadata<-object$annot_pb
rownames(metadata)<-as.character(metadata$group_sample)

data<-as.data.frame(t(data))
data<-data[,!duplicated(colnames(data))]

data<-voom(data)$E

results<-merge(data,metadata,by='row.names')
results$diagnosis<-factor(results$Disease,levels = c('Control','AD'))

results<-results[results$Region=='EC',]

# results$diagnosis<-factor(results$diagnosis,levels = c('Control','MS'))
# results$diagnosis<-factor(results$diagnosis,levels = c('CNTRL','AD'))


#########

















####### AUCell on  micro

library(AUCell)

cells_rankings<-AUCell_buildRankings(object$sumDat)

Important_GeneSets <- readRDS("~/AD/Gene_sets/Important_GeneSets.RDS")

list<-list()
list[['DAM']]<-Important_GeneSets$DAM_logFC0.25_padj0.05
list[['Homeostatic']]<-Important_GeneSets$Rangaraju_micro$Homeostatic
list[['Homeostatic']]<-Important_GeneSets$Rangaraju_micro$Homeostatic
list[['Proinflammatory']]<-Important_GeneSets$Rangaraju_micro$Pro_inflammatory
list[['Core_microglia']]<-Important_GeneSets$Butovsky_micro




cells_AUC<-AUCell_calcAUC(list,cells_rankings)

saveRDS(cells_AUC,'AUCell_gene_sets_.RDS')






##### limma comparisons
library(limma)
sce<-readRDS('sce.RDS')
sce<-subset(sce, subset=new_celltype=='Microglia')
sce<-subset(sce,subset=Gliasubcluster!='9')
gc()


cells_AUC <- readRDS("AUCell_gene_sets.RDS")




# ### quantification using limma
#
#
#
library(limma)

cluster<-'Micro1'

a<-as.data.frame(as.character(sce$Gliasubcluster))
colnames(a)<-'Gliasubcluster'

a$clusters<-NA
a$clusters[a$Gliasubcluster==cluster]<-'test'
a$clusters[a$Gliasubcluster!=cluster]<-'control'
sce$clusters<-as.character(a$clusters)


diagnosis<-sce$Disease
diagnosis<-factor(diagnosis,levels = c('Control','AD'))
cluster<-sce$Gliasubcluster
cluster<-factor(cluster)
clusters<-factor(sce$clusters, levels = c('control','test'))

region<-factor(sce$Region, levels=c('SSC','EC'))
individual<-sce$Case.Region
mm <- model.matrix(~0+cluster)


dc <- duplicateCorrelation(log2(cells_AUC@assays@data@listData$AUC[,colnames(sce)]), design=mm, block=individual)
fit <- lmFit(log2(cells_AUC@assays@data@listData$AUC[,colnames(sce)])*10, design=mm, block=individual, correlation=dc$consensus)

#fit <- lmFit(log2(cells_AUC@assays@data@listData$AUC[,colnames(sce)]), mm)

contr <- makeContrasts(contrasts = 'clusterMicro6-clusterMicro1', levels = colnames(fit$coefficients))
tmp <- contrasts.fit(fit,contrasts = contr) #coefficients = 3)

tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = 'none', n = Inf)

write.table(top.table,'Micro6_vs_Micro1_gene_sets.txt',sep='\t',col.names = NA)


## vlnplots
library(ggplot2)
sce<-readRDS('sce.RDS')
sce<-subset(sce, subset=new_celltype=='Microglia')
sce<-subset(sce,subset=Gliasubcluster!='9')
sce@meta.data<-cbind(sce@meta.data,t(getAUC(cells_AUC)))


sce<-subset(sce,subset=Gliasubcluster %in% c('Micro1','Micro6'))
gc()




VlnPlot(object = sce, features = c("TREM2"), pt.size =0.0, group.by = 'Gliasubcluster')+  scale_x_discrete(labels=c('Homeostatic','Activated'))+NoLegend()+
  labs(x='Microglia subset')#,y='Gene set enrichment')+ggtitle(label = 'Core microglial signature')






ggsave('TREM2_vln_Micro6_vs_homeost.tiff',height=5,width=3)

####### pseudobulk boxplots
#################

sce<-qs::qread('micro1_6_pseudobulk_cluster_sce.qs')

rownames(sce)<-as.character(rowData(sce)$external_gene_name)

source('~/Basic_R_scRNAseq/make_pseudobulk.R')

sce<-make_pseudobulk(sce, pseudobulk_ID = 'manifest', pb_columns = c('Region','Disease','Sex','Case.Region','Gliasubcluster'))
gc()







#### boxplots
library(ggpubr)

data<-voom(sce$sumDat)$E

metadata<-sce$annot_pb
rownames(metadata)<-as.character(metadata$group_sample)

metadata<-metadata[colnames(data),]
data<-t(data)

results<-merge(data,metadata,by='row.names')


ggboxplot(results, x = "Gliasubcluster",
          y = 'C1QB',
          combine = TRUE,
          color = "Gliasubcluster", palette = c(  '#FF0000','#00A08A'),
          ylab = 'Gene expression (normalized counts)',
          xlab='Microglia subset',
          add = c("jitter"),    # Add jittered points
          add.params = list(size = 0.1, jitter = 0.2) )+theme(legend.position = 'right')+theme(plot.title = element_text(hjust = 0.5,face = 'bold'))+
  scale_x_discrete(labels=c('Homeostatic','Activated'))+NoLegend()


ggsave('C1QB_Micro6_vs_homeost_boxplot.tiff',height=4,width=4)







