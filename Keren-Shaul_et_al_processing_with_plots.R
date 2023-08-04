
library(SingleCellExperiment)
library(Seurat)

integr_object<-qs::qread(file ='sce.qs')
integr_object<-as.Seurat(integr_object,data = NULL)
integr_object<-NormalizeData(integr_object)



####### comparison between DAM cells and other microglial cells
 load("dataset.rda")

 DefaultAssay(seu.integrated)<-'RNA'
 
  DimPlot(seu.integrated,reduction = 'tsne',split.by = 'genotype')
library(ggplot2)
 
  pdf('feature_dim_plots.pdf')
  
  FeaturePlot(seu.integrated, features = c('ENSMUSG00000030789'),reduction = 'tsne')+
  labs(title='Itgax')
  FeaturePlot(seu.integrated, features = c('ENSMUSG00000002985'),reduction = 'tsne')+
    labs(title='Apoe')
  FeaturePlot(seu.integrated, features = c('ENSMUSG00000007891'),reduction = 'tsne')+
    labs(title='Ctsd') 
  FeaturePlot(seu.integrated, features = c('ENSMUSG00000015568'),reduction = 'tsne')+
    labs(title='Lpl') 
  FeaturePlot(seu.integrated, features = c('ENSMUSG00000030579'),reduction = 'tsne')+
    labs(title='Tyrobp') 
  FeaturePlot(seu.integrated, features = c('ENSMUSG00000068129'),reduction = 'tsne')+
    labs(title='Cst7')
  FeaturePlot(seu.integrated, features = c('ENSMUSG00000023992'),reduction = 'tsne')+
    labs(title='Trem2')
  FeaturePlot(seu.integrated, features = c('ENSMUSG00000079293'),reduction = 'tsne')+
    labs(title='Clec7a')
  FeaturePlot(seu.integrated, features = c('ENSMUSG00000036353'),reduction = 'tsne')+
    labs(title='P2ry12')
  FeaturePlot(seu.integrated, features = c('ENSMUSG00000024610'),reduction = 'tsne')+
    labs(title='Cd74')
  FeaturePlot(seu.integrated, features = c('ENSMUSG00000052336'),reduction = 'tsne')+
    labs(title='Cx3cr1')
  FeaturePlot(seu.integrated, features = c('ENSMUSG00000061311'),reduction = 'tsne')+
    labs(title='Rag1')
  FeaturePlot(seu.integrated, features = c('ENSMUSG00000021665'),reduction = 'tsne')+
    labs(title='Hexb')
  FeaturePlot(seu.integrated, features = c('ENSMUSG00000008845'),reduction = 'tsne')+
    labs(title='Cd163')
  FeaturePlot(seu.integrated, features = c('ENSMUSG00000004612'),reduction = 'tsne')+
    labs(title='Nkg7')
  
  DimPlot(seu.integrated,reduction = 'tsne',split.by = 'genotype',label = T,group.by = 'cluster_original')
  
dev.off()

seu.integrated<-subset(seu.integrated, subset=genotype=='5XFAD')

seu.integrated<-SetIdent(seu.integrated,value = 'cluster_original')

a<-FindMarkers(seu.integrated,ident.1 = 3, ident.2 = c(0,1,2),test.use = 'MAST',logfc.threshold = 0.01,min.pct = 0.01)
write.table(a, 'DAM_vs_otherMicro.txt',sep='\t',col.names = NA)












# plot
library(ggplot2)
#tspo, p2ry12, apoe, trem2,tyrobp, Cst7


seu.integrated$DAM<-factor(seu.integrated$DAM,levels = c('non-DAM','DAM'))



VlnPlot(object = seu.integrated, features = c("ENSMUSG00000068129"), pt.size =0.0, group.by = 'DAM')+  scale_x_discrete(labels=c('non-DAM','DAM'))+NoLegend()+
  labs(x='Microglia subset')+ggtitle(label = 'Cst7')






ggsave('Cst7_vln_DAM_vs_nonFAM.tiff',height=5,width=3)




####### boxplots

genes_of_interest<-as.data.frame(t(as.data.frame(seu.integrated@assays$RNA@data[c("ENSMUSG00000041736",'ENSMUSG00000036353','ENSMUSG00000002985',
                                                                                  'ENSMUSG00000023992','ENSMUSG00000030579',
                                                                                  'ENSMUSG00000068129'),])))
genes_of_interest<-as.data.frame(cbind(genes_of_interest,seu.integrated$cluster_original))
colnames(genes_of_interest)[7]<-c('subcluster')
cl<-as.data.frame(genes_of_interest$subcluster)
cl$subcluster<-'Other microglia'
cl$subcluster[cl$`genes_of_interest$subcluster`==3]<-'DAM'
genes_of_interest$subcluster<-cl$subcluster

data<-genes_of_interest
data$subcluster<-factor(data$subcluster,levels = c('DAM','Other microglia'))
# ggplot(data, aes(x=factor(data$group, levels = c('CTRL','LPS')), y=Tspo)) +
#    geom_boxplot(width=0.1)+geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1), dotsize = .1)+
#    theme_classic()+ 
#    labs( x='Treatment', y='Normalized expression')+
#    theme(axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15))+
#    theme(axis.title = element_text(size = 20))


library(ggpubr)
gene_ens='ENSMUSG00000041736'
gene='Tspo'


#data<-genes_of_interest[genes_of_interest$Tspo>0,]
ggboxplot(data, x = "subcluster",
          y = gene_ens,
          combine = TRUE,
          color = "subcluster", palette = "jco",
          ylab = "Normalized counts", 
          add = "jitter",    # Add jittered points
          add.params = list(size = 0.1, jitter = 0.2)  # Point size and the amount of jittering
)+theme(legend.position = 'none')+theme(plot.title = element_text(hjust = 0.5,face = 'bold'))+ggtitle(gene)+scale_x_discrete(labels=c('DAM','non-\nDAM'))
#theme(axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15))+
#theme(axis.title = element_text(size =15))

ggsave(paste(gene,'.tiff',sep = '_'),height=4,width=2)


                      


