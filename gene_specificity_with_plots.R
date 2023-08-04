# Gene specificity

library(dplyr)
library(SingleCellExperiment)
library(Seurat)
sce <- qs::qread("dataset.qs")


###### specificity using MAST (fixed effects)

sce<-as.Seurat(sce, data=NULL)
sce<-SetIdent(sce,value='cluster_celltype')



####### pairwise specificity test
pairwise_comparisons<-list()
cell_types<-unique(sce$cluster_celltype)
cell_types<-cell_types[cell_types!='Micro']

for (i in cell_types){

pairwise_comparisons[[paste('Micro_vs_',i,sep = '')]]<-FindMarkers(sce,ident.1 = 'Micro',ident.2 = i,slot='counts',only.pos=F,logfc.threshold = 0,min.pct = 0,test.use = 'MAST',features = c('HK2','TFEC','LCP2','TSPO'))

}

saveRDS(pairwise_comparisons,'pairwise_comparisons.RDS')

comparisons<-do.call(rbind,lapply(pairwise_comparisons,function(x)data.frame(x)))
write.table(comparisons,'pairwise_comparisons.txt',sep='\t',col.names = NA)




###### test specificity of HK2, LCP2,TFEC
library(Seurat)
library(ggplot2)

sce<-subset(sce, subset=barcode %in% sample(sce$barcode,10000))


gc()
sce<-NormalizeData(sce)



VlnPlot(sce, features = c('TSPO'),pt.size = 0)+NoLegend()+xlab('Cell type')

ggsave('VlnPlot_TSPO.tiff',height=4,width=4)

markers<-FindMarkers(seu, ident.1 = 'Micro',features = c('HK2','LCP2','TFEC'),test.use = 'MAST')
write.table(markers,'comparison_Micro_vs_allOther.txt',sep='\t',col.names = NA)

markers_astro<-FindMarkers(seu, ident.1 = 'Astro',features = c('HK2','LCP2','TFEC'),test.use = 'MAST',min.pct = 0,logfc.threshold = 0)
write.table(markers_astro,'comparison_Astro_vs_allOther.txt',sep='\t',col.names = NA)





DotPlot(sce,features =c("HK2",'LCP2','TFEC' )  ,group.by = 'cluster_celltype',cols = RColorBrewer::brewer.pal(8,'RdBu')[c(8,1)], dot.min = .01,scale = T,scale.by = 'radius')+
  ylab('Gene')+coord_flip()+theme(axis.text.x = element_text(angle = 45,hjust=1))
ggsave('Dotplot_TFEC_LCP2_HK2.tiff',height = 3,width = 6)


DotPlot(sce,features =  c("TSPO","HK2",'LCP2','TFEC'),group.by = 'cluster_celltype',scale.by = 'size',cols = RColorBrewer::brewer.pal(8,'RdBu')[c(8,1)], dot.min = .01,assay = 'originalexp',scale = F,idents = c('CNTRL','AD'))+
  coord_flip()+ylab('group')+theme(axis.text.x = element_text(angle = 45,hjust=1))+ ylab('Gene')+scale_size(limits=c(1,50))
#+scale_y_discrete(labels=c('CTR','AD'))
ggsave('Dotplot_TSPO_HK2_TFEC_LCP2.tiff',height = 5,width = 5)


