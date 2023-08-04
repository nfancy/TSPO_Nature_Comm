




########## ANALYSIS with all cells (mixed model)



####### lme4
library(lme4)
library(lmerTest)
library(glmmTMB)
library(DHARMa)

results <- readxl::read_excel('Prox_Processed_Data.xlsx', 
                              sheet = "Plaque_prox_all")



results$diagnosis<-factor(results$diagnosis,levels = c('Non-plaque','Plaque'))

results$SLIDE_ID<-factor(results$SLIDE_ID)
####################################################################
results$logHLADR<-log10(results$HLADR+1)




               #
model1<-glmmTMB(formula = logHLADR~diagnosis+(1|SLIDE_ID),
                data =results,family = ziGamma(link = 'log'),ziformula =~. )
                





car::Anova(model1)
summary(model1)
plot(simulateResiduals(model1))




stats<-list()
stats[['medianCond']]<-median(results$HLADR[results$diagnosis=='Plaque'])
stats[['medianCtrl']]<-median(results$HLADR[results$diagnosis=='Non-plaque'])
stats[['IQRCond']]<-IQR(results$HLADR[results$diagnosis=='Plaque'])
stats[['IQRCtrl']]<-IQR(results$HLADR[results$diagnosis=='Non-plaque'])

stats[['FC']]<-stats$medianCond/stats$medianCtrl
stats[['log2FC']]<-log2(stats$FC)
stats[['anova']]<-as.data.frame(car::Anova(model1))

write.table(as.data.frame(stats),'stats_HLADR_plaque_non_plaque_microglia.txt',sep='\t',row.names = F)

###### boxplot
param='HLADR'
library(ggpubr)
library(tibble)
stats<-read.delim(paste("put_path_here",param,"_plaque_non_plaque_microglia.txt",sep = ''))
if (stats$anova.Pr..Chisq.>=0.05){
  stats$stars<-'NS'
}
if (stats$anova.Pr..Chisq.<.05&stats$anova.Pr..Chisq.>=0.01){
  stats$stars<-'*'
    
}
if (stats$anova.Pr..Chisq.<.01&stats$anova.Pr..Chisq.>=0.001){
  stats$stars<-'**'
  
}
if (stats$anova.Pr..Chisq.<.001&stats$anova.Pr..Chisq.>=0.0001){
  stats$stars<-'***'
  
}
if (stats$anova.Pr..Chisq.<.0001){
  stats$stars<-'****'
  
}

stats_pval<-tribble(~group1,~group2,~p,
             levels(results$diagnosis)[1],levels(results$diagnosis)[2],stats$stars)


ggboxplot(results[], x = "diagnosis",
          y = paste('log',param,sep=''),
          combine = TRUE,
          #color = "diagnosis",
          color = "diagnosis", palette = c( '#00A08A', '#FF0000'),
          ylab =substitute(expression(parameter*" density in IBA1"^"+"~"cells (log"[10]*")"),list(parameter=param)),
          xlab = '',
          fill='grey',
          add = c("jitter"),    # Add jittered points
          bxp.errorbar = T,
          bxp.errorbar.width = .2,
          add.params = list(size = 01, jitter = 0.2)  # Point size and the amount of jittering
)+theme(legend.position = 'none')+theme(plot.title = element_text(hjust = 0.5,face = 'bold'))+
  theme(axis.text.x = element_text(angle = 45,hjust=1))+
  stat_pvalue_manual(stats_pval,y.position =1.05*max(results[,paste('log',param,sep='')]),bracket.size = 1,size = 4)

ggsave(paste(param,'_density_plaque_prox.tiff',sep=''),height=5,width=2)




