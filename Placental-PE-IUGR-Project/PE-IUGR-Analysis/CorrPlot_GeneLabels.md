# CorrPlot_GeneLabels


```r
library(limma)
library(gplots)
library(ggplot2)
library(grid)
library(gridExtra)
library(lumi)

memory.limit(10000000000000)
setwd('Z:/ROBLAB1 coredata-databases/1 Samantha DATA Folder/PROJECTS/PE_IUGR_Array/Cox Validation Cohort')

##Read in phenotype data
des<-read.csv('des.matrix.csv',header=T)
for (i in 1:nrow(des)){
  des$Row[i]<-paste(substr(des[i,"Sentrix_Position"], start=1, stop=3))
}
des$Row<- as.factor(des$Row)
str(des)
des$Sentrix_ID<-as.factor(des$Sentrix_ID)
des$GA<-as.numeric(des$GA)

##As IUGR is not fully represented in all 6 row positions, I group rows into levels 1,2,and 3 on the microarray
library(plyr)
des$row_grouped<-revalue(des$Row,c("R01"="1","R02"="1","R03"="2","R04"="2","R05"="3",
                                   "R06"="3"))
des$row_grouped<-as.numeric(des$row_grouped)
rownames(des)<-des$ParticipantID

##Loading in DNAm Data- DNAm measured on the Illumina 450K array

load('CoxCohort.fnorm_Jan2016.RData')
Data<-exprs(PROJECT.fun)
dim(Data)

all(colnames(Data)==rownames(des))##FALSE
Data<-Data[,rownames(des)]
all(colnames(Data)==rownames(des))##TRUE
```


```r
##Run the linear model looking for changes in DNAm between pathological groups, correcting for fetal sex
Des= model.matrix(~0+GRP +Sex, data = des)
#head(Des)
fit1 = lmFit(Data, Des)
fit1= eBayes(fit1)

cont.matrix = makeContrasts(PreTvsEOPE=GRPPreT-GRPEOPE,PreTvsLOPE=GRPPreT-GRPLOPE, TermvsEOPE=GRPTerm-GRPEOPE, TermvsLOPE=GRPTerm-GRPLOPE,TermvsPreT=GRPTerm-GRPPreT,Sex=SexM,levels = Des)
fitCont = contrasts.fit(fit1, cont.matrix)###getting row names of contrasts don't match col names of coefficients- not sure why-look into
EbfitCont = eBayes(fitCont)

##Pulling our differentially methylated sites for pathological groups compared to controls and fetal sex (Male vs Female)
tt_sex = topTable(EbfitCont, coef = "Sex", n = Inf)
Sex_Pval<-data.frame(CpG=rownames(tt_sex), Nominal_P=tt_sex$P.Value)
Sex<-ggplot(Sex_Pval, aes(Nominal_P))+geom_histogram(fill="grey90", color="black")+theme_bw()+xlab("Nominal P Value")+ggtitle("Sex")+
  ylim(0,50000)

tt_EOPEvsTerm = topTable(EbfitCont, coef = "TermvsEOPE", n = Inf)
EOPEvsTerm_Pval<-data.frame(CpG=rownames(tt_EOPEvsTerm), Nominal_P=tt_EOPEvsTerm$P.Value)
EOPEvsTerm<-ggplot(EOPEvsTerm_Pval, aes(Nominal_P))+geom_histogram(fill="grey90", color="black")+theme_bw()+xlab("Nominal P Value")+ggtitle("EOPEvsTerm")+
  ylim(0,50000)

tt_LOPEvsTerm = topTable(EbfitCont, coef = "TermvsLOPE", n = Inf)
LOPEvsTerm_Pval<-data.frame(CpG=rownames(tt_LOPEvsTerm), Nominal_P=tt_LOPEvsTerm$P.Value)
LOPEvsTerm<-ggplot(LOPEvsTerm_Pval, aes(Nominal_P))+geom_histogram(fill="grey90", color="black")+theme_bw()+xlab("Nominal P Value")+ggtitle("LOPEvsTerm")+
  ylim(0,50000)

tt_PreTvsTerm = topTable(EbfitCont, coef = "TermvsPreT", n = Inf)
PreTvsTerm_Pval<-data.frame(CpG=rownames(tt_PreTvsTerm), Nominal_P=tt_PreTvsTerm$P.Value)
PreTvsTerm<-ggplot(PreTvsTerm_Pval, aes(Nominal_P))+geom_histogram(fill="grey90", color="black")+theme_bw()+xlab("Nominal P Value")+ggtitle("PreTvsTerm")+
  ylim(0,50000)

tt_EOPEvsPreT = topTable(EbfitCont, coef = "PreTvsEOPE", n = Inf)
EOPEvsPreT_Pval<-data.frame(CpG=rownames(tt_EOPEvsPreT), Nominal_P=tt_EOPEvsPreT$P.Value)
EOPEvsPreT<-ggplot(EOPEvsPreT_Pval, aes(Nominal_P))+geom_histogram(fill="grey90", color="black")+theme_bw()+xlab("Nominal P Value")+ggtitle("EOPEvsPreT")+
  ylim(0,50000)
##save(tt_EOPEvsPreT,file='tt_EOPEvsPreT_Cox_Aug2016.RData')

tt_LOPEvsPreT = topTable(EbfitCont, coef = "PreTvsLOPE", n = Inf)
LOPEvsPreT_Pval<-data.frame(CpG=rownames(tt_LOPEvsPreT), Nominal_P=tt_LOPEvsPreT$P.Value)
LOPEvsPreT<-ggplot(LOPEvsPreT_Pval, aes(Nominal_P))+geom_histogram(fill="grey90", color="black")+theme_bw()+xlab("Nominal P Value")+ggtitle("LOPEvsPreT")+
  ylim(0,50000)

##p-value distribution for each comparison
#grid.arrange(Sex,EOPEvsTerm,LOPEvsTerm,PreTvsTerm,EOPEvsPreT,LOPEvsPreT,nrow=4)
```


```r
cutoff.p = 0.05  # FDR threshold

##EOPE vs PreT hits
hits_EOPEvsPreT = as.data.frame(rownames(tt_EOPEvsPreT)[tt_EOPEvsPreT$adj.P.Val < cutoff.p])
rownames(hits_EOPEvsPreT)<-hits_EOPEvsPreT$`rownames(tt_EOPEvsPreT)[tt_EOPEvsPreT$adj.P.Val < cutoff.p]`
#str(hits_EOPEvsPreT)##0

##LOPEvsTerm hits
hits_LOPEvsTerm = as.data.frame(rownames(tt_LOPEvsTerm)[tt_LOPEvsTerm$adj.P.Val < cutoff.p])
rownames(hits_LOPEvsTerm)<-hits_LOPEvsTerm$`rownames(tt_LOPEvsTerm)[tt_LOPEvsTerm$adj.P.Val < cutoff.p]`
#str(hits_LOPEvsTerm)##148
```


```r
Data_beta<-m2beta(Data)
##head(Data_beta)

fit_beta = lmFit(Data_beta, Des)
Ebfit_beta = eBayes(fit_beta)

cont.matrix = makeContrasts(PreTvsEOPE=GRPPreT-GRPEOPE,PreTvsLOPE=GRPPreT-GRPLOPE, TermvsEOPE=GRPTerm-GRPEOPE, TermvsLOPE=GRPTerm-GRPLOPE,TermvsPreT=GRPTerm-GRPPreT,Sex=SexM,levels = Des)
fitCont_beta = contrasts.fit(Ebfit_beta, cont.matrix)
EbfitCont_beta = eBayes(fitCont_beta)
ttbeta_all = topTable(EbfitCont_beta, n = Inf)
##save(ttbeta_all,file="ttbeta_all_Cox_LM.RData")
```


```r
##Do the roblab EOPE hits FDR<0.05 correlate in Cox Data?
##Read in hits in the discovery cohort that meet FDR<0.05
Roblab_EOPE0.05<-read.table('Z:/ROBLAB1 coredata-databases/1 Samantha DATA Folder/PROJECTS/PE_IUGR_Array/Robinson Cohort/Roblab_EOPE_FDR0.05.txt',header=T)

EOPEvsterm_betas<-ttbeta_all[which(rownames(ttbeta_all) %in% rownames(Roblab_EOPE0.05)),]
#head(EOPEvsterm_betas)
EOPEvsterm_betas2<-as.data.frame(EOPEvsterm_betas[,c("PreTvsEOPE")])
rownames(EOPEvsterm_betas2)<-rownames(EOPEvsterm_betas)
EOPEvsterm_betas2$BVal_EOPE<-EOPEvsterm_betas2$`EOPEvsterm_betas[, c("PreTvsEOPE")]`
EOPEvsterm_betas2$`EOPEvsterm_betas[, c("PreTvsEOPE")]`<-NULL

##
Roblab_Cox_EOPEBetas<-merge(Roblab_EOPE0.05,EOPEvsterm_betas2,by='row.names')##difference in probe number due to bad quality probes removed from Cox data
rownames(Roblab_Cox_EOPEBetas)<-Roblab_Cox_EOPEBetas$Row.names

#head(Roblab_Cox_EOPEBetas)
cor.test(Roblab_Cox_EOPEBetas$PreTvsEOPE,Roblab_Cox_EOPEBetas$BVal_EOPE)##p<2.2e-16,R=0.62
```

```
## 
## 	Pearson's product-moment correlation
## 
## data:  Roblab_Cox_EOPEBetas$PreTvsEOPE and Roblab_Cox_EOPEBetas$BVal_EOPE
## t = 212.96, df = 74014, p-value < 2.2e-16
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  0.6119035 0.6208378
## sample estimates:
##       cor 
## 0.6163905
```

```r
##Annotation
library(grid)
Corr<-grobTree(textGrob("R=0.62,p<2.2e-16", x=0.01,  y=0.97, hjust=0,
  gp=gpar(col="black", fontsize=20,fontface="bold")))

##Load Annotation to get gene name
anno<-read.table('Uber annotation.txt',header=T)
##Merge with EOPE Data
Roblab_Cox_EOPEBetas_Gene<-merge(Roblab_Cox_EOPEBetas,anno,by='row.names')##74016 probes
rownames(Roblab_Cox_EOPEBetas_Gene)<-Roblab_Cox_EOPEBetas_Gene$IlmnID
##Use only the data needed (Betas,IlmnID,GeneName)
Roblab_Cox_EOPEBetas_Gene<-Roblab_Cox_EOPEBetas_Gene[,c("PreTvsEOPE","BVal_EOPE","IlmnID","Closest_TSS_gene_name")]

##Creating a label that I would like the graph to show- with both the gene name and the probe ID
Roblab_Cox_EOPEBetas_Gene$Label<-paste(Roblab_Cox_EOPEBetas_Gene$Closest_TSS_gene_name,"_",Roblab_Cox_EOPEBetas_Gene$IlmnID)

labels<-Roblab_Cox_EOPEBetas_Gene[which(Roblab_Cox_EOPEBetas_Gene$PreTvsEOPE>0.20 & Roblab_Cox_EOPEBetas_Gene$BVal_EOPE> 0.10|Roblab_Cox_EOPEBetas_Gene$PreTvsEOPE< -0.10 & Roblab_Cox_EOPEBetas_Gene$BVal_EOPE< -0.10|Roblab_Cox_EOPEBetas_Gene$IlmnID=="cg11079619"|Roblab_Cox_EOPEBetas_Gene$IlmnID=="cg20971407"|Roblab_Cox_EOPEBetas_Gene$IlmnID=="cg01924561"|Roblab_Cox_EOPEBetas_Gene$IlmnID=="cg02494582"|Roblab_Cox_EOPEBetas_Gene$IlmnID=="cg12436772"),]

#write.table(labels,file='Tophits_intprobes.txt')

##Plot
library(ggrepel)##Package to avoid overlap of labels
EOPE_RoblabvsCox<-ggplot(Roblab_Cox_EOPEBetas_Gene, aes(x = PreTvsEOPE, y = BVal_EOPE)) +
  ylim(-0.3,0.3)+
  xlim(-0.3,0.3)+
  geom_point(alpha=0.05)+
  geom_point(data=labels, colour="red")+
  geom_text_repel(aes(label=ifelse(rownames(Roblab_Cox_EOPEBetas_Gene) %in% rownames(labels),as.character(Roblab_Cox_EOPEBetas_Gene$Closest_TSS_gene_name),'')),size=1.8,colour = "red") +
 stat_smooth(method="lm",se=FALSE) +
 theme_bw() + 
   xlab("Delta Beta_EOPE vs PreT-Robinson Lab") + 
   ylab("Delta Beta_EOPE vs PreT-Cox Lab") + 
   ggtitle("Correlation between Robinson Lab betas and Cox Lab betas:EOPEvsPreT Hits")+
  theme(
    axis.text = element_text(size = 14),
      axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  plot.title=element_text(size = 16))

EOPE_RoblabvsCox
```

![](CorrPlot_GeneLabels_files/figure-html/Labelling the top hits and other important hits on graph-1.png)<!-- -->

