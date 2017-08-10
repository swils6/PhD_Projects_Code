# PMD_Function
SLW  
February 17, 2017  
Using the annotation from Price et al (2013), I want to pull all CpG sites that are located in Schroder et al (2013) placental specific partially methylated domain (PMD) list.


```r
setwd('Z:/ROBLAB1 coredata-databases/1 Samantha DATA Folder/PROJECTS')

##Illumina 450K annotation from Price et al (2013)
anno<-read.table('Uber annotation.txt',header=T)

##The PMD List does not have chromosome Y, so I will remove it here to avoid having the chromosome factor having differing levels between the annotation and the PMD list
rm<-anno[c(anno$CHR=="chrY"),]##416 probes
anno<-anno[-which(rownames(anno) %in% rownames(rm)),]
anno$CHR<-droplevels(anno$CHR)
anno$MAPINFO<-as.numeric(anno$MAPINFO)

##Subsetting anno to only the columns we will need: CHR, MAPINFO, IlmnID
anno_sub<-anno[,c("IlmnID","CHR","MAPINFO")]

##Placental PMD list from Schroder et al (2013)
PMD_List<-read.table('PMD_list.txt',header=T)
```


```r
##This function will be designed to pull probes that are located in PMD based on the Schroder et al (2013) list. This should also work on the EPIC array

##Split both the list of PMDs with the start and end genomic locations, and the Illumina probe annotation by chromosome
PMD_Chr<-split(PMD_List,PMD_List$CHR)
##head(PMD_Chr)
CpG_anno<-split(anno_sub,anno_sub$CHR)
##head(CpG_anno)

##This will result in a list of chromosomes, each with a number of dataframes for each PMD and the probes in that PMD
PMD_Probes_List<-lapply(1:23,function(CHR){
 lapply (1:nrow(PMD_Chr[[CHR]]), function(PMD){
  PMD_probes<-CpG_anno[[CHR]][which(CpG_anno[[CHR]]$MAPINFO>=PMD_Chr[[CHR]]$Start[PMD] & CpG_anno[[CHR]]$MAPINFO<=PMD_Chr[[CHR]]$End[PMD]),]
  PMD_probes$PMD<-if(nrow(PMD_probes)>0){PMD_probes$PMD<-paste(CHR,PMD_Chr[[CHR]]$Start[PMD],PMD_Chr[[CHR]]$End[PMD])}else{}
  PMD_probes
  })})

##Now I will restructure the data to have data frame that is organized by chromosome with the PMD and probes present in each PMD
PMDProbe_list <- unlist(PMD_Probes_List, recursive = FALSE)
df<-do.call("rbind",PMDProbe_list )

##write.table(df,file="PMDProbes.txt")
```

