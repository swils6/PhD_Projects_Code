# RandomSampling
SLW  
March 20, 2017  

Load Packages



```r
memory.limit(10000000000000)
```

```
## [1] 1e+13
```

```r
setwd('Z:/ROBLAB1 coredata-databases/1 Samantha DATA Folder/PROJECTS/PE_IUGR_Array/Cox Validation Cohort')

##Read in phenotype data
des<-read.csv('des.matrix.csv',header=T)
for (i in 1:nrow(des)){
  des$Row[i]<-paste(substr(des[i,"Sentrix_Position"], start=1, stop=3))
}
des$Row<- as.factor(des$Row)
str(des)
```

```
## 'data.frame':	48 obs. of  14 variables:
##  $ ParticipantID   : Factor w/ 48 levels "COX_10134","COX_10170",..: 1 2 3 4 5 6 7 8 9 10 ...
##  $ GRP             : Factor w/ 4 levels "EOPE","LOPE",..: 2 1 2 1 1 2 1 1 3 1 ...
##  $ GRP.2           : Factor w/ 5 levels "CHR.AB","EOPE",..: 3 2 3 2 2 3 2 2 4 2 ...
##  $ Plate           : Factor w/ 1 level "WG0022048-MSA4": 1 1 1 1 1 1 1 1 1 1 ...
##  $ IUGR            : Factor w/ 2 levels "NO","YES": 1 2 2 2 1 1 2 1 1 1 ...
##  $ Sex             : Factor w/ 2 levels "F","M": 2 1 1 1 1 1 2 1 1 2 ...
##  $ GA              : int  37 29 35 31 29 37 30 34 30 30 ...
##  $ BW              : int  2710 860 1800 890 1080 2620 1040 1673 1380 1320 ...
##  $ BW_SD           : num  -0.85 -1.38 -1.67 -2.07 -0.62 -0.81 -1.46 -1.45 -0.13 -0.55 ...
##  $ CH              : Factor w/ 2 levels "NO","YES": 1 2 1 1 1 1 1 1 1 1 ...
##  $ Sentrix_ID      : num  1.00e+10 1.00e+10 9.98e+09 1.00e+10 9.98e+09 ...
##  $ Sentrix_Position: Factor w/ 12 levels "R01C01","R01C02",..: 3 12 9 4 5 3 11 12 10 5 ...
##  $ Ethnicity       : Factor w/ 5 levels "Asian","Black",..: 3 3 3 2 1 3 3 5 1 2 ...
##  $ Row             : Factor w/ 6 levels "R01","R02","R03",..: 2 6 5 2 3 2 6 6 5 3 ...
```

```r
des$Sentrix_ID<-as.factor(des$Sentrix_ID)
des$GA<-as.numeric(des$GA)

##As IUGR is not fully represented in all 6 row positions, I group rows into levels 1,2,and 3 on the microarray
des$row_grouped<-revalue(des$Row,c("R01"="1","R02"="1","R03"="2","R04"="2","R05"="3",
                                   "R06"="3"))
des$row_grouped<-as.numeric(des$row_grouped)
rownames(des)<-des$ParticipantID

##Loading in DNAm Data- DNAm measured on the Illumina 450K array
load('CoxCohort.fnorm_Jan2016.RData')
Data<-exprs(PROJECT.fun)
dim(Data)
```

```
## [1] 442028     48
```

```r
all(colnames(Data)==rownames(des))##FALSE
```

```
## [1] FALSE
```

```r
Data<-Data[,rownames(des)]
all(colnames(Data)==rownames(des))##TRUE
```

```
## [1] TRUE
```

Question: Of the 1703 EOPE differentially methylated sites in the discvoery cohort, 539 (31.6%) validated in the validation cohort. Is this more than we would expect by chance?

To look into this I will sample 1703 random CpG sites from the validation cohort, run the linear model, and record how many validate, 1000 times. Then I will compare this average to the actual number of hits that validated to determine if it is more than we would expect by chance.


```r
##Take a rondom 1703 probes from the Validation cohort data, run a linear model correcting for fetal sex
random<-function(A){
  B<-A[sample(1703,replace=TRUE),]
  lapply(1:nrow(B), function(CpG){
  pheno<-des
  pheno$Mval<-B[CpG,]
  model_hits<-lm(Mval~GRP+Sex, data=pheno)
  coef(summary(model_hits))[2,4]}) ##pulls row 2, column 4 corresponds to EOPEvsPreT (Row 2), p-values (Column 4)
}

##Perform 1000 iterations of random sampling and linear modelling
p.vals<-as.vector(replicate(1000,random(Data))) ##This will run the random function 1000 times

##For each column in p.vals, how many meet nominal p<0.05
Sig<-as.vector(colSums(p.vals<0.05))

##Number of hits that were actually validated in the data
PH<-as.numeric(539)

##Compare the number of "validated" hits in the random sampling the the actual number of validated hits. Are they significantly different.
#Source code- Adapted from Rachel Edgar
perm_num<-1000 ##Number of iterations
sum(Sig>PH) + 1/(perm_num + 1)##p=0.000999001
```

```
## [1] 0.000999001
```

```r
##Yes the number of validated hits we obtain in the analysis is singificantly more than we expect by chance. 
```

