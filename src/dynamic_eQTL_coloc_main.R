library(tidyverse) ; library(glue) ; library(magrittr)
library(vroom) ; library(arrow) ; library(coloc)
library(TwoSampleMR) ; library(MRutils)
source('dynamic_eQTL_coloc_functions.R')

setwd("/scratch/richards/julian.willett/eQTLpQTLDisconnect")

#coloc work for soskic data
# soskic.major.cells.coloc.rel7 = soskicColoc(vector.specific=NULL,rel6=F,allCells=F,returnCommonQTLs=F)
soskic.major.cells.coloc.rel7 = readRDS('soskic_rel7_majorcells_allcoloc.rds')
soskic.all.cells.coloc.rel7 = soskicColoc(vector.specific=NULL,rel6=F,allCells=T,returnCommonQTLs=F)
saveRDS(soskic.all.cells.coloc.rel7,'soskic_rel7_allcells_allcoloc.rds')

#coloc work for bqc19 at same significant loci as soskic
bqc.coloc.rel7.eQTL = bqcColoc(type.of.qtl='eQTL',vector.specific=NULL,rel6=F,returnCommonQTLs=F,mr=F)
saveRDS(bqc.coloc.rel7.eQTL,'bqc_rel7_eQTL_colocall.rds')
bqc.coloc.rel7.sQTL = bqcColoc(type.of.qtl='sQTL',vector.specific=NULL,rel6=F,mr=F)
saveRDS(bqc.coloc.rel7.sQTL,'bqc_rel7_sQTL_colocall.rds')
tmp = bqcColoc('eQTL',vector.specific = c('A2',6,31182658,'noninf','ENSG00000204472'),rel6=F)
tmp = bqcColoc('sQTL',vector.specific = c('A2',11,34482745,'noninf','ENSG00000121691'),rel6=F)

#coloc work for gtex at same sig loci
gtex.coloc.eQTL = gtexColoc('eQTL',NULL,returnCommonQTLs=F,mr=F,onlySoskicLoci=T)
sens.napsa.a2 = gtexColoc('eQTL',c('A2',19,50363849,NA,'ENSG00000131400'),returnCommonQTLs=F,mr=F,onlySoskicLoci=T)
sens.napsa.b2 = gtexColoc('eQTL',c('B2',19,50362278,NA,'ENSG00000131400'),returnCommonQTLs=F,mr=F,onlySoskicLoci=T)
sens.napsa.c2 = gtexColoc('eQTL',c('C2',19,50362278,NA,'ENSG00000131400'),returnCommonQTLs=F,mr=F,onlySoskicLoci=T)

#H4.PP figures
#OUTCOME A2 SOSKIC
cd4.naive.a2 = plotDynamicPP(soskic.major.cells.coloc.rel7,'CD4_Naive','A2',c(6,11,19),c('AIF1'))
TN.a2 = plotDynamicPP(soskic.major.cells.coloc.rel7,'TN_','A2',c(6,11,19),c('AIF1'))
cd4.memory.a2 = plotDynamicPP(soskic.major.cells.coloc.rel7,'CD4_Memory','A2',c(6,11,16),c('AIF1'))
TCM.a2 = plotDynamicPP(soskic.major.cells.coloc.rel7,'TCM_','A2',c(6,11),c('AIF1'))
TEM.a2 = plotDynamicPP(soskic.major.cells.coloc.rel7,'TEM_','A2',c(6,9,19),c('HLA-A','AIF1'))

#OUTCOME B2
cd4.naive.b2 = plotDynamicPP(soskic.major.cells.coloc.rel7,'CD4_Naive','B2',c(8,11,19),c('HSD17B14'))
tn.b2 = plotDynamicPP(soskic.major.cells.coloc.rel7,'TN_','B2',c(8,11,19),c('HSD17B14'))
cd4.memory.b2 = plotDynamicPP(soskic.major.cells.coloc.rel7,'CD4_Memory','B2',c(8,11,19),c('HSD17B14'))
tcm.b2 = plotDynamicPP(soskic.major.cells.coloc.rel7,'TCM_','B2',c(8,11,19),c())
tem.b2 = plotDynamicPP(soskic.major.cells.coloc.rel7,'TEM_','B2',c(9,19),c())

#OUTCOME C2
cd4.naive.c2 = plotDynamicPP(soskic.major.cells.coloc.rel7,'CD4_Naive','C2',c(11,19),c())
TN.c2 = plotDynamicPP(soskic.major.cells.coloc.rel7,'TN_','C2',c(1,11,19),c('GBA'))
cd4.memory.c2 = plotDynamicPP(soskic.major.cells.coloc.rel7,'CD4_Memory','C2',c(11),c())
tcm.c2 = plotDynamicPP(soskic.major.cells.coloc.rel7,'TCM_','C2',c(1,11),c())
TEM.c2 = plotDynamicPP(soskic.major.cells.coloc.rel7,'TEM_','C2',c(9,19),c())

#bqc figs
bqc.a2.coloc = plotDynamicPP.BQC(bqc.coloc.rel7.eQTL,'A2',c(9,11,16,19),c('ABO'))
bqc.b2.coloc = plotDynamicPP.BQC(bqc.coloc.rel7.eQTL,'B2',c(8,9,11,19),c('ABO','NTN5'))
bqc.c2.coloc = plotDynamicPP.BQC(bqc.coloc.rel7.eQTL,'C2',c(1,9,11,19),c('NTN5','GBAP1','NAPSB','TYK2','THBS3'))

#sensitivity analyses
#OUTCOME A2 by cell
tmp = soskicColoc(vector.specific=c('A2',6,31182658,'CD4_Naive','_40h_','ENSG00000204472'),rel6=F,allCells=F,returnCommonQTLs = F,mr=F)
tmp = soskicColoc(vector.specific=c('A2',6,31182658,'CD4_Naive','_5d_','ENSG00000204472'),rel6=F,allCells=F,returnCommonQTLs = F)
tmp = soskicColoc(vector.specific=c('A2',11,34482745,'CD4_Naive','_0h_','ENSG00000121691'),rel6=F,allCells=F,returnCommonQTLs = F)
tmp = soskicColoc(vector.specific=c('A2',19,50363849,'CD4_Naive','_40h_','ENSG00000131400'),rel6=F,allCells=F,returnCommonQTLs = F)
tmp = soskicColoc(vector.specific=c('A2',6,31182658,'CD4_Memory','_40h_','ENSG00000204472'),rel6=F,allCells=F,returnCommonQTLs = F,mr=F)
tmp = soskicColoc(vector.specific=c('A2',6,31182658,'CD4_Memory','_5d_','ENSG00000204472'),rel6=F,allCells=F,returnCommonQTLs = F)
tmp = soskicColoc(vector.specific=c('A2',11,34482745,'CD4_Memory','_5d_','ENSG00000121691'),rel6=F,allCells=F,returnCommonQTLs = F)
tmp = soskicColoc(vector.specific=c('A2',16,89196249,'CD4_Memory','_5d_','ENSG00000176715'),rel6=F,allCells=F,returnCommonQTLs = F)
tmp = soskicColoc(vector.specific=c('A2',6,29946723,'TEM_40h','_40h_','ENSG00000206503'),rel6=F,allCells=T,returnCommonQTLs = F)
tmp = soskicColoc(vector.specific=c('A2',6,31182658,'TEM_40h','_40h_','ENSG00000204472'),rel6=F,allCells=T,returnCommonQTLs = F,mr=F)
tmp = soskicColoc(vector.specific=c('A2',6,31182658,'TEM_5d','_5d_','ENSG00000204472'),rel6=F,allCells=T,returnCommonQTLs = F)
tmp = soskicColoc(vector.specific=c('A2',9,133271182,'TEM_16h','_16h_','ENSG00000160271'),rel6=F,allCells=T,returnCommonQTLs = F)
tmp = soskicColoc(vector.specific=c('A2',19,50363849,'TEM_40h','_40h_','ENSG00000131400'),rel6=F,allCells=T,returnCommonQTLs = F)
tmp = soskicColoc(vector.specific=c('A2',6,31182658,'TCM_','_40h_','ENSG00000204472'),rel6=F,allCells=T,returnCommonQTLs = F,mr=F)
tmp = soskicColoc(vector.specific=c('A2',6,31182658,'TCM_','_5d_','ENSG00000204472'),rel6=F,allCells=T,returnCommonQTLs = F)
tmp = soskicColoc(vector.specific=c('A2',11,34482745,'TCM_','_5d_','ENSG00000121691'),rel6=F,allCells=T,returnCommonQTLs = F)
tmp = soskicColoc(vector.specific=c('A2',6,31182658,'TN_40h','_40h_','ENSG00000204472'),rel6=F,allCells=T,returnCommonQTLs = F,mr=F)
tmp = soskicColoc(vector.specific=c('A2',6,31182658,'TN_5d','_5d_','ENSG00000204472'),rel6=F,allCells=T,returnCommonQTLs = F)
tmp = soskicColoc(vector.specific=c('A2',11,34482745,'TN_0h','_0h_','ENSG00000121691'),rel6=F,allCells=T,returnCommonQTLs = F)
tmp = soskicColoc(vector.specific=c('A2',19,50363849,'TN_40h','_40h_','ENSG00000131400'),rel6=F,allCells=T,returnCommonQTLs = F)

#OUTCOME B2
tmp = soskicColoc(vector.specific=c('B2',8,60515641,'CD4_Naive','_16h_','ENSG00000104388'),rel6=F,allCells=F,returnCommonQTLs = F)
tmp = soskicColoc(vector.specific=c('B2',8,60515641,'CD4_Naive','_40h_','ENSG00000104388'),rel6=F,allCells=F,returnCommonQTLs = F)
tmp = soskicColoc(vector.specific=c('B2',11,34482745,'CD4_Naive','_0h_','ENSG00000121691'),rel6=F,allCells=F,returnCommonQTLs = F)
tmp = soskicColoc(vector.specific=c('B2',19,48703160,'CD4_Naive','_0h_','ENSG00000087076'),rel6=F,allCells=F,returnCommonQTLs = F,mr=F)
tmp = soskicColoc(vector.specific=c('B2',19,48703160,'CD4_Naive','_5d_','ENSG00000087076'),rel6=F,allCells=F,returnCommonQTLs = F,mr=F)
tmp = soskicColoc(vector.specific=c('B2',19,50362278,'CD4_Naive','_40h_','ENSG00000131400'),rel6=F,allCells=F,returnCommonQTLs = F)
tmp = soskicColoc(vector.specific=c('B2',8,60515641,'CD4_Memory','_16h_','ENSG00000104388'),rel6=F,allCells=F,returnCommonQTLs = F)
tmp = soskicColoc(vector.specific=c('B2',8,60515641,'CD4_Memory','_40h_','ENSG00000104388'),rel6=F,allCells=F,returnCommonQTLs = F)
tmp = soskicColoc(vector.specific=c('B2',11,34482745,'CD4_Memory','_5d_','ENSG00000121691'),rel6=F,allCells=F,returnCommonQTLs = F)
tmp = soskicColoc(vector.specific=c('B2',19,48703160,'CD4_Memory','_0h_','ENSG00000087076'),rel6=F,allCells=F,returnCommonQTLs = F,mr=F)
tmp = soskicColoc(vector.specific=c('B2',8,60515641,'TCM_','_16h_','ENSG00000104388'),rel6=F,allCells=T,returnCommonQTLs = F)
tmp = soskicColoc(vector.specific=c('B2',11,34482745,'TCM_','_5d_','ENSG00000121691'),rel6=F,allCells=T,returnCommonQTLs = F)
tmp = soskicColoc(vector.specific=c('B2',8,60515641,'TN_16h','_16h_','ENSG00000104388'),rel6=F,allCells=T,returnCommonQTLs = F)
tmp = soskicColoc(vector.specific=c('B2',8,60515641,'TN_40h','_40h_','ENSG00000104388'),rel6=F,allCells=T,returnCommonQTLs = F)
tmp = soskicColoc(vector.specific=c('B2',11,34482745,'TN_0h','_0h_','ENSG00000121691'),rel6=F,allCells=T,returnCommonQTLs = F)
tmp = soskicColoc(vector.specific=c('B2',19,48703160,'TN_0h','_0h_','ENSG00000087076'),rel6=F,allCells=T,returnCommonQTLs = F)
tmp = soskicColoc(vector.specific=c('B2',19,48703160,'TN_5d','_5d_','ENSG00000087076'),rel6=F,allCells=T,returnCommonQTLs = F)
tmp = soskicColoc(vector.specific=c('B2',19,50362278,'TN_40h','_40h_','ENSG00000131400'),rel6=F,allCells=T,returnCommonQTLs = F)
tmp = soskicColoc(vector.specific=c('B2',9,133263862,'TEM_16h','_16h_','ENSG00000160271'),rel6=F,allCells=T,returnCommonQTLs = F)
tmp = soskicColoc(vector.specific=c('B2',19,50362278,'TEM_40h','_40h_','ENSG00000131400'),rel6=F,allCells=T,returnCommonQTLs = F)

#OUTCOME C2
tmp = soskicColoc(vector.specific=c('C2',11,34432600,'CD4_Naive','_0h_','ENSG00000121691'),rel6=F,allCells=F,returnCommonQTLs = F,mr=F)
tmp = soskicColoc(vector.specific=c('C2',19,50362278,'CD4_Naive','_40h_','ENSG00000131400'),rel6=F,allCells=F,returnCommonQTLs = F)
tmp = soskicColoc(vector.specific=c('C2',11,34432600,'CD4_Memory','_5d_','ENSG00000121691'),rel6=F,allCells=F,returnCommonQTLs = F,mr=F)
tmp = soskicColoc(vector.specific=c('C2',1,155162859,'TCM_','_5d_','ENSG00000143537'),rel6=F,allCells=T,returnCommonQTLs = F)
tmp = soskicColoc(vector.specific=c('C2',11,34432600,'TCM_','_5d_','ENSG00000121691'),rel6=F,allCells=F,returnCommonQTLs = F,mr=F)
tmp = soskicColoc(vector.specific=c('C2',1,155162859,'TN_5d','_5d_','ENSG00000177628'),rel6=F,allCells=T,returnCommonQTLs = F)
tmp = soskicColoc(vector.specific=c('C2',11,34432600,'TN_0h','_0h_','ENSG00000121691'),rel6=F,allCells=F,returnCommonQTLs = F,mr=F)
tmp = soskicColoc(vector.specific=c('C2',19,50362278,'TN_40h','_40h_','ENSG00000131400'),rel6=F,allCells=T,returnCommonQTLs = F)
tmp = soskicColoc(vector.specific=c('C2',9,133271745,'TEM_16h','_16h_','ENSG00000160271'),rel6=F,allCells=T,returnCommonQTLs = F)
tmp = soskicColoc(vector.specific=c('C2',19,50362278,'TEM_40h','_40h_','ENSG00000131400'),rel6=F,allCells=T,returnCommonQTLs = F)

#BQC A2
tmp = bqcColoc('eQTL',vector.specific = c('A2',19,50363849,'noninf','ENSG00000131400'),rel6=F,returnCommonQTLs=F,mr=F)
#BQC B2
tmp = bqcColoc('eQTL',vector.specific = c('B2',8,60515641,'noninf','ENSG00000104388'),rel6=F,returnCommonQTLs=F,mr=F)
tmp = bqcColoc('eQTL',vector.specific = c('B2',8,60515641,'inf','ENSG00000104388'),rel6=F,returnCommonQTLs=F,mr=F)
tmp = bqcColoc('eQTL',vector.specific = c('B2',19,50362278,'noninf','ENSG00000131400'),rel6=F,returnCommonQTLs=F,mr=F)
#BQC C2
tmp = bqcColoc('eQTL',vector.specific = c('C2',19,50362278,'noninf','ENSG00000131400'),rel6=F,returnCommonQTLs=F,mr=F)

#conduct MR
#OUTCOME A2
mr.boskic.cat.a2.naive.0h = conductMRBoskic('A2','11:34482745','CD4_Naive','_0h_','ENSG00000121691') #ivw sig
mr.boskic.napsa.a2.naive.40h = conductMRBoskic('A2','19:50363849','CD4_Naive','_40h_','ENSG00000131400') #not ivw sig
mr.boskic.cat.a2.memory.5d = conductMRBoskic('A2','11:34482745','CD4_Memory','_5d_','ENSG00000121691') #ivw sig
mr.boskic.acsf3.a2.memory.5d = conductMRBoskic('A2','16:89196249','CD4_Memory','_5d_','ENSG00000176715') #ivw sig
mr.boskic.ralgds.a2.TEM.16h = conductMRBoskic('A2','9:133271182','TEM_16h','_16h_','ENSG00000160271')
mr.boskic.napsa.a2.TEM.40h = conductMRBoskic('A2','19:50363849','TEM_40h','_40h_','ENSG00000131400')
mr.boskic.cat.a2.TCM.5d = conductMRBoskic('A2','11:34482745','TCM_','_5d_','ENSG00000121691')
mr.boskic.cat.a2.TN.0h = conductMRBoskic('A2','11:34482745','TN_0h','_0h_','ENSG00000121691')
mr.boskic.napsa.a2.TN.40h = conductMRBoskic('A2','19:50363849','TN_40h','_40h_','ENSG00000131400')

#OUTCOME B2
mr.boskic.rab2a.b2.naive.16h = conductMRBoskic('B2','8:60515641','CD4_Naive','_16h_','ENSG00000104388') #ivw sig
mr.boskic.rab2a.b2.naive.40h = conductMRBoskic('B2','8:60515641','CD4_Naive','_40h_','ENSG00000104388') #ivw sig
mr.boskic.cat.b2.naive.0h = conductMRBoskic('B2','11:34482745','CD4_Naive','_0h_','ENSG00000121691') #ivw sig
mr.boskic.napsa.b2.naive.40h = conductMRBoskic('B2','19:50362278','CD4_Naive','_40h_','ENSG00000131400') #ivw sig
mr.boskic.rab2a.b2.memory.16h = conductMRBoskic('B2','8:60515641','CD4_Memory','_16h_','ENSG00000104388') #ivw sig
mr.boskic.rab2a.b2.memory.40h = conductMRBoskic('B2','8:60515641','CD4_Memory','_40h_','ENSG00000104388') #ivw sig
mr.boskic.cat.b2.memory.5d = conductMRBoskic('B2','11:34482745','CD4_Memory','_5d_','ENSG00000121691') #ivw sig
mr.boskic.rab2a.b2.tcm.16h = conductMRBoskic('B2','8:60515641','TCM_','_16h_','ENSG00000104388') 
mr.boskic.cat.b2.tcm.5d = conductMRBoskic('B2','11:34482745','TCM_','_5d_','ENSG00000121691') 
mr.boskic.rab2a.b2.tn.16h = conductMRBoskic('B2','8:60515641','TN_16h','_16h_','ENSG00000104388') 
mr.boskic.rab2a.b2.tn.40h = conductMRBoskic('B2','8:60515641','TN_40h','_40h_','ENSG00000104388') 
mr.boskic.cat.b2.tn.0h = conductMRBoskic('B2','11:34482745','TN_0h','_0h_','ENSG00000121691') 
mr.boskic.napsa.b2.tn.40h = conductMRBoskic('B2','19:50362278','TN_40h','_40h_','ENSG00000131400') 
mr.boskic.ralgds.b2.tem.16h = conductMRBoskic('B2','9:133263862','TEM_16h','_16h_','ENSG00000160271') 
mr.boskic.napsa.b2.tem.40h = conductMRBoskic('B2','19:50362278','TEM_40h','_40h_','ENSG00000131400') 

#OUTCOME C2
mr.boskic.cat.c2.naive.0h = conductMRBoskic('C2','11:34432600','CD4_Naive','_0h_','ENSG00000121691') #ivw sig
mr.boskic.napsa.c2.naive.40h = conductMRBoskic('C2','19:50362278','CD4_Naive','_40h_','ENSG00000131400') #ivw sig
mr.boskic.cat.c2.memory.5d = conductMRBoskic('C2','11:34432600','CD4_Memory','_5d_','ENSG00000121691') #ivw sig
mr.boskic.adam15.c2.tcm.5d = conductMRBoskic('C2','1:155162859','TCM_','_5d_','ENSG00000143537') #ivw sig
mr.boskic.cat.c2.tcm.5d = conductMRBoskic('C2','11:34432600','TCM_','_5d_','ENSG00000121691') #ivw sig
mr.boskic.cat.c2.tn.0h = conductMRBoskic('C2','11:34432600','TN_0h','_0h_','ENSG00000121691') #ivw sig
mr.boskic.napsa.c2.tn.40h = conductMRBoskic('C2','19:50362278','TN_40h','_40h_','ENSG00000131400') #ivw sig
mr.boskic.ralgds.c2.tem.16h = conductMRBoskic('C2','9:133271745','TEM_16h','_16h_','ENSG00000160271') #ivw sig
mr.boskic.napsa.c2.tem.40h = conductMRBoskic('C2','19:50362278','TEM_40h','_40h_','ENSG00000131400') #ivw sig

#BQC MR
#a2
mr.bqc.napsa.a2.noninf = conductMRBQC('eQTL','A2','19:50363849','noninf','ENSG00000131400') #ivw sig
#bqc b2
mr.bqc.rab2a.b2.noninf = conductMRBQC('eQTL','B2','8:60515641','noninf','ENSG00000104388') #not ivw sig
mr.bqc.rab2a.b2.inf = conductMRBQC('eQTL','B2','8:60515641','inf','ENSG00000104388') #IVW sig, sens tests suggest issue
mr.bqc.napsa.b2.noninf = conductMRBQC('eQTL','B2','19:50362278','noninf','ENSG00000131400') #IVW sig
#bqc c2
mr.bqc.napsa.b2.noninf = conductMRBQC('eQTL','C2','19:50362278','noninf','ENSG00000131400') #not IVW sig

#GTEx MR
mr.gtex.napsa.a2 = conductMRGTEx('eQTL','A2','19:50363849',NA,'ENSG00000131400') #concerning sens
mr.gtex.napsa.b2 = conductMRGTEx('eQTL','B2','19:50362278',NA,'ENSG00000131400')
mr.gtex.napsa.c2 = conductMRGTEx('eQTL','C2','19:50362278',NA,'ENSG00000131400')

#plot MR results
a2.soskic = plotMRResults(list(mr.boskic.cat.a2.naive.0h,mr.boskic.napsa.a2.naive.40h,
                               mr.boskic.cat.a2.memory.5d,mr.boskic.acsf3.a2.memory.5d,
                               mr.boskic.ralgds.a2.TEM.16h,mr.boskic.napsa.a2.TEM.40h,
                               mr.boskic.cat.a2.TCM.5d,mr.boskic.cat.a2.TN.0h,
                               mr.boskic.napsa.a2.TN.40h),
                          c(1,0,1,1,1,1,1,1,1),outcome='A2')
b2.soskic = plotMRResults(list(mr.boskic.rab2a.b2.naive.16h,mr.boskic.rab2a.b2.naive.40h,
                               mr.boskic.cat.b2.naive.0h,mr.boskic.napsa.b2.naive.40h,
                               mr.boskic.rab2a.b2.memory.16h,mr.boskic.rab2a.b2.memory.40h,
                               mr.boskic.cat.b2.memory.5d,mr.boskic.rab2a.b2.tcm.16h,
                               mr.boskic.cat.b2.tcm.5d,mr.boskic.rab2a.b2.tn.16h,
                               mr.boskic.rab2a.b2.tn.40h,mr.boskic.cat.b2.tn.0h,
                               mr.boskic.napsa.b2.tn.40h,mr.boskic.ralgds.b2.tem.16h,
                               mr.boskic.napsa.b2.tem.40h),
                          c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),'B2')
c2.soskic = plotMRResults(list(mr.boskic.cat.c2.naive.0h,mr.boskic.napsa.c2.naive.40h,
                               mr.boskic.cat.c2.memory.5d,mr.boskic.adam15.c2.tcm.5d,
                               mr.boskic.cat.c2.tcm.5d,mr.boskic.cat.c2.tn.0h,
                               mr.boskic.napsa.c2.tn.40h,mr.boskic.ralgds.c2.tem.16h,
                               mr.boskic.napsa.c2.tem.40h),
                          c(0,1,0,1,0,0,1,1,1),'C2')

#plot MR for BQC data
bqc.mr.plot = plotMRResults.BQC(list(mr.bqc.napsa.a2.noninf,mr.bqc.rab2a.b2.noninf,
                                 mr.bqc.rab2a.b2.inf,mr.bqc.napsa.b2.noninf,
                                 mr.bqc.napsa.b2.noninf),c(1,0,1,0,0),
                            c('A2','B2','B2','B2','C2'))

#plot MR for GTEx data
bqc.mr.plot = plotMRResults.GTEx(list(mr.gtex.napsa.a2,mr.gtex.napsa.b2,
                                     mr.gtex.napsa.c2),c(0,1,1),
                                c('A2','B2','C2'))
