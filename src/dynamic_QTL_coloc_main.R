library(tidyverse) ; library(glue) ; library(magrittr)
library(vroom) ; library(arrow) ; library(coloc)
library(SomaDataIO) ; library(TwoSampleMR) ; library(MRutils)
library(data.table)
source('dynamic_QTL_coloc_functions.R')
source('dynamic_QTL_coloc_minorfunctions.R')
source('dynamic_QTL_MRthenColoc_functions.R')

setwd("/scratch/richards/julian.willett/eQTLpQTLDisconnect")
`%notin%` <- Negate(`%in%`)
tok = 'a9d29e4a815d'

#############################################
#CODE FOR MR THEN COLOCALIZATION
# getAllMRDataReady.MRthenColoc('soskic','eQTL') #time intensive step (several days)
# getAllMRDataReady.MRthenColoc('bqc','eQTL') #time intensive step (around a day)
# getAllMRDataReady.MRthenColoc('gtex','eQTL') # time intensive step (around 8 hr)
all.mr.data.soskic = conductMR.MRthenColoc('soskic','eQTL',rev=F,recollect=F) ; saveRDS(all.mr.data.soskic,'MRthenColoc_MRData/mr_results_all.rds')
all.mr.data.bqc.eqtl = conductMR.MRthenColoc('bqc','eQTL',rev=F) ; saveRDS(all.mr.data.bqc.eqtl,'MRthenColoc_MRData_BQC_eQTL/mr_results_all.rds')
all.mr.data.gtex.eqtl = conductMR.MRthenColoc('gtex','eQTL',rev=F) ; saveRDS(all.mr.data.gtex.eqtl,'MRthenColoc_MRData_GTEx_eQTL/mr_results_all.rds')

soskic.coloc.mrthencoloc = colocFromMR('soskic','eQTL',all.mr.data.soskic,sens.data=NULL) ; saveRDS(soskic.coloc.mrthencoloc,'MRthenColoc_MRData/coloc_results_all.rds')
bqc.coloc.mrthencoloc = colocFromMR('bqc','eQTL',all.mr.data.bqc.eqtl,sens.data = NULL) ; saveRDS(bqc.coloc.mrthencoloc,'MRthenColoc_MRData_BQC_eQTL/coloc_results_all.rds')
gtex.coloc.mrthencoloc = colocFromMR('gtex','eQTL',all.mr.data.gtex.eqtl,sens.data = NULL) ; saveRDS(gtex.coloc.mrthencoloc,'MRthenColoc_MRData_GTEx_eQTL/coloc_results_all.rds')

soskic.coloc.mrthencoloc.sens = colocFromMR('soskic','eQTL',all.mr.data.soskic,sens.data = soskic.coloc.mrthencoloc) ; saveRDS(soskic.coloc.mrthencoloc.sens,'MRthenColoc_MRData/coloc_sens_all.rds')
bqc.coloc.mrthencoloc.sens = colocFromMR('bqc','eQTL',all.mr.data.bqc.eqtl,sens.data = bqc.coloc.mrthencoloc) ; saveRDS(bqc.coloc.mrthencoloc.sens,'MRthenColoc_MRData_BQC_eQTL/coloc_sens_all.rds')
gtex.coloc.mrthencoloc.sens = colocFromMR('gtex','eQTL',all.mr.data.gtex.eqtl,sens.data = gtex.coloc.mrthencoloc) ; saveRDS(gtex.coloc.mrthencoloc.sens,'MRthenColoc_MRData_GTEx_eQTL/coloc_sens_all.rds')

doSensTesting(soskic.coloc.mrthencoloc,soskic.coloc.mrthencoloc.sens[[1]],1)
doSensTesting(bqc.coloc.mrthencoloc,bqc.coloc.mrthencoloc.sens[[1]],1)
doSensTesting(gtex.coloc.mrthencoloc,gtex.coloc.mrthencoloc.sens[[1]],1)

soskic.mr.plot = plotSignificantMR(data.source='soskic',outcome='A2',
                  mr.coloc.data=soskic.coloc.mrthencoloc.sens[[3]],
                  drop.indices=c(1,4,5,6,13,15))
bqc.mr.plot = plotSignificantMR(data.source='bqc',outcome=c('A2','B2','C2'),
                                   mr.coloc.data=bqc.coloc.mrthencoloc.sens[[3]],
                                   drop.indices=c(1,5))
gtex.mr.plot = plotSignificantMR(data.source='gtex',outcome=c('A2','B2','C2'),
                                mr.coloc.data=gtex.coloc.mrthencoloc.sens[[3]],
                                drop.indices=c(1,3,4))

#make heatmap for coloc
# heatmap.soskic = getHeatmapData(data.source='soskic',
#                              data=soskic.coloc.mrthencoloc.sens[[3]][-c(1,4,5,6,13,15),])
# saveRDS(heatmap.soskic,'MRthenColoc_MRData/heatmap_data.rds')
# heatmap.bqc = getHeatmapData(data.source='bqc',
#                                 data=bqc.coloc.mrthencoloc.sens[[3]][-c(1,5),])
# saveRDS(heatmap.bqc,'MRthenColoc_MRData_BQC_eQTL/heatmap_data.rds')

soskic.heatmap.plot = plotColocHeatmap('soskic',soskic.mr.plot,coloc.data=heatmap.soskic)
bqc.heatmap.plot = plotColocHeatmap('bqc',bqc.mr.plot,coloc.data=heatmap.bqc)

#############################################
#Downstream analysis

#logistic regression
# il10rb.a2 = doLogisticAnalysis('A2','IL10RB')
# oas1.a2 = doLogisticAnalysis('A2','OAS1')
# rab2a.b2 = doLogisticAnalysis('B2','RAB2A')
# abo.b2 = doLogisticAnalysis('B2','ABO')
# oas1.b2 = doLogisticAnalysis('B2','OAS1')
# adam15.c2 = doLogisticAnalysis('C2','ADAM15')
# oas1.c2 = doLogisticAnalysis('C2','OAS1')
# tyk2.c2 = doLogisticAnalysis('C2','TYK2')
# protein.analysis = list(il10rb.a2,oas1.a2,rab2a.b2,abo.b2,oas1.b2,adam15.c2,oas1.c2,tyk2.c2)
saveRDS(protein.analysis,'somalogic_protein_analysis.rds')
fig.data = plotOddsRatios(protein.analysis,c('A2 IL10RB','A2 OAS1','B2 RAB2A','B2 ABO',
                                             'B2 OAS1','C2 ADAM15','C2 OAS1','C2 TYK2'),
                          c(3,4,6))

#############################################
#Code for doing colocalization for all transcripts

#coloc work for soskic data
soskic.major.cells.coloc.rel7 = soskicColoc(vector.specific=NULL,allCells=F,returnCommonQTLs=F,mr=F,specific.getdf=F)
saveRDS(soskic.major.cells.coloc.rel7,'soskic_rel7_majorcells_allcoloc.rds')
soskic.major.cells.coloc.rel7 = readRDS('soskic_rel7_majorcells_allcoloc.rds')
# soskic.all.cells.coloc.rel7 = soskicColoc(vector.specific=NULL,allCells=T,returnCommonQTLs=F)
# saveRDS(soskic.all.cells.coloc.rel7,'soskic_rel7_allcells_allcoloc.rds')

#coloc work for bqc19
# bqc.coloc.rel7.eQTL = bqcColoc(type.of.qtl='eQTL',vector.specific=NULL,returnCommonQTLs=F,mr=F,specific.getdf=F)
# saveRDS(bqc.coloc.rel7.eQTL,'bqc_rel7_eQTL_colocall.rds')
bqc.coloc.rel7.eQTL = readRDS('bqc_rel7_eQTL_colocall.rds')

#coloc work for gtex at same sig loci
# gtex.coloc.eQTL = gtexColoc('eQTL',NULL,returnCommonQTLs=F,mr=F)
# saveRDS(gtex.coloc.eQTL,'gtex_rel8_eQTL_colocall.rds')
gtex.coloc.eQTL = readRDS('gtex_rel8_eQTL_colocall.rds')

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

#BQC H4.PP FIGS eQTL
bqc.a2.coloc = plotDynamicPP.BQC(bqc.coloc.rel7.eQTL,'A2','eQTL',c(7,9,10,17,19,21),c('NSF','KANSL1','NUTM2B','RP11-259G18.3'))
bqc.b2.coloc = plotDynamicPP.BQC(bqc.coloc.rel7.eQTL,'B2','eQTL',c(8,9,17,19,21),c('RP11-259G18.3','NSF','NTN5'))
bqc.c2.coloc = plotDynamicPP.BQC(bqc.coloc.rel7.eQTL,'C2','eQTL',c(1,19,21),c('NTN5'))
rm(bqc.a2.coloc,bqc.b2.coloc,bqc.c2.coloc)

#GTEX H4.PP FIGS eQTL
gtex.a2.coloc = plotDynamicPP.GTEx(gtex.coloc.eQTL,'A2','eQTL',c(5,7,9,17,19),exclude.genes=c('LINC02863','KANSL1-AS1'))
gtex.b2.coloc = plotDynamicPP.GTEx(gtex.coloc.eQTL,'B2','eQTL',c(5,9,17,19),c('LINC02863','KANSL1-AS1'))
gtex.c2.coloc = plotDynamicPP.GTEx(gtex.coloc.eQTL,'C2','eQTL',c(1,12,19),c('NTN5'))
rm(gtex.a2.coloc,gtex.b2.coloc,gtex.c2.coloc)

###sens testing
soskic.eqtl.sens = doAllSensTests('soskic','eQTL',soskic.major.cells.coloc.rel7) ; saveRDS(soskic.eqtl.sens,'soskic_eqtl_sens.rds')
# bqc.eqtl.sens = doAllSensTests('bqc','eQTL',bqc.coloc.rel7.eQTL) ; saveRDS(bqc.eqtl.sens,'bqc_eqtl_sens.rds')
# bqc.sqtl.sens = doAllSensTests('bqc','sQTL',bqc.coloc.rel7.sQTL) ; saveRDS(bqc.sqtl.sens,'bqc_sqtl_sens.rds')
gtex.eqtl.sens = doAllSensTests('gtex','eQTL',gtex.coloc.eQTL) ; saveRDS(gtex.eqtl.sens,'gtex_eqtl_sens.rds')
# gtex.sqtl.sens = doAllSensTests('gtex','sQTL',gtex.coloc.sQTL) ; saveRDS(gtex.sqtl.sens,'gtex_sqtl_sens.rds')
# item=48
# sensitivity(gtex.sqtl.sens[[1]][[item]],rule='H4 > 0.5')
# print(gtex.sqtl.sens[[2]][[item]])
