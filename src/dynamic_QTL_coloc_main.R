library(tidyverse) ; library(glue) ; library(magrittr)
library(vroom) ; library(arrow) ; library(coloc)
library(SomaDataIO) ; library(TwoSampleMR) ; library(MRutils)
library(data.table)
source('dynamic_QTL_coloc_functions.R')
source('dynamic_QTL_coloc_minorfunctions.R')
source('dynamic_QTL_MRthenColoc_functions.R')

setwd("/scratch/richards/julian.willett/eQTLpQTLDisconnect")
`%notin%` <- Negate(`%in%`)

#############################################
#CODE FOR MR THEN COLOCALIZATION
# getAllMRDataReady.MRthenColoc('soskic','eQTL') #time intensive step (several days)
# getAllMRDataReady.MRthenColoc('bqc','eQTL') #time intensive step (around a day)
# getAllMRDataReady.MRthenColoc('gtex','eQTL') # time intensive step (around 8 hr)
all.mr.data.soskic = conductMR.MRthenColoc('soskic','eQTL',rev=F) ; saveRDS(all.mr.data.soskic,'MRthenColoc_MRData/mr_results_all.rds')
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
