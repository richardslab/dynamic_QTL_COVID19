library(tidyverse) ; library(glue) ; library(magrittr)
library(vroom) ; library(arrow) ; library(coloc)
library(SomaDataIO) ; library(TwoSampleMR) ; library(MRutils)
library(data.table)
source('dynamic_eQTL_coloc_functions.R')
source('dynamic_eQTL_coloc_minorfunctions.R')

setwd("/scratch/richards/julian.willett/eQTLpQTLDisconnect")
`%notin%` <- Negate(`%in%`)

#coloc work for soskic data
# soskic.major.cells.coloc.rel7 = soskicColoc(vector.specific=NULL,allCells=F,returnCommonQTLs=F)
soskic.major.cells.coloc.rel7 = readRDS('soskic_rel7_majorcells_allcoloc.rds')
# soskic.all.cells.coloc.rel7 = soskicColoc(vector.specific=NULL,allCells=T,returnCommonQTLs=F)
# saveRDS(soskic.all.cells.coloc.rel7,'soskic_rel7_allcells_allcoloc.rds')

#coloc work for bqc19
# bqc.coloc.rel7.eQTL = bqcColoc(type.of.qtl='eQTL',vector.specific=NULL,rel6=F,returnCommonQTLs=F,mr=F)
# saveRDS(bqc.coloc.rel7.eQTL,'bqc_rel7_eQTL_colocall.rds')
# bqc.coloc.rel7.sQTL = bqcColoc(type.of.qtl='sQTL',vector.specific=NULL,rel6=F,mr=F)
# saveRDS(bqc.coloc.rel7.sQTL,'bqc_rel7_sQTL_colocall.rds')
bqc.coloc.rel7.eQTL = readRDS('bqc_rel7_eQTL_colocall.rds')
bqc.coloc.rel7.sQTL = readRDS('bqc_rel7_sQTL_colocall.rds')

#coloc work for gtex at same sig loci
# gtex.coloc.eQTL = gtexColoc('eQTL',NULL,returnCommonQTLs=F,mr=F,onlySoskicLoci=F)
# saveRDS(gtex.coloc.eQTL,'gtex_rel8_eQTL_colocall.rds')
# gtex.coloc.sQTL = gtexColoc('sQTL',NULL,returnCommonQTLs=F,mr=F,onlySoskicLoci=F)
# saveRDS(gtex.coloc.sQTL,'gtex_rel8_sQTL_colocall.rds')
gtex.coloc.eQTL = readRDS('gtex_rel8_eQTL_colocall.rds')
gtex.coloc.sQTL = readRDS('gtex_rel8_sQTL_colocall.rds')

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

#BQC H4.PP FIGS sQTL
bqc.a2.coloc = plotDynamicPP.BQC(bqc.coloc.rel7.sQTL,'A2','sQTL',c(12,21),c('chr12:112911235:112916509:clu_9748_+:ENSG00000089127.15','chr12:112911258:112916509:clu_9748_+:ENSG00000089127.15','chr21:33230216:33241840:clu_28170_+:IFNAR2','chr21:33230216:33241886:clu_28170_+:IFNAR2'))
bqc.b2.coloc = plotDynamicPP.BQC(bqc.coloc.rel7.sQTL,'B2','sQTL',c(8,12,21),c('chr12:112911235:112916509:clu_9748_+:ENSG00000089127.15'))
bqc.c2.coloc = plotDynamicPP.BQC(bqc.coloc.rel7.sQTL,'C2','sQTL',c(12,21),c('chr12:112911235:112916509:clu_9748_+:ENSG00000089127.15','chr21:33230216:33241840:clu_28170_+:IFNAR2','OAS1 Splice Variant G'))
rm(bqc.a2.coloc,bqc.b2.coloc,bqc.c2.coloc)

#GTEX H4.PP FIGS eQTL
gtex.a2.coloc = plotDynamicPP.GTEx(gtex.coloc.eQTL,'A2','eQTL',c(5,7,9,17,19),exclude.genes=c('LINC02863','KANSL1-AS1'))
gtex.b2.coloc = plotDynamicPP.GTEx(gtex.coloc.eQTL,'B2','eQTL',c(5,9,17,19),c('LINC02863','KANSL1-AS1'))
gtex.c2.coloc = plotDynamicPP.GTEx(gtex.coloc.eQTL,'C2','eQTL',c(1,12,19),c('NTN5'))
rm(gtex.a2.coloc,gtex.b2.coloc,gtex.c2.coloc)

#GTEX H4.PP FIGS sQTL
gtex.a2.coloc = plotDynamicPP.GTEx(gtex.coloc.sQTL,'A2','sQTL',c(12,19,21),c('chr19:10351162:10352434:clu_19126:TYK2','chr21:33230216:33241840:clu_24781:ENSG00000159110.19','chr21:33252830:33260597:clu_24783:ENSG00000159110.19','IFNAR2 Splice Variant E'))
gtex.b2.coloc = plotDynamicPP.GTEx(gtex.coloc.sQTL,'B2','sQTL',c(12,19,21),c('chr19:10351162:10352434:clu_19126:TYK2','chr21:33230216:33241840:clu_24781:ENSG00000159110.19','chr21:33252830:33260597:clu_24783:ENSG00000159110.19','IFNAR2 Splice Variant E'))
gtex.c2.coloc = plotDynamicPP.GTEx(gtex.coloc.sQTL,'C2','sQTL',c(1,12,19,21),c('chr19:10353098:10354042:clu_19127:TYK2','chr21:33230216:33241840:clu_24781:ENSG00000159110.19','chr21:33252830:33260597:clu_24783:ENSG00000159110.19','IFNAR2 Splice Variant E'))

###sens testing
soskic.eqtl.sens = doAllSensTests('soskic','eQTL',soskic.major.cells.coloc.rel7) ; saveRDS(soskic.eqtl.sens,'soskic_eqtl_sens.rds')
# bqc.eqtl.sens = doAllSensTests('bqc','eQTL',bqc.coloc.rel7.eQTL) ; saveRDS(bqc.eqtl.sens,'bqc_eqtl_sens.rds')
# bqc.sqtl.sens = doAllSensTests('bqc','sQTL',bqc.coloc.rel7.sQTL) ; saveRDS(bqc.sqtl.sens,'bqc_sqtl_sens.rds')
gtex.eqtl.sens = doAllSensTests('gtex','eQTL',gtex.coloc.eQTL) ; saveRDS(gtex.eqtl.sens,'gtex_eqtl_sens.rds')
# gtex.sqtl.sens = doAllSensTests('gtex','sQTL',gtex.coloc.sQTL) ; saveRDS(gtex.sqtl.sens,'gtex_sqtl_sens.rds')
# item=48
# sensitivity(gtex.sqtl.sens[[1]][[item]],rule='H4 > 0.5')
# print(gtex.sqtl.sens[[2]][[item]])

#Get all data for MR input
mr.data.soskic.eqtl = getAllMRDataReady('soskic','eQTL',soskic.major.cells.coloc.rel7 %>% dplyr::filter(CHR != 6)) ; saveRDS(mr.data.soskic.eqtl,'mr_data_soskic_eqtl.rds')
# mr.data.bqc.eqtl = getAllMRDataReady('bqc','eQTL',bqc.coloc.rel7.eQTL %>% dplyr::filter(CHR != 6)) ; saveRDS(mr.data.bqc.eqtl,'mr_data_bqc_eqtl.rds')
# mr.data.bqc.sqtl = getAllMRDataReady('bqc','sQTL',bqc.coloc.rel7.sQTL %>% dplyr::filter(CHR != 6)) ; saveRDS(mr.data.bqc.sqtl,'mr_data_bqc_sqtl.rds')
# mr.data.gtex.eqtl = getAllMRDataReady('gtex','eQTL',gtex.coloc.eQTL %>% dplyr::filter(CHR != 6)) ; saveRDS(mr.data.gtex.eqtl,'mr_data_gtex_eqtl.rds')
# mr.data.gtex.sqtl = getAllMRDataReady('gtex','sQTL',gtex.coloc.sQTL %>% dplyr::filter(CHR != 6)) ; saveRDS(mr.data.gtex.sqtl,'mr_data_gtex_sqtl.rds')

#Do MR
mr.computed.boskic.eqtl = conductMR(mr.data.soskic.eqtl)
# mr.computed.bqc.eqtl = conductMR(mr.data.bqc.eqtl) ; saveRDS(mr.computed.bqc.eqtl,'mr_computed_bqc_eqtl.rds')
# mr.computed.bqc.sqtl = conductMR(mr.data.bqc.sqtl) ; saveRDS(mr.computed.bqc.sqtl,'mr_computed_bqc_sqtl.rds')
# mr.computed.gtex.eqtl = conductMR(mr.data.gtex.eqtl) ; saveRDS(mr.computed.gtex.eqtl,'mr_computed_gtex_eqtl.rds')
# mr.computed.gtex.sqtl = conductMR(mr.data.gtex.sqtl) ; saveRDS(mr.computed.gtex.sqtl,'mr_computed_gtex_sqtl.rds')

#plot MR results
mr.plot.soskic.eqtl.a2 = plotMRResults('soskic','eQTL',mr.computed.boskic.eqtl,'A2',33,c())
mr.plot.soskic.eqtl.b2 = plotMRResults('soskic','eQTL',mr.computed.boskic.eqtl,'B2',33,c())
mr.plot.soskic.eqtl.c2 = plotMRResults('soskic','eQTL',mr.computed.boskic.eqtl,'C2',33,c())

#plot MR for GTEx data
#Plot MR results eQTL
mr.plot.bqc.eqtl.a2 = plotMRResults('bqc','eQTL',mr.computed.bqc.eqtl,'A2',22,c('ENSG00000188199.10','ENSG00000120071.15','ENSG00000262539.1','ENSG00000073969.18'))
mr.plot.bqc.eqtl.b2 = plotMRResults('bqc','eQTL',mr.computed.bqc.eqtl,'B2',22,c('ENSG00000262539.1','ENSG00000073969.18','ENSG00000142233.14'))
mr.plot.bqc.eqtl.c2 = plotMRResults('bqc','eQTL',mr.computed.bqc.eqtl,'C2',22,c('ENSG00000142233.14'))
mr.plot.gtex.eqtl.a2 = plotMRResults('gtex','eQTL',mr.computed.gtex.eqtl,'A2',13,c('ENSG00000238160.1','ENSG00000214401.4'))
mr.plot.gtex.eqtl.b2 = plotMRResults('gtex','eQTL',mr.computed.gtex.eqtl,'B2',13,c('ENSG00000238160.1','ENSG00000214401.4'))
mr.plot.gtex.eqtl.c2 = plotMRResults('gtex','eQTL',mr.computed.gtex.eqtl,'C2',13,c('ENSG00000142233.11'))

#Plot MR results sQTL
mr.plot.bqc.sqtl.a2 = plotMRResults('bqc','sQTL',mr.computed.bqc.sqtl,'A2',39,c('chr12:112911258:112916509:clu_9748_+:ENSG00000089127.15','chr12:112911235:112916509:clu_9748_+:ENSG00000089127.15','chr17:46038686:46039702:clu_17314_-:ENSG00000120071.15','chr17:46094701:46170855:clu_18029_-:ENSG00000120071.15','chr19:48662041:48663463:clu_21139_-:ENSG00000142233.14','chr19:49670488:49673698:clu_22852_+:ENSG00000126453.10','chr21:33230216:33241840:clu_28170_+:ENSG00000159110.21','chr21:33230216:33241886:clu_28170_+:ENSG00000159110.21'))
mr.plot.bqc.sqtl.b2 = plotMRResults('bqc','sQTL',mr.computed.bqc.sqtl,'B2',39,c('chr9:132935095:132944543:clu_41810_-:ENSG00000165699.15','chr12:112911235:112916509:clu_9748_+:ENSG00000089127.15','chr17:46038686:46039702:clu_17314_-:ENSG00000120071.15','chr17:46094701:46170855:clu_18029_-:ENSG00000120071.15','chr19:48662041:48663463:clu_21139_-:ENSG00000142233.14'))
mr.plot.bqc.sqtl.c2 = plotMRResults('bqc','sQTL',mr.computed.bqc.sqtl,'C2',39,c('chr9:132935095:132944543:clu_41810_-:ENSG00000165699.15','chr12:112911258:112916509:clu_9748_+:ENSG00000089127.15','chr12:112911235:112916509:clu_9748_+:ENSG00000089127.15','chr19:48662041:48663463:clu_21139_-:ENSG00000142233.14','chr21:33230216:33241840:clu_28170_+:ENSG00000159110.21'))
mr.plot.gtex.sqtl.a2 = plotMRResults('gtex','sQTL',mr.computed.gtex.sqtl,'A2',16,c('chr1:155732033:155736946:clu_43109:ENSG00000132676.15','chr1:155737063:155738157:clu_43109:ENSG00000132676.15','chr5:131988226:131988805:clu_31067:ENSG00000164398.12','chr10:79808832:79814613:clu_8678:ENSG00000225484.6','chr10:79766324:79804485:clu_8677:ENSG00000225484.6','chr17:46172232:46192823:clu_11356:ENSG00000120071.13','chr17:46094701:46170855:clu_11354:ENSG00000120071.13','chr17:46038686:46039702:clu_11351:ENSG00000120071.13','chr19:10351162:10352434:clu_19126:ENSG00000105397.13','chr21:33230216:33241840:clu_24781:ENSG00000159110.19','chr21:33230216:33241886:clu_24781:ENSG00000159110.19','chr21:33252830:33260597:clu_24783:ENSG00000159110.19'))
mr.plot.gtex.sqtl.b2 = plotMRResults('gtex','sQTL',mr.computed.gtex.sqtl,'B2',16,c('chr1:155732033:155736946:clu_43109:ENSG00000132676.15','chr1:155737063:155738157:clu_43109:ENSG00000132676.15','chr17:46094701:46170855:clu_11354:ENSG00000120071.13','chr17:46038686:46039702:clu_11351:ENSG00000120071.13','chr19:10351162:10352434:clu_19126:ENSG00000105397.13','chr21:33230216:33241840:clu_24781:ENSG00000159110.19','chr21:33230216:33241886:clu_24781:ENSG00000159110.19','chr21:33252830:33260597:clu_24783:ENSG00000159110.19'))
mr.plot.gtex.sqtl.c2 = plotMRResults('gtex','sQTL',mr.computed.gtex.sqtl,'C2',16,c('chr19:10353098:10354042:clu_19127:ENSG00000105397.13','chr21:33230216:33241840:clu_24781:ENSG00000159110.19','chr21:33230216:33241886:clu_24781:ENSG00000159110.19','chr21:33252830:33260597:clu_24783:ENSG00000159110.19'))

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
                                             'B2 OAS1','C2 ADAM15','C2 OAS1','C2 TYK2'))
