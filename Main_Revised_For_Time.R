library(tidyverse) ; library(stringr) ; library(magrittr) ; library(ggplot2)
library(TwoSampleMR) ; library(glue) ; library(vroom) ; library(openxlsx)
# library(ggpubr) ; library(gridExtra)
# library(edgeR) ; library(coloc) ; library(ieugwasr)
setwd("/scratch/richards/julian.willett/7.eQTLpQTLDisconnect")
source('Functions_Revised_for_Time.R') ; source('dynamic_QTL_coloc_minorfunctions.R')
`%notin%` = Negate(`%in%`)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# remember to give 60 GB Ram for clumping, 30 Gb for MR
args = commandArgs(trailingOnly=TRUE)
curr.index = as.numeric(args[[1]])
stop()

#first, need to make smaller files that we can do MR/coloc on to speed things up
# produceGTExSubpositionData() #grch38, done, slow
# produceSingleCellWholeBloodSubpositionData() # slow
# produceHGISubpositionData() #grch38, done, fast

order.index=c(1,87) #to facilitate running time-intensive steps in parallel
# getPreHarmonizedData('GTEx',order.index) #done
# getPreHarmonizedData('BQC19',c(curr.index,curr.index))
# getPreHarmonizedData('SC_Brain',order.index) #done

# doMR('GTEx',order.index) #done
# doMR('BQC19',c(curr.index,curr.index))
# doMR('SC_Brain',order.index) #done

gtex.coloc = doColoc('GTEx',NA)# ; makeColocSensitivityPlots(gtex.coloc,'gtex')
# scbrain.coloc = doColoc('SC_Brain',NA) #; makeColocSensitivityPlots(scbrain.coloc,'scbrain')
bqc.coloc = doColoc('BQC19',NA) #; makeColocSensitivityPlots(bqc.coloc,'bqc')

bqc.hits = getFinalHits('BQC',bqc.coloc,to.cut=c(1:4,7:8,10))
gtex.hits = getFinalHits('GTEx',gtex.coloc,to.cut=c(2:5,8:15,20:21,23,25,27:31,33:37,
                                             40:42,45:61,64:117,121,126:129,
                                             131,137:142,145,148,150:162,164,170,
                                             174:181,183,190:196,198,200,205:211,
                                             215:217,220:226,231:232,234:243,246,
                                             248,251:253,255:256,260:265,275,277:286,
                                             288:308,312:316,319,323:326,331:333,
                                             336:338,340:342,344:345,350:351,353:354,
                                             356:357,360:374,376:377))
soskic.hits = readRDS('MRthenColoc_MRData/coloc_sens_all.rds')[[3]][-c(1,4,5,6,13,15),]
# scbrain.hits = getFinalHits('SC_Brain',scbrain.coloc,to.cut=c(2:6,8:9,13:18,22:26,28:29))

# num that passed MR sens AND colocalized confidently
(nrow(bqc.hits)+nrow(gtex.hits)+nrow(soskic.hits)) / 5015

mr.plot.bulk = plotNewMR('Bulk',NA,NA,bqc.hits,gtex.hits,16)
mr.plot.sc =  plotNewMR('SC',soskic.hits,NA,NA,NA,20)

plotNewGTExHeatmap(gtex.hits)
# plotNewGTExHeatmap(scbrain.hits)

####################
# Make percent colocalized fig
makePercentColocFig() # fig 3


# Produce locus zoom datasets
makeDataLocusZoom('IFNAR2')
makeDataLocusZoom('RALGDS')
makeDataLocusZoom('IL10RB')

# loci for cutting hla regions that were included in original analysis
length(which(!str_detect(bqc.mr.data$Label,'chr6_29946723') & 
               !str_detect(bqc.mr.data$Label,'chr6_31182658') & 
               !str_detect(bqc.mr.data$Label,'chr6_31276554') & 
               !str_detect(bqc.mr.data$Label,'chr6_31306250') & 
               !str_detect(bqc.mr.data$Label,'chr6_32707868')))

makePercentViolateSCVFig()

makeSupplementalTables()