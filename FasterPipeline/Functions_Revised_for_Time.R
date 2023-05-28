produceGTExSubpositionData = function(index) { #gtex grch38
  loci = vroom('subpositions_chrpos.txt',col_names=F) %>%
    dplyr::arrange(X1) %>% dplyr::rename(Chr=X1,Pos=X2,Outcome=X3)
  
  files = list.files('/scratch/richards/satoshi.yoshiji/database/GTEx_v8/eQTL',
                     recursive=T,full.names=T,all.files=T,include.dirs=T,
                     pattern='tsv')
  data.list = list()
  for (loc in 1:nrow(loci)) {
    if (loc != index) next
    curr.chr = loci$Chr[[loc]]
    print(glue('On locus {loc} of {nrow(loci)}. Time: {Sys.time()}'))
    curr.files = files[which(str_detect(files,glue('\\.{loci$Chr[[loc]]}.eQTL')))]
    file.names = gsub('.*\\/','',curr.files) %>% gsub('\\..*','',.)
    
    #limit repeat time-intensive operations.
    for (index in 1:length(curr.tissue.indices)) {
      print(glue('Reading in file {index} of {length(curr.tissue.indices)}'))
      data.list[[index]] = vroom(curr.files[[index]],show_col_types=F)
      gc()
    }
    
    for (df in 1:length(data.list)) {
      print(glue('Subsetting df {df} of {length(curr.tissue.indices)}. Time: {Sys.time()}'))
      vroom_write(data.list[[df]] %>% dplyr::filter(POS >= loci$Pos[[loc]]-500000,
                                                    POS <= loci$Pos[[loc]]+500000),
                  glue('Fast_GTEx_eQTL/{file.names[[df]]}_chr{loci$Chr[[loc]]}_{loci$Pos[[loc]]}.txt'))
      gc()
    }
  }
}
produceHGISubpositionData = function(window,index) {
  loci = vroom('subpositions_chrpos.txt',col_names=F) %>%
    dplyr::arrange(X3) %>% dplyr::rename(Chr=X1,Pos=X2,Outcome=X3)
  
  hgi.a2 = vroom('HGI_Data/COVID19_HGI_A2_ALL_eur_leave23andme_20220403.tsv.gz',show_col_types = F)
  hgi.b2 = vroom('HGI_Data/COVID19_HGI_B2_ALL_eur_leave23andme_20220403.tsv.gz',show_col_types = F)
  hgi.c2 = vroom('HGI_Data/COVID19_HGI_C2_ALL_eur_leave23andme_20220403.tsv.gz',show_col_types = F)
  
  for (l in 1:nrow(loci)) {
    print(glue('On locus {l} of {nrow(loci)}'))
    if (l %notin% index) next
    if (file.exists(glue('Fast_HGI_Outcomes_{window}/{loci$Chr[[l]]}_{loci$Pos[[l]]}_{loci$Outcome[[l]]}.txt'))) 
      next
    if (loci$Outcome[[l]] == 'A2') sub = hgi.a2 
    else if (loci$Outcome[[l]] == 'B2') sub = hgi.b2 
    else if (loci$Outcome[[l]] == 'C2') sub = hgi.c2 
    
    sub %<>% dplyr::filter(`#CHR`==loci$Chr[[l]],POS>=loci$Pos[[l]]-window,
                           POS<=loci$Pos[[l]]+window)
    vroom_write(sub,glue('Fast_HGI_Outcomes_{window}/{loci$Chr[[l]]}_{loci$Pos[[l]]}_{loci$Outcome[[l]]}.txt'))
  }
  gc()
}
getPreHarmonizedData = function(dataset,oi) {
  # harmonizing takes a lot of memory when vectorized (crashing R),
  # and clumping takes time, so separating the processes
  loci = vroom('subpositions_chrpos.txt',col_names=F) %>%
    dplyr::arrange(X3) %>% dplyr::rename(Chr=X1,Pos=X2,Outcome=X3)
  
  if (dataset == 'Soskic') { #run in dynamic_QTL_MRthenColoc_functions.R, uses 500 kb window
    
  }else if (dataset == 'BQC19') {
    for (l in oi[[1]]:oi[[2]]) { #locus index
      print(glue('On locus {l} of {nrow(loci)}. Time: {Sys.time()}'))
      
      #first load in the data and split as necessary by gene
      outcome.df = vroom(glue('Fast_HGI_Outcomes_500000/{loci$Chr[[l]]}_{loci$Pos[[l]]}_{loci$Outcome[[l]]}.txt'),show_col_types=F) %>%
        dplyr::mutate(SNP=glue('chr{SNP}'))
      curr.loc = loci[l,]
      curr.files = (list.files('BQC_eQTL_Data/Subpositions',full.names=T))
      curr.files = curr.files[which(str_detect(curr.files,glue('chr{curr.loc$Chr}_{curr.loc$Pos}')))]
      
      for (f in curr.files) { #inf vs noninf
        print(glue('On file {f}. Locus {l}'))
        f.name = gsub('.*\\/','',f) %>% gsub('\\.txt','',.)
        harmonize.saved.names = paste0(f.name,'_',curr.loc$Outcome,'.rds')
        print(harmonize.saved.names)
        if (harmonize.saved.names %in% list.files('Fast_BQC_eQTL/pre_harmonized')) next
        
        exp.out.df = vroom(f,show_col_types=F)
        if (nrow(exp.out.df)<1) next
        exp.out.df %<>% dplyr::rename(SNP=variant_id) %>% 
          dplyr::filter(SNP %in% outcome.df$SNP,abs(tss_distance)<=500000) %>%
          merge(.,outcome.df,by='SNP') %>% drop_na(rsid) %>% 
          dplyr::mutate(SNP=rsid,P_exp=pchisq((slope/slope_se)^2,df=1,lower.tail=F),
                        Outcome_Total_N=all_inv_var_meta_cases + all_inv_var_meta_controls)
        split.by.gene = split(exp.out.df,f=exp.out.df$gene_id,drop=T) 
        split.by.gene = split.by.gene[sapply(split.by.gene, nrow)>=50]
        
        #get mr data ready too
        #exposure
        curr.exp = lapply(split.by.gene,function(x) {
          format_data(x,type='exposure',chr_col='CHR',pos_col='POS.x',
                      beta_col='slope',se_col='slope_se',snp_col='rsid',
                      effect_allele_col='ALT.x',other_allele_col='REF.x',
                      eaf_col='maf',samplesize_col='ma_samples',pval_col='P_exp')})
        
        curr.exp = lapply(curr.exp,function(x) { 
          x %<>% dplyr::rename(rsid=SNP,pval=pval.exposure)
          x = ld_clump(dat=x,plink_bin='/home/richards/julian.willett/plink',
                       bfile='/scratch/richards/satoshi.yoshiji/database/1KG/EUR/EUR')
          x %<>% dplyr::rename(SNP=rsid,pval.exposure=pval)
          gc()
          x
        })
        gc()
        
        #outcome
        curr.out = lapply(split.by.gene,function(x) {
          format_data(x,type='outcome',chr_col='CHR',pos_col='POS.y',
                      beta_col='all_inv_var_meta_beta',se_col='all_inv_var_meta_sebeta',
                      snp_col='rsid',
                      effect_allele_col='ALT.y',other_allele_col='REF.y',
                      eaf_col='all_meta_AF',samplesize_col='Outcome_Total_N',
                      pval_col='all_inv_var_meta_p')})
        gc()
        exp.out.list = lapply(1:length(curr.exp),function(x) {
          list(curr.exp[[x]],curr.out[[x]])
        })
        names(exp.out.list) = names(split.by.gene)
        saveRDS(exp.out.list,glue('Fast_BQC_eQTL/pre_harmonized/{f.name}_{curr.loc$Outcome}.rds'))
      }
      gc()
    }
  }else if (dataset == 'GTEx') {
    for (l in oi[[1]]:oi[[2]]) {
      print(glue('On locus {l} of {nrow(loci)}. Time: {Sys.time()}'))
      
      #first load in the data and split as necessary by gene
      outcome.df = vroom(glue('Fast_HGI_Outcomes/{loci$Chr[[l]]}_{loci$Pos[[l]]}_{loci$Outcome[[l]]}.txt'),show_col_types=F) %>%
        dplyr::mutate(SNP=str_replace_all(SNP,':','_')) %>%
        dplyr::mutate(SNP=glue('chr{SNP}_b38'))
      curr.loc = loci[l,]
      curr.files = (list.files('Fast_GTEx_eQTL',full.names=T))
      curr.files = curr.files[which(str_detect(curr.files,glue('chr{curr.loc$Chr}_{curr.loc$Pos}')))]

      for (f in curr.files) { 
        print(glue('On file {f}. Locus {l}'))
        f.name = gsub('.*\\/','',f) %>% gsub('\\.txt','',.)
        harmonize.saved.names = paste0(f.name,'_',curr.loc$Outcome,'.rds')
        print(harmonize.saved.names)
        if (harmonize.saved.names %in% list.files('Fast_GTEx_eQTL/pre_harmonized')) next
        
        exp.out.df = vroom(f,show_col_types=F)
        if (nrow(exp.out.df)<1) next
        exp.out.df %<>% dplyr::rename(SNP=variant_id) %>% 
          dplyr::filter(SNP %in% outcome.df$SNP,abs(tss_distance)<=500000) %>%
          merge(.,outcome.df,by='SNP') %>% drop_na(rsid) %>% 
          dplyr::mutate(SNP=rsid,P_exp=pchisq((slope/slope_se)^2,df=1,lower.tail=F),
                        Outcome_Total_N=all_inv_var_meta_cases + all_inv_var_meta_controls)
        split.by.gene = split(exp.out.df,f=exp.out.df$phenotype_id,drop=T) 
        split.by.gene = split.by.gene[sapply(split.by.gene, nrow)>=50]
        
        #get mr data ready too
        #exposure
        curr.exp = lapply(split.by.gene,function(x) {
          format_data(x,type='exposure',chr_col='CHR',pos_col='POS.x',
                      beta_col='slope',se_col='slope_se',snp_col='rsid',
                      effect_allele_col='ALT.x',other_allele_col='REF.x',
                      eaf_col='maf',samplesize_col='ma_samples',pval_col='P_exp')})
        
        curr.exp = lapply(curr.exp,function(x) { 
          x %<>% dplyr::rename(rsid=SNP,pval=pval.exposure)
          x = ld_clump(dat=x,plink_bin='/home/richards/julian.willett/plink',
                       bfile='/scratch/richards/satoshi.yoshiji/database/1KG/EUR/EUR')
          x %<>% dplyr::rename(SNP=rsid,pval.exposure=pval)
          gc()
          x
        })
        gc()
        
        #outcome
        curr.out = lapply(split.by.gene,function(x) {
          format_data(x,type='outcome',chr_col='CHR',pos_col='POS.y',
                      beta_col='all_inv_var_meta_beta',se_col='all_inv_var_meta_sebeta',
                      snp_col='rsid',
                      effect_allele_col='ALT.y',other_allele_col='REF.y',
                      eaf_col='all_meta_AF',samplesize_col='Outcome_Total_N',
                      pval_col='all_inv_var_meta_p')})
        gc()
        exp.out.list = lapply(1:length(curr.exp),function(x) {
          list(curr.exp[[x]],curr.out[[x]])
        })
        names(exp.out.list) = names(split.by.gene)
        saveRDS(exp.out.list,glue('Fast_GTEx_eQTL/pre_harmonized/{f.name}_{curr.loc$Outcome}.rds'))
      }
      gc()
    }
  }else if (dataset == 'SC_Brain') {
    files = list.files('Bryois_Brain_SingleCell',
                       recursive=T,full.names=T,all.files=T,include.dirs=T,
                       pattern='gz')
    for (l in oi[[1]]:oi[[2]]) {
      print(glue('On locus {l} of {nrow(loci)}. Time: {Sys.time()}'))
      curr.loc = loci[l,]
      curr.files = files[which(str_detect(files,glue('\\.{curr.loc$Chr}\\.')))]
      file.names = gsub('.*\\/','',curr.files) %>% gsub('\\..*','',.)
      outcome.df = vroom(glue('Fast_HGI_Outcomes/{loci$Chr[[l]]}_{loci$Pos[[l]]}_{loci$Outcome[[l]]}.txt'),show_col_types=F)
      
      for (f in curr.files) { 
        print(glue('On file {f}. Locus {l}. Time: {Sys.time()}'))
        f.name = gsub('.*\\/','',f) %>% gsub('\\..*','',.)
        harmonize.saved.names = glue('{f.name}_{curr.loc$Chr}_{curr.loc$Pos}_{curr.loc$Outcome}')
        if (paste0(harmonize.saved.names,'.rds') %in% list.files('Fast_Brain_sceQTL/pre_harmonized')) next
        
        exp.out.df = vroom(f,show_col_types=F,col_names=c('GENE','rsid','Distance_TSS','P','BETA')) %>%
          dplyr::filter(abs(Distance_TSS) <= 500000,rsid %in% outcome.df$rsid)
        
        if (nrow(exp.out.df)<1) next
        exp.out.df %<>% merge(.,outcome.df,by='rsid') %>% drop_na(rsid) %>% 
          dplyr::mutate(SNP=rsid,P_exp=P,Exp_N=192,
                        Outcome_Total_N=(all_inv_var_meta_cases + all_inv_var_meta_controls))
        split.by.gene = split(exp.out.df,f=exp.out.df$GENE,drop=T) 
        split.by.gene = split.by.gene[sapply(split.by.gene, nrow)>=50]
        
        #exposure
        curr.exp = lapply(split.by.gene,function(x) {
          format_data(x,type='exposure',chr_col='#CHR',pos_col='POS',
                      beta_col='BETA',#se_col='slope_se',
                      snp_col='rsid',
                      effect_allele_col='ALT',other_allele_col='REF',
                      eaf_col='all_meta_AF',samplesize_col='Exp_N',pval_col='P_exp')})
        curr.exp = lapply(curr.exp,function(x) { clump_data(x,pop='EUR') })
        gc()
        
        #outcome
        curr.out = lapply(split.by.gene,function(x) {
          format_data(x,type='outcome',chr_col='#CHR',pos_col='POS',
                      beta_col='all_inv_var_meta_beta',se_col='all_inv_var_meta_sebeta',
                      snp_col='rsid',
                      effect_allele_col='ALT',other_allele_col='REF',
                      eaf_col='all_meta_AF',samplesize_col='Outcome_Total_N',
                      pval_col='all_inv_var_meta_p')})
        gc()
        
        exp.out.list = lapply(1:length(curr.exp),function(x) {
          list(curr.exp[[x]],curr.out[[x]])
        })
        names(exp.out.list) = names(split.by.gene)
        saveRDS(exp.out.list,glue('Fast_Brain_sceQTL/pre_harmonized/{harmonize.saved.names}.rds'))
      }
      gc()
    }
  }
}
MRFunction = function(h) {
  steig = directionality_test(h)
  curr.mr = data.frame()
  
  if (sum(h$mr_keep)==0) {
    print('No SNP left post-harmonization')
    curr.mr = data.frame(Label=NA)
  }else if (sum(h$mr_keep)==1) {
    m = mr(h, method_list = "mr_wald_ratio")
    curr.mr = data.frame(Label=NA,nSnps=m$nsnp,Beta=m$b,SE=m$se,
                         Pval=m$pval,SteigerP=steig$steiger_pval)
  }else if (sum(h$mr_keep)>1) {
    m = mr(h, method_list = c("mr_ivw","mr_two_sample_ml","mr_weighted_median","mr_penalised_weighted_median","mr_weighted_mode"))
    if (sum(h$mr_keep)>2) {
      egger = mr_egger_regression(h$beta.exposure[h$mr_keep],
                                  h$beta.outcome[h$mr_keep],
                                  h$se.exposure[h$mr_keep],
                                  h$se.outcome[h$mr_keep])
      curr.mr = data.frame(Label=NA,nSnps=m$nsnp[1],Beta=m$b[1],
                           SE=m$se[1],Pval=m$pval[1],WMedianBeta=m$b[3],
                           WModeBeta=m$b[5],EggerBeta=egger$b,
                           EggerInterceptP=egger$pval_i,SteigerP=steig$steiger_pval)
    } else {
      curr.mr = data.frame(Label=NA,nSnps=m$nsnp[1],Beta=m$b[1],
                           SE=m$se[1],Pval=m$pval[1],SteigerP=steig$steiger_pval)
    }
  }
  return(curr.mr)
}
doMR = function(dataset,oi) {
  loci = vroom('subpositions_chrpos.txt',col_names=F) %>%
    dplyr::arrange(X3) %>% dplyr::rename(Chr=X1,Pos=X2,Outcome=X3)
  
  if (dataset == 'BQC19') {
    for (l in oi[[1]]:oi[[2]]) {
      print(glue('On locus {l} of {nrow(loci)}. Time: {Sys.time()}'))
      
      curr.loc = loci[l,]
      
      curr.files = (list.files('Fast_BQC_eQTL/pre_harmonized',full.names=T))
      curr.files = curr.files[which(str_detect(curr.files,glue('chr{curr.loc$Chr}_{curr.loc$Pos}')))]
      
      for (f in curr.files) { 
        print(glue('Curr file: {f}. Locus {l}'))
        label = gsub('.*\\/','',f) %>% gsub('\\.rds','',.)
        if (glue('{label}.txt') %in% list.files('Fast_BQC_eQTL/mr')) next
        
        pre.harm.list = readRDS(f)
        site = (str_split(label,'_')[[1]]) ; site = site[c(1:(length(site)-1))]
        if (site[[length(site)]] == 'inf') curr.N = 112
        else curr.N = 166
        
        curr.mr = lapply(pre.harm.list,function(x) {
          if (nrow(x[[1]])>0) {
            x[[1]]$samplesize.exposure = curr.N #needing to correct N size
            harm = harmonise_data(x[[1]],x[[2]]) 
            mr.out = MRFunction(harm)
          }})
        curr.mr.collapsed = dplyr::bind_rows(curr.mr)
        curr.mr.collapsed$Label = glue('{label}_{names(curr.mr)[!sapply(curr.mr,is.null)]}')
        if (nrow(curr.mr.collapsed)>0) 
          vroom_write(curr.mr.collapsed %>% drop_na(Beta),glue('Fast_BQC_eQTL/mr/{label}.txt'))
      }
    }
  }else if (dataset == 'GTEx') {
    sample.database = vroom('/project/richards/restricted/dbGap/prj_32756/Documents/linking_file_GTEx.txt')
    for (l in oi[[1]]:oi[[2]]) {
      print(glue('On locus {l} of {nrow(loci)}. Time: {Sys.time()}'))
      
      curr.loc = loci[l,]
      
      curr.files = (list.files('Fast_GTEx_eQTL/pre_harmonized',full.names=T))
      curr.files = curr.files[which(str_detect(curr.files,glue('chr{curr.loc$Chr}_{curr.loc$Pos}')))]
      
      for (f in curr.files) { 
        print(glue('Curr file: {f}. Locus {l}'))
        label = gsub('.*\\/','',f) %>% gsub('\\.rds','',.)
        if (glue('{label}.txt') %in% list.files('Fast_GTEx_eQTL/mr')) next

        pre.harm.list = readRDS(f)
        site = (str_split(label,'_')[[1]]) ; site = site[c(1:(length(site)-3))]
        curr.N = sum(Reduce(`&`, lapply(site, grepl, x=sample.database$body_site)))
        if (curr.N == 0) stop('Unknown N for exposure')
        
        curr.mr = lapply(pre.harm.list,function(x) {
          if (nrow(x[[1]])>0) {
            x[[1]]$samplesize.exposure = curr.N #needing to correct N size
            harm = harmonise_data(x[[1]],x[[2]]) 
            mr.out = MRFunction(harm)
          }})
        curr.mr.collapsed = dplyr::bind_rows(curr.mr)
        curr.mr.collapsed$Label = glue('{label}_{names(curr.mr)[!sapply(curr.mr,is.null)]}')
        if (nrow(curr.mr.collapsed)>0) vroom_write(curr.mr.collapsed %>% drop_na(Beta),glue('Fast_GTEx_eQTL/mr/{label}.txt'))
      }
    }
  }else if (dataset == 'SC_Brain') {
    for (l in oi[[1]]:oi[[2]]) {
      print(glue('On locus {l} of {nrow(loci)}. Time: {Sys.time()}'))
      
      curr.loc = loci[l,]
      
      curr.files = (list.files('Fast_Brain_sceQTL/pre_harmonized',full.names=T))
      curr.files = curr.files[which(str_detect(curr.files,glue('{curr.loc$Chr}_{curr.loc$Pos}')))]
      
      for (f in curr.files) { 
        print(glue('Curr file: {f}. Locus {l}'))
        label = gsub('.*\\/','',f) %>% gsub('\\.rds','',.)
        if (glue('{label}.txt') %in% list.files('Fast_Brain_sceQTL/mr')) next
        
        pre.harm.list = readRDS(f)
        site = (str_split(label,'_')[[1]]) ; site = site[c(1:(length(site)-3))]
        
        curr.mr = lapply(pre.harm.list,function(x) {
          if (nrow(x[[1]])>0) {
            x[[1]]$se.exposure = abs(x[[1]]$beta.exposure/sqrt(qchisq(x[[1]]$pval.exposure,df=1,lower.tail=F)))
            harm = harmonise_data(x[[1]],x[[2]]) 
            mr.out = MRFunction(harm)
          }})
        curr.mr.collapsed = bind_rows(curr.mr)
        label.names = names(curr.mr[!sapply(curr.mr,is.null)])
        curr.mr.collapsed$Label = glue('{label}_{label.names}')
        if (nrow(curr.mr.collapsed)>0)
          vroom_write(curr.mr.collapsed,glue('Fast_Brain_sceQTL/mr/{label}.txt'))
        gc()
      }
    }
  }
}
readMRData = function(dataset) {
  if (dataset == 'BQC19') {
    # mr.files = list.files('Fast_BQC_eQTL/mr')
    # df = lapply(mr.files,function(x) {
    #   if (which(mr.files == x) %% 20 == 0) { print(which(mr.files == x)) ; gc()}
    #   vroom(glue('Fast_BQC_eQTL/mr/{x}'),show_col_types=F)
    # })
    # mr.df = dplyr::bind_rows(df)
    # mr.pass.sens = mr.df %>%
    #   dplyr::mutate(PvalAdj = p.adjust(Pval,method='BH'),PassSens=NA)
    # 
    # for (row in 1:nrow(mr.pass.sens)) {
    #   if (row %% 1000 == 0) print(row)
    #   curr.row = mr.pass.sens[row,]
    #   if (is.na(curr.row$WMedianBeta))
    #     mr.pass.sens$PassSens[[row]] = (curr.row$PvalAdj<=0.05 & curr.row$SteigerP <= 0.05)
    #   else
    #     mr.pass.sens$PassSens[[row]] = (curr.row$PvalAdj<=0.05 & curr.row$SteigerP <= 0.05 &
    #                                    sign(curr.row$Beta)==sign(curr.row$WMedianBeta) &
    #                                    sign(curr.row$Beta)==sign(curr.row$WModeBeta) &
    #                                    sign(curr.row$Beta)==sign(curr.row$EggerBeta) &
    #                                    curr.row$EggerInterceptP > 0.05)
    # }
    # saveRDS(mr.pass.sens,'Fast_BQC_eQTL/mr/collective_mr.results')
    return(readRDS('Fast_BQC_eQTL/mr/collective_mr.results'))
  }else if (dataset == 'GTEx') {
    # mr.files = list.files('Fast_GTEx_eQTL/mr')
    # df = lapply(mr.files,function(x) {
    #   if (which(mr.files == x) %% 20 == 0) { print(which(mr.files == x)) ; gc()}
    #   vroom(glue('Fast_GTEx_eQTL/mr/{x}'),show_col_types=F)
    # })
    # mr.df = dplyr::bind_rows(df)
    # mr.pass.sens = mr.df %>%
    #   dplyr::mutate(PvalAdj = p.adjust(Pval,method='BH'),PassSens=NA)
    # 
    # for (row in 1:nrow(mr.pass.sens)) {
    #   if (row %% 1000 == 0) print(row)
    #   curr.row = mr.pass.sens[row,]
    #   if (is.na(curr.row$WMedianBeta))
    #     mr.pass.sens$PassSens[[row]] = (curr.row$PvalAdj<=0.05 & curr.row$SteigerP <= 0.05)
    #   else
    #     mr.pass.sens$PassSens[[row]] = (curr.row$PvalAdj<=0.05 & curr.row$SteigerP <= 0.05 &
    #                                    sign(curr.row$Beta)==sign(curr.row$WMedianBeta) &
    #                                    sign(curr.row$Beta)==sign(curr.row$WModeBeta) &
    #                                    sign(curr.row$Beta)==sign(curr.row$EggerBeta) &
    #                                    curr.row$EggerInterceptP > 0.05)
    # }
    # saveRDS(mr.pass.sens,'Fast_GTEx_eQTL/mr/collective_mr.results')
    return(readRDS('Fast_GTEx_eQTL/mr/collective_mr.results'))
  }else if (dataset == 'SC_Brain') {
    # mr.files = list.files('Fast_Brain_sceQTL/mr')
    # df = lapply(mr.files,function(x) {
    #   if (which(mr.files == x) %% 20 == 0) { print(which(mr.files == x)) ; gc()}
    #   vroom(glue('Fast_Brain_sceQTL/mr/{x}'),show_col_types=F)
    # })
    # mr.df = dplyr::bind_rows(df)
    # mr.pass.sens = mr.df %>%
    #   dplyr::mutate(PvalAdj = p.adjust(Pval,method='BH'),PassSens=NA)
    # 
    # for (row in 1:nrow(mr.pass.sens)) {
    #   if (row %% 1000 == 0) print(row)
    #   curr.row = mr.pass.sens[row,]
    #   if (is.na(curr.row$WMedianBeta))
    #     mr.pass.sens$PassSens[[row]] = (curr.row$PvalAdj<=0.05 & curr.row$SteigerP <= 0.05)
    #   else
    #     mr.pass.sens$PassSens[[row]] = (curr.row$PvalAdj<=0.05 & curr.row$SteigerP <= 0.05 &
    #                                    sign(curr.row$Beta)==sign(curr.row$WMedianBeta) &
    #                                    sign(curr.row$Beta)==sign(curr.row$WModeBeta) &
    #                                    sign(curr.row$Beta)==sign(curr.row$EggerBeta) &
    #                                    curr.row$EggerInterceptP > 0.05)
    # }
    # saveRDS(mr.pass.sens,'Fast_Brain_sceQTL/mr/collective_mr.results')
    return(readRDS('Fast_Brain_sceQTL/mr/collective_mr.results'))
  }
}
getColocInputs = function(dataset,label) {
  if (dataset == 'BQC19') {
    split.label = str_split(label,'_')[[1]] %>% str_replace(.,'chr','')
    outcome.df = vroom(glue('Fast_HGI_Outcomes_500000/{paste0(split.label[c(1,2,5)],collapse="_")}.txt'),show_col_types=F) %>%
      dplyr::mutate(SNP = glue('chr{SNP}'))
    
    exp.out.df = vroom(glue('BQC_eQTL_Data/Subpositions/chr{paste0(split.label[1:(length(split.label)-2)],collapse="_")}.txt'),show_col_types=F) %>%
      dplyr::rename(SNP=variant_id) %>% 
      dplyr::filter(SNP %in% outcome.df$SNP,abs(tss_distance)<=500000,
                    gene_id == split.label[[length(split.label)]]) %>%
      merge(.,outcome.df,by='SNP') %>% drop_na(rsid) %>% 
      dplyr::mutate(SNP=rsid,P_exp=pval_nominal,
                    Outcome_Total_N=all_inv_var_meta_cases + all_inv_var_meta_controls,
                    Exp_N=NA) 
    if (split.label[[4]] == 'inf') exp.out.df$Exp_N = 112
    else if (split.label[[4]] == 'noninf') exp.out.df$Exp_N = 166
    
    # get sub df, and return
    eqtl.df = data.frame(pvalues=exp.out.df$pval_nominal,N=exp.out.df$Exp_N,
                         MAF=exp.out.df$maf,beta=exp.out.df$slope,varbeta=exp.out.df$slope_se^2,
                         type='quant',snp=exp.out.df$rsid,position=exp.out.df$POS.x)
    out.df = data.frame(pvalues=exp.out.df$all_inv_var_meta_p,N=exp.out.df$Outcome_Total_N,
                        MAF=exp.out.df$all_meta_AF,beta=exp.out.df$all_inv_var_meta_beta,
                        varbeta=exp.out.df$all_inv_var_meta_sebeta^2,type='cc',
                        s=(exp.out.df$all_inv_var_meta_cases/exp.out.df$Outcome_Total_N),
                        snp=exp.out.df$rsid,position=exp.out.df$POS.y) %>% drop_na(MAF)
    return(list(eqtl.df,out.df))
  }else if (dataset == 'GTEx') {
    # get master df
    split.label = str_split(label,'_')[[1]] %>% str_replace(.,'chr','')
    outcome.df = vroom(glue('Fast_HGI_Outcomes/{paste0(split.label[(length(split.label)-3):(length(split.label)-1)],collapse="_")}.txt'),show_col_types=F) %>%
      dplyr::mutate(SNP=str_replace_all(SNP,':','_')) %>%
      dplyr::mutate(SNP=glue('chr{SNP}_b38'))
    
    split.label = str_split(label,'_')[[1]]
    exp.out.df = vroom(glue('Fast_GTEx_eQTL/{paste0(split.label[1:(length(split.label)-2)],collapse="_")}.txt'),show_col_types=F) %>%
      dplyr::rename(SNP=variant_id) %>% 
      dplyr::filter(SNP %in% outcome.df$SNP,abs(tss_distance)<=500000) %>%
      merge(.,outcome.df,by='SNP') %>% drop_na(rsid) %>% 
      dplyr::mutate(SNP=rsid,P_exp=pchisq((slope/slope_se)^2,df=1,lower.tail=F),
                    Outcome_Total_N=all_inv_var_meta_cases + all_inv_var_meta_controls,
                    Exp_N=round(ma_count/maf*(1/2))) %>%
      dplyr::filter(phenotype_id == split.label[[length(split.label)]])
    
    # get sub df, and return
    eqtl.df = data.frame(pvalues=exp.out.df$pval_nominal,N=exp.out.df$Exp_N,
                         MAF=exp.out.df$maf,beta=exp.out.df$slope,varbeta=exp.out.df$slope_se^2,
                         type='quant',snp=exp.out.df$rsid,position=exp.out.df$POS.x)
    out.df = data.frame(pvalues=exp.out.df$all_inv_var_meta_p,N=exp.out.df$Outcome_Total_N,
                        MAF=exp.out.df$all_meta_AF,beta=exp.out.df$all_inv_var_meta_beta,
                        varbeta=exp.out.df$all_inv_var_meta_sebeta^2,type='cc',
                        s=(exp.out.df$all_inv_var_meta_cases/exp.out.df$Outcome_Total_N),
                        snp=exp.out.df$rsid,position=exp.out.df$POS.y) %>% drop_na(MAF)
    return(list(eqtl.df,out.df))
  }else if (dataset == 'SC_Brain') {
    split.label = str_split(label,'_')[[1]]
    outcome.df = vroom(glue('Fast_HGI_Outcomes/{paste0(split.label[(length(split.label)-4):(length(split.label)-2)],collapse="_")}.txt'),show_col_types=F) %>%
      dplyr::mutate(SNP=str_replace_all(SNP,':','_')) %>%
      dplyr::mutate(SNP=glue('chr{SNP}_b38'))
    exp.out.df = vroom(glue('Bryois_Brain_SingleCell/{split.label[[1]]}.{split.label[[2]]}.gz'),
                       show_col_types = F,col_names=c('GENE','rsid','Distance_TSS','P','BETA')) %>%
      dplyr::filter(abs(Distance_TSS) <= 500000,rsid %in% outcome.df$rsid,
                    GENE==paste(split.label[(length(split.label)-1):length(split.label)],collapse='_')) %>%
      merge(.,outcome.df,by='rsid') %>% drop_na(rsid) %>% 
      dplyr::mutate(P_exp=P,Exp_N=192,Outcome_Total_N=(all_inv_var_meta_cases + all_inv_var_meta_controls)) %>%
      dplyr::mutate(SE_exp=abs(BETA/sqrt(qchisq(P_exp,df=1,lower.tail=F))))
    
    # get sub df, and return
    eqtl.df = data.frame(pvalues=exp.out.df$P_exp,N=exp.out.df$Exp_N,
                         MAF=exp.out.df$all_meta_AF,beta=exp.out.df$BETA,varbeta=exp.out.df$SE_exp^2,
                         type='quant',snp=exp.out.df$rsid,position=exp.out.df$POS)
    out.df = data.frame(pvalues=exp.out.df$all_inv_var_meta_p,N=exp.out.df$Outcome_Total_N,
                        MAF=exp.out.df$all_meta_AF,beta=exp.out.df$all_inv_var_meta_beta,
                        varbeta=exp.out.df$all_inv_var_meta_sebeta^2,type='cc',
                        s=(exp.out.df$all_inv_var_meta_cases/exp.out.df$Outcome_Total_N),
                        snp=exp.out.df$rsid,position=exp.out.df$POS) %>% drop_na(MAF)
    return(list(eqtl.df,out.df))
  }
}
doColoc = function(dataset,label.index) { #label index for getting data for sensitivity testing
  if (dataset == 'BQC19') {
    if (!file.exists('Fast_BQC_Coloc_Results.rds') | !is.na(label.index)) {
      mr.data = readMRData(dataset)
      sig.mr = mr.data %>% dplyr::filter(PassSens)
      if (!is.na(label.index)) sig.mr = sig.mr[label.index,]

      results = lapply(sig.mr$Label,function (x) { #get coloc for each sig element
        print(glue('On {x} {which(sig.mr$Label == x)} of {nrow(sig.mr)}'))
        gc()
        ins = getColocInputs(dataset,x)
        coloc.res = coloc.abf(ins[[1]] %>% drop_na(pvalues) %>% dplyr::filter(N>0,pvalues<1),
                              ins[[2]] %>% drop_na(pvalues))$summary
        coloc.df = data.frame(Label=x,H0.PP=coloc.res['PP.H0.abf'],H1.PP=coloc.res['PP.H1.abf'],
                              H2.PP=coloc.res['PP.H2.abf'],H3.PP=coloc.res['PP.H3.abf'],
                              H4.PP=coloc.res['PP.H4.abf'])
        if (coloc.res['PP.H4.abf'] >= 0.8) list(ins,coloc.df)
        else list(NA,coloc.df)
      })
      saveRDS(results,'Fast_BQC_Coloc_Results.rds')
    }else{
      results = readRDS('Fast_BQC_Coloc_Results.rds')
    }
    non.null = lapply(results,function(x) {
      if (length(x[[1]])>1) x
    })
    non.null = Filter(Negate(is.null),non.null)
  }else if (dataset == 'GTEx') { 
    if (!file.exists('Fast_GTEx_Coloc_Results.rds') | !is.na(label.index)) {
      mr.data = readMRData(dataset)
      sig.mr = mr.data %>% dplyr::filter(PassSens)
      if (!is.na(label.index)) sig.mr = sig.mr[label.index,]

      coloc.out = lapply(sig.mr$Label,function (x) { #get coloc for each sig element
        print(glue('On {x} {which(sig.mr$Label == x)} of {nrow(sig.mr)}'))
        gc()
        ins = getColocInputs(dataset,x)
        coloc.res = coloc.abf(ins[[1]] %>% drop_na(pvalues) %>% dplyr::filter(N>0),
                              ins[[2]] %>% drop_na(pvalues))$summary
        coloc.df = data.frame(Label=x,H0.PP=coloc.res['PP.H0.abf'],H1.PP=coloc.res['PP.H1.abf'],
                        H2.PP=coloc.res['PP.H2.abf'],H3.PP=coloc.res['PP.H3.abf'],
                        H4.PP=coloc.res['PP.H4.abf'])
        if (coloc.res['PP.H4.abf'] >= 0.8) list(ins,coloc.df)
        else list(NA,coloc.df)
      })
      saveRDS(coloc.out,'Fast_GTEx_Coloc_Results.rds')
    }else{
      results = readRDS('Fast_GTEx_Coloc_Results.rds')
    }
    non.null = lapply(results,function(x) {
      if (length(x[[1]])>1) x
    })
    non.null = Filter(Negate(is.null),non.null)
  }
  return(non.null)
}
makeColocSensitivityPlots = function(df,name) {
  plts = list()
  for (plt in 1:length(df)) {
    plt.name = df[[plt]][[2]]$Label
    sens.one = ggplot(df[[plt]][[1]][[1]],aes(x=position,y=-log10(pvalues))) +
      geom_point() 
    sens.two = ggplot(df[[plt]][[1]][[2]],aes(x=position,y=-log10(pvalues))) +
      geom_point() 
    plts[[plt]] = list(sens.one,sens.two,plt.name)
  }
  
  pdf(glue('{name}_sens_plots.pdf'))
  for (plt in plts) {
    print(plt[[3]])
    print(grid.arrange(plt[[1]],plt[[2]],ncol=2,top=plt[[3]],heights=c(2,2),widths=c(1,1)))
  }
  dev.off()
}
getBroaderWindow = function(dataset,df,label) {
  split.label = str_split(label,'_')[[1]]
  
  if (dataset == 'GTEx') {
    chr = str_replace(split.label[[length(split.label)-3]],'chr','')
    pos = as.numeric(split.label[[length(split.label)-2]])
    out = split.label[[length(split.label)-1]]
    file = list.files('/scratch/richards/satoshi.yoshiji/database/GTEx_v8/eQTL',
                      recursive=T,full.names=T,all.files=T,include.dirs=T,
                      pattern='tsv') %>%
      grep(pattern=paste(split.label[1],collapse='_'),value=T) %>%
      grep(pattern=glue('\\.{chr}\\.eQTL'),value=T)
    print(paste('File:',file))
    
    outcome.df = vroom(glue('Fast_HGI_Outcomes_1e+06/{chr}_{pos}_{out}.txt'),show_col_types=F)
    eqtl = vroom(file) %>% 
      dplyr::filter(tss_distance <= 500000,POS>=pos-1000000,POS<=pos+1000000,
                    phenotype_id==split.label[[length(split.label)]])
    plt.one = ggplot(eqtl,aes(x=POS,y=-log10(pval_nominal))) + geom_point() + 
            geom_vline(xintercept=c(pos,pos-1000000,pos+1000000))
    plt.two = ggplot(outcome.df,aes(x=POS,y=-log10(all_inv_var_meta_p))) + geom_point() + 
      geom_vline(xintercept=c(pos,pos-1000000,pos+1000000))
    grid.arrange(plt.one,plt.two)
  }else if (dataset == 'SC_Brain') {
    chr = split.label[[length(split.label)-4]]
    pos = split.label[[length(split.label)-3]]
    out = split.label[[length(split.label)-2]]
    gene = paste(split.label[(length(split.label)-1):length(split.label)],collapse='_')
    file = list.files('Bryois_Brain_SingleCell',
                      recursive=T,full.names=T,all.files=T,include.dirs=T,
                      pattern='gz') %>%
      grep(pattern=split.label[[1]],value=T) %>%
      grep(pattern=glue('\\.{chr}\\.gz'),value=T)
    
    #read in required files
    outcome.df = vroom(glue('Fast_HGI_Outcomes_1e+06/{chr}_{pos}_{out}.txt'),show_col_types=F)
    pos = as.numeric(pos)
    print(file)
    eqtl = vroom(file,show_col_types=F,col_names=c('GENE','rsid','Distance_TSS','P','BETA')) %>% 
      dplyr::filter(Distance_TSS <= 500000,rsid %in% outcome.df$rsid,
                    GENE==gene) %>% 
      merge(.,outcome.df,by='rsid') %>% dplyr::filter(POS>=pos-1000000,POS<=pos+1000000)
    plt.one = ggplot(eqtl,aes(x=POS,y=-log10(P))) + geom_point() + 
            geom_vline(xintercept=c(pos,pos-1000000,pos+1000000))
    plt.two = ggplot(eqtl,aes(x=POS,y=-log10(all_inv_var_meta_p))) + geom_point() + 
      geom_vline(xintercept=c(pos,pos-1000000,pos+1000000))
    grid.arrange(plt.one,plt.two)
  }
}
getFinalHits = function(dataset,df,to.cut) {
  if (is.na(to.cut))
    subset.list = df
  else
    subset.list = df[-to.cut]
  
  out.df = lapply(subset.list,function(x) {
    spl = str_split(x[[2]]$Label,'_')[[1]]
    len = length(spl)
    if (dataset == 'GTEx') {
      tissue = paste(spl[1:(len-4)],collapse='_')
      outcome = spl[[len-1]]
      gene = spl[[len]] 
      chr = spl[[len-3]] %>% str_replace('chr','')
      pos = spl[[len-2]]
      data.frame(Tissue=tissue,Gene=gene,Outcome=outcome,H4=x[[2]]$H4.PP,Chr=chr,Pos=pos)
    }else if (dataset == 'SC_Brain') {
      tissue = paste(spl[1:(len-5)],collapse='_')
      outcome = spl[[len-2]]
      gene = paste(spl[(len-1):len])
      chr = spl[[len-4]]
      pos = spl[[len-3]]
      data.frame(Tissue=tissue,Gene=gene,Outcome=outcome,H4=x[[2]]$H4.PP,Chr=chr,Pos=pos)
    }else if (dataset == 'BQC') {
      tissue = paste(spl[3:4],collapse='_')
      outcome = spl[[5]]
      gene = spl[[6]]
      chr = as.numeric(str_replace(spl[[1]],'chr',''))
      pos = spl[[2]]
      data.frame(Tissue=tissue,Gene=gene,Outcome=outcome,H4=x[[2]]$H4.PP,Chr=chr,Pos=pos)
    }
  })
  out.df = dplyr::bind_rows(out.df) %>% 
    #dplyr::mutate(Gene = replaceGeneNames(Gene),Gene = sub('\\..*','',Gene)) %>%
    dplyr::filter(Gene %notin% c('ENSG00000204525','ENSG00000204622')) %>%
    dplyr::filter(!(Chr == 6 & Pos >= 28510120 & Pos <= 33480577))
  
  return(out.df)
}
plotNewMR = function(dataset,soskic,brain.sc,bqc,gtex,fs) { #fs = font size
  if (!file.exists('mr_plot_final_sc_plot.rds') & dataset == 'SC') {
    # Note that brain hits duplicated entries, putting simple and long gene names on different lines vs same line
    # mr.plt.df = brain.sc %>% dplyr::mutate(B=NA,SE=NA,P=NA) %>%
    #   dplyr::mutate(Tissue=ifelse(Tissue == 'Excitatory.neurons','Excitatory',Tissue)) %>%
    #   dplyr::slice(which(dplyr::row_number() %% 2 == 1)) #select long gene name row
    # mr.data = readMRData('SC_Brain')
    # 
    # for (r in 1:nrow(mr.plt.df)) {
    #   print(glue('Row: {r} of {nrow(mr.plt.df)}'))
    #   curr.row = mr.plt.df[r,]
    #   match = mr.data %>% dplyr::filter(str_detect(Label,curr.row$Tissue),
    #                                     str_detect(Label,curr.row$Outcome),
    #                                     str_detect(Label,curr.row$Gene),
    #                                     str_detect(Label,curr.row$Pos))
    #   if (nrow(match)>1) stop(match)
    #   mr.plt.df[r,'B'] = match$Beta
    #   mr.plt.df[r,'SE'] = match$SE
    #   mr.plt.df[r,'P'] = match$PvalAdj
    # }
    # mr.plt.df %<>% dplyr::mutate(Tissue = str_replace_all(Tissue,'pb','Pan Brain'),
    #                              Tissue = factor(Tissue,levels=c('Astrocytes','Excitatory',
    #                                                 'Microglia','Pan Brain')))
    
    mr.plt.df = soskic %>% dplyr::mutate(Tissue=CellState,B=Beta)
    # mr.plt.df %<>% dplyr::bind_rows(soskic.df)
    mr.plt.df %<>% dplyr::mutate(PosOR=as.numeric(B>=0),PosOR=ifelse(PosOR==0,-1,PosOR)) %>%
      dplyr::mutate(Gene = replaceGeneNames(Gene),Gene = sub('\\..*','',Gene))
      
    
    final.plt = mr.plt.df
    for (o in c('A2','B2','C2')) { #add missing entries
      tmp = mr.plt.df %>% dplyr::filter(Outcome == o)
      for (g in tmp$Gene) { #add missing entries
        curr.set = tmp %>% dplyr::filter(Gene==g)
        missing.organs = unique((mr.plt.df %>% dplyr::filter(Tissue %notin% curr.set$Tissue))$Tissue)
        final.plt %<>% dplyr::add_row(Tissue=missing.organs,Gene=g,PosOR=0,Outcome=o)
      }
    }
    final.plt %<>% dplyr::mutate(Organ = NA) %>%
      dplyr::mutate(Organ = ifelse(Tissue %in% soskic$CellState,'Whole Blood','Brain')) %>%
      dplyr::mutate(Tissue = str_replace_all(Tissue,'_',' ')) 
      
    saveRDS(final.plt,'mr_plot_final_sc_plot.rds')
  }else if (!file.exists('mr_plot_final_bulk_plot.rds') & dataset == 'Bulk') {
    mr.plt.df = gtex %>% dplyr::mutate(B=NA,SE=NA)
    bqc.plt.df = bqc %>% dplyr::mutate(B=NA,SE=NA)
    
    gtex.mr.data = readMRData('GTEx')
    bqc.mr.data = readMRData('BQC19')

    for (r in 1:nrow(mr.plt.df)) { #organize gtex data
      print(glue('Row: {r} of {nrow(mr.plt.df)}'))
      curr.row = mr.plt.df[r,]
      match = gtex.mr.data %>% dplyr::filter(str_detect(Label,curr.row$Tissue),
                                        str_detect(Label,curr.row$Outcome),
                                        str_detect(Label,curr.row$Gene),
                                        str_detect(Label,curr.row$Pos))
      if (nrow(match)>1) stop(match)
      mr.plt.df[r,'B'] = match$Beta
      mr.plt.df[r,'SE'] = match$SE
    }
    mr.plt.df %<>% dplyr::mutate(Tissue = ifelse(Tissue=='Whole_Blood',glue('{Tissue} COVID- Symptoms-'),Tissue))

    # next get data for bqc results
    for (r in 1:nrow(bqc.plt.df)) {
      print(glue('Row: {r} of {nrow(bqc.plt.df)}'))
      curr.row = bqc.plt.df[r,]
      match = bqc.mr.data %>% dplyr::filter(str_detect(Label,curr.row$Tissue),
                                        str_detect(Label,curr.row$Outcome),
                                        str_detect(Label,curr.row$Gene),
                                        str_detect(Label,curr.row$Pos))
      if (nrow(match)>1) stop(match)
      bqc.plt.df[r,'B'] = match$Beta
      bqc.plt.df[r,'SE'] = match$SE
    }
    #edit bqc df to go together with gtex (all bulk on one figure)
    bqc.plt.df %<>% dplyr::mutate(Label=Tissue) %>%
      dplyr::mutate(Tissue=ifelse(str_detect(Label,'_inf'),glue('Whole Blood COVID- Symptoms+'),Tissue),
                    Tissue=ifelse(str_detect(Label,'_noninf'),glue('Whole Blood COVID+ Symptoms+'),Tissue))

    mr.plt.df %<>% dplyr::bind_rows(bqc.plt.df) #merge data
    mr.plt.df %<>% dplyr::mutate(PosOR=as.numeric(B>=0),PosOR=ifelse(PosOR==0,-1,PosOR)) %>%
      dplyr::mutate(Gene = replaceGeneNames(Gene),Gene = sub('\\..*','',Gene)) %>%
      dplyr::mutate(Tissue = str_replace_all(Tissue,'_',' '))

    final.plt = mr.plt.df
    for (o in c('A2','B2','C2')) { #add missing entries
      tmp = mr.plt.df %>% dplyr::filter(Outcome == o)
      for (g in tmp$Gene) { #add missing entries
        curr.set = tmp %>% dplyr::filter(Gene==g)
        missing.organs = unique((mr.plt.df %>% dplyr::filter(Tissue %notin% curr.set$Tissue))$Tissue)
        final.plt %<>% dplyr::add_row(Tissue=missing.organs,Gene=g,PosOR=0,Outcome=o)
      }
    }
    saveRDS(final.plt,'mr_plot_final_bulk_plot.rds')
  }
  
  if (dataset == 'Bulk') final.plt = readRDS('mr_plot_final_bulk_plot.rds')
  else if (dataset == 'SC') final.plt = readRDS('mr_plot_final_sc_plot.rds')
  final.plt %<>% dplyr::mutate(Outcome = str_replace_all(Outcome,'A2','Severe'),
                               Outcome = str_replace_all(Outcome,'B2','Hospitalized'),
                               Outcome = str_replace_all(Outcome,'C2','Susceptibility'),
                               Outcome = factor(Outcome, levels = c('Severe','Hospitalized','Susceptibility'))) %>%
    dplyr::filter(Gene != 'ENSG00000204525')
  
  final.plt %<>% dplyr::filter(PosOR != 0)
  plt = ggplot(final.plt,aes(x=Gene,y=Tissue,fill=PosOR)) + geom_tile(color='black') + theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),text=element_text(size=fs),
          legend.position='none') +
    ylab('') + scale_y_discrete(limits = rev) + xlab('') +
    scale_fill_gradient2(low = "#075AFF",mid = "#FFFFCC",high = "#FF0000")
  
  if (dataset == 'SC') plt = plt + facet_grid(Organ ~ Outcome,scales='free')
  else plt = plt + facet_grid(. ~ Outcome,scales='free')
  
  # print(plt)
  print('Red = OR > 1. Blue = OR < 1. White = did not colocalize')
  
  ggsave(glue('mr_plot_{dataset}.png'),plt,width=15,height=7,dpi=1200)
  
  return(final.plt)
}
plotNewHeatmap = function(dataset,df) {
  for (out in c('A2','B2','C2')) {
    plt.df = df %>% dplyr::filter(Outcome == out) %>%
      dplyr::mutate(Tissue=as.factor(Tissue),Coloc=1)
    if (nrow(plt.df)==0) next
    
    final.df = plt.df
    for (g in plt.df$Gene) { #add missing entries
      curr.set = plt.df %>% dplyr::filter(Gene==g) 
      missing.organs = unique((df %>% dplyr::filter(Tissue %notin% curr.set$Tissue))$Tissue)
      final.df %<>% dplyr::add_row(Tissue=missing.organs,Gene=g,Outcome=out,
                            Coloc=0)
    }

    plt = ggplot(final.df,aes(x=Gene,y=Tissue,fill=Coloc)) +
      geom_tile(aes(linetype=(H4>=0.8)),show.legend=F) + theme_bw() +
      ylab('') + scale_fill_continuous(limits=c(0,1)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
            text=element_text(size=14)) +
      ggtitle(out) + scale_y_discrete(limits = rev(levels(plt.df$Tissue)))
    print(plt)
  }
}
makeDataLocusZoom = function(gene) {
  #GRCh38 build
  if (gene == 'IFNAR2') {
    out.list = list(vroom('BQC_eQTL_Data/Subpositions/chr21_33237639_window_inf.txt'),
                    vroom('BQC_eQTL_Data/Subpositions/chr21_33237639_window_noninf.txt'),
                    vroom('GTEx_eQTL_Data/Subpositions/chr21_33237639_window.txt'))
    names(out.list) = c('COVID+_Sx+','COVID-_Sx+','COVID-_Sx-')
    out.list[1:2] = lapply(out.list[1:2],function(x) {
      x %>% dplyr::filter(str_detect(gene_id,'ENSG00000159110'),
                          tss_distance < 500000,maf >= 0.01) %>%
        dplyr::select(CHR,POS,REF,ALT,pval_nominal) %>%
        dplyr::arrange(CHR,POS)
    })
    out.list[[3]] %<>% dplyr::filter(str_detect(phenotype_id,'ENSG00000159110'),
                                     tss_distance < 500000,maf >= 0.01) %>%
      dplyr::select(CHR,POS,REF,ALT,pval_nominal) %>%
      dplyr::arrange(CHR,POS)
    outcome.out.a = vroom('Fast_HGI_Outcomes_500000/21_33237639_A2.txt') %>%
      dplyr::filter(all_meta_AF >= 0.01) %>%
      dplyr::select(`#CHR`,POS,REF,ALT,all_inv_var_meta_p)
    outcome.out.b = vroom('Fast_HGI_Outcomes_500000/21_33237639_B2.txt') %>%
      dplyr::filter(all_meta_AF >= 0.01) %>%
      dplyr::select(`#CHR`,POS,REF,ALT,all_inv_var_meta_p)
    
    vroom_write(out.list[[1]],glue('LocusZoom_Data/{gene}_{names(out.list)[[1]]}.txt'))
    vroom_write(out.list[[2]],glue('LocusZoom_Data/{gene}_{names(out.list)[[2]]}.txt'))
    vroom_write(out.list[[3]],glue('LocusZoom_Data/{gene}_{names(out.list)[[3]]}.txt'))
    vroom_write(outcome.out.a,glue('LocusZoom_Data/{gene}_Outcome_A2.txt'))
    vroom_write(outcome.out.b,glue('LocusZoom_Data/{gene}_Outcome_B2.txt'))
  }else if (gene == 'RALGDS') { #A2
    out.list = readRDS('MRthenColoc_MRData/19_A2_9_133271182.rds')
    names(out.list[[1]]) = out.list[[2]]
    out.list = out.list[[1]][which(str_detect(out.list[[2]],'ENSG00000160271'))] #RALGDS
    out.list = out.list[c(17,18,10)]
    names(out.list) = c('TEM_0h','TEM_16h','TN_16h')
    
    outcome.df = out.list[[1]][[2]] %>% dplyr::select(CHR,position,REF,ALT,pvalues,rsid)
    out.list = lapply(out.list,function(x) {
      x[[1]] %>% dplyr::select(chr,position,REF,ALT,pvalues)
    })
    outcome.out.b = vroom('Fast_HGI_Outcomes_500000/9_133263862_B2.txt') %>%
      dplyr::filter(all_meta_AF >= 0.01) %>%
      dplyr::select(`#CHR`,POS,REF,ALT,all_inv_var_meta_p)
    outcome.out.c = vroom('Fast_HGI_Outcomes_500000/9_133271745_C2.txt') %>%
      dplyr::filter(all_meta_AF >= 0.01) %>%
      dplyr::select(`#CHR`,POS,REF,ALT,all_inv_var_meta_p)
    
    vroom_write(out.list[[1]],glue('LocusZoom_Data/{gene}_{names(out.list)[[1]]}.txt'))
    vroom_write(out.list[[2]],glue('LocusZoom_Data/{gene}_{names(out.list)[[2]]}.txt'))
    vroom_write(out.list[[3]],glue('LocusZoom_Data/{gene}_{names(out.list)[[3]]}.txt'))
    vroom_write(outcome.df,glue('LocusZoom_Data/{gene}_Outcome_A2.txt'))
    vroom_write(outcome.out.b,glue('LocusZoom_Data/{gene}_Outcome_B2.txt'))
    vroom_write(outcome.out.c,glue('LocusZoom_Data/{gene}_Outcome_C2.txt'))
  }else if (gene == 'IL10RB') { #A2
    outcome.df = vroom('Fast_HGI_Outcomes_500000/21_33237639_A2.txt',show_col_types=F)
    out.list = list(vroom('Bryois_Brain_SingleCell/Excitatory.neurons.21.gz',show_col_types = F,col_names=c('GENE','rsid','Distance_TSS','P','BETA')),
                    vroom('Bryois_Brain_SingleCell/Microglia.21.gz',show_col_types = F,col_names=c('GENE','rsid','Distance_TSS','P','BETA')),
                    vroom('Bryois_Brain_SingleCell/pb.21.gz',show_col_types = F,col_names=c('GENE','rsid','Distance_TSS','P','BETA')))
    names(out.list) = c('Excitatory','Microglia','BulkBrain')
    out.list = lapply(out.list,function(x) {
      x %>% dplyr::filter(GENE=='IL10RB_ENSG00000243646',Distance_TSS<=500000) %>% 
        merge(.,outcome.df,by='rsid') %>% dplyr::filter(all_meta_AF >= 0.01) %>%
        dplyr::select(`#CHR`,POS,REF,ALT,P) %>%
        dplyr::arrange(`#CHR`,POS)
    })
    brain.cortex = vroom('Fast_GTEx_eQTL/Brain_Cortex_chr21_33237639.txt')
    brain.cortex %<>% dplyr::filter(str_detect(phenotype_id,'ENSG00000243646'),
                                    tss_distance < 500000,maf >= 0.01) %>%
      dplyr::select(CHR,POS,REF,ALT,pval_nominal) %>%
      dplyr::arrange(CHR,POS)
    
    outcome.out.a = outcome.df %>%
      dplyr::filter(all_meta_AF >= 0.01) %>%
      dplyr::select(`#CHR`,POS,REF,ALT,all_inv_var_meta_p)
    
    vroom_write(out.list[[1]],glue('LocusZoom_Data/{gene}_{names(out.list)[[1]]}.txt'))
    vroom_write(out.list[[2]],glue('LocusZoom_Data/{gene}_{names(out.list)[[2]]}.txt'))
    vroom_write(out.list[[3]],glue('LocusZoom_Data/{gene}_{names(out.list)[[3]]}.txt'))
    vroom_write(brain.cortex,glue('LocusZoom_Data/{gene}_brain_cortex.txt'))
    vroom_write(outcome.out.a,glue('LocusZoom_Data/{gene}_Outcome_A2.txt'))
  }else if (gene == 'MUC5B') {
    out.list = list(vroom('Fast_GTEx_eQTL/Lung_chr11_1219991.txt'),
                    vroom('Fast_GTEx_eQTL/Brain_Cortex_chr11_1219991.txt'),
                    vroom('Fast_GTEx_eQTL/Whole_Blood_chr11_1219991.txt'))
    names(out.list) = c('Lung','Brain_Cortex','WholeBlood')
    out.list = lapply(out.list,function(x) {
      x %>% dplyr::filter(str_detect(phenotype_id,'ENSG00000117983'),
                          tss_distance < 500000,maf >= 0.01) %>%
        dplyr::select(CHR,POS,REF,ALT,pval_nominal) %>%
        dplyr::arrange(CHR,POS)
    })
    outcome.out.a = vroom('Fast_HGI_Outcomes_500000/11_1219991_A2.txt') %>%
      dplyr::filter(all_meta_AF >= 0.01) %>%
      dplyr::select(`#CHR`,POS,REF,ALT,all_inv_var_meta_p)
    
    vroom_write(out.list[[1]],glue('LocusZoom_Data/{gene}_{names(out.list)[[1]]}.txt'))
    vroom_write(out.list[[2]],glue('LocusZoom_Data/{gene}_{names(out.list)[[2]]}.txt'))
    vroom_write(out.list[[3]],glue('LocusZoom_Data/{gene}_{names(out.list)[[3]]}.txt'))
    vroom_write(outcome.out.a,glue('LocusZoom_Data/{gene}_Outcome_A2.txt'))
  }
}
makePercentColocFig = function() {
  df = data.frame(Dataset=c('CD4+ Whole Blood','BQC19 Whole Blood','GTEx Whole Blood','GTEx All'),
                  PercentColoc=c(0.011,0.050,0.027,0.031),
                  SequencingModality=c('SC','Bulk','Bulk','Bulk'))
  df$Dataset = factor(df$Dataset,levels=c('CD4+ Whole Blood','BQC19 Whole Blood','GTEx Whole Blood','GTEx All'))
  plt = ggplot(df,aes(x=Dataset,y=PercentColoc,fill=SequencingModality)) +
    geom_bar(stat='identity') + theme_bw() + xlab('') +
    theme(text=element_text(size=22),axis.text.x=element_text(angle=90),
          legend.position='none') + scale_fill_manual(values=cbPalette[1:2]) +
    ylab('% Colocalizing')
  print(plt)
  ggsave('percent_coloc.png',plt,width=10,height=6,dpi=1200)
}
makePercentViolateSCVFig = function() {
  df = data.frame(Dataset=c('CD4+ Whole Blood','BQC19 Whole Blood','GTEx Whole Blood','GTEx All'),
                  PercentViol=c(0.316,0.667,0.500,0.701),
                  SequencingModality=c('SC','Bulk','Bulk','Bulk'))
  df$Dataset = factor(df$Dataset,levels=c('CD4+ Whole Blood','BQC19 Whole Blood','GTEx Whole Blood','GTEx All'))
  plt = ggplot(df,aes(x=Dataset,y=PercentViol,fill=SequencingModality)) +
    geom_bar(stat='identity') + theme_bw() + xlab('') +
    theme(text=element_text(size=22),axis.text.x=element_text(angle=90),
          legend.position='none') + scale_fill_manual(values=cbPalette[1:2]) +
    ylab('% Violating Assumption')
  print(plt)
  ggsave('percent_violate.png',plt,width=10,height=6,dpi=1200)
}
makeSupplementalTables = function() {
  leadVariants = vroom('subpositions_chrpos.txt',col_names=F) %>%
    mutate(X3 = ifelse(X3 == 'A2','Severe',X3),X3 = ifelse(X3 == 'B2','Hospitalized',X3),
           X3 = ifelse(X3 == 'C2','Susceptibility',X3))
  names(leadVariants) = c('CHR','POS','Outcome')
  leadVariants = leadVariants[-c(26:31),]
  
  soskic.mr = readRDS('MRthenColoc_MRData/mr_results_all.rds')
  bqc.mr = readMRData('BQC19')
  gtex.mr = readMRData('GTEx')
  
  gtex.coloc = doColoc('GTEx',NA)
  bqc.coloc = doColoc('BQC19',NA) 
  
  soskic.coloc = readRDS('MRthenColoc_MRData/coloc_sens_all.rds')[[3]] 
  soskic.coloc$PassedSensitivity = ifelse(row.names(soskic.coloc) %notin% c(1,4,5,6,13,15),TRUE,FALSE)
  bqc.hits = getFinalHits('BQC',bqc.coloc,to.cut=NA)
  bqc.hits$PassedSensitivity = ifelse(row.names(bqc.hits) %notin% c(1:4,7:8,10),TRUE,FALSE)
  gtex.hits = getFinalHits('GTEx',gtex.coloc,to.cut=NA)
  gtex.hits$PassedSensitivity = ifelse(row.names(gtex.hits) %notin% c(2:5,8:15,20:21,23,25,27:31,33:37,
                                                                      40:42,45:61,64:117,121,126:129,
                                                                      131,137:142,145,148,150:162,164,170,
                                                                      174:181,183,190:196,198,200,205:211,
                                                                      215:217,220:226,231:232,234:243,246,
                                                                      248,251:253,255:256,260:265,275,277:286,
                                                                      288:308,312:316,319,323:326,331:333,
                                                                      336:338,340:342,344:345,350:351,353:354,
                                                                      356:357,360:374,376:377),
                                       TRUE,FALSE)
  
  wb = createWorkbook()
  addWorksheet(wb,'S1. Lead variants')
  addWorksheet(wb,'S3. MR CD4 Cells')
  addWorksheet(wb,'S4. MR BQC19')
  addWorksheet(wb,'S5. MR GTEx')
  addWorksheet(wb,'S6. Coloc CD4 Cells')
  addWorksheet(wb,'S7. Coloc BQC19')
  addWorksheet(wb,'S8. Coloc GTEx')
  
  writeData(wb,'S1. Lead variants',leadVariants)
  writeData(wb,'S3. MR CD4 Cells',soskic.mr)
  writeData(wb,'S4. MR BQC19',bqc.mr)
  writeData(wb,'S5. MR GTEx',gtex.mr)
  writeData(wb,'S6. Coloc CD4 Cells',soskic.coloc)
  writeData(wb,'S7. Coloc BQC19',bqc.hits)
  writeData(wb,'S8. Coloc GTEx',gtex.hits)
  
  saveWorkbook(wb,'Supplemental_Tables.xlsx')
}


