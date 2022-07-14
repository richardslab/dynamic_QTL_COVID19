getHGIGwasData = function(out) { #take in string for outcome
  if (out == 'A2') hgi = vroom(glue('HGI_Data/COVID19_HGI_A2_ALL_eur_leave23andme_20220403.tsv.gz'),show_col_types = F)
  else if (out == 'B2') hgi = vroom(glue('HGI_Data/COVID19_HGI_B2_ALL_eur_leave23andme_20220403.tsv.gz'),show_col_types = F)
  else if (out == 'C2') hgi = vroom(glue('HGI_Data/COVID19_HGI_C2_ALL_eur_leave23andme_20220403.tsv.gz'),show_col_types = F)
  return(hgi)
}
getAllMRDataReady.MRthenColoc = function(dataset,type.of.qtl) {
  loci = read.table('subpositions_chrpos.txt') ; names(loci) = c('Chr','Pos','Outcome')
  loci %<>% dplyr::arrange(Outcome)
  df.mr = data.frame(Label=character(),nSnps=numeric(),Beta=numeric(),SE=numeric(),
                     Pval=numeric(),WMedianBeta=numeric(),WModeBeta=numeric(),
                     EggerBeta=numeric(),EggerInterceptP=numeric(),SteigerP=numeric())
  for (locus in 1:nrow(loci)) {
    mr.data.list = mr.data.label = list()
    outcome = loci$Outcome[[locus]]
    
    print(glue('On locus {locus} of {nrow(loci)}'))
    hgi.for.outcome = getHGIGwasData(outcome)
    curr.chr = loci$Chr[[locus]] ; curr.pos = loci$Pos[[locus]]
    hgiQTL = getHGIData(hgi.for.outcome,chr=curr.chr,lead.snp.position=curr.pos,MR=T) 
    
    if (dataset == 'soskic') {
      cells = c('CD4_Memory_uns_0h','CD4_Memory_stim_16h','CD4_Memory_stim_40h',
                'CD4_Memory_stim_5d','CD4_Naive_uns_0h','CD4_Naive_stim_16h',
                'CD4_Naive_stim_40h','CD4_Naive_stim_5d','TN_0h','TN_16h',
                'TN_40h','TN_5d','TCM_0h','TCM_16h','TCM_40h','TCM_5d',
                'TEM_0h','TEM_16h','TEM_40h',
                'TEM_5d','TEMRA_0h','TEMRA_16h','TEMRA_40h','TEMRA_5d',
                'nTreg_0h','nTreg_16h','nTreg_40h')
      for (cell in cells) {
        print(glue('On cell {cell}'))
        files = list.files('Soskic_eQTL_summ_stats',pattern=cell,full.names = T)
        files = files[which(str_detect(files,'awk.addedcolumns.txt'))]
        
        for (file in files) {
          full.mQTL = vroom(file,show_col_types=F) %>% dplyr::filter(CHR==curr.chr) %>% dplyr::filter(POS >= curr.pos-500000,POS <= curr.pos+500000)
          all.genes = unique(full.mQTL[['phenotype_id']])
          
          for (gene in all.genes) {
            mQTL = full.mQTL %>% dplyr::filter(phenotype_id == gene) %>% 
              mutate(N=round(ma_count/maf*(1/2)))
            mQTL %<>% dplyr::rename(beta=slope,se=slope_se,pvalues=pval_nominal,MAF=maf,
                                    chr=CHR,position=POS) %>% mutate(varbeta=se^2) %>%
              dplyr::select(chr,position,REF,ALT,beta,varbeta,pvalues,MAF,variant_id,phenotype_id,ma_samples,N,se) %>%
              dplyr::filter(!is.na(variant_id),!duplicated(variant_id),!is.na(varbeta),
                            pvalues < 1) #coloc freaks out if pvalue = 1
            mQTL %<>% mutate(variant_id = glue('chr{mQTL$chr}:{mQTL$position}:{mQTL$REF}:{mQTL$ALT}'))
            print(paste(nrow(mQTL),nrow(hgiQTL)))
            plt.data = list(mQTL,hgiQTL)
            label = glue('{outcome}_{cell}_{curr.chr}_{curr.pos}_{gene}')
            
            mr.data.list[[length(mr.data.list)+1]] = plt.data
            mr.data.label %<>% append(label)
          }
        }
      }
      return(list(mr.data.list,mr.data.label))
      saveRDS(list(mr.data.list,mr.data.label),glue('MRthenColoc_MRData/{locus}_{outcome}_{curr.chr}_{curr.pos}.rds'))
    } else if (dataset == 'bqc') {
      states = c('_inf_No_','_inf_Yes_')
      for (state in states) {
        if (type.of.qtl == 'eQTL') files = list.files('BQC_eQTL_Data',pattern=state,full.names = T)
        files = files[which(str_detect(files,'awk.addedcolumns.txt'))]
        files = files[which(str_detect(files,glue('\\.chr{curr.chr}\\.')))]
        for (file in files) {
          print(glue('Curr file: {file}'))
          str.chr = glue('chr{curr.chr}')
          full.mQTL = vroom(file,show_col_types=F) %>% dplyr::filter(CHR==str.chr) %>% dplyr::filter(POS >= curr.pos-500000,POS <= curr.pos+500000)
          all.genes = unique(full.mQTL[['gene_id']])
          for (gene in all.genes) {
            mQTL = full.mQTL %>% dplyr::filter(gene_id == gene) %>% 
              mutate(N=round(ma_count/maf*(1/2)))
            mQTL %<>% dplyr::rename(beta=slope,se=slope_se,pvalues=pval_nominal,MAF=maf,
                                    chr=CHR,position=POS) %>% mutate(varbeta=se^2) %>%
              dplyr::select(chr,position,REF,ALT,beta,varbeta,pvalues,MAF,variant_id,gene_id,ma_samples,N,se) %>%
              dplyr::filter(!is.na(variant_id),!duplicated(variant_id),!is.na(varbeta),
                            pvalues < 1) #coloc freaks out if pvalue = 1
            mQTL %<>% mutate(variant_id = glue('{mQTL$chr}:{mQTL$position}:{mQTL$REF}:{mQTL$ALT}'))
            plt.data = list(mQTL,hgiQTL)
            label = glue('{outcome}_{state}_{curr.chr}_{curr.pos}_{gene}')
            
            mr.data.list[[length(mr.data.list)+1]] = plt.data
            mr.data.label %<>% append(label)
          }
        }
      }
      saveRDS(list(mr.data.list,mr.data.label),glue('MRthenColoc_MRData_BQC/{locus}_{outcome}_{curr.chr}_{curr.pos}.rds'))
    } else if (dataset == 'gtex') {
      if (type.of.qtl == 'eQTL') files = list.files('GTEx_eQTL_Data',full.names = T)
      files = files[which(str_detect(files,'awk.addedcolumns.txt'))]
      files = files[which(str_detect(files,glue('\\.chr{curr.chr}\\.')))]
      for (file in files) {
        print(glue('Curr file: {file}'))
        str.chr = glue('chr{curr.chr}')
        full.mQTL = vroom(file,show_col_types=F) %>% dplyr::filter(CHR==str.chr) %>% dplyr::filter(POS >= curr.pos-500000,POS <= curr.pos+500000)
        all.genes = unique(full.mQTL[['phenotype_id']])
        for (gene in all.genes) {
          mQTL = full.mQTL %>% dplyr::filter(phenotype_id == gene) %>% 
            mutate(N=round(ma_count/maf*(1/2)))
          mQTL %<>% dplyr::rename(beta=slope,se=slope_se,pvalues=pval_nominal,MAF=maf,
                                  chr=CHR,position=POS) %>% mutate(varbeta=se^2) %>%
            dplyr::select(chr,position,REF,ALT,beta,varbeta,pvalues,MAF,variant_id,phenotype_id,ma_samples,N,se) %>%
            dplyr::filter(!is.na(variant_id),!duplicated(variant_id),!is.na(varbeta),
                          pvalues < 1) #coloc freaks out if pvalue = 1
          mQTL %<>% mutate(variant_id = glue('{mQTL$chr}:{mQTL$position}:{mQTL$REF}:{mQTL$ALT}'))
          plt.data = list(mQTL,hgiQTL)
          label = glue('{outcome}_{curr.chr}_{curr.pos}_{gene}')
          
          mr.data.list[[length(mr.data.list)+1]] = plt.data
          mr.data.label %<>% append(label)
        }
      }
      saveRDS(list(mr.data.list,mr.data.label),glue('MRthenColoc_MRData_GTEx/{locus}_{outcome}_{curr.chr}_{curr.pos}.rds'))
    }
  }
}
conductMR.MRthenColoc = function(dataset,type.of.qtl,rev,recollect) {
  to.return = data.frame(Label=character(),nSnps=numeric(),Beta=numeric(),SE=numeric(),
                         Pval=numeric(),WMedianBeta=numeric(),WModeBeta=numeric(),
                         EggerBeta=numeric(),EggerInterceptP=numeric(),SteigerP=numeric(),
                         LeaveOneOutP=numeric())
  if (dataset == 'soskic') path = 'MRthenColoc_MRData'
  else if (dataset == 'bqc' & type.of.qtl == 'eQTL') path = 'MRthenColoc_MRData_BQC_eQTL'
  else if (dataset == 'gtex' & type.of.qtl == 'eQTL') path = 'MRthenColoc_MRData_GTEx_eQTL'
  
  if (!rev) files = sort(list.files(path,pattern='rds',full.names=T),decreasing=F)
  else files = sort(list.files(path,pattern='rds',full.names=T),decreasing=T)
  
  for (file in 1:length(files)) { #already removed data for loci in MHC
    missing.rsids = vroom('missing_rsids.txt',show_col_types = FALSE)
    if (str_detect(files[[file]],'mr_results_all.rds')) next
    else if (str_detect(files[[file]],'coloc_results_all.rds')) next
    else if (str_detect(files[[file]],'coloc_sens_all.rds')) next
    else if (str_detect(files[[file]],'heatmap_data.rds')) next
    print(glue('Curr file: {files[[file]]}'))
    f.check = sub('.*\\/','',files[[file]])
    f.check = sub('\\.rds','',f.check)
    f.check = glue('{f.check}_mrresults.rds')
    file.processed = (f.check %in% list.files(glue('{path}/MR_Results')))
    if (file.processed & !recollect) {
      tmp = readRDS(glue('{path}/MR_Results/{f.check}'))
      to.return %<>% add_row(tmp)
      next
    }
    
    curr.results = readRDS(files[[file]])
    tmp.curr.mr = data.frame(Label=character(),nSnps=numeric(),Beta=numeric(),SE=numeric(),
                             Pval=numeric(),WMedianBeta=numeric(),WModeBeta=numeric(),
                             EggerBeta=numeric(),EggerInterceptP=numeric(),SteigerP=numeric(),
                             LeaveOneOutP=numeric())
    
    for (item in 1:length(curr.results[[1]])) {
      print(glue('Doing MR {item} of {length(curr.results[[1]])}, file: {files[[file]]}'))
      common.qtls = curr.results[[1]][[item]] ; rsids = common.qtls[[2]][['rsid']]
      
      #prep exposure data. 
      exp.rsid.df = data.frame(VID=common.qtls[[1]][['variant_id']],Rsid=NA)
      num.overlap = 0
      for (item2 in 1:nrow(exp.rsid.df)) {
        if (exp.rsid.df$VID[[item2]] %in% common.qtls[[2]][['variant_id']]) {
          exp.rsid.df$Rsid[[item2]] = rsids[which(common.qtls[[2]][['variant_id']]==exp.rsid.df$VID[[item2]])]
          num.overlap = num.overlap + 1
        }else if (exp.rsid.df$VID[[item2]] %in% missing.rsids$VID) {
          exp.rsid.df$Rsid[[item2]] = missing.rsids$RSID[which(missing.rsids$VID == exp.rsid.df$VID[[item2]])]
        }else{ #ie missing rsid
          curr.row = common.qtls[[1]][item2,]
          rs.error = F
          rs = tryCatch(get_rsid_from_position(chrom=curr.row$chr,pos=curr.row$position,ref=curr.row$REF,
                                      alt=curr.row$ALT,assembly='hg38'),
                        error = function(c) {print(glue('Missing rsid for: {curr.row$variant_id}')) ; rs.error <<- T})
          if (!rs.error) missing.rsids %<>% add_row(VID = curr.row$variant_id,RSID = rs)
          exp.rsid.df$Rsid[[item2]] = rs
        }
      }
      print(glue('Num remaining missing rsids in exposure: {length(which(is.na(exp.rsid.df$Rsid)))}'))
      print(glue('On test {nrow(tmp.curr.mr)} in {files[[file]]}'))
      vroom_write(missing.rsids,'missing_rsids.txt')

      if (num.overlap < 50) {print('Fewer than 50 variants shared') ; next} #follow reasoning by soskic
      
      exposure.data = data.frame(chr=common.qtls[[1]][['chr']],position=common.qtls[[1]][['position']],
                                 beta=common.qtls[[1]][['beta']],se=common.qtls[[1]][['se']],
                                 SNP=exp.rsid.df$Rsid,effect_allele=common.qtls[[1]][['ALT']],
                                 other_allele=common.qtls[[1]][['REF']],eaf=common.qtls[[1]][['MAF']],
                                 samplesize=common.qtls[[1]][['N']],pval=common.qtls[[1]][['pvalues']])
      exposure.data = format_data(exposure.data,type='exposure')
      exposure.data = clump_data(exposure.data,pop='EUR')
      if (length(which(is.na(exposure.data$SNP))>0)) print('Missing an rsid for exposure data, need proxy')
      
      #prep outcome data
      outcome.data = data.frame(chr=common.qtls[[2]][['CHR']],position=common.qtls[[2]][['position']],
                                beta=common.qtls[[2]][['beta']],se=common.qtls[[2]][['se']],
                                SNP=rsids,effect_allele=common.qtls[[2]][['ALT']],
                                other_allele=common.qtls[[2]][['REF']],eaf=common.qtls[[2]][['MAF']],
                                samplesize=common.qtls[[2]][['N']],pval=common.qtls[[2]][['pvalues']])
      outcome.data = format_data(outcome.data,type='outcome')
      
      #harmonize and Steiger
      harmonized = harmonise_data(exposure.data,outcome.data)
      if (nrow(harmonized) < nrow(exposure.data)) {
        print('Lost SNPs when harmonized, need to get proxy')
        rsid.non.overlap = exposure.data[['SNP']][which(exposure.data[['SNP']] %notin% outcome.data[['SNP']])]
        proxy = getProxy(rsid.non.overlap)
        if (proxy == 'MissingRs') {
          print(glue('Missing proxy for {rsid.non.overlap}'))
          return(outcome.data)
        }
      }
      steiger = directionality_test(harmonized)
      
      #do mr functions
      if (sum(harmonized$mr_keep)==0) print('No SNP remaining after harmonization')
      else if (sum(harmonized$mr_keep)==1) {
        result = mr(harmonized, method_list = "mr_wald_ratio")
        tmp.curr.mr %<>% add_row(Label=curr.results[[2]][[item]],nSnps=result$nsnp,Beta=result$b,
                                 SE=result$se,Pval=result$pval,WMedianBeta=NA,WModeBeta=NA,
                                 EggerBeta=NA,EggerInterceptP=NA,SteigerP=steiger$steiger_pval,
                                 LeaveOneOutP=NA)
      }else if (sum(harmonized$mr_keep)>1) {
        result <- mr(harmonized, method_list = c("mr_ivw","mr_two_sample_ml","mr_weighted_median","mr_penalised_weighted_median","mr_weighted_mode"))
        leaveoneout = mr_leaveoneout(harmonized)
        if (sum(harmonized$mr_keep)>2) {
          egger_result <- mr_egger_regression(harmonized$beta.exposure[harmonized$mr_keep],
                                              harmonized$beta.outcome[harmonized$mr_keep],
                                              harmonized$se.exposure[harmonized$mr_keep],
                                              harmonized$se.outcome[harmonized$mr_keep])
          
          tmp.curr.mr %<>% add_row(Label=curr.results[[2]][[item]],nSnps=result$nsnp[1],Beta=result$b[1],
                                   SE=result$se[1],Pval=result$pval[1],WMedianBeta=result$b[3],
                                   WModeBeta=result$b[5],EggerBeta=egger_result$b,
                                   EggerInterceptP=egger_result$pval_i,SteigerP=steiger$steiger_pval,
                                   LeaveOneOutP=leaveoneout$p)
        } else {
          tmp.curr.mr %<>% add_row(Label=curr.results[[2]][[item]],nSnps=result$nsnp[1],Beta=result$b[1],
                                   SE=result$se[1],Pval=result$pval[1],WMedianBeta=NA,WModeBeta=NA,
                                   EggerBeta=NA,EggerInterceptP=NA,SteigerP=steiger$steiger_pval,
                                   LeaveOneOutP=leaveoneout$p)
        }
        if (result$pval[1] <= 0.05) return(harmonized)
      }
    }
    to.return %<>% add_row(tmp.curr.mr)
    f.str = sub("\\.rds","",files[[file]])
    f.str = sub('.*\\/','',f.str)
    saveRDS(tmp.curr.mr,glue('{path}/MR_Results/{f.str}_mrresults.rds'))
  }
  p.adj = p.adjust(to.return$Pval,method='BH')
  to.return %<>% mutate(PvalAdj = p.adj,PassSens = NA)
  to.return %<>% dplyr::mutate(Outcome=NA,CHR=NA,POS=NA,Gene=NA,.after=Label)
  if (dataset != 'gtex') to.return %<>% dplyr::mutate(CellState=NA,.after=Outcome)

  for (row in 1:nrow(to.return)) {
    curr.row = to.return[row,]
    if (is.na(curr.row$WMedianBeta)) 
      to.return$PassSens[[row]] = (curr.row$PvalAdj<=0.05 & curr.row$SteigerP <= 0.05)
    else
      to.return$PassSens[[row]] = (curr.row$PvalAdj<=0.05 & curr.row$SteigerP <= 0.05 &
                                     sign(curr.row$Beta)==sign(curr.row$WMedianBeta) &
                                     sign(curr.row$Beta)==sign(curr.row$WModeBeta) &
                                     sign(curr.row$Beta)==sign(curr.row$EggerBeta) &
                                     curr.row$EggerInterceptP > 0.05)
    if (str_detect(curr.row$Label,'A2')) to.return$Outcome[[row]] = 'A2'
    else if (str_detect(curr.row$Label,'B2')) to.return$Outcome[[row]] = 'B2'
    else if (str_detect(curr.row$Label,'C2')) to.return$Outcome[[row]] = 'C2'
    
    to.return$Gene[[row]] = sub('.*ENSG','ENSG',curr.row$Label)
    s = sub('.*__','',curr.row$Label) ; s = sub('_ENSG.*','',s)
    if (dataset == 'bqc') {
      if (str_detect(curr.row$Label,'inf_No')) to.return$CellState[[row]] = 'inf_No'
      else if (str_detect(curr.row$Label,'inf_Yes')) to.return$CellState[[row]] = 'inf_Yes'
      to.return$CHR[[row]] = sub('_.*','',s)
      to.return$POS[[row]] = sub('.*_','',s)
    }else if (dataset == 'soskic') {
      spl = str_split(s,'_')[[1]]
      if (str_detect(spl[[2]],'CD4')) {
        to.return$CellState[[row]] = paste(spl[2:5],collapse='_')
        to.return$CHR[[row]] = spl[[6]]
        to.return$POS[[row]] = spl[[7]]
      } else {
        to.return$CellState[[row]] = paste(spl[2:3],collapse='_')
        to.return$CHR[[row]] = spl[[4]]
        to.return$POS[[row]] = spl[[5]]
      }
    }else if (dataset == 'gtex') {
      spl = str_split(s,'_')[[1]]
      to.return$CHR[[row]] = spl[[2]]
      to.return$POS[[row]] = spl[[3]]
    }
  }
  
  return(to.return %>% dplyr::arrange(Outcome,CHR,POS,Gene))
}
colocFromMR = function(dataset,type.of.qtl,mr.data,sens.data) {
  sens.tests.list = list()
  coloc.mr.df = mr.data[1,] %>% dplyr::filter(Gene == '') #just get columns
  if (dataset == 'soskic') {
    path = 'MRthenColoc_MRData'
    df = data.frame(Cell=character(),Outcome=character(),CHR=numeric(),POS=numeric(),
                    H0.PP=numeric(),
                    H1.PP=numeric(),H2.PP=numeric(),H3.PP=numeric(),
                    H4.PP=numeric(),Gene=character())
  }else if (dataset == 'bqc' & type.of.qtl == 'eQTL') {
    path = 'MRthenColoc_MRData_BQC_eQTL'
    df = data.frame(InfState=character(),Outcome=character(),CHR=numeric(),POS=numeric(),
                    Gene=character(),H0.PP=numeric(),
                    H1.PP=numeric(),H2.PP=numeric(),H3.PP=numeric(),
                    H4.PP=numeric())
  }else if (dataset == 'gtex' & type.of.qtl == 'eQTL') {
    path = 'MRthenColoc_MRData_GTEx_eQTL'
    df = data.frame(Outcome=character(),CHR=numeric(),POS=numeric(),
                    Gene=character(),H0.PP=numeric(),
                    H1.PP=numeric(),H2.PP=numeric(),H3.PP=numeric(),
                    H4.PP=numeric())
  }
  files = list.files(path)
  
  sig.mr = mr.data %>% dplyr::filter(PassSens)
  for (row in 1:nrow(sig.mr)) {
    print(glue('On row {row} of {nrow(sig.mr)}'))
    curr.row = sig.mr[row,]
    match.sens = data.frame()
    if (!is.null(sens.data)) { #if checking sens
      if (dataset == 'soskic') 
        match.sens = sens.data %>% dplyr::filter(Cell==curr.row$CellState,
                                                 Outcome==curr.row$Outcome,
                                                 CHR==curr.row$CHR,POS==curr.row$POS,
                                                 Gene==curr.row$Gene,H4.PP>=0.8)
      else if (dataset=='bqc') 
        match.sens = sens.data %>% dplyr::filter(InfState==curr.row$CellState,
                                                 Outcome==curr.row$Outcome,
                                                 CHR==curr.row$CHR,POS==curr.row$POS,
                                                 Gene==curr.row$Gene,H4.PP>=0.8)
      else if (dataset == 'gtex')
        match.sens = sens.data %>% dplyr::filter(Outcome==curr.row$Outcome,
                                                 CHR==curr.row$CHR,POS==curr.row$POS,
                                                 Gene==curr.row$Gene,H4.PP>=0.8)
      if (nrow(match.sens)==0) {next}
    }
    file = files[which(str_detect(files,curr.row$Outcome) & 
                         str_detect(files,glue('_{curr.row$CHR}_')) &
                         str_detect(files,glue('_{curr.row$POS}')))]
    curr.results = readRDS(glue('{path}/{file}'))

    if (dataset == 'soskic') 
      label.index = which(str_detect(curr.results[[2]],curr.row$Outcome) &
                            str_detect(curr.results[[2]],curr.row$CellState) &
                            str_detect(curr.results[[2]],curr.row$CHR) &
                            str_detect(curr.results[[2]],curr.row$POS) &
                            str_detect(curr.results[[2]],curr.row$Gene))
    else if (dataset == 'bqc') 
      label.index = which(str_detect(curr.results[[2]],curr.row$Outcome) &
                            str_detect(curr.results[[2]],curr.row$CellState) &
                            str_detect(curr.results[[2]],curr.row$CHR) &
                            str_detect(curr.results[[2]],curr.row$POS) &
                            str_detect(curr.results[[2]],curr.row$Gene))
    else if (dataset == 'gtex')
      label.index = which(str_detect(curr.results[[2]],curr.row$Outcome) &
                            str_detect(curr.results[[2]],curr.row$CHR) &
                            str_detect(curr.results[[2]],curr.row$POS) &
                            str_detect(curr.results[[2]],curr.row$Gene))
    if (length(label.index) > 1) stop('Two matches for data')
    
    curr.data = curr.results[[1]][[label.index]]
    common.qtls.list = getCommonQTLLists(curr.data[[1]],curr.data[[2]] %>% mutate(rsid=NULL)) #get rid of rsid as messes up colocalization
    if (length(common.qtls.list[[1]][['pvalues']]) == 0) { print(glue('No overlapping rows {file}')) ; break }
    if (length(common.qtls.list[[1]][['pvalues']]) < 50) {
      if (!is.null(vector.specific)) print(glue('Not enough variant On locus {item} of {nrow(loci)}, cell {cell}'))
      next #follow criteria per Soskic
    }
    coloc.res = coloc.abf(common.qtls.list[[1]],common.qtls.list[[2]])
    if (!is.null(sens.data)) { sens.tests.list[[length(sens.tests.list)+1]] = coloc.res }
    coloc.res = coloc.res$summary
    if (dataset == 'soskic') df %<>% add_row(Cell=curr.row$CellState,Outcome=curr.row$Outcome,
                    CHR=as.numeric(curr.row$CHR),POS=as.numeric(curr.row$POS),
                    H0.PP=coloc.res['PP.H0.abf'],H1.PP=coloc.res['PP.H1.abf'],
                    H2.PP=coloc.res['PP.H2.abf'],H3.PP=coloc.res['PP.H3.abf'],
                    H4.PP=coloc.res['PP.H4.abf'],Gene=curr.row$Gene)
    else if (dataset == 'bqc') 
      df %<>% add_row(InfState=curr.row$CellState,Outcome=curr.row$Outcome,
                      CHR=as.numeric(curr.row$CHR),POS=as.numeric(curr.row$POS),
                      Gene=curr.row$Gene,
                      H0.PP=coloc.res['PP.H0.abf'],H1.PP=coloc.res['PP.H1.abf'],
                      H2.PP=coloc.res['PP.H2.abf'],H3.PP=coloc.res['PP.H3.abf'],
                      H4.PP=coloc.res['PP.H4.abf'])
    else if (dataset == 'gtex')
      df %<>% add_row(Outcome=curr.row$Outcome,
                      CHR=as.numeric(curr.row$CHR),POS=as.numeric(curr.row$POS),
                      Gene=curr.row$Gene,
                      H0.PP=coloc.res['PP.H0.abf'],H1.PP=coloc.res['PP.H1.abf'],
                      H2.PP=coloc.res['PP.H2.abf'],H3.PP=coloc.res['PP.H3.abf'],
                      H4.PP=coloc.res['PP.H4.abf'])
    if (coloc.res['PP.H4.abf'] >= 0.8) coloc.mr.df %<>% add_row(curr.row)
  }
  if (is.null(sens.data)) return(df)
  else return(list(sens.tests.list,df,coloc.mr.df))
}
doSensTesting = function(data,sens.data,index) {
  print(glue('Num cases strong colocalization: {nrow(data %>% subset(H4.PP>=0.8))}'))
  strong.coloc = (data %>% dplyr::filter(H4.PP >= 0.8))[index,]
  print(paste(strong.coloc$Outcome,strong.coloc$Cell,strong.coloc$InfState,strong.coloc$Gene))
  sensitivity(sens.data[[index]],rule='H4 > 0.5')
}
plotSignificantMR = function(data.source,outcome,mr.coloc.data,drop.indices) { #drop MR results failing coloc sensitivity testing
  if (data.source == 'soskic') {
    df = mr.coloc.data[-drop.indices,] %>% dplyr::filter(Outcome == outcome) %>% 
      mutate(Gene = replaceGeneNames(Gene),CellTimeGene=paste(CellState,Gene)) %>%
      mutate(CellTimeGene = str_replace(CellTimeGene,'stim_','')) %>%
      mutate(CellTimeGene = str_replace_all(CellTimeGene,'_',' '))
  }else if (data.source == 'bqc') {
    mr.coloc.data$CellState = str_replace(mr.coloc.data$CellState,'inf_No','Control')
    mr.coloc.data$CellState = str_replace(mr.coloc.data$CellState,'inf_Yes','Infectious')
    df = mr.coloc.data[-drop.indices,] %>% dplyr::filter(Outcome %in% outcome) %>% 
      mutate(Gene = replaceGeneNames(Gene),OutcomeStateGene=paste(Outcome,CellState,Gene))
  }else if (data.source == 'gtex') {
    df = mr.coloc.data[-drop.indices,] %>% dplyr::filter(Outcome %in% outcome) %>% 
      mutate(Gene = replaceGeneNames(Gene),OutcomeGene=paste(Outcome,Gene))
  }
  
  if (data.source == 'soskic') {
    plt = ggplot(df,aes(y=CellTimeGene,x=exp(Beta),xmin=exp(Beta-1.96*SE),xmax=exp(Beta+1.96*SE))) +
      scale_y_discrete(limits=rev) + xlim(0.7,2)
  }else if (data.source == 'bqc') {
    df %<>% mutate(OutcomeStateGene=str_replace(OutcomeStateGene,'A2','Severe'))
    df %<>% mutate(OutcomeStateGene=str_replace(OutcomeStateGene,'B2','Hospitalized'))
    df %<>% mutate(OutcomeStateGene=str_replace(OutcomeStateGene,'C2','Susceptibility'))
    plt = ggplot(df,aes(y=OutcomeStateGene,x=exp(Beta),xmin=exp(Beta-1.96*SE),xmax=exp(Beta+1.96*SE))) +
      scale_y_discrete(limits=rev(c('Severe Control HIP1','Severe Control IFNAR2','Severe Control JD275616',
                                    'Hospitalized Control ABO','Hospitalized Infectious IFNAR2',
                                    'Susceptibility Control NAPSB')))
  }else if (data.source == 'gtex') {
    df %<>% mutate(OutcomeGene=str_replace(OutcomeGene,'A2','Severe'))
    df %<>% mutate(OutcomeGene=str_replace(OutcomeGene,'B2','Hospitalized'))
    df %<>% mutate(OutcomeGene=str_replace(OutcomeGene,'C2','Susceptibility'))
    plt = ggplot(df,aes(y=OutcomeGene,x=exp(Beta),xmin=exp(Beta-1.96*SE),xmax=exp(Beta+1.96*SE))) +
      scale_y_discrete(limits=rev) + scale_x_continuous(breaks=c(0.7,0.8,0.9,1.0))
  }
  
  plt = plt + geom_point(size=4) + geom_errorbar(width=0.2) + xlab('Odds ratio') +
    theme_bw() + theme(text=element_text(size=28),legend.position="none") +
    geom_vline(xintercept=1,linetype='dotted') + ylab('')

  print(plt)
  return(df)
}
getCellandTime = function(cell.state) {
  spl = str_split(cell.state,'_')[[1]]
  if (str_detect(cell.state,'CD4_Memory')) { cell = 'CD4_Memory' ; time = spl[[4]] }
  else if (str_detect(cell.state,'CD4_Naive')) { cell = 'CD4_Naive' ; time = spl[[4]] }
  else if (str_detect(cell.state,'TN_')) { cell = 'TN' ; time = spl[[2]] }
  else if (str_detect(cell.state,'TCM_')) { cell = 'TCM_' ; time = spl[[2]] }
  else if (str_detect(cell.state,'TEM_')) { cell = 'TEM' ; time = spl[[2]] }
  return(c(cell,time))
}
getHeatmapData = function(data.source,data) { #get coloc data for all timepoints for hits
  if (data.source=='soskic')
    df = data.frame(Cell=character(),Outcome=character(),CHR=numeric(),POS=numeric(),
                    H0.PP=numeric(),
                    H1.PP=numeric(),H2.PP=numeric(),H3.PP=numeric(),
                    H4.PP=numeric(),Gene=character())
  else if (data.source == 'bqc')
    df = data.frame(InfState=character(),Outcome=character(),CHR=numeric(),POS=numeric(),
                    H0.PP=numeric(),
                    H1.PP=numeric(),H2.PP=numeric(),H3.PP=numeric(),
                    H4.PP=numeric(),Gene=character())
  
  for (item in 1:nrow(data)) {
    curr.row = data[item,]
    if (data.source == 'soskic') {
      cell.time = getCellandTime(curr.row$CellState)
      times = c('0h','16h','40h','5d')
      for (time in times) {
        print(glue('On data row {item} of {nrow(data)}, time {time}'))
        if (cell.time[[1]] == 'TN' | cell.time[[1]] == 'TEM') 
          curr.cell = glue('{cell.time[[1]]}_{time}') #get underscore from fxn
        else curr.cell = cell.time[[1]]

        coloc.df = soskicColoc(c(curr.row$Outcome,curr.row$CHR,curr.row$POS,
                                 curr.cell,glue('_{time}_'),curr.row$Gene),
                               allCells=F,returnCommonQTLs=F,mr=F,specific.getdf=T)
        df %<>% add_row(coloc.df %>% mutate('H3+H4'=NULL))
      }
    } else if (data.source == 'bqc') {
      states = c('inf','noninf')
      for (state in states) {
        print(glue('On data row {item} of {nrow(data)}, state {state}'))
        coloc.df = bqcColoc('eQTL',c(curr.row$Outcome,curr.row$CHR,curr.row$POS,
                                     state,curr.row$Gene),
                            returnCommonQTLs=F,mr=F,specific.getdf=T)
        df %<>% add_row(coloc.df %>% mutate('H3+H4'=NULL))
      }
    }
  }
  return(df)
}
plotColocHeatmap = function(dataset,mr.plot.data,coloc.data) {
  new.plot.data = mr.plot.data  
  new.coloc.data = coloc.data %>% mutate(Gene=replaceGeneNames(Gene))
  if (dataset=='soskic') {
    new.plot.data %<>% mutate(H4.0h=NA,H4.16h=NA,H4.40h=NA,H4.5d=NA)
    plot.df = data.frame(CellGene=character(),Time=character(),H4.PP=numeric())
  } else if (dataset == 'bqc') {
    new.plot.data %<>% mutate(H4.Noninf=NA,H4.Inf=NA)
    plot.df = data.frame(OutcomeGene=character(),State=character(),H4.PP=numeric())
  }
  
  for (row in 1:nrow(new.plot.data)) {
    curr.row = new.plot.data[row,]
    if (row==1) print(glue('Outcome for data: {curr.row$Outcome}')) #avoid losing track of data
    if (dataset=='soskic') {
      curr.cell = getCellandTime(curr.row$CellState)[[1]]
      matched.coloc = new.coloc.data %>% 
        dplyr::filter(Outcome==curr.row$Outcome,Gene==curr.row$Gene,
                      str_detect(Cell,curr.cell)) %>% distinct()
      H4_0 = which(str_detect(matched.coloc$Cell,'_0h_'))
      H4_16 = which(str_detect(matched.coloc$Cell,'_16h_'))
      H4_40 = which(str_detect(matched.coloc$Cell,'_40h_'))
      H4_5 = which(str_detect(matched.coloc$Cell,'_5d_'))

      if (length(H4_0) > 0) {
        new.plot.data[[row,'H4.0h']] = matched.coloc[H4_0,'H4.PP']
        plot.df %<>% add_row(CellGene=paste(curr.cell,curr.row$Gene),
                             Time='0h',H4.PP=matched.coloc[H4_0,'H4.PP'])
      }
      if (length(H4_16) > 0) {
        new.plot.data[[row,'H4.16h']] = matched.coloc[H4_16,'H4.PP']
        plot.df %<>% add_row(CellGene=paste(curr.cell,curr.row$Gene),
                             Time='16h',H4.PP=matched.coloc[H4_16,'H4.PP'])
      }
      if (length(H4_40) > 0) { 
        new.plot.data[[row,'H4.40h']] = matched.coloc[H4_40,'H4.PP']
        plot.df %<>% add_row(CellGene=paste(curr.cell,curr.row$Gene),
                             Time='40h',H4.PP=matched.coloc[H4_40,'H4.PP'])
      }
      if (length(H4_5) > 0) {
        new.plot.data[[row,'H4.5d']] = matched.coloc[H4_5,'H4.PP']
        plot.df %<>% add_row(CellGene=paste(curr.cell,curr.row$Gene),
                             Time='5d',H4.PP=matched.coloc[H4_5,'H4.PP'])
      }
    }else if (dataset == 'bqc') {
      matched.coloc = new.coloc.data %>% 
        dplyr::filter(Outcome==curr.row$Outcome,Gene==curr.row$Gene) %>% 
        distinct()
      H4_noninf = which(matched.coloc$InfState == 'noninf')
      H4_inf = which(matched.coloc$InfState == 'inf')

      if (length(H4_noninf) > 0) {
        new.plot.data[[row,'H4.Noninf']] = matched.coloc[H4_noninf,'H4.PP']
        plot.df %<>% add_row(OutcomeGene=paste(curr.row$Outcome,curr.row$Gene),
                             State='Noninf',H4.PP=matched.coloc[H4_noninf,'H4.PP'])
      }
      if (length(H4_inf) > 0) {
        new.plot.data[[row,'H4.Inf']] = matched.coloc[H4_inf,'H4.PP']
        plot.df %<>% add_row(OutcomeGene=paste(curr.row$Outcome,curr.row$Gene),
                             State='Inf',H4.PP=matched.coloc[H4_inf,'H4.PP'])
      }
    }
  }
  if (dataset=='soskic') {
    plot.df %<>% mutate(CellGene = str_replace_all(CellGene,'_',' '))
    plt = ggplot(plot.df,aes(x=Time,y=CellGene,fill=H4.PP)) + scale_y_discrete(limits=rev)
  }else if (dataset == 'bqc') {
    plot.df %<>% mutate(OutcomeGene = str_replace(OutcomeGene,'A2','Severe'))
    plot.df %<>% mutate(OutcomeGene = str_replace(OutcomeGene,'B2','Hospitalized'))
    plot.df %<>% mutate(OutcomeGene = str_replace(OutcomeGene,'C2','Susceptibility'))
    plt = ggplot(plot.df,aes(x=State,y=OutcomeGene,fill=H4.PP)) +
      scale_y_discrete(limits=rev(c('Severe HIP1','Severe IFNAR2','Severe JD275616',
                                    'Hospitalized ABO','Hospitalized IFNAR2','Susceptibility NAPSB'))) +
      scale_x_discrete(limits=c('Noninf','Inf'))
  }
  plt = plt + geom_tile(aes(linetype=(H4.PP>=0.8)),show.legend=F) + theme_bw() + theme(text=element_text(size=28)) +
            ylab('') + scale_fill_continuous(limits=c(0,1)) +
    geom_text(aes(label=round(H4.PP,2)),color='white',size=10)
  print(plt)

  return(new.plot.data)
}