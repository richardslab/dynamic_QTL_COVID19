getHGIData = function(HGI,chr,lead.snp.position,MR) { #process HGI data for coloc
  print('Reading hgi data')
  names(HGI)[1] = 'CHR'

  hgi = HGI %>% subset(CHR == chr) %>% mutate(varbeta=all_inv_var_meta_sebeta^2) %>% 
    dplyr::filter(!duplicated(SNP),!is.na(varbeta),varbeta!=0) %>%
    dplyr::rename(position=POS,beta=all_inv_var_meta_beta,MAF=all_meta_AF,
                  pvalues=all_inv_var_meta_p,se=all_inv_var_meta_sebeta) %>%
    mutate(variant_id=paste0('chr',SNP))
  hgi %<>% mutate(N=all_inv_var_meta_cases + all_inv_var_meta_controls)
  if (!MR) hgi %<>% mutate(rsid=NULL)
  
  if (!is.null(lead.snp.position)) hgi = hgi %>% 
    dplyr::filter(position >= lead.snp.position-500000,position<=lead.snp.position+500000)
  
  return(hgi %>% dplyr::filter(!is.na(MAF)))
}
getCommonQTLLists = function(mqtl,hgi) { #process data for coloc
  common.snps = inner_join(mqtl,hgi,by='variant_id') %>% dplyr::select(variant_id)
  mqtl.common = inner_join(mqtl,common.snps,by='variant_id') %>% arrange(position,pvalues) %>%
    dplyr::filter(varbeta!=0)
  hgi.common = inner_join(hgi,common.snps,by='variant_id') %>% arrange(position,pvalues) %>%
    dplyr::filter(varbeta!=0)

  mqtl.list = as.list(mqtl.common)
  mqtl.list$snp=mqtl.list$variant_id
  mqtl.list$type = 'quant'
  
  hgi.list = as.list(hgi.common)
  hgi.list$type = 'cc'
  hgi.list$snp = hgi.list$variant_id
  hgi.list$s = hgi.list$all_inv_var_meta_cases/(hgi.list$all_inv_var_meta_cases + hgi.list$all_inv_var_meta_controls)

  return(list(mqtl.list,hgi.list))
}
makeTxtAllSoskicParquets = function() { #process soskic data for my pipe
  files = list.files('Soskic_eQTL_summ_stats',pattern='.parquet',full.names=F)
  files = str_replace(files,'.parquet','')
  for (item in files) {
    print(glue('On item {which(files == item)} of {length(files)}'))
    curr.file = read_parquet(glue('Soskic_eQTL_summ_stats/{item}.parquet'))
    vroom_write(curr.file,glue('Soskic_eQTL_summ_stats/{item}.txt'),col_names=T,
                quote="none")
  }
}
soskicColoc = function(vector.specific,allCells,returnCommonQTLs,mr,specific.getdf) { #coloc for Soskic et al data, permitting looking at specific entry
  df = data.frame(Cell=character(),Outcome=character(),CHR=numeric(),POS=numeric(),
                  H0.PP=numeric(),
                  H1.PP=numeric(),H2.PP=numeric(),H3.PP=numeric(),
                  H4.PP=numeric(),Gene=character())
  loci = read.table('subpositions_chrpos.txt') ; names(loci) = c('chr','pos','outcome')

  if (is.null(vector.specific)) {
    hgi.a2 = vroom(glue('HGI_Data/COVID19_HGI_A2_ALL_eur_leave23andme_20220403.tsv.gz'),show_col_types = F)
    hgi.b2 = vroom(glue('HGI_Data/COVID19_HGI_B2_ALL_eur_leave23andme_20220403.tsv.gz'),show_col_types = F)
    hgi.c2 = vroom(glue('HGI_Data/COVID19_HGI_C2_ALL_eur_leave23andme_20220403.tsv.gz'),show_col_types = F)
  }else{
    if (vector.specific[[1]] == 'A2') hgi.a2 = vroom(glue('HGI_Data/COVID19_HGI_A2_ALL_eur_leave23andme_20220403.tsv.gz'),show_col_types = F)
    else if (vector.specific[[1]] == 'B2') hgi.b2 = vroom(glue('HGI_Data/COVID19_HGI_B2_ALL_eur_leave23andme_20220403.tsv.gz'),show_col_types = F)
    else if (vector.specific[[1]] == 'C2') hgi.c2 = vroom(glue('HGI_Data/COVID19_HGI_C2_ALL_eur_leave23andme_20220403.tsv.gz'),show_col_types = F)
  }
  for (item in 1:nrow(loci)) {
    if (is.null(vector.specific)) print(glue('On row {item} of {nrow(loci)}'))
    curr.chr = loci[['chr']][[item]] ; targ.pos = loci[['pos']][[item]]
    outcome = loci[['outcome']][[item]]
    if (!is.null(vector.specific)) {
      if (outcome != vector.specific[[1]]) next
      else if (curr.chr != vector.specific[[2]] | targ.pos != vector.specific[[3]]) next
    }
    print(glue('Investigating outcome {outcome} chr{curr.chr} pos{targ.pos}'))
    
    if (outcome == 'A2') hgiQTL = getHGIData(hgi.a2,chr=curr.chr,lead.snp.position=targ.pos,MR=mr) #grch38
    else if (outcome == 'B2') hgiQTL = getHGIData(hgi.b2,chr=curr.chr,lead.snp.position=targ.pos,MR=mr) #grch38
    else if (outcome == 'C2') hgiQTL = getHGIData(hgi.c2,chr=curr.chr,lead.snp.position=targ.pos,MR=mr) #grch38
    
    if (!allCells) cells = c('CD4_Memory','CD4_Naive','TN_0h','TN_16h','TN_40h','TN_5d',
                             'TCM_','TEM_0h','TEM_16h','TEM_40h','TEM_5d','TEMRA_',
                             'nTreg_')
    else cells = c('CD4_Memory','CD4_Naive','HSP_','nTreg_','T_ER-stress_','TCM_','TEM_HLApositive_',
                   'TEM_LA_','TEMRA_','TM_cycling_','TM_ER-stress_','TN_cycling_','TN_HSP_','TN_IFN_',
                   'TN_LA_','TN_NFKB_','TN2_','TN_0h','TN_5d','TN_16h','TN_40h','TEM_0h','TEM_5d',
                   'TEM_16h','TEM_40h')
    for (cell in cells) { #go by cell
      if (!is.null(vector.specific)) { #select for cell
        if (cell != vector.specific[[4]]) next
      }
      files = list.files('Soskic_eQTL_summ_stats',pattern=cell,full.names = T)
      files = files[which(str_detect(files,'awk.addedcolumns.txt'))]
      for (file in files) {
        if (!is.null(vector.specific)) { #select for time
          if (!str_detect(file,vector.specific[[5]])) next
        }
        #get mqtl
        full.mQTL = vroom(file,show_col_types=F) %>% dplyr::filter(CHR==curr.chr) %>% dplyr::filter(POS >= targ.pos-500000,POS <= targ.pos+500000)
        if (nrow(full.mQTL) < 1) { print('No rows in full mQTL') ; next }
        all.genes = unique(full.mQTL[['phenotype_id']])
        if (!is.null(vector.specific)) { #select for gene
          if (vector.specific[[6]] %in% all.genes) print('Gene in the given mQTL data')
          else print('Gene not in mQTL data, ie not expressed')
        }
        for (gene in all.genes) {
          if (!is.null(vector.specific)) {
            if (gene != vector.specific[[6]]) { next }
          }
          mQTL = full.mQTL %>% dplyr::filter(phenotype_id == gene) %>% 
            mutate(N=round(ma_count/maf*(1/2)))
          mQTL %<>% dplyr::rename(beta=slope,se=slope_se,pvalues=pval_nominal,MAF=maf,
                                  chr=CHR,position=POS) %>% mutate(varbeta=se^2) %>%
            dplyr::select(chr,position,REF,ALT,beta,varbeta,pvalues,MAF,variant_id,phenotype_id,ma_samples,N,se) %>%
            dplyr::filter(!is.na(variant_id),!duplicated(variant_id),!is.na(varbeta),
                          pvalues < 1) #coloc freaks out if pvalue = 1
          mQTL %<>% mutate(variant_id = glue('chr{mQTL$chr}:{mQTL$position}:{mQTL$REF}:{mQTL$ALT}'))
          if (nrow(mQTL) == 0) { print(glue('No rows left {gene}')) ; next}
          
          common.qtls.list = getCommonQTLLists(mQTL,hgiQTL)
          if (length(common.qtls.list[[1]][['pvalues']]) == 0) { print(glue('No overlapping rows {file}')) ; break }
          if (length(common.qtls.list[[1]][['pvalues']]) < 50) {
            if (!is.null(vector.specific)) print(glue('Not enough variant On locus {item} of {nrow(loci)}, cell {cell}'))
            next #follow criteria per Soskic
          }

          if (returnCommonQTLs & !is.null(vector.specific)) return(list(mQTL,hgiQTL))
          coloc.res = coloc.abf(common.qtls.list[[1]],common.qtls.list[[2]])
          if (!is.null(vector.specific) & !specific.getdf) {
            print(glue('Overlapped rows: {nrow(common.qtls.list[[1]])}'))
            sens = sensitivity(coloc.res,rule='H4 > 0.5')
            return(coloc.res)
          }
          coloc.res = coloc.res$summary
          df %<>% add_row(Cell=file,Outcome=outcome,CHR=curr.chr,POS=targ.pos,
                          H0.PP=coloc.res['PP.H0.abf'],H1.PP=coloc.res['PP.H1.abf'],
                          H2.PP=coloc.res['PP.H2.abf'],H3.PP=coloc.res['PP.H3.abf'],
                          H4.PP=coloc.res['PP.H4.abf'],Gene=gene)
          if (!is.null(vector.specific) & specific.getdf) return(df)
          if (is.null(vector.specific)) print(glue('On locus {item} of {nrow(loci)}, cell {cell}'))
        }
      }
    }
  }
  return(df %>% add_column('H3+H4'=df$H3.PP+df$H4.PP))
}
bqcColoc = function(type.of.qtl,vector.specific,returnCommonQTLs,mr,specific.getdf) { #bqc coloc, used for single entries or all entries
  df = data.frame(InfState=character(),Outcome=character(),CHR=numeric(),POS=numeric(),
                  H0.PP=numeric(),
                  H1.PP=numeric(),H2.PP=numeric(),H3.PP=numeric(),
                  H4.PP=numeric(),Gene=character())
  loci = read.table('subpositions_chrpos.txt') #rel7
  
  names(loci) = c('chr','pos','outcome')
  print('Reading hgiQTL gwas')
  
  if (is.null(vector.specific)) {
    hgi.a2 = vroom(glue('HGI_Data/COVID19_HGI_A2_ALL_eur_leave23andme_20220403.tsv.gz'),show_col_types = F)
    hgi.b2 = vroom(glue('HGI_Data/COVID19_HGI_B2_ALL_eur_leave23andme_20220403.tsv.gz'),show_col_types = F)
    hgi.c2 = vroom(glue('HGI_Data/COVID19_HGI_C2_ALL_eur_leave23andme_20220403.tsv.gz'),show_col_types = F)
  }else{
    if (vector.specific[[1]] == 'A2') hgi.a2 = vroom(glue('HGI_Data/COVID19_HGI_A2_ALL_eur_leave23andme_20220403.tsv.gz'),show_col_types = F)
    else if (vector.specific[[1]] == 'B2') hgi.b2 = vroom(glue('HGI_Data/COVID19_HGI_B2_ALL_eur_leave23andme_20220403.tsv.gz'),show_col_types = F)
    else if (vector.specific[[1]] == 'C2') hgi.c2 = vroom(glue('HGI_Data/COVID19_HGI_C2_ALL_eur_leave23andme_20220403.tsv.gz'),show_col_types = F)
  }
  
  for (item in 1:nrow(loci)) {
    if (is.null(vector.specific)) print(glue('On row {item} of {nrow(loci)}'))
    curr.chr = loci[['chr']][[item]] ; targ.pos = loci[['pos']][[item]]
    outcome = loci[['outcome']][[item]]
    if (!is.null(vector.specific)) { #match for outcome when using vector specific
      if (outcome != vector.specific[[1]]) next
      else if (curr.chr != vector.specific[[2]] | targ.pos != vector.specific[[3]]) next
    }
    print(glue('Investigating outcome {outcome} chr{curr.chr} pos{targ.pos}'))
    
    if (outcome == 'A2') hgiQTL = getHGIData(hgi.a2,chr=curr.chr,lead.snp.position=targ.pos,MR=mr) #grch38
    else if (outcome == 'B2') hgiQTL = getHGIData(hgi.b2,chr=curr.chr,lead.snp.position=targ.pos,MR=mr) #grch38
    else if (outcome == 'C2') hgiQTL = getHGIData(hgi.c2,chr=curr.chr,lead.snp.position=targ.pos,MR=mr) #grch38
    
    inf.state = c('inf','noninf')
    
    for (state in inf.state) {
      if (!is.null(vector.specific)) { #select for state
        if (state != vector.specific[[4]]) next
      }
      if (type.of.qtl == 'eQTL') full.mQTL = vroom(glue('BQC_eQTL_Data/Subpositions/chr{curr.chr}_{targ.pos}_window_{state}.txt'),show_col_types=F)
      else if (type.of.qtl == 'sQTL') full.mQTL = vroom(glue('BQC_sQTL_Data/Subpositions/chr{curr.chr}_{targ.pos}_window_{state}.txt'),show_col_types=F)
      
      #get mqtl
      full.mQTL %<>% dplyr::filter(POS >= targ.pos-500000,POS <= targ.pos+500000)
      if (nrow(full.mQTL) < 1) { print('No remaining rows in full mQTL') ; next }
      all.genes = unique(full.mQTL[['gene_id']])
      for (gene in all.genes) {
        if (!is.null(vector.specific)) {
          if (type.of.qtl == 'eQTL' & !str_detect(gene,vector.specific[[5]])) {next} #given bqc using transcript ids as decimal
          else if (type.of.qtl == 'sQTL' & gene != vector.specific[[5]]) next
        }
        print('getting mQTL')
        mQTL = full.mQTL %>% dplyr::filter(gene_id == gene) %>% 
          mutate(N=round(ma_count/maf*(1/2)))
        if (nrow(mQTL) == 0) { print(glue('No rows {gene}')) ; next}
        mQTL %<>% dplyr::rename(beta=slope,se=slope_se,pvalues=pval_nominal,MAF=maf,
                                chr=CHR,position=POS,phenotype_id=gene_id) %>% mutate(varbeta=se^2) %>%
          dplyr::select(chr,position,REF,ALT,beta,varbeta,pvalues,MAF,variant_id,phenotype_id,ma_samples,N,se) %>%
          dplyr::filter(!is.na(variant_id),!duplicated(variant_id),!is.na(varbeta),
                        pvalues < 1,varbeta != 0) #coloc freaks out if pvalue = 1
        mQTL %<>% mutate(variant_id = glue('{mQTL$chr}:{mQTL$position}:{mQTL$REF}:{mQTL$ALT}'))
        if (nrow(mQTL) == 0) { print(glue('No rows {gene}')) ; next}
        
        common.qtls.list = getCommonQTLLists(mQTL,hgiQTL)
        if (returnCommonQTLs & !is.null(vector.specific)) return(list(mQTL,hgiQTL))
        
        if (length(common.qtls.list[[1]][['pvalues']]) == 0) { print(glue('No overlapping rows')) ; break }
        if (length(common.qtls.list[[1]][['pvalues']]) < 50) {
          print(glue('Not enough variant On locus {item} of {nrow(loci)}, infstate {state}'))
          next #follow criteria per Soskic
        }
        error.next=F
        coloc.res = tryCatch(coloc.abf(common.qtls.list[[1]],common.qtls.list[[2]]),
                             error = function(c) {print(glue('Error for {gene} for loci {item} - likely no expression')) ; error.next <<- T})
        if (error.next) next
        if (!is.null(vector.specific) & !specific.getdf) {
          print(glue('Overlapped rows: {nrow(common.qtls.list[[1]])}'))
          sens = sensitivity(coloc.res,rule='H4 > 0.5')
          return(coloc.res)
        }
        coloc.res = coloc.res$summary
        df %<>% add_row(InfState=state,Outcome=outcome,CHR=curr.chr,POS=targ.pos,
                        H0.PP=coloc.res['PP.H0.abf'],H1.PP=coloc.res['PP.H1.abf'],
                        H2.PP=coloc.res['PP.H2.abf'],H3.PP=coloc.res['PP.H3.abf'],
                        H4.PP=coloc.res['PP.H4.abf'],Gene=gene)
        print(glue('On locus {item} of {nrow(loci)}, inf.state {state}'))
      }
    }
  }
  return(df %>% add_column('H3+H4'=df$H3.PP+df$H4.PP))
}
gtexColoc = function(type.of.qtl,vector.specific,returnCommonQTLs,mr) { #gtex coloc, permitting specific coloc or do all
  df = data.frame(Outcome=character(),CHR=numeric(),POS=numeric(),
                  H0.PP=numeric(),
                  H1.PP=numeric(),H2.PP=numeric(),H3.PP=numeric(),
                  H4.PP=numeric(),LeadSNP_eQTL=character(),LeadSNP_HGI=character(),Gene=character(),
                  Lead_eQTL_Beta=numeric(),Lead_eQTL_P=numeric(),
                  Lead_HGI_Beta=numeric(),Lead_HGI_P=numeric(),
                  HGI_MAF=numeric())
  loci = read.table('subpositions_chrpos.txt') %>% dplyr::arrange(V3) #rel7
  
  names(loci) = c('chr','pos','outcome')
  print('Reading hgiQTL gwas')
  
  if (is.null(vector.specific)) {
    hgi.a2 = vroom(glue('HGI_Data/COVID19_HGI_A2_ALL_eur_leave23andme_20220403.tsv.gz'),show_col_types = F)
    hgi.b2 = vroom(glue('HGI_Data/COVID19_HGI_B2_ALL_eur_leave23andme_20220403.tsv.gz'),show_col_types = F)
    hgi.c2 = vroom(glue('HGI_Data/COVID19_HGI_C2_ALL_eur_leave23andme_20220403.tsv.gz'),show_col_types = F)
  }else{
    if (vector.specific[[1]] == 'A2') hgi.a2 = vroom(glue('HGI_Data/COVID19_HGI_A2_ALL_eur_leave23andme_20220403.tsv.gz'),show_col_types = F)
    else if (vector.specific[[1]] == 'B2') hgi.b2 = vroom(glue('HGI_Data/COVID19_HGI_B2_ALL_eur_leave23andme_20220403.tsv.gz'),show_col_types = F)
    else if (vector.specific[[1]] == 'C2') hgi.c2 = vroom(glue('HGI_Data/COVID19_HGI_C2_ALL_eur_leave23andme_20220403.tsv.gz'),show_col_types = F)
  }
  
  for (item in 1:nrow(loci)) {
    if (is.null(vector.specific)) print(glue('On row {item} of {nrow(loci)}'))
    curr.chr = loci[['chr']][[item]] ; targ.pos = loci[['pos']][[item]]
    outcome = loci[['outcome']][[item]]
    if (!is.null(vector.specific)) {
      if (outcome != vector.specific[[1]]) next
      else if (curr.chr != vector.specific[[2]] | targ.pos != vector.specific[[3]]) next
    }
    print(glue('Investigating outcome {outcome} chr{curr.chr} pos{targ.pos}'))
    
    if (outcome == 'A2') hgiQTL = getHGIData(hgi.a2,chr=curr.chr,lead.snp.position=targ.pos,outcome,MR=mr) #grch38
    else if (outcome == 'B2') hgiQTL = getHGIData(hgi.b2,chr=curr.chr,lead.snp.position=targ.pos,outcome,MR=mr) #grch38
    else if (outcome == 'C2') hgiQTL = getHGIData(hgi.c2,chr=curr.chr,lead.snp.position=targ.pos,outcome,MR=mr) #grch38
    
    if (type.of.qtl == 'eQTL') full.mQTL = vroom(glue('GTEx_eQTL_Data/Subpositions/chr{curr.chr}_{targ.pos}_window.txt'),show_col_types=F)
    else if (type.of.qtl == 'sQTL') full.mQTL = vroom(glue('GTEx_sQTL_Data/Subpositions/chr{curr.chr}_{targ.pos}_window.txt'),show_col_types=F)
    #get mqtl
    # full.mQTL %<>% dplyr::filter(POS >= targ.pos-500000,POS <= targ.pos+500000) #files already filtered
    if (nrow(full.mQTL) < 1) { print('No remaining rows in full mQTL') ; next }
    all.genes = unique(full.mQTL[['phenotype_id']])
    for (gene in all.genes) {
      if (!is.null(vector.specific)) {
        if (!str_detect(gene,vector.specific[[5]])) next #given bqc using transcript ids as decimal
      }
      mQTL = full.mQTL %>% dplyr::filter(phenotype_id == gene) %>% 
        mutate(N=round(ma_count/maf*(1/2)))
      if (nrow(mQTL) == 0) { print(glue('No rows {gene}')) ; next}
      mQTL %<>% dplyr::rename(beta=slope,se=slope_se,pvalues=pval_nominal,MAF=maf,
                              chr=CHR,position=POS) %>% mutate(varbeta=se^2) %>%
        dplyr::select(chr,position,REF,ALT,beta,varbeta,pvalues,MAF,variant_id,phenotype_id,ma_samples,N,se) %>%
        dplyr::filter(!is.na(variant_id),!duplicated(variant_id),!is.na(varbeta),
                      pvalues < 1) #coloc freaks out if pvalue = 1
      mQTL %<>% mutate(variant_id = glue('{mQTL$chr}:{mQTL$position}:{mQTL$REF}:{mQTL$ALT}'))
      if (nrow(mQTL) == 0) { print(glue('No rows {gene}')) ; next}
      
      common.qtls.list = getCommonQTLLists(mQTL,hgiQTL,F,outcome)
      if (returnCommonQTLs & !is.null(vector.specific)) return(list(mQTL,hgiQTL))
      
      if (length(common.qtls.list[[1]][['pvalues']]) == 0) { print(glue('No overlapping rows')) ; break }
      if (length(common.qtls.list[[1]][['pvalues']]) < 50) {
        print(glue('Not enough variant On locus {item} of {nrow(loci)}'))
        next #follow criteria per Soskic
      }
      error.next=F
      coloc.res = tryCatch(coloc.abf(common.qtls.list[[1]],common.qtls.list[[2]]),
                           error = function(c) {print(glue('Error for {gene} for loci {item} - likely no expression')) ; error.next <<- T})
      if (error.next) next
      if (!is.null(vector.specific)) {
        print(glue('Overlapped rows: {nrow(common.qtls.list[[1]])}'))
        sens = sensitivity(coloc.res,rule='H4 > 0.5')
        return(coloc.res)
      }
      lead.eqtl.index = which(common.qtls.list[[1]][['pvalues']] == min(common.qtls.list[[1]][['pvalues']]))[[1]]
      lead.hgi.index = which(common.qtls.list[[2]][['pvalues']] == min(common.qtls.list[[2]][['pvalues']]))[[1]]
      lead.snp.eqtl = common.qtls.list[[1]][['variant_id']][[lead.eqtl.index]]
      lead.snp.hgi = common.qtls.list[[2]][['variant_id']][[lead.hgi.index]]
      eQTL.beta = common.qtls.list[[1]][['beta']][[lead.eqtl.index]]
      eQTL.P = common.qtls.list[[1]][['pvalues']][[lead.eqtl.index]]
      HGI.beta = common.qtls.list[[2]][['beta']][[lead.hgi.index]]
      HGI.P = common.qtls.list[[2]][['pvalues']][[lead.hgi.index]]
      HGI.maf = common.qtls.list[[2]][['MAF']][[lead.hgi.index]]
      coloc.res = coloc.res$summary
      df %<>% add_row(Outcome=outcome,CHR=curr.chr,POS=targ.pos,
                      H0.PP=coloc.res['PP.H0.abf'],H1.PP=coloc.res['PP.H1.abf'],
                      H2.PP=coloc.res['PP.H2.abf'],H3.PP=coloc.res['PP.H3.abf'],
                      H4.PP=coloc.res['PP.H4.abf'],LeadSNP_eQTL=lead.snp.eqtl,LeadSNP_HGI=lead.snp.hgi,
                      Gene=gene,
                      Lead_eQTL_Beta=eQTL.beta,Lead_eQTL_P=eQTL.P,
                      Lead_HGI_Beta=HGI.beta,Lead_HGI_P=HGI.P,
                      HGI_MAF=HGI.maf)
      print(glue('On locus {item} of {nrow(loci)}'))
    }
  }
  return(df %>% add_column('H3+H4'=df$H3.PP+df$H4.PP))
}
doAllSensTests = function(data.source,type.of.data,data) { # do all sens tests when did coloc on all entries rather than a single one
  #collect the data so sensitivity tests can be done more quickly
  #skip chr6 because all loci are in MHC, that is not done in MR.
  #as an aside, even when chr6 is done, single causal variant assumption violated
  arranged.data = data %>% dplyr::arrange(Outcome,CHR,POS) %>% 
    dplyr::filter(H4.PP>=0.8,CHR != 6)
  plot.list = list()
  plot.label = list()
  if (data.source == 'soskic') {
    for (item in 1:nrow(arranged.data)) {
      print(glue('Working on plot {item} of {nrow(arranged.data %>% dplyr::filter(H4.PP>=0.8))}'))
      curr.row = arranged.data[item,]
      cell = gsub('.*/','',curr.row$Cell) ; cell = gsub('_500.*','',cell)
      cell.split = str_split(cell,'_')[[1]]
      if (length(cell.split) > 2) curr.cell = glue('{cell.split[[1]]}_{cell.split[[2]]}')
      else curr.cell = cell.split[[1]]
      time = cell.split[[length(cell.split)]] ; 
      if (curr.cell == 'TEM') curr.cell = glue('{curr.cell}_{time}')
      else if (curr.cell == 'TN') curr.cell = glue('{curr.cell}_{time}')
      else if (curr.cell == 'TCM') curr.cell = 'TCM_'
      time = glue('_{time}_')
      # tmp = soskicColoc(vector.specific=c('A2',9,133271182,'TEM_16h','_16h_','ENSG00000160271'),allCells=T,returnCommonQTLs = F,mr=F)
      print(paste(curr.row$Outcome,curr.row$CHR,curr.row$POS,curr.cell,time,curr.row$Gene))
      plt.data = soskicColoc(vector.specific=c(curr.row$Outcome,curr.row$CHR,curr.row$POS,curr.cell,time,curr.row$Gene),allCells=F,returnCommonQTLs = F,mr=F)
      plot.list[[item]] = plt.data
      plot.label[[item]] = glue('{curr.row$Outcome}_{curr.cell}_{time}_{curr.row$CHR}_{curr.row$POS}_{curr.row$Gene}')
    }
  }else if (data.source == 'bqc') {
    for (item in 1:nrow(arranged.data)) {
      print(glue('Working on plot {item} of {nrow(arranged.data %>% dplyr::filter(H4.PP>=0.8))}'))
      curr.row = arranged.data[item,]
      plt.data = bqcColoc(type.of.data,vector.specific = c(curr.row$Outcome,curr.row$CHR,curr.row$POS,curr.row$InfState,curr.row$Gene),returnCommonQTLs=F,mr=F)
      plot.list[[item]] = plt.data
      plot.label[[item]] = glue('{curr.row$Outcome}_{curr.row$CHR}_{curr.row$POS}_{curr.row$Gene}')
    }
  }else if (data.source == 'gtex') {
    for (item in 1:nrow(arranged.data)) {
      print(glue('Working on plot {item} of {nrow(arranged.data %>% dplyr::filter(H4.PP>=0.8))}'))
      curr.row = arranged.data[item,]
      plt.data = gtexColoc(type.of.data,c(curr.row$Outcome,curr.row$CHR,curr.row$POS,NA,curr.row$Gene),returnCommonQTLs=F,mr=F,onlySoskicLoci=F)
      plot.list[[item]] = plt.data
      plot.label[[item]] = glue('{curr.row$Outcome}_{curr.row$CHR}_{curr.row$POS}_{curr.row$Gene}')
    }
  }
  return(list(plot.list,plot.label))
}
plotDynamicPP = function(data,cell,outcome,chromosomes,exclude.genes) { #plotting H4.PP when doing all coloc
  plot.data = data %>% dplyr::filter(str_detect(data$Cell,cell),Outcome==outcome,CHR %in% chromosomes)
  plot.data$Gene = replaceGeneNames(plot.data$Gene)
  
  #ensure focused on exact sample
  if (cell == 'TEM_') plot.data %<>% dplyr::filter(!str_detect(plot.data$Cell,'TEM_HLA'),
                                                  !str_detect(plot.data$Cell,'TEM_LA'))
  if (cell == 'TCM_') plot.data %<>% dplyr::filter(!str_detect(plot.data$Cell,'TCM_LA'))
  if (cell == 'TN_') plot.data %<>% dplyr::filter(!str_detect(plot.data$Cell,'TN_IFN'),
                                                   !str_detect(plot.data$Cell,'TN_cycling'),
                                                  !str_detect(plot.data$Cell,'TN_HSP'),
                                                  !str_detect(plot.data$Cell,'TN_LA'),
                                                  !str_detect(plot.data$Cell,'TN_NFKB'))
  # colocalizing.genes = unique(plot.data[['Gene']][which(plot.data$H4.PP >= 0.8)])
  colocalizing.genes = c('CAT','NAPSA','ACSF3','RALGDS','RAB2A','ADAM15') 
  plot.data %<>% dplyr::filter(Gene %in% colocalizing.genes)
  times = c('0h','16h','40h','5d')
  plot.df = data.frame(Time=character(),PP=numeric(),Gene=character(),Position=numeric())
  for (time in times) {
    time.rows = plot.data %>% dplyr::filter(str_detect(plot.data$Cell,glue('_{time}')))
    plot.df %<>% add_row(Time = time,PP=time.rows$H4.PP,Gene=time.rows$Gene,
                         Position=time.rows$POS)
  }
  if (cell == 'CD4_Naive') cell.title = 'CD4 Antigen Naive'
  else if (cell == 'CD4_Memory') cell.title = 'CD4 Memory'
  else if (cell == 'T_ER-stress') cell.title = 'T_ER-stress'
  else if (cell == 'TEM_') cell.title = 'TEM'
  else if (cell == 'TN_HSP') cell.title = 'TN_HSP'
  else if (cell == 'TCM_') cell.title = 'TCM'
  else if (cell == 'TN_cycling') cell.title = 'TN_cycling'
  else if (cell == 'TN_IFN') cell.title = 'TN_IFN'
  else if (cell == 'TN_') cell.title = 'TN'
  else cell.title = 'Missing Cell Title'
  if (outcome == 'A2') text.outcome = 'Severe'
  else if (outcome == 'B2') text.outcome = 'Hospitalized'
  else if (outcome == 'C2') text.outcome = 'Susceptible'
  
  plot.df %<>% dplyr::filter(!(Gene %in% exclude.genes))
  
  group.colors = c(CAT = '#F8766D',NAPSA = '#00BFC4',ACSF3='#00BA38',
                   RALGDS = '#B79F00',RAB2A = '#619CFF',ADAM15='#F564E3')
  group.shapes = c(CAT = 16,NAPSA = 17,ACSF3=15,RALGDS=3,RAB2A=7,ADAM15=8)
  # groups.to.remove = numeric() ; genes.of.interest = c('CAT','NAPSA','ACSF3','RALGDS','RAB2A','ADAM15')
  # for (g in 1:length(genes.of.interest)) {
  #   if (genes.of.interest[[g]] %notin% plot.df$Gene) groups.to.remove %<>% append(g)
  # }
  # group.colors = group.colors[-groups.to.remove]
  # group.shapes = group.shapes[-groups.to.remove]

  plot = ggplot(plot.df,aes(x=Time,y=PP,group=Gene)) + geom_line(aes(color=Gene)) + theme_bw() +
          geom_point(aes(shape=Gene),size=6) + theme(text=element_text(size=40)) + geom_hline(yintercept=0.8,linetype='dotted') +
          #ggtitle(glue('Cell: {cell.title}. Outcome: {text.outcome}')) +
          scale_x_discrete(limits=c('0h','16h','40h','5d')) + ylim(c(0,1))
  plot = plot + scale_color_manual(values=group.colors) +
    scale_shape_manual(values=group.shapes)
  print(plot)
  return(plot.df)
}
plotDynamicPP.BQC = function(data,outcome,type.of.qtl,chromosomes,exclude.genes) {
  plot.data = data %>% dplyr::filter(Outcome==outcome,CHR %in% chromosomes)

  colocalizing.genes = unique(plot.data[['Gene']][which(plot.data$H4.PP >= 0.8)])
  plot.data %<>% dplyr::filter(Gene %in% colocalizing.genes)
  plot.df = data.frame(State=character(),PP=numeric(),Gene=character(),Position=numeric())

  plot.df %<>% add_row(State = plot.data$InfState,PP=plot.data$H4.PP,Gene=plot.data$Gene,
                       Position=plot.data$POS)
  
  if (outcome == 'A2') text.outcome = 'Severe'
  else if (outcome == 'B2') text.outcome = 'Hospitalized'
  else if (outcome == 'C2') text.outcome = 'Susceptible'
  
  plot.df$Gene = replaceGeneNames(plot.df$Gene)
  
  if (type.of.qtl == 'eQTL') { #not done for sQTLs because too many variants - made figures hard to appreciate
    colocalizing.genes = c('ABO','HIP1','IFNAR2','JD275616','NAPSA','RAB2A',
                           'GBAP1','NAPSB','THBS3','TYK2','ENSG00000279996.1')
    plot.df %<>% dplyr::filter(Gene %in% colocalizing.genes)
    hues = hue_pal()(length(colocalizing.genes))
    group.colors = c(ABO=hues[1],HIP1=hues[2],IFNAR2=hues[3],
                     JD275616=hues[4],NAPSA=hues[5],RAB2A=hues[6],
                     GBAP1=hues[7],NAPSB=hues[8],THBS3=hues[9],TYK2=hues[10],
                     ENSG00000279996.1=hues[11])
    group.shapes = c(ABO=16,HIP1=17,IFNAR2=15,JD275616=3,NAPSA=7,
                     RAB2A=8,GBAP1=1,NAPSB=2,THBS3=5,TYK2=6,ENSG00000279996.1=9)
  }
    
  plot.df %<>% dplyr::filter(Gene %notin% exclude.genes)
  if (outcome=='A2') plot.df %<>% dplyr::filter(!(Gene=='IFNAR2' & Position==33969937))
  if (outcome=='A2') plot.df %<>% dplyr::filter(!(Gene=='IFNAR2 Splice Variant A' & Position==33969937))
  if (outcome=='A2') plot.df %<>% dplyr::filter(!(Gene=='IFNAR2 Splice Variant B' & Position==33969937))
  if (outcome=='A2') plot.df %<>% dplyr::filter(!(Gene=='IFNAR2 Splice Variant E' & Position==33969937))
  if (outcome=='A2') plot.df %<>% dplyr::filter(!(Gene=='IL10RB Splice Variant A' & Position==33969937))
  if (outcome=='B2') plot.df %<>% dplyr::filter(!(Gene=='IFNAR2' & Position==33940612))
  if (outcome=='B2') plot.df %<>% dplyr::filter(!(Gene=='IFNAR2 Splice Variant A' & Position==33940612))
  if (outcome=='B2') plot.df %<>% dplyr::filter(!(Gene=='IFNAR2 Splice Variant B' & Position==33940612))
  if (outcome=='B2') plot.df %<>% dplyr::filter(!(Gene=='IL10RB Splice Variant A' & Position==33940612))
  
  plot = ggplot(plot.df,aes(x=State,y=PP,group=Gene)) + theme_bw() +
          geom_line(aes(color=Gene)) + geom_point(aes(shape=Gene),size=3) +
          theme(text=element_text(size=24)) + geom_hline(yintercept=0.8,linetype='dotted') +
          # ggtitle(glue('Outcome: {text.outcome}')) +
          ylim(c(0,1)) + scale_x_discrete(limits=c('noninf','inf'))
  if (type.of.qtl == 'eQTL') plot = plot + scale_color_manual(values=group.colors) +
    scale_shape_manual(values=group.shapes)
  print(plot)
  return(plot.df)
}
plotDynamicPP.GTEx = function(data,outcome,type.of.qtl,chromosomes,exclude.genes) {
  plot.data = data %>% dplyr::filter(Outcome==outcome,CHR %in% chromosomes)
  
  colocalizing.genes = unique(plot.data[['Gene']][which(plot.data$H4.PP >= 0.8)])
  plot.data %<>% dplyr::filter(Gene %in% colocalizing.genes)
  plot.df = data.frame(State=character(),PP=numeric(),Gene=character(),Position=numeric())
  
  plot.df %<>% add_row(State='Resting',PP=plot.data$H4.PP,Gene=plot.data$Gene,
                       Position=plot.data$POS)
  
  if (outcome == 'A2') text.outcome = 'Severe'
  else if (outcome == 'B2') text.outcome = 'Hospitalized'
  else if (outcome == 'C2') text.outcome = 'Susceptible'
  
  plot.df$Gene = replaceGeneNames(plot.df$Gene)
  
  plot.df %<>% dplyr::filter(!(Gene %in% exclude.genes))
  
  if (outcome=='A2') plot.df %<>% dplyr::filter(!(Gene=='IFNAR2 Splice Variant E' & Position==33969937))
  if (outcome=='A2') plot.df %<>% dplyr::filter(!(Gene=='IFNAR2 Splice Variant F' & Position==33969937))
  if (outcome=='B2') plot.df %<>% dplyr::filter(!(Gene=='IFNAR2 Splice Variant F' & Position==33940612))
  
  if (type.of.qtl == 'eQTL') { #not done for sQTLs because there are too many - makes a hard to interpret plot
    colocalizing.genes = c('ABO','HIP1','NAPSA','NTN5','ENSG00000236263',
                           'ICAM5','OAS3','RAVER1','TYK2')
    plot.df %<>% dplyr::filter(Gene %in% colocalizing.genes)
    hues = hue_pal()(length(colocalizing.genes))
    group.colors = c(ABO=hues[1],HIP1=hues[2],NAPSA=hues[3],
                     NTN5=hues[4],ENSG00000236263=hues[5],ICAM5=hues[6],
                     OAS3=hues[7],RAVER1=hues[8],TYK2=hues[9])
    group.shapes = c(ABO=16,HIP1=17,NAPSA=15,NTN5=3,ENSG00000236263=7,
                     ICAM5=8,OAS3=1,RAVER1=2,TYK2=5)
  }
  
  plot = ggplot(plot.df,aes(x=State,y=PP,group=Gene)) + theme_bw() +
          geom_point(aes(shape=Gene),size=3) +
          theme(text=element_text(size=24)) + geom_hline(yintercept=0.8,linetype='dotted') +
          # ggtitle(glue('Outcome: {text.outcome}')) +
          ylim(c(0,1))
  if (type.of.qtl == 'eQTL') plot = plot + scale_color_manual(values=group.colors) +
    scale_shape_manual(values=group.shapes)
  print(plot)
  return(plot.df)
}
#########################################
#proteomics functions
getSomaData = function(analyte) { 
  #get somascan column name of interest, use medNormRef file
  seq.id = getSEQID(analyte)
  
  soma.batch1 = read.adat('~/SS-200150_v4_ACDPlasma.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.medNormRefSMP.adat') %>% 
    dplyr::select(SubjectID,seq.id,PlateId)
  soma.batch2 = read.adat('~/SS-215174.hybNorm.medNormInt.plateScale.medNormRefSMP.adat') %>%
    dplyr::select(SubjectID,seq.id,PlateId)
  all.soma = soma.batch1 %>% add_row(soma.batch2) %>% drop_na(SubjectID) %>%
    distinct(SubjectID,.keep_all=T)
  
  #take care of misnaming
  misnamed.rows = which(str_detect(all.soma[['SubjectID']],'CHUM'))
  for (row in misnamed.rows) { all.soma[['SubjectID']][[row]] = glue('VAP{str_replace(all.soma[["SubjectID"]][[row]],"_CHUM","")}') }
  
  #and get the bqcid
  vcode.bqcid.dict = readxl::read_xlsx('~/Vcode_mapping_Release6.xlsx',sheet='All') %>%
    dplyr::select(bqc_id,plasma_vcode...21) %>% dplyr::rename(SubjectID = plasma_vcode...21,BQCID=bqc_id)
  all.soma = inner_join(all.soma,vcode.bqcid.dict,by='SubjectID') 
  names(all.soma)[[2]] = 'Analyte'
  
  return(all.soma)
}
getRelevantRedcap = function() { #get patient information
  sample.genotype.dictionary = vroom(file='~/pt_age_sex.tsv') %>%
    dplyr::select(Alias.BQCid,genotypeID) %>% dplyr::rename(BQCID=Alias.BQCid)
  redcap.data = vroom(file='~/redcap_clinical_data_raw_2022-05-12.csv') %>%
    dplyr::select(BQCID,copy_female,copy_age) %>% distinct(BQCID,.keep_all=T)
  common.data = inner_join(redcap.data,sample.genotype.dictionary,by='BQCID')
  outcome.data = readRDS('~/BQC19_A2B2C2.rds') %>%
    dplyr::rename(BQCID = BQC.identifier..public.) %>% dplyr::select(BQCID,A2,B2,C2)
  common.data = inner_join(common.data,outcome.data,by='BQCID')
  pc.data = vroom('~/EUR.pc') %>%
    dplyr::rename(genotypeID = FID) %>% dplyr::select(-IID)
  common.data = inner_join(common.data,pc.data,by='genotypeID') %>%
    dplyr::rename(Female=copy_female,Age=copy_age)
  testing.site = vroom('~/center_mapping_BQCID_20220126.csv')
  for (item in 1:nrow(testing.site)) {
    if (str_detect(testing.site$individual_id[[item]],'JGH')) testing.site$individual_id[[item]] = 'JGH'
    else if (str_detect(testing.site$individual_id[[item]],'CRCHUM')) testing.site$individual_id[[item]] = 'CRCHUM'
    else if (str_detect(testing.site$individual_id[[item]],'SLSJ')) testing.site$individual_id[[item]] = 'SLSJ'
    else if (str_detect(testing.site$individual_id[[item]],'MUHC')) testing.site$individual_id[[item]] = 'MUHC'
    else if (str_detect(testing.site$individual_id[[item]],'HSCM')) testing.site$individual_id[[item]] = 'HSCM'
    else if (str_detect(testing.site$individual_id[[item]],'CHUS')) testing.site$individual_id[[item]] = 'CHUS'
    else if (str_detect(testing.site$individual_id[[item]],'CHUQ')) testing.site$individual_id[[item]] = 'CHUQ'
  }
  testing.site %<>% dplyr::select(BQCID,individual_id) %>% dplyr::rename(Site=individual_id)
  common.data = inner_join(common.data,testing.site,by='BQCID')
  return(common.data)
}
doLogisticAnalysis = function(outcome,analyte) {
  soma = getSomaData(analyte)
  redcap = getRelevantRedcap()
  merged.data = inner_join(soma,redcap,by='BQCID')
  
  base.form = "Age + Female + Analyte + PlateId + Site + PC1 + PC2 + PC3 + 
                 PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"
  if (outcome == 'A2') {
    modl = glm(as.formula(glue('A2 ~ {base.form}')),family=binomial,
               data=merged.data)
  } else if (outcome == 'B2') {
    modl = glm(as.formula(glue('B2 ~ {base.form}')),family=binomial,
               data=merged.data)
  } else if (outcome == 'C2')
    modl = glm(as.formula(glue('C2 ~ {base.form}')),family=binomial,
               data=merged.data)
  #need centre of recruitment as covariate
  return(modl)
}
plotOddsRatios = function(modls,proteins,keep.indices) {
  df = data.frame(Product=character(),Beta=numeric(),SE=numeric())
  for (modl in 1:length(modls)) {
    if (modl %notin% keep.indices) next
    curr.modl = modls[[modl]]

    df %<>% add_row(Product=proteins[[modl]],Beta=summary(curr.modl)$coefficients[4,1],
                    SE=summary(curr.modl)$coefficients[4,2]) %>% 
      mutate(Significant=summary(curr.modl)$coefficients[4,4] <= 0.05/length(proteins))
  }
  df %<>% mutate(Product = str_replace(Product,'B2','Hospitalized'))
  df %<>% mutate(Product = str_replace(Product,'C2','Susceptibility'))
  print(ggplot(df,aes(y=Product,x=exp(Beta),xmin=exp(Beta-1.96*SE),xmax=exp(Beta+1.96*SE),color=Significant)) +
          geom_point() + geom_errorbar(width=0.2) + xlab('Odds ratio') +
          scale_y_discrete(limits=rev) + theme_bw() + theme(text=element_text(size=24),legend.position='None') +
          geom_vline(xintercept = 1.00,linetype='dotted'))
  return(df)
}
