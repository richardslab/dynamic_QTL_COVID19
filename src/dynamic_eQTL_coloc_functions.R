getPQTLData = function(PQTL,chr,lead.snp.position,curr.outcome,MR) {
  print('Reading pQTL data')
  # pqtl = vroom(glue('pQTL_Data/COVID19_HGI_{curr.outcome}_ALL_eur_leave23andme_20220403.tsv.gz'),show_col_types = F)
  names(PQTL)[1] = 'CHR'
  # print('Read in data, now subsetting/cleaning up data')
  
  pqtl = PQTL %>% subset(CHR == chr) %>% mutate(varbeta=all_inv_var_meta_sebeta^2) %>% 
    dplyr::filter(!duplicated(SNP),!is.na(varbeta)) %>%
    dplyr::rename(position=POS,beta=all_inv_var_meta_beta,MAF=all_meta_AF,
                  pvalues=all_inv_var_meta_p,se=all_inv_var_meta_sebeta) %>%
    mutate(variant_id=paste0('chr',SNP))
  pqtl %<>% mutate(N=all_inv_var_meta_cases + all_inv_var_meta_controls)
  if (!MR) pqtl %<>% mutate(rsid=NULL)
  
  if (!is.null(lead.snp.position)) pqtl = pqtl %>% 
    dplyr::filter(position >= lead.snp.position-500000,position<=lead.snp.position+500000)
  
  return(pqtl %>% dplyr::filter(!is.na(MAF)))
}
getCommonQTLLists = function(mqtl,pqtl,infectious,curr.outcome) {
  common.snps = inner_join(mqtl,pqtl,by='variant_id') %>% dplyr::select(variant_id)
  mqtl.common = inner_join(mqtl,common.snps,by='variant_id') %>% arrange(position,pvalues)
  pqtl.common = inner_join(pqtl,common.snps,by='variant_id') %>% arrange(position,pvalues)

  mqtl.list = as.list(mqtl.common)
  mqtl.list$snp=mqtl.list$variant_id
  mqtl.list$type = 'quant'
  
  pqtl.list = as.list(pqtl.common)
  pqtl.list$type = 'cc'
  pqtl.list$snp = pqtl.list$variant_id
  pqtl.list$s = pqtl.list$all_inv_var_meta_cases/(pqtl.list$all_inv_var_meta_cases + pqtl.list$all_inv_var_meta_controls)

  return(list(mqtl.list,pqtl.list))
}
makeTxtAllSoskicParquets = function() {
  files = list.files('Soskic_eQTL_summ_stats',pattern='.parquet',full.names=F)
  files = str_replace(files,'.parquet','')
  for (item in files) {
    print(glue('On item {which(files == item)} of {length(files)}'))
    curr.file = read_parquet(glue('Soskic_eQTL_summ_stats/{item}.parquet'))
    vroom_write(curr.file,glue('Soskic_eQTL_summ_stats/{item}.txt'),col_names=T,
                quote="none")
  }
}
soskicColoc = function(vector.specific,rel6,allCells,returnCommonQTLs,mr) {
  df = data.frame(Cell=character(),Outcome=character(),CHR=numeric(),POS=numeric(),
                  H0.PP=numeric(),
                  H1.PP=numeric(),H2.PP=numeric(),H3.PP=numeric(),
                  H4.PP=numeric(),LeadSNP_eQTL=character(),LeadSNP_HGI=character(),Gene=character(),
                  Lead_eQTL_Beta=numeric(),Lead_eQTL_P=numeric(),
                  Lead_HGI_Beta=numeric(),Lead_HGI_P=numeric(),
                  HGI_MAF=numeric())
  if (rel6) loci = read.table('subpositions_chrpos_rel6.txt')
  else loci = read.table('subpositions_chrpos.txt') #rel7
  
  names(loci) = c('chr','pos','outcome')
  print('Reading pQTL gwas')
  
  if (is.null(vector.specific)) {
    if (!rel6) {
      pqtl.a2 = vroom(glue('pQTL_Data/COVID19_HGI_A2_ALL_eur_leave23andme_20220403.tsv.gz'),show_col_types = F)
      pqtl.b2 = vroom(glue('pQTL_Data/COVID19_HGI_B2_ALL_eur_leave23andme_20220403.tsv.gz'),show_col_types = F)
      pqtl.c2 = vroom(glue('pQTL_Data/COVID19_HGI_C2_ALL_eur_leave23andme_20220403.tsv.gz'),show_col_types = F)
    }else{
      pqtl.b2 = vroom('/scratch/richards/julian.willett/eQTLpQTLDisconnect/pQTL_Data/COVID19_HGI_B2_ALL_eur_leave_ukbb_23andme_20210622.txt.gz',show_col_types = F)
      pqtl.c2 = vroom('/scratch/richards/julian.willett/eQTLpQTLDisconnect/pQTL_Data/COVID19_HGI_C2_ALL_eur_leave_23andme_20210622.txt.gz',show_col_types = F)
    }
  }else{
    if (vector.specific[[1]] == 'A2') pqtl.a2 = vroom(glue('pQTL_Data/COVID19_HGI_A2_ALL_eur_leave23andme_20220403.tsv.gz'),show_col_types = F)
    else if (vector.specific[[1]] == 'B2') pqtl.b2 = vroom(glue('pQTL_Data/COVID19_HGI_B2_ALL_eur_leave23andme_20220403.tsv.gz'),show_col_types = F)
    else if (vector.specific[[1]] == 'C2') pqtl.c2 = vroom(glue('pQTL_Data/COVID19_HGI_C2_ALL_eur_leave23andme_20220403.tsv.gz'),show_col_types = F)
  }
  for (item in 1:nrow(loci)) {
    if (is.null(vector.specific)) print(glue('On row {item} of {nrow(loci)}'))
    curr.chr = loci[['chr']][[item]] ; targ.pos = loci[['pos']][[item]]
    outcome = loci[['outcome']][[item]]
    if (outcome == 'A2' & rel6) next
    if (!is.null(vector.specific)) {
      if (outcome != vector.specific[[1]]) next
      else if (curr.chr != vector.specific[[2]] | targ.pos != vector.specific[[3]]) next
    }
    print(glue('Investigating outcome {outcome} chr{curr.chr} pos{targ.pos}'))
    
    if (outcome == 'A2') pQTL = getPQTLData(pqtl.a2,chr=curr.chr,lead.snp.position=targ.pos,outcome,MR=mr) #grch38
    else if (outcome == 'B2') pQTL = getPQTLData(pqtl.b2,chr=curr.chr,lead.snp.position=targ.pos,outcome,MR=mr) #grch38
    else if (outcome == 'C2') pQTL = getPQTLData(pqtl.c2,chr=curr.chr,lead.snp.position=targ.pos,outcome,MR=mr) #grch38
    
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
        if (!is.null(vector.specific)) {
          if (vector.specific[[6]] %in% all.genes) print('Gene in the given mQTL data')
          else print('Gene not in mQTL data, ie not expressed')
        }
        for (gene in all.genes) {
          if (!is.null(vector.specific)) {
            if (gene != vector.specific[[6]]) { next }
          }
          mQTL = full.mQTL %>% dplyr::filter(phenotype_id == gene) %>% 
            mutate(N=round(ma_count/maf*(1/2)))
          if (nrow(mQTL) == 0) { print(glue('No rows {gene}')) ; next}
          mQTL %<>% dplyr::rename(beta=slope,se=slope_se,pvalues=pval_nominal,MAF=maf,
                                  chr=CHR,position=POS) %>% mutate(varbeta=se^2) %>%
            dplyr::select(chr,position,REF,ALT,beta,varbeta,pvalues,MAF,variant_id,phenotype_id,ma_samples,N,se) %>%
            dplyr::filter(!is.na(variant_id),!duplicated(variant_id),!is.na(varbeta),
                          pvalues < 1) #coloc freaks out if pvalue = 1
          mQTL %<>% mutate(variant_id = glue('chr{mQTL$chr}:{mQTL$position}:{mQTL$REF}:{mQTL$ALT}'))
          if (nrow(mQTL) == 0) { print(glue('No rows {gene}')) ; next}
          
          common.qtls.list = getCommonQTLLists(mQTL,pQTL,group.infectious,outcome)
          if (length(common.qtls.list[[1]][['pvalues']]) == 0) { print(glue('No overlapping rows {file}')) ; break }
          if (length(common.qtls.list[[1]][['pvalues']]) < 50) {
            if (!is.null(vector.specific)) print(glue('Not enough variant On locus {item} of {nrow(loci)}, cell {cell}'))
            next #follow criteria per Soskic
          }

          if (returnCommonQTLs & !is.null(vector.specific)) return(list(mQTL,pQTL))
          coloc.res = coloc.abf(common.qtls.list[[1]],common.qtls.list[[2]])
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
          df %<>% add_row(Cell=file,Outcome=outcome,CHR=curr.chr,POS=targ.pos,
                          H0.PP=coloc.res['PP.H0.abf'],H1.PP=coloc.res['PP.H1.abf'],
                          H2.PP=coloc.res['PP.H2.abf'],H3.PP=coloc.res['PP.H3.abf'],
                          H4.PP=coloc.res['PP.H4.abf'],LeadSNP_eQTL=lead.snp.eqtl,LeadSNP_HGI=lead.snp.hgi,
                          Gene=gene,
                          Lead_eQTL_Beta=eQTL.beta,Lead_eQTL_P=eQTL.P,
                          Lead_HGI_Beta=HGI.beta,Lead_HGI_P=HGI.P,
                          HGI_MAF=HGI.maf)
          if (!is.null(vector.specific)) print(glue('On locus {item} of {nrow(loci)}, cell {cell}'))
        }
      }
    }
  }
  return(df %>% add_column('H3+H4'=df$H3.PP+df$H4.PP))
}
bqcColoc = function(type.of.qtl,vector.specific,rel6,returnCommonQTLs,mr) {
  df = data.frame(InfState=character(),Outcome=character(),CHR=numeric(),POS=numeric(),
                  H0.PP=numeric(),
                  H1.PP=numeric(),H2.PP=numeric(),H3.PP=numeric(),
                  H4.PP=numeric(),LeadSNP_eQTL=character(),LeadSNP_HGI=character(),Gene=character(),
                  Lead_eQTL_Beta=numeric(),Lead_eQTL_P=numeric(),
                  Lead_HGI_Beta=numeric(),Lead_HGI_P=numeric(),
                  HGI_MAF=numeric())
  if (rel6) loci = read.table('subpositions_chrpos_rel6.txt')
  else loci = read.table('subpositions_chrpos.txt') #rel7
  
  names(loci) = c('chr','pos','outcome')
  print('Reading pQTL gwas')
  
  if (is.null(vector.specific)) {
    if (!rel6) {
      pqtl.a2 = vroom(glue('pQTL_Data/COVID19_HGI_A2_ALL_eur_leave23andme_20220403.tsv.gz'),show_col_types = F)
      pqtl.b2 = vroom(glue('pQTL_Data/COVID19_HGI_B2_ALL_eur_leave23andme_20220403.tsv.gz'),show_col_types = F)
      pqtl.c2 = vroom(glue('pQTL_Data/COVID19_HGI_C2_ALL_eur_leave23andme_20220403.tsv.gz'),show_col_types = F)
    }else{
      pqtl.b2 = vroom('/scratch/richards/julian.willett/eQTLpQTLDisconnect/pQTL_Data/COVID19_HGI_B2_ALL_eur_leave_ukbb_23andme_20210622.txt.gz',show_col_types = F)
      pqtl.c2 = vroom('/scratch/richards/julian.willett/eQTLpQTLDisconnect/pQTL_Data/COVID19_HGI_C2_ALL_eur_leave_23andme_20210622.txt.gz',show_col_types = F)
    }
  }else{
    if (vector.specific[[1]] == 'A2') pqtl.a2 = vroom(glue('pQTL_Data/COVID19_HGI_A2_ALL_eur_leave23andme_20220403.tsv.gz'),show_col_types = F)
    else if (vector.specific[[1]] == 'B2') pqtl.b2 = vroom(glue('pQTL_Data/COVID19_HGI_B2_ALL_eur_leave23andme_20220403.tsv.gz'),show_col_types = F)
    else if (vector.specific[[1]] == 'C2') pqtl.c2 = vroom(glue('pQTL_Data/COVID19_HGI_C2_ALL_eur_leave23andme_20220403.tsv.gz'),show_col_types = F)
  }
  
  for (item in 1:nrow(loci)) {
    if (is.null(vector.specific)) print(glue('On row {item} of {nrow(loci)}'))
    curr.chr = loci[['chr']][[item]] ; targ.pos = loci[['pos']][[item]]
    outcome = loci[['outcome']][[item]]
    if (outcome == 'A2' & rel6) next
    if (!is.null(vector.specific)) {
      if (outcome != vector.specific[[1]]) next
      else if (curr.chr != vector.specific[[2]] | targ.pos != vector.specific[[3]]) next
    }
    print(glue('Investigating outcome {outcome} chr{curr.chr} pos{targ.pos}'))
    
    if (outcome == 'A2') pQTL = getPQTLData(pqtl.a2,chr=curr.chr,lead.snp.position=targ.pos,outcome,MR=mr) #grch38
    else if (outcome == 'B2') pQTL = getPQTLData(pqtl.b2,chr=curr.chr,lead.snp.position=targ.pos,outcome,MR=mr) #grch38
    else if (outcome == 'C2') pQTL = getPQTLData(pqtl.c2,chr=curr.chr,lead.snp.position=targ.pos,outcome,MR=mr) #grch38
    
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
        # print(paste('On: ',gene,item,state))
        if (!is.null(vector.specific)) {
          if (!str_detect(gene,vector.specific[[5]])) next #given bqc using transcript ids as decimal
        }
        mQTL = full.mQTL %>% dplyr::filter(gene_id == gene) %>% 
          mutate(N=round(ma_count/maf*(1/2)))
        if (nrow(mQTL) == 0) { print(glue('No rows {gene}')) ; next}
        mQTL %<>% dplyr::rename(beta=slope,se=slope_se,pvalues=pval_nominal,MAF=maf,
                                chr=CHR,position=POS,phenotype_id=gene_id) %>% mutate(varbeta=se^2) %>%
          dplyr::select(chr,position,REF,ALT,beta,varbeta,pvalues,MAF,variant_id,phenotype_id,ma_samples,N,se) %>%
          dplyr::filter(!is.na(variant_id),!duplicated(variant_id),!is.na(varbeta),
                        pvalues < 1) #coloc freaks out if pvalue = 1
        mQTL %<>% mutate(variant_id = glue('{mQTL$chr}:{mQTL$position}:{mQTL$REF}:{mQTL$ALT}'))
        if (nrow(mQTL) == 0) { print(glue('No rows {gene}')) ; next}
        
        common.qtls.list = getCommonQTLLists(mQTL,pQTL,group.infectious,outcome)
        if (returnCommonQTLs) return(common.qtls.list)
        
        if (length(common.qtls.list[[1]][['pvalues']]) == 0) { print(glue('No overlapping rows')) ; break }
        if (length(common.qtls.list[[1]][['pvalues']]) < 50) {
          print(glue('Not enough variant On locus {item} of {nrow(loci)}, infstate {state}'))
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
        df %<>% add_row(InfState=state,Outcome=outcome,CHR=curr.chr,POS=targ.pos,
                        H0.PP=coloc.res['PP.H0.abf'],H1.PP=coloc.res['PP.H1.abf'],
                        H2.PP=coloc.res['PP.H2.abf'],H3.PP=coloc.res['PP.H3.abf'],
                        H4.PP=coloc.res['PP.H4.abf'],LeadSNP_eQTL=lead.snp.eqtl,LeadSNP_HGI=lead.snp.hgi,
                        Gene=gene,
                        Lead_eQTL_Beta=eQTL.beta,Lead_eQTL_P=eQTL.P,
                        Lead_HGI_Beta=HGI.beta,Lead_HGI_P=HGI.P,
                        HGI_MAF=HGI.maf)
        print(glue('On locus {item} of {nrow(loci)}, inf.state {state}'))
      }
    }
  }
  return(df %>% add_column('H3+H4'=df$H3.PP+df$H4.PP))
}
gtexColoc = function(type.of.qtl,vector.specific,returnCommonQTLs,mr,onlySoskicLoci) {
  df = data.frame(Outcome=character(),CHR=numeric(),POS=numeric(),
                  H0.PP=numeric(),
                  H1.PP=numeric(),H2.PP=numeric(),H3.PP=numeric(),
                  H4.PP=numeric(),LeadSNP_eQTL=character(),LeadSNP_HGI=character(),Gene=character(),
                  Lead_eQTL_Beta=numeric(),Lead_eQTL_P=numeric(),
                  Lead_HGI_Beta=numeric(),Lead_HGI_P=numeric(),
                  HGI_MAF=numeric())
  loci = read.table('subpositions_chrpos.txt') %>% dplyr::arrange(V3) #rel7
  
  names(loci) = c('chr','pos','outcome')
  print('Reading pQTL gwas')
  
  if (is.null(vector.specific)) {
    pqtl.a2 = vroom(glue('pQTL_Data/COVID19_HGI_A2_ALL_eur_leave23andme_20220403.tsv.gz'),show_col_types = F)
    pqtl.b2 = vroom(glue('pQTL_Data/COVID19_HGI_B2_ALL_eur_leave23andme_20220403.tsv.gz'),show_col_types = F)
    pqtl.c2 = vroom(glue('pQTL_Data/COVID19_HGI_C2_ALL_eur_leave23andme_20220403.tsv.gz'),show_col_types = F)
  }else{
    if (vector.specific[[1]] == 'A2') pqtl.a2 = vroom(glue('pQTL_Data/COVID19_HGI_A2_ALL_eur_leave23andme_20220403.tsv.gz'),show_col_types = F)
    else if (vector.specific[[1]] == 'B2') pqtl.b2 = vroom(glue('pQTL_Data/COVID19_HGI_B2_ALL_eur_leave23andme_20220403.tsv.gz'),show_col_types = F)
    else if (vector.specific[[1]] == 'C2') pqtl.c2 = vroom(glue('pQTL_Data/COVID19_HGI_C2_ALL_eur_leave23andme_20220403.tsv.gz'),show_col_types = F)
  }
  
  for (item in 1:nrow(loci)) {
    if (is.null(vector.specific)) print(glue('On row {item} of {nrow(loci)}'))
    if (onlySoskicLoci & !(item %in% c(23,35,28,19,52,57,70,54,80,86,74,79))) next
    curr.chr = loci[['chr']][[item]] ; targ.pos = loci[['pos']][[item]]
    outcome = loci[['outcome']][[item]]
    if (!is.null(vector.specific)) {
      if (outcome != vector.specific[[1]]) next
      else if (curr.chr != vector.specific[[2]] | targ.pos != vector.specific[[3]]) next
    }
    print(glue('Investigating outcome {outcome} chr{curr.chr} pos{targ.pos}'))
    
    if (outcome == 'A2') pQTL = getPQTLData(pqtl.a2,chr=curr.chr,lead.snp.position=targ.pos,outcome,MR=mr) #grch38
    else if (outcome == 'B2') pQTL = getPQTLData(pqtl.b2,chr=curr.chr,lead.snp.position=targ.pos,outcome,MR=mr) #grch38
    else if (outcome == 'C2') pQTL = getPQTLData(pqtl.c2,chr=curr.chr,lead.snp.position=targ.pos,outcome,MR=mr) #grch38
    
    if (type.of.qtl == 'eQTL') full.mQTL = vroom(glue('GTEx_eQTL_Data/Subpositions/chr{curr.chr}_{targ.pos}_window.txt'),show_col_types=F)
    else if (type.of.qtl == 'sQTL') full.mQTL = vroom(glue('GTEx_eQTL_Data/Subpositions/chr{curr.chr}_{targ.pos}_window.txt'),show_col_types=F)
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
      
      common.qtls.list = getCommonQTLLists(mQTL,pQTL,F,outcome)
      if (returnCommonQTLs) return(common.qtls.list)
      
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
plotDynamicPP = function(data,cell,outcome,chromosomes,exclude.genes) {
  plot.data = data %>% dplyr::filter(str_detect(data$Cell,cell),Outcome==outcome,CHR %in% chromosomes)
  #ensure focused on exact sample
  if (cell == 'TEM_') plot.data %<>% dplyr::filter(!str_detect(plot.data$Cell,'TEM_HLA'),
                                                  !str_detect(plot.data$Cell,'TEM_LA'))
  if (cell == 'TCM_') plot.data %<>% dplyr::filter(!str_detect(plot.data$Cell,'TCM_LA'))
  if (cell == 'TN_') plot.data %<>% dplyr::filter(!str_detect(plot.data$Cell,'TN_IFN'),
                                                   !str_detect(plot.data$Cell,'TN_cycling'),
                                                  !str_detect(plot.data$Cell,'TN_HSP'),
                                                  !str_detect(plot.data$Cell,'TN_LA'),
                                                  !str_detect(plot.data$Cell,'TN_NFKB'))
  colocalizing.genes = unique(plot.data[['Gene']][which(plot.data$H4.PP >= 0.8)])
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
  if (outcome == 'A2') text.outcome = 'Very severe'
  else if (outcome == 'B2') text.outcome = 'Severe'
  else if (outcome == 'C2') text.outcome = 'Susceptible'
  
  plot.df$Gene = str_replace(plot.df$Gene,'ENSG00000104388','RAB2A')
  plot.df$Gene = str_replace(plot.df$Gene,'ENSG00000121691','CAT')
  plot.df$Gene = str_replace(plot.df$Gene,'ENSG00000130816','DNMT1')
  plot.df$Gene = str_replace(plot.df$Gene,'ENSG00000131400','NAPSA')
  plot.df$Gene = str_replace(plot.df$Gene,'ENSG00000144802','NFKBIZ')
  plot.df$Gene = str_replace(plot.df$Gene,'ENSG00000166997','CNPY4')
  plot.df$Gene = str_replace(plot.df$Gene,'ENSG00000176715','ACSF3')
  plot.df$Gene = str_replace(plot.df$Gene,'ENSG00000213420','GPC2')
  plot.df$Gene = str_replace(plot.df$Gene,'ENSG00000232300','FAM215B')
  plot.df$Gene = str_replace(plot.df$Gene,'ENSG00000204472','AIF1')
  plot.df$Gene = str_replace(plot.df$Gene,'ENSG00000087076','HSD17B14')
  plot.df$Gene = str_replace(plot.df$Gene,'ENSG00000206503','HLA-A')
  plot.df$Gene = str_replace(plot.df$Gene,'ENSG00000160271','RALGDS')
  plot.df$Gene = str_replace(plot.df$Gene,'ENSG00000143537','ADAM15')
  plot.df$Gene = str_replace(plot.df$Gene,'ENSG00000177628','GBA')
  
  plot.df %<>% dplyr::filter(!(Gene %in% exclude.genes))
  
  print(ggplot(plot.df,aes(x=Time,y=PP,group=Gene)) + geom_line(aes(color=Gene)) +
          geom_point(aes(shape=Gene),size=4) + theme(text=element_text(size=20)) + geom_hline(yintercept=0.8) +
          ggtitle(glue('Cell: {cell.title}. Outcome: {text.outcome}')) +
          scale_x_discrete(limits=c('0h','16h','40h','5d')) + ylim(c(0,1)))
  return(plot.df)
}
plotDynamicPP.BQC = function(data,outcome,chromosomes,exclude.genes) {
  plot.data = data %>% dplyr::filter(Outcome==outcome,CHR %in% chromosomes)

  colocalizing.genes = unique(plot.data[['Gene']][which(plot.data$H4.PP >= 0.8)])
  plot.data %<>% dplyr::filter(Gene %in% colocalizing.genes)
  plot.df = data.frame(State=character(),PP=numeric(),Gene=character(),Position=numeric())

  plot.df %<>% add_row(State = plot.data$InfState,PP=plot.data$H4.PP,Gene=plot.data$Gene,
                       Position=plot.data$POS)
  
  if (outcome == 'A2') text.outcome = 'Very severe'
  else if (outcome == 'B2') text.outcome = 'Severe'
  else if (outcome == 'C2') text.outcome = 'Susceptible'
  
  plot.df$Gene = str_replace(plot.df$Gene,'ENSG00000104388.15','RAB2A')
  plot.df$Gene = str_replace(plot.df$Gene,'ENSG00000121691','CAT')
  plot.df$Gene = str_replace(plot.df$Gene,'ENSG00000130816','DNMT1')
  plot.df$Gene = str_replace(plot.df$Gene,'ENSG00000131400.8','NAPSA')
  plot.df$Gene = str_replace(plot.df$Gene,'ENSG00000144802','NFKBIZ')
  plot.df$Gene = str_replace(plot.df$Gene,'ENSG00000166997','CNPY4')
  plot.df$Gene = str_replace(plot.df$Gene,'ENSG00000176715','ACSF3')
  plot.df$Gene = str_replace(plot.df$Gene,'ENSG00000213420','GPC2')
  plot.df$Gene = str_replace(plot.df$Gene,'ENSG00000232300','FAM215B')
  plot.df$Gene = str_replace(plot.df$Gene,'ENSG00000204472','AIF1')
  plot.df$Gene = str_replace(plot.df$Gene,'ENSG00000087076','HSD17B14')
  plot.df$Gene = str_replace(plot.df$Gene,'ENSG00000206503','HLA-A')
  plot.df$Gene = str_replace(plot.df$Gene,'ENSG00000160271','RALGDS')
  plot.df$Gene = str_replace(plot.df$Gene,'ENSG00000143537','ADAM15')
  plot.df$Gene = str_replace(plot.df$Gene,'ENSG00000177628','GBA')
  plot.df$Gene = str_replace(plot.df$Gene,'ENSG00000175164.16','ABO')
  plot.df$Gene = str_replace(plot.df$Gene,'ENSG00000142233.14','NTN5')
  plot.df$Gene = str_replace(plot.df$Gene,'ENSG00000169231.13','THBS3')
  plot.df$Gene = str_replace(plot.df$Gene,'ENSG00000160766.14','GBAP1')
  plot.df$Gene = str_replace(plot.df$Gene,'ENSG00000105397.14','TYK2')
  plot.df$Gene = str_replace(plot.df$Gene,'ENSG00000131401.11','NAPSB')
   
  plot.df %<>% dplyr::filter(!(Gene %in% exclude.genes))
  
  print(ggplot(plot.df,aes(x=State,y=PP,group=Gene)) + 
          geom_bar(aes(fill=Gene),stat='identity',position=position_dodge()) +
          theme(text=element_text(size=20)) + geom_hline(yintercept=0.8) +
          ggtitle(glue('Outcome: {text.outcome}')) +
          ylim(c(0,1)))
  return(plot.df)
}
conductMRBoskic = function(outcome,variant_id,cell,time,gene) {
  split.id = str_split(variant_id,':')[[1]]
  curr.chr = split.id[[1]] ; targ.pos = split.id[[2]]
  common.qtls = soskicColoc(vector.specific=c(outcome,curr.chr,targ.pos,cell,time,gene),rel6=F,allCells=T,returnCommonQTLs = T,mr=T)
  rsids = common.qtls[[2]][['rsid']]

  #prep exposure data. 
  exp.rsid.df = data.frame(VID=common.qtls[[1]][['variant_id']],Rsid=NA)
  for (item in 1:nrow(exp.rsid.df)) {
    if (exp.rsid.df$VID[[item]] %in% common.qtls[[2]][['variant_id']])
      exp.rsid.df$Rsid[[item]] = rsids[which(common.qtls[[2]][['variant_id']]==exp.rsid.df$VID[[item]])]
  }
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
  
  #harmonize
  harmonized = harmonise_data(exposure.data,outcome.data)

  #do mr functions
  to.return = data.frame()
  if (sum(harmonized$mr_keep)==0) print('No SNP remaining after harmonization')
  else if (sum(harmonized$mr_keep)==1) {
    result = mr(harmonized, method_list = "mr_wald_ratio")
    to.return = result[1,]
    cat("eBMD | Wald-Ratio | MR nsnp",result$nsnp,"| b:",result$b,"| se",result$se,"| p-value",result$pval, collapse="\n")
  }else if (sum(harmonized$mr_keep)>1) {
    result <- mr(harmonized, method_list = c("mr_ivw","mr_two_sample_ml","mr_weighted_median","mr_penalised_weighted_median","mr_weighted_mode"))
    to.return = result[1,]
    cat("IVW | MR nsnp",result$nsnp[1],"| b:",result$b[1],"| se",result$se[1],"| p-value",result$pval[1], collapse="\n")
    if (sum(harmonized$mr_keep)>2) {
      cat("WMedian | MR nsnp",result$nsnp[3],"| b:",result$b[3],"| se",result$se[3],"| p-value",result$pval[3], collapse="\n")
      cat("Wmode | MR nsnp",result$nsnp[5],"| b:",result$b[5],"| se",result$se[5],"| p-value",result$pval[5], collapse="\n")
      egger_result <- mr_egger_regression(harmonized$beta.exposure[harmonized$mr_keep],
                                          harmonized$beta.outcome[harmonized$mr_keep],
                                          harmonized$se.exposure[harmonized$mr_keep],
                                          harmonized$se.outcome[harmonized$mr_keep])
      cat("Egger | MR nsnp",egger_result$nsnp,"| b:",egger_result$b,"| se",egger_result$se,"| p-value",egger_result$pval, collapse="\n")
      cat("Egger-intercept | MR nsnp",egger_result$nsnp,"| b:",egger_result$b_i,"| se",egger_result$se_i,"| p-value",egger_result$pval_i, collapse="\n")
    }
  }
  
  #do mr
  mr.results = directionality_test(harmonized)
  cat("Directionality | p-value",mr.results$steiger_pval, collapse="\n")
  return(to.return %>% mutate(Gene=gene,Time=time,Cell=cell,Outcome=outcome))  
}
conductMRBQC = function(qtl.type,outcome,variant_id,state,gene) {
  split.id = str_split(variant_id,':')[[1]]
  curr.chr = split.id[[1]] ; targ.pos = split.id[[2]]
  common.qtls = bqcColoc(qtl.type,vector.specific=c(outcome,curr.chr,targ.pos,state,gene),rel6=F,returnCommonQTLs = T,mr=T)
  rsids = common.qtls[[2]][['rsid']]

  #prep exposure data. 
  exposure.data = data.frame(chr=common.qtls[[1]][['chr']],position=common.qtls[[1]][['position']],
                             beta=common.qtls[[1]][['beta']],se=common.qtls[[1]][['se']],
                             SNP=rsids,effect_allele=common.qtls[[1]][['ALT']],
                             other_allele=common.qtls[[1]][['REF']],eaf=common.qtls[[1]][['MAF']],
                             samplesize=common.qtls[[1]][['N']],pval=common.qtls[[1]][['pvalues']])
  exposure.data = format_data(exposure.data,type='exposure')
  exposure.data = clump_data(exposure.data,pop='EUR')
  
  #prep outcome data
  outcome.data = data.frame(chr=common.qtls[[2]][['CHR']],position=common.qtls[[2]][['position']],
                            beta=common.qtls[[2]][['beta']],se=common.qtls[[2]][['se']],
                            SNP=rsids,effect_allele=common.qtls[[2]][['ALT']],
                            other_allele=common.qtls[[2]][['REF']],eaf=common.qtls[[2]][['MAF']],
                            samplesize=common.qtls[[2]][['N']],pval=common.qtls[[2]][['pvalues']])
  outcome.data = format_data(outcome.data,type='outcome')
  
  #harmonize
  harmonized = harmonise_data(exposure.data,outcome.data)
  
  #do mr functions
  to.return=data.frame()
  if (sum(harmonized$mr_keep)==0) print('No SNP remaining after harmonization')
  else if (sum(harmonized$mr_keep)==1) {
    result = mr(harmonized, method_list = "mr_wald_ratio")
    to.return = result[1,]
    cat("eBMD | Wald-Ratio | MR nsnp",result$nsnp,"| b:",result$b,"| se",result$se,"| p-value",result$pval, collapse="\n")
  }else if (sum(harmonized$mr_keep)>1) {
    result <- mr(harmonized, method_list = c("mr_ivw","mr_two_sample_ml","mr_weighted_median","mr_penalised_weighted_median","mr_weighted_mode"))
    to.return = result[1,]
    cat("IVW | MR nsnp",result$nsnp[1],"| b:",result$b[1],"| se",result$se[1],"| p-value",result$pval[1], collapse="\n")
    if (sum(harmonized$mr_keep)>2) {
      cat("WMedian | MR nsnp",result$nsnp[3],"| b:",result$b[3],"| se",result$se[3],"| p-value",result$pval[3], collapse="\n")
      cat("Wmode | MR nsnp",result$nsnp[5],"| b:",result$b[5],"| se",result$se[5],"| p-value",result$pval[5], collapse="\n")
      egger_result <- mr_egger_regression(harmonized$beta.exposure[harmonized$mr_keep],
                                          harmonized$beta.outcome[harmonized$mr_keep],
                                          harmonized$se.exposure[harmonized$mr_keep],
                                          harmonized$se.outcome[harmonized$mr_keep])
      cat("Egger | MR nsnp",egger_result$nsnp,"| b:",egger_result$b,"| se",egger_result$se,"| p-value",egger_result$pval, collapse="\n")
      cat("Egger-intercept | MR nsnp",egger_result$nsnp,"| b:",egger_result$b_i,"| se",egger_result$se_i,"| p-value",egger_result$pval_i, collapse="\n")
    }
  }
  
  #do mr
  mr.results = directionality_test(harmonized)
  cat("Directionality | p-value",mr.results$steiger_pval, collapse="\n")
  return(to.return %>% mutate(State=state,Outcome=outcome,Gene=gene))  
}
conductMRGTEx = function(qtl.type,outcome,variant_id,state,gene) {
  split.id = str_split(variant_id,':')[[1]]
  curr.chr = split.id[[1]] ; targ.pos = split.id[[2]]
  common.qtls = gtexColoc(qtl.type,vector.specific=c(outcome,curr.chr,targ.pos,NA,gene),returnCommonQTLs = T,mr=T,onlySoskicLoci = T)
  rsids = common.qtls[[2]][['rsid']]
  
  #prep exposure data. 
  exposure.data = data.frame(chr=common.qtls[[1]][['chr']],position=common.qtls[[1]][['position']],
                             beta=common.qtls[[1]][['beta']],se=common.qtls[[1]][['se']],
                             SNP=rsids,effect_allele=common.qtls[[1]][['ALT']],
                             other_allele=common.qtls[[1]][['REF']],eaf=common.qtls[[1]][['MAF']],
                             samplesize=common.qtls[[1]][['N']],pval=common.qtls[[1]][['pvalues']])
  exposure.data = format_data(exposure.data,type='exposure')
  exposure.data = clump_data(exposure.data,pop='EUR')
  
  #prep outcome data
  outcome.data = data.frame(chr=common.qtls[[2]][['CHR']],position=common.qtls[[2]][['position']],
                            beta=common.qtls[[2]][['beta']],se=common.qtls[[2]][['se']],
                            SNP=rsids,effect_allele=common.qtls[[2]][['ALT']],
                            other_allele=common.qtls[[2]][['REF']],eaf=common.qtls[[2]][['MAF']],
                            samplesize=common.qtls[[2]][['N']],pval=common.qtls[[2]][['pvalues']])
  outcome.data = format_data(outcome.data,type='outcome')
  
  #harmonize
  harmonized = harmonise_data(exposure.data,outcome.data)
  
  #do mr functions
  to.return=data.frame()
  if (sum(harmonized$mr_keep)==0) print('No SNP remaining after harmonization')
  else if (sum(harmonized$mr_keep)==1) {
    result = mr(harmonized, method_list = "mr_wald_ratio")
    to.return = result[1,]
    cat("eBMD | Wald-Ratio | MR nsnp",result$nsnp,"| b:",result$b,"| se",result$se,"| p-value",result$pval, collapse="\n")
  }else if (sum(harmonized$mr_keep)>1) {
    result <- mr(harmonized, method_list = c("mr_ivw","mr_two_sample_ml","mr_weighted_median","mr_penalised_weighted_median","mr_weighted_mode"))
    to.return = result[1,]
    cat("IVW | MR nsnp",result$nsnp[1],"| b:",result$b[1],"| se",result$se[1],"| p-value",result$pval[1], collapse="\n")
    if (sum(harmonized$mr_keep)>2) {
      cat("WMedian | MR nsnp",result$nsnp[3],"| b:",result$b[3],"| se",result$se[3],"| p-value",result$pval[3], collapse="\n")
      cat("Wmode | MR nsnp",result$nsnp[5],"| b:",result$b[5],"| se",result$se[5],"| p-value",result$pval[5], collapse="\n")
      egger_result <- mr_egger_regression(harmonized$beta.exposure[harmonized$mr_keep],
                                          harmonized$beta.outcome[harmonized$mr_keep],
                                          harmonized$se.exposure[harmonized$mr_keep],
                                          harmonized$se.outcome[harmonized$mr_keep])
      cat("Egger | MR nsnp",egger_result$nsnp,"| b:",egger_result$b,"| se",egger_result$se,"| p-value",egger_result$pval, collapse="\n")
      cat("Egger-intercept | MR nsnp",egger_result$nsnp,"| b:",egger_result$b_i,"| se",egger_result$se_i,"| p-value",egger_result$pval_i, collapse="\n")
    }
  }
  
  #do mr
  mr.results = directionality_test(harmonized)
  cat("Directionality | p-value",mr.results$steiger_pval, collapse="\n")
  return(to.return %>% mutate(State=state,Outcome=outcome,Gene=gene))  
}
plotMRResults = function(dfs,solid.hits,outcome) { #take in list of IVW dfs
  df = dfs[[1]]
  if (length(dfs)>1) { for (item in 2:length(dfs)) df %<>% add_row(dfs[[item]]) }
  df$Gene = str_replace(df$Gene,'ENSG00000204472','AIF1')
  df$Gene = str_replace(df$Gene,'ENSG00000121691','CAT')
  df$Gene = str_replace(df$Gene,'ENSG00000176715','ACSF3')
  df$Gene = str_replace(df$Gene,'ENSG00000131400','NAPSA')
  df$Gene = str_replace(df$Gene,'ENSG00000104388','RAB2A')
  df$Gene = str_replace(df$Gene,'ENSG00000087076','HSD17B14')
  df$Gene = str_replace(df$Gene,'ENSG00000160271','RALGDS')
  df$Gene = str_replace(df$Gene,'ENSG00000143537','ADAM15')
  df$Cell = str_replace(df$Cell,'CD4_Naive','CD4 Naive')
  df$Cell = str_replace(df$Cell,'TN_0h','TN')
  df$Cell = str_replace(df$Cell,'TN_16h','TN')
  df$Cell = str_replace(df$Cell,'TN_40h','TN')
  df$Cell = str_replace(df$Cell,'CD4_Memory','CD4 Memory')
  df$Cell = str_replace(df$Cell,'TEM_16h','TEM')
  df$Cell = str_replace(df$Cell,'TEM_40h','TEM')
  df$Cell = str_replace(df$Cell,'TCM_','TCM')
  
  out=outcome
  df %<>% mutate(CellGeneTime = paste(df$Cell,df$Gene,str_replace_all(df$Time,"_",'')),Solid=solid.hits)

  plot = ggplot(df,aes(y=CellGeneTime,x=exp(b),xmin=exp(b-1.96*se),xmax=exp(b+1.96*se),color=Solid)) +
          geom_point() + geom_errorbar() + xlab('Odds ratio') +
          scale_y_discrete() + theme(text=element_text(size=20),legend.position="none") +
          geom_vline(xintercept=1) + ylab('')
  if (outcome == 'A2')
    plot = plot + scale_y_discrete(limits=rev(c('CD4 Naive CAT 0h','CD4 Naive NAPSA 40h',
                                            'TN CAT 0h','TN NAPSA 40h',
                                            'CD4 Memory ACSF3 5d','CD4 Memory CAT 5d',
                                            'TCM CAT 5d','TEM NAPSA 40h',
                                            'TEM RALGDS 16h')))
  else if (outcome == 'B2')
    plot = plot + scale_y_discrete(limits=rev(c('CD4 Naive CAT 0h','CD4 Naive NAPSA 40h',
                                                'CD4 Naive RAB2A 16h','CD4 Naive RAB2A 40h',
                                                'TN CAT 0h','TN NAPSA 40h','TN RAB2A 16h',
                                                'TN RAB2A 40h','CD4 Memory CAT 5d',
                                                'CD4 Memory RAB2A 16h','CD4 Memory RAB2A 40h',
                                                'TCM CAT 5d','TCM RAB2A 16h','TEM NAPSA 40h',
                                                'TEM RALGDS 16h')))
  else if (outcome == 'C2')
    plot = plot + scale_y_discrete(limits=rev(c('CD4 Naive CAT 0h','CD4 Naive NAPSA 40h',
                                                'TN CAT 0h','TN NAPSA 40h',
                                                'CD4 Memory CAT 5d','TCM ADAM15 5d',
                                                'TCM CAT 5d','TEM NAPSA 40h','TEM RALGDS 16h')))
  
  print(plot)
  
  return(df %>% mutate(OR=exp(b),Lower_OR=exp(b-1.96*se),Upper_OR=exp(b+1.96*se)))
}
plotMRResults.BQC = function(dfs,solid.hits,outcome) { #take in list of IVW dfs
  df = dfs[[1]]
  if (length(dfs)>1) { for (item in 2:length(dfs)) df %<>% add_row(dfs[[item]]) }
  df$Gene = str_replace(df$Gene,'ENSG00000204472','AIF1')
  df$Gene = str_replace(df$Gene,'ENSG00000121691','CAT')
  df$Gene = str_replace(df$Gene,'ENSG00000176715','ACSF3')
  df$Gene = str_replace(df$Gene,'ENSG00000131400','NAPSA')
  df$Gene = str_replace(df$Gene,'ENSG00000104388','RAB2A')
  df$Gene = str_replace(df$Gene,'ENSG00000087076','HSD17B14')
  df$Gene = str_replace(df$Gene,'ENSG00000160271','RALGDS')
  
  out=outcome
  df %<>% mutate(OutcomeGeneState = paste(out,df$Gene,df$State),Solid=solid.hits)
  
  print(ggplot(df,aes(y=OutcomeGeneState,x=exp(b),xmin=exp(b-1.96*se),xmax=exp(b+1.96*se),color=Solid)) +
          geom_point() + geom_errorbar() + xlab('Odds ratio') +
          scale_y_discrete(limits=rev) + theme(text=element_text(size=20),legend.position="none") +
          geom_vline(xintercept=1) + ylab(''))
  return(df %>% mutate(OR=exp(b),Lower_OR=exp(b-1.96*se),Upper_OR=exp(b+1.96*se)))
}
plotMRResults.GTEx = function(dfs,solid.hits,outcome) { #take in list of IVW dfs
  df = dfs[[1]]
  if (length(dfs)>1) { for (item in 2:length(dfs)) df %<>% add_row(dfs[[item]]) }
  df$Gene = str_replace(df$Gene,'ENSG00000131400','NAPSA')
  
  out=outcome
  df %<>% mutate(OutcomeGene = paste(out,df$Gene),Solid=solid.hits)
  
  print(ggplot(df,aes(y=OutcomeGene,x=exp(b),xmin=exp(b-1.96*se),xmax=exp(b+1.96*se),color=Solid)) +
          geom_point() + geom_errorbar() + xlab('Odds ratio') +
          scale_y_discrete(limits=rev) + theme(text=element_text(size=16),legend.position="none") +
          geom_vline(xintercept=1) + ylab(''))
  return(df %>% mutate(OR=exp(b),Lower_OR=exp(b-1.96*se),Upper_OR=exp(b+1.96*se)))
}
