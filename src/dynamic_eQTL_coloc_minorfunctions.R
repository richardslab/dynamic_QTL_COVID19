replaceGeneNames = function(n) {
  new.names = n
  
  #rename sQTLs
  new.names[which(str_detect(new.names,'chr1:155210627:155212090') & str_detect(new.names,'ENSG00000173171'))] = 'MTX1 Splice Variant A'
  new.names[which(str_detect(new.names,'chr1:155921527:155921861') & str_detect(new.names,'ENSG00000132680'))] = 'KHDC4 Splice Variant A'
  new.names[which(str_detect(new.names,'chr8:60572113:60576234') & str_detect(new.names,'ENSG00000104388'))] = 'RAB2A Splice Variant A'
  new.names[which(str_detect(new.names,'chr8:60576297:60584208') & str_detect(new.names,'ENSG00000104388'))] = 'RAB2A Splice Variant B'
  new.names[which(str_detect(new.names,'chr8:60572113:60584208') & str_detect(new.names,'ENSG00000104388'))] = 'RAB2A Splice Variant C'
  new.names[which(str_detect(new.names,'chr12:112917700:112931878') & str_detect(new.names,'ENSG00000089127'))] = 'OAS1 Splice Variant A'
  new.names[which(str_detect(new.names,'chr12:112917700:112919487') & str_detect(new.names,'ENSG00000089127'))] = 'OAS1 Splice Variant B'
  new.names[which(str_detect(new.names,'chr12:112917700:112918597') & str_detect(new.names,'ENSG00000089127'))] = 'OAS1 Splice Variant C'
  new.names[which(str_detect(new.names,'chr12:112917700:112919389') & str_detect(new.names,'ENSG00000089127'))] = 'OAS1 Splice Variant D'
  new.names[which(str_detect(new.names,'chr12:112908824:112911098') & str_detect(new.names,'ENSG00000089127'))] = 'OAS1 Splice Variant E'
  new.names[which(str_detect(new.names,'chr12:112908824:112911051') & str_detect(new.names,'ENSG00000089127'))] = 'OAS1 Splice Variant F'
  new.names[which(str_detect(new.names,'chr12:112911258:112916509') & str_detect(new.names,'ENSG00000089127'))] = 'OAS1 Splice Variant G'
  new.names[which(str_detect(new.names,'chr12:132581824:132582062') & str_detect(new.names,'ENSG00000112787'))] = 'FBRSL1 Splice Variant A'
  new.names[which(str_detect(new.names,'chr19:10351162:10352926') & str_detect(new.names,'ENSG00000105397'))] = 'TYK2 Splice Variant A'
  new.names[which(str_detect(new.names,'chr19:10359302:10361511') & str_detect(new.names,'ENSG00000105397'))] = 'TYK2 Splice Variant B'
  new.names[which(str_detect(new.names,'chr19:10359302:10360962') & str_detect(new.names,'ENSG00000105397'))] = 'TYK2 Splice Variant C'
  new.names[which(str_detect(new.names,'chr21:33252830:33268394') & str_detect(new.names,'ENSG00000159110'))] = 'IFNAR2 Splice Variant A'
  new.names[which(str_detect(new.names,'chr21:33245074:33246254') & str_detect(new.names,'ENSG00000159110'))] = 'IFNAR2 Splice Variant B'
  new.names[which(str_detect(new.names,'chr21:33245074:33246718') & str_detect(new.names,'ENSG00000159110'))] = 'IFNAR2 Splice Variant C'
  new.names[which(str_detect(new.names,'chr21:33246482:33246718') & str_detect(new.names,'ENSG00000159110'))] = 'IFNAR2 Splice Variant D'
  new.names[which(str_detect(new.names,'chr21:33230216:33241886') & str_detect(new.names,'ENSG00000159110'))] = 'IFNAR2 Splice Variant E'
  new.names[which(str_detect(new.names,'chr21:33252830:33262561') & str_detect(new.names,'ENSG00000159110'))] = 'IFNAR2 Splice Variant F'
  new.names[which(str_detect(new.names,'chr21:33252830:33268394') & str_detect(new.names,'ENSG00000243646'))] = 'IL10RB Splice Variant A'
  
  #rename eQTLs
  new.names = str_replace(new.names,'ENSG00000104388','RAB2A')
  new.names = str_replace(new.names,'ENSG00000121691','CAT')
  new.names = str_replace(new.names,'ENSG00000130816','DNMT1')
  new.names = str_replace(new.names,'ENSG00000131400.7','NAPSA')
  new.names = str_replace(new.names,'ENSG00000131400','NAPSA')
  new.names = str_replace(new.names,'ENSG00000144802','NFKBIZ')
  new.names = str_replace(new.names,'ENSG00000166997','CNPY4')
  new.names = str_replace(new.names,'ENSG00000176715','ACSF3')
  new.names = str_replace(new.names,'ENSG00000213420','GPC2')
  new.names = str_replace(new.names,'ENSG00000232300','FAM215B')
  new.names = str_replace(new.names,'ENSG00000204472','AIF1')
  new.names = str_replace(new.names,'ENSG00000087076','HSD17B14')
  new.names = str_replace(new.names,'ENSG00000206503','HLA-A')
  new.names = str_replace(new.names,'ENSG00000160271','RALGDS')
  new.names = str_replace(new.names,'ENSG00000143537','ADAM15')
  new.names = str_replace(new.names,'ENSG00000177628','GBA')
  new.names = str_replace(new.names,'ENSG00000104388.15','RAB2A')
  new.names = str_replace(new.names,'ENSG00000121691','CAT')
  new.names = str_replace(new.names,'ENSG00000130816','DNMT1')
  new.names = str_replace(new.names,'ENSG00000131400.8','NAPSA')
  new.names = str_replace(new.names,'ENSG00000144802','NFKBIZ')
  new.names = str_replace(new.names,'ENSG00000166997','CNPY4')
  new.names = str_replace(new.names,'ENSG00000176715','ACSF3')
  new.names = str_replace(new.names,'ENSG00000213420','GPC2')
  new.names = str_replace(new.names,'ENSG00000232300','FAM215B')
  new.names = str_replace(new.names,'ENSG00000204472','AIF1')
  new.names = str_replace(new.names,'ENSG00000087076','HSD17B14')
  new.names = str_replace(new.names,'ENSG00000206503','HLA-A')
  new.names = str_replace(new.names,'ENSG00000160271','RALGDS')
  new.names = str_replace(new.names,'ENSG00000143537','ADAM15')
  new.names = str_replace(new.names,'ENSG00000177628','GBA')
  new.names = str_replace(new.names,'ENSG00000175164.16','ABO')
  new.names = str_replace(new.names,'ENSG00000142233.14','NTN5')
  new.names = str_replace(new.names,'ENSG00000169231.13','THBS3')
  new.names = str_replace(new.names,'ENSG00000160766.14','GBAP1')
  new.names = str_replace(new.names,'ENSG00000105397.14','TYK2')
  new.names = str_replace(new.names,'ENSG00000131401.11','NAPSB')
  new.names = str_replace(new.names,'ENSG00000120071.15','KANSL1')
  new.names = str_replace(new.names,'ENSG00000073969.18','NSF')
  new.names = str_replace(new.names,'ENSG00000127946.17','HIP1')
  new.names = str_replace(new.names,'ENSG00000159110.21','IFNAR2')
  new.names = str_replace(new.names,'ENSG00000188199.10','NUTM2B')
  new.names = str_replace(new.names,'ENSG00000262539.1','RP11-259G18.3')
  new.names = str_replace(new.names,'ENSG00000279996.1','JD275616')
  new.names = str_replace(new.names,'ENSG00000238160.1','LINC02863')
  new.names = str_replace(new.names,'ENSG00000214401.4','KANSL1-AS1')
  new.names = str_replace(new.names,'ENSG00000175164.13','ABO')
  new.names = str_replace(new.names,'ENSG00000142233.11','NTN5')
  new.names = str_replace(new.names,'ENSG00000127946.16','HIP1')
  new.names = str_replace(new.names,'ENSG00000105376.4','ICAM5')
  new.names = str_replace(new.names,'ENSG00000105397.13','TYK2')
  new.names = str_replace(new.names,'ENSG00000111331.12','OAS3')
  new.names = str_replace(new.names,'ENSG00000161847.13','RAVER1')
  new.names = str_replace(new.names,'ENSG00000236263.1','ENSG00000236263')
  return(new.names)
}
getSEQID = function(id) {
  seq.id = ''
  if (id == 'IL-4') seq.id = 'seq.2906.55'
  else if (id == 'IL-18') seq.id = 'seq.5661.15'
  else if (id == 'IL-23') seq.id = 'seq.10365.132'
  else if (id == 'IL-27') seq.id = 'seq.2829.19'
  else if (id == 'STAT1') seq.id = 'seq.10370.21'
  else if (id == 'CNPY4') seq.id = 'seq.15465.79'
  else if (id == 'Catalase') seq.id = 'seq.3488.64'
  else if (id == 'Glypican-2') seq.id = 'seq.3315.15'
  else if (id == 'RAB2A') seq.id = 'seq.16857.2'
  else if (id == 'AIF1') seq.id = 'seq.2849.49'
  else if (id == 'HSD17B14') seq.id = 'seq.13972.4'
  else if (id == 'ABO') seq.id = 'seq.9253.52'
  else if (id == 'ADAM15') seq.id = 'seq.15455.40'
  else if (id == 'IL10RB') seq.id = 'seq.2631.50'
  else if (id == 'OAS1') seq.id = 'seq.10361.25'
  else if (id == 'TYK2') seq.id = 'seq.5260.80'
  
  if (seq.id == '') stop('Missing seq.id')
  
  return(seq.id)
}