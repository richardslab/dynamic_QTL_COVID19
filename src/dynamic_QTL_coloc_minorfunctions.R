getProxy = function(rs) { #minimum R2 for proxy is 0.8
  if (rs == 'rs10479009') return('No proxies') #best rs62385261 R2 0.7964
  else if (rs == 'rs708455') return('No proxies')
  else if (rs == 'rs564898521') return('No proxies') #not in reference panel LDLink
  else if (rs == 'rs200591918') return('No proxies') #not in ref panel
  else if (rs == 'rs9687004') return('No proxies') #not in ref panel
  else return('MissingRs')
}
replaceGeneNames = function(n) {
  new.names = n
  
  #rename eQTLs
  new.names = str_replace(new.names,'ENSG00000104388.15','RAB2A')
  new.names = str_replace(new.names,'ENSG00000104388','RAB2A')
  new.names = str_replace(new.names,'ENSG00000121691','CAT')
  new.names = str_replace(new.names,'ENSG00000130816','DNMT1')
  new.names = str_replace(new.names,'ENSG00000131400.7','NAPSA')
  new.names = str_replace(new.names,'ENSG00000131400.8','NAPSA')
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