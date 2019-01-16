##### Make processed species-reaction mappings for both AGORA and PICRUSt
library(data.table)
library(mimosa) #, lib = "/data/shiny-server/r-packages/")

args = commandArgs(trailingOnly = T)
database = args[1]

if(database == "PICRUSt"){
  ## Get list of all OTUs
  picrust_ko_table_directory = "data/picrustGenomeData/indivGenomes/"
  picrust_ko_table_suffix = "_genomic_content.tab"
  all_otus = gsub(picrust_ko_table_suffix, "", list.files(picrust_ko_table_directory))

  picrust_norm_file = "data/picrustGenomeData/16S_13_5_precalculated.tab"
  #Get normalization data
  picrust_normalization_table = fread(picrust_norm_file, header = T)#fread(paste("gunzip -c ", picrust_norm_file, sep=""), header=T)
  colnames(picrust_normalization_table) = c("OTU", "norm_factor")
  picrust_normalization_table[,OTU:= as.character(OTU)]
  
  #Merge with reaction info
  kegg_paths = c("data/KEGGfiles/reaction_mapformula.lst", "data/KEGGfiles/reaction_ko.list", "data/KEGGfiles/reaction")
  if(length(list.files(path = gsub("/reaction.*", "", kegg_paths[3]), pattern = "network_template.txt")) > 0){
    network_template = fread(paste0(gsub("/reaction.*", "", kegg_paths[3]), "/network_template.txt"))
  } else {
    all_kegg = get_kegg_reaction_info(kegg_paths[2], reaction_info_file = kegg_paths[3], save_out = F)
    network_template = generate_network_template_kegg(kegg_paths[1], all_kegg = all_kegg, write_out = F) 
  }
  for(x in all_otus){
      genomic_content = get_genomic_content_from_picrust_table(x, picrust_ko_table_directory, picrust_ko_table_suffix)
      spec_mod = generate_genomic_network(genomic_content[,unique(Gene)], keggSource = "KeggTemplate", degree_filter = 0, rxn_table = network_template, return_mats = F) #Not going to filter anymore
      #Also get copy number info
      spec_mod = merge(spec_mod, genomic_content, by.x = "KO", by.y = "Gene", all.x=T)
      spec_mod = merge(spec_mod, picrust_normalization_table, by = "OTU", all.x=T, all.y=F)
      spec_mod[,normalized_copy_number:=copy_number/norm_factor]
      spec_mod = spec_mod[,list(OTU, KO, Reac, Prod, stoichReac, stoichProd, normalized_copy_number)]
      write.table(spec_mod, file = paste0("data/picrustGenomeData/indivModels/", x, "_rxns.txt"), quote=F, row.names=F, sep = "\t")
  }
} else { #AGORA - type
  
  dat_path = paste0("data/", database, "/") #could be AGORA, AGORA_EuropeanConstrained, etc
  print(dat_path)
  otu_list = list.files(path = dat_path, pattern = ".mat$")
  otu_list = gsub(".mat$", "", otu_list)
  print(otu_list)
  genome_info = fread("data/blastDB/AGORA_full_genome_info.txt")
  for(spec in otu_list){
    mod1 = load_agora_models(spec, agora_path = dat_path)
    mod1 = get_S_mats(mod1, spec, edge_list = T)
    ##Read back in and add copy number info????
    print(spec)
    mod1[,copy_number:=1]
    mod1 = merge(mod1, genome_info[,list(ModelAGORA, CopyNum)], all.x = T, by.x = "Species", by.y = "ModelAGORA")
    mod1[,normalized_copy_number:=ifelse(CopyNum==0|is.na(CopyNum), 1, copy_number/CopyNum)]
    mod1 = mod1[,list(Species, KO, Reac, Prod, stoichReac, stoichProd, normalized_copy_number, LB, UB, Rev)]
    write.table(mod1, file = paste0(dat_path, spec, "_rxns.txt"), quote=F, row.names=F, sep = "\t")
  }

}


##Debugging this - issue was 2 different versions of emm_to_edge_list function...
# spec1 = "Actinomyces_cardiffensis_F0333"
# comp1 = "12dgr180[e]"
# foo = load_agora_models(spec1, agora_path = "../Cecilia_server/MIMOSA2shiny/data/AGORA/")
# foo1 = get_S_mats(foo, spec1, edge_list  = T)
# foo2 = get_S_mats(foo, spec1, edge_list = F)
# emm1 = data.table(foo2[[1]], Compound = row.names(foo2[[1]]))
# edge_list2 = emm_to_edge_list(emm1)
# #Issue is with FBA_functions emm_to_edge_list
# all_S_mat = foo2
# all_S_mat = rbindlist(lapply(1:length(all_S_mat), function(x){
#   foo = emm_to_edge_list(data.table(all_S_mat[[x]], Compound = row.names(all_S_mat[[x]])))
#   foo[,Species:=spec1]
#   return(foo)
# }))
# compare_lists = merge(foo1, edge_list2, by = c("Reac", "KO", "Prod"), all = T)
# missing_edges = compare_lists[is.na(stoichReac.x) & is.na(stoichProd.x)]
# 
# foo4 = emm_to_edge_list(data.table(foo2[[1]], Compound = row.names(foo2[[1]])))
# net_melted = melt(emm1, id.var = "Compound")
# net_melted[Compound==comp1 & value != 0]
# net_melted = net_melted[value != 0]
# net_melted[,Prod:=ifelse(value > 0, Compound,0)]
# net_melted[,Reac:=ifelse(value < 0, Compound,0)]
# all_rxn_ids = net_melted[,unique(as.character(variable))]
# edge_list = data.table()
# for(k in 1:length(all_rxn_ids)){
#   rxn_sub = net_melted[variable==all_rxn_ids[k]]
#   if(nrow(rxn_sub[Prod !=0 ]) > 0 & nrow(rxn_sub[Reac != 0]) > 0){
#     edge_list_sub = data.table(expand.grid(rxn_sub[,unique(Reac[Reac !=0])], rxn_sub[,unique(Prod[Prod != 0])]))
#     setnames(edge_list_sub, c("Reac", "Prod"))
#   } else {
#     edge_list_sub = rxn_sub[,list(Reac, Prod)]
#     edge_list_sub[Reac==0, Reac:=NA]
#     edge_list_sub[Prod==0, Prod:=NA]
#   }
#   edge_list_sub[,KO:=all_rxn_ids[k]]
#   edge_list_sub[,stoichReac:=sapply(Reac, function(x){ return(rxn_sub[Compound==x,abs(value)])})]
#   edge_list_sub[,stoichProd:=sapply(Prod, function(x){ return(rxn_sub[Compound==x,value])})]
#   edge_list = rbind(edge_list, edge_list_sub)
# }

## DEbugging European constraints
# spec1 = "Abiotrophia_defectiva_ATCC_49176"
# otu_list = list.files(path = "data/AGORA/", pattern = ".txt$")
# otu_list = gsub("_rxns.txt$", "", otu_list)
# spec2 = sample(otu_list, 1)
# 
# mod1 = load_agora_models(spec2, agora_path = "../Cecilia_server/MIMOSA2shiny/data/AGORA/")[[1]]
# mod1c = load_agora_models(spec2, agora_path = "data/AGORA_EuropeanConstrained/")[[1]]
# table(mod1$lb)
# table(mod1c$lb < 0) #Ok, slightly more informative
# table(mod1c$lb)
# table(mod1$ub)
# table(mod1c$ub)
# edge_list = get_S_mats(list(mod1), spec1, edge_list = T)
# edge_list2 = get_S_mats(list(mod1c), spec1, edge_list = T) #Ok, identical
#We have just been assuming everything is not reversible, I guess
#Need to use lb and ub
#think about other ways to use constraints - specifically, lower bounds 
#Rules, GR rules: rules linking genes to reactions
