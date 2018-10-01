##### Make processed species-reaction mappings for both AGORA and PICRUSt
library(data.table)
library(mimosa, lib = "/data/shiny-server/r-packages/")

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
} else { #AGORA
  otu_list = list.files(path = "data/AGORA/", pattern = ".mat$")
  otu_list = gsub(".mat$", "", otu_list)
  print(otu_list)
  all_mods = load_agora_models(otu_list)
  all_mods = get_S_mats(all_mods, otu_list, edge_list = T)
  ##Read back in and add copy number info????
  genome_info = fread("data/blastDB/agora_NCBItax_processed_nodups.txt")
  all_mods[,copy_number:=1]
  all_mods = merge(all_mods, genome_info[,list(AGORA_ID, CopyNum)], all.x = T, by.x = "Species", by.y = "AGORA_ID")
  all_mods[,normalized_copy_number:=ifelse(is.na(CopyNum), 1, copy_number/CopyNum)]
  all_mods = all_mods[,list(Species, KO, Reac, Prod, stoichReac, stoichProd, normalized_copy_number)]
  for(spec in otu_list){
    write.table(all_mods[Species==spec], file = paste0("data/AGORA/", spec, "_rxns.txt"), quote=F, row.names=F, sep = "\t")
  }

}