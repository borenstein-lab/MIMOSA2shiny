##### Make processed species-reaction mappings for both AGORA and PICRUSt
library(data.table)
library(mimosa)

args = commandArgs(trailingOnly = T)
database = args[1]

if(database == "PICRUSt"){
  ## Get list of all OTUs
  picrust_ko_table_directory = "data/picrustGenomeData/indivGenomes/"
  picrust_ko_table_suffix = "_genomic_content.tab"
  all_otus = gsub(picrust_ko_table_suffix, "", list.files(picrust_ko_table_directory))
  # picrust_norm_file = "data/picrustGenomeData/16S_13_5_precalculated.tab.gz"
  # Get normalization data
  # picrust_normalization_table = fread(paste("gunzip -c ", picrust_norm_file, sep=""), header=T)
  # colnames(picrust_normalization_table) = c("OTU", "norm_factor")
  # picrust_normalization_table[,OTU:= as.character(OTU)]
  
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
      spec_mod = generate_genomic_network(genomic_content[,unique(Gene)], keggSource = "KeggTemplate", degree_filter = 30, rxn_table = network_template, return_mats = F)
      #Also get copy number info
      spec_mod = merge(spec_mod, genomic_content, by.x = "KO", by.y = "Gene", all.x=T)
      write.table(spec_mod, file = paste0("data/picrustGenomeData/indivModels/", x, "_rxns.txt"), quote=F, row.names=F, sep = "\t")
  }
} else { #AGORA
  otu_list = list.files(path = "data/AGORA/", pattern = ".mat")
  otu_list = gsub(".mat$", "", otu_list)
  all_mods = load_agora_models(otu_list)
  all_mods = get_S_mats(all_mods, otu_list, edge_list = T)
  #Should we split up? Depends how big it is.
  for(spec in otu_list){
    write.table(all_mods[Species==spec], file = paste0("data/AGORA/", spec, "_rxns.txt"), quote=F, row.names=F, sep = "\t")
  }
  ##Read back in and add copy number info????

}