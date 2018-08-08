##### Make processed species-reaction mappings for both AGORA and PICRUSt
library(data.table)
library(mimosa)

args = commandArgs(trailingOnly = T)
database = args[1]

if(database == "PICRUSt"){
  ## Get list of all OTUs
  
  # Get normalization data
  picrust_normalization_table = fread(paste("gunzip -c ", picrust_norm_file, sep=""), header=T)
  colnames(picrust_normalization_table) = c("OTU", "norm_factor")
  picrust_normalization_table[,OTU:= as.character(OTU)]
  
  all_genomic_content = lapply(all_otus, function(x){
    return(get_genomic_content_from_picrust_table(otu, picrust_ko_table_directory, picrust_ko_table_suffix))
  })
  if(length(list.files(path = gsub("/reaction.*", "", kegg_paths[3]), pattern = "network_template.txt")) > 0){
    network_template = fread(paste0(gsub("/reaction.*", "", kegg_paths[3]), "/network_template.txt"))
  } else {
    all_kegg = get_kegg_reaction_info(kegg_paths[2], reaction_info_file = kegg_paths[3], save_out = F)
    network_template = generate_network_template_kegg(kegg_paths[1], all_kegg = all_kegg, write_out = F) #We should really speed this thing up
  }
  spec_table = rbindlist(lapply(all_otus, function(x){
    spec_mod = generate_genomic_network(contribution_table[OTU==x, unique(Gene)], keggSource = "KeggTemplate", degree_filter = 30, rxn_table = network_template, return_mats = F)
    spec_mod[,OTU:=x]
  }))
} else { #AGORA
  otu_list = list.files(path = "data/AGORA/", pattern = ".mat")
  otu_list = gsub(".mat$", "", otu_list)
  all_mods = load_agora_models(otu_list)
  all_mods = get_S_mats(all_mods, otu_list, edge_list = T)
  #Should we split up? Depends how big it is.
  for(spec in otu_list){
    write.table(all_mods[OTU==spec], file = paste0("data/AGORA/", spec, "_rxns.txt"), quote=F, row.names=F, sep = "\t")
  }

}