#AGORA-KEGG compound mapping

all_comps = unique(unlist(lapply(all_mods, function(x){ return(unlist(x$mets))})))

comp_table = rbindlist(lapply(agora_mods, function(x){
  return(data.table(met = unlist(x$mets), KEGG = sapply(x$metKEGGID, function(y){
    return(y[[1]][1])
  })))
}))
comp_table[KEGG=="[]", KEGG:=NA]
comp_table = unique(comp_table)
write.table(comp_table, file = "../MIMOSA2shiny/data/KEGGfiles/AGORA_KEGG_met_mappings.txt", quote=F, row.names=F, sep = "\t")


### AGORA Sequence database processing
seq_data = fread("../MIMOSA2shiny/data/blastDB/agora_NCBItax_mergedAssemblies.txt")
seq_data[,databaseID:=gsub(".*\\/", "", RefSeq_ftp_path)]
seq_data_nice = seq_data[,list(databaseID, Organism)]
seq_data_nice = seq_data_nice[!is.na(databaseID)]
write.table(seq_data_nice, file = "../MIMOSA2shiny/data/blastDB/agora_NCBItax_processed.txt", quote=F, row.names=F, sep = "\t")

all_agora_mods = list.files("../MIMOSA2shiny/data/AGORA/")
all_agora_mods = gsub("\\.mat$", "", all_agora_mods)
#all_agora_mods = gsub("_$", "", all_agora_mods) #This is a problem for reading things in
seq_data_nice[,AGORA_ID:=gsub(" ", "_",Organism)]
seq_data_nice[,AGORA_ID:=gsub("-", "_", AGORA_ID)]
seq_data_nice[,AGORA_ID:=gsub("\\.", "_", AGORA_ID)]
seq_data_nice[,AGORA_ID:=gsub("\\/", "_", AGORA_ID)]
seq_data_nice[,AGORA_ID:=gsub(":", "_", AGORA_ID)]
seq_data_nice[,AGORA_ID:=gsub("_=_.*", "", AGORA_ID)] #Remove extra strain IDs
seq_data_nice[,AGORA_ID:=gsub(",", "", AGORA_ID)]
seq_data_nice[,AGORA_ID:=gsub("\\'", "", AGORA_ID)]
seq_data_nice[,AGORA_ID:=gsub("\\)", "", AGORA_ID)]
seq_data_nice[,AGORA_ID:=gsub("\\(", "", AGORA_ID)]
seq_data_nice[,AGORA_ID:=gsub("__", "_", AGORA_ID)]

seq_data_nice[grepl("Bacillus_timonensis", AGORA_ID), AGORA_ID:="Peptoniphilus_timonensis_JC401"]
seq_data_nice[grepl("Gemella_moribillum", AGORA_ID), AGORA_ID:="Gemella_morbillorum_M424"]
seq_data_nice[grepl("Clostridium_saccharoperbutylacetonicum", AGORA_ID), AGORA_ID:="Clostridium_saccharoperbutylacetonicum_N1_4_HMT"]
seq_data_nice[grepl("Bifidobacterium_stercoris", AGORA_ID), AGORA_ID:="Bifidobacterium_stercoris_DSM_24849"]

seq_data_nice[!AGORA_ID %in% all_agora_mods, AGORA_ID:=paste0(AGORA_ID, "_")]
seq_data_nice[!AGORA_ID %in% all_agora_mods]
write.table(seq_data_nice, file = "../MIMOSA2shiny/data/blastDB/agora_NCBItax_processed.txt", quote=F, row.names=F, sep = "\t")
