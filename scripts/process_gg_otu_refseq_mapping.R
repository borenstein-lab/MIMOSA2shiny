library(data.table)
refdir = "../data/rep_seqs/"
#vsearch_results_files = paste0(refdir, c("gg_13_8_99_to_agora_97id_all.txt", "silva_132_99_16S_to_agora_97id_all.txt"))
vsearch_results_files = paste0(refdir, c("gg_13_8_99_to_refseq_97id_all.txt", "silva_132_99_to_refseq_97id_all.txt"))
genome_info_file = "../data/blastDB/AGORA_full_genome_info.txt"
genome_info_file = "../data/embl_gems/model_list_processed.txt"
#out_files = paste0(refdir, c("gg_13_8_99_toAGORA_97_map.txt", "silva_132_99_toAGORA_97_map.txt"))
out_files = paste0(refdir, c("gg_13_8_99_toRefSeq_97_map.txt", "silva_132_99_toRefSeq_97_map.txt"))

genome_info = fread(genome_info_file)

for(j in 1:length(vsearch_results_files)){
	file_id = vsearch_results_files[j]
	foo = fread(file_id, header = F)
	foo[,assemblyID:=gsub("_lcl.*$", "", V2)]
	foo = unique(foo)
	map_all = merge(foo, genome_info, by.x = "assemblyID", by.y = "assembly_accession", all.x = T)
	map_all = map_all[,list(ModelID, assemblyID, V1, V2, V3, V4, CopyNum16S)]
	setnames(map_all, c("ModelID", "Ref_ID", "OTU", "fullGeneID", "matchPerc", "matchLen", "CopyNum"))
	write.table(map_all, file = out_files[j], quote=F, row.names=F, sep = "\t")
}

