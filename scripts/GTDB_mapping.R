## PRocess taxonomy mapping - from GTDB
# 6/1/2019
library(data.table)
tax_table = fread("data/rep_seqs/bac_metadata_r86.tsv")
tax_table = fread("data/rep_seqs/bac120_metadata_r86.2.tsv")
names(tax_table)

met_names = fread("data/rep_seqs/metadata_field_desc.tsv")

table_sub = tax_table[,list(ssu_gg_blast_subject_id, ssu_silva_blast_subject_id)]

ref_seq_table = fread("data/embl_gems/model_list_processed.txt")
agora_table = fread("data/blastDB/agora_NCBItax_processed_nodups.txt")

tax_table[accession=="RS_GCF_000003135.1", ssu_gg_blast_subject_id]
tax_table[,accession2:=gsub("^RS_", "", accession)]
tax_table[,accession2:=gsub("^GB_", "", accession2)]

prev_agora_gg_mapping = fread("data/rep_seqs/gg_13_8_99_toAGORA_97_map.txt")

agora_table[,accession2:=gsub("\\.1_.*", "\\.1", databaseID)]
agora_table[,accession2:=gsub("\\.2_.*", "\\.2", accession2)]
agora_table[,accession2:=gsub("\\.3_.*", "\\.3", accession2)]
agora_table[,table(accession2 %in% tax_table[,accession2])]
agora_table[!accession2 %in% tax_table[,accession2]]
prev_agora_gg_mapping[,table(databaseID %in% tax_table[,accession2])]

tax_table_sub = tax_table[,list(accession, accession2, ncbi_organism_name, ssu_count, ssu_gg_blast_subject_id, ssu_silva_blast_subject_id, gtdb_taxonomy)]
tax_table_sub = merge(tax_table_sub, ref_seq_table, by.x = "accession2", by.y = "assembly_accession", all = T)
tax_table_sub = merge(tax_table_sub, agora_table, by.x = "accession2", by.y = "accession2", all = T)
tax_table_sub = tax_table_sub[!is.na(AGORA_ID)|!is.na(ModelID)]
tax_table_sub[,table(ssu_gg_blast_subject_id == "none")]
tax_table_sub[,table(ssu_silva_blast_subject_id == "none")]

setnames(tax_table_sub, c("accession", "fullAccession", "ncbi_organism_name", "ssu_count_gtdb", "ggID_blast", "silvaID_blast", "gtdb_taxonomy", "ncbi_taxid", "organism_name", "strain_name", "file_path", "RefSeqCopyNum16S", "RefSeqModelID", "AGORA_databaseID", "AGORA_ID", "AGORA_CopyNum"))
write.table(tax_table_sub, file = "data/rep_seqs/GTDB_RefSeq_AGORA_SILVA_GG_mapping.txt", quote=F, row.names = F, sep = "\t")
