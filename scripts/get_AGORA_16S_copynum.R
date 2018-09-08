## Get AGORA 16S copy number
library(data.table)
gene_table = fread("data/blastDB/agora_NCBItax_processed_nodups.txt")
gene_list = fread("data/blastDB/agora_NCBI_gene_list.txt", header = F)
gene_list[,TaxID:=gsub("^>", "", V1)]
gene_list[,TaxID:=gsub(" lcl$", "", TaxID)]
gene_list[,TaxID:=gsub("_rna_from_genomic.fna.gz.1","", TaxID)]
num_copies = gene_list[,length(V2), by=TaxID]
setnames(num_copies, "V1", "CopyNum")
gene_table = merge(gene_table, num_copies, by.x = "databaseID", by.y = "TaxID", all.x = T)
## NAs had no annotated 16S genes? I guess
write.table(gene_table, file = "data/blastDB/agora_NCBItax_processed_nodups.txt", quote=F, row.names=F, sep = "\t")
