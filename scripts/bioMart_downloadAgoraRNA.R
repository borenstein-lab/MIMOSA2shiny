### Use Biomart to download AGORA RNA
library(biomartr)
library(data.table)
genome_list = fread("data/AGORA_genomes/agora_genome_list.txt", header = F, sep = "\t")
#Take 2
all_agora_mods = list.files(path = "data/AGORA/", pattern = "_rxns.txt$")
all_agora_mods = gsub("_rxns.txt$", "", all_agora_mods)
#all_agora_mods = gsub("_$", "", all_agora_mods) #This is a problem for reading things in

org_mapping = fread("data/blastDB/agora1_StrainInfo.txt")
additional_orgs = fread("data/blastDB/agora2_newStrainInfo.txt")
setnames(additional_orgs, c("AGORA Model ID", "Strain"), c("ModelAGORA", "Organism"))
all_orgs = rbind(org_mapping, additional_orgs, fill = T)
all_orgs[!ModelAGORA %in% all_agora_mods]
all_agora_mods[!all_agora_mods %in% all_orgs[,ModelAGORA]]
all_orgs[ModelAGORA=="Leuconostoc_lactis_KACC 91922", ModelAGORA:="Leuconostoc_lactis_KACC_91922"]
all_orgs[ModelAGORA=="Gemella_moribillum_M424", ModelAGORA:="Gemella_morbillorum_M424"]
all_orgs[ModelAGORA=="Bacillus_timonensis_JC401", ModelAGORA:="Bacillus_timonensis_10403023"]
all_orgs[ModelAGORA=="Bifidobacterium_stercoris_ATCC_43183", ModelAGORA:="Bifidobacterium_stercoris_DSM_24849"]

all_orgs[,isAvail:=sapply(Organism, is.genome.available, db = "refseq")]
all_orgs[isAvail==F & !is.na(NCBI_ID),isAvailNCBI:=sapply(NCBI_ID, is.genome.available, db = "refseq")]
all_orgs[,Name:=gsub("\\.", "", Organism)]
all_orgs[,Name:=gsub("\\,", "", Name)]
all_orgs[isAvail==F & isAvailNCBI==F & Name != Organism, isAvailName2:=sapply(Name, is.genome.available, db = "refseq")]
all_orgs[isAvail==F & (isAvailNCBI==F|is.na(isAvailNCBI)),isAvailEnsembl:=sapply(Organism, is.genome.available, db = "ensemblgenomes")]
all_orgs[isAvail==F & (isAvailNCBI==F|is.na(isAvailNCBI)) & isAvailEnsembl==F, isAvailGB:=sapply(Organism, is.genome.available, db = "genbank")]
all_orgs[isAvailGB==F | (is.na(isAvailNCBI) & isAvail==F & isAvailEnsembl==F) & is.na(Name3), Name3:=c("187327", "Acetomicrobium hydrogeniformans ATCC BAA-1850", "Anaerostipes caccae DSM 14662", "Bacillus endophyticus", "Bacteroides massiliensis Timone 84634", "GCA_000159875.2", "GCA_000513195.1", "Blautia hansenii DSM 20583", "GCA_000168755.1", "GCA_002215605.1", "GCA_000027105.1", "GCA_000156055.1", "GCA_000158655.1", "GCA_000196455.1", "GCA_000186525.1", "GCA_000210115.1", "GCA_000165065.1", "GCA_000157115.2", "GCA_900209925.1", "GCA_000209955.1", "GCA_000209915.1", "GCA_000016305.1", "GCA_000163075.1", "GCA_000209755.1", "GCA_000612865.1", "Proteus penneri ATCC 35198", "GCA_000165085.1", "Blautia obeum A2-162", "GCA_000210035.1", "GCA_000006945.2", "GCA_900095765.1", 
                                "GCA_000196735.1", "GCA_000330805.1", "GCA_000011645.1", "GCA_000025805.1", "GCA_000009045.1", "Klebsiella michiganensis KCTC 1686", "GCA_000709265.1", "GCA_000012445.1", "GCA_000513215.1", "GCA_000632085.1", "GCA_000026585.1")]
#Manual search

all_orgs[is.na(isAvailName2), isAvailName2:=F]
all_orgs[Organism=="Melainabacterium MEL.A1", Name3:="GCA_001765415.1"]

get_all_seqs = sapply(796:nrow(all_orgs), function(x){
  f_path = paste0("data/AGORA_genomes/", all_orgs[x,ModelAGORA])
  system(paste0("mkdir ", f_path))
  if(all_orgs[x,isAvail==T]){
    getRNA(all_orgs[x,Organism], db = "refseq", path = f_path, reference = F)
  } else if(all_orgs[x, identical(isAvailNCBI, T)]){
    getRNA(all_orgs[x,NCBI_ID], db = "refseq", path = f_path, reference = F)
  } else if(all_orgs[x,isAvailName2==T]){
    getRNA(all_orgs[x,Name], db = "refseq", path = f_path, reference = F)
  } else if(all_orgs[x,identical(isAvailEnsembl, T)]){
    getRNA(all_orgs[x,Name], db = "ensemblgenomes", path = f_path, reference = F)
  } else if(all_orgs[x,identical(isAvailGB, T)]){
    getRNA(all_orgs[x,Name], db = "genbank", path = f_path, reference = F)
  } else if(is.genome.available(all_orgs[x,Name3], db = "refseq")){
    getRNA(all_orgs[x,Name3], db = "refseq", path = f_path, reference = F)
  } else {
    getRNA(all_orgs[x,Name3], db = "genbank", path = f_path, reference = F)
  }
})

all_orgs[grepl("Nitrososphaera", Name), Name3:="GCA_000303155.1"]
f_path = paste0("data/AGORA_genomes/", all_orgs[217,ModelAGORA])
getRNA(all_orgs[217,Name3], db = "genbank", path = f_path, reference = F)
# genome_list[,isAvail:=sapply(V1, is.genome.available, db = "refseq")]
# genome_list[,isAvailDB:=sapply(V1)]
# 
# genome_list[isAvail==T, sapply(V1, getRNA, db = "refseq", path = "data/AGORA_genomes/", reference = F)]
# 
# genomes_to_check = genome_list[isAvail==F]
# genomes_to_check[,isAvailEnsembl:=sapply(V1, is.genome.available, db = "ensemblgenomes")]
# genomes_to_check[isAvailEnsembl==T, sapply(V1, getRNA, db = "ensemblgenomes")]
# 
# genomes_continued = genomes_to_check[V1 > "Serratia fontincola"]
# genomes_continued[isAvailEnsembl==T, sapply(V1, getRNA, db = "ensemblgenomes")]
# 
# ensembl_list = list.files(path= "_ncbi_downloads/RNA/", pattern = ".gz")
# ensembl_list = sort(gsub("^_", "", ensembl_list))
# foo = cbind(ensembl_list, genomes_to_check[isAvailEnsembl==T, V1])
# weird_ones = genomes_to_check[isAvailEnsembl==T][c(3, 17, 37, 63, 145)]
# weird_ones[,sapply(V1, getRNA, db = "ensemblgenomes")]
# getRNA("Brevibacillus agri BAB-2500", db = "refseq", reference = F)
# getRNA("1219013", db = "refseq", reference = F) #NCBI tax ID for Rhodococcus equi
# weird_ones
# 
# genomes_to_check2 = genomes_to_check[isAvailEnsembl==F]
# genomes_to_check2[,Name2:=gsub(" [A-Za-z0-9]+$", "", V1)]
# genomes_to_check2[,isAvail2:=sapply(Name2, is.genome.available, db = "refseq")]
# genomes_to_check2[isAvail2==T, sapply(Name2, getRNA, db = "refseq", path = "data/AGORA_genomes/", reference = F)]
# 
# genomes_to_check3 = genomes_to_check2[isAvail2==F]
# genomes_to_check3[,Name3:=gsub(" [A-Za-z0-9]+ $", "",Name2)]
# genomes_to_check3[1,Name3:="Acinetobacter pittii PHEA-2"]
# genomes_to_check3[,isAvail3:=sapply(Name3, is.genome.available, db = "refseq")]
# 



all_gene_info = fread("data/AGORA_genomes/all_16S_counts.txt")
all_gene_info = merge(all_gene_info, all_orgs[,list(ModelAGORA, Organism)], by.x = "V1", by.y = "ModelAGORA", all = T)
setnames(all_gene_info, c("ModelAGORA", "CopyNum", "Organism"))
all_gene_info2 = all_gene_info[,max(CopyNum), by=list(ModelAGORA, Organism)]
setnames(all_gene_info2, "V1", "CopyNum")
write.table(all_gene_info2, file = "data/blastDB/AGORA_full_genome_info.txt", quote=F, row.names = F, sep = "\t")
