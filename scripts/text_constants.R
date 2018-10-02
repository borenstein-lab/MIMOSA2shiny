###### Text constants in UI
### NEED to be used to update mimosa package sysdata, which is where the app actually gets them
#### Save to package sysdata
#load("../mimosa2/R/sysdata.rda")
path_to_mimosa2 = commandArgs(trailingOnly = T)[1]
print(path_to_mimosa2)
file_out = paste0(path_to_mimosa2, "/R/sysdata.rda")

load(file_out)
microbiome_header = "Microbiome data"
database_title = "16S rRNA taxonomic data format"
database_choices = c("Sequence variants", "Greengenes 13_5 or 13_8", "SILVA", "No 16S rRNA data; use metagenome functional data only")
microbiome_input_title = "Upload 16S rRNA abundance file"
metagenome_option = "Metagenome KO abundances"
metagenome_input_title = "Upload file of metagenomic KO abundances"
#metagenome_use_option = "Use metagenome for core analysis instead of 16S rRNA data (see documentation)"

metabolome_header = "Metabolome data"
met_type_title = "Compound ID type"
met_type_choices = c("KEGG Compound IDs", "Metabolite names (search for matching ID)") #"MetaCyc Compound IDs",
selected_met_type = "KEGG Compound IDs"
metabolome_upload_title = "Upload metabolite abundance file"

network_header = "Model settings"
source_title = "Gene content and metabolic model source"
source_choices = c("PICRUSt KO genomes and KEGG metabolic model", "AGORA genomes and models")
net_mod_option = "Upload modifications to metabolic models"
net_mod_input_title = "Upload file of modifications at the species, gene, and/or reaction levels (optional)"
#closest_title = ""
#closest_options = c("Use closest AGORA species", "Use AGORA models for species within a % similarity threshold")
sim_title = "Minimum similarity threshold for mapping sequence variants"
gapfill_option = "Gap-fill model for each species using x program"
algorithm_header = "Algorithm settings"
stat_title = "Metabolite statistic to analyze:"
stat_choices = c("Variance (analytically calculated)") #, "Differential abundance (Wilcoxon rank-sum, permutation-based)", "Paired-sample differential abundance (paired Wilcoxon rank-sum, permutation based)")

### Replace previous sysdata
all_dat = ls()
all_dat = all_dat[all_dat != "file_out"]
save(list = all_dat, file = file_out)
