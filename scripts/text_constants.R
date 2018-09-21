###### Text constants in UI
### NEED to be used to update mimosa package sysdata, which is where the app actually gets them
#### Save to package sysdata
#load("../mimosa2/R/sysdata.rda")
path_to_mimosa2 = commandArgs(trailingOnly = T)[1]
print(path_to_mimosa2)

load(paste0(path_to_mimosa2, "/R/sysdata.rda"))
microbiome_header = "Microbiome data upload"
database_title = "16S rRNA taxonomic data format"
database_choices = c("Sequence variants", "Greengenes 13_5 or 13_8", "SILVA", "No 16S rRNA data; use metagenome functional data only (upload below)")
microbiome_input_title = "Upload 16S rRNA abundance file (example format linked here)"
metagenome_option = "Metagenome KO abundances"
metagenome_input_title = "Upload file of metagenomic KO abundances (example format linked here)"
#metagenome_use_option = "Use metagenome for core analysis instead of 16S rRNA data (see documentation)"

metabolome_header = "Metabolome data upload"
met_type_title = "Compound IDs:"
met_type_choices = c("KEGG Compound IDs", "Metabolite names (search for matching ID)") #"MetaCyc Compound IDs",
selected_met_type = "KEGG Compound IDs"
metabolome_upload_title = "Upload metabolite file (example format linked here)"

network_header = "Model settings"
source_title = "Gene content and metabolic model source"
source_choices = c("Assign KOs with PICRUSt and use KEGG metabolic model", "Map sequences to AGORA genomes and models")
net_mod_option = "Upload modifications to metabolic models"
net_mod_input_title = "Upload file of modifications at the species, gene, and/or reaction levels (see examples and documentation here)"
#closest_title = ""
#closest_options = c("Use closest AGORA species", "Use AGORA models for species within a % similarity threshold")
sim_title = "Similarity threshold:"
gapfill_option = "Gap-fill model for each species using x program"
algorithm_header = "Algorithm settings"
stat_title = "Metabolite statistic to analyze:"
stat_choices = c("Variance (analytically calculated)") #, "Differential abundance (Wilcoxon rank-sum, permutation-based)", "Paired-sample differential abundance (paired Wilcoxon rank-sum, permutation based)")

### Replace previous sysdata
save(list = ls(), file = paste0(path_to_mimosa2, "/R/sysdata.rda"))
