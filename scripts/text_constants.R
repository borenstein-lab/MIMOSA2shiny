###### Text constants in UI
### NEED to be used to update mimosa package sysdata, which is where the app actually gets them
#### Save to package sysdata
#load("../mimosa2/R/sysdata.rda")
path_to_mimosa2 = commandArgs(trailingOnly = T)[1]
print(path_to_mimosa2)
file_out = paste0(path_to_mimosa2, "/R/sysdata.rda")

load(file_out)
microbiome_header = "Microbiome data"
database_title = "Microbiome data format"
database_choices = c("Sequence variants (ASVs)", "Greengenes 13_5 or 13_8 OTUs", "SILVA 132 OTUs", "No taxonomic abundance data; use functional abundance data only")
microbiome_input_title = "Upload taxonomic abundance file (16S rRNA)"
metagenome_title = "Functional abundance data format"
metagenome_options = c("Total KO abundances", "Taxon-stratified KO abundances (HUMAnN2 or PICRUSt/PICRUSt2)")
metagenome_input_title = "Upload file of KEGG Ortholog functional abundances"
#metagenome_use_option = "Use metagenome for core analysis instead of 16S rRNA data (see documentation)"
microbiome_tooltip = "Sequence variants are preferred. MIMOSA2 will map each taxon to any of the source choices described below."
microbiome_description = "Upload a dataset of taxonomic and/or functional abundances and select the corresponding format."
metagenome_tooltip = "A table of functional (KEGG KO) abundances is optional. It can either replace the taxonomic table or be used in comparison with the taxa-based analysis."

metabolome_header = "Metabolome data"
met_type_title = "Compound ID type"
met_type_choices = c("KEGG Compound IDs", "Metabolite names (search for matching ID)") #"MetaCyc Compound IDs",
selected_met_type = "KEGG Compound IDs"
metabolome_upload_title = "Upload metabolite abundance file"
metabolome_tooltip = "If KEGG IDs are not provided, MIMOSA2 will try to map compound names to KEGG IDs using MetaboAnalystR"
metabolome_description = "Upload a metabolite concentration table and select its metabolite ID format."
metabolome_norm_description = "Log transform metabolite values"
metabolome_transform_tooltip = "See the MIMOSA documentation for recommendations on whether to transform your data."

network_header = "Metabolic model settings"
source_title = "Gene content and metabolic model source"
source_choices = c("PICRUSt KO genomes and KEGG metabolic model", "AGORA genomes and models", "RefSeq/EMBL_GEMs genomes and models")
net_mod_option = "Upload modifications to metabolic models"
net_mod_input_title = "Upload file of modifications at the species, gene, and/or reaction levels (optional)"
network_tooltip = "See the documentation for descriptions and citations for each source option."
network_description = "Choose a source metabolic model template. If providing 16S rRNA sequence variants, choose a mapping threshold. Optionally, provide modifications to the metabolic model template - see the documentation for more details on formatting of modification files."
network_mapping_tooltip = "A higher threshold will result in fewer taxa more likely to be present. A lower threshold will include more species models but with more potential for error."

#closest_title = ""
#closest_options = c("Use closest AGORA species", "Use AGORA models for species within a % similarity threshold")
sim_title = "Minimum similarity threshold (for mapping ASVs only)"
gapfill_option = "Gap-fill model for each species using x program"
algorithm_header = "Algorithm settings"
algorithm_description = "Select the regression method for comparing metabolite concentrations with community metabolic potential (CMP) scores."
algorithm_tooltip = "Rank-based regression is recommended and is more robust and sensitive; OLS regression is much faster. See the documentation for more details."
stat_title = "Metabolite statistic to analyze"
stat_choices = c("Variance (analytically calculated)") #, "Differential abundance (Wilcoxon rank-sum, permutation-based)", "Paired-sample differential abundance (paired Wilcoxon rank-sum, permutation based)")
regression_title = "CMP-Metabolite model type"
regression_choices = c("Rank-based regression", "Least-squares (OLS) regression")

skip_contribs_option = "Skip taxonomic contribution analysis; only perform CMP-metabolite comparisons"
skip_contribs_tooltip = "Taxonomic contribution analysis can be slow when rank-based regression is selected."

result_table_description = "Each row of the table summarizes the MIMOSA2 results for a given metabolite, including the overall concordance between metabolic potential and concentration, and the contributing taxa."

### Replace previous sysdata
all_dat = ls()
all_dat = all_dat[all_dat != "file_out"]
save(list = all_dat, file = file_out)
