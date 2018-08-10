###### Text constants in UI
microbiome_header = "Upload a tab-delimited file with a row for each species and a column for each sample."
database_title = "16S format:"
database_choices = c("Sequence variants (recommended for AGORA)", "Greengenes 13_5 or 13_8", "SILVA")
microbiome_input_title = "Upload 16S rRNA abundance file (example format linked here)"
metagenome_option = "Also upload metagenome KO abundances"
metagenome_input_title = "Upload file of metagenomic KO abundances (example format linked here)"

metabolome_header = "Upload a tab-delimited file with a row for each species and a column for each sample."
met_type_title = "Compound IDs:"
met_type_choices = c("KEGG Compound IDs", "MetaCyc Compound IDs", "Metabolite names (search for matching ID)")
selected_met_type = "KEGG Compound IDs"
metabolome_upload_title = "Upload metabolite file (example format linked here)"

network_header = "Model settings"
source_title = "Gene content and metabolic model source"
source_choices = c("Assign KOs with PICRUSt and use KEGG metabolic model", "Map sequences to AGORA genomes and models")
gene_mod_option = "Upload modifications to taxon-gene mapping"
gene_mod_input_title = "Upload file of species-gene modifications (see examples and documentation here)"
net_mod_option = "Upload modifications to taxon-reaction mapping"
net_mod_input_title = "Upload file of species-reaction modifications (see examples and documentation here)"
closest_title = ""
closest_options = c("Use closest AGORA species", "Use AGORA models for species within a % similarity threshold")
sim_title = "Similarity threshold:"
gapfill_option = "Gap-fill model for each species using x program"
algorithm_header = "Algorithm settings"
stat_title = "Metabolite statistic to analyze:"
stat_choices = c("Variance (analytically calculated)", "Differential abundance (Wilcoxon rank-sum, permutation-based)", "Paired-sample differential abundance (paired Wilcoxon rank-sum, permutation based)")

