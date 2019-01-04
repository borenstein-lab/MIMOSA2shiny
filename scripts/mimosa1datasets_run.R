#Run MIMOSA2 on MIMOSA1 datasets

library(data.table)
datadir = "data/testData/mimosa1data/"

run_mimosa1 = function(genefile, metfile, contribfile, file_prefix, kegg_file){
  datasets = read_files(genefile, metfile)
  genes = datasets[[1]]
  mets = datasets[[2]]
  write_net = T
  net_method = 'KeggTemplate'
  degree_filter = 30
  cor_method = "spearman"
  rxn_table = fread(kegg_file)
  num_permute = 10000
  nonzero_filt = 3
  tax_file = ""
  sum_to_genus = F
  run_all_metabolites(genes, mets, file_prefix = file_prefix, id_met = F, met_id_file = "",
                      net_method = net_method, net_file = net_file, rxn_table_source = rxn_table,
                      correction = "fdr", degree_filter = degree_filter, minpath_file = "", cor_method = cor_method, nperm = num_permute, nonzero_filter = nonzero_filt)
  if(spec_contribs_all){
    spec_contribs = get_spec_contribs(contribs_file, data_dir = getwd(), results_file = paste0(file_prefix, "_out.rda"), out_prefix = file_prefix, otu_id = "all", valueVar = "singleMusicc", sum_to_genus = sum_to_genus, write_out = T, taxonomy_file = tax_file, comparison = "cmps") #will also save to file
  }
  load(paste0(file_prefix, "_out.rda"))
  good_mets = node_data[,compound]
  subjects = names(mets)[names(mets) != "KEGG"]
  if(gene_contribs_all){
    cmps = get_cmp_scores(ko_net[[1]], norm_kos)
    cmps_sub_good = cmps[compound %in% good_mets]
    all_rxns = lapply(good_mets, function(x){ return(get_non_rev_rxns(ko_net[[3]][Reac==x | Prod==x]))})
    all_gene_contribs = lapply(1:length(good_mets), gene_contributions, cmps_sub_good = cmps_sub_good, all_rxns = all_rxns,
                               subjects=subjects, norm_kos = norm_kos, ko_net = ko_net)
    all_ko_cors = rbindlist(lapply(all_gene_contribs, function(x){ return(x[[1]])}))
    all_comp_info = rbindlist(lapply(all_gene_contribs, function(x){ return(x[[2]])}))
  }
  
}

compare_mimosa_1_2_corrs = function(species_file, met_file, kos_file = "", simulated = F, config_table_list = NULL){ #Set up with default then add options
  if(!simulated){
    #Set up
    if(kos_file==""){
      filelist1 =  c(species_file, met_file)
      names(filelist1) = c("file1", "file2")
    } else {
      filelist1 = c(species_file, met_file, kos_file)
      names(filelist1) = c("file1", "file2", "metagenome")
    }
    if(is.null(config_table)){
      config_table = data.table(V1 = c("database", "genomeChoices","metType", "kegg_prefix", "data_prefix", "vsearch_path"), 
                                V2 = c("Greengenes 13_5 or 13_8", "PICRUSt KO genomes and KEGG metabolic model", "KEGG Compound IDs", "data/KEGGfiles/", "data/", "bin/vsearch"))
    } else {
      config_table = config_table_list[[1]]
    }
    datafiles = read_mimosa2_files(file_list = filelist1, configTable = config_table, app =F)
    species = datafiles[[1]]
    mets = datafiles[[2]]
    
    ### MIMOSA2 
    results_kegg = run_mimosa2(config_table)
    if(is.null(config_table)){
      config_table[V1=="genomeChoices", V2:=get_text("source_choices")[2]]
    } else {
      config_table = config_table_list[[2]]
    }
    results_agora = run_mimosa2(config_table)
    
    ### MIMOSA 1
    contrib_table = generate_contribution_table_using_picrust(species, "data/picrustGenomeData/16S_13_5_precalculated.tab.gz", "data/picrustGenomeData/indivGenomes/", "_genomic_content.tab")
    genes = contrib_table[,sum(contribution), by=list(Sample, Gene)]
    genes = dcast(genes, Gene~Sample, value.var = "V1")
    setnames(genes, "Gene", "KO")
    setnames(mets, "compound", "KEGG")
    setkey(mets, "KEGG")
    file_prefix = "mimosa1"
    rxn_table = fread("data/KEGGfiles/full_rxn_table.txt")
    run_all_metabolites(genes, mets, file_prefix = file_prefix, id_met = F, met_id_file = NULL,
                        net_method = "KeggTemplate", net_file = NULL, rxn_table_source = rxn_table,
                        correction = "fdr", degree_filter = 30, minpath_file = NULL, cor_method = "spearman", nperm = 1000, nonzero_filter = 4)
    load(paste0(file_prefix, "_out.rda"))
    spec_contribs = get_spec_contribs(contribs_table, data_dir = getwd(), results_file = paste0(file_prefix, "_out.rda"), out_prefix = file_prefix, otu_id = "all", valueVar = "singleMusicc", sum_to_genus = F, write_out = F, taxonomy_file = NULL, comparison = "cmps") 
    spec_contribs2 = get_spec_contribs(contribs_table, data_dir = getwd(), results_file = paste0(file_prefix, "_out.rda"), out_prefix = file_prefix, otu_id = "all", valueVar = "singleMusicc", sum_to_genus = F, write_out = F, taxonomy_file = NULL, comparison = "mets") 
    
        #Gene contribs
    load(paste0(file_prefix, "_out.rda"))
    good_mets = node_data[,compound]
    subjects = names(mets)[names(mets) != "KEGG"]
    cmps = get_cmp_scores(ko_net[[1]], norm_kos)
    cmps_sub_good = cmps[compound %in% good_mets]
    all_rxns = lapply(good_mets, function(x){ return(get_non_rev_rxns(ko_net[[3]][Reac==x | Prod==x]))})
    all_gene_contribs = lapply(1:length(good_mets), gene_contributions, cmps_sub_good = cmps_sub_good, all_rxns = all_rxns,
                               subjects=subjects, norm_kos = norm_kos, ko_net = ko_net)
    all_ko_cors = rbindlist(lapply(all_gene_contribs, function(x){ return(x[[1]])}))
    all_comp_info = rbindlist(lapply(all_gene_contribs, function(x){ return(x[[2]])}))
    
    ## Correlations
    species_melt = melt(species, id.var = "OTU", variable.name = "Sample")
    setnames(species_melt, "OTU", "Species")
    mets_melt = melt(mets, id.var = "KEGG", variable.name = "Sample")
    setnames(mets_melt, "KEGG", "compound")
    spec_met_corrs = basic_correlation_matrix(species_melt, mets_melt, method="spearman")
    
    ##Make venn diagrams of shared consistent metabolites and key contributors?? I think this would be good
  } else {
    
  }
}

#Also try comparing Mantel test with species var share contributions
species_file1 = paste0(datadir, "Dataset2_otu_table.txt")
mets_file1 = paste0(datadir, "Dataset2_mets.txt")
kos_file1 = paste0(datadir, "validation_picrust_good.txt")
#setnames(species, "Species", "OTU")

configTable = data.table(V1 = c("data_prefix", "database", "genomeChoices"), V2 = c("data/", "Greengenes 13_5 or 13_8", "PICRUSt KO genomes and KEGG metabolic model"))
network_results = build_metabolic_model(species1, configTable) #, input_data$netAdd) #input_data$geneAdd, 
network = network_results[[1]]
species1 = network_results[[2]] #Species should be unchanged in this case
rm(network_results)
species1 = dcast(species1, OTU~Sample, value.var="value")
indiv_cmps = get_species_cmp_scores(species1, network)


mets_melt = melt(mets1, id.var = "KEGG", variable.name = "Sample", value.var = "value")
setnames(mets_melt, "KEGG", "compound")
mets_melt[is.na(value), value:=0]
mets_melt[,value:=value/1000]
cmp_mods = fit_cmp_mods(indiv_cmps, mets_melt[,list(Sample, compound, value)])
indiv_cmps = add_residuals(indiv_cmps, cmp_mods[[1]], cmp_mods[[2]])

var_shares = calculate_var_shares(indiv_cmps)


species2 = fread(paste0(datadir, "BV_kos_qpcr_good.txt"))
mets2 = fread(paste0(datadir, "BV_mets_good.txt"))
#setnames(species, "Species", "OTU")

configTable = data.table(V1 = c("data_prefix", "database", "genomeChoices", "kegg_prefix"), V2 = c("data/", get_text("database_choices")[4], "PICRUSt KO genomes and KEGG metabolic model", "data/KEGGfiles/"))
network_results = build_metabolic_model(species2, configTable) #, input_data$netAdd) #input_data$geneAdd, 
network = network_results[[1]]
species2 = network_results[[2]] #Species should be unchanged in this case
indiv_cmps2 = get_cmp_scores_kos(species2, network)
mets_melt2 = melt(mets2, id.var = "KEGG", variable.name = "Sample")
setnames(mets_melt2, "KEGG", "compound")
cmp_mods2 = fit_cmp_mods(indiv_cmps2, mets_melt2)
indiv_cmps2 = add_residuals(indiv_cmps2, cmp_mods2[[1]], cmp_mods2[[2]])
var_shares_metagenome = calculate_var_shares(indiv_cmps2)
var_shares_metagenome
#shinyjs::logjs(devtools::session_info())
#Order dataset for plotting