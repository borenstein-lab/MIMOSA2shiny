#Run MIMOSA2 on MIMOSA1 datasets

library(data.table)
library(seqsetvis)
library(mimosa)
library(cowplot)
library(pROC)

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

clean_spec_contribs = function(spec_contribs, node_data, threshold = 0.1, varShare = F, threshold2 = 0.2, return_all = F){
  if(!"QValPos" %in% names(spec_contribs)){
    spec_contribs = merge(spec_contribs, node_data, by ="compound")
  }
  if(varShare == T){
    spec_contribs[,Pass:=ifelse(!is.na(VarShare) & VarShare > threshold2, 1, 0)]
  }
  if(return_all == F){
    spec_contribs_good = spec_contribs[Pass==1 & QValPos < threshold]
  } else {
    spec_contribs_good = spec_contribs
    if("VarShare" %in% names(spec_contribs_good)){
      spec_contribs_good[,VarShare2:=ifelse(QValPos < threshold, 0, VarShare)]
      spec_contribs_good[,PosVarShare:=ifelse(V1 > 0, V1/sum(V1[V1 > 0], na.rm=T),0), by=compound]
      
    } else {
      spec_contribs_good[,Cor2:=ifelse(QValPos < threshold, 0, Cor)]
      spec_contribs_good[is.na(Cor2), Cor2:=0]
    }
  }
  spec_contribs_good[,ContribPair:=paste0(compound, "_", Species)]
  return(spec_contribs_good)
}

get_overlaps = function(set_list){
  num_types = length(set_list)
  overlaps = data.table(expand.grid(names(set_list), names(set_list)))
  overlaps = overlaps[Var1 != Var2]
  overlaps[,NumShared:=sapply(1:nrow(overlaps), function(x){
    length(intersect(set_list[[overlaps[x,Var1]]], set_list[[overlaps[x,Var2]]]))
  })]
  overlaps[,NumVar1:=sapply(Var1, function(x){
    length(set_list[[x]])
  })]
  overlaps[,NumVar2:=sapply(Var2, function(x){
    length(set_list[[x]])
  })]
  return(overlaps)
}

get_overlaps_grid = function(set_list, all_features = NULL, melt = T){
  if(is.null(all_features)){
    all_features = unique(unlist(set_list)) #all mets or contribs
  }
  overlaps = data.table(ID = all_features)
  for(j in 1:length(names(set_list))){
    overlaps[,names(set_list)[j] := ifelse(ID %in% set_list[[j]], 1, 0)]
  }
  if(melt){
    overlaps = melt(overlaps, id.var = "ID", variable.name = "Method")
    ##Remove ones that are 0 for everything
    bad_mets = overlaps[,sum(value), by=ID][V1==0, ID]
    overlaps = overlaps[!ID %in% bad_mets]
  } else {
    overlaps = overlaps[rowSums(overlaps[,2:ncol(overlaps), with=F]) != 0]
  }
  return(overlaps)
}


compare_all_continuous = function(contrib_list){
  true_contribs = contrib_list[[1]]
  all_contribs = merge(true_contribs)
}

compare_mimosa_1_2_corrs = function(species_file, met_file, kos_file = NULL, fluxes_file = NULL, simulated = F, config_table_list = NULL, consist_threshold = 0.05, corr_threshold = 0.01, varShare_threshold = 0.2, varShare_threshold_met = 0.1){ #Set up with default then add options
  if(!simulated){
    #Set up
    if(is.null(kos_file)){
      filelist1 =  c(species_file, met_file)
      names(filelist1) = c("file1", "file2")
    } else {
      filelist1 = c(species_file, met_file, kos_file)
      names(filelist1) = c("file1", "file2", "metagenome")
    }
    if(is.null(config_table_list)){ #IF specific set of configs is not supplied, will follow a set ordering of main options
      config_table = data.table(V1 = c("database", "genomeChoices","metType", "kegg_prefix", "data_prefix", "vsearch_path", "metagenome_format"), 
                                V2 = c("Greengenes 13_5 or 13_8", "PICRUSt KO genomes and KEGG metabolic model", "KEGG Compound IDs", "data/KEGGfiles/", "data/", "bin/vsearch", F))
    } else { #Otherwise can also do custom config
      config_table = config_table_list[[1]] 
    }
    datafiles = read_mimosa2_files(file_list = filelist1, configTable = config_table, app =F)
    species = datafiles[[1]]
    mets = datafiles[[2]]
    
    ### MIMOSA2 
    results_kegg = run_mimosa2(config_table, species = species, mets = mets)
    if(is.null(config_table_list)){
      config_table[V1=="genomeChoices", V2:=get_text("source_choices")[2]]
    } else {
      config_table = config_table_list[[2]]
    }
    results_agora_noRev = run_mimosa2(config_table, species = species, mets = mets) 
    
    config_table = rbind(config_table, data.table(V1 = "revRxns", V2 = T))
    
    results_agora = run_mimosa2(config_table, species = species, mets = mets) 
      
    config_table = rbind(config_table, data.table(V1 = c( "rxnEdit"), V2 = T))
    mimosa2_results_edit_agora = run_mimosa2(config_table, species = species, mets = mets)
    
    config_table[V1=="genomeChoices", V2:=get_text("source_choices")[1]]
    config_table = config_table[V1 != "revRxns"] #Don't add rev rxns for KEGG, already accounted for
    mimosa2_results_edit_kegg = run_mimosa2(config_table, species = species, mets = mets)
    
    ### MIMOSA 1#../vol2_server/MIMOSA2shiny/
    contrib_table = generate_contribution_table_using_picrust(species, "data/picrustGenomeData/16S_13_5_precalculated.tab.gz", "data/picrustGenomeData/indivGenomes/", "_genomic_content.tab", copynum_column = T)
    genes = contrib_table[,sum(contribution), by=list(Sample, Gene)]
    genes = dcast(genes, Gene~Sample, value.var = "V1", fill = 0)
    setnames(genes, "Gene", "KO")
    setnames(mets, "compound", "KEGG")
    setkey(mets, "KEGG")
    setnames(contrib_table, c("contribution", "copy_number"), c("CountContributedByOTU", "GeneCountPerGenome"))
    file_id = gsub(".txt", "", basename(species_file))

    file_prefix = file_id
    rxn_table = fread("data/KEGGfiles/full_rxn_table.txt")
    run_all_metabolites(genes, mets, file_prefix = file_prefix, id_met = F, met_id_file = NULL,
                        net_method = "KeggTemplate", net_file = NULL, rxn_table_source = rxn_table,
                        correction = "fdr", degree_filter = 30, minpath_file = NULL, cor_method = "spearman", nperm = 1000, nonzero_filter = 4)
    load(paste0(file_prefix, "_out.rda"))
    spec_contribs = get_spec_contribs(contrib_table, data_dir = getwd(), results_file = paste0(file_prefix, "_out.rda"), out_prefix = file_prefix, otu_id = "all", valueVar = "singleMusicc", sum_to_genus = F, write_out = F, taxonomy_file = NULL, comparison = "cmps") 
    spec_contribs2 = get_spec_contribs(contrib_table, data_dir = getwd(), results_file = paste0(file_prefix, "_out.rda"), out_prefix = file_prefix, otu_id = "all", valueVar = "singleMusicc", sum_to_genus = F, write_out = F, taxonomy_file = NULL, comparison = "mets", met_data = mets) 
    spec_contribs3 = get_spec_contribs(contrib_table, data_dir = getwd(), results_file = paste0(file_prefix, "_out.rda"), out_prefix = file_prefix, otu_id = "all", valueVar = "singleMusicc", sum_to_genus = F, write_out = F, taxonomy_file = NULL, var_share = T)
    spec_contribs4 = get_spec_contribs(contrib_table, data_dir = getwd(), results_file = paste0(file_prefix, "_out.rda"), out_prefix = file_prefix, otu_id = "all", valueVar = "singleMusicc", sum_to_genus = F, write_out = F, taxonomy_file = NULL, cov_share = T, met_data = mets)
    #Varshares of cov
    
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
    if("KEGG" %in% names(mets)) setnames(mets, "KEGG", "compound")
    mets_melt = melt(mets, id.var = "compound", variable.name = "Sample")
    spec_met_corrs = basic_correlation_matrix(species_melt, mets_melt, method="spearman")
    spec_met_corrs[,Qval:=p.adjust(p.value, method="fdr")]
    ##Also do correlation plus gene presence-absence
    #Need species-specific net
    setkey(contrib_table, NULL)
    spec_met_links = unique(merge(contrib_table, ko_net[[3]], by.x = "Gene",by.y = "KO", all.x = T, allow.cartesian = T)[,list(OTU, Reac, Prod)])
    network_links = melt(spec_met_links[,list(OTU, Reac, Prod)], id.var = "OTU", value.name = "compound")
    spec_met_corrs[,GeneNum:=sapply(1:nrow(spec_met_corrs), function(x){
      return(nrow(network_links[OTU==spec_met_corrs[x,Species] & compound==spec_met_corrs[x,compound]]))
    })]
    spec_met_corrs[,hasGene:=ifelse(GeneNum > 0, 1, 0)]
    
    ##Make venn diagrams of shared consistent metabolites and key contributors?? I think this would be good
    consist_mets = node_data[QValPos < consist_threshold, compound]
    contrast_mets = node_data[QValNeg < consist_threshold, compound]
    m2kegg_mets = results_kegg[[2]][Rsq > varShare_threshold_met & Slope > 0, compound]
    m2kegg_contrast_mets= results_kegg[[2]][Rsq > varShare_threshold_met & Slope < 0, compound]
    m2agora_mets = results_agora[[2]][Rsq > varShare_threshold_met & Slope > 0, compound]
    m2agora_contrast_mets  = results_agora[[2]][Rsq > varShare_threshold_met & Slope < 0, compound]
    m2agora_noRev_mets = results_agora_noRev[[2]][Rsq > varShare_threshold_met & Slope > 0, compound]
    m2agora_noRev_contrast_mets = results_agora_noRev[[2]][Rsq > varShare_threshold_met & Slope < 0, compound]
    m2edit_agora_mets = mimosa2_results_edit_agora[[2]][Rsq > varShare_threshold_met & Slope > 0, compound]
    m2edit_kegg_mets = mimosa2_results_edit_kegg[[2]][Rsq > varShare_threshold_met & Slope > 0, compound]
    m2edit_agora_contrast_mets = mimosa2_results_edit_agora[[2]][Rsq > varShare_threshold_met & Slope < 0, compound]
    m2edit_kegg_contrast_mets = mimosa2_results_edit_kegg[[2]][Rsq > varShare_threshold_met & Slope < 0, compound]
    corr_mets = spec_met_corrs[,sum(Qval < corr_threshold), by=compound][V1 > 0, compound] # At least one correlated species
    corr_mets_gene = spec_met_corrs[,sum(Qval < corr_threshold & hasGene==T), by=compound][V1 > 0, compound ]
    consist_mets_compare = list(mimosa1 = consist_mets, 
                                m2agora = m2agora_mets,  
                                m2agora_noRev = m2agora_noRev_mets, m2kegg = m2kegg_mets, 
                                m2edit_agora = m2edit_agora_mets, m2edit_kegg = m2edit_kegg_mets,
                                correlation = corr_mets, correlation_gene = corr_mets_gene)
    contrast_mets_compare =  list(mimosa1 = contrast_mets, 
                                  m2agora = m2agora_contrast_mets,  
                                  m2agora_noRev = m2agora_noRev_contrast_mets, 
                                  m2kegg = m2kegg_contrast_mets, 
                                  m2edit_agora = m2edit_agora_contrast_mets, m2edit_kegg = m2edit_kegg_contrast_mets,
                                  correlation = corr_mets, correlation_gene = corr_mets_gene)
    
    plot1 = ssvFeatureVenn(consist_mets_compare[c("mimosa1", "m2agora", "m2kegg")])
    plot2 = ssvFeatureVenn(consist_mets_compare[c("mimosa1", "m2edit_kegg", "correlation")])
    plot3 = ssvFeatureEuler(consist_mets_compare) + annotate("text", x = rep(1.75, length(consist_mets_compare)), y = seq(3.5, 4.9, by = 1.4/(length(consist_mets_compare)-1)), label = paste0(names(consist_mets_compare), " ", sapply(consist_mets_compare, length)))
    compare_met_dat = get_overlaps(consist_mets_compare)
    #plot4 = ssvFeatureEuler(consist_mets_compare[])
    plot5 = ssvFeatureEuler(contrast_mets_compare) + annotate("text", x = rep(1.75, length(contrast_mets_compare)), y = seq(3.5, 4.9, by = 1.4/(length(contrast_mets_compare)-1)), label = paste0(names(contrast_mets_compare), " ", sapply(contrast_mets_compare, length)))
    
    contrib_list_cmps = clean_spec_contribs(spec_contribs, node_data, threshold = consist_threshold)[,ContribPair]
    contrib_list_mets = clean_spec_contribs(spec_contribs2, node_data, threshold = consist_threshold)[,ContribPair]
    contrib_list_varshares = clean_spec_contribs(spec_contribs3, node_data, varShare = T)[,ContribPair]
    contrib_list_covshares = clean_spec_contribs(spec_contribs4, node_data, varShare = T)[,ContribPair]
    mimosa2kegg = results_kegg$varShares[VarShare > varShare_threshold_met & Slope > 0 & Species != "Residual", paste0(compound, "_", Species)]
    mimosa2agora = results_agora$varShares[VarShare > varShare_threshold_met &  Slope > 0 & Species != "Residual", paste0(compound, "_", Species)]
    mimosa2kegg_edit = mimosa2_results_edit_kegg$varShares[VarShare > varShare_threshold_met &  Slope > 0 & Species != "Residual", paste0(compound, "_", Species)]
    mimosa2agora_edit = mimosa2_results_edit_agora$varShares[VarShare > varShare_threshold_met &  Slope > 0 & Species != "Residual", paste0(compound, "_", Species)]
    mimosa2agora_noRev = results_agora_noRev$varShares[VarShare > varShare_threshold_met &  Slope > 0 & Species != "Residual", paste0(compound, "_", Species)]
    
    corr_list = spec_met_corrs[Qval < corr_threshold, paste0(compound, "_", Species)]
    corr_list_gene = spec_met_corrs[Qval < corr_threshold & hasGene==1, paste0(compound, "_", Species)]

    
    contribs_compare = list(mimosa1 = contrib_list_cmps, mimosa1_mets = contrib_list_mets,
                            #mimosa1_varshares = contrib_list_varshares, mimosa1_covshares = contrib_list_covshares, 
                            mimosa2_agora_norev = mimosa2agora_noRev,
                            mimosa2kegg = mimosa2kegg, mimosa2agora = mimosa2agora, m2kegg_edit = mimosa2kegg_edit, m2agora_edit = mimosa2agora_edit,
                            correlation = corr_list, correlation_gene = corr_list_gene)
    plot6 = ssvFeatureEuler(contribs_compare) + annotate("text", x = rep(1.75, length(contribs_compare)), y = seq(3.5, 4.9, by = 1.4/(length(contribs_compare)-1)), label = paste0(names(contribs_compare), " ", sapply(contribs_compare, length)))

    compare_dat = get_overlaps(contribs_compare)
    save(list = ls(), file = paste0("results/", file_id, "_runAll.rda"))
    save_plot(plot3, file = paste0("results/", file_id, "_consistMets.png"), base_width = 6, base_height = 5.5)
    #save_plot(plot4, file = paste0("results/", file_id, "_contribs.png"), base_width = 6, base_height = 5.5)
    save_plot(plot5, file = paste0("results/", file_id, "_contrastMets.png"), base_width = 6, base_height = 5.5)
    return(list(plot1, plot2, plot3, plot4, compare_met_dat, compare_dat))
    
  } else {
    source(paste0(datadir2, "FBA_functions.R"))
    #Run on simulation data
    #Get S matrix
    species = fread(species_file)
    if(!"Sample" %in% names(species)){
      species[,Sample:=paste0("run", SimRun, "__TP", TimePoint)]
      if("noiseLevel" %in% names(species)){
        species[,Sample:=paste0(Sample, "_", noiseLevel*10)]
      }
    }
    species = dcast(species, Species~Sample, value.var = "value", fill = 0)
    setnames(species, "Species", "OTU")
    species = spec_table_fix(species)
    mets = fread(met_file)
    if("medium" %in% names(mets)){
      setnames(mets, "medium", "compound")
    }
    if(!"Sample" %in% names(mets)){
      mets[,Sample:=paste0("run", SimRun, "__TP", TimePoint)]
      if("noiseLevel" %in% names(mets)){
        mets[,Sample:=paste0(Sample, "_", noiseLevel*10)]
      }
    }
    
    mets = dcast(mets, compound~Sample, value.var = "value",fill = 0)
    mets = met_table_fix(mets)
    mets[,compound:=gsub("[env]", "[e]", compound, fixed = T)]
    
    if(!is.null(fluxes_file)) met_fluxes = fread(fluxes_file)
    contribs = getContributions(met_fluxes)
    contribs[,compound:=gsub("[env]", "[e]", compound, fixed = T)]
    #rm(met_fluxes)
    
    config_table = data.table(V1 = c("database", "genomeChoices","metType", "kegg_prefix", "data_prefix", "vsearch_path", "revRxns"), 
                              V2 = c("Greengenes 13_5 or 13_8", "AGORA genomes and models", "KEGG Compound IDs", "data/KEGGfiles/", "data/", "bin/vsearch", T))
    config_table = rbind(config_table, data.table(V1 = "manualAGORA", V2 = T))
    ##MIMOSA2
    mimosa2_results = run_mimosa2(config_table, species = species, mets = mets)
    #Get network to use for MIMOSA1
    network = build_metabolic_model(species, config_table = config_table, manual_agora = T, degree_filt = 0)[[1]]
    network = add_rev_rxns(network)
    
    config_table = rbind(config_table, data.table(V1 = c("rxnEdit"), V2 = T))
    mimosa2_results_edit = run_mimosa2(config_table, species = species, mets = mets)
    
    config_table = config_table[!V1 %in% c("rxnEdit", "revRxns")]
    m2_noRev = run_mimosa2(config_table, species = species, mets = mets)
    
    network_noRev = build_metabolic_model(species, config_table = config_table, manual_agora = T, degree_filt = 0)[[1]]
    
    mimosa1_results = run_all_metabolites_FBA2("HMPsims", fake_spec= species, fake_mets = mets, network = network, nperm = 1000, spec_codes = spec_codes)
    cmp_mat = dcast(mimosa1_results[[4]], compound~Sample, value.var = "value", fill = 0)
    cmp_mat = cmp_mat[compound %in% mets[,compound]]
    setkey(cmp_mat, compound)
    comps = cmp_mat[,compound]
    all_rxns = lapply(comps, function(x){ 
      if(nrow(network[Reac==x|Prod==x]) > 0){
        return(data.table(network[Reac==x|Prod==x], Reversible = 0))
      } else { return(NA) }
    })
    single_spec_cmps = get_species_cmp_scores(species, network, manual_agora = T)
    single_spec_cmps = dcast(single_spec_cmps, compound+Species~Sample, value.var = "CMP", fill = 0)
    single_spec_cmps = single_spec_cmps[compound %in% comps]
    setkey(single_spec_cmps, compound)
    spec_codes = data.table(Species = sort(unique(species[,OTU])), Code = sort(unique(species[,OTU])))
    single_spec_cmps = lapply(spec_codes[,Code], function(x){ return(single_spec_cmps[Species==x])})
    subjects = intersect(names(species), names(mets))
    
    spec_contribs1 = rbindlist(lapply(1:length(comps), cmp_species_contributions_picrust, cmps_sub_good = cmp_mat, all_rxns = all_rxns, subjects = subjects, norm_kos = "", ko_net = "", all_taxa = spec_codes[,Code], single_spec_cmps = single_spec_cmps, comparison = "cmps", met_data = mets))
    spec_contribs2 = rbindlist(lapply(1:length(comps), cmp_species_contributions_picrust, cmps_sub_good = cmp_mat, all_rxns = all_rxns, subjects = subjects, norm_kos = "", ko_net = "", all_taxa = spec_codes[,Code], single_spec_cmps = single_spec_cmps, comparison = "mets", met_data = mets))
    spec_contribs3 = var_shares_cmps(cmp_mat, all_rxns, subjects, norm_kos = "", ko_net = "", all_taxa = spec_codes[,Code], single_spec_cmps = single_spec_cmps)
    spec_contribs4 = var_shares_cmps(cmp_mat, all_rxns, subjects, norm_kos = "", ko_net = "", all_taxa = spec_codes[,Code], single_spec_cmps = single_spec_cmps, cov_shares = T, met_data = mets)
    
    mimosa1_results_noRev = run_all_metabolites_FBA2("HMPsims", fake_spec= species, fake_mets = mets, network = network_noRev, nperm = 1000, spec_codes = spec_codes)
    cmp_mat_nr = dcast(mimosa1_results_noRev[[4]], compound~Sample, value.var = "value", fill = 0)
    cmp_mat_nr = cmp_mat_nr[compound %in% mets[,compound]]
    setkey(cmp_mat_nr, compound)
    comps_nr = cmp_mat_nr[,compound]
    all_rxns_nr = lapply(comps_nr, function(x){ 
      if(nrow(network_noRev[Reac==x|Prod==x]) > 0){
        return(data.table(network_noRev[Reac==x|Prod==x], Reversible = 0))
      } else { return(NA) }
    })
    single_spec_cmps_nr = get_species_cmp_scores(species, network_noRev, manual_agora = T)
    single_spec_cmps_nr = dcast(single_spec_cmps_nr, compound+Species~Sample, value.var = "CMP", fill = 0)
    single_spec_cmps_nr = single_spec_cmps_nr[compound %in% comps_nr]
    setkey(single_spec_cmps_nr, compound)
    single_spec_cmps_nr = lapply(spec_codes[,Code], function(x){ return(single_spec_cmps_nr[Species==x])})

    spec_contribs1_nr = rbindlist(lapply(1:length(comps_nr), cmp_species_contributions_picrust, cmps_sub_good = cmp_mat_nr, all_rxns = all_rxns_nr, subjects = subjects, norm_kos = "", ko_net = "", all_taxa = spec_codes[,Code], single_spec_cmps = single_spec_cmps_nr, comparison = "cmps", met_data = mets))
    #spec_contribs2_nr = rbindlist(lapply(1:length(comps_nr), cmp_species_contributions_picrust, cmps_sub_good = cmp_mat_nr, all_rxns = all_rxns_nr, subjects = subjects, norm_kos = "", ko_net = "", all_taxa = spec_codes[,Code], single_spec_cmps = single_spec_cmps_nr, comparison = "mets", met_data = mets))
    spec_contribs3_nr = var_shares_cmps(cmp_mat_nr, all_rxns_nr, subjects, norm_kos = "", ko_net = "", all_taxa = spec_codes[,Code], single_spec_cmps = single_spec_cmps_nr)
    spec_contribs4_nr = var_shares_cmps(cmp_mat_nr, all_rxns_nr, subjects, norm_kos = "", ko_net = "", all_taxa = spec_codes[,Code], single_spec_cmps = single_spec_cmps_nr, cov_shares = T, met_data = mets)
    
    #Correlations
    species_melt = melt(species, id.var = "OTU", variable.name = "Sample")
    setnames(species_melt, "OTU", "Species")
    mets_melt = melt(mets, id.var = "compound", variable.name = "Sample")
    spec_met_corrs = basic_correlation_matrix(species_melt, mets_melt, method="spearman")
    spec_met_corrs[,Qval:=p.adjust(p.value, method="fdr")]
    network_links = melt(network[,list(OTU, Reac, Prod)], id.var = "OTU", value.name = "compound")
    #This was an issue
    network_links = unique(network_links[compound %in% mets[,compound],list(OTU, compound)]) #Don't need to keep reac/prod
    spec_met_corrs[,GeneNum:=sapply(1:nrow(spec_met_corrs), function(x){
      return(nrow(network_links[OTU==spec_met_corrs[x,Species] & compound==spec_met_corrs[x,compound]]))
      })]
    spec_met_corrs[,hasGene:=ifelse(GeneNum > 0, 1, 0)]
    contribs = merge(contribs, network_links, by.x = c("compound", "Species"), by.y = c("compound", "OTU"), all.x = T)
    #contribs[is.na(variable) & VarShare != 0]
    
    node_data = mimosa1_results[[2]]
    consist_mets = node_data[QValPos < consist_threshold, compound]
    contrast_mets = node_data[QValNeg < consist_threshold, compound]
    node_data2 = mimosa1_results_noRev[[2]]
    consist_mets2 = node_data2[QValPos < consist_threshold, compound]
    contrast_mets2 = node_data2[QValNeg < consist_threshold, compound]
    
    #m2kegg_mets = results_kegg[[2]][Rsq > 0.1, compound]
    m2agora_mets = mimosa2_results[[2]][Rsq > varShare_threshold_met & Slope > 0, compound]
    m2agora_contrast_mets = mimosa2_results[[2]][Rsq > varShare_threshold_met & Slope < 0, compound]
    m2edit_mets = mimosa2_results_edit[[2]][Rsq > varShare_threshold_met & Slope >0, compound]
    m2edit_contrast_mets = mimosa2_results_edit[[2]][Rsq > varShare_threshold_met & Slope < 0, compound]
    m2fwd_mets = m2_noRev[[2]][Rsq > varShare_threshold_met & Slope > 0, compound]
    m2fwd_contrast_mets = m2_noRev[[2]][Rsq > varShare_threshold_met & Slope < 0, compound]
    corr_mets = spec_met_corrs[,sum(Qval < corr_threshold), by=compound][V1 > 0, compound] # At least one correlated species
    corr_mets_gene = spec_met_corrs[,sum(Qval < corr_threshold & hasGene==T), by=compound][V1 > 0, compound ]
    
    contribs[,PosVarShare:=ifelse(V1 > 0, V1/sum(V1[V1 > 0], na.rm=T),0), by=compound]
    contribs[,VarShareMagnitude:=abs(V1)/sum(abs(V1)), by=compound]
    
    true_consist_mets = contribs[Species=="Inflow" & PosVarShare < (1-varShare_threshold_met), compound]
    true_contrast_mets = contribs[Species=="Inflow" & PosVarShare > (1-varShare_threshold_met), compound]
    
    contribs[Species=="Inflow", table(VarShare <(1-varShare_threshold))]
    contribs[Species=="Inflow", table(PosVarShare <(1-varShare_threshold))]
    contribs[Species=="Inflow", table(PosVarShare <(1-varShare_threshold_met))]
    
    m2_noRev[[2]][,table(Rsq > varShare_threshold)]
    
    consist_mets_compare = list(mimosa1 = consist_mets,mimosa1_noRev = consist_mets2, 
                                m2agora = m2agora_mets, m2agora_edit = m2edit_mets, m2agora_fwd = m2fwd_mets, correlation = corr_mets, correlation_gene = corr_mets_gene, 
                                true_consist = true_consist_mets)
    
    contrast_mets_compare = list(mimosa1 = contrast_mets,  mimosa1_noRev = contrast_mets2, 
                                m2agora = m2agora_contrast_mets, m2agora_edit = m2edit_contrast_mets, m2agora_fwd = m2fwd_contrast_mets, correlation = corr_mets, correlation_gene = corr_mets_gene, 
                                true_contrast = true_contrast_mets)
    
    plot1 = ssvFeatureVenn(consist_mets_compare[c("mimosa1", "m2agora_edit", "true_consist")])
    plot2 = ssvFeatureVenn(consist_mets_compare[c("true_consist", "mimosa1", "correlation")])
    plot3 = ssvFeatureEuler(consist_mets_compare) + annotate("text", x = rep(1.75, length(consist_mets_compare)), y = seq(3.5, 4.9, by = 1.4/(length(consist_mets_compare)-1)), label = paste0(names(consist_mets_compare), " ", sapply(consist_mets_compare, length)))
    plot5 = ssvFeatureEuler(contrast_mets_compare) + annotate("text", x = rep(1.75, length(contrast_mets_compare)), y = seq(3.5, 4.9, by = 1.4/(length(contrast_mets_compare)-1)), label = paste0(names(contrast_mets_compare), " ", sapply(contrast_mets_compare, length)))
    
    compare_met_dat = get_overlaps(consist_mets_compare)
    plot4 = ssvFeatureEuler(consist_mets_compare[c("mimosa1", "m2agora_edit", "correlation", "correlation_gene", "true_consist")])
    
    compare_grid_met = get_overlaps_grid(consist_mets_compare, all_features = mets[,compound], melt = T)
    compare_grid_met2 = get_overlaps_grid(consist_mets_compare, all_features = mets[,compound], melt = F)
    compare_grid_met2[,table(true_consist, m2agora_fwd)]
    compare_grid_met2[,table(true_consist, m2agora_edit)]
    compare_grid_met2[,table(true_consist, correlation)]
    compare_grid_met2[,table(true_consist, correlation_gene)]
    
    
    feature_order = compare_grid_met[Method=="true_consist"][order(value), ID]
    compare_grid_met[,ID:=factor(ID, levels = feature_order)]
    ggplot(compare_grid_met, aes(x=ID, y = Method, fill = factor(value))) + geom_tile()
    
    contrib_list_cmps = clean_spec_contribs(spec_contribs1, node_data, threshold = consist_threshold)[,ContribPair]
    #contrib_list_mets = clean_spec_contribs(spec_contribs2, node_data, threshold = consist_threshold)[,ContribPair]
    contrib_list_varshares = clean_spec_contribs(spec_contribs3, node_data, threshold = consist_threshold, varShare = T)[,ContribPair]
    contrib_list_covshares = clean_spec_contribs(spec_contribs4, node_data, threshold = consist_threshold, varShare = T)[,ContribPair]
    contrib_list_cmps2 = clean_spec_contribs(spec_contribs1_nr, node_data2, threshold = consist_threshold)[,ContribPair]
    #contrib_list_mets2 = clean_spec_contribs(spec_contribs2_nr, node_data2, threshold = consist_threshold)[,ContribPair]
    
    
    
    mimosa2_results$varShares[,ContribPair:=paste0(compound, "_", Species)]
    mimosa2_results_edit$varShares[,ContribPair:=paste0(compound, "_", Species)]
    m2_noRev$varShares[,ContribPair:=paste0(compound, "_", Species)]
    mimosa2_results_edit$varShares[,PosVarShare:=ifelse(V1 > 0, V1/sum(V1[V1 > 0], na.rm=T),0), by=compound]
    mimosa2_results$varShares[,PosVarShare:=ifelse(V1 > 0, V1/sum(V1[V1 > 0], na.rm=T),0), by=compound]
    m2_noRev$varShares[,PosVarShare:=ifelse(V1 > 0, V1/sum(V1[V1 > 0], na.rm=T),0), by=compound]
    
    
    mimosa2_contribs = mimosa2_results$varShares[VarShare > varShare_threshold & Slope > 0 & Species != "Residual", ContribPair]
    mimosa2_contribs_edit = mimosa2_results_edit$varShares[VarShare > varShare_threshold & Slope > 0 & Species != "Residual", ContribPair]
    mimosa2_contribs_edit_scaled = mimosa2_results_edit$varShares[PosVarShare > varShare_threshold_met & Slope > 0 & Species != "Residual", ContribPair]
    
    mimosa2_contribs_noRev = m2_noRev$varShares[VarShare > varShare_threshold & Slope > 0 & Species != "Residual", ContribPair]
    mimosa2_contribs_noRev_scaled = m2_noRev$varShares[PosVarShare > varShare_threshold_met & Slope > 0 & Species != "Residual", ContribPair]
    corr_list = spec_met_corrs[Qval < corr_threshold, paste0(compound, "_", Species)]
    corr_list_gene = spec_met_corrs[Qval < corr_threshold & hasGene==1, paste0(compound, "_", Species)]
    contribs[,ContribPair:=paste0(compound, "_", Species)]
    true_contribs = contribs[VarShare > varShare_threshold & Species != "Inflow", ContribPair]
    #try stricter true def also
    
    spec_met_corrs[,ContribPair:=paste0(compound, "_", Species)]
    spec_met_corrs[,AbsCor:=abs(estimate)]
    spec_met_corrs[,AbsCorGene:=ifelse(hasGene==1, AbsCor, 0)]
    
    contribs[,PosVarShare:=ifelse(V1 > 0, V1/sum(V1[V1 > 0], na.rm=T),0), by=compound]
    contribs[,VarShareMagnitude:=abs(V1)/sum(abs(V1)), by=compound]
    true_contribs_scaled = contribs[PosVarShare > varShare_threshold_met & Species != "Inflow", ContribPair]
    
    contribs_compare = list(mimosa1 = contrib_list_cmps,# mimosa1_mets = contrib_list_mets,
                            mimosa1_varshares = contrib_list_varshares, mimosa1_covshares = contrib_list_covshares, 
                            mimosa1_nr = contrib_list_cmps2, #mimosa1_nr_mets = contrib_list_mets2,
                            mimosa2 = mimosa2_contribs, mimosa2_edit = mimosa2_contribs_edit, mimosa2_edit_scale = mimosa2_contribs_edit_scaled, 
                            mimosa2_noRev = mimosa2_contribs_noRev, mimosa2_noRev_scaled = mimosa2_contribs_noRev_scaled, correlation = corr_list, correlation_gene = corr_list_gene, 
                            trueContribs = true_contribs_scaled)
    compare_grid2 = get_overlaps_grid(contribs_compare, all_features = contribs[,ContribPair], melt = F)
    compare_grid2[trueContribs==1 & rowSums(compare_grid2[,2:ncol(compare_grid2), with=F])==1]
    
    plot6 = ssvFeatureEuler(contribs_compare)+ annotate("text", x = rep(1.75, length(contribs_compare)), y = seq(3.5, 4.9, by = 1.4/(length(contribs_compare)-1)), label = paste0(names(contribs_compare), " ", sapply(contribs_compare, length)))
    plot5a = ssvFeatureEuler(contribs_compare[c("mimosa1", "mimosa2", "mimosa2_edit_scale", "correlation", "correlation_gene", "trueContribs")])
    plot5a = ssvFeatureEuler(contribs_compare[c("mimosa1_nr", "mimosa2_noRev", "mimosa2_edit_scale", "correlation", "correlation_gene", "trueContribs")])
    
    compare_dat = get_overlaps(contribs_compare)
    compare_grid = get_overlaps_grid(contribs_compare, all_features = contribs[,ContribPair], melt = T)
    feature_order = compare_grid[Method=="trueContribs"][order(value), ID]
    compare_grid[,ID:=factor(ID, levels = feature_order)]
    plot_grid = ggplot(compare_grid, aes(x=ID, y = Method, fill = factor(value))) + geom_tile()
    compare_grid2 = get_overlaps_grid(contribs_compare, all_features = contribs[,ContribPair], melt = F)
    compare_grid2[trueContribs==1 & rowSums(compare_grid2[,2:ncol(compare_grid2), with=F])==1]
    
    compare_grid2[,table(trueContribs, mimosa2)]
    compare_grid2[,table(trueContribs, mimosa2_noRev_scaled)]
    compare_grid2[,table(trueContribs, mimosa2_edit_scale)]
    compare_grid2[,prop.table(table(trueContribs, mimosa2_edit_scale), 2)]
    compare_grid2[,prop.table(table(trueContribs, mimosa2_noRev_scaled), 2)]
    compare_grid2[,prop.table(table(trueContribs, mimosa1_nr), 2)]
    compare_grid2[,prop.table(table(trueContribs, correlation), 2)]
    compare_grid2[,prop.table(table(trueContribs, correlation_gene), 2)]
    sens_spec_table = data.table()
    for(j in c("mimosa2_noRev_scaled","mimosa2_edit_scale", "mimosa1_nr", "correlation", "correlation_gene")){
      results_tab = data.table(Metric = j, Sensitivity = compare_grid2[trueContribs==1, sum(get(j)==1)/length(get(j))], 
                               Specificity = compare_grid2[trueContribs==0, sum(get(j)==0)/length(get(j))], 
                               Precision = compare_grid2[get(j)==1, sum(trueContribs==1)/length(trueContribs)])
      sens_spec_table = rbind(sens_spec_table, results_tab)
    }
    
    
    spec_contribs3 = clean_spec_contribs(spec_contribs3, node_data, threshold = consist_threshold, varShare = T, return_all = T)
    spec_contribs2 = clean_spec_contribs(spec_contribs2, node_data, threshold = consist_threshold, varShare = F, return_all = T)
    #spec_contribs2_nr = clean_spec_contribs(spec_contribs2_nr, node_data, threshold = consist_threshold, varShare = F, return_all = T)
      
    compare_grid3 = merge(contribs[,list(compound, Species, ContribPair, VarShare, PosVarShare)], mimosa2_results$varShares[,list(ContribPair, VarShare, PosVarShare)], by = "ContribPair", all = T)
    setnames(compare_grid3, c("VarShare.x", "PosVarShare.x", "VarShare.y", "PosVarShare.y"), c("TrueVarShare", "TruePosVarShare", "M2VarShare", "M2PosVarShare"))
    compare_grid3 = merge(compare_grid3, mimosa2_results_edit$varShares[,list(ContribPair, VarShare, PosVarShare)], by = "ContribPair", all = T)
    setnames(compare_grid3, c("VarShare", "PosVarShare"), c("M2Edit_VarShare", "M2Edit_PosVarShare"))
    compare_grid3 = merge(compare_grid3, spec_contribs3[,list(ContribPair, VarShare2, PosVarShare)], by = "ContribPair", all = T)
    setnames(compare_grid3, c("VarShare2", "PosVarShare"), c("M1_VarShare", "M1_PosVarShare"))
    compare_grid3 = merge(compare_grid3, m2_noRev$varShares[,list(ContribPair, VarShare, PosVarShare)], by = "ContribPair", all = T)
    setnames(compare_grid3, c("VarShare", "PosVarShare"), c("M2noRev_VarShare", "M2noRev_PosVarShare"))
    compare_grid3 = compare_grid3[Species != "Inflow"]
    compare_grid3[,BinaryContribPos:=ifelse(TruePosVarShare > varShare_threshold_met, 1, 0)]
    compare_grid3 = merge(compare_grid3, spec_met_corrs[,list(ContribPair, AbsCor, AbsCorGene)], by = "ContribPair", all = T)
    compare_grid3 = merge(compare_grid3, spec_contribs2[,list(ContribPair, Cor2)], by="ContribPair", all = T)
    setnames(compare_grid3, "Cor2", "M1_Corr")
    #compare_grid3 = merge(compare_grid3, spec_contribs2_nr[,list(ContribPair, Cor2)], by="ContribPair", all = T)
    #setnames(compare_grid3, "Cor2", "M1_noRev_Corr")
    
    roc_tab = data.table()
    roc_plots = data.table()
    roc_results = list()
    prec_recall = data.table()
    for(j in 1:(length(names(compare_grid3))-5)){
      if(names(compare_grid3)[j+5] != "BinaryContribPos"){
        roc_results[[j]] = roc(compare_grid3[,BinaryContribPos], compare_grid3[,get(names(compare_grid3)[j+5])])
        roc_tab = rbind(roc_tab, data.table(Metric = names(compare_grid3)[j+5], AUC = roc_results[[j]]$auc))
        roc_plots = rbind(roc_plots, data.table(Sens = roc_results[[j]]$sensitivities, Spec = roc_results[[j]]$specificities, Metric = names(compare_grid3)[j+5]))
        recall_dat = coords(roc_results[[j]], x = "all", ret = c("recall", "precision"))
        prec_recall = rbind(prec_recall, data.table(Recall = recall_dat["recall",], Precision = recall_dat["precision",], Metric = names(compare_grid3)[j+5]))
      }
    }
    plot6 = ggplot(roc_tab, aes(x=Metric, y = AUC)) + geom_bar(stat = "identity") + scale_y_continuous(expand = c(0,0), limits = c(0,1)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    plot7 = ggplot(roc_plots, aes(x=1-Spec, y = Sens, color = Metric)) + geom_line() + xlim(0, 0.25)
    plot8 = ggplot(prec_recall, aes(x=Recall, y = Precision, color = Metric)) + geom_line() + scale_x_continuous(limits = c(0,1), expand = c(0,0))
    plot8a = ggplot(prec_recall[Metric %in% c("AbsCor", "AbsCorGene", "M1_noRev_Corr", "M2noRev_PosVarShare")], aes(x=Recall, y = Precision, color = Metric)) + geom_line() + scale_x_continuous(limits = c(0,0.75), expand = c(0,0)) + scale_color_brewer(palette = "Set1", labels = c("Correlation", "Correlation w/ reaction presence", "MIMOSA1", "MIMOSA2"))
    
    plot8b = ggplot(prec_recall[Metric %in% c("AbsCor", "AbsCorGene", "M1_Corr", "M2Edit_PosVarShare")], aes(x=Recall, y = Precision, color = Metric)) + geom_line() + scale_x_continuous(limits = c(0,0.75), expand = c(0,0)) + scale_color_brewer(palette = "Set1", labels = c("Correlation", "Correlation w/ reaction presence", "MIMOSA1", "MIMOSA2"))
    
    compare_grid3[,table(BinaryContribPos, M2Edit_PosVarShare > 0.1)]
    file_id = gsub(".txt", "", basename(species_file))
    sens_spec_table[,Dataset:=file_id]
    write.table(sens_spec_table, file = paste0("results/fixed/", file_id, "sensSpec.txt"), quote=F, row.names = F, sep = "\t")
    
    save(list = ls(), file = paste0("results/fixed/", file_id, "_runAll.rda"))
    save_plot(plot3, file = paste0("results/fixed/", file_id, "_consistMets.png"), base_width = 6, base_height = 5.5)
    save_plot(plot5, file = paste0("results/fixed/", file_id, "_contrastMets.png"), base_width = 6, base_height = 5.5)
    save_plot(plot4, file = paste0("results/fixed/", file_id, "_contribs.png"), base_width = 6, base_height = 5.5)
    save_plot(plot6, file = paste0("results/fixed/", file_id, "_AUCs.png"), base_width = 6, base_height = 5.5)
    save_plot(plot7, file = paste0("results/fixed/", file_id, "_ROCs.png"), base_width = 6.5, base_height = 5.5)
    save_plot(plot8a, file = paste0("results/fixed/", file_id, "_precRecall.png"), base_width = 7, base_height = 3.8)

    ###Sensitivity/etc across different p-value thresholds
    # thresholds = c(0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.00000001, 10e-10, 10e-11, 10e-12, 10e-13)
    # sens_spec = data.table()
    # for(j in thresholds){
    #   sens_spec = rbind(sens_spec, data.table(spec_met_corrs[BinaryContribPos==0, sum(p.value > j)/length(p.value)], Stat="Specificity", Threshold = j), data.table(spec_met_corrs[BinaryContribPos==1, sum(p.value < j)/length(p.value)], Stat = "Sensitivity", Threshold = j), data.table(spec_met_corrs[ p.value < j, sum(BinaryContribPos)/length(p.value)], Stat = "PPV", Threshold = j), data.table(spec_met_corrs[, sum((p.value < j & BinaryContribPos==1)|(p.value > j & BinaryContribPos==0))/length(p.value)], Stat = "Accuracy", Threshold = j))
    # }
    
    return(list(plot1, plot2, plot3, plot4, compare_met_dat, compare_dat))
  }
}
 
 
#We should have this make linked bar plots of overlap to see total # of metabolites as well

  # consistent_mets = list(node_data[QValPos < 0.01 & PValPos < 0.01, compound], )
  # venn_consistent_mets = ssvFeatureVenn
  # 
  # spec_contrib_compare_set = list()
#}

datadir = "data/testData/mimosa1data/"
datadir2 = "data/testData/sim_data/"

args = commandArgs(trailingOnly = T)
spec_file = args[1]
met_file = args[2]
if("sim" %in% args){
  flux_file = args[3]
  sim = T
  ko_file = NULL
} else if(length(args) > 2){
  ko_file = args[3]
  sim = F
  flux_file = NULL
} else {
  ko_file = NULL
  flux_file = NULL
  sim = F
}

compare_mimosa_1_2_corrs(spec_file, met_file, kos_file = ko_file, fluxes_file = flux_file, simulated = sim)



sens_spec_table_plot = sens_spec_table
sens_spec_table2 = fread("~/Google Drive File Stream/My Drive/Genome_Sciences/MIMOSA2Data/HMP_test577_HMP_noise_specAbundssensSpec.txt")
sens_spec_table_plot = rbind(sens_spec_table_plot, sens_spec_table2)
sens_spec_table_plot = melt(sens_spec_table_plot, id.var = c("Metric", "Dataset"), variable.name = "Result")
sens_spec_table_plot[,Dataset2:=ifelse(Dataset == "HMP", "Dataset 2", "Dataset 1")]
sens_spec_plot = ggplot(sens_spec_table_plot[Metric %in% c("correlation", "correlation_gene", "mimosa1_nr", "mimosa2_noRev_scaled")], aes(x=Metric, y = value, shape = Dataset, fill = Metric)) + 
  geom_bar(stat = "identity", position = "dodge") + facet_grid(Dataset2~Result)+ 
  scale_fill_brewer(palette = "Set1", name = "Method", labels = c("Correlation", "Correlation w/ reaction presence", "MIMOSA1", "MIMOSA2"),
                    guide = guide_legend(nrow = 2)) +
  scale_x_discrete(name = "") + scale_y_continuous(expand = c(0,0)) + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
                                                                            strip.background = element_blank(), legend.position = "bottom")

sens_spec_plot2 = ggplot(sens_spec_table_plot[Metric %in% c("correlation", "mimosa1_nr", "mimosa2_noRev_scaled") & Result=="Precision"], aes(x=Metric, y = value, shape = Dataset, fill = Metric)) + 
  geom_bar(stat = "identity", position = "dodge") + facet_grid(~Dataset2)+ 
  scale_fill_manual(values = brewer.pal(4, "Set1")[c(1,2,4)], name = "Method", labels = c("Correlation", "MIMOSA1", "MIMOSA2"),
                    guide = guide_legend()) +
  scale_x_discrete(name = "") + scale_y_continuous(expand = c(0,0)) + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
                                                                            strip.background = element_blank(), legend.position = "right", legend.title = element_blank()) + ylab("Precision")

save_plot(sens_spec_plot, file = "simDatasets_compare_sensSpec.png", base_width = 6.5, base_height = 5)
