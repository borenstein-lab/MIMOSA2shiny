## RFit contributions benchmarking/cleaup
# 6/14/2019

library(data.table)
library(mimosa)

source("scripts/mimosa2_dev_functions.R")
orig_config_table = data.table(V1 = c("database", "genomeChoices","metType", "kegg_prefix", "data_prefix", "vsearch_path", "metagenome_format"),
                               V2 = c(get_text("database_choices")[2], get_text("source_choices")[1], "KEGG Compound IDs", "data/KEGGfiles/", "data/", "bin/vsearch", F))
agora_config_table = data.table(V1 = c("database", "genomeChoices","metType", "kegg_prefix", "data_prefix", "vsearch_path", "metagenome_format"), 
                                V2 = c(get_text("database_choices")[2], get_text("source_choices")[2], "KEGG Compound IDs", "data/KEGGfiles/", "data/", "/usr/local/bin/vsearch", F))

base_config_table = rbind(orig_config_table, data.table(V1 = "met_transform", V2 = "zscore"))
base_config_table_a = rbind(agora_config_table, data.table(V1 = "met_transform", V2 = "zscore"))
base_config_table_sim = copy(base_config_table)
base_config_table_sim[V1 == "genomeChoices", V2:=get_text("source_choices")[2]]
base_config_table_sim = rbind(base_config_table_sim, data.table(V1 = "manualAGORA", V2 = T))

#everything rank-based
base_config_table = rbind(base_config_table, data.table(V1 = "rankBased", V2 = T))
base_config_table_a = rbind(base_config_table_a, data.table(V1 = "rankBased", V2 = T))
base_config_table_sim = rbind(base_config_table_sim, data.table(V1 = "rankBased", V2 = T))

datasetIDs = c("bv", "asd", "abx", "metz", "ten_spec", "hmp", "ten_spec_noNoise")
ref_options = c("KEGG", "agora")
bv_datasets = process_abundances("data/testData/mimosa1data/Dataset2_otu_table.txt", "data/testData/mimosa1data/Dataset2_mets.txt")
abx_datasets = process_abundances("data/testData/mimosa1data/mice_otu_table.txt", "data/testData/mimosa1data/mice_mets_good.txt")
metz_datasets = process_abundances("data/testData/mimosa1data/metz_unc_otus_good.txt", "data/testData/mimosa1data/metz_unc_mets_good.txt")
asd_datasets = process_abundances("data/testData/miceASD_otus.txt", "data/testData/metsNMR_kegg.txt")

#simulated
ten_spec_dat = process_abundances("data/testData/sim_data/allSpeciesEnv3.txt", "data/testData/sim_data/allMetabolitesEnv3.txt", fluxes_file = "data/testData/sim_data/allMetFluxesEnv3.txt", simulated = T)
hmp_data = process_abundances("data/testData/sim_data/HMP_test577_HMP_noise_specAbunds.txt", "data/testData/sim_data/HMP_test577_HMP_noise_metAbunds.txt", fluxes_file = "data/testData/sim_data/HMP_test577_HMP_noise_metFluxesFinal.txt", simulated = T)
ten_spec_noNoise = process_abundances("data/testData/sim_data/allSpeciesFinalMain.txt", "data/testData/sim_data/allMetabolitesFinalMain.txt", fluxes_file = "data/testData/sim_data/allMetFluxesFinalMain.txt", simulated = T)

#contribs_long = rank_based_rsq_contribs(indiv_cmps, mets_melt, config_table, nperm = 2500, return_perm = T) 
run_contribs = function(dataset, configTable){
  species = dataset[[1]]
  mets = dataset[[2]]
  network_results = build_metabolic_model(species, config_table = configTable) #, input_data$netAdd) #input_data$geneAdd, 
  network = network_results[[1]]
  species = network_results[[2]] #Allow for modifying this for AGORA
  if("rxnEdit" %in% configTable[,V1]){
    rxn_param = T
    cat("Will refine reaction network\n")
  } else rxn_param = F
  rank_based = T
  cat("Will use rank-based/robust regression\n")
  rank_type = "rfit"
  if(configTable[V1=="database", V2==get_text("database_choices")[4]] & configTable[V1=="metagenome_format", V2==get_text("metagenome_options")[1]]){
    no_spec_param = T
    humann2_param = F
  } else if(configTable[V1=="database", V2==get_text("database_choices")[4]] & configTable[V1=="metagenome_format", V2==get_text("metagenome_options")[2]]){
    no_spec_param = F
    humann2_param = T
  } else {
    no_spec_param = F
    humann2_param = F
  }
  if("met_transform" %in% configTable[,V1]){
    met_transform = configTable[V1=="met_transform", V2]
    cat(paste0("Will transform metabolite values, transform is ", met_transform))
  } else met_transform = ""
  mets_melt = melt(mets, id.var = "compound", variable.name = "Sample")
  if(met_transform != ""){
    mets_melt = transform_mets(mets_melt, met_transform)
  }
  if(any(mets[,grep("[e]", compound, fixed = T)])){
    manual_agora1 = T
  } else manual_agora1 = F
  indiv_cmps = get_species_cmp_scores(species, network, normalize = T, leave_rxns = F, manual_agora = manual_agora1, kos_only = no_spec_param, humann2 = humann2_param)
  indiv_cmps = indiv_cmps[compound %in% mets[,compound]]
  cmp_mods = fit_cmp_mods(indiv_cmps, mets_melt, rank_based = rank_based, rank_type = rank_type)
  all_var_shares_dat = data.table()
  spec_dat = melt(dataset[[1]], id.var = "OTU", variable.name = "Sample")[,list(value/sum(value), OTU), by=Sample] #convert to relative abundance
  bad_spec = spec_dat[,list(length(V1[V1 != 0])/length(V1), max(V1)), by=OTU]
  bad_spec = bad_spec[V1 < 0.2 & V1 < 0.1, OTU] #Never higher than 10% and absent in at least 80% of samples
  time_dat = data.table(Iter = 1:10, Runtime = rep(0))
  for(m in 1:10){ # run shap with default # 10 times
    print(m)
    time1 = Sys.time()
    var_shares_dat = rank_based_rsq_contribs(indiv_cmps, mets_melt = mets_melt, config_table = configTable, cmp_mods = cmp_mods, 
                                             merge_low_abund = bad_spec, signif_threshold = 0.05)
    time2 = Sys.time()
    var_shares_dat[,PermuteRound:=m]
    all_var_shares_dat = rbind(all_var_shares_dat, var_shares_dat, fill = T)
    time_dat[m,Runtime:=difftime(time2, time1, units = "secs")]
    print(time_dat)
  }
  return(list(varShares = all_var_shares_dat, timeDat = time_dat))
}


config_table1 = base_config_table
config_table1_sim = base_config_table_sim
datasets = list(bv_datasets, asd_datasets, abx_datasets, metz_datasets, ten_spec_dat, hmp_data, ten_spec_noNoise)
contribs_list = lapply(datasets, function(x){ #list of empty data tables
  return(data.table())
})
timeDat = data.table()
for(i in 1:length(ref_options)){
    #Get correct base option
    print(ref_options[i])
    if(ref_options[i] == "KEGG"){
      config_table1 = copy(base_config_table)
    } else {
      config_table1 = copy(base_config_table_a)
    }
    for(k in 1:length(datasetIDs)){
      print(datasetIDs[k])
      if(!datasetIDs[k] %in% c("ten_spec", "hmp", "ten_spec_noNoise")){ # If not simulated, run
        print(config_table1)
        run_contribs_results = run_contribs(dataset = datasets[[k]], configTable = config_table1)
        run_contribs_dat = run_contribs_results$varShares
        run_time = run_contribs_results$timeDat
      } else if(ref_options[i] != "agora"){ 
        print(config_table1_sim)
        run_contribs_results = run_contribs(datasets[[k]], config_table1_sim)
        run_contribs_dat = run_contribs_results$varShares
        run_time = run_contribs_results$timeDat
      }
      run_contribs_dat[,Dataset:=datasetIDs[k]]
      run_contribs_dat[,Reference:=ref_options[i]]
      run_time[,Dataset:=datasetIDs[k]]
      run_time[,Reference:=ref_options[i]]
      run_time[,NumMets:=run_contribs_dat[,length(unique(compound))]]
      run_time[,NumSpec:=run_contribs_dat[,length(unique(Species))]]
      run_time[,NumSamps:=ncol(datasets[[k]][[1]])-1]
      print(run_time)
      timeDat = rbind(timeDat, run_time, fill = T)
      contribs_list[[k]] = rbind(contribs_list[[k]], run_contribs_dat, fill = T)
    }
}


#for each dataset
#Number of taxa
#Number of metabolites
#Sample size
#50*taxa orderings
#Time of analysis
#SD and cv of 5*t subsets

#Then repeat with merging low abundance taxa and compare

contribs_list[[1]][VarShare > 0.01,sd(value)/mean(value), by=list(compound, Species, Reference)] #2%

contribs_list[[1]][VarShare > 0.005,sd(value)/mean(value), by=list(compound, Species, Reference)][,summary(V1)] #2%
#< 5% for 86% of compounds and metabolites
contribs_list[[2]][VarShare > 0.005,sd(value)/mean(value), by=list(compound, Species, Reference)][,summary(V1)] #2%

contribs_list[[4]][VarShare > 0.005,sd(value)/mean(value), by=list(compound, Species, Reference)][,summary(V1)] #2%
 #2%
#Ok this did not run for the simulated datasets actually
for(j in 1:length(contribs_list)){
  print(contribs_list[[j]][VarShare > 0.005,sd(value)/mean(value), by=list(compound, Species, Reference)][,summary(V1)])
}