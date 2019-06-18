## Shapley analysis test
library(mimosa)
source("scripts/mimosa2_dev_functions.R")
orig_config_table = data.table(V1 = c("database", "genomeChoices","metType", "kegg_prefix", "data_prefix", "vsearch_path", "metagenome_format"),
                               V2 = c("Greengenes 13_5 or 13_8", "PICRUSt KO genomes and KEGG metabolic model", "KEGG Compound IDs", "data/KEGGfiles/", "data/", "bin/vsearch", F))
agora_config_table = data.table(V1 = c("database", "genomeChoices","metType", "kegg_prefix", "data_prefix", "vsearch_path", "metagenome_format"), 
                                V2 = c("Greengenes 13_5 or 13_8", "AGORA genomes and models", "KEGG Compound IDs", "data/KEGGfiles/", "data/", "/usr/local/bin/vsearch", F))

base_config_table = rbind(orig_config_table, data.table(V1 = "met_transform", V2 = "zscore"))
base_config_table_a = rbind(agora_config_table, data.table(V1 = "met_transform", V2 = "zscore"))
base_config_table_sim = copy(base_config_table)
base_config_table_sim[V1 == "genomeChoices", V2:="AGORA genomes and models"]
base_config_table_sim = rbind(base_config_table_sim, data.table(V1 = "manualAGORA", V2 = T))

#everything rank-based
base_config_table = rbind(base_config_table, data.table(V1 = "rankBased", V2 = T))
base_config_table_a = rbind(base_config_table_a, data.table(V1 = "rankBased", V2 = T))
base_config_table_sim = rbind(base_config_table_sim, data.table(V1 = "rankBased", V2 = T))

ten_spec_dat = process_abundances("data/testData/sim_data/allSpeciesEnv3.txt", "data/testData/sim_data/allMetabolitesEnv3.txt", 
                                  fluxes_file = "data/testData/sim_data/allMetFluxesEnv3.txt", simulated = T)

all_datasets = ten_spec_dat
config_table = base_config_table_sim
if(!"manualAGORA" %in% config_table[,V1]){
  network_results = build_metabolic_model(all_datasets[[1]], config_table, degree_filt = 0)
  network = network_results[[1]]
  species = network_results[[2]]
} else {
  network_results = build_metabolic_model(all_datasets[[1]], config_table, manual_agora = T, degree_filt = 0)
  network = network_results[[1]]
  species = network_results[[2]]
}
if("met_transform" %in% config_table[,V1]){
  met_transform = config_table[V1=="met_transform", V2]
} else met_transform = ""

mets_melt = melt(all_datasets[[2]], id.var = "compound", variable.name = "Sample")
if(met_transform != ""){
  mets_melt = transform_mets(mets_melt, met_transform)
}
indiv_cmps = get_species_cmp_scores(species, network, normalize = T, leave_rxns = F, manual_agora = ifelse("manualAGORA" %in% config_table[,V1], T, F), kos_only = F)
indiv_cmps = indiv_cmps[compound %in% mets_melt[,compound]]

nperm = 20
if("rankBased" %in% config_table[,V1]){
  rank_based = T
  if("rank_type" %in% config_table[,V1]){
    rank_type = config_table[V1=="rank_type", V2]
  } else rank_type = "rfit"
} else rank_based = F

species_cmps = indiv_cmps
cmp_mods = fit_cmp_mods(species_cmps, mets_melt, rank_based = rank_based, rank_type = rank_type)
mod_fit_true = cmp_mods[[1]][!is.na(Rsq) & PVal < signif_threshold,list(compound, Rsq)]
species_cmps = species_cmps[compound %in% mod_fit_true[,compound]]
#Fill in 0s
species_cmps = melt(dcast(species_cmps, compound+Sample~Species, value.var = "CMP", fill = 0), id.vars = c("compound", "Sample"), variable.name = "Species", value.name = "CMP")

spec_list = species_cmps[,sort(unique(as.character(Species)))]
nspec = length(spec_list)
R1 = sapply(1:nperm, function(x){
  sample.int(nspec)
}) # Matrix of permutations #fread(perm_file, header = F)[,get(paste0("V",perm_id))]

allContribs = data.table()
for(perm_id in 1:nperm){
  cumulMetVars = copy(mod_fit_true)
  setnames(cumulMetVars, "Rsq", "TrueRsq")
  spec_order = R1[,perm_id]
  for(j in 1:nspec){
    if(j < nspec){
      perm_dat = species_cmps[!Species %in% spec_list[spec_order[1:j]]] #Remove species incrementally
      #Fill in 0s for all species even if they don't do anything for that compound
      #When they are equal this works fine
      #fit model under permutation
      mod_fit1 = fit_cmp_mods(perm_dat, mets_melt, rank_based = rank_based, rank_type = rank_type)
      cumulMetVar = mod_fit1[[1]][,list(compound, Rsq)]
      cumulMetVars = merge(cumulMetVars, cumulMetVar, by = "compound", all.x = T)
      cumulMetVars[is.na(Rsq), Rsq:=0]
      setnames(cumulMetVars, "Rsq", spec_list[spec_order[j]])
    } else {
      cumulMetVars[,(spec_list[spec_order[j]]):=0]
    }
    if(j > 1){
      cumulMetVars[,paste0("Marg_", spec_list[spec_order[j]]):=get(spec_list[spec_order[j-1]]) - get(spec_list[spec_order[j]])]
    } else {
      cumulMetVars[,paste0("Marg_", spec_list[spec_order[j]]):=TrueRsq - get(spec_list[spec_order[j]])]
    }
  }
  cumulMetVars[,OrderID:=perm_id]
  cumulMetVars = cumulMetVars[,c("compound","TrueRsq", sort(names(cumulMetVars)[3:ncol(cumulMetVars)])), with=F]
  allContribs = rbind(allContribs, cumulMetVars, fill = T)
}
allContribs_mean = allContribs[,lapply(.SD, mean), by=compound, .SDcols = paste0("Marg_", spec_list)]
setnames(allContribs_mean, gsub("Marg_", "", names(allContribs_mean)))
allContribs_mean = melt(allContribs_mean, variable.name = "Species")
allContribs_mean = merge(allContribs_mean, mod_fit_true, by = "compound", all = T)

zero_cmps = species_cmps[,length(CMP[CMP != 0]), by=list(compound, Species)]
allContribs_mean = merge(allContribs_mean, zero_cmps, by = c("compound", "Species"), all = T)
allContribs_mean[V1==0 & value != 0]

shap_contribs = run_shapley_analysis(indiv_cmps, mets_melt, base_config_table_sim, nperm = 100)


ten_spec_dat = process_abundances("data/testData/sim_data/allSpeciesEnv3.txt", "data/testData/sim_data/allMetabolitesEnv3.txt", 
                                  fluxes_file = "data/testData/sim_data/allMetFluxesEnv3.txt", simulated = T)
if(!"manualAGORA" %in% config_table[,V1]){
  network_results = build_metabolic_model(all_datasets[[1]], config_table, degree_filt = 0)
  network = network_results[[1]]
  species = network_results[[2]]
} else {
  network_results = build_metabolic_model(all_datasets[[1]], config_table, manual_agora = T, degree_filt = 0)
  network = network_results[[1]]
  species = network_results[[2]]
}
if("met_transform" %in% config_table[,V1]){
  met_transform = config_table[V1=="met_transform", V2]
} else met_transform = ""

mets_melt = melt(all_datasets[[2]], id.var = "compound", variable.name = "Sample")
if(met_transform != ""){
  mets_melt = transform_mets(mets_melt, met_transform)
}
indiv_cmps = get_species_cmp_scores(species, network, normalize = T, leave_rxns = F, manual_agora = ifelse("manualAGORA" %in% config_table[,V1], T, F), kos_only = F)
indiv_cmps = indiv_cmps[compound %in% mets_melt[,compound]]
datasets2 = run_shapley_analysis()


foo = fit_cmp_mods(species_cmps = species_cmps, mets_melt = mets_melt, rank_based = T, rank_type = "svm")
tot_cmps = species_cmps[,sum(CMP), by=list(compound, Sample)]
tot_cmps = merge(tot_cmps, mets_melt[,list(compound, Sample, value)], by = c("compound", "Sample"))
all_comps = tot_cmps[,unique(compound)]
model_dat = data.table(compound = all_comps, Intercept = 0, Slope = 0, Rsq = 0, PVal = 0) #Make all other columns numeric
resid_dat = data.table(expand.grid(compound = all_comps, Sample = tot_cmps[,unique(Sample)]))
for(x in 1:length(all_comps)){
  if(!rank_based){
    scaling_mod = tryCatch(tot_cmps[compound==all_comps[x], lm(value~V1)], error=function(e){ NA})
  } else {
    if(rank_type == "mblm"){
      scaling_mod = tryCatch(tot_cmps[compound==all_comps[x], mblm::mblm(value~V1)], error=function(e){ NA})
    } else if (rank_type == "svm"){
      #Nu = 15% of sample size by default
      scaling_mod = tryCatch(tot_cmps[compound==all_comps[x], e1071::svm(V1, value, type = "nu-regression", nu = 0.15)], error=function(e){ NA})
    } else {
      scaling_mod = tryCatch(tot_cmps[compound==all_comps[x], Rfit::rfit(value~V1)], error=function(e){ NA})
    }
  }
  if(!identical(scaling_mod, NA)){
    scaling_coefs = coef(scaling_mod)
    scaling_resids = resid(scaling_mod)
    model_dat[x,Intercept:=scaling_coefs[1]]
    model_dat[x,Slope:=scaling_coefs[2]]
    if(!rank_based){
      model_dat[x,Rsq:=summary(scaling_mod)$r.squared]
      model_dat[x,PVal:=anova(scaling_mod)[["Pr(>F)"]][1]]
    } else {
      if(rank_type == "mblm"){
        model_dat[x,Rsq:=0]
        model_dat[x,PVal:=tryCatch(mblm::summary.mblm(scaling_mod)$coefficients[2,4], error=function(e){ NA})]
      } else { 
        mod_rsq = tryCatch(Rfit::summary.rfit(scaling_mod, overall.test = "drop")$R2, error=function(e){ NA})
        mod_pval = tryCatch(Rfit::drop.test(scaling_mod)$p.value, error=function(e){ NA})
        model_dat[x,Rsq:=mod_rsq]
        model_dat[x,PVal:=mod_pval]
      }
    }
    if(length(scaling_resids) != nrow(resid_dat[compound==all_comps[x]])) stop("Missing residuals")
    resid_dat[compound==all_comps[x], Resid:=scaling_resids]
  }
}


scaling_mod = tryCatch(tot_cmps[compound==all_comps[x], e1071::svm(V1, value, kernel = "linear")], error=function(e){ NA})
plot_dat = tot_cmps[compound==all_comps[x], list(V1, value)]
plot_dat[,SVMPred:=scaling_mod$fitted]
ggplot(plot_dat, aes(x=V1, y = value))+ geom_point(color = "red") + geom_point(aes(y=SVMPred), color = "black")
rank_reg = tot_cmps[compound==all_comps[x], Rfit::rfit(value~V1)]
ggplot(plot_dat, aes(x=V1, y = value))+ geom_point(color = "red") + geom_point(aes(y=SVMPred), color = "black") + geom_abline(aes(x=V1, y = value), slope = coef(rank_reg)[2], intercept = coef(rank_reg)[1]) +
  geom_smooth(method = "lm", se = F, alpha = 0.5)


scaling_mod2 = tot_cmps[compound==all_comps[x], liquidSVM::qtSVM(V1, value)]
plot_dat = tot_cmps[compound==all_comps[x], list(V1, value)]
plot_dat[,SVMPred:=predict(scaling_mod2, tot_cmps[compound==all_comps[x], V1])]
ggplot(plot_dat, aes(x=V1, y = value))+ geom_point(color = "red") + geom_point(aes(y=SVMPred), color = "black")
errors(scaling_mod2)
