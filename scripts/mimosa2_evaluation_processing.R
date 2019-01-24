# ##HMP sim specific
# spec_file = paste0(datadir2, "HMP_test577_HMP_noise_specAbunds.txt")
# met_file = paste0(datadir2, "HMP_test577_HMP_noise_metAbunds.txt")
# flux_file = paste0(datadir2, "HMP_test577_HMP_noise_metFluxesFinal.txt")
# #results1 = compare_mimosa_1_2_corrs(spec_file, met_file, kos_file = NULL, fluxes_file = flux_file, simulated = T)
# load("results/HMP_test577_HMP_noise_specAbunds_runAll.rda")
# mimosa2_results$modelData[,hist(Rsq, breaks = 30)] #Ouch
# 
# #Try adding posvarshare
# contribs[,PosVarShare:=ifelse(V1 > 0, V1/sum(V1[V1 > 0]),0), by=compound]
# contribs[,VarShareMagnitude:=abs(V1)/sum(abs(V1)), by=compound]
# 
# env_contribs = contribs[Species=="Inflow"]
# compare_consist = merge(mimosa2_results$modelData, env_contribs, by = "compound", all = T)
# 
# ggplot(compare_consist, aes(x=PosVarShare, y = Rsq)) + geom_point() + xlim(0,1)+ geom_smooth(method = lm)
# ggplot(compare_consist, aes(x=1-PosVarShare, y = Rsq, color = ifelse(Rsq > 0.2, T, F))) + geom_point() + xlim(0,1) 
# compare_consist[,table(1-PosVarShare > 0.2, Rsq > 0.2)]
# compare_consist[,MicrobialContrib:=1-PosVarShare]
# lm(Rsq~MicrobialContrib, data = compare_consist) ##Hmm
# compare_consist[Rsq < 0.2 & PosVarShare < 0.2, compound]
# 
# consist_threshold = 0.01 #Try strciter here
# consist_mets = node_data[QValPos < consist_threshold, compound]
# contrast_mets = node_data[QValNeg < consist_threshold, compound]
# #m2kegg_mets = results_kegg[[2]][Rsq > 0.1, compound]
# varShare_threshold_met = 0.1 #Looser here
# m2agora_mets = mimosa2_results[[2]][Rsq > varShare_threshold_met, compound]
# corr_threshold = 0.01
# corr_mets = spec_met_corrs[,sum(Qval < corr_threshold), by=compound][V1 > 0, compound] # At least one correlated species
# corr_mets_gene = spec_met_corrs[,sum(Qval < corr_threshold & hasGene==T), by=compound][V1 > 0, compound ]
# spec_met_corrs[,list(sum(Qval < corr_threshold), sum(Qval < corr_threshold & hasGene==T)), by=compound][V1 != 0 & V2 ==0] #16 don't meet the second criteria
# 
# true_threshold_met = 0.2
# true_consist_mets = contribs[Species=="Inflow" & VarShare < (1-true_threshold_met), compound]
# consist_mets_compare = list(mimosa1 = consist_mets, mimosa1contrast = contrast_mets, 
#                             m2agora = m2agora_mets, correlation = corr_mets, correlation_gene = corr_mets_gene, 
#                             true_consist = true_consist_mets)
# plot3 = ssvFeatureEuler(consist_mets_compare)
# compare_met_dat = get_overlaps(consist_mets_compare)
# contrib_compare[VarShare != 0 & hasGene==0, unique(compound)] #Ok now it's just a specific compound issue. weird. degree filter? Ok, this is fixed.
# 
# contrib_compare[Qval < corr_threshold & hasGene==F & VarShare > true_threshold_met]
# contrib_compare[Qval < corr_threshold & hasGene==T & VarShare > true_threshold_met]
# 
#What do we know about true contribs that aren't detected by anything?
contribs[,PosVarShare:=ifelse(V1 > 0, V1/sum(V1[V1 > 0], na.rm=T),0), by=compound]
contribs[,VarShareMagnitude:=abs(V1)/sum(abs(V1)), by=compound]

contrib_compare = merge(spec_met_corrs, contribs[Species != "Inflow"], by = c("Species", "compound"), all= T)

contrib_compare = merge(contrib_compare, mimosa2_results[[1]], by = c("compound", "Species"), all.x = T)
# 
contrib_compare[VarShare.x > 0.45 & Qval > 0.05 & VarShare.y < 0.1, table(as.character(compound))]
contrib_compare[VarShare.x > 0.45 & Qval > 0.05 & VarShare.y < 0.1, length(unique(compound))] #It is a subset of compounds
contrib_compare[VarShare.x > 0.45 & Qval > 0.05 & VarShare.y < 0.1, length(unique(Species))]
contrib_compare[VarShare.x > 0.45 & Qval > 0.05 & VarShare.y < 0.1, table(as.character(Species))] #A few top species as well - wonder what the deal is
contrib_compare[VarShare.x > 0.45 & Qval > 0.05 & VarShare.y < 0.1, list(table(as.character(Species)), names(table(as.character(Species))))][order(V1)]
# setkey(contrib_compare, NULL)
# dim(unique(contrib_compare)) #why are there mysterious duplicates again? 
# contrib_compare = unique()
contrib_compare[VarShare.x > 0.45 & Qval > 0.05 & VarShare.y < 0.1][order(VarShare.y)] #negative contributins
contrib_compare[VarShare.x > 0.45 & Qval > 0.05 & VarShare.y < 0.1][,table(VarShare.y < 0)] #negative contributins
contrib_compare[VarShare.x > 0.45 & Qval > 0.05 & VarShare.y < 0.1][order(abs(VarShare.y))] #negative contributins

bad_comps = contrib_compare[VarShare.x > 0.45 & Qval > 0.05 & VarShare.y < 0.1][,unique(compound)]
contrib_compare[,length(V1.x[V1.x != 0]), by=compound][,summary(V1), by=compound %in% bad_comps] #Ok it's not dramatically more species or anything

contrib_compare[VarShare.x > 0.45]
contrib_compare[VarShare.x > 0.45 ,summary(PosVarShare), by=p.value > 0.05 & VarShare.y < 0.1] #Ok yes this is a big part of it
contrib_compare[VarShare.x > 0.45 & PosVarShare > 0.45 & Qval > 0.05 & VarShare.y < 0.1]

cmps_w_rxns = get_species_cmp_scores(species, network, normalize = F, manual_agora = T, leave_rxns = T)
cmps_w_rxns = cmps_w_rxns[CMP != 0]
cmps_w_rxns = cmps_w_rxns[compound %in% mets[,compound]]
setkey(cmps_w_rxns, NULL)
cmps_w_rxns = unique(cmps_w_rxns) #Ok good

gly_cmps = cmps_w_rxns[compound=="gly[e]"]
gly_cmps = merge(gly_cmps, mets_melt[compound=="gly[e]"], by=c("Sample", "compound"))
gly_contribs = contribs[compound=="gly[e]" & VarShare > 0.2, unique(Species)]
ggplot(gly_cmps[Species %in% gly_contribs], aes(x=CMP, y = value.y , color = Species, shape = KO)) + geom_point() + facet_wrap(~Species, scales = "free")
ggplot(gly_cmps[Species %in% gly_contribs], aes(x=value.x, y = value.y , color = Species)) + geom_point() + facet_wrap(~Species, scales = "free")
#Cmps are just not representative of contribs
tot_gly_cmps = gly_cmps[,sum(CMP), by=list(Species, compound, Sample, value.x, value.y)]
ggplot(tot_gly_cmps[Species %in% gly_contribs], aes(x=V1, y = value.y , color = Species)) + geom_point() + facet_wrap(~Species, scales = "free")

ggplot(contrib_compare[compound=="gly[e]"], aes(x=VarShare.x, y = VarShare.y)) + geom_point()
ggplot(contrib_compare[compound=="gly[e]"], aes(x=PosVarShare, y = VarShare.y, color = ifelse(VarShare.x > 0.45, T, F))) + geom_point()
##Who knows. Back to subselection algorihtm
true_contribs_scaled = contribs[PosVarShare > varShare_threshold & Species != "Inflow", ContribPair]

contribs_compare = list(mimosa1 = contrib_list_cmps, mimosa1_mets = contrib_list_mets,
                        mimosa1_varshares = contrib_list_varshares, mimosa1_covshares = contrib_list_covshares, 
                        mimosa2 = mimosa2_contribs, correlation = corr_list, correlation_gene = corr_list_gene, 
                        trueContribs = true_contribs, trueContribs_scaled = true_contribs_scaled)
ssvFeatureEuler(contribs_compare)
ggplot(cmps_w_rxns[compound=="gly[e]"], aes(x=Sample, y = CMP, color = Species)) + geom_point()

# #Make wide model table for a single metabolite
# #Stepwise fit adding species?
#Try with 10 species
load("results/fixed/allSpeciesEnv5_runAll.rda")
#Ok network has bounds now



#For every reversible reaction, test scaling mod in each direction or both (no effect)
#For every reaction, test scaling mod with and without

met1 = "gly[e]"
network_orig = copy(network)
get_best_rxn_subset = function(met1, network, species, mets_melt){
  network = add_rev_rxns(network)
  cmps_w_rxns = get_species_cmp_scores(species, network, normalize = F, manual_agora = T, leave_rxns = T)[CMP != 0 & compound %in% mets[,compound]]
  cmps_w_rxns[,SpecRxn:=paste0(Species, "_", KO)]
  network[,SpecRxn:=paste0(OTU, "_", KO)]
  
  cmp_dat = cmps_w_rxns[compound==met1]
  met_net = network[Reac==met1|Prod==met1]
  
  uniq_rxns = cmp_dat[,unique(SpecRxn)]
  orig_scaling_mod = fit_scaling_mod(met1, cmp_dat, mets_melt)
  new_scaling_mod = orig_scaling_mod
  new_rsq = 1
  rsqs = rep(1, length(uniq_rxns))
  rxns_to_keep = c()
  rxns_to_remove = c()
  min_rxns = 3
  while(any(!uniq_rxns %in% rxns_to_keep) & length(uniq_rxns) > min_rxns & new_rsq > 1.1*summary(orig_scaling_mod)$r.squared){
    #If new one isn't that much better htan old one we're done
    orig_scaling_mod = new_scaling_mod
    met_net = met_net[!SpecRxn %in% rxns_to_remove]
    uniq_rxns = uniq_rxns[!uniq_rxns %in% c(rxns_to_remove, rxns_to_keep)]
    #update full model
    print(uniq_rxns)
    
    all_scaling_mods = list()
    for(j in 1:length(uniq_rxns)){
      cmp_dat1 = cmp_dat[SpecRxn != uniq_rxns[j]]
      all_scaling_mods[[j]] = fit_scaling_mod(met1, cmp_dat1, mets_melt)
    }
    slopes = sapply(all_scaling_mods, function(x){ return(x$coefficients[2])}) #Rule out changes that make this negative
    rxns_to_keep = c(rxns_to_keep, uniq_rxns[slopes < 0])
    rsqs = sapply(all_scaling_mods, function(x){ return(summary(x)$r.squared)}) #Really?
    rxns_to_remove = c(rxns_to_remove, sample(uniq_rxns[which(rsqs==max(rsqs))], 1)) #randomly pick 1 if multiple
    print(rxns_to_remove)
    cmp_dat = cmp_dat[!SpecRxn %in% rxns_to_remove]
    new_scaling_mod = fit_scaling_mod(met1, cmp_dat, mets_melt)
    print(summary(new_scaling_mod))
    new_rsq = summary(new_scaling_mod)$r.squared
  }
  
  #Ok we removed 5 reactions
  #Ok now if we get true contribs does it work better?
  #Also think about the idea that this is fixed for the rest of the network/metabolites - order of metabolites will matter a lot
  #Will have to report nice summary of rxns removed, rxns direction switched, etc
  scaling_coefs = coef(orig_scaling_mod)
  scaling_resids = resid(orig_scaling_mod)
  model_dat = data.table(compound = met1)
  resid_dat = data.table(expand.grid(compound = met1, Sample = cmp_dat[,unique(Sample)]))
  model_dat[compound==met1 ,Intercept:=scaling_coefs[1]]
  model_dat[compound==met1,Slope:=scaling_coefs[2]]
  model_dat[compound==met1,Rsq:=summary(orig_scaling_mod)$r.squared]
  if(length(scaling_resids) != nrow(resid_dat[compound==met1])) stop("Missing residuals")
  resid_dat[compound==met1, Resid:=scaling_resids]
  
  new_cmps = get_species_cmp_scores(species, network, normalize = F, manual_agora = T)
  new_cmps_test = new_cmps[compound==met1]
  new_cmps_test = add_residuals(new_cmps_test, model_dat, resid_dat) 
  var_shares = calculate_var_shares(new_cmps_test)
  #What does this look like?
  tot_cmps = new_cmps_test[,sum(newValue), by=list(compound, Sample)]
  tot_cmps = merge(tot_cmps, mets_melt, by = c("Sample", "compound"))
  ggplot(tot_cmps, aes(x=V1, y = value)) + geom_point()  
  flux_dat = fread(fluxes_file)
  flux_dat[,Sample:=paste0("run", SimRun, "__TP", TimePoint, "_", noiseLevel*10)]
  flux_dat[,compound:=gsub("[env]", "[e]", compound, fixed = T)]
  flux_dat1 = flux_dat[compound==met1]
  flux_dat1 = merge(flux_dat1, mets_melt, by = c("Sample", "compound"))
  ggplot(flux_dat1, aes(x=cumulFlux, y = value, color = Species)) + geom_point() + facet_wrap(~Species)
  #Ok let's see how it works for all compounds, and try with HMP dataset
  
}
#Somehow we need to fix a reaction if its direction gets set for one metabolite - must be the same for all
#If it gets removed for one metabolite - must get removed for all?
cmps_w_rxns[,table(grepl("_1", KO))] #Some
tot_cmps = cmps_w_rxns[,sum(CMP), by=list(Sample, compound)]
tot_cmps = merge(tot_cmps, mets_melt, by = c("Sample", "compound"))
#tot_cmps1 = merge(tot_cmps, )
scaling_mod = tot_cmps[compound==met1, lm(value~V1)]

model_dat = dcast(cmps_w_rxns[compound==met1], Sample~Species+KO, value.var = "CMP")
model_dat = merge(model_dat, mets_melt[compound==met1, list(Sample, value)], by="Sample")
lm(value~., data = model_dat[,2:ncol(model_dat), with=F])

# cmps_w_rxns

# contrib_list_cmps = clean_spec_contribs(spec_contribs1, node_data, threshold = consist_threshold)[,ContribPair]
# contrib_list_mets = clean_spec_contribs(spec_contribs2, node_data, threshold = consist_threshold)[,ContribPair]
# 
# #try lower var share spec threshold
# spec_contribs3[,Pass:=ifelse(VarShare > 0.1, 1, 0)]
# spec_contribs4[,Pass:=ifelse(VarShare > 0.1, 1, 0)]
# contrib_list_varshares = clean_spec_contribs(spec_contribs3, node_data, threshold = consist_threshold)[,ContribPair]
# contrib_list_covshares = clean_spec_contribs(spec_contribs4, node_data, threshold = consist_threshold)[,ContribPair]
# 
# m2_varShare_threshold = 0.05
# mimosa2_contribs = mimosa2_results$varShares[VarShare > m2_varShare_threshold & Species != "Residual", paste0(compound, "_", Species)]
# corr_list = spec_met_corrs[Qval < corr_threshold, paste0(compound, "_", Species)]
# corr_list_gene = spec_met_corrs[Qval < corr_threshold & hasGene==1, paste0(compound, "_", Species)]
# contribs[,ContribPair:=paste0(compound, "_", Species)]
# true_contribs = contribs[VarShare > varShare_threshold & Species != "Inflow", ContribPair]
# 
# 
# contribs_compare = list(mimosa1 = contrib_list_cmps, mimosa1_mets = contrib_list_mets,
#                         mimosa1_varshares = contrib_list_varshares, mimosa1_covshares = contrib_list_covshares, 
#                         mimosa2 = mimosa2_contribs, correlation = corr_list, correlation_gene = corr_list_gene, 
#                         trueContribs = true_contribs)
# plot4 = ssvFeatureEuler(contribs_compare)
# compare_dat = get_overlaps(contribs_compare)
# 
# 
# #Ok, trying smarter model fitting algorithm
# network = build_metabolic_model(species, config_table = config_table, manual_agora = T, degree_filt = 0)[[1]]
# cmps_w_rxns = get_species_cmp_scores(species_table = species, network = network, manual_agora = T, normalize = F, leave_rxns = T)
# cmps_w_rxns = cmps_w_rxns[CMP != 0]
# cmps_w_rxns = cmps_w_rxns[compound %in% mets[,compound]]
# setkey(cmps_w_rxns, NULL)
# cmps_w_rxns = unique(cmps_w_rxns)
# #Make wide model table for a single metabolite
# #Stepwise fit adding species?
# met1 = cmps_w_rxns[1,compound]
# model_dat = dcast(cmps_w_rxns[compound==met1], Species+KO~Sample, value.var = "CMP")
# 
#   

# species_file1 = paste0(datadir, "Dataset2_otu_table.txt")
# mets_file1 = paste0(datadir, "Dataset2_mets.txt")
#kos_file1 = paste0(datadir, "validation_picrust_good.txt")
#setnames(species, "Species", "OTU")
# 
# configTable = data.table(V1 = c("data_prefix", "database", "genomeChoices"), V2 = c("data/", "Greengenes 13_5 or 13_8", "PICRUSt KO genomes and KEGG metabolic model"))
# network_results = build_metabolic_model(species1, configTable) #, input_data$netAdd) #input_data$geneAdd, 
# network = network_results[[1]]
# species1 = network_results[[2]] #Species should be unchanged in this case
# rm(network_results)
# species1 = dcast(species1, OTU~Sample, value.var="value")
# indiv_cmps = get_species_cmp_scores(species1, network)
# 
# 
# mets_melt = melt(mets1, id.var = "KEGG", variable.name = "Sample", value.var = "value")
# setnames(mets_melt, "KEGG", "compound")
# mets_melt[is.na(value), value:=0]
# mets_melt[,value:=value/1000]
# cmp_mods = fit_cmp_mods(indiv_cmps, mets_melt[,list(Sample, compound, value)])
# indiv_cmps = add_residuals(indiv_cmps, cmp_mods[[1]], cmp_mods[[2]])
# 
# var_shares = calculate_var_shares(indiv_cmps)
# 
# 
# species2 = fread(paste0(datadir, "BV_kos_qpcr_good.txt"))
# mets2 = fread(paste0(datadir, "BV_mets_good.txt"))
# #setnames(species, "Species", "OTU")
# 
# configTable = data.table(V1 = c("data_prefix", "database", "genomeChoices", "kegg_prefix"), V2 = c("data/", get_text("database_choices")[4], "PICRUSt KO genomes and KEGG metabolic model", "data/KEGGfiles/"))
# network_results = build_metabolic_model(species2, configTable) #, input_data$netAdd) #input_data$geneAdd, 
# network = network_results[[1]]
# species2 = network_results[[2]] #Species should be unchanged in this case
# indiv_cmps2 = get_cmp_scores_kos(species2, network)
# mets_melt2 = melt(mets2, id.var = "KEGG", variable.name = "Sample")
# setnames(mets_melt2, "KEGG", "compound")
# cmp_mods2 = fit_cmp_mods(indiv_cmps2, mets_melt2)
# indiv_cmps2 = add_residuals(indiv_cmps2, cmp_mods2[[1]], cmp_mods2[[2]])
# var_shares_metagenome = calculate_var_shares(indiv_cmps2)
# var_shares_metagenome
#shinyjs::logjs(devtools::session_info())
#Order dataset for plotting


## Testing new functions 1/17
load("results/fixed/HMP_test577_HMP_noise_specAbunds_runAll.rda")
load("results/fixed/allSpeciesEnv5_runAll.rda")
config_table = data.table(V1 = c("database", "genomeChoices","metType", "kegg_prefix", "data_prefix", "vsearch_path"),
                          V2 = c("Greengenes 13_5 or 13_8", "AGORA genomes and models", "KEGG Compound IDs", "data/KEGGfiles/", "data/", "bin/vsearch"))
config_table = rbind(config_table, data.table(V1 = "manualAGORA", V2 = T))
config_table = rbind(config_table, data.table(V1 = c("revRxns", "rxnEdit"), V2 = c(T, T)))
config_table = check_config_table(config_table, app = T)
data_inputs = list(species = species, mets = mets)
network_results = build_metabolic_model(species, config_table, manual_agora = T)
network = network_results[[1]]
species = network_results[[2]]
network = add_rev_rxns(network)
mets_melt = melt(mets, id.var = "compound", variable.name = "Sample")


#cmp_mods =  fit_cmp_net_edit(network, species, mets_melt, manual_agora = T)
#new_network = cmp_mods[[3]]
#indiv_cmps = cmp_mods[[4]]
#break it down

species_cmps = get_species_cmp_scores(species, network, normalize = F, manual_agora = T, leave_rxns = T)[compound %in% mets_melt[,unique(compound)]]
species_cmps[,SpecRxn:=paste0(Species, "_", KO)]
network[,SpecRxn:=paste0(OTU, "_", KO)]
comp_order = mets_melt[,sd(value)/mean(value), by=compound][order(V1, decreasing = T), compound]
model_dat = data.table(compound = comp_order)
resid_dat = data.table(expand.grid(compound = comp_order, Sample = species_cmps[,unique(Sample)]))
all_rxns_removed = data.table()

#Go through in order of largest to smallest coefficient of variation
for(j in 1:length(comp_order)){
  met1 = comp_order[j]
  print(met1)
  cmp_dat = species_cmps[compound==met1]
  met_net = network[Reac==met1|Prod==met1]
  
  if(nrow(cmp_dat)==0){ #If not actually in network (environmental metabolite or whatever)
    next
  }
  uniq_rxns = cmp_dat[,unique(SpecRxn)]
  orig_scaling_mod = fit_single_scaling_mod(met1, cmp_dat, mets_melt)
  new_scaling_mod = orig_scaling_mod
  new_rsq = 1
  rsqs = rep(1, length(uniq_rxns))
  rxns_to_keep = c()
  rxns_to_remove = c()
  min_rxns = 3
  while(any(!uniq_rxns %in% rxns_to_keep) & length(uniq_rxns) > min_rxns & new_rsq > 1.1*summary(orig_scaling_mod)$r.squared){
    #If new one isn't much better htan old one we're done
    orig_scaling_mod = new_scaling_mod
    #update full model
    #print(uniq_rxns)
    
    all_scaling_mods = list()
    for(k in 1:length(uniq_rxns)){
      cmp_dat1 = cmp_dat[SpecRxn != uniq_rxns[k]]
      all_scaling_mods[[k]] = try(fit_single_scaling_mod(met1, cmp_dat1, mets_melt))
    }
    slopes = sapply(all_scaling_mods, function(x){ 
      if(class(x) != "try-error"){
        return(x$coefficients[2])
      } else return(0)
    }) #Rule out changes that make this negative
    rxns_to_keep = c(rxns_to_keep, uniq_rxns[slopes <= 0]) #Keep in rxns taht would break it if removed
    rsqs = sapply(all_scaling_mods, function(x){ return(summary(x)$r.squared)}) #Really?
    rxns_to_remove = c(rxns_to_remove, sample(uniq_rxns[which(rsqs==max(rsqs))], 1)) #randomly pick 1 if multiple
    print(rxns_to_remove)
    cmp_dat = cmp_dat[!SpecRxn %in% rxns_to_remove]
    new_scaling_mod = fit_single_scaling_mod(met1, cmp_dat, mets_melt)
    #print(summary(new_scaling_mod))
    new_rsq = summary(new_scaling_mod)$r.squared
    #Update stuff for next loop
    met_net = met_net[!SpecRxn %in% rxns_to_remove]
    uniq_rxns = uniq_rxns[!uniq_rxns %in% c(rxns_to_remove, rxns_to_keep)]
    
  }
  #Do we recalculate CMPs now? I guess so?
  if(length(rxns_to_remove) > 0){
    if(any(comp_order[(j+1):length(comp_order)] %in% c(network[SpecRxn %in% rxns_to_remove,unique(Reac)],network[SpecRxn %in% rxns_to_remove,unique(Prod)]))){
      #Recalculate only if necessary
      species_cmps = get_species_cmp_scores(species, network[!SpecRxn %in% rxns_to_remove], normalize = F, manual_agora = T, leave_rxns = T)[compound %in% mets_melt[,unique(compound)]]
      species_cmps[,SpecRxn:=paste0(Species, "_", KO)]
    }
    #Remove from full network
    network = network[!SpecRxn %in% rxns_to_remove] 
    all_rxns_removed = rbind(all_rxns_removed, data.table(compound = met1, Rxn = rxns_to_remove))
  }
  print(all_rxns_removed)
  #Model data
  scaling_coefs = coef(orig_scaling_mod)
  scaling_resids = resid(orig_scaling_mod)
  model_dat[compound==met1 ,Intercept:=scaling_coefs[1]]
  model_dat[compound==met1, Slope:=scaling_coefs[2]]
  model_dat[compound==met1, Rsq:=summary(orig_scaling_mod)$r.squared]
  if(length(scaling_resids) != nrow(resid_dat[compound==met1])) stop("Missing residuals")
  resid_dat[compound==met1, Resid:=scaling_resids]
  
}
new_cmps = get_species_cmp_scores(species, network, normalize = F, manual_agora = T)
new_cmps = add_residuals(new_cmps, model_dat, resid_dat)
var_shares = calculate_var_shares(new_cmps)

spec_table_wide = dcast(new_cmps, Sample+compound~Species, value.var = "newValue", fill = 0, fun.aggregate=sum)
spec_list = new_cmps[,unique(Species)]
var_shares = rbindlist(lapply(spec_list, function(y){
  all1 = rbindlist(lapply(spec_list, function(x){
    foo = spec_table_wide[,cov(get(x), get(y), use="complete.obs"), by=compound]
    foo[,Species:=x]
    return(foo)
  }))
  all1[,Species2:=y]
}))
var_shares = var_shares[,sum(V1),by=list(compound, Species)]
tot_vals = new_cmps[,sum(newValue), by = list(compound, Sample)]
true_met_var = tot_vals[,list(var(V1), mean(V1)), by = compound]
setnames(true_met_var, c("V1", "V2"), c("Var", "Mean"))
var_shares = merge(var_shares, true_met_var, by="compound")
var_shares[,VarShare:=V1/Var]
#Ok so why did this fail?
#are more variable mets predicted better? seems like an issue
setkey(model_dat, compound)
model_dat[comp_order]
model_dat[,CompRank:=match(compound, comp_order)]
model_dat[,cor.test(Rsq, CompRank, method = "spearman")] ## #Ooops

#Should have this be customizable to just do it for specific metabolites

orig_cmps = get_species_cmp_scores(species, network_results[[1]], normalize = F, manual_agora = T)
comp_fit = fit_cmp_mods(orig_cmps, mets_melt)
orig_cmps = add_residuals(orig_cmps, comp_fit[[1]], comp_fit[[2]])
var_shares_orig = calculate_var_shares(orig_cmps)
comp_fit[[1]][,CompRank:=match(compound, comp_order)]
comp_fit[[1]][,cor.test(Rsq, CompRank, method = "spearman")] ## Ok thsi was a strong effect before
model_dat_comp = merge(model_dat, comp_fit[[1]], by="compound")
ggplot(model_dat_comp, aes(x=Rsq.y, y = Rsq.x)) + geom_point() 
model_dat_comp[Rsq.x < Rsq.y & Rsq.y > 0.5] #Ok nothing really gets worse for 10 species - less true with HMP

model_dat_comp[,summary(CompRank.x), by=Rsq.x < Rsq.y] #Ok nothing really gets worse for 10 species - less true with HMP
#Lower ranks for these, aka higher CVs?? How/why? hmm
#Ok less ideal. and weird.
model_dat_comp[,wilcox.test(Rsq.x, Rsq.y)]
model_dat_comp[,list(median(Rsq.x), median(Rsq.y))]

contribs[,Species2:=ifelse(Species=="Inflow", "Residual", Species)]
compare_results = merge(contribs, var_shares_orig, by.x = c("Species2", "compound"), by.y = c("Species", "compound"), all = T)
compare_results = merge(compare_results, var_shares, by.x = c("Species2", "compound"), by.y = c("Species", "compound"),all = T)
compare_results[!is.na(V1.x),PosVarShare:=ifelse(V1.x > 0, V1.x/sum(V1.x[V1.x > 0], na.rm=T),0), by=compound]
compare_results[!is.na(V1.x),VarShareMagnitude:=abs(V1.x)/sum(abs(V1.x)), by=compound]
#true = varshare.x
#VarShareOrig = VarShare.y
#VarShareNew = VarShare
setnames(compare_results, c("VarShare.x", "VarShare.y", "VarShare"), c("TrueVarShare", "OrigVarShare", "NewVarShare"))
ggplot(compare_results[Species=="Inflow"], aes(x=PosVarShare, y = NewVarShare)) + geom_point() + xlim(-0.1, 1.5)
ggplot(compare_results[Species=="Inflow"], aes(x=PosVarShare, y = OrigVarShare)) + geom_point()+ xlim(-0.1, 1.5)
ggplot(compare_results[Species=="Inflow"], aes(x=VarShare.x, y = VarShare.y)) + geom_point()+ xlim(-0.1, 1.5)

library(GGally)
ggscatmat(compare_results[Species=="Inflow",list(PosVarShare, OrigVarShare, NewVarShare)], corMethod = "spearman")
ggscatmat(compare_results[Species=="Inflow",list(TrueVarShare, OrigVarShare, NewVarShare)], corMethod = "spearman")
compare_results[Species=="Inflow", table(PosVarShare > 0.1)]
compare_results[Species=="Inflow", table(PosVarShare > 0.1, OrigVarShare > 0.3)] #Models are just not very good
compare_results[Species=="Inflow", table(PosVarShare > 0.1, NewVarShare > 0.3)]
compare_results[Species=="Inflow", cor.test(PosVarShare, NewVarShare, method = "spearman", use = "complete.obs")]
compare_results[Species=="Inflow", cor.test(PosVarShare, OrigVarShare, method = "spearman", use = "complete.obs")]

##Bad Models: Lots of cases with estimated large microbial effect but not true
compare_results[Species=="Inflow" & PosVarShare > 0.1, cor.test(PosVarShare, OrigVarShare, method = "spearman", use = "complete.obs")]
compare_results[Species=="Inflow" & PosVarShare > 0.1, cor.test(PosVarShare, NewVarShare, method = "spearman", use = "complete.obs")]
compare_results[Species=="Inflow" & PosVarShare > 0.05, cor.test(PosVarShare, OrigVarShare, method = "spearman", use = "complete.obs")]
compare_results[Species=="Inflow" & PosVarShare > 0.05, cor.test(PosVarShare, NewVarShare, method = "spearman", use = "complete.obs")]
compare_results[Species=="Inflow" & PosVarShare < 0.05 & NewVarShare > 0.2]

ggscatmat(compare_results[Species!="Inflow",list(TrueVarShare, OrigVarShare, NewVarShare)], corMethod = "spearman")
ggscatmat(compare_results[Species!="Inflow",list(PosVarShare, OrigVarShare, NewVarShare)], corMethod = "spearman")
#compare_results[abs(VarShare.y) > 100]
#mets_melt[,list(var(value), sd(value)/mean(value)), by=compound]
#ggscatmat(compare_results[Species!="Inflow",list(PosVarShare, OrigVarShare, NewVarShare)], corMethod = "spearman")
ggscatmat(compare_results[Species!="Inflow",list(VarShareMagnitude, abs(OrigVarShare), abs(NewVarShare))], corMethod = "spearman")
compare_results[Species != "Inflow",table(PosVarShare > 0.1, OrigVarShare > 0.1)]
compare_results[Species != "Inflow",table(PosVarShare > 0.1, NewVarShare > 0.1)]
compare_results[Species != "Inflow" & PosVarShare > 0.1, table(NewVarShare > 0.1)]
compare_results[Species != "Inflow" & PosVarShare > 0.1, table(OrigVarShare > 0.1)]


consist_mets_new = model_dat_comp[Rsq.x > 0.2, compound]
consist_mets_orig = model_dat_comp[Rsq.y > 0.2, compound] #We get more consistnet mets, that's good
ggscatmat(compare_results[Species!="Inflow" & compound %in% consist_mets_orig,list(PosVarShare, OrigVarShare, NewVarShare)], corMethod = "spearman")
ggscatmat(compare_results[Species!="Inflow" & compound %in% consist_mets_new,list(PosVarShare, OrigVarShare, NewVarShare)], corMethod = "spearman")
ggscatmat(compare_results[Species!="Inflow" & compound %in% consist_mets_orig,list(VarShareMagnitude, abs(OrigVarShare), abs(NewVarShare))], corMethod = "spearman")
ggscatmat(compare_results[Species!="Inflow" & compound %in% consist_mets_new,list(VarShareMagnitude, abs(OrigVarShare), abs(NewVarShare))], corMethod = "spearman")
compare_results[Species != "Inflow" & compound %in% consist_mets_new, table(VarShareMagnitude > 0.2, OrigVarShare > 0.1)]
compare_results[Species != "Inflow" & compound %in% consist_mets_new, table(VarShareMagnitude > 0.2, NewVarShare > 0.1)]
compare_results[Species != "Inflow" & compound %in% consist_mets_orig, table(VarShareMagnitude > 0.2, OrigVarShare > 0.1)]
compare_results[Species != "Inflow" & compound %in% consist_mets_orig, table(VarShareMagnitude > 0.2, NewVarShare > 0.1)]
library(pROC)

compare_results[Species != "Inflow" & compound %in% consist_mets_new, roc(response = VarShareMagnitude > 0.2, scale(abs(NewVarShare)))]
compare_results[Species != "Inflow" & compound %in% consist_mets_new, roc(response = VarShareMagnitude > 0.2, scale(abs(OrigVarShare)))]
compare_results[Species != "Inflow" & compound %in% consist_mets_new, roc(response = PosVarShare > 0.1, scale(NewVarShare))]
compare_results[Species != "Inflow" & compound %in% consist_mets_new, roc(response = PosVarShare > 0.1, scale(OrigVarShare))]

compare_results[Species != "Inflow" & compound %in% consist_mets_orig, roc(response = PosVarShare > 0.1, scale(NewVarShare))]
compare_results[Species != "Inflow" & compound %in% consist_mets_orig, roc(response = PosVarShare > 0.1, scale(OrigVarShare))] #Wow


compare_results[Species != "Inflow" & compound %in% consist_mets_new & PosVarShare > 0.1 & OrigVarShare > 0.1 & NewVarShare < 0.1]
#Wait this is only a few?
#I'm so confused by this

#Reduce species
plot_summary_contributions(var_shares[compound %in% consist_mets_new & VarShare > 0.05], include_zeros = F, met_id_col = "compound")

#Consistent_mets - this is a good result (for 10-spec data)
compare_results[Species == "Inflow", table(PosVarShare < 0.9, NewVarShare < 0.9)]
compare_results[Species == "Inflow", table(PosVarShare < 0.9, OrigVarShare < 0.9)]

compare_results[VarShareMagnitude > 0.2, table(abs(OrigVarShare) > 0.2)]
compare_results[VarShareMagnitude > 0.2, table(abs(NewVarShare) > 0.2)]
#HMMM
# compare_results[VarShareMagnitude > 0.1, table(abs(OrigVarShare) > 0.1)]
# compare_results[VarShareMagnitude > 0.1, table(abs(NewVarShare) > 0.1)]
#Implement adding this to comparison

config_table = data.table(V1 = c("database", "genomeChoices","metType", "kegg_prefix", "data_prefix", "vsearch_path"), 
                          V2 = c("Greengenes 13_5 or 13_8", "AGORA genomes and models", "KEGG Compound IDs", "data/KEGGfiles/", "data/", "bin/vsearch"))
config_table = rbind(config_table, data.table(V1 = "manualAGORA", V2 = T))
config_table = rbind(config_table, data.table(V1 = c("revRxns", "rxnEdit"), V2 = T))

results_adj = run_mimosa2(config_table, species = species, mets = mets)
results_adj$modelData[,hist(Rsq, breaks = 30)]
results_adj$modelData[,table(Rsq > 0.2)]

compare_results[Species == "Inflow", table(TrueVarShare > 0.8)] #Ok
