#Run MIMOSA on HMP simulation data
#skip mapping step? 1st
#process arguments
library(data.table)
datadir = "~/Google Drive File Stream/My Drive/Genome_Sciences/MetabolSpeciesContribs/microbiome-metabolome-evaluation/"
species = fread(paste0(datadir, "HMP_test577_HMP_noise_specAbunds.txt"))
mets = fread(paste0(datadir, "HMP_test577_HMP_noise_metAbunds.txt"))
setnames(species, "Species", "OTU")

configTable = data.table(V1 = c("data_prefix", "database"), V2 = c("data/", "AGORA_manual"))
network_results = build_metabolic_model(species, configTable, manual_agora = T) #, input_data$netAdd) #input_data$geneAdd, 
network = network_results[[1]]
species = network_results[[2]] #Species should be unchanged in this case
rm(network_results)
species = dcast(species, OTU~Sample, value.var="value")
indiv_cmps = get_species_cmp_scores(species, network, manual_agora = T)

setnames(mets, 'medium', "compound")
mets[,compound:=gsub("[env]", "[e]", compound, fixed = T)]

#mets_melt = melt(mets, id.var = "compound", variable.name = "Sample", value.var = "value")
cmp_mods = fit_cmp_mods(indiv_cmps, mets[,list(Sample, compound, value)])
indiv_cmps = add_residuals(indiv_cmps, cmp_mods[[1]], cmp_mods[[2]])

var_shares = calculate_var_shares(indiv_cmps)
#shinyjs::logjs(devtools::session_info())
#Order dataset for plotting

met_fluxes = fread(paste0(datadir, "HMP_test577_HMP_noise_metFluxesFinal.txt"))
source(paste0(datadir, "FBA_functions.R"))
shap_contribs = getContributions(met_fluxes)
rm(met_fluxes)
shap_contribs[,PosVarShare:=ifelse(V1 > 0, V1/sum(V1[V1 > 0]),0), by=compound]
shap_contribs[,VarShareMagnitude:=abs(V1)/sum(abs(V1)), by=compound]
contrib_threshold_pos = 0.1
contrib_threshold_mag = 0.2
shap_contribs[,BinaryContribPos:=ifelse(PosVarShare > contrib_threshold_pos, 1, 0)]
shap_contribs[,BinaryContribMag:=ifelse(VarShareMagnitude > contrib_threshold_mag, 1, 0)]
shap_contribs[,compound:=gsub("[env]", "[e]", compound, fixed = T)]

var_shares[Species=="Residual", Species:="Inflow"]
var_shares_compare = merge(var_shares, shap_contribs, all = T, by = c("compound", "Species"))
var_shares_compare = var_shares_compare[compound %in% mets[,unique(compound)]]

met_summary = mets[,list(mean(value,na.rm=T), var(value,na.rm=T), sd(value, na.rm=T)/mean(value, na.rm=T)),by=compound]
setnames(met_summary, c("V1", "V2", "V3"), c("Mean", "Variance", "CoefVar"))

var_mets = met_summary[CoefVar > 3, compound]
ggplot(var_shares_compare[compound %in% var_mets & Species != "Inflow" & Species != "Residual"], aes(x=VarShare.x, y = VarShare.y)) + geom_point() + facet_wrap(~compound, scales = "free")

good_mets = var_shares_compare[!is.na(VarShare.x), unique(compound)]
good_mets = good_mets[good_mets %in% shap_contribs[,unique(compound)]]
var_shares_compare = var_shares_compare[compound %in% good_mets]
var_shares_compare[Species != "Inflow",table(BinaryContribPos, useNA = "ifany")]

ggplot(var_shares_compare[Species != "Inflow" & !is.na(BinaryContribPos)], aes(x=VarShare.x, color = factor(BinaryContribPos))) + 
  geom_density() + xlim(-2,2)
ggplot(var_shares_compare[!is.na(BinaryContribPos)], aes(x=VarShare.x, color = factor(BinaryContribPos))) + 
  geom_density() + xlim(-2,2)


var_shares_compare[,PosVarShare.x:=ifelse(V1.x > 0, V1.x/sum(V1.x[V1.x > 0], na.rm=T),0), by=compound]

var_shares_compare[,table(PosVarShare.x > 0.1, BinaryContribPos)] #Ok
library(pROC)
roc(var_shares_compare[,BinaryContribPos], var_shares_compare[,PosVarShare.x]) #eh, correlations are better?
roc(var_shares_compare[!Species %in% c("Inflow", "Residual"),BinaryContribPos], var_shares_compare[!Species %in% c("Inflow", "Residual"),PosVarShare.x]) #eh, correlations are better?

#Ok well if you do this we're good
consistent_mets = var_shares_compare[Species=='Inflow' & VarShare.x < 0.75, unique(compound)]
var_shares_compare[,length(unique(compound)), by=compound %in% consistent_mets]
roc(var_shares_compare[compound %in% consistent_mets,BinaryContribPos], var_shares_compare[compound %in% consistent_mets,PosVarShare.x]) #eh, correlations are better?

consistent_mets = var_shares_compare[Species=='Inflow' & VarShare.x < 0.85, unique(compound)]
var_shares_compare[,length(unique(compound)), by=compound %in% consistent_mets]
roc(var_shares_compare[compound %in% consistent_mets,BinaryContribPos], var_shares_compare[compound %in% consistent_mets,PosVarShare.x]) #eh, correlations are better?

consistent_mets = var_shares_compare[Species=='Inflow' & VarShare.x < 0.75, unique(compound)]
roc(var_shares_compare[compound %in% consistent_mets,BinaryContribPos], var_shares_compare[compound %in% consistent_mets,PosVarShare.x]) #eh, correlations are better?

met_level_outcomes = var_shares_compare[,list(sum(PosVarShare.x > 0.1 & BinaryContribPos==1)/sum(BinaryContribPos==1), 
                                              sum(PosVarShare.x < 0.1 & BinaryContribPos==0)/sum(BinaryContribPos==0),
                                              sum(PosVarShare.x > 0.1 & BinaryContribPos==1)/sum(PosVarShare.x > 0.1),
                                              sum((PosVarShare.x > 0.1 & BinaryContribPos==1)|(PosVarShare.x < 0.1 & BinaryContribPos==0))/length(BinaryContribPos),
                                              sum(BinaryContribPos),
                                              PosVarShare[Species=="Inflow"],
                                              PosVarShare.x[Species=="Inflow"]
                                              ),
                                        by = compound]
setnames(met_level_outcomes, c("compound", "Sens", "Spec", "PPV", "Accuracy", "NumContribs", "EnvContrib", "EstEnvContrib"))
env_contribs = ggplot(met_level_outcomes, aes(x=EnvContrib, y = EstEnvContrib)) + geom_point() + xlab("Environmental contribution") +
  ylab("Residual contribution")
save_plot(env_contribs, file = "results/HMP_mimosa_env_contrib.png")

ggplot(met_level_outcomes[EstEnvContrib < 0.6], aes(x=PPV)) + geom_histogram() #OK
ggplot(met_level_outcomes[EstEnvContrib < 0.6], aes(x=Sens)) + geom_histogram() #OK
ggplot(met_level_outcomes[EstEnvContrib < 0.6], aes(x=Spec)) + geom_histogram() #OK

ggplot(met_level_outcomes, aes(x=PPV)) + geom_histogram() #OK
ggplot(met_level_outcomes, aes(x=Sens)) + geom_histogram() #OK
ggplot(met_level_outcomes, aes(x=Spec)) + geom_histogram() #OK


ggplot(var_shares_compare, aes(x=VarShare.y, y = VarShare.x)) + geom_point(alpha = 0.7) + xlim(-2, 2) + ylim(-2,2) 

near_1_correct = var_shares_compare[abs(VarShare.x-1) < 0.1 & abs(VarShare.y - 1) < 0.1]

var_shares_compare_plot = ggplot(var_shares_compare[Species != "Inflow"], aes(x=VarShare.y, y = VarShare.x)) + geom_point(alpha = 0.7) + xlim(-2, 2) + ylim(-2,2)+ facet_wrap(~ifelse(compound %in% consistent_mets, "Consistent metabolites (Rsq > 0.25)", "Other")) +
  xlab("True relative contribution") + ylab("Estimated relative contribution")
save_plot(var_shares_compare_plot, file = "results/HMP_true_est_contribs.png", base_width = 7)


bad_mets = var_shares[]

#Output plots!
### COmpare with correlations
spec_met_corrs = get_correlation_contrib_comparison(species, fake_mets_melt = mets, shap_contribs, outcome_var = "BinaryContribPos")

