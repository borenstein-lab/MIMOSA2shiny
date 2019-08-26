### Compare more limited versions of mimosa 1 and 2 in more detail
#4/18/2019

library(data.table)
library(mimosa)
library(ggplot2)
library(cowplot)

source("scripts/mimosa2_dev_functions.R")
contrib_comparison = function(contribs1, contribs2){
  
}

bv_datasets = process_abundances("data/testData/mimosa1data/Dataset2_otu_table.txt", "data/testData/mimosa1data/Dataset2_mets.txt")

config_table = data.table(V1 = c("database", "genomeChoices","metType", "kegg_prefix", "data_prefix", "vsearch_path", "metagenome_format"), 
                          V2 = c("Greengenes 13_5 or 13_8", "PICRUSt KO genomes and KEGG metabolic model", "KEGG Compound IDs", "data/KEGGfiles/", "data/", "bin/vsearch", F))
bv2_kegg = run_mimosa2(config_table, species = bv_datasets[[1]], mets = bv_datasets[[2]])
bv1_kegg = run_mimosa1(bv_datasets[[1]], bv_datasets[[2]])
bv_compare_dat = compare_mimosa12(bv1_kegg[[5]], bv2_kegg[[2]])
bv_comparison = compare_mimosa12_plot(bv_compare_dat)
bv_met_plot = cmp_met_compare(bv_datasets, config_table)
bv_met_plot2 = cmp_met_compare(bv_datasets, config_table, compare_dat = bv_compare_dat)
bv_met_plot2a = cmp_met_compare(bv_datasets, config_table, compare_dat = bv_compare_dat, rank_based = T)

config_table = rbind(config_table, data.table(V1 = "rankBased", V2 = T))
bv2_kegg_rank = run_mimosa2(config_table, species = bv_datasets[[1]], mets = bv_datasets[[2]])

bv_compare_rank = compare_mimosa12(bv2_kegg[[2]], bv2_kegg_rank[[2]])
bv_rank_plot = compare_mimosa12_plot(bv_compare_rank)

asd_datasets = process_abundances("data/testData/miceASD_otus.txt", "data/testData/metsNMR_kegg.txt")
config_table = data.table(V1 = c("database", "genomeChoices","metType", "kegg_prefix", "data_prefix", "vsearch_path", "metagenome_format"), 
                          V2 = c("Greengenes 13_5 or 13_8", "PICRUSt KO genomes and KEGG metabolic model", "KEGG Compound IDs", "data/KEGGfiles/", "data/", "bin/vsearch", F))
asd2_kegg = run_mimosa2(config_table, species = asd_datasets[[1]], mets = asd_datasets[[2]])
asd1_kegg = run_mimosa1(asd_datasets[[1]], asd_datasets[[2]])
asd_compare_dat = compare_mimosa12(asd1_kegg, asd2_kegg)
asd_comparison = compare_mimosa12_plot(asd_compare_dat)
asd_met_plot = cmp_met_compare(asd_datasets, config_table)
asd_met_plot2 = cmp_met_compare(asd_datasets, config_table, compare_dat = asd_compare_dat)
asd_met_plot2a = cmp_met_compare(asd_datasets, config_table, compare_dat = asd_compare_dat, rank_based = T)

config_table = rbind(config_table, data.table(V1 = "rankBased", V2 = T))
asd2_kegg_rank = run_mimosa2(config_table, species = asd_datasets[[1]], mets = asd_datasets[[2]])

asd_compare_rank = compare_mimosa12(asd2_kegg[[2]], asd2_kegg_rank[[2]])
asd_rank_plot = compare_mimosa12_plot(asd_compare_rank)

abx_datasets = process_abundances("data/testData/mimosa1data/mice_otu_table.txt", "data/testData/mimosa1data/mice_mets_good.txt")
config_table = data.table(V1 = c("database", "genomeChoices","metType", "kegg_prefix", "data_prefix", "vsearch_path", "metagenome_format"), 
                          V2 = c("Greengenes 13_5 or 13_8", "PICRUSt KO genomes and KEGG metabolic model", "KEGG Compound IDs", "data/KEGGfiles/", "data/", "bin/vsearch", F))
abx2_kegg = run_mimosa2(config_table, species = abx_datasets[[1]], mets = abx_datasets[[2]])
abx1_kegg = run_mimosa1(abx_datasets[[1]], abx_datasets[[2]])
abx_compare_dat = compare_mimosa12(abx1_kegg, abx2_kegg)
abx_comparison = compare_mimosa12_plot(abx_compare_dat)
abx_met_plot = cmp_met_compare(abx_datasets, config_table)
abx_met_plot2 = cmp_met_compare(abx_datasets, config_table, compare_dat = abx_compare_dat)
abx_met_plot2a = cmp_met_compare(abx_datasets, config_table, compare_dat = abx_compare_dat, rank_based = T)

config_table = rbind(config_table, data.table(V1 = "rankBased", V2 = T))
abx2_kegg_rank = run_mimosa2(config_table, species = abx_datasets[[1]], mets = abx_datasets[[2]])

abx_compare_rank = compare_mimosa12(abx2_kegg[[2]], abx2_kegg_rank[[2]])
abx_rank_plot = compare_mimosa12_plot(abx_compare_rank)

sim_data = process_abundances("data/testData/sim_data/HMP_test577_HMP_noise_specAbunds.txt", "data/testData/sim_data/HMP_test577_HMP_noise_metAbunds.txt", 
                              fluxes_file = "data/testData/sim_data/HMP_test577_HMP_noise_metFluxesFinal.txt", simulated = T)
config_table = data.table(V1 = c("database", "genomeChoices","metType", "kegg_prefix", "data_prefix", "vsearch_path", "revRxns"), 
                          V2 = c("Greengenes 13_5 or 13_8", "AGORA genomes and models", "KEGG Compound IDs", "data/KEGGfiles/", "data/", "bin/vsearch", T))
config_table = rbind(config_table, data.table(V1 = "manualAGORA", V2 = T))
config_table = config_table[!V1 %in% c("rxnEdit", "revRxns")]

hmp2 = run_mimosa2(config_table, species = sim_data[[1]], mets = sim_data[[2]])
hmp1 = run_mimosa1(species = sim_data[[1]], mets = sim_data[[2]], simulated = T, config_table = config_table)

config_table = rbind(config_table, data.table(V1 = "rankBased", V2 = T))
hmp2_rank = run_mimosa2(config_table, species = sim_data[[1]], mets = sim_data[[2]])

hmp_compare_rank = compare_mimosa12(hmp2[[2]], hmp2_rank[[2]], contribs = sim_data[[3]])
hmp_rank_plot = compare_mimosa12_plot(hmp_compare_rank, simulated = T)
# Does rank id microbial metabolites better?
hmp_compare_rank[,cor(1-EnvVarShare, Rsq.x, use = "complete.obs", method = "spearman")]
hmp_compare_rank[,cor(1-EnvVarShare, Rsq.y, use = "complete.obs", method = "spearman")] 
hmp_compare_rank[!is.na(Rsq.x) & !is.na(EnvVarShare), sum((1 -EnvVarShare > 0.2) == (Rsq.x > 0.2))/length(Rsq.x)]
hmp_compare_rank[!is.na(Rsq.y) & !is.na(EnvVarShare), sum((1 -EnvVarShare > 0.2) == (Rsq.y > 0.2))/length(Rsq.y)] #Hm
hmp_compare_rank[!is.na(Rsq.x) & !is.na(EnvVarShare), sum((1 -EnvVarShare > 0.5) == (Rsq.x > 0.5))/length(Rsq.x)]
hmp_compare_rank[!is.na(Rsq.y) & !is.na(EnvVarShare), sum((1 -EnvVarShare > 0.5) == (Rsq.y > 0.5))/length(Rsq.y)] #Hm

hmp_compare_dat = compare_mimosa12(hmp1, hmp2, contribs = sim_data[[3]])
hmp_comparison = compare_mimosa12_plot(hmp_compare_dat, simulated = T)
hmp_met_plot = cmp_met_compare(sim_data, config_table, simulated = T)
hmp_met_plot2 = cmp_met_compare(sim_data, config_table, simulated = T, compare_dat = hmp_compare_dat)
hmp_met_plot2a = cmp_met_compare(sim_data, config_table, simulated = T, compare_dat = hmp_compare_dat, rank_based = T)


#Do another analysis with a higher level of env noise also
ten_spec_dat = process_abundances("data/testData/sim_data/allSpeciesEnv3.txt", "data/testData/sim_data/allMetabolitesEnv3.txt", fluxes_file = "data/testData/sim_data/allMetFluxesEnv3.txt", simulated = T)
config_table = data.table(V1 = c("database", "genomeChoices","metType", "kegg_prefix", "data_prefix", "vsearch_path", "revRxns"), 
                          V2 = c("Greengenes 13_5 or 13_8", "AGORA genomes and models", "KEGG Compound IDs", "data/KEGGfiles/", "data/", "bin/vsearch", T))
config_table = rbind(config_table, data.table(V1 = "manualAGORA", V2 = T))
config_table = config_table[!V1 %in% c("rxnEdit", "revRxns")]

ten_spec2 = run_mimosa2(config_table, species = ten_spec_dat[[1]], mets = ten_spec_dat[[2]])
ten_spec1 = run_mimosa1(species = ten_spec_dat[[1]], mets = ten_spec_dat[[2]], simulated = T, config_table = config_table)

config_table = rbind(config_table, data.table(V1 = "rankBased", V2 = T))
ten_spec2_rank = run_mimosa2(config_table, species = ten_spec_dat[[1]], mets = ten_spec_dat[[2]])

ten_spec_compare_rank = compare_mimosa12(ten_spec2[[2]], ten_spec2_rank[[2]], contribs = ten_spec_dat[[3]])
ten_spec_rank_plot = compare_mimosa12_plot(ten_spec_compare_rank, simulated = T)
#Does rank id true contrib mets more accurately?
ten_spec_compare_rank[,cor(1-EnvVarShare, Rsq.x, use = "complete.obs")]
ten_spec_compare_rank[,cor(1-EnvVarShare, Rsq.y, use = "complete.obs")] #Yes
ten_spec_compare_rank[!is.na(Rsq.x) & !is.na(EnvVarShare), sum((1 -EnvVarShare > 0.2) == (Rsq.x > 0.2))/length(Rsq.x)]
ten_spec_compare_rank[!is.na(Rsq.y) & !is.na(EnvVarShare), sum((1 -EnvVarShare > 0.2) == (Rsq.y > 0.2))/length(Rsq.y)] #Hm
ten_spec_compare_rank[!is.na(Rsq.x) & !is.na(EnvVarShare), sum((1 -EnvVarShare > 0.5) == (Rsq.x > 0.5))/length(Rsq.x)]
ten_spec_compare_rank[!is.na(Rsq.y) & !is.na(EnvVarShare), sum((1 -EnvVarShare > 0.5) == (Rsq.y > 0.5))/length(Rsq.y)] #Hm

ggplot(ten_spec_compare_rank, aes(x=factor(1-EnvVarShare > 0.2), y = Rsq.x)) + geom_boxplot()
ggplot(ten_spec_compare_rank, aes(x=factor(1-EnvVarShare > 0.2), y = Rsq.y)) + geom_boxplot()

ten_spec_compare_dat = compare_mimosa12(ten_spec1, ten_spec2, contribs = ten_spec_dat[[3]])
ten_spec_comparison = compare_mimosa12_plot(ten_spec_compare_dat, simulated = T)
ten_spec_met_plot = cmp_met_compare(ten_spec_dat, config_table, simulated = T)
ten_spec_met_plot2 = cmp_met_compare(ten_spec_dat, config_table, simulated = T, compare_dat = ten_spec_compare_dat)
ten_spec_met_plot2a = cmp_met_compare(ten_spec_dat, config_table, simulated = T, compare_dat = ten_spec_compare_dat, rank_based = T)

save_plot(asd_comparison, file = "results/cmpFixed/ASD_met_comparison.png", base_width = 9, base_height = 3.5)
save_plot(bv_comparison, file = "results/cmpFixed/BV_met_comparison.png", base_width = 9, base_height = 3.5)
save_plot(abx_comparison, file = "results/cmpFixed/Abx_met_comparison.png", base_width = 9, base_height = 3.5)
save_plot(hmp_comparison, file = "results/cmpFixed/HMP_sim_met_comparison.png", base_width = 9, base_height = 6)
save_plot(ten_spec_comparison, file = "results/cmpFixed/tenSpec_sim_met_comparison.png", base_width = 9, base_height = 6)

save_plot(asd_met_plot, file = "results/cmpFixed/ASD_cmp_met.png", base_width = 15, base_height = 10)
save_plot(abx_met_plot, file = "results/cmpFixed/Abx_cmp_met.png", base_width = 18, base_height = 12)
save_plot(bv_met_plot, file = "results/cmpFixed/BV_cmp_met.png", base_width = 15, base_height = 10)
save_plot(hmp_met_plot, file = "results/cmpFixed/HMP_cmp_met.png", base_width = 18, base_height = 12)
save_plot(ten_spec_met_plot, file = "results/cmpFixed/TenSpec_cmp_met.png", base_width = 15, base_height = 10)

save_plot(asd_met_plot2, file = "results/cmpFixed/ASD_cmp_metMimosa.png", base_width = 15, base_height = 10)
save_plot(abx_met_plot2, file = "results/cmpFixed/Abx_cmp_metMimosa.png", base_width = 18, base_height = 12)
save_plot(bv_met_plot2, file = "results/cmpFixed/BV_cmp_metMimosa.png", base_width = 15, base_height = 10)
save_plot(hmp_met_plot2, file = "results/cmpFixed/HMP_cmp_metMimosa.png", base_width = 20, base_height = 13)
save_plot(ten_spec_met_plot2, file = "results/cmpFixed/TenSpec_cmp_metMimosa.png", base_width = 15, base_height = 10)

save_plot(asd_met_plot2a, file = "results/cmpFixed/ASD_cmp_metMimosaRank.png", base_width = 15, base_height = 10)
save_plot(abx_met_plot2a, file = "results/cmpFixed/Abx_cmp_metMimosaRank.png", base_width = 18, base_height = 12)
save_plot(bv_met_plot2a, file = "results/cmpFixed/BV_cmp_metMimosaRank.png", base_width = 15, base_height = 10)
save_plot(hmp_met_plot2a, file = "results/cmpFixed/HMP_cmp_metMimosaRank.png", base_width = 20, base_height = 13)
save_plot(ten_spec_met_plot2a, file = "results/cmpFixed/TenSpec_cmp_metMimosaRank.png", base_width = 15, base_height = 10)

save_plot(asd_rank_plot, file = "results/cmpFixed/ASD_rankReg_results.png", base_width = 9, base_height = 3.5)
save_plot(bv_rank_plot, file = "results/cmpFixed/BV_rankReg_results.png", base_width = 9, base_height = 3.5)
save_plot(abx_rank_plot, file = "results/cmpFixed/Abx_rankReg_results.png", base_width = 9, base_height = 3.5)
save_plot(hmp_rank_plot, file = "results/cmpFixed/HMP_rankReg_results.png", base_width = 9, base_height = 6)
save_plot(ten_spec_rank_plot, file = "results/cmpFixed/tenSpec_rankReg_results.png", base_width = 9, base_height = 6)


##Making sure I understand contributions with rank regression
rsq_compare = bv2_kegg_rank[[1]][Species=="Residual"]
rsq_compare[,ContribRsq:=1-VarShare]
rsq_compare[,hist(Rsq-ContribRsq, breaks = 20)]
rsq_compare[Rsq-ContribRsq > 0.2]

## If we ran each species' cmp scores on their own, what would we get? Different scaling.
# Can we just scale var explained within a different residual contribution?
#Bit of a philosophical question I think
# But still where does the Rfit Rsq come from?
#Residual sum of squared Jaeckel dispersion? - read this
library(Rfit)
?rfit

#Just to double-check they are equal for normal regression
bv2_kegg[[1]][,ContribRsq:=1-VarShare] 
bv2_kegg[[1]][Species=="Residual", hist(Rsq-ContribRsq, breaks = 20)]

### example divergent metabolites
bv_compare_dat = compare_mimosa12(bv1_kegg, bv2_kegg)
bv_compare_dat[PosNegCompare == "TRUE_FALSE" | PosNegCompare=="FALSE_TRUE"]
bv_config_table = data.table(V1 = c("database", "genomeChoices","metType", "kegg_prefix", "data_prefix", "vsearch_path", "metagenome_format"), 
                            V2 = c("Greengenes 13_5 or 13_8", "PICRUSt KO genomes and KEGG metabolic model", "KEGG Compound IDs", "data/KEGGfiles/", "data/", "bin/vsearch", F))

bv_network = build_metabolic_model(bv_datasets[[1]], bv_config_table)
bv_spec_cmp_scores = get_species_cmp_scores(species_table = bv_datasets[[1]], bv_network[[1]], normalize = T)
bv_cmp_scores = bv_spec_cmp_scores[,sum(CMP), by=list(Sample, compound)]


bv_cmp_scores1 = melt(bv1_kegg[[6]], id.var = "compound", variable.name = "Sample")
bv_contrib_table = generate_contribution_table_using_picrust(bv_datasets[[1]], "data/picrustGenomeData/16S_13_5_precalculated.tab.gz", "data/picrustGenomeData/indivGenomes/", "_genomic_content.tab", copynum_column = T)
setnames(bv_contrib_table, c("contribution", "copy_number"), c("CountContributedByOTU", "GeneCountPerGenome"))
bv_contrib_table = bv_contrib_table[CountContributedByOTU != 0]
all_otus = sort(unique(bv_contrib_table[,OTU]))
bv_contrib_table = single_spec_musicc(bv_contrib_table)
all_koAbunds_byOTU = contribs_by_species_list(bv_contrib_table, valueVar = "singleMusicc", "foo", write_out = F)
rxn_table = fread("data/KEGGfiles/full_rxn_table.txt")
ko_net = generate_genomic_network(bv_contrib_table[,unique(Gene)], keggSource = "KeggTemplate", normalize = T, rxn_table = rxn_table, degree_filter = 0)
bv1_networks_alone = lapply(1:length(all_koAbunds_byOTU), function(x){
  sub_ko_net = generate_genomic_network(all_koAbunds_byOTU[[x]][,KO], keggSource = "KeggTemplate", degree_filter = 0, rxn_table = rxn_table, normalize = T)
  sub_ko_net[[3]][,OTU:=all_otus[x]]
  return(sub_ko_net)
})
  
bv1_network = rbindlist(lapply(bv1_networks_alone, function(x){
  return(x[[3]])
}))
bv1_network[,stoichRatio:=stoichProd/stoichReac]
bv_network[[1]][,stoichRatio:=stoichProd/stoichReac]
bv_network_compare = merge(bv1_network, bv_network[[1]], by=c("OTU", "KO", "Reac", "Prod", "stoichRatio"), all = T)
bv_network_compare[is.na(normalized_copy_number)|is.na(stoichProd.x)|is.na(stoichProd.y)|is.na(stoichReac.x)|is.na(stoichReac.y)] #Hmm how did we get 116 more rows? 
#Ok, so the edge list networks are the same

x=1
x=18
x = 3
testNet = bv1_networks_alone[[x]][[1]][row.names(bv1_networks_alone[[x]][[1]]) %in% bv_datasets[[2]][,KEGG],]
testNet1 = testNet[,apply(testNet, 2, function(x){ any(x != 0)})]
testNet2 = bv1_networks_alone[[x]][[3]][Prod %in% bv_datasets[[2]][,KEGG]|Reac %in% bv_datasets[[2]][,KEGG]]
testNet_m2 = bv_network[[1]][OTU==all_otus[x] & (Prod %in% bv_datasets[[2]][,KEGG]|Reac %in% bv_datasets[[2]][,KEGG])]

testDat = get_cmp_scores(bv1_networks_alone[[x]][[1]][row.names(bv1_networks_alone[[x]][[1]]) %in% bv_datasets[[2]][,KEGG],], 
                         norm_kos = all_koAbunds_byOTU[[x]])
testDat = melt(testDat, id.var = "compound", variable.name = "Sample")

testDat2 = get_species_cmp_scores(bv_datasets[[1]], network = testNet_m2)# Hmm
testDat2a = bv_spec_cmp_scores[Species==all_otus[x]]
testDat2 = merge(testDat2, testDat2a, by = c("Species", "Sample", "compound", "value"), all = T)[compound %in% bv_datasets[[2]][,KEGG]]
testDat2[abs(CMP.x - CMP.y) > 10e-12]
setnames(testDat2, "value", "RelAbund")
testDat2[is.na(CMP.x)|is.na(CMP.y)]

testDat = merge(testDat, testDat2, by=c("compound", "Sample"), all = T)
#testDat[,value:=value*100]
testDat[value==0 & CMP.x != 0]
testDat[value!=0 & CMP.x == 0] #OK I'm not worrying about this
testDat[value != 0 & CMP.x != 0, list(median(value/CMP.x), min(value/CMP.x), max(value/CMP.x)), by=compound][abs(V1-1) > 10e-10]
testDat[abs(value/CMP.x-1) > 10e-16 , unique(compound)]
testDat[abs(value/CMP.x-1) > 10e-16]
testDat[is.na(value)|is.na(CMP.x)]

testDat_sub = testDat[compound=="C00078"]
otu_abunds = melt(bv_datasets[[1]], variable.name = "Sample")
otu_abunds[,OTURelAbund:=value/sum(value), by=Sample]
setnames(otu_abunds, "value", "ReadCount")
testDat_sub = merge(testDat_sub, otu_abunds[OTU==227000], by="Sample", all = T)
testDat_sub[is.na(OTURelAbund), OTURelAbund:=0]
testDat_sub[,Factor1:=value/(OTURelAbund*100000)]
testDat_sub[,Factor2:=CMP.x/(OTURelAbund*100000)]


testDat_sub = testDat[compound=="C00037"]
otu_abunds = melt(bv_datasets[[1]], variable.name = "Sample")
otu_abunds[,OTURelAbund:=value/sum(value), by=Sample]
setnames(otu_abunds, "value", "ReadCount")
testDat_sub = merge(testDat_sub, otu_abunds[OTU==227000], by="Sample", all = T)
testDat_sub[is.na(OTURelAbund), OTURelAbund:=0]
testDat_sub[,Factor1:=value/(OTURelAbund*100000)]
testDat_sub[,Factor2:=CMP.x/(OTURelAbund*100000)]
testNet_m2[Reversible==0 & (Reac=="C00037"|Prod=="C00037")]

network = testNet_m2[Reversible==0]
network_reacs = network[,list(OTU, KO, Reac, stoichReac, normalized_copy_number)]
network_prods = network[,list(OTU, KO, Prod, stoichProd, normalized_copy_number)]
network_reacs[,stoichReac:=-1*stoichReac]
setnames(network_reacs, c("Reac", "stoichReac"), c("compound", "stoich"))
setnames(network_prods, c("Prod", "stoichProd"), c("compound", "stoich"))
#Remove multiple-encoded things
setkey(network_reacs, NULL)
setkey(network_prods, NULL)
network_reacs = unique(network_reacs)
network_prods = unique(network_prods)
network_reacs[,stoich:=as.double(stoich)]
network_prods[,stoich:=as.double(stoich)]
#network_prods = network_prods[compound %in% c("C00031", "C00078")]
#network_reacs = network_reacs[compound %in% c("C00031", "C00078")]
network_reacs[,stoich2:=stoich/abs(sum(stoich)), by=list(OTU, compound)]
network_prods[,stoich2:=stoich/sum(stoich), by=list(OTU, compound)]
net2 = rbind(network_reacs, network_prods, fill = T)
net2 = net2[!is.na(compound)] #remove NAs
net2[,stoich3:=stoich2*normalized_copy_number]
net2[,sum(stoich2), by=list(OTU, KO, compound)][abs(V1) < 10e-12]

netmat = melt(data.table(testNet[row.names(testNet) %in% bv_datasets[[2]][,KEGG],], compound = row.names(testNet[row.names(testNet) %in% bv_datasets[[2]][,KEGG],])), variable.name = "KO")
netmat = netmat[value != 0]
net2 = net2[compound %in% bv_datasets[[2]][,KEGG]]
net2_comp = merge(net2, netmat, by=c("KO", "compound"), all = T)
net2_comp[value != stoich, unique(compound)]
net2_comp[is.na(value), unique(compound)]
#net2_comp[is.na(stoich)]
bad_comps = net2_comp[,sum(value==stoich2, na.rm = T)/length(value), by=compound][V1 != 1, compound] 

net2_comp1 = net2_comp[compound %in% bad_comps]
#Ok, pick up here tomorrow and figure out how generate_genomic_network is adjusting the numbers
net2_comp1[order(compound, KO)]
#Something with duplicated KOs?

spec_cmps_kos = merge(otu_abunds[OTU==227000], net2, by = "OTU", allow.cartesian = T)
spec_cmps_kos[OTURelAbund != 0,table(abs(stoich/OTURelAbund))]
spec_cmps_kos2 = spec_cmps_kos[,sum(stoich3*OTURelAbund), by=list(OTU, Sample, compound, OTURelAbund)]
spec_cmps_kos2[OTURelAbund != 0, table(abs(V1/OTURelAbund))]
#Ok, try working through generate_genomic_network and figure out what's up
#Is it because this version accts for 16S and the other one doesn't? Can't be entirely because is constant by species

testDat_sub = testDat[compound=="C00031"]
otu_abunds = melt(bv_datasets[[1]], variable.name = "Sample")
otu_abunds[,OTURelAbund:=value/sum(value), by=Sample]
setnames(otu_abunds, "value", "ReadCount")
testDat_sub = merge(testDat_sub, otu_abunds[OTU==227000], by="Sample", all = T)
testDat_sub[is.na(OTURelAbund), OTURelAbund:=0]
testDat_sub[,Factor1:=value/(OTURelAbund*100000)]
testDat_sub[,Factor2:=CMP.x/(OTURelAbund*100000)]


##What is going on with these compounds????? C00022, C00031, C00078, 65, 49, 37
testNet_m2[Reversible==0 & (Reac=="C00078"|Prod=="C00078")]


testDat[]
test_kos = bv_contrib_table[OTU==all_otus[1]]
test_combine = merge(test_kos, testNet_m2, by.x = c("OTU", "Gene"), by.y = c("OTU", "KO"))
test_combine[abs(normalized_copy_number-GeneCountPerGenome) > 10e-16] #Ok so they are the same. that's good.

testDat[,cor(value, CMP.x), by=compound][V1 != 1]

testDat[compound=="C00031"]
testDat[compound=="C00022"]


bv1_networks_alone[[1]][[1]]["C00031",bv1_networks_alone[[1]][[1]]["C00031",] != 0]
sort(names(bv1_networks_alone[[1]][[1]]["C00031",bv1_networks_alone[[1]][[1]]["C00031",] != 0]))
"K00845" %in% names(bv1_networks_alone[[1]][[1]])
testNet2[Prod=="C00031"|Reac=="C00031"]
testNet2[Prod=="C00031"|Reac=="C00031", sort(unique(KO))]
testNet2[KO=="K00845"]

bv1_networks_alone[[1]][[1]]["C00022",bv1_networks_alone[[1]][[1]]["C00022",] != 0]
kos22 = sort(names(bv1_networks_alone[[1]][[1]]["C00022",bv1_networks_alone[[1]][[1]]["C00022",] != 0]))
testNet2[Prod=="C00022"|Reac=="C00022"]
testNet2[Prod=="C00022"|Reac=="C00022", sort(unique(KO))]
abs_kos = testNet2[Prod=="C00022"|Reac=="C00022", sort(unique(KO[!KO %in% kos22]))]
testNet_m2[KO %in% abs_kos & Reversible==0] #None of these have C00022???? 
testDat[compound=="C00022" & value != CMP.x]
ko_abunds = melt(all_koAbunds_byOTU[[1]][KO %in% testNet_m2[(Prod=="C00022"|Reac=="C00022") & Reversible==0, KO]], id.var = "KO", variable.name = "Sample")
bad_samps = testDat[compound=="C00022" & value != CMP.x, Sample]
ko_abunds[Sample %in% bad_samps]
#C00022  j31tr55

test_mat = bv1_networks_alone[[1]][[1]]
good_kos = names(test_mat[,apply(testNet, 2, function(x){ any(x != 0)})])
table(good_kos %in% network[Reversible==0, KO])
table(network[Reversible==0, unique(KO)] %in% good_kos) #Ok the same
test_mat = data.table(test_mat, compound = row.names(test_mat))

get_cmp_scores(test_mat)


cmps_alone = get_all_singleSpec_cmps(all_otus, all_koAbunds_byOTU, valueVar = "singleMusicc", out_prefix = "foo", rxn_table = ko_net[[3]], write_out = F, degree_filter = 0) #ko_net from results file
cmps_tot_alone = rbindlist(lapply(1:length(cmps_alone), function(x){
  newDat = data.table(cmps_alone[[x]], OTU = all_otus[x])
  return(newDat)
}))
cmps_tot_alone = melt(cmps_tot_alone, id.var = c("compound", "OTU"), variable.name = "Sample")
cmps_compare = merge(cmps_tot_alone, bv_spec_cmp_scores, by.x = c("OTU", "Sample", "compound"), by.y = c("Species", "Sample", "compound"), all = T)
cmps_compare = cmps_compare[compound %in% bv_datasets[[2]][,KEGG]] #narrow down
#cmps_compare[,value:=100*value]
cmps_compare[,cor(CMP, value.x), by=list(compound,OTU)][!is.na(V1)][order(V1)][1:100] #Ok this is closer than before
cmps_compare[,cor(CMP, value.x), by=list(compound,OTU)][,hist(V1, breaks = 40)]
cmps_compare[,hist(CMP-value.x, breaks = 40)]
cmps_compare[order(abs(CMP-value), decreasing = T)][1:100]
cmps_compare[CMP*value.x < 0] #One pos other neg - all v small numbers

###Ok we're going to walk through this with C00031 and 227000
which(all_otus==227000)
testNet = bv1_networks_alone[[18]][[1]][row.names(bv1_networks_alone[[18]][[1]])=="C00031",]
testNet1 = testNet[,apply(testNet, 2, function(x){ any(x != 0)})]
testNet2 = bv1_networks_alone[[18]][[3]][(Prod=="C00031"|Reac=="C00031")]
testNet_m2 = bv_network[[1]][OTU==all_otus[18] & (Prod=="C00031"|Reac=="C00031")]
testDat = get_cmp_scores(bv1_networks_alone[[18]][[1]], norm_kos = all_koAbunds_byOTU[[18]])
testDat = melt(testDat, id.var = "compound", variable.name = "Sample")
testDat2 = get_species_cmp_scores(bv_datasets[[1]], network = testNet_m2)# Hmm
testDat2a = bv_spec_cmp_scores[Species==all_otus[18]]
testDat2 = merge(testDat2, testDat2a, by = c("Species", "Sample", "compound"), all = T)[compound=="C00031"]
testDat2[CMP.x != CMP.y]

testDat = merge(testDat, testDat2, by=c("compound", "Sample"))
#testDat[,value:=value*100]
testDat[value==0 & CMP.x != 0]
testDat[value!=0 & CMP.x == 0]
testDat[value != 0 & CMP.x != 0, summary(value/CMP.x)]
testDat[abs(value/CMP.x-1) > 10e-16 & value > 10e-7, unique(compound)]



tot_cmps_compare = cmps_compare[,list(sum(value), sum(CMP)), by=list(Sample, compound)]
tot_cmps_compare[,cor(V1, V2), by=compound][,hist(V1, breaks = 40)] #yeow


bv_cmp_scores1a 

mets_melt = melt(bv_datasets[[2]], id.var = "KEGG", variable.name = "Sample")
bv_compare_mets = merge(bv_cmp_scores, mets_melt, by.x=c("Sample", "compound"), by.y = c("Sample", "KEGG"))

ggplot(bv_compare_mets[compound=="C00082"], aes(x=V1, y = value)) + geom_point()
ggplot(bv_compare_mets[compound=="C00123"], aes(x=V1, y = value)) + geom_point()
bad_mets = bv_compare_mets[,var(V1), by=compound][V1==0, compound]
ggplot(bv_compare_mets[!compound %in% bad_mets], aes(x=V1, y = value)) + geom_point() + geom_smooth() + facet_wrap(~compound, scales = "free")
#Separate out by synth only, deg only, both
cmp_met_plots_bv = ggplot(bv_compare_mets[!compound %in% bad_mets], aes(x=V1, y = log1p(value))) + geom_point() + geom_smooth() + facet_wrap(~compound, scales = "free")
##HMM

abx_mets_melt = melt(abx_datasets[[2]], id.var = "KEGG", variable.name = "Sample")
abx_network = build_metabolic_model(abx_datasets[[1]], bv_config_table)
abx_spec_cmp_scores = get_species_cmp_scores(species_table = abx_datasets[[1]], abx_network[[1]], normalize = T)
abx_cmp_scores = abx_spec_cmp_scores[,sum(CMP), by=list(Sample, compound)]
abx_cmp_scores1 = melt(abx1_kegg[[6]], id.var = "compound", variable.name = "Sample")

abx_compare_mets = merge(abx_cmp_scores, abx_mets_melt, by.x=c("Sample", "compound"), by.y = c("Sample", "KEGG"))

bad_mets = abx_compare_mets[,var(V1), by=compound][V1==0, compound]
ggplot(abx_compare_mets[!compound %in% bad_mets], aes(x=V1, y = value)) + geom_point() + geom_smooth() + facet_wrap(~compound, scales = "free")
#Separate out by synth only, deg only, both
cmp_met_plots_abx= ggplot(abx_compare_mets[!compound %in% bad_mets], aes(x=V1, y = log1p(value))) + geom_point() + geom_smooth() + facet_wrap(~compound, scales = "free")
##HMM
cmp_met_plots_abx2= ggplot(abx_compare_mets[!compound %in% bad_mets & grepl("^No", Sample)], aes(x=V1, y = log1p(value))) + geom_point() + geom_smooth() + facet_wrap(~compound, scales = "free")
cmp_met_plots_abx2a= ggplot(abx_compare_mets[!compound %in% bad_mets & grepl("^No", Sample)], aes(x=V1, y = value)) + geom_point() + geom_smooth() + facet_wrap(~compound, scales = "free")
#account for covariates/study desgin?
synth_deg = abx_compare_mets[,list(sum(V1 > 0), sum(V1 < 0)), by=compound]
setnames(synth_deg, c("V1", "V2"), c("SynthSamples", "DegSamples"))
synth_deg[,table(SynthSamples > 0, DegSamples > 0)]
abx_compare_mets = merge(abx_compare_mets, synth_deg, by="compound")
abx_compare_mets[,CMPType:=ifelse(SynthSamples > 0 & DegSamples==0, "Synth", "Both")]
abx_compare_mets[,CMPType:=ifelse(DegSamples > 0 & SynthSamples==0, "Deg", CMPType)]
ggplot(abx_compare_mets[!compound %in% bad_mets & CMPType=="Synth"], aes(x=V1, y = value)) + geom_point() + geom_smooth() + facet_wrap(~compound, scales = "free")
ggplot(abx_compare_mets[!compound %in% bad_mets & CMPType=="Both"], aes(x=V1, y = value)) + geom_point() + geom_smooth() + facet_wrap(~compound, scales = "free")


head(bv_compare_mets)
bv_compare_mets = merge(bv_compare_mets, bv_cmp_scores1, by = c("compound", "Sample"), all.x = T)
setnames(bv_compare_mets, c("V1", "value.x", "value.y"), c("CMP2", "Met", "CMP1"))
bv_compare_mets[,cor(CMP1, CMP2), by=compound]

#### Get mimosa1 single species to be able to compare for sure
# I feel like we've done this previously.
spec_rel_abunds = melt(bv_datasets[[1]], id.var = "OTU", variable.name = "Sample")
spec_rel_abunds[,RelAbund:=value/sum(value), by=Sample]# 0.0167597765
bv_contrib_table[GeneCountPerGenome==1 & OTU==1038074 & Sample=="p6z1tr12"] #0.02287166 #different
bv_contrib_table[,sum(GeneCountPerGenome==1), by=OTU] 
foo = unique(bv_contrib_table[GeneCountPerGenome==1, list(Sample, OTU, CountContributedByOTU)])
##Oh if there are species with no single-copy genes? #Yep a couple - definitely better ways to handle this

get_cmp_scores
get_species_cmp_scores 
bv_contrib_table = single_spec
all_koAbunds_byOTU = contribs_by_species_list(contribs, valueVar = valueVar, out_prefix, write_out)

