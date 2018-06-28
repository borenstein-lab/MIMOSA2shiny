#data from mimosaMainRun.rmd
library(data.table)
library(ggplot2)
library(mimosa)
library(RColorBrewer)
library(pROC)
datadir = "~/Documents/GS_server/PROJECTS/MORE_METABOLOMICS_CN/RESULTS/AGORA/"
homedir =  "~/Documents/GS_server/PROJECTS/MORE_METABOLOMICS_CN/"
datadir2 = "../../../MetabolSpeciesContribs/microbiome-metabolome-evaluation/"
source("../../../MetabolSpeciesContribs/microbiome-metabolome-evaluation/FBA_functions.R")
load(paste0(datadir, "TenSpecVary_relAbund_chemostat_fixRevRxns_mimosa_pieces.rda"))
media_file = paste0(homedir, "DATA/GutMedia/FaithMedia_AGORA_F_final.csv")
dictionary_file = paste0(homedir, "DATA/Dictionary_AGORA_complete.csv")
met_key_file = paste0(homedir, "/hmdb_brite_metabolite_classes.txt")
spec_file = paste0(datadir2, "allSpeciesFinalMain.txt")
met_file = paste0(datadir2, "allMetabolitesFinalMain.txt")
all_fluxes_file = paste0(datadir2, "allMetFluxesFinalMain.txt")

media = fread(media_file)
media[,V1:=gsub("[e]", "[env]", V1, fixed = T)]
setnames(media, "V1", "medium")

spec_codes = make_spec_codes()

dictionary = fread(dictionary_file)
dictionary[,medium:=gsub("[e]", "[env]", medium, fixed = T)]
dictionary[Primary=='L-lysinium(1+)', Primary:="L-lysine"]
dictionary[,Primary:=gsub(" \\(.*\\)$","", Primary)]
dictionary[,Primary:=gsub("\\(.*\\)$","", Primary)]
dictionary[Primary=='L-argininium', Primary:="L-arginine"]
dictionary[Primary=="proton", Primary:="H+"]
dictionary[Primary=="hydrogenphosphate", Primary:="Orthophosphate"]
dictionary[grepl("Amylopectin", Primary), Primary:="Amylopectin"]
setnames(dictionary, "Primary", "Metabolite")
kegg_translate = dictionary[,list(Metabolite, medium)]
dictionary[,niceLab:=tolower(Metabolite)]

path_key = fread(met_key_file)
path_key = merge(path_key, dictionary[,list(niceLab, medium)], by = "niceLab", all = T)
path_key[,Metabolite:=NULL]

##Color scale for plotting species
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
col_spec = c(getPalette(10), "grey50", "grey20", "grey80")
names(col_spec) = c(spec_codes[1:10][order(SpeciesName), Code], "in", "out", "Resid")
col_spec2 = col_spec
names(col_spec2) = c(spec_codes[1:10][order(SpeciesName), SpeciesName], "Inflow", "Outflow")


timePointFinal = 577

fake_spec = fread(spec_file)
fake_spec[,SampleID:=SimRun]
fake_spec = merge(fake_spec, spec_codes, by = c("Species", "Code"))
fake_spec[,SpeciesName:=factor(SpeciesName, levels = rev(spec_codes[,SpeciesName]))] #Order by
fake_spec[,Sample:=paste0("run", SimRun, "_", paste0(spec_codes[1:10,Code], collapse=""), "_TP", TimePoint)]

fake_mets_melt = fread(met_file) #melt(fake_mets, id.vars=c("Metabolite","medium"), variable.name = "Sample")
fake_mets_melt = merge(fake_mets_melt, path_key, by.x = c("compound", "niceLab"), by.y = c("medium", "niceLab"), all.x=T)
fake_mets_melt[,Sample:=paste0("run", SimRun, "_", paste0(spec_codes[1:10,Code], collapse=""), "_TP", TimePoint)]

met_order = fake_mets_melt[,var(value), by=niceLab][order(V1, decreasing=T), niceLab]
met_summary = fake_mets_melt[,list(mean(value,na.rm=T), var(value,na.rm=T), sd(value, na.rm=T)/mean(value, na.rm=T)),by=list(niceLab,compound)]
setnames(met_summary, c("V1", "V2", "V3"), c("Mean", "Variance", "CoefVar"))

all_met_fluxes_final1 = fread(all_fluxes_file)
all_met_fluxes_final1[,niceLab:=factor(niceLab, levels = met_order)]

fake_spec_wide = dcast(fake_spec, Species~Sample, value.var = "value")
contrib_table_genes = FBA_picrust_contribs(fake_spec_wide, mimosa_pieces_revFix[[4]])
contrib_table_genes[,Sample:=gsub("revFix", "run", Sample)]
setnames(contrib_table_genes, c("Species", "RxnID"), c("OTU","Gene"))
spec_codes = spec_codes[1:10]
all_koAbunds_byOTU_genes = contribs_by_species_list(contrib_table_genes, "singleMusicc", "FBAout", write_out = F)[match(spec_codes[,Code], sort(spec_codes[,Code]))]
emms_revFix = mimosa_pieces_revFix[[1]]
cmps_alone = singleSpecCMP_FBA(spec_codes[,Code], all_koAbunds_byOTU_genes, emms_revFix)


shap_contribs1 = getContributions(all_met_fluxes_final1, spec_codes[Code != "out"], kegg_translate, path_key)
shap_contribs1[,niceLab:=factor(tolower(niceLab), levels=met_order)]

var_cutoff = met_summary[Variance !=0, quantile(Variance, 0.25)*1.000001] #handle numerical issues

shap_contribs1[TrueVar > var_cutoff, unique(niceLab)]
shap_contribs = shap_contribs1[TrueVar > var_cutoff]
orig_mets = shap_contribs[,unique(compound)]
contrib_threshold_pos = 0.1
contrib_threshold_mag = 0.2
fake_mets_melt_good = fake_mets_melt[niceLab %in% shap_contribs[,niceLab]]
spec_met_corrs = get_correlation_contrib_comparison(fake_spec, fake_mets_melt_good, shap_contribs)

#New stuff
species_cmps = rbindlist(lapply(1:length(cmps_alone), function(x){
  foo = melt(cmps_alone[[x]], id.var = "compound", variable.name = "Sample")
  foo[,Species:=spec_codes[x,Species]]
  return(foo)
}))
species_cmps[,compound:=gsub("[e]", "[env]", compound, fixed = T)]
species_cmps = species_cmps[compound %in% orig_mets]

cmp_mods = fit_cmp_mods(species_cmps, fake_mets_melt)
species_cmps = add_residuals(species_cmps, cmp_mods[[1]], cmp_mods[[2]])
var_shares = calculate_var_shares(species_cmps)

var_shares_compare = merge(var_shares, spec_met_corrs, by = c("compound", "Species"), all.x = T)

var_shares_compare[,cor(VarShare.x, VarShare.y, use = "complete.obs")]
var_shares_compare[,summary(VarShare.x[Species=="Residual"])]
var_shares_compare[VarShare.x==1 & Species=="Residual"]
var_shares_compare = var_shares_compare[Var > var_cutoff]
bad_mets = var_shares_compare[VarShare.x > 0.5 & Species=="Residual",compound]
bad_mets2 = var_shares_compare[VarShare.x > 0.3 & Species=="Residual",compound]
var_shares_compare[compound %in% bad_mets]
#what to do with these...

cbind(var_shares_compare[Species=="Residual", list(compound,VarShare.x, 1-VarShare.x)], cmp_mods[[1]]) #Confirm that residual var share=1-rsq from model

var_shares_compare[,PredContrib:=ifelse(VarShare.x > contrib_threshold, 1, 0)]
var_shares_compare[,PosVarShare.x:=ifelse(V1.x > 0, V1.x/sum(V1.x[V1.x > 0]),0), by=compound]
var_shares_compare[,VarShareMagnitude.x:=abs(V1.x)/sum(abs(V1.x)), by=compound]
var_shares_compare[,PosPredContrib:=ifelse(PosVarShare.x > contrib_threshold_pos, 1, 0)]
resid_vals = var_shares_compare[Species=="Residual", list(V1.x, VarShare.x, compound)]

var_shares_compare[,prop.table(table(PosPredContrib, BinaryContribPos,compound %in% bad_mets))] #Ok this kind of works
var_shares_compare[,table(PosPredContrib, p.value < 0.01)] #Not the same

var_shares_compare
roc(predictor = var_shares_compare[,PosVarShare.x], response =  var_shares_compare[,BinaryContribPos], ci = T, ci.method = "bootstrap", boot.n = 500)
roc(predictor = var_shares_compare[!compound %in% bad_mets & Species != "Residual",PosVarShare.x], response =  var_shares_compare[!compound %in% bad_mets & Species != "Residual",BinaryContribPos], ci = T, ci.method = "bootstrap", boot.n = 500)
#Is this better than MIMOSA??

roc(predictor = var_shares_compare[!compound %in% bad_mets2 & Species != "Residual",PosVarShare.x], response =  var_shares_compare[!compound %in% bad_mets & Species != "Residual",BinaryContribPos], ci = T, ci.method = "bootstrap", boot.n = 500)

#Can we fit a model where every species has its own scaling coefficient but they are constrained somehow? prior question again. Or maybe every reaction?

var_share_accuracy = var_shares_compare[Species != "Residual",list(sum(PosPredContrib==BinaryContribPos)/length(PosPredContrib), sum(PosPredContrib[BinaryContrib==1])/sum(BinaryContrib==1), sum((PosPredContrib==0)[BinaryContrib==0])/sum(BinaryContrib==0), sum(PosPredContrib[BinaryContrib==1])/sum(PosPredContrib==1), sum(BinaryContribPos), sum(PosPredContrib)), by=compound]
setnames(var_share_accuracy, c("compound", "Accuracy", "Sensitivity", "Specificity", "PPV", "NumTruePos", "NumFalsePos"))
var_share_accuracy = merge(var_share_accuracy, cmp_mods[[1]], by="compound")
ggplot(var_share_accuracy, aes(x=Rsq, y = Sensitivity)) + geom_jitter()
ggplot(var_share_accuracy, aes(x=Rsq, y = PPV)) + geom_jitter()
ggplot(var_share_accuracy, aes(x=Sensitivity==1, y = Rsq)) + geom_violin() + geom_jitter()
ggplot(var_share_accuracy, aes(x=Accuracy==1, y = Rsq)) + geom_violin() + geom_jitter()
ggplot(var_share_accuracy, aes(x=PPV==1, y = Rsq)) + geom_violin()+geom_jitter()

var_share_accuracy[Rsq < 0.4]
ggplot(species_cmps[compound %in% bad_mets], aes(x=Sample, y= value, color=Species)) + geom_point() +
  facet_wrap(~compound, scales = "free") + theme(axis.text.x = element_blank())
species_cmps = merge(species_cmps, fake_mets_melt, by = c("compound", "Sample"), all.x =T)
species_cmps = merge(species_cmps, spec_codes, by = "Species", all.x = T)
species_cmps[Species=="Residual", Code:="Resid"]

ggplot(species_cmps[compound %in% bad_mets], aes(x=value.x, y= value.y, color=Code)) + geom_point() +
  facet_wrap(~compound, scales = "free") + scale_color_manual(values= col_spec)

tot_cmps = species_cmps[Species != "Residual",sum(value.x), by=list(Sample, compound, value.y)]


tot_cmp_conc_result = ggplot(tot_cmps[compound %in% bad_mets], aes(x=V1, y = value.y)) + geom_point() + facet_wrap(~compound, scales = "free")


#Some sort of iterative fitting where we flip reaction directions to deal with this?
var_shares_compare[,cor(V1.x, V1.y, use = "complete.obs", method="spearman"), by=compound][,table(abs(V1) > 0.9)]


met_order = var_share_accuracy[order(Rsq, decreasing=T), compound]
tot_cmps = merge(tot_cmps, var_share_accuracy, by = "compound")
tot_cmps[,compound:=factor(compound, levels = met_order)]
cmp_scaling_plot = ggplot(tot_cmps, aes(x=V1, y = value.y, col = Slope < 0)) + geom_point() + facet_wrap(~compound, scales = "free") +
  xlab("CMP Score") + ylab("Metabolite conc") + theme(axis.text = element_text(size=7), strip.text = element_text(size=7))
save_plot(cmp_scaling_plot, file = "CMP_scaling_simulation.png", base_width = 12, base_height = 9)
#Steps
#Set up nice pipeline to run these algorithms on simulation data
#Then a real data
# then can try some iterative model improvements
# # Plotting function for var shares


#Take a look at # of reactions for each compound and species
orig_mets_v2 = gsub("[env]", "[e]", orig_mets, fixed = T)
num_rxns = data.table(rbindlist(lapply(1:length(emms_revFix), function(x){
  foo = rbindlist(lapply(orig_mets_v2, function(y){
    if(y %in% emms_revFix[[x]][,Compound]){
      dat = emms_revFix[[x]][,lapply(.SD, function(z){ if(z[Compound==y] != 0) return(z[Compound==y])})]
      dat = melt(dat, id.var = "Compound", variable.name = "Rxn")
      return(dat)
    }
  }), fill = T)
  foo[,Species:=spec_codes[x,Species]]
  return(foo)
}), fill = T))
num_rxns[,length(Rxn), by=list(Compound, Species, value)]
rxn_count = num_rxns[,list(length(Rxn[value < 0]), length(Rxn[value > 0])), by=list(Compound, Species)]
setnames(rxn_count, c("Compound", "V1", "V2"), c("compound", "DegRxns", "SynthRxns"))

rxn_count[,compound:=gsub("[e]", "[env]", compound, fixed = T)]
rxn_count[compound %in% bad_mets]

rxn_count[,list(mean(DegRxns), mean(SynthRxns)), by=compound %in% bad_mets]
rxn_count_tot = rxn_count[,list(sum(DegRxns), sum(SynthRxns)), by=compound]

rxn_count_tot[,fisher.test(table(V2 > 0 & V1>0, compound %in% bad_mets))] #meh


