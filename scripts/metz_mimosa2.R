### Re-analysis of Snijders data for MIMOSA2 paper
library(data.table)
library(seqsetvis)
library(mimosa)
library(cowplot)
library(pROC)

load("~/Documents/MIMOSA2shiny/results/fixed/metz_unc_otus_good_all.rda")

diet_content = fread("~/Google Drive File Stream/My Drive/Genome_Sciences/metzMiceMetabolites/diet_mets_processed_autocl.txt")
diet_mets = diet_content[meanUNC_A != 0, KEGG]
diet_mets = diet_mets[diet_mets %in% mets[,compound]]
non_diet_mets = mets[!compound %in% diet_mets, compound]

consist_mets_compare = list(mimosa1 = consist_mets, mimosa1contrast = contrast_mets, 
                            m2agora_noRev = m2agora_noRev_mets, m2kegg = m2kegg_mets, 
                           m2edit_kegg = m2edit_kegg_mets, # m2edit_agora = m2edit_agora_mets, 
                            correlation = corr_mets, correlation_gene = corr_mets_gene, in_diet = diet_mets, not_in_diet = non_diet_mets)

consist_mets_compare2 = list( # m2edit_agora = m2edit_agora_mets, 
  in_diet = diet_mets, not_in_diet = non_diet_mets,
                            correlation = corr_mets, correlation_gene = corr_mets_gene, 
                            mimosa1 = consist_mets, 
                            mimosa2 = m2edit_kegg_mets)

plot1 = ssvFeatureVenn(consist_mets_compare[c("mimosa1", "m2agora", "m2kegg")])
plot2 = ssvFeatureVenn(consist_mets_compare[c("mimosa1", "m2edit_kegg", "correlation")])
plot3 = ssvFeatureEuler(consist_mets_compare)
compare_met_dat = get_overlaps(consist_mets_compare)
plot4 = ssvFeatureEuler(consist_mets_compare[c("in_diet", "not_in_diet", "m2edit_kegg", "correlation", "mimosa1")])
method_colors = brewer.pal(4, "Set1")
method_colors = c("gray50", "gray80", method_colors)
#names(method_colors) = c("correlation", "correlation_gene", "mimosa1", "mimosa2", "in_diet", "not_in_diet")

names(consist_mets_compare2) = c("In diet", "Not in diet", "Correlation", "Correlation w/ reaction presence", 
                                 "MIMOSA1", "MIMOSA2")
names(method_colors) = c("In diet", "Not in diet", "Correlation", "Correlation w/ reaction presence", 
                         "MIMOSA1", "MIMOSA2")
plot4dat = ssvFeatureEuler(consist_mets_compare2, return_data = T)
plot4dat$group_names = factor(plot4dat$group_names, levels = c("In diet", "Not in diet", "Correlation", "Correlation w/ reaction presence", 
                                                               "MIMOSA1", "MIMOSA2"))
cn = colnames(object)
eu = eulerr::euler(object, shape = shape)
dd = eu$ellipses
h <- dd$h
k <- dd$k
a <- dd$a
b <- dd$b
phi <- dd$phi
p = ggellipse(xcentres = plot4dat$xcentres, ycentres = plot4dat$ycentres, r = plot4dat$r, r2 = plot4dat$r2, 
              phi = plot4dat$phi, circle_colors = method_colors, group_names = plot4dat$group_names, 
              line_alpha = 1, fill_alpha = 0.3, line_width = 2, 
              n_points = 200)
p

method_colors2 = alpha(method_colors, alpha = 0.3)
names(method_colors2) = c("In diet", "Not in diet", "Correlation", "Correlation w/ reaction presence", 
                          "MIMOSA1", "MIMOSA2")
plot4 = ssvFeatureEuler(consist_mets_compare2, circle_colors = method_colors) + theme(legend.position = "bottom") +
  scale_fill_manual(values = method_colors2, breaks = c("In diet", "Not in diet", "Correlation", "Correlation w/ reaction presence","MIMOSA1", "MIMOSA2")) +
  scale_color_manual(values = method_colors, breaks = c("In diet", "Not in diet", "Correlation", "Correlation w/ reaction presence", 
                                                                                                                  "MIMOSA1", "MIMOSA2")) 

#+ scale_alpha_manual(values = c(0.3, 1))
#+ 
  #guides(color = guide_legend(order = c("correlation", "correlation_gene", "mimosa1", "mimosa2", "in_diet", "not_in_diet")))

compare_met_grid = get_overlaps_grid(consist_mets_compare, melt = F)
compare_met_grid[,table(m2edit_kegg, mimosa1)]
compare_met_grid[,fisher.test(table(in_diet, mimosa1))]
compare_met_grid[,fisher.test(table(in_diet, m2edit_kegg))]
compare_met_grid[,table(not_in_diet, mimosa1)]

varShare_threshold_met2 = 0.05
mimosa2kegg = results_kegg$varShares[VarShare > varShare_threshold_met2 & Species != "Residual", paste0(compound, "_", Species)]
mimosa2agora = results_agora$varShares[VarShare > varShare_threshold_met2 & Species != "Residual", paste0(compound, "_", Species)]
mimosa2kegg_edit = mimosa2_results_edit_kegg$varShares[VarShare > varShare_threshold_met2 & Species != "Residual", paste0(compound, "_", Species)]
mimosa2agora_edit = mimosa2_results_edit_agora$varShares[VarShare > varShare_threshold_met2 & Species != "Residual", paste0(compound, "_", Species)]
mimosa2agora_noRev = results_agora_noRev$varShares[VarShare > varShare_threshold_met2 & Species != "Residual", paste0(compound, "_", Species)]

corr_list = spec_met_corrs[Qval < corr_threshold, paste0(compound, "_", Species)]
corr_list_gene = spec_met_corrs[Qval < corr_threshold & hasGene==1, paste0(compound, "_", Species)]


contribs_compare = list(mimosa1 = contrib_list_cmps, mimosa1_mets = contrib_list_mets,
                        #mimosa1_varshares = contrib_list_varshares, mimosa1_covshares = contrib_list_covshares, 
                        mimosa2_agora_norev = mimosa2agora_noRev,
                        mimosa2kegg = mimosa2kegg, mimosa2agora = mimosa2agora, m2kegg_edit = mimosa2kegg_edit, m2agora_edit = mimosa2agora_edit,
                        correlation = corr_list, correlation_gene = corr_list_gene)
plot5 = ssvFeatureEuler(contribs_compare)
plot5a = ssvFeatureEuler(contribs_compare[c("mimosa1_mets", "m2kegg_edit", "correlation")])
compare_dat = get_overlaps(contribs_compare)
signif_taxa = mimosa2_results_edit_kegg$varShares[VarShare > 0.02, unique(Species)]

taxonomy = fread("~/Google Drive File Stream/My Drive/Genome_Sciences/KEGGfiles/97_otu_taxonomy.txt", header = F)
taxonomy[,OTU:=gsub("\t.*", "", V1)]
taxonomy[,Kingdom:=gsub(".*\t", "", V1)]
setnames(taxonomy, c("V2", "V3", "V4", "V5", "V6", "V7"), c("Phylum", "Class", "Order", "Family", "Genus", "Species"))
var_dat = mimosa2_results_edit_kegg$varShares
var_dat = merge(var_dat, taxonomy, by.x = "Species", by.y = "OTU", all.x = T)
var_dat[,TaxonName:=ifelse(Species.y=="s__", Species, paste0(Species.y, "_", Species))]
var_dat[,TaxonName:=ifelse(Genus=="g__", TaxonName, paste0(Genus, "_", TaxonName))]
var_dat[,TaxonName:=ifelse(Family=="f__", TaxonName, paste0(Family, "_", TaxonName))]
var_dat[,TaxonName:=ifelse(Order=="o__", TaxonName, paste0(Order, "_", TaxonName))]

plot_summary_contributions(mimosa2_results_edit_kegg$varShares[compound %in% m2edit_kegg_mets & Species %in% signif_taxa], include_zeros = T)
var_dat[,metID:=met_names(as.character(compound))]
met_order = var_dat[Species=="Residual"][order(VarShare, decreasing = F), metID]
var_dat[,metID:=factor(metID, levels = met_order)]
resid_dat = var_dat[Species == "Residual"]
varShares = var_dat[Species != "Residual"]
spec_order = varShares[,length(VarShare[abs(VarShare) > 0.05]), by=TaxonName][order(V1), TaxonName]
varShares[,Species2:=factor(as.character(TaxonName), levels = spec_order)]
plot_var = "VarShare"
color_lab = "Contribution to variance"
resid_plot = ggplot(resid_dat[compound %in% m2edit_kegg_mets], aes(x=metID, y = 1-VarShare)) + geom_bar(stat = "identity") + scale_y_continuous(expand = c(0,0))+ theme_minimal() +
  theme(axis.text.x = element_text(angle=90, hjust=0, vjust =0.5), axis.line = element_blank(), axis.title.x = element_blank()) + ylab("Model R-squared") + theme(axis.title.y = element_blank())
plot1 = ggplot(varShares[compound %in% m2edit_kegg_mets & Species %in% signif_taxa], aes(x=metID, y = Species2)) + geom_tile(aes_string(fill = plot_var)) + theme_minimal() +
  theme(axis.text.x = element_blank(), axis.line = element_blank(), legend.position = "bottom") +
  scale_fill_gradient2(low = brewer.pal(9,"Reds")[9], mid = "white",  high = brewer.pal(9, "Blues")[9], midpoint = 0, name = color_lab) +
  ylab("Taxon")+xlab("Metabolite")  + theme(axis.title.y = element_blank())


plot_all = plot_grid(resid_plot, plot1, nrow = 2, align = "v", axis = "r", rel_heights = c(1, 2.5))

paper_fig= plot_grid(plot4,plot1 + theme(axis.text.x = element_text(angle=90, hjust=0, vjust =0.5), axis.title.x = element_blank()) + scale_x_discrete(position = "top"), nrow = 2, rel_heights = c(1,2.5), labels = c("A", "B"))
save_plot(paper_fig, file = "results/metz_unc_results_figure.png", base_width = 9, base_height = 8.5)
