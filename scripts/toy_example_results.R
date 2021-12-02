## Toy example for paper/docs
library(mimosa)
library(data.table)
library(cowplot)


## From evaluation project
toy_colors = c("Taxon 1" = "#05D673", "Taxon 2" = "#ED7D31" , "Taxon 3" = "#FFC000", "Taxon 4" = "#5B9BD5")
toy_colors = c("Taxon 1" = "#FFC000", "Taxon 2" = "#FF0000", "Taxon 3" = "#00B0F0" )

toy_cmps = data.table(Taxon = c(rep("Taxon 1", 6), rep("Taxon 2", 6), rep("Taxon 3", 6)), Sample = rep(c("Sample C", "Sample F", "Sample D", "Sample A", "Sample E", "Sample B"), 3), CMP = c(7.3,6, 10, 5,7,8, 3,4,6,4.3,5,6,  -5, -3, -8, -7, -1, -8.5)) 

toy_cmps = data.table(Taxon = c(rep("Taxon 1", 6), rep("Taxon 2", 6), rep("Taxon 3", 6)), 
                      Sample = rep(c("Sample C", "Sample F", "Sample D", "Sample A", "Sample E", "Sample B"), 3), 
                      CMP = c(7.3,6, 10, 5,7,8, 3,4,4,4.3,5,6,  -5, -2, -8, -7, -1, -8.5)) 

#, rep("Residual Error", 6)

toy_conc = toy_cmps[,sum(CMP), by=Sample]
toy_conc[,Conc:=c(7, 13, 9, 4, 10, 5.5)]
toy_conc[,Conc:=c(6, 14, 12, 4, 10, 5.5)]
toy_model = toy_conc[,lm(Conc~V1)]
toy_model_rank = toy_conc[,Rfit::rfit(Conc~V1)]

1-toy_model_rank$D1/toy_model_rank$D0
summary(toy_model)

toy_conc[,Taxon:="Total"]
setnames(toy_conc, "V1", "CMP")
toy_conc[,Sample2:=gsub("Sample ", "", Sample)]
plot_cmp_met = ggplot(toy_conc, aes(x=CMP, y = Conc, label = Sample2)) + geom_abline(slope = coef(toy_model)[2], intercept = coef(toy_model)[1], linetype = 2, color = "springgreen4", size = 1.15) +
  geom_abline(slope = coef(toy_model_rank)[2], intercept = coef(toy_model_rank)[1], linetype = 2, color = "slateblue3", size = 1.15) + 
  annotate(geom = "text", x=9.3, y = 6, label = paste0("OLS ", "~R^2", ": 0.43"), color = "springgreen4", parse = T, size = 4)+ geom_point(size = 2) + geom_text(nudge_y = -0.4, size = 3.5) + annotate(geom = "text", x=9.3, y = 5.2, label = paste0("Rank  ", "~R^2", ": 0.31"), color = "slateblue3", parse = T, size = 4) +
  xlab("Total CMP") + ylab("Concentration [M]")

toy_conc2 = toy_conc[,list(Sample, Conc)]
toy_conc2[,Taxon:="Concentration [M]"]
setnames(toy_conc2, "Conc", "CMP")
toy_cmps_conc = rbind(toy_cmps, toy_conc, toy_conc2, fill = T)

# ggplot(toy_fluxes_conc[Sample=="A"], aes(x=Sample, y = Flux, col=factor(Taxon), shape=factor(Taxon))) + geom_point(size=4) + geom_hline(yintercept=0, linetype=2) + scale_shape_manual(values = c(8,15:18)) + scale_color_manual(values = c(toy_colors, "Concentration" = "black")) + guides(color = guide_legend(title = ""), shape = guide_legend(title = ""))+ ylim(-5.5, 16.5) + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle=90)) + ylab("[M] Flux")
# 
# ggplot(toy_fluxes, aes(x=Sample, y = Flux, col=factor(Taxon), shape=factor(Taxon))) + geom_point(size=3) + geom_hline(yintercept=0, linetype=2) + scale_shape_manual(values = c(15:18)) + scale_color_manual(values = toy_colors) + guides(color = guide_legend(title = ""), shape = guide_legend(title = ""))+ ylim(-9.5, 10.5) + theme(axis.text.x = element_blank())

toy_cmps_conc[,Taxon2:=ifelse(grepl("Concentration", Taxon), Taxon, paste0(Taxon, " CMP"))]
toy_color2 = toy_colors
names(toy_color2) = paste0(names(toy_colors), " CMP")
toy_color2 = c(toy_color2,  "Total CMP" = "gray40", "Concentration [M]" = "black")
toy_cmps_conc[,Sample2:=gsub("Sample ", "", Sample)]
# toy_alpha = c(rep(0.7, 3), 1, 1)
# names(toy_alpha) = names(toy_color2)
# toy_size = c(rep(3, 4), 9)
# names(toy_size) = names(toy_color2) #Shapes c(42, 15:18)+ scale_alpha_manual(values = toy_alpha) + scale_size_manual(values = toy_size) , alpha = factor(Taxon2), size = factor(Taxon2)
plot1 = ggplot(toy_cmps_conc, aes(x=Sample2, y = CMP, col=factor(Taxon2), shape=factor(Taxon2))) + geom_point(size = 3, stroke = 1.15) + geom_hline(yintercept=0, linetype=2) + scale_shape_manual(values = c(19, 0,1,2, 18)) +  
  scale_color_manual(values = toy_color2) +
  guides(color = guide_legend(title = ""), shape = guide_legend(title = ""), alpha = guide_legend(title = ""), size = guide_legend(title = ""))+ theme(axis.text = element_text(size = 9)) + 
  ylab("CMP/Concentration") + xlab("Sample")#, axis.text.x = element_text(angle=90, vjust=0.5))#+ ylim(-5.5, 16.5)


#Get contributions
species_contrib_table = toy_cmps
species_contrib_table[,compound:="M"]
species_contrib_table[,Species:=Taxon]
met_table = toy_conc[,list(Sample, Conc)]
met_table[,compound:="M"]
setnames(met_table, "Conc", "value")
mods1 = fit_cmp_mods(species_contrib_table, met_table, rank_based = F)
mods2 = fit_cmp_mods(species_contrib_table, met_table, rank_based = T)
var_shares1 = calculate_var_shares(species_contrib_table, met_table, mods1, config_table = data.table(V1 = ""))
var_shares2 = calculate_var_shares(species_contrib_table, met_table, mods2, config_table = data.table(V1 = "rankBased"))
var_shares1[,Method:="OLS"]
var_shares2[,Method:="Rank-based"]
var_shares_all = rbind(var_shares1, var_shares2, fill = T)
contrib_plot = plot_contributions(var_shares_all, metabolite = "M", metIDcol = "compound", color_palette = toy_colors, order_spec = F) + 
  facet_wrap(Method~., nrow = 2) + theme(strip.text = element_text(), plot.title = element_blank())
contrib_plot_rank = plot_contributions(var_shares2, metabolite = "M", metIDcol = "compound", color_palette = toy_colors, order_spec = F) + theme(plot.title = element_blank())

toy_example = plot_grid(plot1, plot_grid(plot_cmp_met, contrib_plot + theme(legend.position = "right"), nrow = 1, labels = c("B", "C")), nrow = 2, labels = c("A", ""))


library(GGally)
library(network)
library(sna)
library(ggnetwork)
net = matrix(c(0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0), nrow = 4, ncol = 4)
net = network(net, directed = T)
net %v% "label" = c("M", "", "", "")
net %v% "TaxonN" = c("", "Taxon 1", "Taxon 2", "Taxon 3")
net %e% "Taxon" = c("Taxon 1", "Taxon 2", "Taxon 3")
ggnet(net, arrow.size = 12)
net2 = ggnetwork(net)
net2[net2$TaxonN == "Taxon 3", "x"] = 1
net2[net2$TaxonN == "Taxon 3", "y"] = 0.5
net2[net2$label == "M", "x"] = 0.52
net2[net2$label == "M", "y"] = 0.50
net2[net2$TaxonN == "Taxon 3", "yend"] = 0.5
net2[6,"yend"] = 0.5
net2[net2$TaxonN %in% c("Taxon 1", "Taxon2"), "yend"] = 0.5
net2[net2$TaxonN %in% c("Taxon 1", "Taxon2"), "xend"] = 0.5
net2[1, "xend"] = 0.5
net2[net2$TaxonN %in% c("Taxon 1", "Taxon 2"), "x"] = 0.1
net2[net2$TaxonN == "Taxon 1", "y"] = 0.8
net2[net2$TaxonN == "Taxon 2", "y"] = 0.2
net2[net2$TaxonN == "Taxon 2", "xend"] = 0.5
toy_diagram = ggplot(net2,  aes(x = x, y = y, xend = xend, yend = yend)) + geom_edges(arrow = arrow(length = unit(14, "pt"), type = "closed"), aes(color = Taxon), size = 1.6)  +  #+ geom_nodes(aes(color = TaxonN))
  theme_blank() + scale_color_manual(values = toy_colors) + guides(color = F) + annotate(geom = "label", x = net2$x[1], y = net2$y[1], label = "M", size = 6) + xlim(0, 1) + ylim(0, 1)


toy_example = plot_grid(plot_grid(toy_diagram + theme(plot.margin = margin(0.45,0.3, 0.45, 0.3, "inches")), plot1, nrow = 1, rel_widths = c(1, 2.3), labels = c("A", "B")), 
                        plot_grid(plot_cmp_met + theme(plot.margin = margin(0.2, 0.35, 0.2, 0.2, "inches")), contrib_plot + theme(legend.position = "right"), nrow = 1, labels = c("C", "D"), rel_widths = c(1.2, 1)), nrow = 2, rel_heights = c(1, 1.17))
save_plot(toy_example, file = "paper/toyExample.png", base_width = 8.7, base_height = 6.5)

toy_example2 = plot_grid(plot_grid(toy_diagram + theme(plot.margin = margin(0.45,0.3, 0.45, 0.3, "inches")), plot1 + theme(plot.margin = margin(0.25, 0.1, 0.1, 0.1, "inches")), nrow = 1, rel_widths = c(1, 2.3), labels = c("A) Reaction network", "B) Taxon-level CMP scores"), hjust=-0.1), 
                        plot_grid(plot_cmp_met + theme(plot.margin = margin(0.25, 0.35, 0.2, 0.2, "inches")), contrib_plot + theme(legend.position = "right", plot.margin = margin(0.25, 0.1, 0.1, 0.1, "inches")), nrow = 1, labels = c("C) CMP-Metabolite comparison", "D) Taxon contributions"), hjust = -0.1, rel_widths = c(1.2, 1)), nrow = 2, rel_heights = c(1, 1.17))
save_plot(toy_example2, file = "paper/toyExampleTitles.png", base_width = 8.7, base_height = 6.5)
toy_example3 = plot_grid(plot_grid(toy_diagram + theme(plot.margin = margin(0.45,0.3, 0.45, 0.3, "inches")), plot1 + theme(plot.margin = margin(0.25, 0.1, 0.1, 0.1, "inches")), nrow = 1, rel_widths = c(1, 2.3), labels = c("A) Reaction network", "B) Taxon-level CMP scores"), hjust=-0.1), 
                         plot_grid(plot_cmp_met + theme(plot.margin = margin(0.25, 0.35, 0.2, 0.2, "inches")), contrib_plot_rank + theme(legend.position = "right", plot.margin = margin(0.5, 0.1, 0.6, 0.1, "inches")), nrow = 1, labels = c("C) CMP-Metabolite comparison", "D) Taxon contributions"), hjust = -0.1, rel_widths = c(1.2, 1)), nrow = 2, rel_heights = c(1, 1.17))
save_plot(toy_example3, file = "paper/toyExampleTitlesRank.png", base_width = 8.7, base_height = 6.5)
save_plot(toy_example3, file = "paper/toyExampleTitlesRank.pdf", base_width = 8.7, base_height = 6.5)

toy_var_shares = rbindlist(lapply(c("Taxon 1", "Taxon 2","Taxon 3"), function(y){
  all1 = rbindlist(lapply(c("Taxon 1", "Taxon 2","Taxon 3"), function(x){
    foo = data.table(toy_fluxes_fill[,cov(get(as.character(x)), get(as.character(y)), use="complete.obs")])
    foo[,Species:=x]
    return(foo)
  }))
  all1[,Species2:=y]
}))
toy_var_shares = toy_var_shares[,sum(V1),by=Species]
toy_var_shares[,TrueVar:=toy_conc[,var(Flux)]]
toy_var_shares[,VarShare:=V1/TrueVar]

ggplot(toy_var_shares, aes(y=VarShare, x = factor(Species, levels = rev(sort(Species))), fill = Species)) + geom_bar(stat = "identity") + scale_fill_manual(values = toy_colors) + geom_abline(intercept = 0, slope = 0, linetype = 2) + theme(strip.background = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(size=10),axis.title = element_text(size=16), axis.text.y = element_blank(), legend.title = element_blank(), legend.text = element_text(size=15), strip.text = element_text(size=17)) + ylab("Contribution to variance") + xlab("") + coord_flip()# 


toy_data = data.table()