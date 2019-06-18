#### MIMOSA2 testing functions
process_abundances = function(species_file, met_file, kos_file = NULL, fluxes_file = NULL, simulated = F, config_table_list = NULL, return_fluxes = F){
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
    return(datafiles)
  } else {
    source("data/testData/sim_data/FBA_functions.R")
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
    if(return_fluxes == F){
      return(list(species, mets, contribs))
    } else {
      return(list(species, mets, contribs, met_fluxes))
    }
  }
}

compare_mimosa12 = function(m1results, m2results, contribs = NULL){
  compare_dat = merge(m1results, m2results, by = "compound", all = T)
  if("PValPos" %in% names(compare_dat)){
    compare_dat[,PValBoth:=ifelse(Correlation > 0, PValPos, PValNeg)]
    compare_dat[,PosNegCompare:=paste(Correlation > 0, Slope > 0, sep = "_")]
  } else {
    compare_dat[,PosNegCompare:=paste(Slope.x > 0, Slope.y > 0, sep = "_")]
  }
  if(!is.null(contribs)){
    compare_dat = merge(compare_dat, contribs[Species=="Inflow"], by = "compound", all.x = T)
    setnames(compare_dat, "VarShare", "EnvVarShare")
    #Compare contrib also
  } 
  return(compare_dat)
}

compare_mimosa12_plot = function(compare_dat, simulated = F){
  slope_color_scale = RColorBrewer::brewer.pal(7, "RdBu")
  slope_color_scale = c(slope_color_scale, "black", "green") #differentiate completely diff cases
  names(slope_color_scale) = c("FALSE_FALSE", "FALSE_NA", "NA_FALSE","NA_NA", "NA_TRUE", "TRUE_NA", "TRUE_TRUE", "FALSE_TRUE", "TRUE_FALSE")
  if("PValBoth" %in% names(compare_dat)){
    plot1 = ggplot(compare_dat, aes(x=Rsq, y = -log(PValBoth), color = PosNegCompare, label = compound)) + geom_point() + geom_text(aes(y=-0.985*log(PValBoth)), size = 2) + scale_color_manual(values = slope_color_scale) + guides(color = F)
    plot2 = ggplot(compare_dat, aes(x=sqrt(Rsq), y = abs(Correlation), color = PosNegCompare, label = compound)) + geom_abline(intercept = 0, slope = 1, linetype = 2) + geom_point() + geom_text(aes(y=0.985*abs(Correlation)), size = 2) + scale_color_manual(values = slope_color_scale) + theme(legend.title = element_blank())
    if(simulated){
      plot3 = ggplot(compare_dat, aes(y=sqrt(Rsq), x = sqrt(1-EnvVarShare), label = compound, color = ifelse((Rsq > 0.2)==(1-EnvVarShare > 0.2), T, F))) + geom_abline(intercept = 0, slope = 1, linetype = 2) + geom_point() + geom_text(aes(y=0.985*sqrt(Rsq)), size = 2) + guides(color = F)
      plot4 = ggplot(compare_dat, aes(y=abs(Correlation), x = sqrt(1-EnvVarShare), label = compound, color = ifelse((Correlation^2 > 0.2)==(1-EnvVarShare > 0.2), T, F))) + geom_abline(intercept = 0, slope = 1, linetype = 2) + geom_point() + geom_text(aes(y=0.985*abs(Correlation)), size = 2) + theme(legend.title = element_blank())
      final_plot = plot_grid(plot_grid(plot1, plot2, rel_widths = c(1, 1.2)), plot_grid(plot3, plot4, rel_widths = c(1, 1.2)), nrow = 2, align = "v", axis = "l")
    } else {
      final_plot = plot_grid(plot1, plot2, rel_widths = c(1, 1.4))
    }
  } else {
    plot1 = ggplot(compare_dat, aes(x=-log10(PVal.x), y = -log10(PVal.y), color = PosNegCompare, label = compound)) +
      geom_hline(yintercept = 2, linetype = 2) + geom_vline(xintercept = 2, linetype = 2) + geom_point() + 
      geom_text(aes(y=-0.985*log10(PVal.y)), size = 2) + scale_color_manual(values = slope_color_scale) + guides(color = F)
    plot2 = ggplot(compare_dat, aes(x=Rsq.x, y = Rsq.y, color = PosNegCompare, label = compound)) + geom_abline(intercept = 0, slope = 1, linetype = 2) + geom_point() + geom_text(aes(y=0.985*Rsq.y), size = 2) + scale_color_manual(values = slope_color_scale) + theme(legend.title = element_blank())
    if(simulated){
      plot3 = ggplot(compare_dat, aes(y=Rsq.x, x = 1-EnvVarShare, label = compound, color = ifelse((Rsq.x > 0.2)==(1-EnvVarShare > 0.2), T, F))) + 
        geom_abline(intercept = 0, slope = 1, linetype = 2) + geom_point() + geom_text(aes(y=0.985*Rsq.x), size = 2) + guides(color = F) +
        xlim(-0.5, 1.5)
      plot4 = ggplot(compare_dat, aes(y=Rsq.y, x = 1-EnvVarShare, label = compound, color = ifelse((Rsq.y > 0.2)==(1-EnvVarShare > 0.2), T, F))) + 
        geom_abline(intercept = 0, slope = 1, linetype = 2) + geom_point() + geom_text(aes(y=0.985*Rsq.y), size = 2) + 
        theme(legend.title = element_blank()) + xlim(-0.5, 1.5)
      final_plot = plot_grid(plot_grid(plot1, plot2, rel_widths = c(1, 1.2)), plot_grid(plot3, plot4, rel_widths = c(1, 1.2)), nrow = 2, align = "v", axis = "l")
    } else {
      final_plot = plot_grid(plot1, plot2, rel_widths = c(1, 1.4))
    }
    
  }
  return(final_plot)
}


cmp_met_compare = function(m_datasets, config_table, simulated = F, compare_dat = NULL, m1cmp = NULL, rank_based = F, 
                           score_transform = "", met_transform = "", usePvals = T, use_smooth = T){
  if("KEGG" %in% names(m_datasets[[2]])){
    met_var = "KEGG"
  } else {
    met_var = "compound"
  }
  mets_melt = melt(m_datasets[[2]], id.var = met_var, variable.name = "Sample")
  if(met_transform != ""){
    if("KEGG" %in% names(m_datasets[[2]])) setnames(mets_melt, "KEGG", "compound")
    met_var = "compound"
    mets_melt = transform_mets(mets_melt, met_transform = met_transform)
  }
  if(!is.null(m1cmp)){
    #Use M1 CMP scores instead
    cmp_scores = melt(m1cmp, id.var = "compound", variable.name = "Sample")
  } else {
    if(!simulated) m_network = build_metabolic_model(m_datasets[[1]], config_table) else m_network = build_metabolic_model(m_datasets[[1]], config_table, manual_agora = T)
    spec_cmp_scores = get_species_cmp_scores(species_table = m_network[[2]], m_network[[1]], 
                                             normalize = T, manual_agora = T)
    if(score_transform != ""){
      spec_cmp_scores = transform_cmps(spec_cmp_scores, score_transform = score_transform)
    }
    cmp_scores = spec_cmp_scores[,sum(CMP), by=list(Sample, compound)]
  }
  
  m_compare_mets = merge(cmp_scores, mets_melt, by.x=c("Sample", "compound"), by.y = c("Sample", met_var))
  
  bad_mets = m_compare_mets[,var(V1), by=compound][V1==0, compound]
  #Separate out by synth only, deg only, both
  synth_deg = m_compare_mets[,list(sum(V1 > 0), sum(V1 < 0)), by=compound]
  setnames(synth_deg, c("V1", "V2"), c("SynthSamples", "DegSamples"))
  synth_deg[,table(SynthSamples > 0, DegSamples > 0)]
  m_compare_mets = merge(m_compare_mets, synth_deg, by="compound")
  m_compare_mets[,CMPType:=ifelse(SynthSamples > 0 & DegSamples==0, "Synth", "Both")]
  m_compare_mets[,CMPType:=ifelse(DegSamples > 0 & SynthSamples==0, "Deg", CMPType)]
  setnames(m_compare_mets, c("V1", "value"), c("CMP", "Met"))
  
  if(!is.null(compare_dat)){ #If including m1/m2 results
    m_compare_mets = merge(m_compare_mets, compare_dat, by = "compound", all = T)
    #m1_order = m_compare_mets[order(Correlation, decreasing = T), unique(compound)]
    #m_compare_mets[,compound:=factor(compound, levels = m1_order)]
    m_compare_mets[,m2_R:=sqrt(Rsq)]
    m_compare_mets[,m2_R:=ifelse(Slope < 0, -1*m2_R, m2_R)]
    
    if(usePvals){
      m_compare_mets[,annotResult:=paste0(compound, "\n", round(PVal, 3), " ", round(PValBoth, 3))] #Just always compare with m1 I guess
    } else {
      m_compare_mets[,annotResult:=paste0(compound, "\n", round(Correlation, 2), " ", round(m2_R, 2))]
    }
    #print(m_compare_mets)
    m1_order = m_compare_mets[order(Correlation, decreasing = T), unique(annotResult)]
    m2_order = m_compare_mets[order(m2_R, decreasing = T), unique(annotResult)]
    m_compare_mets[,annotResult:=factor(annotResult, levels = m2_order)]
    m_compare_mets[,M2Signif:=ifelse(PVal < 0.01, T, F)]
    setnames(m_compare_mets, "PosNegCompare", "M1M2Positive")
    if(rank_based){
      m_compare_mets[,MetValue:=Met]
      m_compare_mets[,CMPValue:=CMP]
      m_compare_mets[,Met:=NULL]
      m_compare_mets[,Met:=rank(MetValue), by=compound]
      m_compare_mets[,CMP:=NULL]
      m_compare_mets[,CMP:=rank(CMPValue), by=compound]
    }
    if(simulated){
      m_compare_mets[,EnvMet:=ifelse(EnvVarShare > 0.5, T, F)]
      met_plot = ggplot(m_compare_mets[!compound %in% bad_mets], aes(x=CMP, y = Met)) + geom_smooth(
        aes(linetype = M1Signif, color = M1M2Positive), span = 0.85) + 
        geom_point(alpha = 0.6, size = 0.8, aes(shape = EnvMet)) + facet_wrap(~annotResult, scales = "free") + 
        theme(strip.text = element_text(size = 9))
    } else {
      met_plot = ggplot(m_compare_mets[!compound %in% bad_mets], aes(x=CMP, y = Met))  + 
        geom_point(alpha = 0.6, size = 0.8) + facet_wrap(~annotResult, scales = "free") + 
        theme(strip.text = element_text(size = 9)) 
      if(use_smooth){
        met_plot = met_plot + geom_smooth(
          aes(linetype = M1Signif, color = M1M2Positive), span = 0.85)
      } else {
        met_plot = met_plot + geom_abline(aes(slope = Slope, intercept = Intercept, linetype = M2Signif, color = M1M2Positive))
      }
    }
  } else {
    #Plot by synth/deg
    comp_order = m_compare_mets[order(CMPType), unique(compound)]
    m_compare_mets[,compound:=factor(compound, levels = comp_order)]
    met_plot = ggplot(m_compare_mets[!compound %in% bad_mets], aes(x=CMP, y = Met)) + geom_smooth(aes(color = CMPType)) + 
      geom_point() + facet_wrap(~compound, scales = "free")
  }
  
  return(met_plot)
}

cmp_met_compare_setup = function(m_datasets, config_table, simulated = F, rank_based_plot = F, 
                           use_smooth = F, contribs = NULL){
  if("KEGG" %in% names(m_datasets[[2]])){
    setnames(m_datasets[[2]], "KEGG", "compound")
  } 
  mets_melt = melt(m_datasets[[2]], id.var = "compound", variable.name = "Sample")
  if(!simulated){
    m_network = build_metabolic_model(m_datasets[[1]], config_table)
    spec_cmp_scores = get_species_cmp_scores(species_table = m_network[[2]], m_network[[1]], normalize = T)
  }  else {
    m_network = build_metabolic_model(m_datasets[[1]], config_table, manual_agora = T)
    spec_cmp_scores = get_species_cmp_scores(species_table = m_network[[2]], m_network[[1]], normalize = T, manual_agora = T)
  }
  if("score_transform" %in% config_table[,V1]){
      spec_cmp_scores = transform_cmps(spec_cmp_scores, score_transform = config_table[V1=="score_transform", V2])
    }
  cmp_scores = spec_cmp_scores[,sum(CMP), by=list(Sample, compound)]

  if("met_transform" %in% config_table[,V1]){
    mets_melt = transform_mets(mets_melt, met_transform = config_table[V1=="met_transform", V2])
  }
  
  m_compare_mets = merge(cmp_scores, mets_melt, by=c("Sample", "compound"))
  
  #Separate out by synth only, deg only, both
  synth_deg = m_compare_mets[,list(sum(V1 > 0), sum(V1 < 0)), by=compound]
  setnames(synth_deg, c("V1", "V2"), c("SynthSamples", "DegSamples"))
  m_compare_mets = merge(m_compare_mets, synth_deg, by="compound")
  m_compare_mets[,CMPType:=ifelse(SynthSamples > 0 & DegSamples==0, "Synth", "Both")]
  m_compare_mets[,CMPType:=ifelse(DegSamples > 0 & SynthSamples==0, "Deg", CMPType)]
  setnames(m_compare_mets, c("V1", "value"), c("CMP", "Met"))
  
  # Run mimosa2
  m2_results = run_mimosa2(config_table, species = m_datasets[[1]], mets = m_datasets[[2]], compare_only = T)

  m_compare_mets = merge(m_compare_mets, m2_results[[2]], by = "compound", all = T)
  m_compare_mets[,m2_R:=sqrt(Rsq)]
  m_compare_mets[,m2_R:=ifelse(Slope < 0, -1*m2_R, m2_R)]
    
  m_compare_mets[,annotResult:=paste0(round(Rsq, 3), " ", round(PVal, 3))] #Just always compare with m1 I guess
  m_compare_mets[,M2Signif:=ifelse(PVal < 0.01, T, F)]
  #setnames(m_compare_mets, "PosNegCompare", "M1M2Positive")
  if(rank_based_plot){
      m_compare_mets[,MetValue:=Met]
      m_compare_mets[,CMPValue:=CMP]
      m_compare_mets[,Met:=NULL]
      m_compare_mets[,Met:=rank(MetValue), by=compound]
      m_compare_mets[,CMP:=NULL]
      m_compare_mets[,CMP:=rank(CMPValue), by=compound]
  }
  if(simulated){
    if(!is.null(contribs)){
      m_compare_mets = merge(m_compare_mets, contribs[Species=="Inflow"], by = "compound", all.x = T)
      setnames(m_compare_mets, "VarShare", "EnvVarShare")
      m_compare_mets[,EnvMet:=ifelse(EnvVarShare > 0.5, T, F)]
      #Compare contrib also
    } 
  }
  return(m_compare_mets)
}

cmp_met_compare_plot = function(m_compare_met_table, use_smooth = F, simulated = F, 
                                all_settings = c("orig", "logMets", "rankBased", "logMetsRank", "mblm")){
  bad_mets = m_compare_met_table[,var(CMP), by=compound][V1==0, compound]
  bad_mets = c(bad_mets, m_compare_met_table[,var(Met), by=compound][V1==0, compound])
  m_compare_met_table = m_compare_met_table[!compound %in% bad_mets]
  if(simulated){
    m_compare_met_table = m_compare_met_table[!is.na(EnvMet)]
  }
  n_comps = m_compare_met_table[,length(unique(compound))]
  #Set up orders
  m2_order = m_compare_met_table[Setting == "orig"][order(m2_R, decreasing = T), unique(compound)]
  m_compare_met_table[,compound:=factor(compound, levels = m2_order)]
  m_compare_met_table[,Setting:=factor(Setting, levels = all_settings)]
  annot_data = m_compare_met_table[,list(compound, Setting, annotResult)]
  setkey(annot_data, NULL)
  annot_data = unique(annot_data)
  if(simulated){
    m_compare_met_table[,EnvOutcome:=paste0(ifelse(EnvMet, "Env", "Micro"), "_", ifelse(M2Signif, "Signif", "Not"))]
    m_compare_met_table[,table(EnvOutcome)]
    met_plot = ggplot(m_compare_met_table[!compound %in% bad_mets], aes(x=CMP, y = Met, color = EnvOutcome)) + 
    geom_point(alpha = 0.6, size = 0.8, aes(shape = EnvMet)) + facet_wrap(~compound+Setting, scales = "free", ncol = length(all_settings), drop = F) + 
    theme(strip.text = element_text(size = 7), axis.text = element_text(size = 5)) + scale_color_brewer(palette = "Set1") + 
      geom_text(data = annot_data, aes(label = annotResult), x = Inf, y = Inf, hjust = 1, vjust = 1, inherit.aes = F, size = 2.5, fontface = "bold")
  } else {
    met_plot = ggplot(m_compare_met_table, aes(x=CMP, y = Met, color = M2Signif))+ 
      geom_point(alpha = 0.6, size = 0.8) + facet_wrap(~compound+Setting, scales = "free", ncol = length(all_settings), drop = F) + 
      theme(strip.text = element_text(size = 7), axis.text = element_text(size = 5))+ scale_color_brewer(palette = "Set1") +
      geom_text(data = annot_data, aes(label = annotResult), x = Inf, y = Inf, hjust = 1, vjust = 1, inherit.aes = F, size = 2.5, fontface = "bold")
  }
  if(use_smooth){
    met_plot = met_plot + geom_smooth(
      aes(linetype = M2Signif), span = 0.85) 
  } else {
    met_plot = met_plot + geom_abline(aes(slope = Slope, intercept = Intercept, linetype = M2Signif, color = M2Signif), size = 1.2, alpha = 0.7, color = "black")
  }
  return(met_plot)

}

# Test transformation options
test_fit_options = function(processed_data, met_transform_options = c("zscore", "logplus", "sqrt", ""), 
                            score_transform_options = c("zscore", "logplus", "sqrt", ""), 
                            reg_options = c("regular", "rankBased"), orig_config_table = NULL){
  if(is.null(orig_config_table)){
    orig_config_table = data.table(V1 = c("database", "genomeChoices","metType", "kegg_prefix", "data_prefix", "vsearch_path", "metagenome_format"), 
                                                V2 = c("Greengenes 13_5 or 13_8", "PICRUSt KO genomes and KEGG metabolic model", "KEGG Compound IDs", "data/KEGGfiles/", "data/", "bin/vsearch", F))
  }
  orig_results = run_mimosa2(orig_config_table, species = processed_data[[1]], mets = processed_data[[2]])
  compare_plot_list = vector("list", length(met_transform_options)*length(score_transform_options)*length(reg_options))
  compare_dat = data.table()
  var_share_compare_dat = data.table()
  var_count = 1
  for(met_option in met_transform_options){
    for(score_option in score_transform_options){
      for(reg_method in reg_options){
        config_table = rbind(orig_config_table, data.table(V1 = c("met_transform", "score_transform"), c(met_option, score_option))) # this is the same...
        if(reg_method == "rankBased"){
          config_table = rbind(orig_config_table,data.table(V1 = "rankBased", V2 = T))
        }
        print(met_option)
        print(score_option)
        print(reg_method)
        transforms_run = run_mimosa2(config_table, species = processed_data[[1]], mets = processed_data[[2]])
        trans_compare = compare_mimosa12(orig_results[[2]], transforms_run[[2]])
        trans_compare[,cor(Rsq.x, Rsq.y, use = "complete.obs")]
        plot1 = ggplot(trans_compare, aes(x=Rsq.x, y = Rsq.y)) + geom_point() + geom_abline(slope = 1, intercept = 0) + scale_x_continuous(limits = c(0, 0.5)) + scale_y_continuous(limits = c(0, 0.5)) + ggtitle(label = paste0(met_option, "_", score_option, "_", reg_method))
        compare_plot_list[[var_count]] = plot1
        trans_compare[,Setting:=paste0(met_option, "_", score_option, "_", reg_method)]
        compare_dat = rbind(compare_dat, trans_compare, fill = T)
        transforms_run$varShares[,Setting:=paste0(met_option, "_", score_option, "_", reg_method)]
        var_share_compare_dat = rbind(var_share_compare_dat, transforms_run$varShares)
        var_count = var_count + 1
      }
    }
  }
  return(list(plotList = compare_plot_list, compareRsq = compare_dat, compareVarShares = var_share_compare_dat))
}

