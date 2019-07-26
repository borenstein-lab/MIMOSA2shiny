#.libPaths(c("/data/shiny-server/r-packages/", "/data/shiny-server/R/x86_64-redhat-linux-gnu-library/3.2/")) #, "/data/shiny-server/app_specific_r_packages/"))
#library(Rcpp, lib.loc = "/data/shiny-server/r-packages/")
#library(Cairo, lib.loc = "/data/shiny-server/r-packages/")
library(shiny)
library(shinyjs)
#logjs(sessionInfo())

library(mimosa) #, lib.loc ="/data/shiny-server/r-packages/")
library(data.table) #, lib.loc ="/data/shiny-server/R/x86_64-redhat-linux-gnu-library/3.2/")
library(readr)
library(ggplot2) #, lib.loc = "/data/shiny-server/r-packages")
library(viridis) #, lib.loc = "/data/shiny-server/r-packages")
library(cowplot)
library(RColorBrewer)
options(datatable.webSafeMode = TRUE, scipen = 20000, stringsAsFactors = F, shiny.usecairo = F, shiny.maxRequestSize=200*1024^2, 
        show.error.locations=TRUE, shiny.trace = F)
theme_set(theme_cowplot() + theme(text = element_text(family = 'Helvetica')))
library(shinyBS)
library(ggpubr)

microbiome_data_upload = function(){
  fluidPage(
    tags$head(tags$style("#file1{color: gray } ")),

    tags$tr(class = "mainContainer", 
      tags$td(
      h4(get_text("microbiome_header"), width='100%'),
  #,tags$a(tags$img(src = "help.png", border = 0), href = "https://www.github.com/borenstein-lab/", target = "_blank"), width = '100%')
      p(get_text("microbiome_description")),
      fluidRow(
      column(radioButtons("database", get_text("database_title"), choices = get_text("database_choices"), width="100%"), width = 8),
      column(tipify(tags$img(src = "help.png", border = 0), title = get_text("microbiome_tooltip"), placement = "right"), width = 4)),
      fluidRow(
      column(fileInput("file1", get_text("microbiome_input_title"),
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv", 
                           ".tsv",
                           "text/tsv"), width = '430px'), width = 8),
        column(tags$a(tags$img(src = "example.png", border = 0), href = "test_seqs.txt", target = "_blank"), width = 4)),
   #   checkboxInput("metagenomeOpt", label=get_text("metagenome_option"), width = '100%'),
    #p(get_text("metagenome_description")),
    fluidRow(
    column(radioButtons("metagenome_format", get_text("metagenome_title"), choices = get_text("metagenome_options"), width="100%"), width = 8),
    column(tipify(tags$img(src = "help.png", border = 0), title = get_text("metagenome_tooltip"), placement = "right"), width = 4)),
   	fluidRow(
   		column(      fileInput("metagenome", get_text("metagenome_input_title"),
                         multiple = FALSE,
                         accept = c("text/csv",
                                    "text/comma-separated-values,text/plain",
                                    ".csv", 
                                    ".tsv",
                                    "text/tsv"), width = '430px'), width = 8), 
   		column(tags$a(tags$img(src = "example.png", border = 0), href = "test_metagenome.txt", target = "_blank"), width = 4))
      #) )#,
      # disabled(checkboxInput("metagenome_use", get_text("metagenome_use_option"),
      #               ))

      )
  ),
  id = "microbiome_section" , tags$style(type = "text/css", "#microbiome_section { horizontal-align: left; width: 100%}" ) #tags$style(".mainContainer { border: 2px black; }"), 
  )
}

metabolome_data_upload = function(){
  fluidPage(
    h4(get_text("metabolome_header")), 
    #tags$a(tags$img(src = "help.png", border = 0), href = "https://www.github.com/borenstein-lab/", target = "_blank"), id = "metabolome", width='100%')
    p(get_text("metabolome_description")),
    fluidRow(
      column(radioButtons("metType", label = get_text("met_type_title"), choices = get_text("met_type_choices"), selected = get_text("selected_met_type"), width="100%"), width = 8),
      column(tipify(tags$img(src = "help.png", border = 0), title = get_text("metabolome_tooltip"), placement = "right"), width = 4), width='100%'),
    fluidRow(
      column(checkboxInput("logTransform", get_text("metabolome_norm_description"), width="100%"), width = 8),
      column(tipify(tags$img(src = "help.png", border = 0), title = get_text("metabolome_transform_tooltip"), placement = "right"), width = 4)
    ),
	fluidRow(
   		column(fileInput("file2", get_text("metabolome_upload_title"),
              multiple = FALSE,
              accept = c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv", 
                         ".tsv",
                         "text/tsv"), width = '430px'), width = 8),
    	column(tags$a(tags$img(src = "example.png", border = 0), href = "test_mets.txt", target = "_blank"), width = 4)),
	id = "metabolome_section", tags$style(type = "text/css", "#metabolome_section { horizontal-align: left; width: 100%}" ))
}

network_settings = function(){
  fluidPage(
    h4(get_text("network_header"), width='100%'),
       #tags$a(tags$img(src = "help.png", border = 0), href = "https://www.github.com/borenstein-lab/", target = "_blank"), id = "genome", width='100%'),
    p(get_text("network_description")),
    fluidRow(
      column(radioButtons("genomeChoices", label = get_text("source_title"), choices = get_text("source_choices"), selected = get_text("source_choices")[2], width="100%"), width = 8),
      column(tipify(tags$img(src = "help.png", border = 0), title = get_text("network_tooltip"), placement = "right"), width = 4)
    ),
    fluidRow(
      column(numericInput("simThreshold", get_text("sim_title"), value = 0.99, min=0.8, max = 1, step = 0.01, width='430px'), width = 8), #, width = '100%'),
      column(tipify(tags$img(src = "help.png", border = 0), title = get_text("network_mapping_tooltip"), placement = "right"), width = 4)
    ),
    fluidRow(
   		column(fileInput("netAdd", get_text("net_mod_input_title"), multiple = FALSE, accept = c("text/csv",
                                 "text/comma-separated-values,text/plain",".csv", 
                                 ".tsv",
                                 "text/tsv"),  width = '550px'), width = 8), #,  #, width = '100%') #,
                                     column(tags$a(tags$img(src = "example.png", border = 0), href = "test_netAdd_species_rxns_KEGG_clean.txt", target = "_blank"), width = 4)
    )

    # radioButtons("modelTemplate", "Metabolic model template", choices = c("Generic KEGG metabolic model", "AGORA metabolic models (recommended)")),
    #disabled(radioButtons("closest", get_text("closest_title"), choices = get_text("closest_options"))),
    #checkboxInput("gapfill", get_text("gapfill_option"), width = '100%')
    , id = "network_section"
  )
}

algorithm_settings = function(){
  fluidPage(
    h4(get_text("algorithm_header"), id = "algorithm"),
    p(get_text("algorithm_description")),
    #radioButtons("contribType", get_text("stat_title"), choices = get_text("stat_choices")),
    fluidRow(
      column( radioButtons("regType", get_text("regression_title"), choices = get_text("regression_choices"), width="100%"), width = 8),
      column(tipify(tags$img(src = "help.png", border = 0), title = get_text("algorithm_tooltip"), placement = "right"), width = 4)
    ),
    fluidRow(
      column( checkboxInput("compare_only", get_text("skip_contribs_option"), width="100%"), width = 8),
      column(tipify(tags$img(src = "help.png", border = 0), title = get_text("skip_contribs_tooltip"), placement = "right"), width = 4)
    ),
    id = "algorithm_section"
  )
}

output_settings = function(){
  fluidPage(
    h4("Output", id = "outputResults"),
    #radioButtons("returnType", "Select desired output:", choices = c("Summarized contribution results", "Full verbose data output")),
    checkboxInput("netVis", "Generate interactive metabolic network visualization of results")
  )
}

###### Core app function
run_pipeline = function(input_data, configTable, analysisID){
  withProgress(message = "Running MIMOSA!", {
    if(is.character(input_data)){ #Error occurred
      stop(input_data) #Print error message
    }
    #process arguments
    species = input_data$species
    mets = input_data$mets
    print(configTable)
    configTable[,V2:=as.character(V2)]
    print(configTable)
	  configTable = check_config_table(configTable, app = T)
	
    incProgress(2/10, detail = "Building metabolic model")
    validate(need(!(configTable[V1=="database", V2==get_text("database_choices")[4]] & 
                                  configTable[V1=="genomeChoices", V2 != get_text("source_choices")[1]]), "Error: Only KEGG metabolic model (network option 1) can be used with KEGG Ortholog data"))
    network_results = build_metabolic_model(species, configTable, netAdd = input_data$netAdd) #, input_data$netAdd) #input_data$geneAdd, 
    network = network_results[[1]]
    species = network_results[[2]] #Allow for modifying this for AGORA
    if(!is.null(input_data$metagenome) & configTable[V1=="database", V2!=get_text("database_choices")[4]]){
      #If we are doing a comparison of the species network and the metagenome network
      #Metagenome data
      metagenome_network = build_metabolic_model(input_data$metagenome, configTable)
    }
    if(configTable[V1=="metType", V2 ==get_text("met_type_choices")[2]]){
      mets = map_to_kegg(mets)
    }
    incProgress(1/10, detail = "Calculating metabolic potential and fitting metabolite concentration model")
    #Get CMP scores
    if("rxnEdit" %in% configTable[,V1]){
      rxn_param = T
      cat("Will refine reaction network\n")
    } else rxn_param = F
    if("rankBased" %in% configTable[,V1]){
      rank_based = T
      cat("Will use rank-based/robust regression\n")
      if("rank_type" %in% configTable[,V1]){
        rank_type = configTable[V1=="rank_type", V2]
      } else {
        rank_type = "rfit"
      }
      cat(paste0("Regression type is ", rank_type, "\n"))
    } else rank_based = F
    if(configTable[V1=="database", V2==get_text("database_choices")[4]] & identical(configTable[V1=="metagenome_format", V2], get_text("metagenome_options")[1])){
      no_spec_param = T
      humann2_param = F
    } else if(configTable[V1=="database", V2==get_text("database_choices")[4]] & identical(configTable[V1=="metagenome_format", V2], get_text("metagenome_options")[2])){
      no_spec_param = F
      humann2_param = T
      cat("Humann2 format\n")
    } else {
      no_spec_param = F
      humann2_param = F
    }
    # if("revRxns" %in% configTable[,V1]){ #Whether to add reverse of reversible-annotated rxns - mainly for agora networks
    #   network = add_rev_rxns(network, sameID = T) # Give reverse the same rxn ID
    #   cat("Will add reverse of reversible reactions\n")
    # }
    if("met_transform" %in% configTable[,V1]){
      met_transform = configTable[V1=="met_transform", V2]
      cat(paste0("Will transform metabolite values, transform is ", met_transform))
    } else met_transform = ""

    #indiv_cmps = get_cmp_scores_kos(species, network) #Use KO abundances instead of species abundances to get cmps
    mets_melt = melt(mets, id.var = "compound", variable.name = "Sample")
    if(met_transform != ""){
      mets_melt = transform_mets(mets_melt, met_transform)
    }
    if(rxn_param){
      cmp_mods =  fit_cmp_net_edit(network, species, mets_melt, manual_agora = agora_param, rank_based = rank_based)
      network = cmp_mods[[3]] #Revised network
      indiv_cmps = cmp_mods[[4]]
      #Will have to report nice summary of rxns removed, rxns direction switched, etc
    } else {
      indiv_cmps = get_species_cmp_scores(species, network, normalize = !rxn_param, leave_rxns = rxn_param, manual_agora = F, kos_only = no_spec_param, humann2 = humann2_param)
      indiv_cmps = indiv_cmps[compound %in% mets[,compound]]
      cmp_mods = fit_cmp_mods(indiv_cmps, mets_melt, rank_based = rank_based, rank_type = rank_type)
    }
    print(cmp_mods[[1]]) 
    if(!configTable[V1 == "compare_only", V2==T]){
      incProgress(2/10, detail = "Calculating microbial contributions")
      #Get low-abundance species to remove
      if(!humann2_param){
        spec_dat = melt(species, id.var = "OTU", variable.name = "Sample")[,list(value/sum(value), OTU), by=Sample] #convert to relative abundance
        bad_spec = spec_dat[,list(length(V1[V1 != 0])/length(V1), max(V1)), by=OTU]
        bad_spec = bad_spec[V1 < 0.2 & V1 < 0.1, OTU] #Never higher than 10% and absent in at least 80% of samples
        print(bad_spec)
      } else bad_spec = NULL
      var_shares = calculate_var_shares(indiv_cmps, met_table = mets_melt, model_results = cmp_mods, config_table = configTable, species_merge = bad_spec)
      ## If nothing significant was analyzed, behave as if compare_only were selected
      if(is.null(var_shares)){
        configTable[V1 == "compare_only", V2:="TRUE"]
      }
      print("Var shares species")
      if(!is.null(var_shares)) print(var_shares[,unique(Species)])
    } else {
      var_shares = NULL
    }
    #Add species/rxn info
    cmp_summary = get_cmp_summary(species, network, normalize = !rxn_param, manual_agora = F, kos_only = no_spec_param, humann2 = humann2_param, 
                                  met_subset = cmp_mods[[1]][!is.na(Rsq) & Rsq != 0,compound], contrib_sizes = var_shares)
    cmp_mods[[1]] = merge(cmp_mods[[1]], cmp_summary, by = "compound", all.x = T)
    print(cmp_mods[[1]])
    cmp_mods[[1]][,compound:=as.character(compound)]
    indiv_cmps[,compound:=as.character(compound)]
    #shinyjs::logjs(devtools::session_info())
    #Order dataset for plotting
    incProgress(1/10, detail = "Making CMP-Metabolite plots")
    print(indiv_cmps[compound %in% mets_melt[,compound]])
    CMP_plots = plot_all_cmp_mets(cmp_table = indiv_cmps, met_table = mets_melt, mod_results = cmp_mods[[1]])
      
    if(configTable[V1 == "compare_only", V2 != TRUE]){
      incProgress(1/10, detail = "Making metabolite contribution plots")
      comp_list = var_shares[!is.na(VarShare), unique(as.character(compound))]
      comp_list = comp_list[!comp_list %in% var_shares[Species == "Residual" & VarShare == 1, as.character(compound)]]
      all_contrib_taxa = var_shares[compound %in% comp_list & !is.na(VarShare) & Species != "Residual" & Species != "Other" & abs(VarShare) >= 0.02, as.character(unique(Species))]
      print(all_contrib_taxa)
      getPalette = colorRampPalette(brewer.pal(12, "Paired"))
      contrib_color_palette = c("gray", getPalette(length(all_contrib_taxa))) #"black", 
      names(contrib_color_palette) = c( "Other", all_contrib_taxa) #"Residual",
      print("Color palette")
      print(contrib_color_palette)
      var_shares[,MetaboliteName:=met_names(as.character(compound))]
      var_shares[is.na(MetaboliteName), MetaboliteName:=compound]
      met_contrib_plots = lapply(comp_list, function(x){
        print(x)
        if(is.na(met_names(x))){
          met_id = x
        } else { met_id = met_names(x)}
        plot_contributions(var_shares, met_id, metIDcol = "MetaboliteName", color_palette = contrib_color_palette, include_residual = F, merge_threshold = 0.01)
      })
    } else {
      met_contrib_plots = NULL
    }
    print(analysisID)
    dir.create(path = analysisID, showWarnings = T)
    print(dir.exists(analysisID))
    for(i in 1:length(CMP_plots)){
      print(paste0(analysisID, "/", analysisID, "_", names(CMP_plots)[i], ".png"))
      print(CMP_plots[[i]])
      if(!identical(CMP_plots[[i]], NA)){
        save_plot(CMP_plots[[i]], file = paste0(analysisID, "/", analysisID, "_", names(CMP_plots)[i], ".png"), base_width = 2, base_height = 2)
      } 
    }
    if(!configTable[V1 == "compare_only", V2==T]){
      for(i in 1:length(met_contrib_plots)){
        print(comp_list[i])
        if(!is.null(met_contrib_plots[[i]])){
          save_plot(met_contrib_plots[[i]] + guides(fill = F), file = paste0(analysisID, "/", analysisID, "_", comp_list[i], "_contribs.png"), base_width = 2, base_height = 2)
        }
      }
      # if(!exists("contrib_legend")){ #Get legend from first non-nulll compound
      #   
      # }
      print("Making legend")
      print(all_contrib_taxa)
      print(contrib_color_palette)
      leg_dat = data.table(V1 = factor(names(contrib_color_palette), levels = c(all_contrib_taxa, "Other"))) #, "Residual"
      setnames(leg_dat, "V1", "Contributing Taxa")
      print(leg_dat)
      legend_plot = ggplot(leg_dat, aes(fill = `Contributing Taxa`, x=`Contributing Taxa`)) + geom_bar() + scale_fill_manual(values = contrib_color_palette, name = "Contributing Taxa")# + theme(legend.text = element_text(size = 10))
      contrib_legend = tryCatch(get_legend(legend_plot), error = function(){ return(NULL)}) 
      #save(contrib_legend, file = "data/exampleData/contrib_legend.rda")
      if(!is.null(contrib_legend)) save_plot(contrib_legend, file = paste0(analysisID, "/", analysisID, "_", "contribLegend.png"), dpi=120, base_width = 4, base_height = 2.5)
    } else {
      contrib_legend = NULL
    }
    if(configTable[V1=="genomeChoices", V2 == get_text("source_choices")[1]]){
      network_sub = network[Prod %in% cmp_mods[[1]][,compound]|Reac %in% cmp_mods[[1]][,compound]] #Network is in KEGG compounds
    } else {
      network[,KEGGReac:=agora_kegg_mets(Reac)]
      network[,KEGGProd:=agora_kegg_mets(Prod)]
      network_sub = network[(KEGGReac %in% cmp_mods[[1]][,compound] & grepl("[e]", Reac, fixed = T))|(KEGGProd %in% cmp_mods[[1]][,compound]& grepl("[e]", Prod, fixed = T))]
    }
    analysis_summary = get_analysis_summary(input_species = input_data[[1]], species = species, mets = mets, network = network, indiv_cmps = indiv_cmps, cmp_mods = cmp_mods, var_shares = var_shares, config_table = configTable)
    return(list(newSpecies = species, varShares = var_shares, modelData = cmp_mods[[1]], configs = configTable[!grepl("prefix", V1) & V1 != "vsearch_path"], 
                networkData = network_sub, CMPScores = indiv_cmps[CMP != 0], CMPplots = CMP_plots, metContribPlots = met_contrib_plots, plotLegend = contrib_legend, 
                analysisSummary = analysis_summary))
    #Send var_shares for download
    #Generate plot of var shares
    #source(other stuff)
  })
}

ui = fluidPage(
  useShinyjs(),
  tags$head(
    #tags$link(rel = "stylesheet", type = "text/css", href = "lab_styles.css"),
    #tags$link(rel = "stylesheet", type = "text/css", href = "shiny-alt.css"),
    #tags$link(rel = "stylesheet", type = "text/css", href = "software.css"),
    tags$link(rel = "stylesheet", type = "text/css", href = "burritostyle.css"), 
    tags$style("#errorMessage{color: red;
                                 font-size: 14px;
                                 }"
    )
  ),
  tags$img(src = "title_lab_800.jpg", border = 0), #Borenstein lab img
  titlePanel("MIMOSA2", windowTitle = "MIMOSA2"),
  mainPanel(fluidPage(
    p("MIMOSA2 is a tool for metabolic model-based evaluation of paired microbiome and metabolomics datasets. MIMOSA2 1) constructs community metabolic models, 
    2) assesses whether metabolite concentrations are consistent with estimated community metabolic potential, and 3) identifies specific taxa and reactions that can
    explain metabolite variation.
      For more information, see the ", 
      tags$a(tags$b("documentation."), href = "https://cnoecker.github.io/MIMOSA2shiny", target = "_blank"))), id="description", width = 12),
  tags$style(type = "text/css", "#title { color: #3CB371; horizontal-align: left; width: 100%}"  ),
  mainPanel(
    uiOutput("uploadPage"),
    tags$style(type = "text/css", "#title { horizontal-align: left; width: 100%}",  ".col-sm-12 {width: 95%}"  ),
    width = 12
  ))



server <- function(input, output, session) {
  #shinyjs::logjs(sessionInfo())
  #shinyjs::logjs(devtools::session_info())
  #shinyjs::logjs(Cairo.capabilities())
  analysisID = randomString()

  output$uploadPage = renderUI({
    if(is.null(input$goButton) & is.null(input$exampleButton)){
      return(fluidPage(
        #h3("Data Input"),
        microbiome_data_upload(),
        network_settings(),
        metabolome_data_upload(),
        algorithm_settings(),
        fluidRow(
          column(
            actionButton("goButton", " Run MIMOSA "),
            actionButton("exampleButton", "Show results for example dataset"),
             br(),
             br(),
            tags$style(type='text/css', "#goButton { vertical-align: middle; horizontal-align: middle; font-size: 22px; background-color: #3CBCDB; padding: 5px;}"), 
        	tags$style(type='text/css', "#exampleButton { vertical-align: middle; horizontal-align: middle; font-size: 14px; background-color: #C3D2D5}"), 
            width = 12, align = "center"
          )), width="100%", tags$style(type = "text/css", "#title { horizontal-align: left; width: 100%}")))
    } else if(input$goButton==0 & input$exampleButton == 0){
      return(fluidPage(
        #h3("Data Input"),
        microbiome_data_upload(),
        network_settings(),
        metabolome_data_upload(),
        #h3("Settings"),
        algorithm_settings(),
        #output_settings(),
        fluidRow(
          column(
            actionButton("goButton", " Run MIMOSA "),
             actionButton("exampleButton", "Show results for example dataset"),
             br(),
             br(),
            tags$style(type='text/css', "#goButton { vertical-align: middle; horizontal-align: middle; font-size: 22px; background-color: #3CBCDB; padding: 5px; margin:5px;}"), 
            tags$style(type='text/css', "#exampleButton { vertical-align: middle; horizontal-align: center; font-size: 14px; background-color: #C3D2D5; margin:5px;}"), 
            width = 12, align = "center"
          )), width="100%", tags$style(type = "text/css", "#title { horizontal-align: left; width: 100%}"  )))
    } else {
      if(input$compare_only == F & !(!is.null(input$metagenome) & input$metagenome_format==get_text("metagenome_options")[1])){
        if(class(datasetInput()) != "character"){
        fluidPage(
            fluidRow(
              column(
                downloadButton("downloadAll", "Download All Results Tables and Plots", class = "downloadButton"),
                downloadButton("downloadSummaryStats", "Analysis Summary Statistics", class = "downloadButton"),
                downloadButton("downloadModelData", "Model Summaries", class = "downloadButton"), 
                downloadButton("downloadData", "Contribution Results", class = "downloadButton"), 
                downloadButton("downloadSpecies", "Mapped Taxa Abundances", class = "downloadButton"),
                downloadButton("downloadNetworkData", "Community Metabolic Network Models", class = "downloadButton"), 
                downloadButton("downloadCMPs", "Community Metabolic Potential Scores", class = "downloadButton"),
                downloadButton("downloadSettings", "Record of Configuration Settings", class = "downloadButton"),
                tags$style(type='text/css', ".downloadButton { float: left; font-size: 14px; margin: 2px; margin-bottom: 3px; }"), width = 12, align = "center")),
            p(get_text("result_table_description")),
            fluidRow( # Big table
              DT::dataTableOutput("allMetaboliteInfo"), width="100%"
            ),
            fluidRow(
              plotOutput("contribLegend")
            ),
            fluidRow(
              downloadButton("downloadComparePlots", "Download selected CMP-Metabolite plots", class = "downloadButton"), 
              downloadButton("downloadContribPlots", "Download selected taxa contribution plots", class = "downloadButton"), 
              downloadButton("downloadContributionHeatmap", "Generate and download contribution heatmap plot", class = "downloadButton")
              #verbatimTextOutput('x4')
              #tableOutput("data2"), width = "100%"
            ), 
            tags$style(
              type="text/css",
              "#image img {max-width: 2in; max-height: 2in; width: auto; height: auto}"
            ), 
            tags$style(
              type="text/css",
              "#downloadAll  {background-color: #3CBCDB}"
            ) 
            

          # ),
          # fluidRow(
          #   plotOutput("contribPlots", height = "650px") #, click = "plot_click", hover = "plot_hover") 
          # ) 
        )
        } else {
          textOutput("errorMessage")
        }
        
      } else {
        fluidPage(
          fluidRow(
            column(
              downloadButton("downloadAll", "Download All Results Tables and Plots", class = "downloadButton"),
              downloadButton("downloadSummaryStats", "Analysis Summary Statistics", class = "downloadButton"),
              downloadButton("downloadModelData", "Model Summaries", class = "downloadButton"), 
              downloadButton("downloadSpecies", "Mapped Taxa Abundances", class = "downloadButton"),
              downloadButton("downloadNetworkData", "Community Metabolic Network Models", class = "downloadButton"), 
              downloadButton("downloadCMPs", "Community Metabolic Potential Scores", class = "downloadButton"),
              downloadButton("downloadSettings", "Record of Configuration Settings", class = "downloadButton"),
              tags$style(type='text/css', ".downloadButton { float: left; font-size: 14px; margin: 2px; margin-bottom: 3px; }"),
              tags$style(
                type="text/css",
                "#downloadAll  {background-color: #3CBCDB}"
              ), width = 12, align = "center")),
          p(get_text("result_table_description")),
          fluidRow( # Big table
            DT::dataTableOutput("allMetaboliteInfo"), width="100%"
          ),
          fluidRow(
            downloadButton("downloadComparePlots", "Download Selected CMP-Metabolite Plots", class = "downloadButton") #, 
            #downloadButton("downloadContribPlots", "Download Selected CMP-Metabolite Plots", class = "downloadButton"), 
            #downloadButton("downloadContributionHeatmap", "Generate and download contribution heatmap plot", class = "downloadButton")
            #verbatimTextOutput('x4')
            #tableOutput("data2"), width = "100%"
          )
          # fluidRow(
          #   plotOutput("contribPlots", height = "650px") #, click = "plot_click", hover = "plot_hover") 
          # ) 
        )
      }
    }

  })
  observeEvent(input$database, {
    if(input$database==get_text("database_choices")[1]){
      #updateRadioButtons(session, "genomeChoices", selected = get_text("source_choices")[2])
    } else {
      #updateRadioButtons(session, "genomeChoices", selected = get_text("source_choices")[1])
    }
    if(input$database==get_text("database_choices")[4]){
      #updateCheckboxInput(session, "metagenomeOpt", value = T)
      disable("file1")
    } else {
      #updateCheckboxInput(session, "metagenomeOpt", value = F)
      enable("file1")
    }
  })
  observeEvent(input$genomeChoices, {
    # if(input$genomeChoices==get_text("source_choices")[2] & input$database==get_text("database_choices")[1]){ #Only enable mapping threshold if providing sequences
    #   enable("simThreshold")
    # } else {
    #   disable("simThreshold")
    #   }
  })

  config_table = reactive({
  	if(is.null(input$exampleButton)|input$exampleButton==0){
  	    initial_inputs = c("database", "metagenome_format", "genomeChoices",  #"geneAdd", 
                       "metType", "simThreshold", "logTransform", "regType", "compare_only") #Check that this is all of them
  	  if(is.null(input$metagenome)){
  	    initial_inputs = initial_inputs[initial_inputs != "metagenome_format"]
  	  }
  	  if(input$database != get_text("database_choices")[1]){
  	    initial_inputs = initial_inputs[initial_inputs != "simThreshold"] #Does nothing if we don't have ASV data
  	  } 
      inputs_provided = initial_inputs[which(sapply(initial_inputs, function(x){ !is.null(input[[x]])}))]
      values_provided = sapply(inputs_provided, function(x){ return(input[[x]])})
  	  inputs_provided = c("analysisID", inputs_provided, "kegg_prefix", "data_prefix", "vsearch_path")
	    values_provided = c(analysisID, values_provided, "data/KEGGfiles/", "data/", "bin/vsearch")    #print(file.exists(input$file1$datapath))
      print(inputs_provided)
      print(values_provided)
      if(input$logTransform == T){
        inputs_provided[inputs_provided == "logTransform"] = "met_transform"
        values_provided[inputs_provided == "met_transform"] = "logplus"
      } 
      if(input$regType == get_text("regression_choices")[1]){
        inputs_provided = c(inputs_provided, "rankBased", "rank_type")
        values_provided = c(values_provided, T, "rfit")
      }
      return(data.table(V1 = inputs_provided, V2 = values_provided))
  	} else {
  		return(fread("data/exampleData/configs_example.txt", header = T))
  	}
  })
  
  datasetInput = reactive({
    if(is.null(input$exampleButton)|input$exampleButton==0){
      validate(need(!is.null(input$file1)|!is.null(input$metagenome), "No microbiome file provided"), 
               need(input$file2, "No metabolite file provided")
               )
      file_list1 = list(input$file1, input$file2, input$metagenome, input$netAdd) # input$geneAddFile,
      names(file_list1) = c("file1","file2", "metagenome","netAdd") # "geneAddFile", 
      #logjs(config_table())
      print(file_list1)
      input_data = tryCatch(read_mimosa2_files(file_list = file_list1, configTable = config_table()), error = function(e){ return(e$message)})
      tryCatch(run_pipeline(input_data, config_table(), analysisID), error=function(e){ 
        print(e$call)
        return(e$message)})
    } else {
    	#logjs("Reading example data")
    	config = config_table()
    	species_dat = fread("data/exampleData/speciesExample.txt")
    	var_shares = fread("data/exampleData/contributionResultsExample.txt")
    	modResults = fread("data/exampleData/modelResultsExample.txt")
    	networkData = fread("data/exampleData/communityNetworkModelsExample.txt")
    	cmpScores = fread("data/exampleData/CMPScoresExample.txt")
    	summaryStats = fread("data/exampleData/StatsExample.txt")
    	#analysisID = "example"
    	
    	return(list(newSpecies = species_dat, varShares = var_shares, modelData = modResults, configs = config, networkData = networkData, 
    	            CMPScores = cmpScores, analysisSummary = summaryStats, example_data = T))
    }
  })
  
  
  
  observeEvent(input$goButton, {
    req(input$file2)
    #if(error_checks){
    config_table()
    datasetInput()
    enable("downloadSettings")
    enable("downloadData")
    enable("downloadModelData")
    enable("downloadNetworkData")
    #var_shares = run_pipeline(input) #Just give it everything and go from there
    
    #}
  })

  observeEvent(input$exampleButton, {
    config_table()
    datasetInput()
    enable("downloadSettings")
    enable("downloadData")
    enable("downloadModelData")
    enable("downloadNetworkData")
  })
  
  output$errorMessage = renderText({
    paste0("Error: ", datasetInput())
  })
  
  output$downloadAll = downloadHandler(
    filename = function(){
      if("example_data" %in% names(datasetInput())){
        "example_allResults.zip"
      } else {
        paste0(analysisID, "_allResults.zip")
      }
    },
    content = function(file){
      if("example_data" %in% names(datasetInput())){
        file_ids = paste0("data/exampleData/", list.files(path = "data/exampleData"))
      } else {
        write.table(datasetInput()$analysisSummary, file = paste0(analysisID, "/summaryStats.txt"), row.names = F, sep = "\t", quote=F)
        write.table(datasetInput()$configs, file = paste0(analysisID, "/configSettings.txt"), row.names = F, sep = "\t", quote=F)
        
        if(!input$compare_only & !(!is.null(input$metagenome) & input$metagenome_format==get_text("metagenome_options")[1])) write.table(datasetInput()$varShares, file = paste0(analysisID, "/contributionResults.txt"), row.names = F, sep = "\t", quote=F)
        write.table(datasetInput()$newSpecies, file = paste0(analysisID, "/mappedTaxaData.txt"), row.names = F, sep = "\t", quote=F)
        write.table(datasetInput()$modelData, file = paste0(analysisID, "/modelResults.txt"), row.names = F, sep = "\t", quote=F)
        write.table(datasetInput()$networkData, file = paste0(analysisID, "/communityNetworkModels.txt"), row.names = F, sep = "\t", quote=F)
        write.table(datasetInput()$CMPScores, file = paste0(analysisID, "/communityMetabolicPotentialScores.txt"), row.names = F, sep = "\t", quote=F)
        
        if(!input$compare_only & !(!is.null(input$metagenome) & input$metagenome_format==get_text("metagenome_options")[1])){
          plotData = datasetInput()$varShares
          
          tableData = datasetInput()$modelData[!is.na(Slope)]
          # plotData = merge(plotData, tableData, by="compound", all.x = T)
          # if("Slope.x" %in% names(plotData)) setnames(plotData, "Slope.x", "Slope")
          save_plot(plot_summary_contributions(plotData, include_zeros = T, remove_resid_rescale = F), filename = paste0(analysisID, "/contributionHeatmapPlotSelected.pdf"), 
                    base_width = 10, base_height = 8)
        }
        print(list.files(path = analysisID))
        file_ids = paste0(analysisID, "/", list.files(path = analysisID))
      }
      print(file_ids)

      zip(zipfile = file, files = file_ids)
    }
  )
  output$downloadSummaryStats <- downloadHandler(
    filename = function() {
      "summaryStats.txt"
      #paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      write.table(datasetInput()$analysisSummary, file, row.names = F, sep = "\t", quote=F)
    }
  )
  output$downloadSettings <- downloadHandler(
    filename = function() {
      "configSettings.txt"
      #paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      write.table(datasetInput()$configs, file, row.names = F, sep = "\t", quote=F)
    }
  )
  output$downloadData <- downloadHandler(
    filename = function() {
      "contributionResults.txt"
      #paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      write.table(datasetInput()$varShares, file, row.names = F, sep = "\t", quote=F)
    }
  )
  output$downloadSpecies = downloadHandler(
    filename = function() {
      "mappedTaxaData.txt"
    },
    content = function(file){
      write.table(datasetInput()$newSpecies, file, row.names = F, sep = "\t", quote=F)
    }
  )
  output$downloadModelData <- downloadHandler(
    filename = function() {
      "modelResults.txt"
      #paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      write.table(datasetInput()$modelData, file, row.names = F, sep = "\t", quote=F)
    }
  )
  output$downloadNetworkData = downloadHandler(
  	filename = function(){
  		"communityNetworkModels.txt"
  	},
  	content = function(file) {
  		write.table(datasetInput()$networkData, file, row.names = F, sep = "\t", quote=F)
  	}
  )
  output$downloadCMPs = downloadHandler(
    filename = function(){
      "communityMetabolicPotentialScores.txt"
    },
    content = function(file) {
      write.table(datasetInput()$CMPScores, file, row.names = F, sep = "\t", quote=F)
    }
  )
  
  output$contribPlots = renderPlot({
    plotData = datasetInput()$varShares
    # plotData = merge(plotData, datasetInput()$modelData, by="compound", all.x = T)
    # if("Slope.x" %in% names(plotData)) setnames(plotData, "Slope.x", "Slope")
    #print(plotData)
    ##met1 = plotData[1,compound]
    #### Make a drop bar to select one metabolite to display plot & data for at a time!!! cool.
    #plotData[,hist(V1)]
    return(plot_summary_contributions(plotData, include_zeros = T, remove_resid_rescale = F))
  }, res = 100)
  
  output$allMetaboliteInfo = DT::renderDT({
    #print(datasetInput())
    print(datasetInput()$modelData)
    print(names(datasetInput()$modelData))
    tableData = datasetInput()$modelData[!is.na(Slope)]
    #Get order before rounding
    tableData[,m2R:=ifelse(Slope < 0, -1*sqrt(Rsq), sqrt(Rsq))]
    compound_order = tableData[order(m2R, decreasing = T), compound] #, tableData[Slope <= 0][order(Rsq, decreasing = T), compound])
    
    print(head(tableData[,Slope]))
    #tableData[,m2R:=ifelse(Slope < 0, -1*sqrt(Rsq), sqrt(Rsq))]
    tableData[,PVal:=round(PVal, 5)]
    tableData[,Rsq:=round(Rsq, 3)]
    tableData[,Slope:=round(Slope, 3)]
    tableData[,Intercept:=round(Intercept, 3)]
    tableData[,metName:=sapply(compound, met_names)]
    tableData2 = tableData[,list(compound, metName, Rsq, PVal, Slope, TopSynthSpecGenes, TopDegSpecGenes, Intercept)]
    tableData2[,compound:=factor(compound, levels = compound_order)]
    print(tableData2)
    #good_comps = list.files(path = getwd(), pattern = analysisID)
    if("example_data" %in% names(datasetInput())){
      analysisID2 = "data/exampleData/example"
    } else {
      analysisID2 = paste0(analysisID, "/", analysisID)
    }
    print(analysisID2)
    tableData2[,Plot:=sapply(paste0(analysisID2, "_", compound, ".png"), function(x){ 
      if(file.exists(x)){
        return(img_uri(x))
      } else {
        return(img_uri("blank_plot.png"))
      }})]
    if(datasetInput()$configs[V1=="compare_only", V2==FALSE]){
      #tableData2 = tableData2[paste0(analysisID, "_",compound,"_contribs.png") %in% good_comps]
      tableData2[,ContribPlot:=sapply(paste0(analysisID2, "_", compound, "_contribs.png"), function(x){
        if(file.exists(x)){
          return(img_uri(x))
        } else {
          return(img_uri("blank_plot.png"))
        }
        })]
      tableData2 = tableData2[,list(compound, metName, Rsq, PVal, Slope, Plot, ContribPlot, TopSynthSpecGenes, TopDegSpecGenes, Intercept)]
      setnames(tableData2, c("Compound ID", "Name", "R-squared", "P-value", "Slope", "Comparison Plot", "Contribution Plot", "Top Producing Taxa and Genes/Rxns", "Top Utilizing Taxa and Genes/Rxns",  "Intercept"))
      tooltip_table = htmltools::withTags(table(
        class = 'display',
        thead(
          tr(
            th('', title = get_text("results_table_titles")[1]),
            th('Compound ID', title = get_text("results_table_titles")[2]),
            th('Name', title = get_text("results_table_titles")[3]),
            th('R-squared', title = get_text("results_table_titles")[4]),
            th('P-value', title = get_text("results_table_titles")[5]),
            th('Slope', title = get_text("results_table_titles")[6]),
            th('Comparison Plot', title = get_text("results_table_titles")[7]),
            th('Contribution Plot', title = get_text("results_table_titles")[8]),
            th("Top Producing Taxa and Genes/Rxns", title = get_text("results_table_titles")[9]),
            th("Top Utilizing Taxa and Genes/Rxns", title = get_text("results_table_titles")[11]),
            th("Intercept", title = get_text("results_table_titles")[13])
          ))))
    } else {
      if(datasetInput()$configs[V1=="metagenome_format", V2==get_text("metagenome_options")[1]]){ #Skip species
        tableData2 = tableData2[,list(compound, metName, Rsq, PVal, Slope, Plot, TopSynthSpecGenes, TopDegSpecGenes, Intercept)]
        setnames(tableData2, c("Compound ID", "Name", "R-squared", "P-value", "Slope", "Comparison Plot", "Top Producing Genes/Rxns","Top Utilizing Genes/Rxns", "Intercept"))
        tooltip_table = htmltools::withTags(table(
          class = 'display',
          thead(
            tr(
              th('', title = get_text("results_table_titles")[1]),
              th('Compound ID', get_text("results_table_titles")[2]),
              th('Name', title = get_text("results_table_titles")[3]),
              th('R-squared', title = get_text("results_table_titles")[4]),
              th('P-value', title = get_text("results_table_titles")[5]),
              th('Slope', title = get_text("results_table_titles")[6]),
              th('Comparison Plot', title = get_text("results_table_titles")[7]),
              th("Top Producing Genes/Rxns", title = get_text("results_table_titles")[9]),
              th("Top Utilizing Genes/Rxns", title = get_text("results_table_titles")[11]),
              th("Intercept", title = get_text("results_table_titles")[13])
            ))))
      } else {
        tableData2 = tableData2[,list(compound, metName, Rsq, PVal, Slope, Plot, TopSynthSpecGenes, TopDegSpecGenes, Intercept)]
        setnames(tableData2, c("Compound ID", "Name", "R-squared", "P-value", "Slope", "Comparison Plot", "Top Producing Taxa and Genes/Rxns", "Top Utilizing Taxa and Genes/Rxns", "Intercept"))
        tooltip_table = htmltools::withTags(table(
          class = 'display',
          thead(
            tr(
              th('', title = get_text("results_table_titles")[1]),
              th('Compound ID', get_text("results_table_titles")[2]),
              th('Name', title = get_text("results_table_titles")[3]),
              th('R-squared', title = get_text("results_table_titles")[4]),
              th('P-value', title = get_text("results_table_titles")[5]),
              th('Slope', title = get_text("results_table_titles")[6]),
              th('Comparison Plot', title = get_text("results_table_titles")[7]),
              th("Top Producing Taxa and Genes/Rxns", title = get_text("results_table_titles")[9]),
              th("Top Utilizing Taxa and Genes/Rxns", title = get_text("results_table_titles")[11]),
              th("Intercept", title = get_text("results_table_titles")[13])
            ))))
      }
    }

    return(DT::datatable(
      tableData2[order(`Compound ID`)], escape = F, options = list(lengthMenu = c(5, 10), pageLength = 5, rowCallback = DT::JS("function(r,d) {$(r).attr('overflow', 'hidden').attr('height', '217px')}")), filter = "top", 
      container = tooltip_table))
    #return(DT::datatable(tableData[order(m2R, decreasing = T),list(compound, Rsq, PVal, Slope, Intercept)], 
     #                    options = list(lengthMenu = c(5, 10), pageLength = 5)))
  })
  output$downloadComparePlots = downloadHandler(
    filename = function() {
      "CMP_Metabolite_Compare.zip"
      #paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      comp_num = input$allMetaboliteInfo_rows_selected
      tableData = datasetInput()$modelData[!is.na(Slope)]
      #Get order before rounding
      tableData[,m2R:=ifelse(Slope < 0, -1*sqrt(Rsq), sqrt(Rsq))]
      compound_order = tableData[order(m2R, decreasing = T), compound]
      #compound_order = c(tableData[Slope > 0][order(Rsq, decreasing = T), compound], tableData[Slope <= 0][order(Rsq, decreasing = T), compound])
      comp_ids = compound_order[comp_num]
      print(comp_ids)
      print(analysisID)
      if("example_data" %in% names(datasetInput())){
        analysisID2 = "data/exampleData/example"
      } else {
        analysisID2 = paste0(analysisID, "/", analysisID)
      }
      print(analysisID2)
      file_ids = paste0(analysisID2, "_", comp_ids, ".png")
      file_ids = file_ids[sapply(file_ids, file.exists)]
      zip(zipfile = file, files = file_ids)
    }
  )
  output$downloadContribPlots = downloadHandler(
    filename = function() {
      "Metabolite_Taxa_Contributions.zip"
      #paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      comp_num = input$allMetaboliteInfo_rows_selected
      tableData = datasetInput()$modelData[!is.na(Slope)]
      #Get order before rounding
      tableData[,m2R:=ifelse(Slope < 0, -1*sqrt(Rsq), sqrt(Rsq))]
      compound_order = tableData[order(m2R, decreasing = T), compound]
      #compound_order = c(tableData[Slope > 0][order(Rsq, decreasing = T), compound], tableData[Slope <= 0][order(Rsq, decreasing = T), compound])
      comp_ids = compound_order[comp_num]
      print(comp_ids)
      if("example_data" %in% names(datasetInput())){
        analysisID2 = "data/exampleData/example"
      } else {
        analysisID2 = paste0(analysisID, "/", analysisID)
      }
      file_ids = paste0(analysisID2, "_", comp_ids, "_contribs.png")
      file_ids = file_ids[sapply(file_ids, file.exists)]
      if(length(file_ids) > 0) zip(zipfile = file, files = file_ids)
    }
  )
  
  output$downloadContributionHeatmap = downloadHandler(
    filename = function() {
      "contributionHeatmapPlotSelected.pdf"
    },
    content = function(file) {
      comp_num = input$allMetaboliteInfo_rows_selected
      tableData = datasetInput()$modelData[!is.na(Slope)]
      #Get order before rounding
      tableData[,m2R:=ifelse(Slope < 0, -1*sqrt(Rsq), sqrt(Rsq))]
      compound_order = tableData[order(m2R, decreasing = T), compound]
      #compound_order = c(tableData[Slope > 0][order(Rsq, decreasing = T), compound], tableData[Slope <= 0][order(Rsq, decreasing = T), compound])
      comp_list = compound_order[comp_num]
      plotData = datasetInput()$varShares[compound %in% comp_list]
      #plotData = merge(plotData, tableData[compound %in% comp_list], by="compound", all.x = T)
      #print(plotData)
      ##met1 = plotData[1,compound]
      #### Make a drop bar to select one metabolite to display plot & data for at a time!!! cool.
      #plotData[,hist(V1)]
      save_plot(plot_summary_contributions(plotData, include_zeros = T, remove_resid_rescale = F), filename = file, base_width = 10, base_height = 8)
      
    }

  )
  output$x4 = renderPrint({
    s = input$allMetaboliteInfo_rows_selected
    if (length(s)) {
      cat('These rows were selected:\n\n')
      cat(s, sep = ', ')
    }
  })
  # output$data2 <- renderTable({
  #   dat <- cars[1:5,]
  #   dat$test <- c("a","b","c","d",
  #                 '<div id="testPlot2" class="shiny-plot-output" style="width: 100px ; height: 100px"></div>')
  #   dat
  # }, sanitize.text.function = function(x) x)
  # output$testPlot2 = renderPlot({
  #   qplot(x = 1:10, y = 2:11, geom = "point")
  # })
  output$indivCellInfo = renderDataTable({
    plotData = datasetInput()$varShares[,list(metID, compound, Species, Var, VarShare)]
    setnames(plotData, c("Metabolite", "KEGG ID", "Taxon", "Total Variance", "Variance Contribution"))
    spec_order = plotData[Taxon != "Residual",length(`Variance Contribution`[abs(`Variance Contribution`) > 0.05]), by=Taxon][order(V1, decreasing = T), Taxon]
    spec_order = c("Residual", spec_order)
    
    if(!is.null(input$plot_click)){
      met_of_interest = plotData[,levels(Metabolite)][round(input$plot_click$x)] 
      print(met_of_interest)
      spec_of_interest = spec_order[round(input$plot_click$y)]
      if(length(spec_of_interest)==1 & length(met_of_interest)==1){
        return(plotData[Taxon==spec_of_interest & Metabolite==met_of_interest])
      }
    } else {
      return(plotData[1])
    }
  })
  met_of_interest = reactive({
  	 plotData = datasetInput()$varShares
  	 if(is.null(input$plot_click$x)){
  	 	return(plotData[1, metID])
  	 } else {
  	 	return(plotData[,levels(metID)][round(input$plot_click$x)]) 
  	 }
  })
  output$indivPlots = renderPlot({
    plotData = datasetInput()$varShares
      met_of_interest = plotData[1,compound]      
   	return(plot_contributions(plotData, metabolite = met_of_interest, include_zeros = F))
  }, res = 100)
  
  output$contribLegend = renderPlot({
    if("example_data" %in% names(datasetInput())){
      load("data/exampleData/example_contrib_legend.rda")
      return(plot_grid(contrib_legend))
    } else {
      if(!is.null(datasetInput()$plotLegend)){
        return(plot_grid(datasetInput()$plotLegend))
      } else { return(NULL)}
    }
    #return(list(src = paste0(analysisID, "_contribLegend.png")))
  })
  
  session$onSessionEnded(function() {
    #remove plots
    plots_made = list.files(path = analysisID)
    print(plots_made)
    file.remove(paste0(analysisID, "/", plots_made))
    file.remove(analysisID)
    #Move files to www dir instead

  })
  
  
}

shinyApp(ui = ui, server = server)
