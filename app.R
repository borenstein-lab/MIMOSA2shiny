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
options(datatable.webSafeMode = TRUE, scipen = 20000, stringsAsFactors = F, shiny.usecairo = F, shiny.maxRequestSize=30*1024^2)
theme_set(theme_get() + theme(text = element_text(family = 'Helvetica')))
library(shinyBS)

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
                           ".csv"), width = '430px'), width = 8),
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
                                    ".csv"), width = '430px'), width = 8), 
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
                         ".csv"), width = '430px'), width = 8),
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
                                 "text/comma-separated-values,text/plain",".csv"),  width = '550px'), width = 8), #,  #, width = '100%') #,
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
    #process arguments
    species = input_data$species
    mets = input_data$mets
    print(configTable)
    print(configTable[V1=="file1", V2])
    configTable[,V2:=as.character(V2)]
    print(file.exists(configTable[V1=="file1", V2]))
    print(configTable)
    print(configTable[V1=="file1", V2])
	  configTable = check_config_table(configTable, app = T)
	
    incProgress(2/10, detail = "Building metabolic model")
    network_results = build_metabolic_model(species, configTable) #, input_data$netAdd) #input_data$geneAdd, 
    network = network_results[[1]]
    species = network_results[[2]] #Allow for modifying this for AGORA
    if(!is.null(input_data$metagenome) & configTable[V1=="database", V2==get_text("database_choices")[4]]){
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
    if("revRxns" %in% configTable[,V1]){ #Whether to add reverse of reversible-annotated rxns - mainly for agora networks
      network = add_rev_rxns(network, sameID = T) # Give reverse the same rxn ID
      cat("Will add reverse of reversible reactions\n")
    }
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
    if(!configTable[V1 == "compare_only", V2==T]){
      incProgress(2/10, detail = "Calculating microbial contributions")
      var_shares = calculate_var_shares(indiv_cmps, met_table = mets_melt, model_results = cmp_mods, config_table = configTable)
    } else {
      var_shares = NULL
    }
    #shinyjs::logjs(devtools::session_info())
    #Order dataset for plotting
    incProgress(1/10, detail = "Making CMP-Metabolite plots")
    CMP_plots = plot_all_cmp_mets(cmp_table = indiv_cmps, met_table = mets_melt, mod_results = cmp_mods[[1]])
    comp_list = cmp_mods[[1]][!identical(Rsq, NA) &  Rsq != 0][order(PVal), compound]
    if(!configTable[V1 == "compare_only", V2==T]){
      incProgress(1/10, detail = "Making metabolite contribution plots")
      met_contrib_plots = lapply(comp_list, function(x){
        plot_contributions(var_shares, x, metIDcol = "compound")
      })
    } else {
      met_contrib_plots = NULL
    }
    for(i in 1:length(CMP_plots)){
      save_plot(CMP_plots[[i]], file = paste0(analysisID, "_", names(CMP_plots)[i], ".png"), dpi = 50, base_width = 2, base_height = 2)
    }
    for(i in 1:length(met_contrib_plots)){
      print(comp_list[i])
      save_plot(met_contrib_plots[[i]], file = paste0(analysisID, "_", comp_list[i], "_contribs.png"), dpi = 50, base_width = 2, base_height = 2)
    }
    return(list(varShares = var_shares, modelData = cmp_mods[[1]], configs = configTable[!grepl("prefix", V1)], networkData = network, CMPplots = CMP_plots, metContribPlots = met_contrib_plots))
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
    tags$link(rel = "stylesheet", type = "text/css", href = "burritostyle.css")
  ),
  tags$img(src = "title_lab_800.jpg", border = 0), #Borenstein lab img
  titlePanel("MIMOSA", tags$style(type = "text/css", "horizontal-align: center;")),
  mainPanel(fluidPage(
    p("MIMOSA is a tool for metabolic model-based evaluation of paired microbiome and metabolomics datasets. For more information, see the ", 
      tags$a("documentation.", href = "https://cnoecker.github.io/MIMOSA2shiny", target = "_blank"))), id="description", width = 12),
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
      fluidPage(
      fluidRow(
        column(
          downloadButton("downloadModelData", "Model Summaries", class = "downloadButton"), 
          downloadButton("downloadData", "Contribution Results", class = "downloadButton"), 
          downloadButton("downloadNetworkData", "Community Metabolic Network Models", class = "downloadButton"), 
      downloadButton("downloadSettings", "Record of Configuration Settings", class = "downloadButton"),
      tags$style(type='text/css', ".downloadButton { float: left; font-size: 14px; margin: 2px; margin-bottom: 3px; }"), width = 12, align = "center")),
      p(get_text("result_table_description")),
      fluidRow( # Big table
        DT::dataTableOutput("allMetaboliteInfo"), width="100%"
      ),
      fluidRow(
        verbatimTextOutput('x4')
        #tableOutput("data2"), width = "100%"
      ),
      fluidRow(
        #plots
        #plotOutput("residPlot", height = "100px", click = "plot_click"),
        plotOutput("contribPlots", height = "650px") #, click = "plot_click", hover = "plot_hover") 
      ) #,
      # fluidRow(
      #   dataTableOutput("indivCellInfo"), width="100%"
      # ),
      # fluidRow(
      #   plotOutput("indivPlots", height = "600px")
      # )
      )
  
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
      inputs_provided = initial_inputs[which(sapply(initial_inputs, function(x){ !is.null(input[[x]])}))]
      values_provided = sapply(inputs_provided, function(x){ return(input[[x]])})
  	  inputs_provided = c(inputs_provided, "kegg_prefix", "data_prefix", "vsearch_path")
	    values_provided = c(values_provided, "data/KEGGfiles/", "data/", "bin/vsearch")    #print(file.exists(input$file1$datapath))
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
  		return(fread("data/exampleData/configs_example_clean.txt", header = T))
  	}
  })
  
  datasetInput <- reactive({
    if(is.null(input$exampleButton)|input$exampleButton==0){
      if(is.null(input$file1) & is.null(input$metagenome)) stop("No microbiome data file provided")
      if(is.null(input$file2)) stop("No metabolite data file provided")
      file_list1 = list(input$file1, input$file2, input$metagenome, input$netAddFile) # input$geneAddFile,
      names(file_list1) = c("file1","file2", "metagenome","netAdd") # "geneAddFile", 
      #logjs(config_table())
      input_data = read_mimosa2_files(file_list = file_list1, configTable = config_table())
      run_pipeline(input_data, config_table(), analysisID)
    } else {
    	#logjs("Reading example data")
    	config = config_table()
    	varShares = fread("data/exampleData/varSharesExample.txt")
    	modResults = fread("data/exampleData/modelResultsExample.txt")
    	networkData = fread("data/exampleData/communityNetworkModelsExample.txt")
    	return(list(varShares = varShares, modelData = modResults, configs = config))
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
  output$contribPlots = renderPlot({
    plotData = datasetInput()$varShares
    plotData = merge(plotData, datasetInput()$modelData, by="compound", all.x = T)
    #print(plotData)
    ##met1 = plotData[1,compound]
    #### Make a drop bar to select one metabolite to display plot & data for at a time!!! cool.
    #plotData[,hist(V1)]
    return(plot_summary_contributions(plotData, include_zeros = T, remove_resid_rescale = F))
  }, res = 100)
  
  output$allMetaboliteInfo = DT::renderDT({
    #print(datasetInput())
    tableData = datasetInput()$modelData[!is.na(Slope)]
    tableData[,m2R:=ifelse(Slope < 0, -1*sqrt(Rsq), sqrt(Rsq))]
    tableData[,PVal:=round(PVal, 3)]
    tableData[,Rsq:=round(Rsq, 3)]
    tableData[,Slope:=round(Slope, 3)]
    tableData[,Intercept:=round(Intercept, 3)]
    tableData[,metName:=met_names(compound)]
    compound_order = tableData[order(m2R, decreasing = T), compound]
    tableData2 = tableData[,list(compound, metName, Rsq, PVal, Slope, Intercept)]
    tableData2[,compound:=factor(compound, levels = compound_order)]
    good_comps = list.files(path = getwd(), pattern = analysisID)
    tableData2 = tableData2[paste0(analysisID, "_", compound, ".png") %in% good_comps]
    tableData2[,Plot:=sapply(paste0(analysisID, "_", compound, ".png"), img_uri)]
    if(input$compare_only != T){
      tableData2 = tableData2[paste0(analysisID, "_",compound,"_contribs.png") %in% good_comps]
      tableData2[,ContribPlot:=sapply(paste0(analysisID, "_", compound, "_contribs.png"), img_uri)]
      setnames(tableData2, c("Compound ID", "Name", "R-squared", "P-value", "Slope", "Intercept", "Comparison Plot", "Contribution Plot"))
    } else {
      setnames(tableData2, c("Compound ID", "Name", "R-squared", "P-value", "Slope", "Intercept", "Comparison Plot"))
    }
    print(tableData2[1])
    
    return(DT::datatable(tableData2[order(`Compound ID`)], escape = F, options = list(lengthMenu = c(5, 10), pageLength = 5)))
    #return(DT::datatable(tableData[order(m2R, decreasing = T),list(compound, Rsq, PVal, Slope, Intercept)], 
     #                    options = list(lengthMenu = c(5, 10), pageLength = 5)))
  })
  
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
  
  session$onSessionEnded(function() {
    #remove plots
    plots_made = list.files(pattern = analysisID)
    print(plots_made)
    file.remove(plots_made)

  })
  
  
}

shinyApp(ui = ui, server = server)
