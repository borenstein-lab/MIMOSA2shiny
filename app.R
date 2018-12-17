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
library(RColorBrewer)
options(datatable.webSafeMode = TRUE, scipen = 20000, stringsAsFactors = F, shiny.usecairo = F)
theme_set(theme_get() + theme(text = element_text(family = 'Helvetica')))
library(shinyBS)

microbiome_data_upload = function(){
  fluidPage(
    tags$head(tags$style("#fileMet{color: gray }")),

    tags$tr(class = "mainContainer", 
      tags$td(
      h4(get_text("microbiome_header"), width='100%'),
  #,tags$a(tags$img(src = "help.png", border = 0), href = "https://www.github.com/borenstein-lab/", target = "_blank"), width = '100%')
      p(get_text("microbiome_description")),
      column(radioButtons("database", get_text("database_title"), choices = get_text("database_choices"), width = '100%'), width = 8),
      column(tipify(tags$img(src = "help.png", border = 0), title = get_text("microbiome_tooltip"), placement = "right"), width = 4),
      fluidRow(
      column(fileInput("file1", get_text("microbiome_input_title"),
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv"), width = '430px'), width = 8),
        column(tags$a(tags$img(src = "example.png", border = 0), href = "test_seqs.txt", target = "_blank"), width = 4)),
   #   checkboxInput("metagenomeOpt", label=get_text("metagenome_option"), width = '100%'),
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
  id = "microbiome_section" #tags$style(".mainContainer { border: 2px black; }"), 
  )
}

metabolome_data_upload = function(){
  fluidPage(
    h4(get_text("metabolome_header"), tipify(tags$img(src = "help.png", border = 0), title = get_text("metabolome_tooltip"), placement = "right"), width='100%'),
    #tags$a(tags$img(src = "help.png", border = 0), href = "https://www.github.com/borenstein-lab/", target = "_blank"), id = "metabolome", width='100%')
    p(get_text("metabolome_description")),
    radioButtons("metType", label = get_text("met_type_title"), choices = get_text("met_type_choices"), selected = get_text("selected_met_type"), width = '100%'),
	fluidRow(
   		column(fileInput("file2", get_text("metabolome_upload_title"),
              multiple = FALSE,
              accept = c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv"), width = '430px'), width = 8),
    	column(tags$a(tags$img(src = "example.png", border = 0), href = "test_mets.txt", target = "_blank"), width = 4)), id = "metabolome_section")
}

network_settings = function(){
  fluidPage(
    h4(get_text("network_header"),  tipify(tags$img(src = "help.png", border = 0), title = get_text("network_tooltip"), placement = "right"), width='100%'),
       #tags$a(tags$img(src = "help.png", border = 0), href = "https://www.github.com/borenstein-lab/", target = "_blank"), id = "genome", width='100%'),
    p(get_text("network_description")),
    radioButtons("genomeChoices", label = get_text("source_title"), choices = get_text("source_choices"), selected = get_text("source_choices")[2], width = '100%'), # width = 7),
	numericInput("simThreshold", get_text("sim_title"), value = 0.99, min=0.8, max = 1, step = 0.01, width='430px'), #, width = '100%'),
    fluidRow(
   		column(fileInput("netAdd", get_text("net_mod_input_title"), multiple = FALSE, accept = c("text/csv",
                                 "text/comma-separated-values,text/plain",".csv"),  width = '430px'), width = 8), #,  #, width = '100%') #,
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
    radioButtons("contribType", get_text("stat_title"), choices = get_text("stat_choices"))
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
run_pipeline = function(input_data, configTable){
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
      #Implement doing stuff with this later
      metagenome_network = build_metabolic_model(input_data$metagenome, configTable)
      # species2 = metagenome_data[[1]]
      # metagenome_network = metagenome_data[[2]]
    }
    if(configTable[V1=="metType", V2 ==get_text("met_type_choices")[2]]){
      mets = map_to_kegg(mets)
    }
    incProgress(1/10, detail = "Calculating metabolic potential and fitting metabolite concentration model")
    if(configTable[V1=="database", V2==get_text("database_choices")[4]]){
      indiv_cmps = get_cmp_scores_kos(species, network) #Use KO abundances instead of species abundances to get cmps
    } else {
      indiv_cmps = get_species_cmp_scores(species, network)
    }
    mets_melt = melt(mets, id.var = "compound", variable.name = "Sample")
    cmp_mods = fit_cmp_mods(indiv_cmps, mets_melt)
    indiv_cmps = add_residuals(indiv_cmps, cmp_mods[[1]], cmp_mods[[2]])
    incProgress(2/10, detail = "Calculating microbial contributions")
    var_shares = calculate_var_shares(indiv_cmps)
    #shinyjs::logjs(devtools::session_info())
    #Order dataset for plotting
    
    return(list(varShares = var_shares, modelData = cmp_mods[[1]], configs = configTable[!grepl("prefix", V1)], networkData = network))
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
  titlePanel("MIMOSA"),
  mainPanel(fluidPage(p("MIMOSA is a tool for metabolic model-based evaluation of paired microbiome and metabolomics datasets. For more information, see the ", tags$a("manual.", href = "https://cnoecker.github.io/MIMOSA2shiny", target = "_blank"))), id="description"),
  tags$style(type = "text/css", "#title { color: #3CB371; horizontal-align: left; }"  ),
  mainPanel(
    uiOutput("uploadPage")
  ), width='100%')

  #sidebarLayout(
#  sidebarPanel(
#    tags$head(tags$script(HTML('
#                               var fakeClick = function(tabName) {
#                               var dropdownList = document.getElementsByTagName("a");
#                               console.log(dropdownList);
#                               for (var i = 0; i < dropdownList.length; i++) {
#                               var link = dropdownList[i];
#                               if(link.getAttribute("id") == tabName) {
#                               link.click();
#                               };
#                               }
#                               };
#                               '))),
 #   fluidPage(
      #h4("Data input"),
 #     fluidRow(HTML("<a href='#microbiome'>Microbiome data upload</a>")),
      #actionLink("gotomicrobiome", "Microbiome data upload", onclick = "fakeClick('microbiome')")),
      #h4("Settings"),
 #     fluidRow(HTML("<a href='#network'>Metabolic model settings</a>")),
 #     fluidRow(HTML("<a href='#metabolome'>Metabolome data upload</a>")) #,
      #fluidRow(HTML("<a href='#algorithm'>Algorithm settings</a>")) #,
      #fluidRow(HTML("<a href='#outputResults'>Output network settings</a>"))
      #fluidRow(actionLink("gotometabolome", "Metabolomics data upload", onclick = "fakeClick('metabolome')"))
      #fluidRow(h4("this 2nd box should lead me to tab2", onclick = "fakeClick('tab2')"))
 #   ),
    #"Data Input",
    # tabPanel("Microbiome data upload", microbiome_data_upload()),
    # tabPanel("Metabolomics data upload", metabolome_data_upload()),
    # "Settings",
    # tabPanel("Metabolic model settings", network_settings()),
    # tabPanel("Algorithm settings", algorithm_settings()),
    # tabPanel("Output settings", output_settings()),
    # widths = c(4,8)
  #  widths = 3.5
  #  ),



server <- function(input, output, session) {
  #shinyjs::logjs(sessionInfo())
  #shinyjs::logjs(devtools::session_info())
  #shinyjs::logjs(Cairo.capabilities())

  output$uploadPage = renderUI({
    if(is.null(input$goButton) & is.null(input$exampleButton)){
      return(fluidPage(
        #h3("Data Input"),
        microbiome_data_upload(),
        network_settings(),
        metabolome_data_upload(),
        fluidRow(
          column(
            actionButton("goButton", " Run MIMOSA "),
            actionButton("exampleButton", "Show results for example dataset"),
             br(),
             br(),
            tags$style(type='text/css', "#goButton { vertical-align: middle; horizontal-align: middle; font-size: 22px; background-color: #3CBCDB; padding: 5px;}"), 
        	tags$style(type='text/css', "#exampleButton { vertical-align: middle; horizontal-align: middle; font-size: 14px; background-color: #C3D2D5}"), 
            width = 12, align = "center"
          ))))
    } else if(input$goButton==0 & input$exampleButton == 0){
      return(fluidPage(
        #h3("Data Input"),
        microbiome_data_upload(),
        network_settings(),
        metabolome_data_upload(),
        #h3("Settings"),
        #algorithm_settings(),
        #output_settings(),
        fluidRow(
          column(
            actionButton("goButton", " Run MIMOSA "),
             actionButton("exampleButton", "Show results for example dataset"),
             br(),
             br(),
            tags$style(type='text/css', "#goButton { vertical-align: middle; horizontal-align: middle; font-size: 22px; background-color: #3CBCDB; padding: 5px;}"), 
            tags$style(type='text/css', "#exampleButton { vertical-align: middle; horizontal-align: center; font-size: 14px; background-color: #C3D2D5}"), 
            width = 12, align = "center"
          ))))
    } else {
      fluidPage(
      fluidRow(
        column(
      downloadButton("downloadSettings", "Download Record of Configuration Settings"),
      downloadButton("downloadNetworkData", "Download Community Metabolic Network Models"),
      downloadButton("downloadData", "Download Variance Contribution Results"),
      downloadButton("downloadModelData", "Download Model Summaries"),
      tags$style(type='text/css', ".downloadButton { vertical-align: middle; horizontal-align: center; font-size: 22px; background-color: #3CB371}"), width = 12, align = "center")),
      fluidRow(
        #plots
        #plotOutput("residPlot", height = "100px", click = "plot_click"),
        plotOutput("contribPlots", height = "650px", click = "plot_click", hover = "plot_hover") #, click = "plot_click", hover = "plot_hover")
      ),
      fluidRow(
        tableOutput("indivCellInfo")
      ),
      fluidRow(
        plotOutput("indivPlots", height = "600px")
      ))
  
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
  	    initial_inputs = c("closest", "contribType", "database","gapfill", "genomeChoices",  #"geneAdd", 
                       "metType", "simThreshold") #Check that this is all of them
    inputs_provided = initial_inputs[which(sapply(initial_inputs, function(x){ !is.null(input[[x]])}))]
    values_provided = sapply(inputs_provided, function(x){ return(input[[x]])})
	inputs_provided = c(inputs_provided, "kegg_prefix", "data_prefix", "vsearch_path")
	values_provided = c(values_provided, "data/KEGGfiles/", "data/", "bin/vsearch")    #print(file.exists(input$file1$datapath))
    print(inputs_provided)
    print(values_provided)
    return(data.table(V1 = inputs_provided, V2 = values_provided))
  	} else {
  		return(fread("data/exampleData/configs_example_clean.txt", header = T))
  	}
  })
  
  datasetInput <- reactive({
    if(is.null(input$exampleButton)|input$exampleButton==0){
    file_list1 = list(input$file1, input$file2, input$metagenome, input$netAddFile) # input$geneAddFile,
    names(file_list1) = c("file1","file2", "metagenome","netAdd") # "geneAddFile", 
    #logjs(config_table())
    input_data = read_mimosa2_files(file_list = file_list1, configTable = config_table())
    run_pipeline(input_data, config_table())
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
    #var_shares = run_pipeline(input) #Just give it everything and go from there
    
    #}
  })

  observeEvent(input$exampleButton, {
    config_table()
    datasetInput()
    enable("downloadSettings")
    enable("downloadData")
    enable("downloadModelData")
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
    #print(plotData)
    ##met1 = plotData[1,compound]
    #### Make a drop bar to select one metabolite to display plot & data for at a time!!! cool.
    #plotData[,hist(V1)]
    return(plot_summary_contributions(plotData, include_zeros = T, remove_resid_rescale = T))
  }, res = 100)
  
  output$indivCellInfo = renderTable({
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
    #print(input$plot_click)
    #print(levels(x)[input$plot_click[x]]) #etc
    if(!is.null(input$plot_click)){
      print(input$plot_click$x)
      print(plotData[,levels(metID)])
      met_of_interest = plotData[,levels(metID)][round(input$plot_click$x)] 
      #print(met_of_interest)
    } else {
      met_of_interest = plotData[1,metID]      
    }
   	return(plot_contributions(plotData, metabolite = met_of_interest(), include_zeros = F))
  }, res = 100)
}

shinyApp(ui = ui, server = server)
