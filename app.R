.libPaths(c("/data/shiny-server/r-packages/", "/data/shiny-server/R/x86_64-redhat-linux-gnu-library/3.2/")) #, "/data/shiny-server/app_specific_r_packages/"))
library(Rcpp, lib.loc = "/data/shiny-server/r-packages/")
library(Cairo, lib.loc = "/data/shiny-server/r-packages/")
library(shiny)
library(shinyjs)
library(mimosa, lib.loc ="/data/shiny-server/r-packages/")
library(data.table) #, lib.loc ="/data/shiny-server/R/x86_64-redhat-linux-gnu-library/3.2/")
library(readr)
library(ggplot2, lib.loc = "/data/shiny-server/r-packages")
library(viridis, lib.loc = "/data/shiny-server/r-packages")
library(RColorBrewer)
options(datatable.webSafeMode = TRUE, scipen = 20000, stringsAsFactors = F, shiny.usecairo = F)
theme_set(theme_get() + theme(text = element_text(family = 'Helvetica')))


microbiome_data_upload = function(){
  fluidPage(
    tags$head(tags$style("#fileMet{color: gray }")),

    tags$tr(class = "mainContainer", 
      tags$td(
      h4(get_text("microbiome_header") ,tags$a(tags$img(src = "help.png", border = 0), href = "https://www.github.com/borenstein-lab/", target = "_blank"), width = '100%'),
      radioButtons("database", get_text("database_title"), choices = get_text("database_choices"), width = '100%'),
      fluidRow(
      column(fileInput("file1", get_text("microbiome_input_title"),
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv"), width = '400px'), width = 8),
        column(tags$a(tags$img(src = "example.png", border = 0), href = "test_seqs.txt", target = "_blank"), width = 4)),
   #   checkboxInput("metagenomeOpt", label=get_text("metagenome_option"), width = '100%'),
   	fluidRow(
   		column(      fileInput("metagenome", get_text("metagenome_input_title"),
                         multiple = FALSE,
                         accept = c("text/csv",
                                    "text/comma-separated-values,text/plain",
                                    ".csv"), width = '400px'), width = 8),
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
    h4(get_text("metabolome_header"), tags$a(tags$img(src = "help.png", border = 0), href = "https://www.github.com/borenstein-lab/", target = "_blank"), id = "metabolome", width='100%'),
    radioButtons("metType", label = get_text("met_type_title"), choices = get_text("met_type_choices"), selected = get_text("selected_met_type"), width = '100%'),
	fluidRow(
   		column(fileInput("file2", get_text("metabolome_upload_title"),
              multiple = FALSE,
              accept = c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv"), width = '400px'), width = 8),
    	column(tags$a(tags$img(src = "example.png", border = 0), href = "test_mets.txt", target = "_blank"), width = 4)), id = "metabolome_section")
}

network_settings = function(){
  fluidPage(
    h4(get_text("network_header"), tags$a(tags$img(src = "help.png", border = 0), href = "https://www.github.com/borenstein-lab/", target = "_blank"), id = "genome", width='100%'),
    radioButtons("genomeChoices", label = get_text("source_title"), choices = get_text("source_choices"), selected = get_text("source_choices")[2], width = '100%'), # width = 7),
	numericInput("simThreshold", get_text("sim_title"), value = 0.99, min=0.8, max = 1, step = 0.01), #, width = '100%'),
    fluidRow(
   		column(fileInput("netAdd", get_text("net_mod_input_title"), multiple = FALSE, accept = c("text/csv",
                                 "text/comma-separated-values,text/plain",".csv"),  width = '400px'), width = 8), #,  #, width = '100%') #,
                                     column(tags$a(tags$img(src = "example.png", border = 0), href = "test_netAdd.txt", target = "_blank"), width = 4)
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
    #cat(input$file1$datapath)
    #cat(input$file2$datapath)
    #run_mimosa2(configTable)

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
    if(configTable[V1=="metType", V2 !=get_text("met_type_choices")[1]]){
      #mets = map_to_kegg(mets)
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
    shinyjs::logjs(devtools::session_info())
    #Order dataset for plotting
    
    return(list(varShares = var_shares, modelData = cmp_mods[[1]]))
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
  shinyjs::logjs(sessionInfo())
  shinyjs::logjs(devtools::session_info())
  shinyjs::logjs(Cairo.capabilities())

  output$uploadPage = renderUI({
    if(is.null(input$goButton)){
      return(fluidPage(
        #h3("Data Input"),
        microbiome_data_upload(),
#        tags$td(
#        tags$img(border="0", src="pixel.gif", width='100%', height="10", hspace="0", vspace="0"),
#        bgcolor="#EBEAE1"),
        network_settings(),
#        tags$img(border="0", src="pixel.gif", width='100%', height="10", hspace="0", vspace="0"),
        metabolome_data_upload(),
        #h3("Settings"),
        #algorithm_settings(),
        #output_settings(),
        fluidRow(
          column(
            actionButton("goButton", "Run MIMOSA"),
            tags$style(type='text/css', "#goButton { vertical-align: middle; horizontal-align: middle; font-size: 22px; background-color: #3CB371}"), 
            width = 12, align = "center"
          ))))
    } else if(input$goButton==0){
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
            actionButton("goButton", "Run MIMOSA"),
            tags$style(type='text/css', "#goButton { vertical-align: middle; horizontal-align: center; font-size: 22px; background-color: #3CB371}"), 
            width = 12, align = "center"
          ))))
    } else {
      fluidPage(
      fluidRow(
        column(
      downloadButton("downloadSettings", "Download Record of Configuration Settings"),
      downloadButton("downloadData", "Download Variance Contribution Results"),
      downloadButton("downloadModelData", "Download Model Summaries"),
      tags$style(type='text/css', ".downloadButton { vertical-align: middle; horizontal-align: center; font-size: 22px; background-color: #3CB371}"), width = 12, align = "center")),
      fluidRow(
        #plots
        plotOutput("contribPlots") #, click = "plot_click", hover = "plot_hover")
      ),
      fluidRow(
        tableOutput("indivCellInfo")
      ),
      fluidRow(
        plotOutput("indivPlots")
      ))
  
    }

  })
  observeEvent(input$database, {
    if(input$database==get_text("database_choices")[1]){
      updateRadioButtons(session, "genomeChoices", selected = get_text("source_choices")[2])
    } else {
      updateRadioButtons(session, "genomeChoices", selected = get_text("source_choices")[1])
    }
    if(input$database==get_text("database_choices")[4]){
      #updateCheckboxInput(session, "metagenomeOpt", value = T)
      disable("file1")
    } else {
      #updateCheckboxInput(session, "metagenomeOpt", value = F)
      enable("file1")
    }
  })
#  observeEvent(input$metagenomeOpt, {
#    if(input$metagenomeOpt==T){
#      enable("metagenome")
#    } else {
#      disable("metagenome")
#    }
#  })
  observeEvent(input$genomeChoices, {
    if(input$genomeChoices==get_text("source_choices")[2] & input$database==get_text("database_choices")[1]){ #Only enable mapping threshold if providing sequences
      enable("simThreshold")
    } else {
      disable("simThreshold")
      }
  })
  # observeEvent(input$closest, {
  #   # Change the following line for more examples
  #   if(input$closest==get_text("closest_options")[2]){
  #     enable("simThreshold")
  #   } else {
  #     disable("simThreshold")
  #   }
  # })
  # observeEvent(input$geneAdd, {
  #   if(input$geneAdd==T){
  #     enable("geneAddFile")
  #   } else {
  #     disable("geneAddFile")
  #   }
  # })
  
#  observeEvent(input$netAdd, {
#    if(input$netAdd==T){
#      enable("netAddFile")
#    } else {
#      disable("netAddFile")
#    }
#  })
  
  # output$contents = renderTable({
  #   
  #   # input$file1 will be NULL initially. After the user selects
  #   # and uploads a file, head of that data file by default,
  #   # or all rows if selected, will be shown.
  #   input$goButton
  #   req(input$file1)
  #   
  #   # when reading semicolon separated files,
  #   # having a comma separator causes `read.csv` to error
  #   tryCatch(
  #     {
  #       df <- read.csv(input$file1$datapath,
  #                      header = input$header,
  #                      sep = input$sep,
  #                      quote = input$quote)
  #     },
  #     error = function(e) {
  #       # return a safeError if a parsing error occurs
  #       stop(safeError(e))
  #     }
  #   )
  #   
  #   if(input$disp == "head") {
  #     return(head(df))
  #   }
  #   else {
  #     return(df)
  #   }
  # })
  # 
  config_table = reactive({
    initial_inputs = c("closest", "contribType", "database","gapfill", "genomeChoices",  #"geneAdd", 
                       "metType", "simThreshold") #Check that this is all of them
    inputs_provided = initial_inputs[which(sapply(initial_inputs, function(x){ !is.null(input[[x]])}))]
    values_provided = sapply(inputs_provided, function(x){ return(input[[x]])})
    inputs_provided = c(inputs_provided, "kegg_prefix")
    values_provided = c(values_provided, "data/KEGGfiles/")
    #print(file.exists(input$file1$datapath))
    print(inputs_provided)
    print(values_provided)
    return(data.table(V1 = inputs_provided, V2 = values_provided))
  })
  
  datasetInput <- reactive({
    file_list1 = list(input$file1, input$file2, input$metagenome, input$netAddFile) # input$geneAddFile,
    names(file_list1) = c("file1","file2", "metagenome","netAdd") # "geneAddFile", 
    print(file_list1$file1)
    print(file_list1$metagenome)
    print(config_table())
    logjs(config_table())
    input_data = read_mimosa2_files(file_list = file_list1, configTable = config_table())
    print(input_data)
    run_pipeline(input_data, config_table())
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
  output$contribPlots = renderPlot({
    plotData = datasetInput()$varShares
    print(plotData)
    ##met1 = plotData[1,compound]
    #### Make a drop bar to select one metabolite to display plot & data for at a time!!! cool.
    #plotData[,hist(V1)]
    return(plot_summary_contributions(plotData, include_zeros = T))
  }, res = 100)
  output$indivCellInfo = renderTable({
    plotData = datasetInput()$varShares
    if(!is.null(input$plot_hover)){
      met_of_interest = plotData[,levels(metID)][round(input$plot_click$x)] #Plot in alphabetical/numeric order
      print(met_of_interest)
      spec_of_interest = plotData[,levels(Species)][round(input$plot_click$y)]
      print(spec_of_interest)
      return(plotData[Species==spec_of_interest & metID==met_of_interest])
    } else {
      return(plotData[1])
    }
  })
  output$indivPlots = renderPlot({
    plotData = datasetInput()$varShares
    print(input$plot_click)
    #levels(x)[input$plot_click[x]] #etc
    if(!is.null(input$plot_click)){
      met_of_interest = plotData[,levels(metID)][round(input$plot_click$x)] #Plot in alphabetical/numeric order
      print(met_of_interest)
    } else {
      met_of_interest = plotData[1,metID]      
    }
   	return(plot_contributions(plotData, metabolite = met_of_interest, include_zeros = F))
  }, res = 100)
}

shinyApp(ui = ui, server = server)
