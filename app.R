library(shiny)
library(shinyjs) #, lib.loc="/data/shiny-server/R/x86_64-redhat-linux-gnu-library/3.2/")
library(mimosa) #, lib.loc="/data/shiny-server/R/x86_64-redhat-linux-gnu-library/3.2/")
library(data.table)
library(ggplot2)
library(viridis)
library(RColorBrewer)
options(datatable.webSafeMode = TRUE, scipen = 20000, stringsAsFactors = F)
#source("scripts/text_constants.R")


microbiome_data_upload = function(){
  fluidPage(
    tags$head(tags$style("#fileMet{color: gray }")),
    h4(get_text("microbiome_header")),
    radioButtons("database", "Microbiome data format:", choices = get_text("database_choices")),
    fileInput("file1", get_text("microbiome_input_title"),
              multiple = FALSE,
              accept = c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")),
    checkboxInput("metagenome", get_text("metagenome_option")),
    #uiOutput("type16S"),
    #uiOutput("typeMetagenome")
    disabled(fileInput("fileMet", get_text("metagenome_input_title"),
                       multiple = FALSE,
                       accept = c("text/csv",
                                  "text/comma-separated-values,text/plain",
                                  ".csv") #,
    ),
    disabled(checkboxInput("metagenome_use", get_text("metagenome_use_option"),
                  ))
    )
  )
  
}
metabolome_data_upload = function(){
  fluidPage(
    h4(get_text("metabolome_header"), id = "metabolome"),
    radioButtons("metType", label = get_text("met_type_title"), choices = get_text("met_type_choices"), selected = get_text("selected_met_type")),
    fileInput("file2", get_text("metabolome_upload_title"),
              multiple = FALSE,
              accept = c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv"))
  )
}

network_settings = function(){
  fluidPage(
    h4(get_text("network_header"), id = "genome"),
    radioButtons("genomeChoices", get_text("source_title"), choices = get_text("source_choices")),
    checkboxInput("geneAdd", get_text("gene_mod_option")),
    disabled(fileInput("geneAddFile", get_text("gene_mod_input_title"),
                       multiple = FALSE,
                       accept = c("text/csv",
                                  "text/comma-separated-values,text/plain",
                                  ".csv"))),
    checkboxInput("netAdd", get_text("net_mod_option")), 
    disabled(fileInput("netAddFile", get_text("net_mod_input_title"), multiple = FALSE, accept = c("text/csv",
                                                                                                                                                         "text/comma-separated-values,text/plain",
                                                                                                                                                         ".csv"))),
    # radioButtons("modelTemplate", "Metabolic model template", choices = c("Generic KEGG metabolic model", "AGORA metabolic models (recommended)")),
    disabled(radioButtons("closest", get_text("closest_title"), choices = get_text("closest_options"))),
    disabled(numericInput("simThreshold", get_text("sim_title"), value = 0.99, min=0.8, max = 1, step = 0.01)),
    p("\n"),
    checkboxInput("gapfill", get_text("gapfill_option"))
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
run_pipeline = function(configTable){
  withProgress(message = "Running MIMOSA!", {
  #process arguments
  print(configTable)
  print(configTable[V1=="file1", V2])
  configTable[,V2:=as.character(V2)]
  print(file.exists(configTable[V1=="file1", V2]))
  print(configTable)
  print(configTable[V1=="file1", V2])
    #cat(input$file1$datapath)
  #cat(input$file2$datapath)
  #run_mimosa2(configTable)
  incProgress(1/10, detail = "Reading data")
  data_inputs = read_mimosa2_files(configTable)
  species = data_inputs[[1]]
  mets = data_inputs[[2]]
  
  if(configTable[V1=="metagenome", V2 != FALSE]){
    metagenome_data = species_network_from_metagenome(configTable[V1=="fileMet", V2])
    species2 = metagenome_data[[1]]
    metagenome_network = metagenome_data[[2]]
  }
  incProgress(2/10, detail = "Building metabolic model")
  network = build_metabolic_model(species, configTable)

  incProgress(1/10, detail = "Calculating metabolic potential")
  indiv_cmps = get_species_cmp_scores(species, network)
  mets_melt = melt(mets, id.var = "compound", variable.name = "Sample")
  cmp_mods = fit_cmp_mods(indiv_cmps, mets_melt)
  indiv_cmps = add_residuals(indiv_cmps, cmp_mods[[1]], cmp_mods[[2]])
  incProgress(1/10, detail = "Calculating microbial contributions")
  var_shares = calculate_var_shares(indiv_cmps)
  return(list(varShares = var_shares, modelData = cmp_mods[[1]]))
  #Send var_shares for download
  #Generate plot of var shares
  #source(other stuff)
  })
}

ui = fluidPage(
  useShinyjs(),
  titlePanel("MIMOSA"),
  #sidebarLayout(
  sidebarPanel(
    tags$head(tags$script(HTML('
                               var fakeClick = function(tabName) {
                               var dropdownList = document.getElementsByTagName("a");
                               console.log(dropdownList);
                               for (var i = 0; i < dropdownList.length; i++) {
                               var link = dropdownList[i];
                               if(link.getAttribute("id") == tabName) {
                               link.click();
                               };
                               }
                               };
                               '))),
    fluidPage(
      h4("Data input"),
      fluidRow(HTML("<a href='#microbiome'>Microbiome data upload</a>")),
      #actionLink("gotomicrobiome", "Microbiome data upload", onclick = "fakeClick('microbiome')")),
      fluidRow(HTML("<a href='#metabolome'>Metabolome data upload</a>")),
      h4("Settings"),
      fluidRow(HTML("<a href='#network'>Metabolic model settings</a>")),
      fluidRow(HTML("<a href='#algorithm'>Algorithm settings</a>")) #,
      #fluidRow(HTML("<a href='#outputResults'>Output network settings</a>"))
      #fluidRow(actionLink("gotometabolome", "Metabolomics data upload", onclick = "fakeClick('metabolome')"))
      #fluidRow(h4("this 2nd box should lead me to tab2", onclick = "fakeClick('tab2')"))
    ),
    #"Data Input",
    # tabPanel("Microbiome data upload", microbiome_data_upload()),
    # tabPanel("Metabolomics data upload", metabolome_data_upload()),
    # "Settings",
    # tabPanel("Metabolic model settings", network_settings()),
    # tabPanel("Algorithm settings", algorithm_settings()),
    # tabPanel("Output settings", output_settings()),
    # widths = c(4,8)
    widths = 3.5
    ),
  mainPanel(
    uiOutput("uploadPage")
  ), fluid = F)



# Define server logic required to draw a histogram ----
server <- function(input, output, session) {
  
  output$uploadPage = renderUI({
    if(is.null(input$goButton)){
      return(fluidPage(
        h3("Data Input"),
        microbiome_data_upload(),
        metabolome_data_upload(),
        h3("Settings"),
        network_settings(),
        algorithm_settings(),
        #output_settings(),
        fluidRow(
          column(
            actionButton("goButton", "Run MIMOSA"),
            tags$style(type='text/css', "#goButton { vertical-align: middle; horizontal-align: center; font-size: 22px; color: #3CB371}"), 
            width = 12, align = "center"
          ))))
    } else if(input$goButton==0){
      return(fluidPage(
        h3("Data Input"),
        microbiome_data_upload(),
        metabolome_data_upload(),
        h3("Settings"),
        network_settings(),
        algorithm_settings(),
        #output_settings(),
        fluidRow(
          column(
            actionButton("goButton", "Run MIMOSA"),
            tags$style(type='text/css', "#goButton { vertical-align: middle; horizontal-align: center; font-size: 22px; color: #3CB371}"), 
            width = 12, align = "center"
          ))))
    } else {
      fluidPage(
      fluidRow(
        column(
      downloadButton("downloadSettings", "Download Record of Configuration Settings"),
      downloadButton("downloadData", "Download Variance Contribution Results"),
      downloadButton("downloadModelData", "Download Model Summaries"),
      tags$style(type='text/css', ".downloadButton { vertical-align: middle; horizontal-align: center; font-size: 22px; color: #3CB371}"), width = 12, align = "center")),
      fluidRow(
        #plots
        plotOutput("contribPlots", click = "plot_click")
      ),
      fluidRow(
        plotOutput("indivPlots")
      ))
  
    }

  })

  observeEvent(input$metagenome, {
    if(input$metagenome==T){
      enable("fileMet")
    } else {
      disable("fileMet")
    }
  })
  observeEvent(input$genomeChoices, {
    if(input$genomeChoices==get_text("source_choices")[2]){
      enable("closest")
    } else {
      disable("closest")
      }
  })
  observeEvent(input$closest, {
    # Change the following line for more examples
    if(input$closest==get_text("closest_options")[2]){
      enable("simThreshold")
    } else {
      disable("simThreshold")
    }
  })
  observeEvent(input$geneAdd, {
    if(input$geneAdd==T){
      enable("geneAddFile")
    } else {
      disable("geneAddFile")
    }
  })
  
  observeEvent(input$netAdd, {
    if(input$netAdd==T){
      enable("netAddFile")
    } else {
      disable("netAddFile")
    }
  })
  
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
    initial_inputs = c("closest", "contribType", "database", "file1", "file2", "fileMet", "gapfill", "geneAdd", "geneAddFile", "genomeChoices", "metagenome", "metagenome_use", 
                       "metType", "netAdd", "netAddFile") #Check that this is all of them
    inputs_provided = initial_inputs[which(sapply(initial_inputs, function(x){ !is.null(input[[x]])}))]
    values_provided = sapply(inputs_provided, function(x){ return(ifelse("datapath" %in% names(x), input[[x]]$datapath, input[[x]]))})
    print(file.exists(input$file1$datapath))
    print(inputs_provided)
    return(data.table(V1 = inputs_provided, V2 = values_provided))
  })
  
  datasetInput <- reactive({
    config_table()
    run_pipeline(config_table())
  })
  
  
  observeEvent(input$goButton, {
    req(input$file1)
    req(input$file2)
    #if(error_checks){
    print(names(input))
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
    plot_summary_contributions(plotData, include_zeros = T)
  })
  output$indivPlots = renderPlot({
    plotData = datasetInput()$varShares
    print(input$plot_click)
    #levels(x)[input$plot_click[x]] #etc
    if(!is.null(input$plot_click)){
      met_of_interest = plotData[,unique(metID)][round(input$plot_click$x)]
      print(met_of_interest)
      plot_contributions(plotData, metabolite = met_of_interest, include_zeros = F)
    }
  })
}

shinyApp(ui = ui, server = server)
