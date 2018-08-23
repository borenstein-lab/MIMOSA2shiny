library(shiny)
library(shinyjs) #, lib.loc="/data/shiny-server/R/x86_64-redhat-linux-gnu-library/3.2/")
library(mimosa) #, lib.loc="/data/shiny-server/R/x86_64-redhat-linux-gnu-library/3.2/")
library(data.table)
library(ggplot2)
library(viridis)
#source("scripts/text_constants.R")

microbiome_data_upload = function(){
  fluidPage(
    tags$head(tags$style("#fileMet{color: gray }")),
    h4(get_text("microbiome_header")),
    radioButtons("database", "16S format:", choices = get_text("database_choices")),
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
    )
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
run_pipeline = function(input){
  withProgress(message = "Running MIMOSA!", {
  #process arguments
  cat(input$file1$datapath)
  cat(input$file2$datapath)
  incProgress(1/10, detail = "Reading data")
  species = fread(input$file1$datapath)
  species = spec_table_fix(species)
  mets = fread(input$file2$datapath)
  met_nonzero_filt = 5
  mets = met_table_fix(mets, met_nonzero_filt)
  #Filter species using default abundance values
  species = filter_species_abunds(species, filter_type = "fracNonzero")
  species = filter_species_abunds(species, filter_type = "mean")
  shared_samps = intersect(names(species), names(mets))
  if(length(shared_samps) < 2) stop("Sample IDs don't match between species and metabolites")
  species = species[,c("OTU", shared_samps), with=F]
  mets = mets[,c("compound", shared_samps), with=F]
  if(input$metagenome){
    #Metagenome data
    #Implement this later
  }
  incProgress(2/10, detail = "Building metabolic model")
  if(input$genomeChoices==get_text("source_choices")[1]){
    if(input$database==get_text("database_choices")[1]){
      seq_list = species[,OTU]
      incProgress(1/10, detail = "Assigning sequences to greengenes OTUs")
      species_table = get_otus_from_seqvar(seq_list, repSeqDir = "~/Documents/MIMOSA2shiny/data/rep_seqs/", repSeqFile = "gg_13_5.fasta.gz", add_agora_names = F, seqID = 0.97) #Run vsearch to get gg OTUs
    } else if(input$database != get_text("database_choices")[2]){
      stop("Only Greengenes currently implemented")
    }
    if(input$geneAdd){
      req(input$geneAddFile)
      contribution_table = generate_contribution_table_using_picrust(species, picrust_norm_file = "data/picrustGenomeData/16S_13_5_precalculated.tab.gz", picrust_ko_table_directory ="data/picrustGenomeData/indivGenomes/", picrust_ko_table_suffix = "_genomic_content.tab")
      contribution_table = add_genes_to_contribution_table(contribution_table, input$geneAddFile)
      network = build_generic_network(contribution_table, kegg_paths = c("data/KEGGfiles/reaction_mapformula.lst", "data/KEGGfiles/reaction_ko.list", "data/KEGGfiles/reaction"))
    } else { #Just load in preprocessed
      network = get_kegg_network(species, net_path = "data/picrustGenomeData/indivModels/")
    }
  } else if(input$genomeChoices==get_text("source_choices")[2]){
    network_results = build_species_networks_w_agora(species, input$database, input$closest, input$simThreshold)
    species = network_results[[1]]
    network = network_results[[2]]
    species = species[OTU %in% network[,OTU]]
  }
  if(input$netAdd){
    network = add_rxns_to_network(network, input$netAddFile)
    #This will need to map between metabolite IDs possibly
  }
  if(input$gapfill){
    #Do stuff
  }
  if(input$metType!="KEGG Compound IDs"){
    #mets = map_to_kegg(mets)
    #Implement this later
  }
  incProgress(1/10, detail = "Calculating metabolic potential")
  indiv_cmps = get_species_cmp_scores(species, network)
  if(input$genomeChoices==get_text("source_choices")[2]){ #Switch to KEGG IDs at this point
    indiv_cmps[,KEGG:=agora_kegg_mets(compound)]
    indiv_cmps = indiv_cmps[,sum(CMP), by=list(Species, KEGG, Sample)] #Check that this makes sense
    #separate internal/external?
    setnames(indiv_cmps, c("KEGG", "V1"), c("compound", "CMP"))
  }
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
      downloadButton("downloadSettings", "Download Configuration/Settings File"),
      downloadButton("downloadData", "Download Variance Contribution Results"),
      downloadButton("downloadModelData", "Download Model Summaries"),
      tags$style(type='text/css', ".downloadButton { vertical-align: middle; horizontal-align: center; font-size: 22px; color: #3CB371}"), width = 12, align = "center")),
      fluidRow(
        #plots
        plotOutput("contribPlots")
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
  datasetInput <- reactive({
    run_pipeline(input)
  })
  
  observeEvent(input$goButton, {
    req(input$file1)
    req(input$file2)
    #if(error_checks){
    print(names(input))
    datasetInput()
    enable("downloadData")
    enable("downloadModelData")
    #var_shares = run_pipeline(input) #Just give it everything and go from there
    
    #}
  })
  
  
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
    met1 = plotData[1,compound]
    #### Make a drop bar to select one metabolite to display plot & data for at a time!!! cool.
    plot_contributions(plotData, metabolite = met1)
  })
}

shinyApp(ui = ui, server = server)
