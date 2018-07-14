.libPaths(c("/data/shiny-server/R/x86_64-redhat-linux-gnu-library/3.2/")) #,"/data/shiny-server/r-packages/", "/data/shiny-server/app_specific_r_packages/"))
library(shiny, lib.loc="/data/shiny-server/R/x86_64-redhat-linux-gnu-library/3.2/")
library(shinyjs, lib.loc="/data/shiny-server/R/x86_64-redhat-linux-gnu-library/3.2/")
library(mimosa, lib.loc="/data/shiny-server/R/x86_64-redhat-linux-gnu-library/3.2/")
library(readr, lib.loc = "/data/shiny-server/R/x86_64-redhat-linux-gnu-library/3.2/")
library(data.table)

microbiome_data_upload = function(){
  fluidPage(
    tags$head(tags$style("#fileMet{color: gray }")),
    h4("Microbiome data upload: provide a tab-delimited file with a row for each species and a column for each sample."),
    radioButtons("database", "16S format:", choices = c("Sequence variants (recommended for AGORA)", "Greengenes 13_5 or 13_8", "SILVA")),
    fileInput("file1", "Upload 16S rRNA abundance file (example format linked here)",
              multiple = FALSE,
              accept = c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")),
    checkboxInput("specType", "Also upload metagenome KO abundances"),
    #uiOutput("type16S"),
    #uiOutput("typeMetagenome")
    disabled(fileInput("fileMet", "Upload file of metagenomic KO abundances (example format linked here)",
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
    h4("Metabolite data upload: provide a tab-delimited file with a row for each metabolite and a column for each sample.", id = "metabolome"),
    radioButtons("metType", label = "Compound IDs:", choices = c("KEGG Compound IDs", "MetaCyc Compound IDs", "Metabolite names (search for matching ID)"), selected = "KEGG Compound IDs"),
    fileInput("file2", "Upload metabolite file (example format linked here)",
              multiple = FALSE,
              accept = c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv"))
  )
}

network_settings = function(){
  fluidPage(
    h4("Metabolic network", id = "network"),
    radioButtons("AGORA", "Metabolic model template", choices = c("Generic PICRUST/KEGG metabolic model", "Map species to AGORA genomes and use AGORA metabolic models (recommended)")),
    disabled(radioButtons("closest", "", choices = c("Use closest AGORA species", "Use AGORA models for species within a % similarity threshold"))),
    disabled(numericInput("simThreshold", "Similarity threshold:", value = 0.99, min=0.8, max = 1, step = 0.01)),
    checkboxGroupInput("netAdd", "", choices = c("Gap-fill metabolic network of each species","Add manual reactions")),
    disabled(fileInput("rxnfile", "Upload file of reactions to add to the model (can be species-specific or generic, example format linked here)",
                       multiple = FALSE,
                       accept = c("text/csv",
                                  "text/comma-separated-values,text/plain",
                                  ".csv")))
  )
}

algorithm_settings = function(){
  fluidPage(
    h4("Algorithm", id = "algorithm"),
    radioButtons("contribType", "Metabolite statistic to analyze:", choices = c("Variance (analytically calculated)", "Differential abundance (Wilcoxon rank-sum, permutation-based)", "Paired-sample differential abundance (paired Wilcoxon rank-sum, permutation based)"))
  )
}

output_settings = function(){
  fluidPage(
    h4("Output", id = "output"),
    radioButtons("returnType", "Select desired output:", choices = c("Summarized contribution results", "Full verbose data output")),
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
  met_col_name = names(mets)[names(mets) %in% c("compound", "KEGG", "Compound", "metabolite", "Metabolite")]
  if(length(met_col_name) != 1) stop("Ambiguous metabolite ID column name, must be one Compound/KEGG/Metabolite")
  setnames(mets, met_col_name, "compound")
  shared_samps = intersect(names(species), names(mets))
  if(length(shared_samps) < 2) stop("Sample IDs don't match between species and metabolites")
  species = species[,c("OTU", shared_samps), with=F]
  mets = mets[,c("compound", shared_samps), with=F]
  if(input$specType){
    #Metagenome data
    #Implement this later
  }
  incProgress(2/10, detail = "Building metabolic model")
  if(input$AGORA=="Generic PICRUST/KEGG metabolic model"){
    network = build_generic_network(species, input$database, picrust_paths = c("data/picrustGenomeData/16S_13_5_precalculated.tab.gz", "data/picrustGenomeData/indivGenomes/",
                                                                               "_genomic_content.tab"), kegg_paths = c("data/KEGGfiles/reaction_mapformula.lst", "data/KEGGfiles/reaction_ko.list", "data/KEGGfiles/reaction"))
    save(network, file = "test_network.rda")
  }else{
    network = build_species_networks_w_agora(species, input$database, input$closest, input$simThreshold)
  }
  if("Add manual reactions" %in% input$netAdd){
    network = add_rxns_to_network(network, input$rxnfile)
    #This will need to map between metabolite IDs possibly
  }
  if("Gap-fill metabolic network of each species" %in% input$netAdd){
    #Do stuff
  }
  if(input$metType!="KEGG Compound IDs"){
    #mets = map_to_kegg(mets)
    #Implement this later
  }
  incProgress(2/10, detail = "Calculating metabolic potential")
  indiv_cmps = get_species_cmp_scores(species, network)
  tot_cmps = indiv_cmps[,sum(CMP), by=list(compound, Sample)]
  setnames(tot_cmps, "V1", "value")
  mets_melt = melt(mets, id.var = "compound", variable.name = "Sample")
  cmp_mods = fit_cmp_mods(tot_cmps, mets_melt)
  indiv_cmps = add_residuals(indiv_cmps, cmp_mods[[1]], cmp_mods[[2]])
  incProgress(2/10, detail = "Calculating microbial contributions")
  var_shares = calculate_var_shares(indiv_cmps)
  return(var_shares)
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
      fluidRow(HTML("<a href='#network'>Metabolic network settings</a>")),
      fluidRow(HTML("<a href='#algorithm'>Algorithm settings</a>")),
      fluidRow(HTML("<a href='#output'>Output network settings</a>"))
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
    h3("Data Input"),
    microbiome_data_upload(),
    metabolome_data_upload(),
    h3("Settings"),
    network_settings(),
    algorithm_settings(),
    output_settings(),
    fluidRow(
      column(
        actionButton("goButton", "Run MIMOSA"),
        disabled(downloadButton("downloadData", "Download Results")),
        tags$style(type='text/css', "#goButton { vertical-align: middle; horizontal-align: center; font-size: 22px; color: #3CB371}"), 
        tags$style(type='text/css', "#downloadData { vertical-align: middle; horizontal-align: center; font-size: 22px; color: #3CB371}"),
        width = 12, align = "center"
      ))
    #)
  ), fluid = F)



# Define server logic required to draw a histogram ----
server <- function(input, output) {
  print(sessionInfo())
  # output$type16S = renderUI({
  #   if (is.null(input$specType)|!"16S rRNA" %in% input$specType) return(
  #     fluidPage(
  #       radioButtons("database", "16S format:", choices = c("Sequence variants (will map to Greengenes)", "Greengenes 13_5 or 13_8 97% OTUs", "SILVA OTUs")),
  #       fileInput("file1", "Upload 16S rRNA abundance file (example format linked here)",
  #                 multiple = FALSE,
  #                 accept = c("text/csv",
  #                            "text/comma-separated-values,text/plain",
  #                            ".csv")),
  #       tags$head(tags$style("#text1{color: gray;
  #                                }"
  #       )
  #     )
  #     )
  #   ) else {
  #     fluidPage(
  #     radioButtons("database", "16S format:", choices = c("Sequence variants (will map to Greengenes)", "Greengenes 13_5 or 13_8 97% OTUs", "SILVA OTUs")),
  #     fileInput("file1", "Upload 16S rRNA abundance file (example format linked here)",
  #               multiple = FALSE,
  #               accept = c("text/csv",
  #                          "text/comma-separated-values,text/plain",
  #                          ".csv"))
  #     )
  #   }
  # })
  observeEvent(input$specType, {
    if(input$specType==T){
      enable("fileMet")
    } else {
      disable("fileMet")
    }
  })
  observeEvent(input$AGORA, {
    if(input$AGORA=="Map species to AGORA genomes and use AGORA metabolic models (recommended)"){
      enable("closest")
    } else {
      disable("closest")
    }
  })
  observeEvent(input$closest, {
    # Change the following line for more examples
    if(input$closest=="Use AGORA models for species within a % similarity threshold"){
      enable("simThreshold")
    } else {
      disable("simThreshold")
    }
  })
  observeEvent(input$netAdd, {
    if("Add manual reactions" %in% input$netAdd){
      enable("rxnFile")
    } else {
      disable("rxnFile")
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
  	print(.libPaths())
    req(input$file1)
    req(input$file2)
    run_pipeline(input)
  })
  
  observeEvent(input$goButton, {
    #shiny::req(input$file1)
    #shiny::req(input$file2)
    #if(error_checks){
    print(names(input))
    
    datasetInput()
    enable("downloadData")
    #var_shares = run_pipeline(input) #Just give it everything and go from there
    
    #}
  })
  
  
  output$downloadData <- downloadHandler(
    filename = function() {
      "results.txt"
      #paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      write.table(datasetInput(), file, row.names = F, sep = "\t", quote=F)
    }
  )
  
}

shinyApp(ui = ui, server = server)
