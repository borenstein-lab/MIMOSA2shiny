.libPaths(c("/data/shiny-server/R/x86_64-redhat-linux-gnu-library/3.2/")) #,"/data/shiny-server/r-packages/", "/data/shiny-server/app_specific_r_packages/"))
library(shiny, lib.loc="/data/shiny-server/R/x86_64-redhat-linux-gnu-library/3.2/")
library(shinyjs, lib.loc="/data/shiny-server/R/x86_64-redhat-linux-gnu-library/3.2/")
library(mimosa, lib.loc="/data/shiny-server/R/x86_64-redhat-linux-gnu-library/3.2/")
library(readr, lib.loc = "/data/shiny-server/R/x86_64-redhat-linux-gnu-library/3.2/")
library(data.table)
library(ggplot2)
library(viridis)

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
    checkboxInput("metagenome", "Also upload metagenome KO abundances"),
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
    h4("Model settings", id = "genome"),
    radioButtons("genomeChoices", "Gene content source", choices = c("Assign KOs with PICRUSt", "Map sequences to AGORA genomes")),
    checkboxInput("geneAdd", "Upload modifications to taxon-gene mapping"),
    disabled(fileInput("geneAddFile", "Upload file of species-gene modifications (see examples and documentation here)",
                       multiple = FALSE,
                       accept = c("text/csv",
                                  "text/comma-separated-values,text/plain",
                                  ".csv"))),
    radioButtons("modelTemplate", "Metabolic model template", choices = c("Generic KEGG metabolic model", "AGORA metabolic models (recommended)")),
    disabled(radioButtons("closest", "", choices = c("Use closest AGORA species", "Use AGORA models for species within a % similarity threshold"))),
    disabled(numericInput("simThreshold", "Similarity threshold:", value = 0.99, min=0.8, max = 1, step = 0.01)),
    p("\n"),
    checkboxInput("netAdd", "Upload modifications to taxon-reaction mapping"), 
    disabled(fileInput("netAddFile", "Upload file of species-reaction modifications (see examples and documentation here)", multiple = FALSE, accept = c("text/csv",
                                  "text/comma-separated-values,text/plain",
                                  ".csv"))),
    checkboxInput("gapfill", "Gap-fill model for each species using x program")
  )
}

algorithm_settings = function(){
  fluidPage(
    h4("Algorithm settings", id = "algorithm"),
    radioButtons("contribType", "Metabolite statistic to analyze:", choices = c("Variance (analytically calculated)", "Differential abundance (Wilcoxon rank-sum, permutation-based)", "Paired-sample differential abundance (paired Wilcoxon rank-sum, permutation based)"))
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
  met_col_name = names(mets)[names(mets) %in% c("compound", "KEGG", "Compound", "metabolite", "Metabolite")]
  if(length(met_col_name) != 1) stop("Ambiguous metabolite ID column name, must be one Compound/KEGG/Metabolite")
  setnames(mets, met_col_name, "compound")
  shared_samps = intersect(names(species), names(mets))
  if(length(shared_samps) < 2) stop("Sample IDs don't match between species and metabolites")
  species = species[,c("OTU", shared_samps), with=F]
  mets = mets[,c("compound", shared_samps), with=F]
  if(input$metagenome){
    #Metagenome data
    #Implement this later
  }
  incProgress(2/10, detail = "Building metabolic model")
  if(input$genomeChoices=="Assign KOs with PICRUSt"){
    if(input$database=="Sequence variants (recommended for AGORA)"){
      seq_list = species[,OTU]
      species_table = get_otus_from_seqvar(seq_list, repSeqDir = "~/Documents/MIMOSA2shiny/data/rep_seqs/", repSeqFile = "gg_13_5.fasta.gz", add_agora_names = F, seqID = simThreshold) #Run vsearch to get gg OTUs
    } else if(input$database != "Greengenes 13_5 or 13_8"){
      stop("Only Greengenes currently implemented")
    }
    contribution_table = generate_contribution_table_using_picrust(species, picrust_norm_file = "data/picrustGenomeData/16S_13_5_precalculated.tab.gz", picrust_ko_table_directory ="data/picrustGenomeData/indivGenomes/", picrust_ko_table_suffix = "_genomic_content.tab")
    contribution_table = contribution_table[contribution != 0]
    if(input$geneAdd){
      req(input$geneAddFile)
      contribution_table = add_genes_to_contribution_table(contribution_table, input$geneAddFile)
    }
  } 
  if(input$modelTemplate=="Generic KEGG metabolic model"){
    network = build_generic_network(contribution_table, kegg_paths = c("data/KEGGfiles/reaction_mapformula.lst", "data/KEGGfiles/reaction_ko.list", "data/KEGGfiles/reaction"))
    save(network, file = "test_network.rda")
  }else{
    network = build_species_networks_w_agora(species, input$database, input$closest, input$simThreshold)
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
  incProgress(2/10, detail = "Calculating metabolic potential")
  indiv_cmps = get_species_cmp_scores(species, network)
  tot_cmps = indiv_cmps[,sum(CMP), by=list(compound, Sample)]
  setnames(tot_cmps, "V1", "value")
  mets_melt = melt(mets, id.var = "compound", variable.name = "Sample")
  cmp_mods = fit_cmp_mods(tot_cmps, mets_melt)
  indiv_cmps = add_residuals(indiv_cmps, cmp_mods[[1]], cmp_mods[[2]])
  incProgress(2/10, detail = "Calculating microbial contributions")
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
      fluidRow(HTML("<a href='#network'>Metabolic network settings</a>")),
      fluidRow(HTML("<a href='#algorithm'>Algorithm settings</a>")),
      fluidRow(HTML("<a href='#outputResults'>Output network settings</a>"))
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
        output_settings(),
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
        output_settings(),
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
    if(input$genomeChoices=="Map sequences to AGORA genomes"){
      enable("closest")
      updateRadioButtons(session, "modelTemplate", selected = "AGORA metabolic models (recommended)")
    } else {
      disable("closest")
      updateRadioButtons(session, "modelTemplate", selected = "Generic KEGG metabolic model")
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
