---
title: Using the mimosa2 R package
layout: default
active: package
---
# Using the mimosa2 R package

The easiest way to run a MIMOSA2 analysis is via the web application. However, you can also install the mimosa2 R package
 to run your own custom analyses and integrate MIMOSA2 into your analysis pipelines.

## Installation

mimosa2 can be easily installed from GitHub using the `devtools` package:

```R
devtools::install_github("borenstein-lab/mimosa2", dependencies = T)
library(mimosa)
``` 

If you see some warnings when loading the package, this is normal. 

If you want to analyze ASV data, you will also need to have the program *vsearch* installed. Visit the [vsearch website](https://github.com/torognes/vsearch) to download and install.

If you want to test your installation (including optionally access to vsearch), you can run: 

```R
test_m2_analysis(test_vsearch = T)
```

This will only test the package installation, not the setup of your reference data (see below).

## Generating or downloading preprocessed reference data

Before running a MIMOSA2 analysis, the reference data you would like to use for the analysis needs to be downloaded and/or generated. MIMOSA2 relies on two separate types of reference data: 

1) Reference data to link ASVs to reference taxa (not necessary if you have metagenomic KO annotation data)

2) Sets of genes and reactions linked to reference taxa

The figure below illustrates all the possible combinations of input and reference data formats. You can set up the reference databases to run all of them or just a subset for a specific analysis.

![reference chart](FigureS1_modelBuilding.png "Reference Flow Chart")

Several of these are available from the [Downloads](download.html) page. These can also be regenerated using scripts provided in the [MIMOSA2 GitHub repository](https://github.com/cnoecker/MIMOSA2shiny/). 

If you would like to run a workflow that uses freely available data, the `download_reference_data` function will obtain the necessary data and format it as expected by the main MIMOSA2 analysis. This function takes two arguments, which correspond to the *file1_type* and *ref_choices* options in the configuration table for a MIMOSA2 analysis (below).

```R
download_reference_data(seq_db = "Sequence variants (ASVs)", target_db = "AGORA genomes and models")

download_reference_data(seq_db = "Greengenes 13_5 or 13_8 OTUs", target_db = "RefSeq/EMBL_GEMs genomes and models")
```
You can use the `save_to` argument to customize where these files are saved, but if you change this you will need to modify the `data_prefix` argument when running your MIMOSA2 analysis (see below). The result of this function should be to produce a directory called "data" containing a sub-directory called either "AGORA" or "embl_gems", which contains mapping data as well as a further subdirectory called "RxnNetworks" containing metabolic reference data for each taxon.

If you would like to run an analysis using KEGG, you need to have a KEGG license and to download 3 files from the KEGG FTP server: annotated pathway reactions (filename *reaction_mapformula.lst*), reaction annotations (filename *reaction*), and reaction-KO links (filename *ko_reaction.list*). Then you can provide those files as input to the `generate_preprocessed_networks` function to set up the reference database for MIMOSA2.

To test that your reference databases are formatted and set up as expected by the package for a particular analysis, you can use `check_ref_data`, for example:

```R
check_ref_data(seq_db = "Sequence variants (ASVs)", target_db = "PICRUSt KO genomes and KEGG metabolic model")
```

## Run a full MIMOSA2 analysis

Once you have downloaded and set up the relevant reference databases, you can run a full MIMOSA2 analysis simply by providing a "configuration table" containing all of the relevant settings for the analysis to the `run_mimosa2` function.
The table below lists the various fields that you can provide in your configuration table. Required fields are in bold. The fields or rows of the table can be specified in any order.

| Field | Description | Possible values |
|------|----------|---------|
|**file1** | Microbiome file path | Valid file path|
|**file2** | Metabolomics file path | Valid file path|
|**file1_type** | Taxonomic abundance file type| One of: "Sequence variants (ASVs)", "Greengenes 13_5 or 13_8 OTUs", "SILVA 132 OTUs", "Metagenome: Total KO abundances", "Metagenome: Taxon-stratified KO abundances (HUMAnN2 or PICRUSt/PICRUSt2)" |
|**ref_choices** | Ref model option | One of: "PICRUSt KO genomes and KEGG metabolic model", "AGORA genomes and models", "RefSeq/EMBL_GEMs genomes and models" |
|**data_prefix** | File path to reference databases (see below for required files)| Valid file path|
|simThreshold | If 16S rRNA ASVs are provided, threshold for mapping them to a reference database | Value from 0 to 1 (default 0.99)|
|netAdd | File path to network modifications file | Valid file path|
|metType | Whether metabolite data is provided as KEGG compound IDs or metabolite names (assumes KEGG if not provided) | One of: "KEGG Compound IDs", "Metabolite names (search for matching ID)" (default "KEGG Compound IDs") |
|signifThreshold | Taxonomic contributors to metabolites will only be evaluated for metabolites with a model fit p-value below this threshold | Value from 0 to 1 (default 0.2)|
|compare_only | Skip the taxonomic contribution analysis, only build the model and compare CMP scores and metabolites | T or F (default F)|
|logTransform | Whether a log transform should be applied to metabolite data| T or F (default F)|
|rankBased | Whether to use rank-based regression for comparing CMP scores and metabolites| T or F (default T)|
|vsearch_path | File path to vsearch executable | Valid path (when not provided, MIMOSA2 assumes vsearch is in the executable path)|

Some example configuration tables are linked below:

- [An ASV-based analysis using EMBL_GEMS models and rank-based regression](config_example1.txt) 
- [A stratified metagenome-based analysis using KEGG and OLS regression, with a custom network modification](config_example2.txt)
- [A GreenGenes OTU table mapped to AGORA models](config_example3.txt)

You can also download the configuration table used to run any analysis on the MIMOSA2 web server, which allows anyone to later reproduce the same analysis in an R session.


Once you have downloaded the necessary reference data, installed MIMOSA2 and vsearch, and created a configuration table, it is easy to run a full MIMOSA2 analysis. Save it as a text document, for example "configuration_table1.txt", and run the following in an R session or script: 

```R
mimosa_results = run_mimosa2("configuration_table1.txt")
```

The run_mimosa2 function returns a list of data tables that is identical to the set of results provided by the web application. More details about the results are provided on the [Results](results.html) page.

If you want to generate plots of metabolic potential and taxonomic contributors for each metabolite, similar to the web app, use the `make_plots` and `save_plots` arguments for run_mimosa2. 

```R
mimosa_results_make_plots = run_mimosa2("configuration_table1.txt", make_plots = T, save_plots = T)
```

In this case lists of plots will also be returned. If `save_plots` is true, the function will save all plots in a folder in the current working directory.

If you want to run MIMOSA2 using microbiome and metabolite data frames in an existing R session, you can use the `species` and `mets` arguments (below). Note that this means that the analysis will skip several initial data filters and quality checks that would otherwise be performed when importing the data.

```R
mimosa_results = run_mimosa2("configuration_table1.txt", species = my_otu_table, mets = my_met_table)
```


<!---

## Run individual components of a MIMOSA2 analysis

## Other utility functions

- `format_humann2_contributions`

- `map_to_kegg`

- `plot_summary_contributions`: Make a heatmap

<h4 id="processRefs">Processing Reference Data for Compatibility with MIMOSA2</h4>

You can also generate your own version of the AGORA or RefSeq databases using the package function `download_ribosomal_ref_seqs`, which uses the [biomartR](https://ropensci.github.io/biomartr/) package to download the relevant list of accessions from NCBI. 

#### Download and reformat metabolic reconstructions
If you do not wish to use the precomputed files provided above (i.e. to use a new version of a database), you can generate a reference database set up for MIMOSA2 yourself. 

- [AGORA](www.vmh.life)
- [embl_gems](www.github.com/cdanielmachado/embl_gems/)
- PICRUSt_KEGG (PICRUSt1 pre-calculated files, plus KEGG FTP downloads): KEGG network with GreenGenes OTUs and PICRUSt 1
- KEGG (KEGG FTP downloads): KEGG network without taxonomic information, to use with KEGG-annotated metagenomic data

Note that downloading the necessary KEGG FTP files requires a license. 

If you use AGORA 1.0.2 and the current (2019) version of embl_gems, information on each model is included in the package data. Otherwise, if you use a different database version,
you may need to create a new model info file, listing each model, its 16S copy number (if known), and an ID that is shared between the model itself and its linked sequence data. You can see an example for AGORA [here](link, put data up). 

Next, use the `generate_preprocessed_networks` function to reformat the model files to be compatible with MIMOSA2. For example, to format the *embl_gems* database models:

```R
generate_preprocessed_networks("embl_gems", dat_path = file_path_to_raw_models, out_path = file_path_for_output)
```

(In the function call above, you would replace `file_path_to_raw_models` and `file_path_for_output` with your corresponding file paths.) 

To use GreenGenes and KEGG, you can download the [GreenGenes representative OTU sequence files](http://greengenes.secondgenome.com/?prefix=downloads/greengenes_database/gg_13_5/).


-->
