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
devtools::install_github("borenstein-lab/mimosa2")
``` 

<!- You can test that the package has been installed correctly by running the test suite:

```R
devtools::test("mimosa2")
```
-->

## Generating preprocessed reference data

In order to run a MIMOSA2 analysis, the first step is to format the reference data you would like to use as needed for the analysis. 

#### Download and reformat metabolic reconstructions
To do so, the first step is to download the reference dataset you would like to use from the appropriate source. The source database options are:

- [AGORA](link)
- [embl_gems](link)
- PICRUSt_KEGG ([PICRUSt1 pre-calculated files], plus KEGG FTP downloads): KEGG network with GreenGenes OTUs and PICRUSt 1
- KEGG (KEGG FTP downloads): KEGG network without taxonomic information, to use with KEGG-annotated metagenomic data

Note that downloading the necessary KEGG FTP files requires a license. 

If you use AGORA 1.0.2 and the current (2019) version of embl_gems, information on each model is included in the package data. Otherwise, if you use a different database version,
you may need to create a new model info file, listing each model, its 16S copy number (if known), and an ID that is shared between the model itself and its linked sequence data. You can see an example for AGORA [here](link, put data up). 

Next, use the `generate_preprocessed_networks` function to reformat the model files to be compatible with MIMOSA2. For example, to format the *embl_gems* database models:

```R
generate_preprocessed_networks("embl_gems", dat_path = file_path_to_raw_models, out_path = file_path_for_output)
```

(In the function call above, you would replace `file_path_to_raw_models` and `file_path_for_output` with your corresponding file paths.) 

#### Download reference ribosomal sequence data (ASVs)
Next, if you intend to use MIMOSA to analyze ASV data, you will need to download ribosomal sequences linked to each model. 

For AGORA v1.0.2 and the 2019 version of embl_gems, you can do so using the package function `download_ribosomal_ref_seqs`, which uses the [biomartR](https://ropensci.github.io/biomartr/) package to download the relevant list of accessions from NCBI. 

```R
download_ribosomal_ref_seqs("AGORA")

```
To use GreenGenes and KEGG, you can simply download the GreenGenes representative OTU sequence file from the [QIIME2 website](https://docs.qiime2.org/2019.4/data-resources/#marker-gene-reference-databases).

Additionally, to use an ASV table as input, you will need to make sure you have [vsearch](https://github.com/torognes/vsearch) installed. 

[Preprocessed reference data for the metabolic reconstruction databases can also be downloaded from the [Downloads page](downloads.html). ]

#### Download OTU-model mapping files
If you would like to link GreenGenes or SILVA reference OTUs to AGORA or embl_gems models, you can skip downloading the sequence data and instead download our pre-generated mappings between these databases from the [Downloads page](downloads.html).


## Run a full MIMOSA2 analysis

Once you have downloaded and set up the relevant reference databases, you can run a full MIMOSA2 analysis simply by providing a "configuration table" containing all of the relevant settings for the analysis to the `run_mimosa2` function.
The table below lists the various fields that you can provide in your configuration table. Required fields are in bold.

| Field | Description |
|------|----------|
|**file1** | Microbiome file path |
|**file2** | Metabolomics file path |
|**database** | Taxonomic abundance file type |
|metagenome_format | Metagenome function abundance file type |
|**genomeChoices** | Ref model option |
|simThreshold | If 16S rRNA ASVs are provided, threshold for mapping them to a reference database |
|netAdd | File path to network modifications file |
|met_type | Whether metabolite data is provided as KEGG compound IDs or metabolite names |
|met_transform | metabolome file path |
|rankBased | Whether to use rank-based regression for comparing CMP scores and metabolites (T or F) |
|**kegg_prefix** | File path to processed generic KEGG network - product of the generate_preprocessed_networks function above database files|
|**data_prefix** | File path to other reference databases |
|**vsearch_path** | File path to vsearch executable |

Notes: 
- **kegg_prefix** is only required when using KEGG-based models.
- **file1** is not required if database option 4 is selected and a metagenome file is provided.

Some example configuration tables:

- [An ASV-based analysis using EMBL_GEMS models and rank-based regression](link) 
- [A metagenome-based analysis using KEGG and OLS regression](link2)

You can also download the contribution table used to run any analysis on the MIMOSA2 web server, which allows anyone to later reproduce the same analysis in an R session.

Once you have set up a configuration table, it is easy to run a full MIMOSA2 analysis: 

```R
mimosa_results = run_mimosa2(config_table)
```

## Run individual components of a MIMOSA2 analysis

## Other utility functions

- `format_humann2_contributions`

- `map_to_kegg`

- plot functions