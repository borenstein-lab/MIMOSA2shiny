---
title: Using the mimosa2 R package
layout: default
---
# Using the mimosa2 R package

The easiest way to run a MIMOSA2 analysis is via the web application. However, you can also install the mimosa2 R package in order
 to run your own custom analyses and integrate MIMOSA2 into analysis pipelines.

## Installation

MIMOSA2 can be easily installed from GitHub using the `devtools` package:

```R
devtools::install_github("borenstein-lab/mimosa2")
``` 

You can test that the package has been installed correctly by running the test suite:

```R
devtools::test("mimosa2")
```

## Downloading reference data

Preprocessed reference data can be downloaded from [this link](http://cnoecker.github.io/MIMOSA2shiny/downloads.html). 

Alternatively, you can run the "generatePreprocessedSpeciesRxnMappings.R" script to obtain MIMOSA2-compatible files from the relevant source database files.

...

## Run a full MIMOSA2 analysis

Once you have downloaded the relevant reference databases, you can run a full MIMOSA2 analysis simply by providing a "configuration table" containing all of the relevant settings for the analysis to the `run_mimosa2` function.
The table below lists the various fields that you can provide in your configuration table. Required fields are in bold.

| Field | Description |
|:------:|:----------:|
|*file1* | Microbiome file path |
|*file2* | Metabolomics file path |
|metagenome | Metagenome file path |
|*database* | Taxonomic abundance file type |
|*metagenome_format* | Taxonomic abundance file type |
|*genomeChoices* | Ref model option |
|simThreshold | If 16S rRNA ASVs are provided, threshold for mapping them to a reference database |
|netAdd | File path to network modifications file |
|met_type | Whether metabolite data is provided as KEGG compound IDs or metabolite names |
|met_transform | metabolome file path |
|rankBased | Whether to use rank-based regression for comparing CMP scores and metabolites (T or F) |
|*kegg_prefix* | metabolome file path |
|*data_prefix* | metabolome file path |
|*vsearch_path* | metabolome file path |

Notes: 
- kegg_prefix is only required when using KEGG-based models.
- file1 is not required if database option 4 is selected and a metagenome file is provided.

Some example configuration tables:

- [An ASV-based analysis using EMBL_GEMS models and rank-based regression](link) 
- [A metagenome-based analysis using KEGG and OLS regression](link2)

You can also download the contribution table used to run any analysis on the MIMOSA2 web server, which allows anyone to later reproduce the same analysis in an R session.

## Run individual components of a MIMOSA2 analysis

## Other utility functions
