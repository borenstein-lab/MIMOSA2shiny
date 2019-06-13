---
title: Using the mimosa2 R package
layout: default
---
# Using the mimosa2 R package

## Installation

MIMOSA2 can be easily installed from GitHub using the devtools package:

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
The table below lists the various fields that you can provide in your configuration table:


| Field | Description |
|:------:|:----------:|
|file1 | microbiome file path |
|file2 | metabolome file path |
| col 2 is      | centered      |
| zebra stripes | are neat      |

Some example configuration tables:

- An ASV-based analysis using EMBL_GEMS models: [here](link) 
- A metagenome-based analysis using KEGG: [here](link2)

You can also download the contribution table used to run any analysis on the MIMOSA2 web server, which allows anyone to later reproduce the same analysis in an R session.

## Run individual components of a MIMOSA2 analysis

## Other utility functions
