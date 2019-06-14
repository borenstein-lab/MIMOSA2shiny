---
title: MIMOSA2 Analysis Settings
layout: default
---
# Running a MIMOSA2 Analysis: Setup and Settings

The MIMOSA2 web application is available at https://elbo-spice.gs.washington.edu/shiny/MIMOSA2shiny

The landing page of the web application is an interface for uploading your data and selecting options for your MIMOSA2 analysis. Each input section and the corresponding
options are described below:

## Microbiome Data Options
Select your microbiome data format. You can provide taxonomic abundances, function abundances, or both, in several formats.

## Metabolic Model Options
Select what metabolic database you would like to use to link your microbiome and metabolite data. Several options are available:
- 

More in-depth information on these resources is available on the [Downloads page](downloads.html). 

If you have provided a table of 16S rRNA sequence variants, 
Mapping of ASVs to ref dbs
For advanced users: Optional modifications to network (include option to provide whole network?)

## Metabolomics Data Options
Identifications: Best to provide KEGG IDs
Log transform - whether to do so - particularly important for OLS regression

## Algorithm Options
OLS or rank-based regression
Rank-based regression is provided using the Rfit package and the methods proposed by Jaeckel, Kloke, etc (citations).