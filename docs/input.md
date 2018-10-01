---
title: Input Data for MIMOSA
layout: default
---
# Input data options for the main MIMOSA analysis

MIMOSA has several options for uploading and viewing datasets. Each data file must be formatted as shown in the examples linked below.

<h2 id="taxonomy">Microbiome Data</h2>

<h3>Taxonomic abundance table</h3>

4 options for providing data on microbiome composition. 

<a href="https://elbo-spice.gs.washington.edu/shiny/MIMOSA2shiny/test_seqs.txt" target="_blank">Example (with sequence variants) </a>

<a href="https://elbo-spice.gs.washington.edu/shiny/MIMOSA2shiny/test_otus.txt" target="_blank">Example (with Greengenes OTUs) </a>
<!--- <a href="https://elbo-spice.gs.washington.edu/shiny/burrito/Data/examples/example_otus.txt" target="_blank">Example </a> -->


<h2 id="function">Metabolomics Data</h2>

<h3>Metabolite concentration table</h3>

MIMOSA will only analyze data with shared sample IDs between the microbiome and metabolome datasets. 

Compound IDs can be provided as KEGG IDs (used directly for analysis), or as metabolite names, which are mapped to KEGG IDs using MetaboAnalystR (link).

<a href="https://elbo-spice.gs.washington.edu/shiny/MIMOSA2shiny/test_seqs.txt" target="_blank">Example (with KEGG IDs) </a>

