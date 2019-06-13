---
title: Input Data for MIMOSA2
layout: default
---

# Input data options for the main MIMOSA analysis

MIMOSA has several options for uploading and viewing datasets. Each data file must be formatted as shown in the examples linked below. 

<h2 id="taxonomy">Microbiome Data</h2>

<h3>Taxonomic abundance table</h3>

4 options for providing data on microbiome composition. 

- Sequence variants (ASVs)
<a href="https://elbo-spice.gs.washington.edu/shiny/MIMOSA2shiny/test_seqs.txt" target="_blank">Example (with sequence variants) </a>

- Greengenes 13_5 or 13_8 OTUs
<a href="https://elbo-spice.gs.washington.edu/shiny/MIMOSA2shiny/test_otus.txt" target="_blank">Example (with Greengenes OTUs) </a>

- SILVA v132 OTUs

- No 16S rRNA data; use metagenome data only

<h3>Metagenome abundance table</h3>

- Total KO abundances
Example: 

- Taxon-stratified KO abundances (HUMAnN2 or PICRUSt/PICRUSt2)
Example: 


<!--- <a href="https://elbo-spice.gs.washington.edu/shiny/burrito/Data/examples/example_otus.txt" target="_blank">Example </a> -->


<h2 id="function">Metabolomics Data</h2>

<h3>Metabolite concentration table</h3>

MIMOSA will only analyze data with shared sample IDs between the microbiome and metabolome datasets. 

Compound IDs can be provided as KEGG IDs (used directly for analysis), or as metabolite names, which are mapped to KEGG IDs using MetaboAnalystR (citation). 
However, this process is generally improved by a manual curation step, so we encourage you to run the MetaboAnalystR utility yourself at their website: 
Then you can curate the annotations and use the resulting table for MIMOSA2 analysis.


<a href="https://elbo-spice.gs.washington.edu/shiny/MIMOSA2shiny/test_seqs.txt" target="_blank">Example (with KEGG IDs) </a>

## For advanced users: Providing modifications to the MIMOSA2 network model

A single table can be provided to specify multiple types of modifications to the model. Modifications include adding and removing KOs and/or reactions, and these modifications
can apply to all community members or to specific taxa.
