---
title: Input Data for MIMOSA2
layout: default
---

# MIMOSA2 Input Data

Generally, MIMOSA2 integrates some form of microbiome dataset with a paired dataset of identified metabolite measurements. It can analyze several different combinations of data types. Each data file must be formatted as shown in the examples linked below. In general,
input data files should be tab-delimited tables with a row for each feature and a column for each sample. There should not be duplicate features or samples. If sample IDs differ between the microbiome and metabolome data files, MIMOSA will only analyze samples that are present in both. 
MIMOSA currently normalizes microbiome data to relative abundances by default. Options for metabolomics data normalization are provided. 


<h2 id="taxonomy">Microbiome Data</h2>

<h3>Taxonomic abundance table</h3>

Microbiome composition data can be provided in 3 different formats: 

1) **Sequence variants (ASVs):** 16S rRNA data processed by denoising tools (i.e. qiime2/DADA2/Deblur). In this case, the first column must contain the DNA sequences themselves. 
If you provide ASV data, you can also control how specifically your sequences are mapped to metabolic models (see [Settings](settings.html)).
<a href="https://elbo-spice.gs.washington.edu/shiny/MIMOSA2shiny/data/formatExamples/test_seqs.txt" target="_blank">Example (with sequence variants) </a>

2) **Greengenes 13_5 or 13_8 OTUs:** 16S rRNA data assigned to closed-reference Greengenes OTUs (e.g. by QIIME1 or vsearch). 
<a href="https://elbo-spice.gs.washington.edu/shiny/MIMOSA2shiny/data/formatExamples/test_gg.txt" target="_blank">Example (with Greengenes OTUs) </a>

3) **SILVA v132 OTUs:** 16S rRNA data assigned to closed-reference SILVA OTUs. **not currently compatible with all metabolic model database options
<a href="https://elbo-spice.gs.washington.edu/shiny/MIMOSA2shiny/data/formatExamples/test_silva.txt" target="_blank">Example (with SILVA OTUs) </a>

4) **No 16S rRNA data; use metagenome data only**: If this option is selected, a metagenome abundance table must be provided (see below). You also must select 
"PICRUSt KO genomes and KEGG metabolic model" as the metabolic model option (see [Running a MIMOSA2 analysis](settings.html) for more details).

<h3>Metagenome abundance table</h3>

A table of functional KEGG Ortholog abundances can be provided instead of *(or in addition to) the taxonomic abundances. This table can take two different forms: 

1) **Total KO abundances**: A table of function abundances, as produced from any metagenomic functional annotation pipeline.
<a href="https://elbo-spice.gs.washington.edu/shiny/MIMOSA2shiny/data/formatExamples/test_kos.txt" target="_blank">Example KO table</a>

2) **Taxon-stratified KO abundances**: A table of taxon-specific function abundances, which can be produced from HUMAnN2, kraken, or PICRUSt/PICRUSt2. This must be formatted in the style
of the "metagenome contribution table" produced by PICRUSt version 1 or 2 (see [the PICRUSt2 documentation](https://github.com/picrust/picrust2/wiki/Full-pipeline-script)). 
You can also use the "format_humann2_contributions" function in [the mimosa2 package](package.html) to correctly format a HuMAnN2 stratified table.
 
<a href="https://elbo-spice.gs.washington.edu/shiny/MIMOSA2shiny/data/formatExamples/test_stratified_kos.txt" target="_blank">Example stratified KO table</a>


<h2 id="function">Metabolomics Data</h2>

<h3>Metabolite concentration table</h3>

MIMOSA only analyzes metabolite features that can be assigned a KEGG compound identification. Compounds can be either be provided directly as KEGG IDs, 
or as metabolite names, which are mapped to KEGG IDs using [MetaboAnalystR](https://www.metaboanalyst.ca/). 
However, this process is generally improved by a manual curation step, so we encourage you to run the MetaboAnalystR utility yourself at [their website](https://www.metaboanalyst.ca/faces/ModuleView.xhtml). 
Then you can manually curate the annotations and use the resulting table for MIMOSA2 analysis.

<a href="https://elbo-spice.gs.washington.edu/shiny/MIMOSA2shiny/data/formatExamples/test_mets.txt" target="_blank">Example metabolites (KEGG IDs) </a>

<a href="https://elbo-spice.gs.washington.edu/shiny/MIMOSA2shiny/data/formatExamples/test_mets_names.txt" target="_blank">Example metabolites (metabolite names) </a>

## For advanced users: Providing modifications to the MIMOSA2 network model

A single table can be provided to specify multiple types of modifications to the model. Modifications include adding and removing KOs and/or reactions, and these modifications
can apply to all community members or to specific taxa. Formatting of this modification file is illustrated in several examples below:

- <a href="https://elbo-spice.gs.washington.edu/shiny/MIMOSA2shiny/data/formatExamples/test_netAdd_rxns_KEGG.txt" target="_blank">Add and remove reactions for all taxa</a>

- <a href="https://elbo-spice.gs.washington.edu/shiny/MIMOSA2shiny/data/formatExamples/test_netAdd_species_rxns_AGORA.txt" target="_blank">Add reactions for specific taxa</a>
 
- <a href="https://elbo-spice.gs.washington.edu/shiny/MIMOSA2shiny/data/formatExamples/test_netAdd_genes_KEGG.txt" target="_blank">Add KEGG orthologs for specific taxa</a>
