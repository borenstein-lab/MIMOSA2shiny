---
title: MIMOSA2
layout: default
---

### Integrative, metabolic model-based analysis of microbiome and metabolomics data

MIMOSA2 relates variation in microbiome composition to paired metabolite measurements, using information from reference databases on the metabolic capabilities of microbial taxa. 
MIMOSA2 identifies microbial features that may underlie differences in microbial metabolite levels between similar communities. It answers questions such as:

- Do the metabolites in a dataset appear to vary depending on microbiome composition? Which ones?
- Can differences in microbiome metabolic capabilities explain metabolite variation?
- Which taxa, genes, and reactions appear to be playing a role in metabolite differences?

MIMOSA2 can be run as either a web application at **[http://elbo-spice.cs.tau.ac.il/shiny/MIMOSA2shiny/](http://elbo-spice.cs.tau.ac.il/shiny/MIMOSA2shiny/)**, or as a standalone R package available from GitHub. 

### Overview

![alt text](schematic_v2.png "MIMOSA2 Flow Chart")

A MIMOSA2 analysis consists of 4 major steps: 

1) Construct a community metabolic model consisting of the set of metabolic reactions that each community member taxon is predicted to be capable of performing. 

2) Use the metabolic model to calculate metabolic potential (MP) scores for each taxon and metabolite, 
describing an approximate relative estimate of the effects of each taxon, on each metabolite, in each sample.

3) Compare total community-level metabolic potential (CMP) scores with metabolite measurements across all samples, 
and use a regression model to assess whether CMP scores are significantly predictive of metabolite levels.

4) Decompose the overall model fit from step 3 into the contributions from each taxon, 
identifying specific taxa that explain variation in each metabolite.

For more on each of these steps, see [What Does A MIMOSA2 Analysis Do?](analysis_description.html).