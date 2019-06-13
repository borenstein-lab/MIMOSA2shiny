---
title: MIMOSA2
layout: default
---

### Integrative, metabolic model-based analysis of microbiome and metabolomics data

MIMOSA2 relates variation in microbiome composition to paired metabolite measurements, using information from reference databases on the metabolic capabilities of microbial taxa. 
MIMOSA2 seeks to identify differences in microbial metabolic potential that may underlie differences in the concentrations of microbial metabolites between similar communities.

MIMOSA2 can be run as either a web application at https://elbo-spice.gs.washington.edu/shiny/MIMOSA2shiny/, or as a standalone R package available from GitHub. 

### Overview

A MIMOSA2 analysis consists of 4 major steps: 

1) Construct a community metabolic model consisting of the set of metabolic reactions that each community member taxon is predicted to be capable of performing. 

2) Use the metabolic model to calculate metabolic potential (MP) scores for each taxon and metabolite, 
describing an approximate relative estimate of the effects of each taxon, on each metabolite, in each sample.

3) Compare total community-level metabolic potential (CMP) scores with metabolite concentrations across all samples, 
and use a regression model to assess whether CMP scores are significantly predictive of concentrations.

4) For significant metabolites, decompose the overall model fit from step 3 into the contributions from each taxon, 
identifying specific taxa that explain variation in each metabolite.