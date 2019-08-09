---
title: What does MIMOSA2 do?
layout: default
active: analysis_description
---
# What does MIMOSA2 do?

MIMOSA2 summarizes paired microbiome-metabolome datasets to support mechanistic interpretation and hypothesis generation. It answers 2 main questions:

1) Do the concentrations of metabolites in my dataset vary in a way that is consistent with what we would expect, based on the varying capabilities of the associated microbiomes?

2) Which specific microbial taxa appear to be responsible for differences in the concentrations of each metabolite across samples?

MIMOSA2 uses a 4-step workflow, described below, to answer these questions. A more technical explanation of the method can be found in the [MIMOSA2 manuscript](link).

**1) Construct a community metabolic model consisting of the set of metabolic reactions that each community member taxon is predicted to be capable of performing.**


2) Use the metabolic model to calculate metabolic potential (MP) scores for each taxon and metabolite, 
describing an approximate relative estimate of the effects of each taxon, on each metabolite, in each sample.

3) Compare total community-level metabolic potential (CMP) scores with metabolite concentrations across all samples, 
and use a regression model to assess whether CMP scores are significantly predictive of concentrations.

4) Decompose the overall model fit from step 3 into the contributions from each taxon, 
identifying specific taxa that explain variation in each metabolite.

