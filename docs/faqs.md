---
title: Frequently asked questions about MIMOSA2
layout: default
active: faqs
---

# Frequently Asked Questions about MIMOSA2

### Input Data

- [What kind of microbiome data can I use? Can I provide transcriptomic data?](faqs.html#transcriptome)

- [Can I calculate metabolic potential scores even if I don't have metabolomics data?](faqs.html#cmpsAlone)

- [What about metabolomics features that have not been assigned a compound identification?](faqs.html#noID)

- [How should I normalize my metabolomics data, or will MIMOSA2 normalize it for me?](faqs.html#normalization)


### Functionality and Analysis Options

- [How is MIMOSA2 different from MIMOSA1?](faqs.html#mimosa1)

- [Why does MIMOSA2 only analyze some of my metabolites? Why does the set of analyzed metabolites depend on the analysis settings?](faqs.html#metabolites)

- [Which reference database should I use? Why do my results change when using different reference databases?](faqs.html#whichReference)

- [I don't have a KEGG license. Can I still use the KEGG-based analysis?](faqs.html#keggLicense)

### Results Interpretation

- [What does it mean for a metabolite to be negatively predicted by metabolic potential (negative slope)?](faqs.html#negatives)

- [What does it mean for a taxon to have a negative contribution value?](faqs.html#negativeTaxa)

- [What cutoff(s) should I use to identify microbial metabolites and contributing taxa?](faqs.html#thresholds)

### Input Data

<h4 id="transcriptome">What kind of microbiome data can I use? Can I provide transcriptomic data?</h4>
You can use MIMOSA2 with multiple different types of microbiome data, and MIMOSA2 does not require any specific normalization or quantification. 
If you can format your data in a way that complies with the input requirements described on the [Input Data](input.html) page, you can use it. Generally, you need 
either abundances of taxa/populations that can be linked with genomic information, or you need genomic information annotated with KEGG Ortholog abundances. We have not 
comprehensively assessed whether using amplicon data, metagenomic data, or metatranscriptomic data consistently produces better results. 

<h4 id="cmpsAlone">Can I calculate metabolic potential scores even if I don't have metabolomics data?</h4>
Technically, you can do this using the `get_species_cmp_scores` function in the R package. However, we do not recommend it. CMP scores are a very approximate method for 
predicting metabolic processes, and they are most useful as a tool for interrogating possible mechanisms that explain observed metabolite trends, rather than
predicting changes in unknown metabolites.

<h4 id="noID">What about metabolomics features that have not been assigned a compound identification?</h4>
Currently, we can't do anything about them. You can assign a putative identification and see if MIMOSA2 is able to significantly link them with microbiome data.

<h4 id="normalization">How should I normalize my metabolomics data, or will MIMOSA2 normalize it for me?</h4>
MIMOSA2 assumes you have already normalized your metabolite data appropriately for the assay and study design with which it was collected. In particular,
MIMOSA2 assumes that you have quantitative measurements for each individual metabolite feature - i.e. differing measurement values or intensities between samples are proportional with underlying differences in metabolite concentrations. MIMOSA2 will optionally log-transform your data, 
which may improve model fit. 

### Functionality and Analysis Options

<h4 id="mimosa1">How is MIMOSA2 different from MIMOSA1?</h4>
MIMOSA2 and MIMOSA1 answer the same general questions: Which metabolites vary consistently with microbial metabolic potential? What taxa and reactions do they appear to depend on?
MIMOSA2 has expanded input and modeling options and uses a more consistent statistical approach for doing so. It also presents analysis results in a more user-friendly format via the web server.

<h4 id="metabolites">Why does MIMOSA2 only analyze some of my metabolites? Why does the set of analyzed metabolites depend on the analysis settings?</h4>
MIMOSA2 assesses whether your observations are consistent with what we know about microbial metabolism from a specific database or set of metabolic reconstructions.
Therefore, it can only analyze metabolites that are present in the specific reference database, and that are also found somewhere in the community metabolic model 
constructed from your microbiome dataset. If you modify the settings for how the community metabolic model is being constructed, you may change whether certain metabolites
are included in the model.

<h4 id="whichReference">Which reference database should I use? Why do my results change when using different reference databases?</h4>
The best choice of reference data depends on the type of microbiome data you have and the environment of your samples. Here are various considerations that can affect the answer to this question: 
- If you have metagenomic or metatranscriptomic data, currently MIMOSA2 can only utilize the KEGG metabolic network to construct a KEGG metabolic model from your data (this may change in the future). 
- The AGORA database is a collection of metabolic reconstructions of human gut microbiome species, so if your dataset is from an environment other than the human gut, analyses using that option may not be ideal (and may result in a smaller share of taxa included in the community metabolic model).
- Similarly, if you provide 16S rRNA ASV data and use the KEGG model option, MIMOSA2 will links ASVs to KOs via the Greengenes database and PICRUSt 1. Greengenes has not been updated in the last few years, so this is another caveat to keep in mind. (Another way to link ASVs to KEGG models would be to make KEGG Ortholog predictions for your data using PICRUSt2,
which will use a more extensive process and larger database to infer KOs linked to each ASV, and upload the resulting contribution table.)

If you have 16S rRNA amplicon data, there is no reason not to try using multiple reference options and assess the effect on your results. Your results may change because the number of taxa that can be mapped to a model and the 
number of metabolites present in the model can differ widely between the different sources.

<h4 id="keggLicense">I don't have a KEGG license. Can I still use the KEGG-based analysis?</h4>
Yes. An advantage of running MIMOSA2 via the web server is that you can run analyses using KEGG on our server. The web server only provides a subset of the community metabolic network to download, so you
will not be able to run any MIMOSA2 analyses locally unless you have access to files from the KEGG FTP database.

### Results Interpretation

<h4 id="negatives">What does it mean for a metabolite to be negatively predicted by metabolic potential (negative slope)?</h4>
"Negatively predicted" metabolites can be interpreted several different ways. Possible reasons for a negative correlation between a metabolite's levels and its metabolic potential include incorrectly annotated or missing reactions,
and effects beyond direct metabolic reactions such as growth promotion or toxicity. In simulations and in a simple validation dataset, we have found that taxonomic contributors to a model with a negative slope are 
slightly less likely to represent a true taxon-metabolite link than contributors to a model with a positive slope. Contributors identified for negatively correlated metabolites *could* therefore 
represent true relationships, but should be interpreted more cautiously than positively correlated metabolites. 

<h4 id="negativeTaxa">What does it mean for a taxon to have a negative contribution value?</h4>
The contribution values produced by MIMOSA2 are a measure of each taxon's importance in explaining variance in measurements of a specific metabolite. The taxa with the largest contribution values, positive or negative, can be considered the largest potential influencers. A positive contribution value indicates that if that taxon were removed, the metabolite would be less variable 
than is actually observed. A negative contribution value indicates that if that taxon were removed, the metabolite is predicted to be more variable than is actually observed. This means that the negatively-contributing taxon is somehow compensating for or mitigating the predicted metabolic effects of other taxa.

<h4 id="thresholds">What cutoff(s) should I use to identify microbial metabolites and contributing taxa?</h4>
This will depend on analysis goals, study design, and microbiome properties. We typically classify "putative microbe-influenced metabolites" as those predicted with a positive slope and a p-value less than 0.1. The largest taxonomic contributors to variation in those metabolites are considered their "potential taxonomic contributors", typically any taxon with a contribution to variance greater than 1%.

