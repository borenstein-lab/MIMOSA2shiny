---
title: Frequently asked questions about MIMOSA2
layout: default
active: faqs
---

# Frequently Asked Questions about MIMOSA2
Make this have its own sidebar menu

- **What kind of microbiome data can I use? Can I provide transcriptomic data?**
You can use MIMOSA2 with multiple different types of microbiome data, and MIMOSA2 does not require any specific normalization or quantification. 
If you can format your data in a way that complies with the input requirements described on the [Input Data](input.html) page, you can use it. Generally, you need 
either abundances of taxa/populations that can be linked with genomic information, or you need genomic information annotated with KEGG Ortholog abundances. We have not 
comprehensively assessed whether using amplicon data, metagenomic data, or metatranscriptomic data consistently produces better results. 

- **How is MIMOSA2 different from MIMOSA1?**
MIMOSA2 and MIMOSA1 answer the same general questions: Which metabolites vary consistently with microbial metabolic potential? What taxa and reactions do they appear to depend on?
MIMOSA2 has expanded input and modeling options and uses a more consistent statistical approach for doing so. It also presents analysis results in a more user-friendly format via the web server.
See [the manuscript](link) for a more detailed comparison of the two.

- **Why does MIMOSA2 only analyze some of my metabolites? Why does the set of analyzed metabolites change depending on the analysis settings?**
MIMOSA2 assesses whether your observations are consistent with what we know about microbial metabolism from a specific database or set of metabolic reconstructions.
Therefore, it can only analyze metabolites that are present in the specific reference database, and that are also found somewhere in the community metabolic model 
constructed from your microbiome dataset. If you modify the settings for how the community metabolic model is being constructed, you may change whether certain metabolites
are included in the model.

- **Which reference database should I use? Why do my results change when using different reference databases?**
The best choice of reference data depends on the type of microbiome data you have and the environment of your samples. If you have metagenomic or metatranscriptomic data, currently
MIMOSA2 can only utilize the KEGG metabolic network to construct a KEGG metabolic model from your data (this may change in the future). If you have 16S rRNA amplicon data,
you may want to try using multiple reference datasets and assess how/whether your results change. The AGORA database is a collection of metabolic reconstructions of human gut microbiome species,
so if your dataset is from an environment other than the human gut, analyses using that option may not be ideal (and may result in a smaller share of taxa included in the community metabolic model).

- **I don't have a KEGG license. Can I still use the KEGG-based analysis?**
Yes. An advantage of running MIMOSA2 via the web server is that you can run analyses using KEGG on our server. The web server only makes a subset of the community metabolic network available for download, so you
will not be able to run any MIMOSA2 analyses locally unless you have access to files from the KEGG FTP database.

- **Can I calculate metabolic potential scores even if I don't have metabolomics data?**
Sort of.

- **What about metabolomics features that have not been assigned a compound identification?**
Currently, we can't do anything about them. You can assign a putative identification and see if MIMOSA is able to link them to anything interesting.

- **What does it mean for a metabolite to be negatively predicted by metabolic potential (negative slope)?**
It's complicated...

