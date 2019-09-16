---
title: MIMOSA2 Analysis Settings
layout: default
active: settings
---
# Running a MIMOSA2 Analysis: Setup and Settings

The MIMOSA2 web application is available at [https://elbo-spice.gs.washington.edu/shiny/MIMOSA2shiny](https://elbo-spice.gs.washington.edu/shiny/MIMOSA2shiny)

The landing page of the web application is an interface for uploading your data and selecting options for your MIMOSA2 analysis. Each input section and the corresponding
options are described below:

## Microbiome Data
Select your microbiome data format, and upload the corresponding data file. You can provide taxonomic abundances, function abundances, or both, in several formats. 
See [Input Data](input.html) for more information.

## Metabolic Model Options
Select what metabolic database you would like to use to link your microbiome and metabolite data. Several options are available:

- KEGG metabolic model: As in the original version of MIMOSA, this option uses pre-computed KEGG Ortholog inferences 
and a generic KEGG metabolic model to link taxa to metabolites. If you provide a taxonomic abundance table, the input taxa will first be mapped to Greengenes OTUs, and the pre-computed PICRUSt KO predictions for each OTU will be used 
as the basis for the community metabolic model. 
If you provide KEGG Ortholog metagenome abundances, these will be directly linked to reactions and metabolites.

- AGORA genomes and models (Magnusdottir et al Nature Biotech 2016): For this option, the 16S rRNA sequences of provided microbiome taxa (either from Greengenes/SILVA, or the ASV sequences themselves) are mapped to 
ribosomal genes from the genomes of the AGORA collection of metabolic reconstructions (currently consisting of 818 reference strains present in human gut microbiota).
The metabolic reconstructions for mapped species are then combined into the community metabolic model and used for the analysis.

- RefSeq/EMBL_GEMs genomes and models ([Machado et al Nucleic Acids Research 2018](https://doi.org/10.1093/nar/gky537)): An alternative collection of metabolic reconstructions, corresponding to 
all 5,587 RefSeq reference genomes and constructed using the CarveMe method. As in the AGORA option, MIMOSA2 maps the 16S rRNA sequences of input taxa 
against a reference ribosomal database to link microbiome features to reference genomes and metabolic reconstructions.


If you have provided a table of 16S rRNA ASVs, you can select how strictly these sequences are aligned against the reference sequences. The default is a relatively strict minimum threshold (99%). 
Mappings are pre-computed for reference OTU sequences at a 97% threshold.

Finally, advanced users can optionally provide a file specifying modifications to the network model, including adding custom reactions and filtering specific reference reactions. See [Input Data](input.html). 

## Metabolomics Data
Select your metabolite data format, and upload the corresponding data file. Metabolite data can be generated using any platform but should consist of identified metabolite abundances (preferably specified as KEGG Compound IDs)

Select whether the metabolite data should be log transformed. We recommend this option if your metabolite measurements tend to have highly skewed distributions (for example, if your study consists of two sample groups and one has greater variability than the other). 
You may find it useful to examine the CMP-Metabolite comparison results with and without log transformation.

## Algorithm Settings

Select whether to compare metabolic potential (CMP) and metabolite levels using ordinary least-squares regression (OLS) or rank-based regression. 
We generally recommend rank-based regression as it can detect metabolite relationships more robustly and sensitively across a wider variety of data distributions. However, 
for rank-based regression, the contributions of individual taxa to the model fit are calculated using a permutation-based approach, which greatly increases the analysis runtime. Therefore, when this option is selected, 
MIMOSA2 will only calculate taxonomic contributors for metabolites with a model p-value less than 0.1 (rather than all metabolites). You also have the option to skip 
the contribution analysis and just compare metabolites with metabolic potential.
 
Rank-based regression is provided using the Rfit package in R. The general idea is to find a regression solution that minimizes a function of both the rank and the size of the model residuals (instead of the sum of squared errors in OLS). You can read more about the statistical approach [here](https://journal.r-project.org/archive/2012-2/RJournal_2012-2_Kloke+McKean.pdf). 
Other non-linear models may be provided as options in the future.

## Run MIMOSA
Once you have selected all options, push this button to run your analysis. Analysis runtime depends on the numbers of samples, microbiome taxa, and metabolites in your dataset. Analyses using OLS regression typically only take a minute or two; analyses using rank-based regression may take ~10 minutes.