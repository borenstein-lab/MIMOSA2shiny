---
title: MIMOSA2 Reference Data
layout: default
active: downloads
---
# MIMOSA2 Reference Data Downloads

MIMOSA2 can use various types of reference data to construct metabolic models from a set of microbiome data. A subset of these are available for download below. 
The methods used to generate these files are described in the MIMOSA2 manuscript. Scripts to re-generate these files can be found in the "scripts" directory of the MIMOSA2 app GitHub repository.

## Download all files needed for a specific workflow
See the [R package tutorial](package.html) for an overview of all workflows. 

- [16S rRNA ASVs to EMBL_GEMS](http://elbo-spice.cs.tau.ac.il/shiny/MIMOSA2shiny/refData/ASV_EMBL.tar.gz)
- [Greengenes or SILVA OTUs to EMBL_GEMs](http://elbo-spice.cs.tau.ac.il/shiny/MIMOSA2shiny/refData/OTU_EMBL.tar.gz)
- [16S rRNA ASVs to AGORA](http://elbo-spice.cs.tau.ac.il/shiny/MIMOSA2shiny/refData/ASV_AGORA.tar.gz)
- [Greengenes or SILVA OTUs to AGORA](http://elbo-spice.cs.tau.ac.il/shiny/MIMOSA2shiny/refData/OTU_AGORA.tar.gz)
- 16S rRNA ASVs to KEGG - see Note
- Greengenes OTUs to KEGG - see Note
- Metagenomes to KEGG - see Note

Note: To use a KEGG-based workflow, you must have a KEGG license to download the necessary files, or you can run your analysis using the [web app.](https://elbo-spice.gs.washington.edu/shiny/MIMOSA2shiny/) 
The specific files and their setup are described in the [R package tutorial](package.html#KEGG). 

You can also download freely available data for a particular workflow with the `download_reference_data` function in MIMOSA2 (see the [R package tutorial](package.html)).

If these links don't work, try these backup links to access the [EMBL_GEMS](https://drive.google.com/drive/folders/1iCuiK1DOujAe9TFLW_SVlskqOJUU_quk?usp=share_link) and [AGORA](https://drive.google.com/drive/folders/1HfcOIo1RsLdFgR_zgaNze8kIKpX-n8sd?usp=share_link) databases.

<!---
Individual processed data files are described and linked below:

### Ribosomal sequence reference datasets (for mapping ASVs)
These are provided in .udb format for efficient download and search with vsearch.

- [GreenGenes 99% OTU representative sequences]() These are the provided Greengenes 13_8 representative 99% OTUs, indexed into a vsearch udb file.

- [Ribosomal sequences from AGORA genomes]()

- [Ribosomal sequences from RefSeq (for EMBL_GEMS)]()

### Processed mappings for linking OTUs to genomes and models:

- [Mapping from Greengenes 13_8 99% OTUs to AGORA-linked genomes]()

- [Mapping from Greengenes 13_8 99% OTUs to RefSeq genomes for EMBL_GEMS]()

- [Mapping from SILVA 132 99% OTUs to AGORA-linked genomes]()

- [Mapping from SILVA 132 99% OTUs to RefSeq genomes]()

The scripts used to generate these mapping files are available at http://github.com/cnoecker/MIMOSA2shiny/scripts.

### Taxon-specific metabolic models, processed and formatted for MIMOSA2 compatibility: 

- [AGORA models (reconstructions from Magnusdottir et al Nature Biotech 2016)](agora_models.tar.gz)

- [CarveMe/EMBL_GEMS models (reconstructions from Machado et al Nucleic Acids Research 2018)](embl_gems_models.tar.gz)


-->
