---
title: MIMOSA2 Reference Data
layout: default
active: downloads
---
# MIMOSA2 Reference Data Downloads

MIMOSA2 can use various types of reference data to construct metabolic models from a set of microbiome data. A subset of these are available for download below. 
The methods used to generate these files are described in the MIMOSA2 manuscript. Scripts to re-generate these files can be found [here](http://github.com/borenstein-lab/MIMOSA2app/scripts/).

## Download all files needed for a specific workflow
See the [R package tutorial](package.html) for an overview of all workflows. 

- 16S rRNA ASVs to EMBL_GEMS
- Greengenes or SILVA OTUs to EMBL_GEMs
- 16S rRNA ASVs to AGORA
- Greengenes or SILVA OTUs to AGORA
- 16S rRNA ASVs to KEGG - see Note
- Greengenes OTUs to KEGG - see Note

To analyze metagenomic data with the KEGG workflow, you must have a KEGG license to download the necessary files, or you can run your analysis using the [web app](https://elbo-spice.gs.washington.edu/shiny/MIMOSA2shiny/)

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


