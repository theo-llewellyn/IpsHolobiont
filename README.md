# A case of invasive holobionts: invasive Ips typographus brings new fungal phytopathogens into the United Kingdom
This repository contains the code associated to the paper:
<br/>

A case of invasive holobionts: invasive Ips typographus brings new fungal phytopathogens into the United Kingdom
<br/>

**Authors**:

Theo Llewellyn<sup>1,2,#,</sup>*, Angelina Ceballos-Escalera<sup>1,2,#</sup>, Emma Handsup>1,2</sup>, John Richards<sup>1,2</sup>, Yuqi Songsup>1,2</sup>,Daegan J. G. Inward<sup>3</sup>, Alfried P. Vogler<sup>1,2</sup>
<br/>

**Affilitions**<br/>
1. Leverhulme Centre for the Holobiont, Department of Life Sciences, Imperial College London, Silwood Park Campus, Ascot, Berkshire, SL5 7PY, UK
2. Department of Life Sciences, Natural History Museum, Cromwell Road, London, SW6 7BD UK
3. Forest Research, Alice Holt Research Station, Farnham, Surrey, GU10 4LH, UK
<br/>
<sup>#</sup>These authors contributed equally

*Corresponding author: t.llewellyn19@imperial.ac.uk

## Data Records

The processed sequence data and metadata are available on the public NCBI SRA under BioProject accession X.  All other data files are available at Figshare https://doi.org/10.6084/m9.figshare.30913712. Figshare contains:
1. The OTU x Sample table
2. Full Taxonomic Identifications for each OTU as determined by dnabarcoder and BLAST against UNITE and subsequent dynamic clustering
3. An interactive KRONA plot showing taxonomy and abundance of OTUs
4. Representative ITS2 sequences for all OTUs in FASTA format
5. Ophiostomatoid ITS2 sequences used for EPA-ng placement
6. Ophiostomatales reference phylogeny with OTUs grafted
7. Microascales reference phylogeny with OTUs grafted

## Analysis scripts
The following scripts contain all the code to produce the results shown in the manuscript. Sequence denoising, and ASV/OTU identification were performed using the code from Ceballos-Escalera et al (2025), which can be found here: https://github.com/theo-llewellyn/UK_Survey

Phylogenetic Placement uses the phylogenetic trees of Llewellyn et al. (2025) under revision in MycoKeys. Sequence and phylogenetic tree files can be found at FigShare repository https://doi.org/10.6084/m9.figshare.30580976.v1

### 1. Contaminant Removal
1. `./isContaminant.R` This script identifies and removes putative contaminant OTUs associated with negative controls, employing a prevalence-based approach. Follows on from dynamic_clustering.R script of https://github.com/theo-llewellyn/UK_Survey

### 2. Phylogenetic Placement
1. `./phylogenetic_placement.sh` pulls ophiostomatoid OTUs from full dataset, aligns to reference ITS MSAs, places phylogenetically and summarises results using GAPPA
2. `./EPA_dnabarcoder_comparison.R` compares IDs from EPA and those from dynamic clustering and dnabarcoder.

### 3. PIME
1. `./PIME.R` uses prevalence filtering and Random Forests to detect core OTUs for Ips typographus
