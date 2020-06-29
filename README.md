# Paeruginosa_acute_infection 

## Patient data
### patient_to_letter.xlsx
Each patient has an alphabetic label (used in paper) and a numeric label (used for sample collection). This file contains conversion between patient letters and numbers.

### all_rx.xlsx
Contains all drugs administered to the patient in the ICU. Patients 45/E and 46/H did not receive any drugs.

### sample_dates.xlsx
Contains the sputum or stool sample day, relative to day 1.

## Data directories
### patient_trees 
Contains maximum parsimony phylogeny for each patient 

### phylogeny_matrix
Contains a matrix of within-patient calls (SNPs and indels) for each patient used to construct the maximum parsimony tree for each patient.

### MIC_data
Contains raw MIC data for all 420 Pseudomonas aeruginosa isolates sequenced on the Illumina platform, measured on the Vitek machine at the clinical microbiology lab of Boston Children's Hospital. 
#### MIC_data/CLINICAL_MICs.csv 
Contains clinical MIC data used for measuring antibiotic resistance during the acute infection, using a single or handful of clinical isolates that were collected by clinical staff (and not sequenced for this study)

### reference_genomes
Contains fasta and genbank file corresponding to each patient-specific reference genome. Patient reference genomes were constructed from PacBio long-read sequencing of a single isolate from day 1. The genbank files in this directory were used for all analyses in the paper. Fasta file and NCBI annotated versions of all reference genomes are found at Sequence Read Archive PRJNA638217 (https://www.ncbi.nlm.nih.gov/bioproject/PRJNA638217).

### Whole genome sequencing of 420 Pseudomonas aeruginosa isolates using Illumina sequencing are available at https://dataview.ncbi.nlm.nih.gov/object/PRJNA622605

## Code 
#### dMRCA_analysis.ipynb 
Code for the d_MRCA analysis (Fig. 3a), including the null model. 

#### plot_phylogenies.ipynb
Code for plotting all phylogenies including MIC information for each isolate (Fig. 2b, 3c, 4a, Extended Data Fig. 4)

### matlab_scripts (dir containing matlab code)
#### plot_MIC_for_paper.m 
Code for plotting all MICs for individual isolates (Fig. 4b, Extended Data Fig. 7)

#### plot_treatment_MIC_correlation.m 
Code for plotting relationship between frequency of antibiotic therapy (by drug class) and mean fold change in antibiotic resistance (Fig. 4d)

#### figure_genes_under_multiple_selection.m 
Code for plotting heatmap of multiply mutated genes (Extended Data Fig 6)
