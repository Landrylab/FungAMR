# FungAMR

Repository containing data analysis scripts and figures for the FungAMR paper. A searchable version of Supplementary Table 1 is available at the [Comprehensive Antibiotic Resistance Database](http://card.mcmaster.ca).

## Pre-Print

bioRxiv: [FungAMR: A comprehensive portrait of antimicrobial resistance mutations in fungi](https://www.biorxiv.org/content/10.1101/2024.10.07.617009v1).

## Description of the confidence scores associated with curated mutations in FungAMR

| Confidence Score | Description |
|-----|-----|
| 1 | The mutation was created in the same susceptible species background and effect on resistance was confirmed by measurements. |
| 2 | The mutation was created in the same gene but expressed in another species with potentially non endogenous level of expression and effect on resistance was confirmed by experimental measurements. |
| 3 | The effect of the mutation was quantified by bulk competition assays such as Deep Mutational Scanning. |
| 4 | The mutation was identified by experimental evolution where its high frequency in independent replicates suggests it causes resistance. |
| 5 | The deletion of the gene causes resistance. |
| 6 | The overexpression or duplication of the gene causes resistance. |
| 7 | Significant association between the mutation and resistance in a population as determined by GWAS (or other population-based association such as QTL or classical genetics). |
| 8 | The mutation was identified in a 'natural' strain (e.g. a clinical isolate) that is resistant but without any further validations. |
| Negative scores (-) | Any confidence score where the mutation was found not to cause resistance using the approach of the corresponding positive score. |
