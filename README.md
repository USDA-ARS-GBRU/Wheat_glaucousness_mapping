# Wheat_glaucousness_mapping

**Daniela's working directory for QTL and GWAS project**

## HilliardxGA06493-13LE6 RIL population
*Linkage Map Construction and QTL Analysis*

Input genotypic SNP data file:

**SunRILs_2021_UX1992_RILs_filt.vcf.gz**


Allele definitions:

'Hilliard' coded as parent A

'GA06493-13LE6' coded as parent B


### Linkage map construction 

A genetic linkage map was constructed from 332 RIL individuals using R packages:

*qtl version 1.48.1*

*magrittr version 2.0.1*

*ASMap version 1.0.4*

as described in R code file:

**HilliardxGA_RIL_GeneticMap-Lite-2023-Final.R**


Linkage map QC

Before finalizing the linkage map, dot plots of each linkage group were constructed for quality control purposes.

semifinal dotplot file name:

UX1992-HilliardxGA-Miller-linkagemap-dotplots-semifinal.jpg or pdf

3 outliers were idenfied visually from these dot plots, on chromosomes 1A and 1B, and were removed from the final map.


Final linkage map 

file name:

**UX1992-HilliardxGA-RIL-Miller-linkage-map.csv**

dot plot of final linkage map:

UX1992-HilliardxGA-Miller-linkagemap-dotplots.jpg and pdf



### Phenotypic Means

HG RIL Population Input Phenotypic Data Files:

**R21_T20-T25_UX1992_final_DMM.csv** # Raleigh 2021

**K21-WaxExp_T6-T10_final_DMM.csv** # Kinston 2021


Calculation of phenotypic means among genotype replicates within an environment was done using a mixed model with R packages:

*lme4 version 1.1.27.1*

*emmeans version 1.6.3*

As described in R code file: 
**lsmeans-GHpop-all-locs-combined.R**

This concatenated all phenotype means in a single file for QTL analysis:
HG_R21_K21_alltrts.csv


### Quantitative Trait Locus (QTL) Mapping

QTL mapping was done using phenotypic data for 205 individuals in Raleigh, NC and 189 individuals in Kinston, NC using R packages:

*qtl version 1.48.1*

*magrittr version 2.0.1*

*ASMap version 1.0.4*

As described in R code file: 
**QTL_mapping_HG_RIL_Miller_final_copy.R**
