# combinedImputeGEX
Code used to impute tissue-specific gene expression level using genotype and whole blood expression as predictors.


## Notice
The purpose of this repository is to share the code used in the study: "Building an optimal predictive model for imputing tissue-specific gene expression by combining genotype and whole blood transcriptome data".
We are not able to share the data used in the study due to the license issue. We instead include a detailed description on how to prepare and use the data to run our code within the README file.


## Data download
In order to run the code, appropriate data must be downloaded from GTEx v7 database. The data include tissue-specific gene expression profile, whole blood transcriptome profile, and genotype information. 


## Code instructions
Code: impute_gex.R
<br/> 
It is coded in R version 4.1.0 (2021-05-18).


### Input data
1. Tissue-specific gene expression profile
* A dataset contains the expression level of all the genes obtained from a tissue of interest (e.g. Muscle Skeletal, Brain Frontal Cortex, etc).
* The input/GTEx_v7 directory must contain the 47 datasets of transciptome profile corresponding to each of the 47 human tissues.
* Each dataset must be stored in the .Rdata format and named after a target tissue (e.g. Brain_FrontalCortex.Rdata contains the transciptome profile for Brain Spinal Cord tissue).
* The rows must represent samples, while the columns must represent genes.
2. Whole blood transcriptome profile
* A dataset contains the expression level of all the genes obtained from whole blood tissue.
* The input/WBE_datasets dictionary must contain the 47 datasets of whole blood transcriptome profile of samples that also have transcriptome profile for the corresponding tissue.
* Each dataset must be stored in the .Rdata format and named after a target tissue (e.g. Brain_FrontalCortex.Rdata contains the whole blood transciptome profile of samples that also have transcriptome profile for Brain Frontal Cortex tissue).
* The rows must represent samples, while the columns must represent genes.
