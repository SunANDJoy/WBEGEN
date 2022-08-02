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


### Input 
1. Tissue-specific gene expression profile (imputation target)
* A dataset contains the expression level of all the genes obtained from a tissue of interest (e.g. Muscle Skeletal, Brain Frontal Cortex, etc).
* The input/GTEx_v7 directory must contain the 47 datasets of transcriptome profile corresponding to each of the 47 human tissues.
* Each dataset must be saved in the .Rdata format and named after the corresponding tissue (e.g. Brain_FrontalCortex.Rdata contains the transcriptome profile for Brain Frontal Cortex tissue).
* The rows must represent samples, while the columns must represent genes.
2. Whole blood transcriptome profile
* A dataset contains the expression level of all the genes obtained from whole blood tissue.
* The input/WBE_datasets dictionary must contain the 47 datasets of whole blood transcriptome profile of samples that also have transcriptome profile for the corresponding tissue.
* Each dataset must be saved in the .Rdata format and named after the corresponding tissue (e.g. Brain_FrontalCortex.Rdata contains the whole blood transcriptome profile of samples that also have transcriptome profile for Brain Frontal Cortex tissue).
* The rows must represent samples, while the columns must represent genes.
3. Genotype information
* A dataset contains the genotype information of cis-SNPS located around a target gene.
* The directory input/WBE_datasets must contain a number of datasets as many as the number of all available genes present in all 47 tissues.
* Each dataset must be saved in the .Rdata format and named after the corresponding gene (e.g. ENSG00000089022.Rdata contains the genotype information of cis-SNPs located around the ENSG00000089022 gene).
* The rows must represent samples, while the columns must represent SNPs.
4. Input information (input_info.Rdata)
* The file contains a list of the 47 human tissues and the sets of sample IDs for the training, validation, and test set.
* To ensure the consistency of the results, we defined the components for the training, validation, and test set in advance.
* We provide the input_info.Rdata file in this repository.
<br/>
The three directories (GTEx_v7, WBE_datasets, and Genotype_datasets) containing the appropriate datasets and one file (input_info.Rdata) must be present in the 'input' directory in order to run the code.


### Output 
1. The summary of the gene expression imputation results across the 47 human tissues (output.Rdata).
* For the metric of imputation accuracy, we used Mean R^2 (mean of R^2 between true expression level and imputed expression level over all available genes belonging to each tissue)
* Imputation accuracy for two standalone models (GEN and WBE) and two combined models (MERGED and IVW) are provided.
* Proportions of genes for which each combined model was applied are also provided.
