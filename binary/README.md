# Binary classification analysis

## 1. Simulation

> __GOAL__: To prove our expectation: ROML predicts better - more accurate and more stable - when inter-study heterogeneity  increases

### Data Generation

1. Subsample balanced samples of LumA & LumB subtypes from TCGA BRCA and MetaBric data
2. Differential expression analysis for MetaBric suσbdataset by limma-trend
3. Set the DE percentage `pDE`
3. Values in training data sample from N(μ, σ) - μ, σ are the means and variances of genes in the MetaBric data; DEG and non-DEG values are generated respectively in two groups
4. Same as step 3, values in testing are based on BRCA data


### Balanced dataset

|             |Training|Testing|
|:-----------:|:------:|:-----:|
|No. of gene  |5000    |5000   |
|No. of sample|200+200 |200+200|
|DEG%         |10%, 30%|10%, 30%|



Compare accuracy (ACC) and Youden index of random forest (RF), ktsp + random forest (tspRF/ROML) and [ktsp](https://academic.oup.com/bioinformatics/article/31/2/273/2365798) methods.


### Imbalanced dataset

Just wait and see if it's worth exploring.



## 2. Real Data

### Data description 

_Maybe adding more datasets afterwards._

#### 1. TCGA BRCA RNA-seq

- Data were download from GEO [GSE62944](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62944), data preprocession and BRCA PAM50 prediction see this [paper](https://www.nature.com/articles/s41598-018-25357-0)

###### P.S. They select one primary tumoral sample from one patient only, thus the ids are all unique, I don't know why they did that.

#### 2. MetaBric microarray

- Data were downloaded from Synapse sofware platform (syn1688369) and preprocessed as mentioned above.

_Update: After considering about the unit issues a lot, I think it's better to use tpm values to do baseline & ROML, since I don't have enough information and it's too much work to start from the very beginning like BAM files or fastq._


### Pipeline

We tried both directions for prediction, i.e. `TCGA -> MB` & `MB -> TCGA`

#### ROML

- Filter out 50% (or 25%) genes with relative low expression using mean ranks
- Log-transformed TPM values from RNA-seq and fluorescent intensities from microarray are used to calculate tsp score
- The ranks of top gene pairs are transformed into 0-1 binary values and inputed into the following machine learning algorithm


#### Baseline

- Filter out 50% (or 25%) genes with relative low expression using mean ranks
- Log-transformed TPM values from RNA-seq and fluorescent intensities from microarray are used to conduct differential expression analysis
- Features are selected based on several thresholds of adjusted P-value and log fold change in the result of DEA
- Perform cross-platform normalization (at least 3 various methods): ----> TO DO!
- Go through machine learning procedure with selected features

### kTSP
--------> TO DO!
- 


---

7/22 To do:

1. Read papers about __cross-platform normalization__, summarize into my wiki page, then include 3 typical ones in my analysis
2. Figure out how to improve binary real data steps
3. Add kTSP for comparison
4. Finish prelimary binary part! (by the end of this week 7/28)



  
