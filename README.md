# Robust Order-based Machine Learning (ROML) framework
> Built binary and multi-class classifiers for cancer subtypes


## Overview
Most algorithms operate under the assumption that the training and test data will be drawn from the same distribution. However, heterogeneity (expression shift) exists in different studies due to various sequencing platforms, protocols, materials, etc. Also, cross-platform normalization methods are not effective enough.

Here we developed a general machine learning framework called Robust Order-based Machine Learning (ROML). An order-based Top Scoring Pairs method is firstly used for transforming features into 0-1 binary data based on gene pair orders. The binarized data are then filtered and input into an existing machine learning method for the final predictive model, where interpretable methods with embedded feature selection such as random forest will be preferred.


### Feature selecion
- [kTSP](https://academic.oup.com/bioinformatics/article/21/20/3896/203010) score calculation
- Gene pairs are converted to 0-1 binary features

### Random Forest Models
- Binary 
- Multi-class
  - Model 1: one-vs-rest
  - Model 2: pairwise
  - Model 3: pairwise binary models + objective function


### Simulation

Results to be interpreted...


### Real data application

- Breast cancer
  - RNA-seq: TCGA-BRCA
  - Microarray: MetaBric
  - Subtypes: Basal, Her2, LumA, LumB

- Colorectal cancer
  - RNA-seq: TCGA-COAD
  - Microarray: [KFSYSCC](https://www.synapse.org/#!Synapse:syn4974668)
  - Subtypes: CMS1, CMS2, CMS3, CMS4

---

TO DO (sep 2019):

Modify all CV and multi-class model 3 functions to parallel computing way
