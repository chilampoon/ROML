# Robust Order-based Machine Learning ROML framework

Built binary and multi-class classifiers for cancer subtypes


### Feature selecion
- kTSP score calculation
- Gene pairs are converted to 0-1 binary features

### Random Forest Models
- Binary 
- Multi-class
  - Model 1: one-vs-rest
  - Model 2: pairwise
  - Model 3: pairwise binary models + objective function


### Simulation

The results were not as we expected...


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
