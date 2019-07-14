# Binary classification analysis


## Real Data
### Data description and preprocession

_Maybe adding more datasets afterwards._

#### TCGA BRCA RNA-seq

1. Data were download from GEO [GSE62944](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62944), 
and were processed as stated in this [paper](https://www.nature.com/articles/s41598-018-25357-0).

###### P.S. They select one primary tumoral sample from one patient only, thus the ids are all unique, I don't know why they did that.
###### Ref: [Alternative preprocessing of RNA-Sequencing data in The Cancer Genome Atlas leads to improved analysis results](https://academic.oup.com/bioinformatics/article/31/22/3666/240143)

2. Use feature counts after median ratio normalization and variance stabilizing transformation (vst) to redo PAM50 breast cancer subtypes assignment, 
to see if the results are different from those using TPM data a lot.

###### Ref: [MLSeq manual](https://bioconductor.org/packages/release/bioc/vignettes/MLSeq/inst/doc/MLSeq.pdf)

3. Use values processed above to apply ML, random forest first, doing cross validation and testing on independent MetaBric data. Then do it again reversely.
For binary classification, LumA and LumB balanced samples were selected.
    - Baseline model: DE analysis to identify DEGs, then do ML procedures
    - ROML: kTSP + a ML model
    - Compare performances, choose the better one as seeds of simulation
    
#### MetaBric microarray

1. Data were downloaded from Synapse sofware platform (syn1688369) and preprocessed as mentioned above.
2. Same step as that in TCGA.
3. Same step as that in TCGA.


---

## Simulation

### Generate RNA-seq read counts



### Generate microarray values




