# clinaml

## Introduction

[CLINAML](https://github.com/clindet/clinaml) is a robust computational method, which can be used to predict the gene expression profiling (GEP)-defined subgroups of AML, as well as the GEP-plus risk scores.

Currently, the prognostic stratification of $ML is mainly based on the integration of clinical phenotype and genomic abnormalities with clinical significance, such as the European LeukemiaNet (ELN) recommendations. Meanwhile, a few studies constructed the prediction models based on survival status or ELN risk stratification using differential expressed genes.

These classification methods provides a valuable reference for the clinical classification and decision of AML. However, existing approaches may ignored the systematic bias of ELN risk while not taking into account the GEP-defined subgroups and the high risk alternative splicing events.

To refined the risk stratification in AML, firstly, we conducted the cross-validated unsupervised hierarchical clustering and biological association analysis and identified six major GEP-defined subgroups. These molecular subgroups characters with distinct clinical and genetic features in multicenter 514 primary AML. To validate these GEP-defined subgroups, we used the TCGA LAML cohort to reproduce the clustering steps and could get a similar GEP classification. Furthermore, we also constructed the supervised predictive model based on Autogluon and XGBoost machine learning methods. The prediction method consists of 320 sub-classifiers based on the data sampling and the methods of RNA-Seq reads counting, normalization, and the batch effect adjust.

Via adopting the combination of high-risk GEP-defined subgroups, we identified the prognostic gene expression and alternative splicing events. Finally, based on multiple feature selection methods, we proposed the GEP-plus risk scores using the age, GEP-defined subgroups, aberrant gene expression, and splicing events. The risk score could independently predict a poor prognosis in AML, as well as an independent APL cohort.

## Requirements

- R (4.0.2)
- Python (3.8)
- XGBoost (1.4)
- Autogluon (0.2.0)
- Model data (https://www.biosino.org/node/search?queryWord=OEZ007857): Restricted now.

## Installation

```r
git clone https://github.com/clindet/clinaml

cd clinaml

conda env create -f conda_env.yaml

conda activate clinaml

R CMD INSTALL .
```

## Usage

```r
# load package and set important options
library(clinaml)
options(clinaml.MODEL_DATA="/env/model/clinaml") # modified required
options(clinaml.condaenv="clinaml") # modifed required

# XGBoost-based model (DESeq2 and TPM normalization)
x <- system.file("extdata", "tcga_log2tpm_test.csv", package = "clinaml")
res <- xgboost_pred(x, exp_type = "log2tpm")

x2 <- system.file("extdata", "tcga_deseq2_test.csv", package = "clinaml")
res2 <- xgboost_pred(x2, exp_type = "deseq2")

# Autogluon-based model (DESeq2 and TPM normalization)
x <- system.file("extdata", "tcga_log2tpm_test.csv", package = "clinaml")
res <- autogluon_pred(x, exp_type = "log2tpm")
x2 <- system.file("extdata", "tcga_deseq2_test.csv", package = "clinaml")
res2 <- autogluon_pred(x2, exp_type = "deseq2")
```

## Maintainer

[Jianfeng Li](https://github.com/Miachol)

## License

R package:

[MIT](https://en.wikipedia.org/wiki/MIT_License)

Related Other Resources:

Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License

