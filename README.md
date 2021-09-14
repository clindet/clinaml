# clinaml

## Introduction

[CLINAML](https://github.com/clindet/clinaml) is a robust computational method, which can be used to predict the gene expression profiling (GEP)-defined subgroups of AML, as well as the GEP-plus risk scores that bear significant prognostic value in AML.

Currently, the risk stratification of AML is mainly based on the integration of clinical features and genomic abnormalities with prognostic significance, such as the European LeukemiaNet (ELN) recommendations. Meanwhile, several studies constructed prediction models based on differentially expressed genes in AML patients with distinct clinical outcomes and ELN risks.

Although these existing categorization approaches provide a valuable evidence for the molecular classification and therapeutic decision-making of AML, they may ignore the systematic bias and inherent limitations of ELN risk, and may not integrate the GEP-defined subgroups and the high risk alternative splicing.

To refined the risk stratification in AML, firstly, we utilized the multicenter RNA-Seq data of 514 primary AML and conducted the cross-validated unsupervised hierarchical clustering and biological association analysis, by which six major GEP-defined subgroups were identified with distinct clinical and genetic features. To validate these GEP-defined subgroups, we used the TCGA LAML cohort to reproduce the clustering steps, which manifested a similar GEP classification. Furthermore, we also constructed the supervised predictive model based on the Autogluon and XGBoost machine learning methods, which consists of 320 sub-classifiers with different data sampling, and different methods of RNA-Seq reads counting, normalization, and the batch effect adjust.

Via combining high-risk GEP-defined subgroups, we identified deregulated gene expression and alternative splicing events with prognostic value. Finally, based on multiple feature selection methods, we proposed the GEP-plus risk scores using age, GEP-defined subgroups, aberrant gene expression, and splicing events, which could independently predict a subset of AML with an extremely poor prognosis, as well as early death in an independent APL cohort.

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
options(clinaml.condaenv="clinaml") # modified required

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

