# clinaml

## Introduction

[CLINAML](https://github.com/clindet/clinaml) is a robust computational method, which can be used to predict the gene expression profiling (GEP)-defined subgroups of AML, as well as the GEP-plus risk scores that bear significant prognostic value in AML.

Currently, the risk stratification of AML is mainly based on the integration of clinical features and genomic abnormalities with prognostic significance, such as the European LeukemiaNet (ELN) recommendations. Meanwhile, several studies constructed prediction models based on differentially expressed genes in AML patients with distinct clinical outcomes and ELN risks.

Although these existing categorization approaches provide a valuable evidence for the molecular classification and therapeutic decision-making of AML, they may ignore the systematic bias and inherent limitations of ELN risk, and may not integrate the GEP-defined subgroups and the high risk alternative splicing.

To refined the risk stratification in AML, firstly, we utilized the multicenter RNA-Seq data of 514 primary AML and conducted the cross-validated unsupervised hierarchical clustering and biological association analysis, by which six major GEP-defined subgroups were identified with distinct clinical and genetic features. To validate these GEP-defined subgroups, we used the TCGA LAML cohort to reproduce the clustering steps, which manifested a similar GEP classification. Furthermore, we also constructed the supervised predictive model based on the Autogluon and XGBoost machine learning methods, which consists of 320 sub-classifiers with different data sampling, and different methods of RNA-Seq reads counting, normalization, and the batch effect adjust.

Via combining high-risk GEP-defined subgroups, we identified deregulated gene expression and alternative splicing events with prognostic value. Finally, based on multiple feature selection methods, we proposed the GEP-plus risk scores using age, GEP-defined subgroups, aberrant gene expression, and splicing events, which could independently predict a subset of AML with an extremely poor prognosis, as well as early death in an independent APL cohort.

## Note

A new online version of the prediction models based on 655 AMLs are provided in http://hiplot.org/clinical-tools/clinaml-gep2, which was compatible with the predictions of this repository. It can be used to get G1-G8 labels of AML:

- G1: PML::RARA
- G2: CBFB::MYH11
- G3: RUNX1::RUNXT1
- G4: biCEBPA/-like
- G5: myelodysplasia-related/-like
- G6: HOX-committed
- G7: HOX-primitive
- G8: HOX-mixed

G1-G8 was associated with distinct prognosis and drug sensitivities, supporting the clinical applicability of this transcriptome-based classification of AML. These molecular subgroups illuminate the complex molecular network of AML, which may promote systematic studies of disease pathogenesis and foster the screening of targeted agents based on omics.

The predicted labels of TCGA LAML and Beat AML cohorts could be directly downloaded from [here](https://www.pnas.org/doi/suppl/10.1073/pnas.2211429119/suppl_file/pnas.2211429119.sd09.xlsx).

We would like to invite you to cite us paper:

- Cheng WY, Li JF, Zhu YM, Lin XJ, Wen LJ, Zhang F, Zhang YL, Zhao M, Fang H, Wang SY, Lin XJ, Qiao N, Yin W, Zhang JN, Dai YT, Jiang L, Sun XJ, Xu Y, Zhang TT, Chen SN, Zhu HH, Chen Z, Jin J, Wu DP, Shen Y, Chen SJ. Transcriptome-based molecular subtypes and differentiation hierarchies improve the classification framework of acute myeloid leukemia. Proc Natl Acad Sci U S A. 2022 Dec 6;119(49):e2211429119. doi: 10.1073/pnas.2211429119. Epub 2022 Nov 28. PMID: 36442087.

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

R CMD INSTALL .

python3 -m pip install -U pip
python3 -m pip install -U setuptools wheel
python3 -m pip install -U "mxnet<2.0.0"
python3 -m pip install autogluon==0.2.0
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

Web interface based on the Hiplot (ORG) platform is available here: https://hiplot.org/clinical-tools/clinaml-gep. Users can get the following output in minutes.

## Maintainer

[Jianfeng Li](https://github.com/Miachol)

## License

R package:

[MIT](https://en.wikipedia.org/wiki/MIT_License)

Related Other Resources:

Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License

