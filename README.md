[![DOI](https://zenodo.org/badge/427084835.svg)](https://zenodo.org/badge/latestdoi/427084835)

# BCG Response Subtype predictor (*BRSpred*)

<img src="/img/flavor.png" alt="BRSpred workflow" title="BRSpred" width="100%"/>

## Introduction

The recommended treatment for high-risk non-muscle invasive bladder cancer (HR-NMIBC) patients is intravesical *Bacillus Calmette-Guérin* (BCG) instillations, yet only 50% of patients benefit from BCG therapy. Through molecular profiling of 283 BCG-naïve, HR-NMIBC patients in two separate patient cohorts, with 44 matched post-BCG recurrences, we delineate three distinct BCG Response Subtypes (BRS1-3). Patients with BRS3 tumors have reduced high-grade recurrence-free and progression-free survival as compared to BRS1/2 tumors. BRS3 tumors expressed high EMT-basal markers and have an immunosuppresive profile, which was confirmed with spatial proteomics. The prognostic relevance of the BRS stratification was validated and improved the current recommended clinical guidelines, which is based on clinicopathological variables.

## Contents

The `BRSpred` R-package (see **citation**) contains a nearest-centroid classifier pipeline that can be used for transcriptomic data of primary HR-NMIBC tumors that were treated with BCG. The BCG Response Subtyper (BRS) identifies 3 distinct molecular classes: BRS1, BRS2 and BRS3 tumors. In addition, the package contains a vignette reproducing the key results from the main manuscript, and comes with two datasets: 

* Cohort A (n=132): Discovery Cohort (pre-BCG), exported data object ```CohortA_pre```
* Cohort A (n=44): Post-BCG samples, exported data object ```CohortA_post```
* Cohort B (n=151): Validation Cohort, exported data object ```CohortB```

## Installation

In order to install ```BRSpred``` directly from GitHub:

```
devtools::install_github("CostelloLab/BRSpred", dependencies=TRUE)
```

The ```install_github``` does not build vignettes for the user by default. In order to have the vignettes, one should install the package with the help of ```Rtools``` while adding:

```
devtools::install_github("CostelloLab/BRSpred", dependencies=TRUE, build_vignettes=TRUE)
```

Alternatively, one can install the package from a local source tarball in the R terminal with:
```
install.packages("BRSpred_VERSION.tar.gz", repos=NULL, type="source")
```

Note that one must install dependencies from a R package repository such as CRAN in order to BRSpred to work.

Load the package to the current workspace by:

```
library(BRSpred)
```

## Dependencies

```BRSpred``` requires the following packages for base functionality:

* pamr
* survminer
* matrixStats
* preprocessCore

In order to fully run the vignette reproducing results and/or build the package, the following R package dependencies are also suggested:

* ConsensusClusterPlus,
* pheatmap,
* RColorBrewer,
* tibble,
* roxygen2,
* knitr,
* rmarkdown,

Use of ```Rtools``` may be needed for compiling the package on Windows.

## Usage

The main BRS predictor function is the ```BRS```-function; for example predicting the validation cohort:
```
BRSpred::BRS(BRSpred::CohortB, scale="independent")
```

For further usage and prior results, see the vignette ```BRS prediction```.

## Citation

de Jong FC, Laajala TD, Hoedemaeker RF, Jordan KR, van der Made ACJ, Boevé ER, van der Schoot DKE, Nieuwkamer B, Janssen EAM, Mahmoudi T, Boormans JL, Theodorescu D, Costello JC, Zuiverloon TCM. _Non-muscle-invasive bladder cancer molecular subtypes predict differential response to intravesical Bacillus Calmette-Guérin._ Sci Transl Med. 2023 May 24;15(697):eabn4118. doi: [10.1126/scitranslmed.abn4118](https://doi.org/10.1126/scitranslmed.abn4118). Epub 2023 May 24. PMID: 37224225.

## Reproducing results from de Jong et al.

See the vignette ```BRS prediction``` listed for the package at ```vignette(package="BRSpred")``` for training and validation using the proposed BRS-methodology. The raw contents of this vignette are located in the folder ```/vignettes/brs.Rmd```, or for the compiled package inside ```/doc/brs.html``` (if compiled to HTML as expected by default).

## Acknowledgements

The BCG Response Subtyper was developed in a collaboration between the Erasmus MC Cancer Institute, the Anschutz Medical Campus, Cedars-Sinai Samuel Oschin Comprehensive Cancer Institute and University of Turku & the Finnish Cancer Institute.

## Quick contact

For quicker answers to technical questions regarding the R parts and implementation, feel free to email teelaa@utu.fi / daniel.laajala@helsinki.fi
