# BCG Response Subtype predictor (*BRSpred*)

<img src="/img/flavor.png" alt="BRSpred workflow" title="BRSpred" width="100%"/>

## Introduction

The recommended treatment for high-risk non-muscle invasive bladder cancer (HR-NMIBC) patients is intravesical *Bacillus Calmette-Guérin* (BCG) instillations, yet only 50% of patients benefit from BCG therapy. Through molecular profiling of 283 BCG-naïve, HR-NMIBC patients in two separate patient cohorts, with 44 matched post-BCG recurrences, we delineate three distinct BCG Response Subtypes (BRS1-3). Patients with BRS3 tumors have reduced high-grade recurrence-free and progression-free survival as compared to BRS1/2 tumors. BRS3 tumors expressed high EMT-basal markers and have an immunosuppresive profile, which was confirmed with spatial proteomics. The prognostic relevance of the BRS stratification was validated and improved the current recommended clinical guidelines, which is based on clinicopathological variables.

## Contents

The `BRSpred` R-package (see **citation**) contains a nearest-centroid classifier that can be used for transcriptomic data of primary HR-NMIBC tumors that were treated with BCG. The BCG response subtyper (BRS) identifies 3 distinct molecular classes: BRS1, BRS2 and BRS3 tumors. In addition, the package contains a vignette reproducing the key results from the main manuscript, and comes with two datasets: 

* Cohort A (n=132): Discovery Cohort
* Cohort B (n=151): Validation Cohort

## Installation

In order to install ```BRSpred``` directly from GitHub:

```
devtools::install_github("CostelloLab/BRSpred", dependencies=TRUE)
library(BRSpred)
```

Alternatively, one can install the package from a local source tarball in the R terminal with:
```
install.packages("BRSpred_VERSION.tar.gz", repos=NULL, type="source")
```

Note that one must install dependencies from a R package repository such as CRAN in order to BRSpred to work.

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

## Usage

The main wrapper for the pre-fitted BRS predictor object is the ```BRS```-function; for example predicting the training cohort:
```
BRSpred::BRS(BRSpred::CohortA.vst)
```

For further usage and prior results, see the vignette ```BRS prediction```.

## Citation

de Jong et al. *Submitted*

## Reproducing results from de Jong et al.

See the vignette ```BRS prediction``` listed for the package at ```vignette(package="BRSpred")``` for training and validation using the proposed BRS-methodology. The raw contents of this vignette are located in the folder ```/vignettes/brs.Rmd```, or for the compiled package inside ```/doc/brs.html``` (if compiled to HTML as expected by default).

## Acknowledgements

The BCG Response Subtyper was developed in a collaboration between the Erasmus MC Cancer Institute, the Anschutz Medical Campus, Cedars-Sinai Samuel Oschin Comprehensive Cancer Institute and University of Turku & the Finnish Cancer Institute.
