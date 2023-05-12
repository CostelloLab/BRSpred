###
#
# Exported .RData objects
#
###

#' Clinical information and gene expression matrices for the Erasmus Medical Center bladder cancer cohort A and B
#' 
#' @name ErasmusMC
#' @rdname ErasmusMC
NULL

#' erasmus_clinical: Clinical information for the Erasmus MC cohorts A (pre- and post-BCG) and B
#'
#' @rdname ErasmusMC
"erasmus_clinical"

#' erasmus_immune: Immune cell estimates (Sup Table S5) for the Erasmus MC cohorts A and B
#'
#' @rdname ErasmusMC
"erasmus_immune"

#' CohortA_pre: ErasmusMC NMIBC Cohort A (training pre-BCG, variance stabilizing transformation)
#'
#' @rdname ErasmusMC
"CohortA_pre"

#' CohortA_post: ErasmusMC NMIBC Cohort A (post-BCG, variance stabilizing transformation)
#'
#' @rdname ErasmusMC
"CohortA_post"

#' CohortB: ErasmusMC NMIBC Cohort B (validation, variance stabilizing transformation)
#'
#' @rdname ErasmusMC
"CohortB"

#' Gene expression matrices of raw aligned read counts in the Erasmus Medical Center bladder cancer cohort A and B
#' 
#' @name CohortAB_counts
#' @rdname CohortAB_counts
NULL

#' CohortA and B aligned raw counts data from RNA-seq
#'
#' @rdname CohortAB_counts
"CohortA_counts"

#' CohortA and B aligned raw counts data from RNA-seq
#'
#' @rdname CohortAB_counts
"CohortB_counts"

## Package documentation

#' BRSpred: BCG Response Subtype Predictor for High-risk Non-Muscle Invasive Bladder Cancer
#'
#' Consensus clustering based molecular predictor was first constructed for characterizing BCG responsiveness in bladder cancer. With these initial classes, a pamr-classifier is constructed for future data. Further, transcriptomics datasets for training and validation are provided along with key clinical variables. See the main BRS vignette for further details.
#' 
#' @import survminer
#' @import ConsensusClusterPlus
#'
#' @docType package
#' @name BRSpred
NULL
