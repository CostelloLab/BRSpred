###
#
# Molecular characterization of non-muscle invasive bladder cancer (pamr transcriptomics predictor)
# de Jong, Laajala et al.
#
###

#' Main BRS subclass predictor
#'
#' @details
#' This is the main BCG response subtype predictor function. It essentially works as a wrapper for the original 3-class pamr-object, with imputation and other convenience functionality.
#'
#' @param x Input data matrix
#' @param train Training data matrix; used for imputing missing gene names (by default original training data matrix)
#' @param pamrobj pamr-package object of class 'pamr', by default the original one trained for the package
#' @param threshold pamr threshold parameter (by default the optimal threshold found by CV in training data)
#' @param type Type of prediction as given to pamr; eligible values 'class', 'posterior', 'centroid', 'nonzero'
#' @param qnormalize Should quantile normalization be applied to 'x' in respect to the training data
#' @param zscale Should z-score scaling be applied to data; by default TRUE
#' @param verb Verbosity
#' @param ... Additional parameters passed on to pamr::pamr.predict
#'
#' @references
#' de Jong F. C., Laajala T. D., et al. Citation
#'
#' @return pamr::pamr.predict-call predictions for the input 'x'
#'
#' @rdname BRS
#'
#' @examples
#' library(BRSpred)
#' BRS(x=CohortB.vst, train=CohortA.vst, threshold=0.6606584)
#'
#' @import pamr
#'
#' @importFrom preprocessCore normalize.quantiles.use.target
#'
#' @export
BRS <- function(
	# Input data matrix
	x,
	# Training data matrix; used for imputing missing gene names (by default original training data matrix)
	train = BRSpred::CohortA.vst,
	# pamr-object used for prediction (trained previously using 'train')
	pamrobj = BRSpred::pamr3cl,
	# pamr threshold parameter (by default the optimal threshold found by CV in training data)
	threshold = 0.1651646,
	# Prediction type; allowed by pamr: 'class', 'posterior', 'centroid', or 'nonzero'
	type = "class",
	# Should quantile normalization be applied to 'x' in respect to the training data
	qnormalize = TRUE,
	# Should z-score scaling be applied to data; by default TRUE
	zscale = TRUE,
	# Verbosity
	verb = TRUE,
	# Additional parameters passed on to pamr.predict
	...
){
	if(qnormalize){
		if(verb) print("Performing quantile normalization")
		dimn <- dimnames(x)
		x <- preprocessCore::normalize.quantiles.use.target(
			x = as.matrix(x), 
			target = c(unlist(BRSpred::CohortA.vst))
		)
		dimnames(x) <- dimn
		if(verb){
			print("Quantiles in x and train post-normalization:")
			print(quantile(as.matrix(x)))
			print(quantile(c(unlist(train))))
		}
	}
	if(zscale){
		if(verb) print("Z-score scaling to 'x' and 'train' (shift to zero mean, unit stdev)")
		x <- t(scale(t(x)))
		train <- t(scale(t(train)))
	}
	if(verb){
		print("input 'x' dimensions prior to BRS_impute:")
		print(dim(x))
	}
	x <- BRSpred:::BRS_impute(x=x, train=train)
	if(verb){
		print("input 'x' dimensions post BRS_impute:")
		print(dim(x))
	}
	x <- x[rownames(pamrobj$centroids),]
	if(verb){
		print("input 'x' dimensions after subsetting with pamr-object centroids:")
		print(dim(x))
	}	
	pamr::pamr.predict(fit=pamrobj, newx=x, threshold=threshold, type=type, ...)
}

#' Imputation based on medians from training data for gene names that are missing from input data
#'
#' @param x New data matrix 'x' with potentially missing rownames
#' @param train Original training data matrix
#' @param verb Verbosity
#'
#' @details
#' Imputation of gene expression values based on median from an existing training data.
#'
#' @rdname BRS
BRS_impute <- function(
	x,
	train,
	verb = TRUE
){
	# Identify row names not present in x
	gs <- setdiff(rownames(train), rownames(x))
	if(!length(gs)==0){
		warning(paste0("Missing gene names imputed: ", paste0(gs, collapse=", ")))
		imputed <- do.call("rbind", lapply(gs, FUN=function(g){
			rep(median(as.numeric(train[g,])), times=ncol(x))
		}))
		rownames(imputed) <- gs
		colnames(imputed) <- colnames(x)
		x <- rbind(x, imputed)
	}
	# Return x in same order as train for rownames
	x[rownames(train),]	
}

#' List top genes from a thresholded pamr-object without forcefully printing it to R stdout
#'
#' @param fit A fitted pamr-object by pamr.train
#' @param data The original data matrix provided for training
#' @param threshold The penalization threshold parameter
#' @param fitcv A pamr cross-validation object
#' @param genenames Inclusion of gene names taken from parameter 'data'
#'
#' @details
#' Original implementation by authors of the 'pamr'-package's function 'pamr.listgenes'; this is merely a small adjustment for user convenience.
#'
#' @export
BRS_pamr.listgenes <- function(
	fit, 
	data, 
	threshold, 
	fitcv = NULL, 
	genenames = FALSE
) {
	x <- data$x[fit$gene.subset, fit$sample.subset]
	if (genenames) {
		gnames <- data$genenames[fit$gene.subset]
	}
	if (!genenames) {
		gnames <- NULL
	}
	geneid <- data$geneid[fit$gene.subset]
	if (!is.null(fit$y)) {
		nc <- length(fit$y)
	}
	if (is.null(fit$y)) {
		nc <- ncol(fit$proby)
	}
	clabs <- colnames(fit$centroids)
	aa <- pamr.predict(fit, x, threshold = threshold, type = "nonzero")
	cen <- pamr.predict(fit, x, threshold = threshold, type = "centroid")
	d <- (cen - fit$centroid.overall)[aa, , drop = FALSE]/fit$sd[aa]
	gene.order <- order(-apply(abs(d), 1, max))
	d <- round(d, 4)
	g <- gnames[aa]
	g1 <- geneid[aa]
	if (is.null(gnames)) {
		gnhdr <- NULL
	}
	if (!is.null(gnames)) {
		gnhdr <- "name"
	}
	if (!is.null(fitcv)) {
		nfold = length(fitcv$cv.objects)
		ind = matrix(F, nrow = nrow(x), ncol = nfold)
		ranks = NULL
		for (ii in 1:nfold) {
			cen = pamr.predict(fitcv$cv.objects[[ii]], x[, -fitcv$folds[[ii]]], 
			threshold = 0, type = "centroid")
			dtemp <- (cen - fitcv$cv.objects[[ii]]$centroid.overall)[, 
			, drop = FALSE]/fitcv$cv.objects[[ii]]$sd
			r <- apply(abs(dtemp), 1, max)
			ranks = cbind(ranks, rank(-abs(r)))
			junk = pamr.predict(fitcv$cv.objects[[ii]], x[, -fitcv$folds[[ii]]], 
			threshold = threshold, type = "nonzero")
			ind[junk, ii] = T
		}
		av.rank = apply(ranks, 1, mean)
		av.rank = round(av.rank[aa], 2)
		prop = apply(ind[aa, , drop = F], 1, sum)/nfold
		}
	# options(width = 500) # TDL: No need for adjusting text output
	schdr <- paste(clabs, "score", sep = "-")
	if (is.null(fitcv)) {
		res <- cbind(as.character(g1), g, d)[gene.order, , drop = F]
		dimnames(res) <- list(NULL, c("id", gnhdr, schdr))
	}
	if (!is.null(fitcv)) {
		res <- cbind(as.character(g1), g, d, av.rank, prop)[gene.order, 
		, drop = F]
		dimnames(res) <- list(NULL, c("id", gnhdr, schdr, 
		"av-rank-in-CV", "prop-selected-in-CV"))
	}
	#print(res, quote = FALSE)
	res # TDL: Instead of printing, return this original result matrix
}
