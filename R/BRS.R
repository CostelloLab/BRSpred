###
#
# Molecular characterization of non-muscle invasive bladder cancer (pamr transcriptomics predictor)
# de Jong, Laajala et al.
#
###


#' Main BRS subclass predictor
#'
#' @details
#' This is the main BCG response subtype predictor function. It essentially works as a wrapper for 3-class pamr-object classifier, with training and optimal threshold determined using cross-validation. Convenience functions, such as imputation for missing genes and quantile normalization, are provided; note however
#'
#' @param x Input data matrix
#' @param train Training data matrix; used for imputing missing gene names (by default original training data matrix)
#' @param pamrobj pamr-package object of class 'pamr', by default the original one trained for the package
#' @param threshold pamr-prediction threshold parameter (if missing, by default the optimal threshold is identified by minimizing CV misclassification rate)
#' @param Type of prediction as given by pamr; eligible values 'class', 'posterior', 'centroid', 'nonzero' (see ?pamr.predict)
#' @param qnormalize Should quantile normalization be applied to 'x' in respect to the training data
#' @param zscale Should z-score scaling be applied to data; by default TRUE
#' @param verb Verbosity
#' @param ... Additional parameters passed on to pamr::pamr.predict
#'
#' @references
#' de Jong F. C., Laajala T. D., et al. Citation
#'
#' @return pamr::pamr.predict-call predictions or a list with the prediction and corresponding data matrices and pamr-object and threshold
#'
#' @rdname BRS
#'
#' @examples
#' library(BRSpred)
#' predict_post <- BRSpred::BRS(newx = BRSpred::CohortA_post, scale = "together")
#' predict_cohortb <- BRSpred::BRS(newx = BRSpred::CohortB, scale = "independent")
#'
#' @import pamr
#' @import matrixStats
#' @import survminer
#'
#' @importFrom preprocessCore normalize.quantiles.use.target
#' @importFrom stats median quantile
#' @importFrom matrixStats rowVars
#'
#' @export
BRS <- function(
	# New data matrix to predict on
	newx,
	# Training data matrix x and classes y
	trainx = as.matrix(BRSpred::CohortA_pre),
	trainy = BRSpred::erasmus_clinical[colnames(BRSpred::CohortA_pre),"Erasmus.BRS"],
	# Should only common genes be considered (prevents need for potential imputation)
	common = TRUE,
	# Subset of genes used in the training; by default top 2000 variable genes in the training data
	genes,
	# pamr-object used for prediction (trained previously using 'train')
	pamrobj,
	# RNG seed, if cross-validation is used for threshold-determination this should be set
	seed = 1234,
	# If provided, will designate number of folds in the CV
	nfold,
	# pamr-prediction threshold parameter (if missing, by default the optimal threshold is identified by minimizing CV misclassification rate)
	threshold,
	# Type of prediction as given by pamr; eligible values 'class', 'posterior', 'centroid', 'nonzero' (see ?pamr.predict)
	type = "class",
	# Scaling of data; none, or together with or independently of the training data
	scale = c("together", "independent", "none"),
	# Should quantile normalization be applied to 'newx' in respect to the training data; TRUE does this prior to pamr.training by default
	qnormalize = TRUE,
	# Should gene imputation allowed via internal function if not pre-processed by user
	impute = FALSE,
	# Should all objects be returned - will instead create a list with predictions, pamr object, newx, trainx, and trainy
	getall = FALSE,
	# Verbosity
	verb = TRUE,
	# Additional parameters passed on to pamr.predict
	...
){
	# If seed is needed; set to NULL if ought to be omitted - essential for CV reproducibility
	if(!is.null(seed)){
		set.seed(seed)
	}
	
	# Should only common genes be considered for training
	if(common){
		newx <- newx[intersect(rownames(newx), rownames(trainx)),]
		trainx <- trainx[rownames(newx),]
	}
	
	# Perform pre-training quantile-normalization	
	if(qnormalize){
		if(verb) cat("\nPerforming pre-pamr quantile normalization\n")
		dimn <- dimnames(newx)
		newx <- preprocessCore::normalize.quantiles.use.target(
			x = as.matrix(newx), 
			target = c(unlist(trainx))
		)
		dimnames(newx) <- dimn
	}

	# Set of genes to use in the pamr object
	# By default 2000 top varying genes
	if(missing(genes)){
		genes <- head(order(matrixStats::rowVars(trainx), decreasing = TRUE), 2000)
	}
	# Subset to genes of interest
	if(missing(pamrobj)){
		trainx <- trainx[genes,]
	}
	newx <- newx[genes,]
	
	# z-scale gene-wise; can be done together or separately from the training data
	if(scale %in% c("together", 1)){
		if(verb) cat("\nRow-wise scaling together with the training data\n")
		# Scale together
		sets <- cbind(newx, trainx)
		sets <- t(scale(t(sets)))
		# Subset back to original
		newx <- sets[,1:ncol(newx)]
		trainx <- sets[,(ncol(newx)+1):ncol(sets)]		
	}else if(scale %in% c("independent", "separately", 2)){
		if(verb) cat("\nScaling newx independently of training data\n")
		newx <- t(scale(t(newx)))
		trainx <- t(scale(t(trainx)))
	}else{
		if(verb) cat("\nNo scaling performed\n")
	}

	# Format a training list for a pamr-object
	train.list <- list(
		x = trainx, # Input data matrix
		y = trainy, # Response vector to train on
		genenames = rownames(trainx), # Names for genes
		sampnames = colnames(trainy), # Sample names
		gene.ids = rownames(trainy) # Identifiers for genes
	)

	# pamrobj missing - train one
	if(missing(pamrobj)){
		if(verb) cat("\npamr-object missing, training one using pamr::pamr.train\n")
		
		pamrobj <- pamr::pamr.train(data = train.list)
	}
	# If threshold is missing, obtain it from cross-validation on the pamr object
	if(missing(threshold)){
		if(verb) cat("\nPre-determined pamr-object threshold not provided, calculating one via CV\n")

		# Run cross-validation
		if(missing(nfold)){
			cv <- pamr::pamr.cv(pamrobj, data = train.list)
		}else{
			cv <- pamr::pamr.cv(pamrobj, data = train.list, nfold = nfold)
		}
		
		# Pick the largest threshold that minimizes the misclassification rate
		threshold <- cv$threshold[max(which(cv$error == min(cv$error)))]
	}

	# If imputation is allowed for missing genes
	if(impute){
		newx <- BRSpred:::BRS_impute(x=newx, train=trainx)
		newx <- newx[rownames(pamrobj$centroids),]
	}


	# Subsetting new data to the training data gene-centroids prior to predictions
	newx <- newx[rownames(pamrobj$centroids),]	

	# Predict using pamr
	pred <- pamr::pamr.predict(fit=pamrobj, newx=newx, threshold=threshold, type=type, ...)
	# Return just predictions
	if(!getall){	
		# Predict classes
		pred
	# Additional output as a list, such as the data matrices and pamr-object
	}else{
		list(
			pred = pred,
			newx = newx,
			trainx = trainx,
			trainy = trainy,
			pamrobj = pamrobj,
			threshold = threshold
		)
	}
}

#' Imputation based on medians from training data for gene names that are missing from input data
#'
#' @param x New data matrix 'x' with potentially missing rownames
#' @param train Original training data matrix
#' @param FUN Function that provides within-row imputations; by default mean, see also e.g. median
#' @param verb Verbosity
#'
#' @details
#' Imputation of gene expression values based on median from an existing training data.
#'
#' @rdname BRS
BRS_impute <- function(
	x,
	train,
	FUN = \(x) { median(x, na.rm=TRUE) },
	verb = TRUE
){
	# Identify row names not present in x
	gs <- setdiff(rownames(train), rownames(x))
	if(!length(gs)==0){
		warning(paste0("Missing gene names imputed: ", paste0(gs, collapse=", ")))
		imputed <- do.call("rbind", lapply(gs, FUN=function(g){
			rep(FUN(as.numeric(train[g,])), times=ncol(x))
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
