#####################################################################
# Package: INSiGHT
# Version: 1.0.0
# Date : 2025-11-14
######################################################################

##########################################################
#                   INSiGHT functions                    #
##########################################################
#' @title An S4 class to store key information.
#'
#' @slot geneExpr A list containing gene expression matrices of all the samples to analyze.
#' The list is of length \eqn{M}, where \eqn{M} is the number of samples.
#' Each gene expression matrix is of dimension \eqn{G \times N}, where \eqn{G} rows represent \eqn{G} genes and
#' \eqn{N} columns represent \eqn{N} spots across this sample.
#' @slot spaCoord A list containing spatial coordinate matrices of all the samples to analyze.
#' The list is of length \eqn{M}, where \eqn{M} is the number of samples.
#' Each spatial coordinate matrix is of dimension \eqn{N \times 2},
#' where \eqn{N} rows represent \eqn{N} spots across this sample,
#' and 2 columns represents x axis and y axis.
#' @slot covariateMat A list containing covariate matrices of all the samples to analyze.
#' The list is of length \eqn{M}, where \eqn{M} is the number of samples.
#' Each covariate matrix is of dimension \eqn{N \times p},
#' where \eqn{N} rows represent \eqn{N} spots across this sample,
#' and \eqn{p} columns represents \eqn{p} covariates.
#' @slot weightMat A list containing weight matrices of all the samples to analyze.
#' The list is of length \eqn{M}, where \eqn{M} is the number of samples.
#' Each covariate matrix is of dimension \eqn{N \times N},
#' where \eqn{N} is the number of spots across this sample.
#' @slot samples A character string vector containing the sample IDs of all the samples to analyze.
#' The vector is of length \eqn{M}, where \eqn{M} is the number of samples.
#' @slot genes A character string vector containing the gene IDs shared acorss all the samples to analyze.
#' The vector is of length \eqn{G}, where \eqn{G} is the number of genes to analyze.
#' @slot barcodes A list containing barcode (spot IDs) vectors of all the samples to analyze.
#' The list is of length \eqn{M}, where \eqn{M} is the number of samples.
#' Each barcode vector contains all the barcodes (spot IDs) of a sample.
#' @slot bandwidths A numeric vector containing bandwidth hyper-parameters for the kernel functions
#' to contruct spatial similarity matrices.
#' @slot W.info A list containing the information of the weighted spatial similarity matrices,
#' including W matrices, the corresponding eigen-values, and etc.
#' @slot rankScores A list containing rank score matrices of all the samples to analyze.
#' The list is of length \eqn{M}, where \eqn{M} is the number of samples.
#' Each rank score matrix is of dimension \eqn{G \times N}, where \eqn{G} rows represent \eqn{G} genes and
#' \eqn{N} columns represent \eqn{N} spots across this sample.
#' @slot testResults A data frame containing the test results of all the genes analyzed.
#' @slot features Contains other features
#' @import methods
#'
#' @export

setClass("INSiGHT", slots = list(
  geneExpr = "ANY",
  spaCoord = "ANY",
  covariateMat = "ANY",
  weightMat = "ANY",
  samples = "character",
  genes = "character",
  barcodes = "ANY",
  bandwidths = "numeric",
  W.info = "ANY",
  rankScores = "ANY",
  testResults = "ANY",
  features = "ANY"
))

# ==============================================================================
# ==============================================================================
check_list_mats <- function(x, nm,
                            allow_sparse     = TRUE,
                            allow_dataframe  = TRUE,
                            df_require_numeric = FALSE) {
  if (!is.list(x)) {
    stop(sprintf('"%s" must be a list.', nm), call. = FALSE)
  }

  is_mat_like <- function(obj) {
    if (is.matrix(obj)) return(TRUE)
    if (allow_sparse && inherits(obj, "Matrix")) return(TRUE)
    if (allow_dataframe && is.data.frame(obj)) {
      if (!df_require_numeric) return(TRUE)
      # require all columns numeric if requested
      return(all(vapply(obj, is.numeric, logical(1L))))
    }
    FALSE
  }

  ok <- vapply(x, is_mat_like, logical(1L))
  if (!all(ok)) {
    bad_idx <- which(!ok)
    bad_cls <- vapply(x[bad_idx], function(obj) paste(class(obj), collapse = "/"), character(1L))
    req <- paste0(
      "base matrix",
      if (allow_sparse) ", Matrix::Matrix" else "",
      if (allow_dataframe) if (df_require_numeric) ", numeric data.frame" else ", data.frame" else ""
    )
    stop(sprintf(
      'All elements of "%s" must be one of: %s. Problem at indices: %s (classes: %s)',
      nm, req, paste(bad_idx, collapse = ", "), paste(bad_cls, collapse = "; ")
    ), call. = FALSE)
  }

  invisible(TRUE)
}

# ==============================================================================
# ==============================================================================
#' @title Create the INSiGHT object.
#' @param data.input A list containing original gene expression matrices
#' of all the samples to analyze.
#' The list is of length \eqn{M}, where \eqn{M} is the number of samples.
#' Each gene expression matrix is of dimension \eqn{G \times N}, where \eqn{G} rows represent \eqn{G} genes and
#' \eqn{N} columns represent \eqn{N} spots across this sample.
#' @param spatial.locs A list containing spatial coordinate matrices of all the samples to analyze.
#' The list is of length \eqn{M}, where \eqn{M} is the number of samples.
#' Each spatial coordinate matrix is of dimension \eqn{N \times 2},
#' where \eqn{N} rows represent \eqn{N} spots across this sample,
#' and 2 columns represents x axis and y axis.
#' @param sample.names A character string vector specifying the sample IDs of all the samples to analyze.
#' The vector is of length \eqn{M}, where \eqn{M} is the number of samples.
#' @param covariates (default NULL) A list containing covariate matrices of all the samples to analyze.
#' The list is of length \eqn{M}, where \eqn{M} is the number of samples.
#' Each covariate matrix is of dimension \eqn{N \times p},
#' where \eqn{N} rows represent \eqn{N} spots across this sample,
#' and \eqn{p} columns represents \eqn{p} covariates.
#' @param weights.list A list containing weight matrices of all the samples to analyze.
#' The list is of length \eqn{M}, where \eqn{M} is the number of samples.
#' Each covariate matrix is of dimension \eqn{N \times N},
#' where \eqn{N} is the number of spots across this sample.
#' @return A INSiGHT object.
#'
#' @import methods
#'
#' @export
#'

createINSiGHTobject <- function(data.input, spatial.locs,
                                sample.names = NULL,
                                covariates = NULL,
                                weights.list = NULL){
  # 1) Every element in each list must be a matrix or data.frame
  check_list_mats(data.input, "data.input",
                  allow_sparse = TRUE,
                  allow_dataframe  = FALSE,
                  df_require_numeric = FALSE)
  check_list_mats(spatial.locs, "spatial.locs",
                  allow_sparse = TRUE,
                  allow_dataframe  = TRUE,
                  df_require_numeric = TRUE)
  if (!is.null(covariates)) {
    check_list_mats(covariates, "covariates",
                    allow_sparse = TRUE,
                    allow_dataframe  = FALSE,
                    df_require_numeric = FALSE)
  }
  if (!is.null(weights.list)) {
    check_list_mats(weights.list, "weights.list",
                    allow_sparse = TRUE,
                    allow_dataframe  = FALSE,
                    df_require_numeric = FALSE)
  }
  # ----------------------------------------------------------------------------
  # 2) sample.names must be a (1D) character vector, if provided
  if (!is.null(sample.names)) {
    if (!is.character(sample.names) || is.matrix(sample.names) || is.array(sample.names)) {
      stop('"sample.names" must be a character vector.', call. = FALSE)
    }
    if (anyNA(sample.names)) {
      stop('"sample.names" contains NA values.', call. = FALSE)
    }
  }
  # ----------------------------------------------------------------------------
  # 3) Lengths of the three lists and sample.names must match
  n.Sample <- length(data.input)

  if (length(spatial.locs) != n.Sample) {
    stop(sprintf(
      "Lengths must match: length(data.input) = %d, length(spatial.locs) = %d",
      n.Sample, length(spatial.locs)
    ), call. = FALSE)
  }

  if (!is.null(covariates) && length(covariates) != n.Sample) {
    stop(sprintf(
      "Lengths must match: length(data.input) = %d, length(covariates) = %d",
      n.Sample, length(covariates)
    ), call. = FALSE)
  }

  if (!is.null(weights.list) && length(weights.list) != n.Sample) {
    stop(sprintf(
      "Lengths must match: length(data.input) = %d, length(weights.list) = %d",
      n.Sample, length(weights.list)
    ), call. = FALSE)
  }

  if (!is.null(sample.names) && length(sample.names) != n.Sample) {
    stop(sprintf(
      'Length(sample.names) must equal %d; got %d',
      n.Sample, length(sample.names)
    ), call. = FALSE)
  }
  # ----------------------------------------------------------------------------
  # Provide default sample.names if missing
  if (is.null(sample.names)) {
    nm <- names(data.input)
    if (is.null(nm) || anyNA(nm) || any(nm == "")) {
      sample.names <- paste0("sample", seq_len(n.Sample))
    } else {
      sample.names <- nm
    }
  }
  # ----------------------------------------------------------------------------
  barcodes.list <- vector("list", n.Sample)
  has_cov <- !is.null(covariates)
  has_weights <- !is.null(weights.list)

  for (iSample in seq_len(n.Sample)) {
    # Helpful sample label (if you have sample.names defined)
    sname <- if (exists("sample.names") && length(sample.names) >= iSample) {
      sample.names[iSample]
    } else {
      paste0("Sample", iSample)
    }

    Xi <- data.input[[iSample]]        # genes x spots
    Si <- spatial.locs[[iSample]]      # spots x 2 (or more)
    Ci <- if (has_cov) covariates[[iSample]] else NULL
    Wi <- if (has_weights) weights.list[[iSample]] else NULL


    ## Check 1: dimension compatibility
    if (ncol(Xi) != nrow(Si)) {
      stop(sprintf("Sample %s: ncol(data.input)=%d must equal nrow(spatial.locs)=%d.",
                   sname, ncol(Xi), nrow(Si)), call. = FALSE)
    }
    if (has_cov && ncol(Xi) != nrow(Ci)) {
      stop(sprintf("Sample %s: ncol(data.input)=%d must equal nrow(covariates)=%d.",
                   sname, ncol(Xi), nrow(Ci)), call. = FALSE)
    }
    if (has_weights && nrow(Wi) != ncol(Wi)) {
      stop(sprintf("Sample %s: nrow(weights.list)=%d must equal ncol(weights.list)=%d.",
                   sname, nrow(Wi), ncol(Wi)), call. = FALSE)
    }
    if (has_weights && ncol(Xi) != nrow(Wi)) {
      stop(sprintf("Sample %s: ncol(data.input)=%d must equal nrow(weights.list)=%d.",
                   sname, ncol(Xi), nrow(Wi)), call. = FALSE)
    }

    ## Check 2: required row/column names exist
    rn_Si <- rownames(Si)
    if (is.null(rn_Si)) {
      stop(sprintf("Sample %s: rownames(spatial.locs[[%d]]) cannot be NULL.", sname, iSample), call. = FALSE)
    }
    cn_Xi <- colnames(Xi)
    if (is.null(cn_Xi)) {
      stop(sprintf("Sample %s: colnames(data.input[[%d]]) cannot be NULL.", sname, iSample), call. = FALSE)
    }
    if (has_cov && is.null(rownames(Ci))) {
      stop(sprintf("Sample %s: rownames(covariates[[%d]]) cannot be NULL.", sname, iSample), call. = FALSE)
    }

    ## (Optional but recommended) ensure barcodes are unique within a sample
    # if (anyDuplicated(rn_Si)) {
    #   stop(sprintf("Sample %s: rownames(spatial.locs) contain duplicates.", sname), call. = FALSE)
    # }

    ## Check 3: exact matching (same order & values)
    if (!identical(cn_Xi, rn_Si)) {
      stop(sprintf("Sample %s: colnames(data.input) must exactly match rownames(spatial.locs).", sname), call. = FALSE)
    }
    if (has_cov && !identical(rownames(Ci), rn_Si)) {
      stop(sprintf("Sample %s: rownames(covariates) must exactly match rownames(spatial.locs).", sname), call. = FALSE)
    }

    ## Store barcodes and strip names to avoid downstream mismatches
    barcodes.list[[iSample]] <- rn_Si
    rownames(spatial.locs[[iSample]]) <- NULL
    colnames(data.input[[iSample]])   <- NULL
    if (has_cov) rownames(covariates[[iSample]]) <- NULL
  }
  # ----------------------------------------------------------------------------
  # Filter genes shared by all the samples
  if (is.null(rownames(data.input[[1]]))) {
    stop("Row names in \'data.input\' cannot be NULL!")
  } else {
    genes.used <- rownames(data.input[[1]])
  }
  if(n.Sample > 1) {
    for (iSample in 2:n.Sample) {
      if (is.null(rownames(data.input[[iSample]]))) {
        stop("Row names in \'data.input\' cannot be NULL!")
      } else {
        genes.used <- intersect(genes.used, rownames(data.input[[iSample]]))
      }
    }
  }
  data.use <- list()
  for (iSample in 1:n.Sample) {
    data.use[[iSample]] <- data.input[[iSample]][genes.used,]
    rownames(data.use[[iSample]]) <- NULL
  }
  rm(data.input)
  cat(paste0('##\t ',
             length(genes.used), " genes shared by ",
             n.Sample, " sample(s) will be examined. \n"))
  W.info <- list()
  W.info[[1]] <- list()
  W.info[[2]] <- list()
  names(W.info) <- c("W.mat", "W.eigenvalues")

  features <- list()

  object <- methods::new(
    Class = "INSiGHT",
    geneExpr = data.use,
    spaCoord = spatial.locs,
    covariateMat = covariates,
    weightMat = weights.list,
    samples =sample.names,
    genes = genes.used,
    barcodes = barcodes.list,
    W.info = W.info
  )
  return(object)
}

# ==============================================================================
# ==============================================================================
#' @title Do quality control for a INSiGHT object
#' Remove spots with low total expression count, low-expressed genes and mt genes.
#'
#' @param object INSiGHT object
#' @param spot.threshold (default 10) filter out the spots whose total expression count
#' lower than the threshold.
#' @param gene.threshold (default 0.05) filter out low-expressed genes
#' who has a percentage of non-zero expression spots lower than the threshold.
#' @return A INSiGHT object.
#'
#' @import methods
#'
#' @export
#'

QualityControl <- function(object, spot.threshold = 10, gene.threshold = 0.05) {
  # spot.threshold = 10; gene.threshold = 0.1

  n.Sample <- length(object@geneExpr)
  has_cov <- !is.null(object@covariateMat)
  has_weights <- !is.null(object@weightMat)

  geneRm.index <- c()
  # Identify mt genes
  if (length(grep("mt-", object@genes, ignore.case = T)) > 0) {
    mt_gene_list <- grep("mt-", object@genes, ignore.case = T)
    geneRm.index <- c(geneRm.index, mt_gene_list)
  }

  for(iSample in 1:n.Sample){
    # iSample = 1
    ## Spot QC
    spots.use.idx <- which(colSums(object@geneExpr[[iSample]]) >= spot.threshold)
    spots.use <- object@barcodes[[iSample]][spots.use.idx]
    pos.use <- object@spaCoord[[iSample]][spots.use.idx, ]
    numSpots.removed <- nrow(object@spaCoord[[iSample]]) - nrow(pos.use)
    object@spaCoord[[iSample]] <- pos.use
    object@barcodes[[iSample]] <- spots.use
    counts.use <- object@geneExpr[[iSample]][,spots.use.idx]
    object@geneExpr[[iSample]] <- counts.use
    if (has_cov) {
      covariates.use <- object@covariateMat[[iSample]][spots.use.idx,]
      object@covariateMat[[iSample]] <- covariates.use
      rm(covariates.use)
    }
    if (has_weights) {
      weights.use <- object@weightMat[[iSample]][spots.use.idx,spots.use.idx]
      object@weightMat[[iSample]] <- weights.use
      rm(weights.use)
    }

    cat(paste0('##\t ',
               numSpots.removed, ' spots with total gene expression less than ', spot.threshold, ' have been removed for the sample ',
               object@samples[iSample], '.\n'))

    rm(spots.use, spots.use.idx, pos.use, counts.use)

    ## Identify low expressed genes
    # The rate of spots with non-zero expression
    NonzeroExpressionRate <- rowSums(object@geneExpr[[iSample]] != 0) / ncol(object@geneExpr[[iSample]])
    gene.rm <- which(NonzeroExpressionRate < gene.threshold)
    geneRm.index <- c(geneRm.index, gene.rm)
  }

  geneRm.index <- sort(unique(geneRm.index), decreasing = FALSE)
  if(length(geneRm.index) != 0){
    genes.use <- object@genes[-geneRm.index]
    object@genes <- genes.use
    ## Gene quality control
    for (iSample in 1:n.Sample) {
      counts.use <- object@geneExpr[[iSample]][-geneRm.index,]
      object@geneExpr[[iSample]] <- counts.use
    }
  }
  cat(paste0('##\t ',
             length(geneRm.index), ' mt genes and genes with non-zero expression spot rate less than ', gene.threshold,
             ' have been removed for all the samples. \n'))
  cat(paste0('##\t ',
             length(object@genes), ' genes will be examined for all the samples after quality control. \n'))

  return(object)
}

# ==============================================================================
# ==============================================================================
#' @title Center a kernel matrix
#'
#' @param K A n*n symmetric matrix
#' @return A centered matrix
#'
#' @export
#'

CenterKernel <- function(K) {
  n <- nrow(K)
  rmean <- rowMeans(K)
  cmean <- colMeans(K)
  grand <- mean(K)
  K - outer(rmean, rep(1, n)) - outer(rep(1, n), cmean) + grand
}

# ==============================================================================
# ==============================================================================
HWH_matvec <- function(x, args) {
  x0 <- x - mean(x)                   # H %*% x
  y  <- as.vector(args$W %*% x0)      # W %*% (Hx)
  y - mean(y)                         # H %*% y
}

# ==============================================================================
# ==============================================================================
#' @title Compute the bandwidth hyper-parameters for a INSiGHT object
#' @description
#' Conpute the bandwidth hyperparameters of the Gaussian kernels for a INSiGHT,
#' using the quantiles of pair-wise distances of spots.
#' Only the pair-wise distances of the k-nearest neighbors for each spot will be calculated,
#' using a kd-tree.
#'
#' @param object INSiGHT object
#' @param k.NN The maximum number of nearest neighbours to compute for the kd-tree.
#' @param kernel.option Charater vector specifying whether a "single" kernel
#' or a "mixture" of multiple kernels will be contructed for spatial similarity.
#' @param quantile.bandwidth (default 0.1) A number of probability with values in \eqn{[0,1]},
#' for computing the quantile of pair-wise distances.
#' Works only when `kernel.option = "single"`.
#' When `kernel.option = "mixture"`, five quantiles (0.1, 0.3, 0.5, 0.7, 0.9) will computed.
#' @return A INSiGHT object.
#'
#' @import RANN
#' @keywords internal
#' @importFrom stats quantile
#'
#' @export
#'

ComputeBandwidth <- function(object, k.NN = 100,
                             kernel.option = c("single", "mixture"),
                             quantile.bandwidth = 0.1) {
  spatial.locs <- object@spaCoord
  n.Sample <- length(spatial.locs)
  kernel.arg <- match.arg(kernel.option)

  if(kernel.arg == "single") {
    for (iSample in 1:n.Sample) {
      n.locs <- nrow(spatial.locs[[iSample]])

      if (k.NN >= n.locs) {
        cat("##\t Warning: 'k.NN' is not less than the number of spots! \n")
        cat("##\t All pair-wise distances will be calculated. It may take a long time... \n")
        k.NN4Kernel <- n.locs
      }

      nn_result <- nn2(data = spatial.locs[[iSample]], k = k.NN + 1)
      if (iSample == 1) {
        Dval <- nn_result$nn.dists[,-1]
      } else {
        Dval <- rbind(Dval, nn_result$nn.dists[,-1])
      }
    }
    Dval.nz <- Dval[Dval>1e-8]

    bandwidths <- quantile(Dval.nz, probs = quantile.bandwidth)
  }else if(kernel.option == "mixture"){
    for (iSample in 1:n.Sample) {
      nn_result <- nn2(data = spatial.locs[[iSample]], k = k.NN + 1)
      if (iSample == 1) {
        Dval <- nn_result$nn.dists[,-1]
      } else {
        Dval <- rbind(Dval, nn_result$nn.dists[,-1])
      }
    }
    Dval.nz <- Dval[Dval>1e-8]

    bandwidths <- quantile(Dval.nz, probs = seq(0.1,0.9,by=0.2))
  }

  object@bandwidths <- bandwidths
  return(object)
}

# ==============================================================================
# ==============================================================================
#' @title Construct a spatial similarity matrix using Gaussian kernel for a INSiGHT object
#' @description
#' Construct a spatial similarity matrix for a INSiGHT object,
#' using a distance-based Gaussian kernel function
#' (i.e. \eqn{K(s_i, s_j) = \exp(-\frac{\lVert s_i - s_j \rVert^2}{2 \phi^2})}),
#' where \eqn{\phi} is the bandwidth parameter.
#'
#' @param object INSiGHT object.
#' @param phi Bandwidth hyper-parameter for the Gaussian kernel.
#' @param NN4Kernel (default FALSE) logical; If true, only the pair-wise distances
#' of k-nearest-neighbor will be computed for each spot, resulting in a sparse matrix.
#' If false, all the pair-wise distances will be computed.
#' @param k.NN4Kernel The maximum number of nearest neighbors to compute
#' for the construction of Gaussian kernel using a kd-tree. Works when `NN4Kernel == TRUE`.
#' @param eigen.less (default FALSE) logical; If true, only the top eigen-values
#' will be computed with `eigs_sym` function to save computing time. Works when `NN4Kernel == TRUE`.
#' If false, all the eigen-values will be computed, which may be time-consuming
#' when the number of spots is large.
#' @param eigen.prop A number between \eqn{[0,1]}
#' specifying the proportion of top eigen-values to compute.
#' Works when `eigen.less == TRUE`.
#' @param weighted (default FALSE)
#' If true, then the each spatial similarity matrix will be weighted
#' by its corresponding weight matrix.
#' @return A INSiGHT object.
#'
#' @import RANN
#' @import Matrix
#' @import matrixStats
#' @import RSpectra
#'
#' @export
#'

constructGaussianW <- function(object, phi,
                               NN4Kernel = FALSE, k.NN4Kernel = 100,
                               eigen.less = FALSE, eigen.prop = 0.3,
                               weighted = FALSE) {
  spatial.locs <- object@spaCoord
  n.Sample <- length(spatial.locs)
  has_cov <- !is.null(object@covariateMat)
  has_weights <- !is.null(object@weightMat)
  W.list <- list()
  W.ev.list <- list()

  if (isTRUE(NN4Kernel)) {
    for (iSample in 1:n.Sample) {
      n.locs <- nrow(spatial.locs[[iSample]])

      if (k.NN4Kernel >= n.locs) {
        cat("##\t Warning: 'k.NN4Kernel' is not less than the number of spots! \n")
        cat("##\t 'NN4Kernel' is not working! It may take a long time... \n")
        k.NN4Kernel <- n.locs
      }
      stopifnot((eigen.prop >= 0) & (eigen.prop <= 1))

      # k nearest neighbors (excluding self)
      nn_result <- nn2(data = spatial.locs[[iSample]], k = k.NN4Kernel + 1)
      idx  <- nn_result$nn.idx[, -1, drop = FALSE]
      dist <- nn_result$nn.dists[, -1, drop = FALSE]
      # Compute Gaussian kernels
      vals <- exp(-dist^2 / (2 * phi^2))
      # Build sparse similarity matrix directly
      i <- rep(1:n.locs, each = k.NN4Kernel)
      j <- as.vector(t(idx))
      x <- as.vector(t(vals))
      Wmat <- sparseMatrix(i = i, j = j, x = x, dims = c(n.locs, n.locs), giveCsparse = TRUE)

      if(has_weights && isTRUE(weighted)) {
        Cweight.sparse <- as(object@weightMat[[iSample]], "dgCMatrix")
        Wmat <- Wmat * Cweight.sparse
      }

      # Symmetrize: (W + t(W)) / 2
      Wmat <- symmpart(Wmat)

      # Threshold in-place and compress structure
      if (length(Wmat@x)) {
        idx <- which(abs(Wmat@x) < 1e-5)
        if (length(idx)) {
          Wmat@x[idx] <- 0
          Wmat <- drop0(Wmat)
        }
      }

      W.list[[iSample]] <- Wmat

      if (has_cov) {
        X <- object@covariateMat[[iSample]]
        XtX <- t(X) %*% X
        XtX_inv <- solve(XtX)
        P <- X %*% XtX_inv %*% t(X)
        H <- diag(1, nrow = nrow(X)) - P
        rm(X, XtX, XtX_inv, P)
      }

      if (eigen.less) {
        diag(Wmat) <- 1
        k <- floor(eigen.prop * n.locs)

        if (has_cov) {
          Wmat <- H %*% Wmat %*% H
          eig <- eigs_sym(Wmat, n = n.locs, k = k, which = "LA")
        } else {
          eig <- eigs_sym(A = HWH_matvec, n = n.locs, k = k, which = "LA",
                          args = list(W = Wmat))
        }
        W.ev <- eig$values
        W.ev <- W.ev[W.ev > 1e-10]
      } else {
        # Center W
        if (has_cov) {
          Wmat.center <- H %*% Wmat %*% H
        } else {
          Wmat.center <- CenterKernel(Wmat)
        }
        rm(Wmat) # zero out small stored values
        Wmat.center@x[Wmat.center@x < 1e-5] <- 0
        # remove explicit stored zeros
        Wmat.center <- drop0(Wmat.center)
        W.ev <- eigen(Wmat.center, only.values = TRUE, symmetric = TRUE)$values
      }
      W.ev.list[[iSample]] <- W.ev

      cat(paste0("##\t Sample ", object@samples[iSample], " is complete. \n"))
    }
  } else {
    for (iSample in 1:n.Sample) {
      n.locs <- nrow(spatial.locs[[iSample]])
      Dmat <- dist(spatial.locs[[iSample]])
      Wmat <- exp(-Dmat^2 / (2 * phi^2))
      Wmat.sparse <- Matrix(as.matrix(Wmat), sparse = TRUE)
      # zero out small stored values
      Wmat.sparse@x[Wmat.sparse@x < 1e-5] <- 0
      # remove explicit stored zeros
      Wmat.sparse <- drop0(Wmat.sparse)

      if(has_weights && isTRUE(weighted)) {
        Cweight.sparse <- as(object@weightMat[[iSample]], "dgCMatrix")
        Wmat.sparse <- Wmat.sparse * Cweight.sparse
      }

      W.list[[iSample]] <- Wmat.sparse

      if (has_cov) {
        X <- object@covariateMat[[iSample]]
        XtX <- t(X) %*% X
        XtX_inv <- solve(XtX)
        P <- X %*% XtX_inv %*% t(X)
        H <- diag(1, nrow = nrow(X)) - P
        Wmat.sparse.center <- H %*% Wmat.sparse %*% H
        rm(X, XtX, XtX_inv, P, H)
      } else {
        Wmat.sparse.center <- CenterKernel(Wmat.sparse)
      }

      rm(Wmat.sparse)
      # zero out small stored values
      Wmat.sparse.center@x[Wmat.sparse.center@x < 1e-5] <- 0
      # remove explicit stored zeros
      Wmat.sparse.center <- drop0(Wmat.sparse.center)

      # Eigen decomposition
      W.ev <- eigen(Wmat.sparse.center, only.values = TRUE, symmetric = TRUE)$values
      W.ev.list[[iSample]] <- W.ev
      cat(paste0("##\t Sample ", object@samples[iSample], " is complete. \n"))
    }
  }

  object@W.info$W.mat <- append(object@W.info$W.mat, list(W.list))
  object@W.info$W.eigenvalues <- append(object@W.info$W.eigenvalues, list(W.ev.list))
  return(object)
}

# ==============================================================================
# ==============================================================================
#' @title Run pipeline to process INSiGHT object
#' @description
#' Process the INSiGHT object before running hypothesis testing, including
#' quality control process, estimation of bandwidth parameter,
#' construction of similarity matrices, and assignment of rank scores.
#'
#' @param object INSiGHT object.
#' @param spot.threshold (default 10) filter out the spots whose total expression count
#' lower than the threshold.
#' @param gene.threshold (default 0.05) filter out low-expressed genes
#' who has a percentage of non-zero expression spots lower than the threshold.
#' @param k.NN4bandwidth (default 100) The maximum number of nearest neighbors
#' to estimate bandwidths.
#' @param quantile.bandwidth (default 0.1) A number of probability with values in \eqn{[0,1]},
#' for computing the quantile of pair-wise distances.
#' Works only when `kernel.option = "single"`.
#' When `kernel.option = "mixture"`, five quantiles (0.1, 0.3, 0.5, 0.7, 0.9) will computed.
#' @param NN4Kernel (default FALSE) logical; If true, only the pair-wise distances
#' of k-nearest-neighbor will be computed for each spot, resulting in a sparse matrix.
#' If false, all the pair-wise distances will be computed.
#' @param k.NN4Kernel The maximum number of nearest neighbors to compute
#' for the construction of Gaussian kernel using a kd-tree. Works when `NN4Kernel == TRUE`.
#' @param eigen.less (default FALSE) logical; If true, only the top eigen-values
#' will be computed with `eigs_sym` function to save computing time. Works when `NN4Kernel == TRUE`.
#' If false, all the eigen-values will be computed, which may be time-consuming
#' when the number of spots is large.
#' @param eigen.prop A number between \eqn{[0,1]}
#' specifying the proportion of top eigen-values to compute.
#' Works when `eigen.less == TRUE`.
#' @param kernel.option Charater vector specifying whether a "single" kernel
#' or a "mixture" of multiple kernels will be contructed for spatial similarity.
#' @param rank.tie A character string specifying how ties are treated.
#' Run `?rank` for more details.
#' @param rank.zero (default FALSE) logical.
#' If true, assign zeros a rank of 0 and rank nonzero values from 1 upward.
#' If false, rank all the values.
#' @param weighted (default FALSE)
#' If true, then the each spatial similarity matrix will be weighted
#' by its corresponding weight matrix.
#'
#' @return A INSiGHT object.
#'
#' @import RANN
#' @import Matrix
#' @import matrixStats
#' @import RSpectra
#'
#' @export
#'

processINSiGHT <- function(object,
                           spot.threshold = 10, gene.threshold = 0.05,
                           k.NN4bandwidth = 100, quantile.bandwidth = 0.1,
                           NN4Kernel = FALSE, k.NN4Kernel = 100,
                           eigen.less = FALSE, eigen.prop = 0.3,
                           kernel.option = c("single","mixture"),
                           rank.tie = c("min", "max", "first", "last", "random", "average"),
                           rank.zero = TRUE,
                           weighted = FALSE) {
  # Run quality control process
  object <- QualityControl(object = object, spot.threshold = spot.threshold, gene.threshold = gene.threshold)

  # Compute bandwidth parameters
  object <- ComputeBandwidth(object = object, k.NN = k.NN4bandwidth, kernel.option = kernel.option, quantile.bandwidth = quantile.bandwidth)

  if(length(object@W.info$W.mat) != 0 | length(object@W.info$W.eigenvalues) != 0) {
    W.info <- list()
    W.info[[1]] <- list()
    W.info[[2]] <- list()
    names(W.info) <- c("W.mat", "W.eigenvalues")
    object@W.info <- W.info
  }

  # ----------------------------------------------------------------------------
  kernel.arg <- match.arg(kernel.option)
  object@features[[1]] <- kernel.arg
  object@features[[2]] <- NN4Kernel
  object@features[[3]] <- eigen.less
  names(object@features) <- c("kernel.option", "NN4Kernel", "eigen.less")

  if(kernel.arg == "mixture"){
    for (iKernel in 1:length(object@bandwidths)) {
      cat(paste0("##\t Constructing Gaussian kernel ", iKernel, "... \n"))
      object <- constructGaussianW(object = object, phi = object@bandwidths[iKernel],
                                   NN4Kernel = NN4Kernel, k.NN4Kernel = k.NN4Kernel,
                                   eigen.less = eigen.less, eigen.prop = eigen.prop,
                                   weighted = weighted)
    }
  } else if (kernel.arg == "single"){
    cat("##\t Constructing the single Gaussian kernel... \n")
    object <- constructGaussianW(object = object, phi = object@bandwidths,
                                 NN4Kernel = NN4Kernel, k.NN4Kernel = k.NN4Kernel,
                                 eigen.less = eigen.less, eigen.prop = eigen.prop,
                                 weighted = weighted)
  }
  # ----------------------------------------------------------------------------
  n.Sample <- length(object@samples)
  rank.tie <- match.arg(rank.tie)
  has_cov <- !is.null(object@covariateMat)
  cat("##\t Working on rank scores...\n")

  rankScore <- lapply(seq_len(n.Sample), function(iSample){
    Y <- as.matrix(object@geneExpr[[iSample]])

    if(isTRUE(rank.zero)) {
      # Row-wise ranks in one shot (ties averaged, same shape as Y)
      R <- rowRanks(Y, ties.method = rank.tie, preserveShape = TRUE)
    } else {
      G <- nrow(Y); N <- ncol(Y)

      # Build matrix of rank scores, zeros left as 0 (excluded from ranking)
      R <- matrix(0, nrow = G, ncol = N)
      # non-zero entries
      nz <- which(Y != 0, arr.ind = TRUE)
      if (nrow(nz) > 0L) {
        # row indices of non-zeros
        i <- nz[, 1]
        # col indices of non-zeros
        j <- nz[, 2]
        # non-zero values
        v <- Y[nz]
        # Rank ONLY within the non-zeros of each row (ties averaged)
        r <- ave(v, i, FUN = function(z) rank(z, ties.method = rank.tie))
        # Fill ranked non-zeros; zeros stay 0
        R[cbind(i, j)] <- r
      }
    }
    if (has_cov) {
      X <- object@covariateMat[[iSample]]
      XtX <- t(X) %*% X
      XtX_inv <- solve(XtX)
      P <- X %*% XtX_inv %*% t(X)
      R <- P %*% R
    } else {
      # Row-wise center & scale without apply()
      mu <- rowMeans2(R)
      sd <- rowSds(R)
      # guard against 0-variance rows
      sd[!is.finite(sd) | sd == 0] <- 1

      R <- sweep(R, 1, mu, FUN = "-")
      R <- sweep(R, 1, sd, FUN = "/")
    }

    cat(paste0("##\t Sample ", object@samples[iSample], " is complete. \n"))
    return(R)
  })
  object@rankScores <- rankScore

  return(object)
}

# ==============================================================================
# ==============================================================================
#' @title Run INSiGHT for multi-sample integative detection of SVGs
#' @description
#' Run INSiGHT hypothesis tests to detect spatially variable genes (SVGs)
#' for multiple samples.
#'
#' @param object INSiGHT object.
#' @param QuadForm.method A character string specifying the method
#' to compute p-values for a quadratic form.
#' @param p.adjust.method A character string specifying the correction methods
#' for FDR control, including "holm", "hochberg", "hommel", "bonferroni", "BH",
#' "BY", "fdr" and "none". See `p.adjust.methods` for more details.
#'
#' @return A INSiGHT object containing test results.
#'
#' @import RANN
#' @import Matrix
#' @import matrixStats
#' @import RSpectra
#' @import CompQuadForm
#' @import ACAT
#'
#' @keywords internal
#' @importFrom stats p.adjust
#' @importFrom stats ave
#'
#' @export
#'

testINSiGHT <- function(object,
                        QuadForm.method = c("davies", "liu"),
                        p.adjust.method = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")) {

  n.Sample <- length(object@samples)
  QuadForm.method <- match.arg(QuadForm.method)
  p.adjust.method <- match.arg(p.adjust.method)

  # ----------------------------------------------------------------------------
  p_values_list <- lapply(1:length(object@bandwidths), function(iKernel) {
    cat(paste0('##\t', "Testing with Gaussian kernel ", iKernel, "...\n"))

    # Compute U = 2 * \sum_{i < j} E_{ij} * W_{ij}
    Ug.list <- vector("list", n.Sample)
    for (iSample in seq_len(n.Sample)) {
      Wc <- summary(object@W.info$W.mat[[iKernel]][[iSample]])
      i_nz <- Wc$i; j_nz <- Wc$j; x_nz <- Wc$x

      Ug.per.gene <- apply(object@rankScores[[iSample]], 1, function(y) {
        2 * sum(y[i_nz] * y[j_nz] * x_nz)
      })
      Ug.list[[iSample]] <- Ug.per.gene
    }
    # Sum over samples (vectorized)
    Ug.allSample <- Reduce(`+`, Ug.list)
    # ----------------------------------------------------------------------------
    # Combine all eigenvalues
    W.ev.combined <- c()
    for (iSample in 1:n.Sample) {
      W.ev.combined <- c(W.ev.combined, object@W.info$W.eigenvalues[[iKernel]][[iSample]])
    }
    W.ev.combined <- sort(W.ev.combined, decreasing = T)
    # ----------------------------------------------------------------------------
    eigen.less <- object@features$eigen.less
    if (isTRUE(eigen.less)) {
      Ug.allSample <- Ug.allSample + sum(W.ev.combined)
    }

    if (QuadForm.method == "davies") {
      p_values.g <- sapply(Ug.allSample, function(u) {
        pv <- CompQuadForm::davies(q = u, lambda = W.ev.combined)$Qq
        if(pv == 0) pv <- 1e-200
        return(pv)
      })
    } else if (QuadForm.method == "liu") {
      p_values.g <- sapply(Ug.allSample, function(u) {
        pv <- CompQuadForm::liu(q = u, lambda = W.ev.combined)
        if(pv == 0) pv <- 1e-200
        return(pv)
      })
    }

    return(p_values.g)

  })

  # ----------------------------------------------------------------------------
  kernel.option <- object@features$kernel.option
  if (kernel.option == "mixture") {
    p_values_matrix <- do.call(cbind, p_values_list)
    colnames(p_values_matrix) <- paste0("p.values_kernel_", 1:length(object@bandwidths))
    p_values_combined <- apply(p_values_matrix, MARGIN = 1, function(pvs){
      pv_comb <- ACAT::ACAT(Pvals = pvs)
    })
    p_values_combined <- as.vector(p_values_combined)
    p_values_adj <- p.adjust(p_values_combined, method = p.adjust.method)

    output <- data.frame(
      ID = object@genes,
      p.values_adjusted = p_values_adj,
      p.values_combined = p_values_combined
    )
    output <- cbind(output, p_values_matrix)
    object@testResults <- output
  } else if (kernel.option == "single") {
    p_values.g <- p_values_list[[1]]
    p_values_adj <- p.adjust(p_values.g, method = p.adjust.method)

    output <- data.frame(
      ID = object@genes,
      p.values_adjusted = p_values_adj,
      p.values_combined = p_values.g
    )
    object@testResults <- output
  }

  return(object)

}
