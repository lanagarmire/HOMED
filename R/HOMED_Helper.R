#' Perform Limma-based Probe Selection for Deconvolution
#'
#' This function uses the `limma` package to select probes based on differential expression between specified cell types.
#'
#' @param beta Methylation beta matrix (CpG by sample).
#' @param pd A data frame containing the cell type annotation, must have a column with name 'CellType'
#' @param celltypes A vector of cell types to be used for probe selection.
#' @param p_val The p-value threshold for significance in the differential expression analysis.
#' @param n_prob The number of probes to select.
#' @param by Specifies the method for probe selection, either by log fold change (`'logfc'`) or adjusted p-value (`'adj.P.Val'`).
#'
#' @return A vector of selected probe names.
#'
#' @importFrom stats relevel model.matrix complete.cases
#' @importFrom limma lmFit eBayes topTable
#' @importFrom lumi beta2m
#'
#' @export
#'
#' @examples
#' \dontrun{
#' selected_probes <- deconv_limma(beta,pd,celltypes,p_val= 1e-11,n_prob= 500,by='logfc')
#' }
#Sys.setenv('_R_CHECK_SYSTEM_CLOCK_' = 0)
deconv_limma <- function(beta,pd,celltypes,p_val= 1e-11,n_prob= 500,by='logfc'){
  trainingProbes <- list()
  for(i in celltypes){
    pd$status <- NA
    pd[pd$CellType==i,]$status <- rep('interest',dim(pd[pd$CellType==i,])[1])
    pd[!pd$CellType==i,]$status <- rep('other',dim(pd[!pd$CellType==i,])[1])
    pd$status = relevel(factor(pd$status), ref="interest")
    design = model.matrix(~pd$status)
    fit = lmFit(beta2m(beta), design)
    fit = eBayes(fit)
    table_limma = topTable(fit, coef=2, number=dim(fit)[1])
    table_limma_complete = table_limma[is.finite(rowSums(table_limma)),]
    sig_limma = table_limma_complete[table_limma_complete$adj.P.Val<p_val,]

    if(by=='logfc'){
      hypo_limma = sig_limma[order(sig_limma[, "logFC"], decreasing = FALSE), ]
      hyper_limma = sig_limma[order(sig_limma[, "logFC"], decreasing = TRUE), ]
      probes = c(rownames(hyper_limma)[seq_len(floor(n_prob/2))],
                 rownames(hypo_limma)[seq_len(floor(n_prob/2))])
    }else{
      probes = rownames(sig_limma[order(sig_limma[, "adj.P.Val"], decreasing = FALSE), ])[1:n_prob]
    }

    trainingProbes = append(trainingProbes,list(probes))
  }
  trainingProbes = unique(unlist(trainingProbes))
  trainingProbes = trainingProbes[complete.cases(trainingProbes)]
  return(trainingProbes)
}




#' Validate Cell Type Proportions
#'
#' This function validates the estimated cell type proportions using linear models or mixed effects models.
#' It computes the F-statistic and associated p-value for each CpG site.
#'
#' @param Y A matrix of methylation values (CpG by sample).
#' @param pheno A data frame containing phenotype data.
#' @param modelFix The fixed effect model formula.
#' @param modelBatch Optional; the batch effect model formula (if any).
#' @param L.forFstat The matrix used for calculating F-statistics.
#' @param verbose Logical; whether to print progress messages.
#'
#' @return A list containing coefficient estimates, variance-covariance matrices, F-statistics, and p-values.
#'
#' @importFrom stats model.matrix lm vcov pf
#' @importFrom nlme lme getVarCov
#'
#' @export
#'
#' @examples
#' \dontrun{
#' validation_results <- validationCellType(Y, pheno, modelFix,
#' modelBatch=NULL,L.forFstat = NULL, verbose = FALSE)
#' }
validationCellType <- function(Y, pheno, modelFix, modelBatch=NULL,
                               L.forFstat = NULL, verbose = FALSE){
  N <- dim(pheno)[1]
  pheno$y <- rep(0, N)
  xTest <- model.matrix(modelFix, pheno)
  sizeModel <- dim(xTest)[2]
  M <- dim(Y)[1]

  if (is.null(L.forFstat)) {
    # NOTE: All non-intercept coefficients
    L.forFstat <- diag(sizeModel)[-1,]
    colnames(L.forFstat) <- colnames(xTest)
    rownames(L.forFstat) <- colnames(xTest)[-1]
  }

  # Initialize various containers
  sigmaResid <- sigmaIcept <- nObserved <- nClusters <- Fstat <- rep(NA, M)
  coefEsts <- matrix(NA, M, sizeModel)
  coefVcovs <- list()

  if (verbose) cat("[validationCellType] ")
  # Loop over each CpG
  for (j in seq_len(M)) {
    # Remove missing methylation values
    ii <- !is.na(Y[j, ])
    nObserved[j] <- sum(ii)
    pheno$y <- Y[j,]

    if (j %% round(M / 10) == 0 && verbose) cat(".") # Report progress

    # Try to fit a mixed model to adjust for plate
    try({
      if (!is.null(modelBatch)) {
        fit <- try(
          lme(modelFix, random = modelBatch, data = pheno[ii, ]))
        # NOTE: If LME can't be fit, just use OLS
        OLS <- inherits(fit, "try-error")
      } else {
        OLS <- TRUE
      }

      if (OLS) {
        fit <- lm(modelFix, data = pheno[ii, ])
        fitCoef <- fit$coef
        sigmaResid[j] <- summary(fit)$sigma
        sigmaIcept[j] <- 0
        nClusters[j] <- 0
      } else {
        fitCoef <- fit$coef$fixed
        sigmaResid[j] <- fit$sigma
        sigmaIcept[j] <- sqrt(getVarCov(fit)[1])
        nClusters[j] <- length(fit$coef$random[[1]])
      }
      coefEsts[j,] <- fitCoef
      coefVcovs[[j]] <- vcov(fit)

      useCoef <- L.forFstat %*% fitCoef
      useV <- L.forFstat %*% coefVcovs[[j]] %*% t(L.forFstat)
      Fstat[j] <- (t(useCoef) %*% solve(useV, useCoef)) / sizeModel
    })
  }
  if (verbose) cat(" done\n")

  # Name the rows so that they can be easily matched to the target data set
  rownames(coefEsts) <- rownames(Y)
  colnames(coefEsts) <- names(fitCoef)
  degFree <- nObserved - nClusters - sizeModel + 1

  # Get P values corresponding to F statistics
  Pval <- 1 - pf(Fstat, sizeModel, degFree)

  list(
    coefEsts = coefEsts,
    coefVcovs = coefVcovs,
    modelFix = modelFix,
    modelBatch = modelBatch,
    sigmaIcept = sigmaIcept,
    sigmaResid = sigmaResid,
    L.forFstat = L.forFstat,
    Pval = Pval,
    orderFstat = order(-Fstat),
    Fstat = Fstat,
    nClusters = nClusters,
    nObserved = nObserved,
    degFree = degFree)
}

#' Select Probes Using Limma for Cell Type-Specific Deconvolution
#'
#' This function selects the most informative probes for deconvolution by performing differential expression analysis
#' using `limma` at different hierarchical levels of cell types.
#'
#' @param p A matrix of methylation data (CpG by sample).
#' @param pd A data frame of phenotype data with cell type information.
#' @param cellTypes The cell types to include in the analysis.
#' @param level The hierarchical level for cell types (1 for primary, 2 for secondary).
#' @param numProbes The number of probes to select.
#' @param p_val The p-value threshold for significance.
#' @param probeSelect The selection criterion: `'logfc'`, `'adj.P.Val'`, or `'both'`.
#' @param mylist Optional list of previously selected probes.
#' @param by The method to use for probe selection (`'logfc'` or `'adj.P.Val'`).
#'
#' @return A vector of selected probe names.
#'
#' @importFrom genefilter rowFtests rowttests
#' @importFrom stats as.formula model.matrix
#' @importFrom matrixStats colMeans2 rowMeans2 rowRanges
#'
#' @export
#'
#' @examples
#' \dontrun{
#' selected_probes <- pickProbes_limma(methylation_data, phenotype_data,
#' cellTypes, level = 1, numProbes = 500)
#' }

pickProbes_limma <- function(p, pd, cellTypes = NULL, level=1, numProbes = 500,
                             probeSelect = 'both',p_val,mylist=NULL,by='logfc') {
  if (level==1) {
    pd$CellType = pd$cellType_level1
    cellTypes = unique(pd$cellType_level1)
  }else if (level==2) {
    pd$CellType = pd$cellType_level2
    cellTypes = unique(pd$cellType_level2)
  }else{
    stop("Input must be integer 1 or 2 indicating the heirchical level.")
  }

  splitit <- function(x) {
    split(seq_along(x), x)
  }

  if (!is.null(cellTypes)) {
    if (!all(cellTypes %in% pd$CellType))
      stop("elements of argument 'cellTypes' is not part of ",
           "'mSet$CellType'")
    keep <- which(pd$CellType %in% cellTypes)
    pd <- pd[keep,]
    p <- p[,keep]
  }
  # NOTE: Make cell type a factor
  pd$CellType <- factor(pd$CellType, levels = cellTypes)
  ffComp <- rowFtests(p, pd$CellType)
  prof <- vapply(
    X = splitit(pd$CellType),
    FUN = function(j) rowMeans2(p, cols = j),
    FUN.VALUE = numeric(nrow(p)))
  r <- rowRanges(p)
  compTable <- cbind(ffComp, prof, r, abs(r[, 1] - r[, 2]))
  names(compTable)[1] <- "Fstat"
  names(compTable)[c(-2, -1, 0) + ncol(compTable)] <-
    c("low", "high", "range")
  tIndexes <- splitit(pd$CellType)
  tstatList <- lapply(tIndexes, function(i) {
    x <- rep(0,ncol(p))
    x[i] <- 1
    return(rowttests(p, factor(x)))
  })

  if(is.null(mylist)){
    trainingProbes <- deconv_limma(p,pd,celltypes=cellTypes,p_val,n_prob=numProbes,by=by)
  }else{
    trainingProbes <- rownames(p)[rownames(p)%in%mylist]
  }

  p <- p[trainingProbes,]

  pMeans <- colMeans2(p)
  names(pMeans) <- pd$CellType

  form <- as.formula(
    sprintf("y ~ %s - 1", paste(levels(pd$CellType), collapse = "+")))
  phenoDF <- as.data.frame(model.matrix(~ pd$CellType - 1))
  colnames(phenoDF) <- sub("^pd\\$CellType", "", colnames(phenoDF))
  if (ncol(phenoDF) == 2) {
    # Two group solution
    X <- as.matrix(phenoDF)
    coefEsts <- t(solve(t(X) %*% X) %*% t(X) %*% t(p))
  } else {
    # > 2 groups solution
    tmp <- validationCellType(Y = p, pheno = phenoDF, modelFix = form)
    coefEsts <- tmp$coefEsts
  }

  list(
    coefEsts = coefEsts,
    compTable = compTable,
    sampleMeans = pMeans)
}

#' Optimize Coefficients Using RMSE for Deconvolution
#'
#' This function optimizes methylation deconvolution coefficients by minimizing the RMSE between predicted and observed proportions.
#'
#' @param candDMRFinderObject A list or object containing candidate probes and coefficient estimates.
#' @param trainingBetas A matrix of methylation values for training samples.
#' @param trainingCovariates A matrix or data frame of the training covariates (e.g., cell proportions).
#' @param libSize The size of the library to use for optimization.
#' @param maxIt The maximum number of iterations for the optimization.
#' @param numCores The number of CPU cores to use for parallel processing.
#' @param IDOLobj Logical; whether to use the IDOL optimization method.
#' @param rmse_improve_thresh The threshold for RMSE improvement to stop the optimization.
#' @param rmse_for The cell type to optimize RMSE for.
#'
#' @return A list containing the optimized probes and coefficients, as well as RMSE and R2 values for each iteration.
#'
#' @importFrom EpiDISH epidish
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom stats cor sd
#' @importFrom parallel makeCluster stopCluster
#' @export
#'
#' @examples
#' \dontrun{
#' idol_results <- IDOL_opt_rmse(candDMRFinderObject,
#' trainingBetas, trainingCovariates, libSize = 500, maxIt = 200)
#' }
IDOL_opt_rmse = function(candDMRFinderObject, trainingBetas, trainingCovariates,
                         libSize = 500, maxIt = 200, numCores = 4, IDOLobj=TRUE,
                         rmse_improve_thresh = 0.0001,rmse_for = c("Hofbauer") ) {

  cl <- makeCluster(numCores)
  registerDoParallel(cl)

  expit = function(w) exp(w)/(1 + exp(w))

  R2_celltype_level = function(obs, pred) {
    pcc = numeric(ncol(obs))
    for (i in 1:ncol(obs)) {
      pcc[i] = cor(obs[, i], pred[, i], method = "pearson")
    }
    mean(pcc, na.rm = TRUE)
  }


  RMSE_celltype_level = function(obs, pred) {
    cols_use = intersect(rmse_for, colnames(obs))
    if (length(cols_use) < 1) stop("Missing cell type in the data.")
    obs = obs[, cols_use, drop = FALSE]
    pred = pred[, cols_use, drop = FALSE]
    rmse = numeric(ncol(obs))
    for (i in 1:ncol(obs)) {
      rmse[i] = sqrt(mean((obs[, i] - pred[, i])^2, na.rm = TRUE))
    }
    mean(rmse, na.rm = TRUE)
  }

  polar = function(x, y, scale = 1) {
    r = sqrt(x^2 + y^2)
    theta = atan2(y, x)
    r * cos(theta - (scale * pi/4))
  }

  if (IDOLobj) {
    trainingProbes1 = rownames(candDMRFinderObject$coefEsts)
    coefEsts = candDMRFinderObject$coefEsts
  } else {
    trainingProbes1 = rownames(candDMRFinderObject)
    coefEsts = candDMRFinderObject
  }

  P = length(trainingProbes1)
  ProbVector = rep(1/P, P)
  V = libSize
  B = maxIt

  R2_best = 0
  RMSE_best = 10000

  R2Vals = numeric(B)
  RMSEVals = numeric(B)

  cellTypes = colnames(coefEsts)
  K = length(cellTypes)

  if (!all(cellTypes %in% colnames(trainingCovariates))) {
    stop("Cell type names in covariates do not match coefEsts columns")
  }

  for (i in 1:B) {
    Probes = sample(1:P, V, prob = ProbVector)
    CpGNames = trainingProbes1[Probes]
    Beta = coefEsts[CpGNames, ]

    ctpred = data.frame(suppressMessages(epidish(trainingBetas[CpGNames, ],
                                                 Beta,
                                                 method = "CP",
                                                 maxit = 50,
                                                 nu.v = c(0.25, 0.5, 0.75),
                                                 constraint = "inequality")$estF))
    omega.tilde = as.matrix(ctpred)
    omega.obs = trainingCovariates[, cellTypes]

    R2_ct = R2_celltype_level(omega.obs, omega.tilde)
    RMSE_ct = RMSE_celltype_level(omega.obs, omega.tilde)

    Perform.q = foreach(j = 1:length(CpGNames)) %dopar% {
      Beta.q = Beta[CpGNames[-j], ]
      ctpred.q = data.frame(suppressMessages(epidish(trainingBetas[CpGNames[-j], ],
                                                     Beta.q,
                                                     method = "CP",
                                                     maxit = 50,
                                                     nu.v = c(0.25, 0.5, 0.75),
                                                     constraint = "inequality")$estF))
      omega.tilde.q = as.matrix(ctpred.q)
      R2.q = R2_celltype_level(omega.obs, omega.tilde.q)
      RMSE.q = RMSE_celltype_level(omega.obs, omega.tilde.q)
      c(R2.q, RMSE.q)
    }

    R2.q = unlist(Perform.q)[seq(1, 2 * V, by = 2)]
    RMSE.q = unlist(Perform.q)[seq(2, 2 * V, by = 2)]

    rmse.dq = (RMSE_ct - RMSE.q) * (-1)
    norm.rmse = rmse.dq / sd(rmse.dq)
    r2.dq = R2_ct - R2.q
    norm.r2 = r2.dq / sd(r2.dq)

    p1 = polar(norm.rmse, norm.r2)
    for (j in 1:length(Probes)) {
      p0 = ProbVector[[Probes[j]]]
      ProbVector[[Probes[j]]] = expit(p1[j]) * p0 + p0 / 2
    }
    ProbVector = ProbVector / sum(ProbVector)

    #if ((RMSE_best - RMSE_ct > rmse_improve_thresh) & R2_ct >= R2_best) {
    if ((RMSE_best - RMSE_ct > rmse_improve_thresh)){#} | (R2_ct > R2_best & RMSE_ct < RMSE_best)){

      RMSE_best = RMSE_ct
      R2_best = R2_ct

      #print(paste("Iteration: ", i, " Cell-RMSE = ", round(RMSE_ct, 4), "; Cell-PCC = ", round(R2_ct, 4), sep = ""))

      IDOL.optim.DMRs = CpGNames
      IDOL.optim.coefEsts = coefEsts[CpGNames, ]
    }

    RMSEVals[i] = RMSE_ct
    R2Vals[i] = R2_ct
  }

  stopCluster(cl)
  IDOLObjects = list(IDOL.optim.DMRs,
                     IDOL.optim.coefEsts,
                     RMSEVals,
                     R2Vals,
                     B,
                     V)
  names(IDOLObjects) = c("IDOL Optimized Library",
                         "IDOL Optimized CoefEsts",
                         "RMSE_CelltypeLevel",
                         "PCC_CelltypeLevel",
                         "Number of Iterations",
                         "LibrarySize")

  print(paste("Final Library - ",rmse_for,": Avg Cell-RMSE = ", round(RMSE_best, 4),"; Cell-PCC = ", round(R2_best, 4), sep = ""))
  return(IDOLObjects)
}

#' Calculate Cell Type Proportions Using Deconvolution Method
#'
#' This function calculates the cell type proportion for a given 2-level structured cell type using a deconvolution method.
#' It uses the `epidish` package to perform deconvolution and returns the estimated cell proportion for a specified cell type.
#'
#' @param beta A matrix of methylation data (CpG by sample).
#' @param reference A reference matrix (cell type-specific methylation profiles).
#' @param celltype A character string specifying the cell type for which the proportion is calculated.
#'
#' @return A numeric vector of predicted cell proportions for the specified cell type.
#' @export
#'
#' @importFrom EpiDISH epidish
#'
#' @examples
#' \dontrun{
#' cell_proportion <- calculate_cell_proportions(beta_data, reference_data)
#' }
calculate_cell_proportions = function(beta, reference, celltype) {
  cell_prop_res = suppressMessages(epidish(beta[rownames(reference), ], reference, method = "CP")$estF[,celltype])
  return(cell_prop_res)
}


#' Summarize Level 1 Cell Type Proportions
#'
#' This function summarizes cell type proportions at level 1 by summing the proportions of subtypes for each cell type.
#' The resulting data frame contains the summed proportions for each cell type across all samples.
#'
#' @param cell_types_list A list where each element contains the names of subtypes corresponding to a level 1 cell type.
#' @param prop_data A data frame or matrix of cell proportions (samples by cell type).
#'
#' @return A data frame with summarized cell proportions for each level 1 cell type.
#' @export
#'
#' @examples
#' \dontrun{
#' summarized_proportions <- summarize_level1(celltypes, prop_data)
#' }
summarize_level1 <- function(cell_types_list, prop_data) {
  summarized_data <- data.frame(matrix(NA, ncol = length(names(cell_types_list)),nrow=dim(prop_data)[2]))  # Initialize the output
  colnames(summarized_data) = names(cell_types_list)
  rownames(summarized_data) = colnames(prop_data)

  for (cell_type in names(cell_types_list)) {
    subtypes <- cell_types_list[[cell_type]]  # Get subtypes for the current Level 1 cell type

    if (is.null(subtypes)) {
      # If no subtypes, keep the original column as is
      summarized_data[,cell_type] <- t(prop_data)[,cell_type]
    } else {
      # Check if subtypes are valid columns in the data frame
      valid_subtypes <- subtypes[subtypes %in% rownames(prop_data)]

      if (length(valid_subtypes) > 0) {
        # Sum the proportions of valid subtypes
        summed_proportion <- rowSums(t(prop_data)[,valid_subtypes])
        summarized_data[,cell_type] <- summed_proportion  # Replace with summed values
      } else {
        # If no valid subtypes, return NA for that cell type
        summarized_data[,cell_type] <- NA
        warning(paste("No valid subtypes found for", cell_type))
      }
    }
  }
  return(summarized_data)
}


#' Merge Cell Type Proportions from Two Layers
#'
#' This function merges cell type proportions from two layers (level 1 and level 2) of a hierarchical cell type model.
#' It normalizes the proportions of the subtypes at level 2 using the proportions at level 1.
#'
#' @param layer1_df A data frame of cell proportions for level 1 cell types.
#' @param layer2_df A data frame of cell proportions for level 2 cell types.
#' @param cell_types A list where names are cell types and values are vectors of subtypes at the respective levels.
#'
#' @return A data frame with merged and normalized cell type proportions from both layers.
#' @export
#'
#' @examples
#' \dontrun{
#' merged_data <- merge_layers(layer1_data, layer2_data, celltypes)
#' }
merge_layers <- function(layer1_df, layer2_df, cell_types) {

  group_1 <- names(cell_types)[sapply(cell_types, is.null)]
  group_2 <- names(cell_types)[sapply(cell_types, function(x) !is.null(x))]

  normalized_level2 = NULL

  for (cell_type in group_2) {
    cell_type = group_2[1]
    subtypes = cell_types[[cell_type]]
    level1_major = layer1_df[,cell_type]
    level2_sub = layer2_df[,subtypes]
    normalized_level2_temp = level2_sub/rowSums(level2_sub)*level1_major

    if (is.null(normalized_level2)) {
      normalized_level2 <- normalized_level2_temp  # Initialize for the first iteration
    } else {
      normalized_level2 <- cbind(normalized_level2, normalized_level2_temp)  # Add subsequent columns
    }
  }
  merged_result = cbind(layer1_df[,group_1],normalized_level2)
  return(merged_result)
}

#' Plot Deconvolution Results
#'
#' This function generates a plot of true vs predicted cell proportions for each cell type and computes
#' the Pearson Correlation Coefficient (PCC) and Concordance Correlation Coefficient (CCC) for each cell type.
#' The plot shows the correlation between observed and predicted proportions with calculated PCC and CCC values.
#'
#' @param beta A matrix of methylation data (CpG by sample).
#' @param ref A reference object containing optimized deconvolution coefficients.
#' @param true_prop A matrix or data frame of true cell proportions.
#' @param title The title for the plot.
#' @param idol Logical; whether to use IDOL-optimized coefficients or not.
#' @param method The deconvolution method to use (e.g., `'RPC'`,`'CP'`).
#' @param con The constraint for deconvolution (e.g., `'equality'`,`'inequality'`).
#' @param pred_prop Optional; a matrix of predicted cell proportions. If not provided, IDOL optimization will be used.
#' @param xylim The limits for the x and y axes of the plot.
#'
#' @return A list containing the PCC and CCC data frames, the plot, and predicted proportions.
#'
#' @importFrom EpiDISH epidish
#' @importFrom ggplot2 ggplot aes geom_point geom_abline facet_wrap labs theme_minimal ylim xlim scale_color_manual geom_text
#' @importFrom dplyr group_by left_join summarise distinct %>% mutate
#' @importFrom stats cor sd reshape
#' @importFrom RColorBrewer brewer.pal
#'
#' @export
#'
#' @examples
#' \dontrun{
#' HOMED_plot(beta_data, reference_data, true_cell_proportions, title = "Deconvolution Results")
#' }
HOMED_plot = function(beta,ref,true_prop,title='Deconvolution Results',
                             idol=TRUE,method='RPC',con='equality',
                             pred_prop=NULL,xylim=0.5){
  cell_type = NULL
  pcc = NULL
  ccc = NULL
  y_true = as.data.frame(t(true_prop))
  if(is.null(pred_prop)){
    if(idol){
      res_idol_epidish2 = epidish(beta, #GSE73377_methyl[common_cpg,],
                                  ref$`IDOL Optimized CoefEsts`,
                                  method = c(method),
                                  maxit = 50,nu.v = c(0.25, 0.5, 0.75),
                                  constraint = c(con))
    }else{
      res_idol_epidish2 = epidish(beta, #GSE73377_methyl[common_cpg,],
                                  ref$coefEsts,
                                  method = c(method),
                                  maxit = 50,nu.v = c(0.25, 0.5, 0.75),
                                  constraint = c(con))
    }
    y_pred = as.data.frame(res_idol_epidish2$estF)[,colnames(y_true)]
  }else{
    common_ct = intersect(colnames(y_true),colnames(pred_prop))
    y_pred = pred_prop[,common_ct]
    y_true = y_true[,common_ct]
  }
  y_true$sample <- rownames(y_true)
  y_pred$sample <- rownames(y_pred)

  # Convert to long format using reshape()
  y_true_long <- reshape(y_true,
                         varying = names(y_true)[-ncol(y_true)],
                         v.names = "y_true",
                         timevar = "cell_type",
                         times = names(y_true)[-ncol(y_true)],
                         direction = "long")

  y_pred_long <- reshape(y_pred,
                         varying = names(y_pred)[-ncol(y_pred)],
                         v.names = "y_pred",
                         timevar = "cell_type",
                         times = names(y_pred)[-ncol(y_pred)],
                         direction = "long")

  # Merge on sample and cell_type
  df <- merge(y_true_long[, c("sample", "cell_type", "y_true")],
              y_pred_long[, c("sample", "cell_type", "y_pred")],
              by = c("sample", "cell_type"))

  df <- df %>%
    mutate(
      cell_type = as.factor(cell_type),  # Ensure `cell_type` is a factor
    )

  # Function to calculate CCC
  ccc <- function(x, y) {
    rho <- cor(x, y, use = "complete.obs")
    x_mean <- mean(x, na.rm = TRUE)
    y_mean <- mean(y, na.rm = TRUE)
    x_sd <- sd(x, na.rm = TRUE)
    y_sd <- sd(y, na.rm = TRUE)

    numerator <- 2 * rho * x_sd * y_sd
    denominator <- x_sd^2 + y_sd^2 + (x_mean - y_mean)^2

    return(numerator / denominator)
  }


  cell_type_pcc_df <- df %>%
    group_by(cell_type) %>%
    summarise(pcc = cor(y_true, y_pred, method = "pearson", use = "complete.obs"))

  cell_type_ccc_df <- df %>%
    group_by(cell_type) %>%
    summarise(ccc = ccc(y_true, y_pred))

  cell_colors <- RColorBrewer::brewer.pal(n = length(unique(df$cell_type)), name = "Set2")

  df <- df %>%
    left_join(cell_type_pcc_df, by = "cell_type") %>%
    left_join(cell_type_ccc_df, by = "cell_type")

  p1 = ggplot(df, aes(x = y_true, y = y_pred, color = cell_type)) +
    geom_point(alpha = 0.8, size = 1.8, show.legend = FALSE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
    facet_wrap(~cell_type, scales = "free", ncol = 4) +
    theme_minimal() +
    ylim(0, xylim) + xlim(0, xylim) +
    scale_color_manual(values = cell_colors) +
    labs(
      x = "True Cell Proportion",
      y = "Predicted Cell Proportion",
      title = title
    )

  p1 = p1 +
    geom_text(data = distinct(df, cell_type, pcc),
              aes(x = 0.05, y = xylim-0.1, label = paste0("PCC = ", round(pcc, 2))),
              color = "black", size = 3, hjust = 0)+
    geom_text(data = distinct(df, cell_type, ccc),
              aes(x = 0.05, y = xylim-0.2, label = paste0("CCC = ", round(ccc, 2))),
              color = "black", size = 3, hjust = 0)

  p1
  return(list(cell_type_pcc_df,cell_type_ccc_df,p1,y_pred))
}

