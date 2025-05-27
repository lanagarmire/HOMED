#' Construct Hierarchical Reference Matrices for Cell Type Deconvolution
#'
#' This function creates hierarchical reference matrices for deconvolution by performing probe selection and optimization
#' of coefficients for both level 1 and level 2 cell types. It uses a combination of `limma` for probe selection and an
#' optimization method (IDOL) to refine the reference matrices based on training data.
#'
#' @param FACS_beta A matrix of methylation data for the FACS samples (CpG by sample).
#' @param FACS_pd A data frame containing phenotype data for the FACS samples, including cell type annotations.
#' @param training_beta A matrix of training methylation data (CpG by sample) used for model optimization.
#' @param training_prop A matrix of cell type proportions (cell type by sample) for the training samples.
#' @param celltypes A named list of cell types, where each cell type is associated with a list of its subtypes at level 2.
#' @param p_val The p-value threshold for probe selection using `limma`.
#' @param numProbes The number of probes to select for each cell type.
#' @param maxIt A vector of two integers specifying the maximum number of iterations for level 1 and level 2 optimization.
#' @param libSize The size of the library to use for optimization.
#' @param rmse_improve_thresh The threshold for improvement in RMSE between iterations to stop the optimization.
#'
#' @return A list containing two elements:
#' \item{layer1_reference_list}{A list of reference matrices for level 1 cell types.}
#' \item{layer2_reference_list}{A list of reference matrices for level 2 cell types.}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' references <- HOMED_Reference(FACS_beta,FACS_pd,
#' training_beta,training_prop,celltypes,p_val = 1e-11,numProbes = 500)
#' }
HOMED_Reference <- function(FACS_beta, FACS_pd, training_beta, training_prop, celltypes,
                            p_val= 0.05, numProbes=500, maxIt = c(20,20), libSize = 500, rmse_improve_thresh = 0.001) {

  probes_limma = pickProbes_limma(FACS_beta, FACS_pd, level = 1, p_val = p_val, numProbes = numProbes)

  layer1_reference_list <- list()
  summarized_proportions <- summarize_level1(celltypes, t(training_prop))

  for (i in colnames(probes_limma$coefEsts)) {
    message(paste0("Calculating level 1 reference for ",i))
    #IDOL optimizing level 1 probes
    layer1_reference <- IDOL_opt_rmse(candDMRFinderObject = probes_limma$coefEsts,
                                      trainingBetas = training_beta,
                                      trainingCovariates = summarized_proportions,
                                      maxIt = maxIt[1], libSize = libSize,
                                      IDOLobj = FALSE, rmse_improve_thresh = rmse_improve_thresh,
                                      rmse_for = i)
    layer1_reference_list[[i]] <- layer1_reference$`IDOL Optimized CoefEsts`
  }

  # Level 2 deconvolution
  layer2_reference_list <- list()

  non_null_layer2_cell_types <- names(celltypes)[sapply(celltypes, function(x) !is.null(x))]

  for (cell_type in non_null_layer2_cell_types) {

    subtypes <- celltypes[[cell_type]]
    sub_ids = rownames(FACS_pd[FACS_pd$cellType_level2%in%subtypes,])

    probes_limma_layer2 = pickProbes_limma(FACS_beta, FACS_pd, level=2, p_val= p_val, numProbes=numProbes)

    for(i in subtypes){
      message(paste0("Calculating level 2 reference for ",i))
      layer2_reference <- IDOL_opt_rmse(candDMRFinderObject = probes_limma_layer2$coefEsts,
                                        trainingBetas = training_beta,
                                        trainingCovariates = training_prop,
                                        maxIt = maxIt[2], libSize = libSize,
                                        IDOLobj = FALSE, rmse_improve_thresh = rmse_improve_thresh,
                                        rmse_for = i)
      layer2_reference_list[[i]] <- layer2_reference$`IDOL Optimized CoefEsts`
    }

  }

  # Return the final references for both layers
  return(list(layer1_reference_list = layer1_reference_list,
              layer2_reference_list = layer2_reference_list))
}
