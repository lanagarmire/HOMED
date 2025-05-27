#' Estimate Cell Type Proportions Using Hierarchical Deconvolution
#'
#' This function performs hierarchical deconvolution by estimating cell type proportions at two levels (level 1 and level 2).
#' It first estimates the proportions for level 1 cell types, and then for level 2 cell types. The results from both levels are merged
#' to generate the final estimated proportions for each sample.
#'
#' @param celltypes A named list where each element contains the names of cell types at level 1 and their corresponding subtypes at level 2.
#'                  For example, `list("Tcell" = c("Treg", "Th1", "Th2"), "Bcell" = c("Naive", "Plasma"))`.
#' @param deconv_beta A matrix of methylation beta values (CpG by sample), which will be used to calculate the cell proportions.
#' @param reference_obj A list containing reference matrices for both level 1 and level 2 deconvolution.
#'                      `reference_obj$layer1_reference_list` and `reference_obj$layer2_reference_list` should each be a list where
#'                      each element corresponds to a reference matrix for a specific cell type.
#'
#' @return A list containing the estimated cell type proportions. The `deconvolution_proportion` element of the list contains
#'         the merged proportions for both level 1 and level 2 cell types.
#'
#' @importFrom dplyr bind_cols
#'
#' @export
#'
#' @examples
#' \dontrun{
#' estimated_proportions <- HOMED_Estimate(celltypes, deconv_beta, reference_obj)
#' }

HOMED_Estimate <- function(celltypes, deconv_beta, reference_obj) {

  layer1_proportion_list <- list()
  for (i in names(celltypes)) {
    message(paste0("Level 1 deconvolution for ",i))
    layer1_pred_prop = calculate_cell_proportions(deconv_beta,reference_obj$layer1_reference_list[[i]],i)
    layer1_proportion_list[[i]] <- as.data.frame(layer1_pred_prop)
  }
  wide_df_layer1<- as.data.frame(bind_cols(lapply(layer1_proportion_list, function(x) x$layer1_pred_prop)))
  rownames(wide_df_layer1) = colnames(deconv_beta)


  layer2_proportion_list <- list()
  non_null_layer2_cell_types <- names(celltypes)[sapply(celltypes, function(x) !is.null(x))]
  for (cell_type in non_null_layer2_cell_types) {
    subtypes <- celltypes[[cell_type]]
    for(i in subtypes){
      message(paste0("Level 2 deconvolution for ",i))
      layer2_pred_prop = calculate_cell_proportions(deconv_beta,reference_obj$layer2_reference_list[[i]], i)
      layer2_proportion_list[[i]] <- as.data.frame(layer2_pred_prop)
    }
  }
  wide_df_layer2 <- as.data.frame(bind_cols(lapply(layer2_proportion_list, function(x) x$layer2_pred_prop)))
  rownames(wide_df_layer2) = colnames(deconv_beta)

  # merge result together
  final_prop <- merge_layers(wide_df_layer1, wide_df_layer2, celltypes)
  return(list(deconvolution_proportion = final_prop))
}
