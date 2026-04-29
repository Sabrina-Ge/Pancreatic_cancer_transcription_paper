library(SingleR)

#' Subtype single cells through SingleR, using bulk tumours as reference.
#'
#' @param mat Count matrix where genes (HGNC format) are rows and cells are columns. LogCPM values are recommended. Can be a sparse matrix.
#' @param reference_file Filename of reference bulk tumour RData containing the objects ref_mat (LogTPM HGNC gene by sample matrix) and ref_labels (named character vector, sample: subtype label). Default "bulk_reference.RData".
#' @param find_mixed Classify cells as  Mixed. Default TRUE.
#' @param ... Additional arguments passed to SingleR
#'
#' @return List with the following:
#' 1) data: data.frame where rows are cells and columns represent the subtype, normalized confidence score and normalized subtype scores
#' 2) singler_pred: SingleR raw output
#'
singler_subtype <- function(mat, reference_file = "references/singler_reference.RData", find_mixed = TRUE, ...) {
  
  load(reference_file)
  if (!all(colnames(ref_mat) == names(ref_labels))) {
    stop("Mismatched sample names in reference file. Please check if you are using the correct reference.")
  }
  if (!("Classical1" %in% ref_labels) | !("Basal2" %in% ref_labels) | "hybrid" %in% ref_labels) {
    stop("Reference subtype labels do not match expected. Please check if you are using the correct reference.")
  }
  
  overlap <- length(intersect(rownames(mat), rownames(ref_mat))) / nrow(ref_mat)
  if (overlap < 0.2) {
    message(
      "Reference and query gene overlap is low (",
      round(overlap * 100, 2),
      "%). Please check if query single cell matrix contains genes (as rownames) in HGNC format."
    )
  }
  
  pred <- SingleR(test = mat, ref = ref_mat, labels = ref_labels, ...)
  
  # Normalize SingleR scores
  score_metadata <- as.data.frame(pred$scores)
  colnames(score_metadata) <- paste0(colnames(score_metadata), "_score")
  rownames(score_metadata) <- rownames(pred)
  score_conf <- rowMeans(score_metadata)
  score_metadata_norm <- score_metadata - score_conf
  
  data <- cbind(subtype = pred$labels,
                conf = score_conf,
                score_metadata_norm)
  
  if (find_mixed) {
    mixed_cutoff <- 0.02
    # Cells identified as Mixed if originally Classical1 or Basal2, and similar in score for both.
    subtype_with_mixed <- ifelse(
      data$subtype %in% c("Classical1", "Basal2") &
        abs(data$Classical1_score - data$Basal2_score) < mixed_cutoff,
      "Mixed",
      data$subtype
    )
    data$subtype <- subtype_with_mixed
  }
  
  return(list(data = data, singler_pred = pred))
  
}
