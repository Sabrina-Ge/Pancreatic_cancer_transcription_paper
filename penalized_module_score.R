library(Seurat)
library(dplyr)
library(Matrix)

penalty_transform <- function(x, k = 20, x0 = 0.625) {
  1 / (1 + exp(-k * (x - x0)))
}

LengthCheck <- function(values) {
  sapply(values, function(x) length(x) > 0)
}

remove_subtype_cols <- function(dataset, subtype_names=c("Basal1", "Basal2", "Classical1", "Classical2")){
  existing_columns_to_remove <- intersect(subtype_names, colnames(dataset@meta.data))
  
  # Remove only the existing columns
  dataset@meta.data <- dataset@meta.data %>%
    dplyr::select(-all_of(existing_columns_to_remove))
  
  return(dataset)
}

AddPenalizedModuleScore <- function(
    object,
    features,
    pool = NULL,
    nbin = 24,
    ctrl = 100,
    assay = NULL,
    name = 'Cluster',
    seed = 1,
    search = FALSE,
    slot = 'data',
    k = 25,
    x0 = 0.25,
    ...
) {
  if (!is.null(x = seed)) {
    set.seed(seed = seed)
  }
  assay.old <- DefaultAssay(object = object)
  assay <- assay %||% assay.old
  DefaultAssay(object = object) <- assay
  assay.data <- GetAssayData(object = object, assay = assay, slot = slot)
  features.old <- features

  if (is.null(x = features)) {
    stop("Missing input feature list")
  }
  features <- lapply(
    X = features,
    FUN = function(x) {
      missing.features <- setdiff(x = x, y = rownames(x = object))
      if (length(x = missing.features) > 0) {
        warning(
          "The following features are not present in the object: ",
          paste(missing.features, collapse = ", "),
          ifelse(
            test = search,
            yes = ", attempting to find updated synonyms",
            no = ", not searching for symbol synonyms"
          ),
          call. = FALSE,
          immediate. = TRUE
        )
        if (search) {
          tryCatch(
            expr = {
              updated.features <- UpdateSymbolList(symbols = missing.features, ...)
              names(x = updated.features) <- missing.features
              for (miss in names(x = updated.features)) {
                index <- which(x == miss)
                x[index] <- updated.features[miss]
              }
            },
            error = function(...) {
              warning(
                "Could not reach HGNC's gene names database",
                call. = FALSE,
                immediate. = TRUE
              )
            }
          )
          missing.features <- setdiff(x = x, y = rownames(x = object))
          if (length(x = missing.features) > 0) {
            warning(
              "The following features are still not present in the object: ",
              paste(missing.features, collapse = ", "),
              call. = FALSE,
              immediate. = TRUE
            )
          }
        }
      }
      return(intersect(x = x, y = rownames(x = object)))
    }
  )
  feature.length <- length(x = features)

  if (!all(LengthCheck(values = features))) {
    warning(paste(
      'Could not find enough features in the object from the following feature lists:',
      paste(names(x = which(x = !LengthCheck(values = features)))),
      'Attempting to match case...'
    ))
    features <- lapply(
      X = features.old,
      FUN = CaseMatch,
      match = rownames(x = object)
    )
  }
  if (!all(LengthCheck(values = features))) {
    stop(paste(
      'The following feature lists do not have enough features present in the object:',
      paste(names(x = which(x = !LengthCheck(values = features)))),
      'exiting...'
    ))
  }
  pool <- pool %||% rownames(x = object)
  data.avg <- Matrix::rowMeans(x = assay.data[pool, ])
  data.avg <- data.avg[order(data.avg)]
  data.cut <- cut_number(x = data.avg + rnorm(n = length(data.avg))/1e30, n = nbin, labels = FALSE, right = FALSE)
  names(x = data.cut) <- names(x = data.avg)
  ctrl.use <- vector(mode = "list", length = feature.length)
  for (i in 1:feature.length) {
    features.use <- features[[i]]
    for (j in 1:length(x = features.use)) {
      ctrl.use[[i]] <- c(
        ctrl.use[[i]],
        names(x = sample(
          x = data.cut[which(x = data.cut == data.cut[features.use[j]])],
          size = ctrl,
          replace = FALSE
        ))
      )
    }
  }
  ctrl.use <- lapply(X = ctrl.use, FUN = unique)
  ctrl.scores <- matrix(
    data = numeric(length = 1L),
    nrow = length(x = ctrl.use),
    ncol = ncol(x = object)
  )
  for (i in 1:length(ctrl.use)) {
    features.use <- ctrl.use[[i]]
    ctrl.scores[i, ] <- Matrix::colMeans(x = assay.data[features.use, ])
  }
  features.scores <- matrix(
    data = numeric(length = 1L),
    nrow = feature.length,
    ncol = ncol(x = object)
  )
  for (i in 1:feature.length) {
    features.use <- features[[i]]
    data.use <- assay.data[features.use, , drop = FALSE]
    features.scores[i, ] <- Matrix::colMeans(x = data.use)
  }
  feature_penalties <- matrix(0, nrow = length(features), ncol = ncol(assay.data))
  rownames(feature_penalties) <- paste0("Module", seq_along(features))
  colnames(feature_penalties) <- colnames(assay.data)
  for (i in seq_along(features)) {
    # Extract the current module and control gene sets
    module_genes <- features[[i]]
    control_genes <- ctrl.use[[i]]

    # Ensure genes are present in the assay data
    module_genes <- intersect(module_genes, rownames(assay.data))
    control_genes <- intersect(control_genes, rownames(assay.data))

    # Get the expression values for module and control genes
    module_expr <- assay.data[module_genes, , drop = FALSE]  # Expression of module genes
    control_expr <- assay.data[control_genes, , drop = FALSE]  # Expression of control genes

    # Calculate the average control expression per cell
    control_avg <- Matrix::colMeans(control_expr)

    # Compare each module gene's expression to the control average for each cell
    higher_than_control <- module_expr > control_avg

    # Calculate the percentage of module genes exceeding control average for each cell
    feature_penalties[i, ] <- penalty_transform(Matrix::colMeans(higher_than_control), k = k, x0 = x0)
  }

  features.scores.use <- (features.scores - ctrl.scores) * feature_penalties
  rownames(x = features.scores.use) <- name
  features.scores.use <- as.data.frame(x = t(x = features.scores.use))
  rownames(x = features.scores.use) <- colnames(x = object)
  object[[colnames(x = features.scores.use)]] <- features.scores.use
  CheckGC()
  DefaultAssay(object = object) <- assay.old
  return(object)
}

add_penalized_module_score <- function(dataset, genelists, k = 25, x0 = 0.25, nbin = 24) {
  # Create a vector of new names based on the names of the genelists
  names_to_add <- make.names(names(genelists))

  # Build a regex pattern that matches any of the names followed by one or more digits (current scores)
  pattern_to_remove <- paste0("(", paste(names_to_add, collapse = "|"), ")\\d+$")

  dataset <- dataset %>% remove_subtype_cols()
  # Add new module scores with updated genelists
  dataset <- dataset %>%
    AddPenalizedModuleScore(genelists, name = names_to_add, k = k, x0 = x0)

  # Remove digits from the end of feature names in the metadata
  dataset@meta.data <- dataset@meta.data %>%
    rename_with(~gsub("\\d+$", "", .), matches(pattern_to_remove))

  return(dataset)
}


