#' Detect edge dryspots in spatial transcriptomics data (Enhanced Version)
#'
#' This function identifies edge dryspots and problem areas in spatial 
#' transcriptomics data by analyzing QC metrics and spatial patterns.
#' Enhanced with parameterized thresholds and improved sensitivity.
#'
#' @param spe A SpatialExperiment object containing spatial transcriptomics data
#' @param qc_metric Character string specifying the QC metric column name to analyze (default: "sum_gene")
#' @param samples Character string specifying the sample ID column name (default: "sample_id")
#' @param mad_threshold Numeric value for MAD threshold for outlier detection (default: 3)
#' @param edge_threshold Numeric threshold for edge detection (default: 0.6, reduced from 0.75)
#' @param min_cluster_size Minimum cluster size for morphological cleaning (default: 20, reduced from 40)
#' @param shifted Logical indicating whether to apply coordinate adjustment for hexagonal arrays (default: FALSE)
#' @param batch_var Character specifying batch variable for outlier detection ("slide", "sample_id", or "both", default: "both")
#' @param name Character string for naming output columns (default: "edge_dryspot")
#'
#' @return A SpatialExperiment object with additional columns in colData:
#'   \item{[name]_edge}{Logical indicating spots identified as edges}
#'   \item{[name]_problem_id}{Character identifying problem area clusters}
#'   \item{[name]_problem_size}{Numeric size of problem area clusters}
#'
#' @examples
#' # Enhanced usage with improved sensitivity
#' \dontrun{
#' spe <- detectEdgeDryspots(spe, edge_threshold = 0.5, min_cluster_size = 15)
#' 
#' # More aggressive detection
#' spe <- detectEdgeDryspots(spe, edge_threshold = 0.4, min_cluster_size = 10)
#'
#' }
#' @importFrom dplyr %>% sym
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom scuttle isOutlier
#' @export
detectEdgeDryspots <- function(
    spe, 
    qc_metric = "sum_gene",
    samples = "sample_id", 
    mad_threshold = 3,
    edge_threshold = 0.75,
    shifted = FALSE,
    batch_var = "both",
    name = "edge_dryspot", 
    min_cluster_size = 40) {
  
  # Input validation
  if (!inherits(spe, "SpatialExperiment")) {
    stop("Input data must be a SpatialExperiment or inherit from SpatialExperiment.")
  }
  if (!qc_metric %in% colnames(colData(spe))) {
    stop("qc_metric must be present in colData.")
  }
  if (!samples %in% colnames(colData(spe))) {
    stop("samples column must be present in colData.")
  }
  required_cols <- c("in_tissue", "array_row", "array_col")
  missing_cols <- required_cols[!required_cols %in% colnames(colData(spe))]
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  if (!is.numeric(mad_threshold) || mad_threshold <= 0) {
    stop("'mad_threshold' must be a positive numeric value.")
  }
  
  # Data preparation
  lg10_metric <- paste0("lg10_", qc_metric)
  colData(spe)[[lg10_metric]] <- log10(colData(spe)[[qc_metric]])
  
  # Outlier detection by batch variable
  if (batch_var %in% c("slide", "both")) {
    if ("slide" %in% colnames(colData(spe))) {
      outlier_slide_col <- paste0(qc_metric, "_3MAD_outlier_slide")
      colData(spe)[[outlier_slide_col]] <- scuttle::isOutlier(
        colData(spe)[[lg10_metric]], 
        subset = colData(spe)$in_tissue, 
        batch = colData(spe)$slide, 
        type = "lower", 
        nmads = mad_threshold
      )
    } else {
      warning("'slide' column not found, skipping slide-level outlier detection")
      colData(spe)[[paste0(qc_metric, "_3MAD_outlier_slide")]] <- FALSE
    }
  }

  if (batch_var %in% c("sample_id", "both")) {
    outlier_sample_col <- paste0(qc_metric, "_3MAD_outlier_sample")
    colData(spe)[[outlier_sample_col]] <- scuttle::isOutlier(
      colData(spe)[[lg10_metric]], 
      subset = colData(spe)$in_tissue, 
      batch = colData(spe)[[samples]], 
      type = "lower", 
      nmads = mad_threshold
    )
  }
  
  # Combine outlier results
  if (batch_var == "both") {
    colData(spe)[[paste0(qc_metric, "_3MAD_outlier_binary")]] <- 
      colData(spe)[[paste0(qc_metric, "_3MAD_outlier_slide")]] | 
      colData(spe)[[paste0(qc_metric, "_3MAD_outlier_sample")]]
  } else if (batch_var == "slide") {
    colData(spe)[[paste0(qc_metric, "_3MAD_outlier_binary")]] <- 
      colData(spe)[[paste0(qc_metric, "_3MAD_outlier_slide")]]
  } else {
    colData(spe)[[paste0(qc_metric, "_3MAD_outlier_binary")]] <- 
      colData(spe)[[paste0(qc_metric, "_3MAD_outlier_sample")]]
  }
  
  # Exclude off-tissue spots from outlier detection
  colData(spe)[[paste0(qc_metric, "_3MAD_outlier_binary")]] <- 
    ifelse(colData(spe)$in_tissue == FALSE, FALSE, 
           colData(spe)[[paste0(qc_metric, "_3MAD_outlier_binary")]])
  
  sampleList <- unique(colData(spe)[[samples]])
  names(sampleList) <- sampleList
  
  # Edge detection
  message("Detecting edges...")
  genes_edges <- lapply(sampleList, function(x) {
    tmp <- colData(spe)[colData(spe)[[samples]] == x, 
                        c("in_tissue", "array_row", "array_col", 
                          paste0(qc_metric, "_3MAD_outlier_binary"))]
    
    # CRITICAL FIX: Convert DataFrame to regular data.frame for raster operations
    tmp_df <- as.data.frame(tmp)
    
    # Ensure numeric columns are properly typed
    tmp_df$array_row <- as.numeric(tmp_df$array_row)
    tmp_df$array_col <- as.numeric(tmp_df$array_col)
    tmp_df[[paste0(qc_metric, "_3MAD_outlier_binary")]] <- as.logical(tmp_df[[paste0(qc_metric, "_3MAD_outlier_binary")]])
    
    result <- clumpEdges(tmp_df[, -1], rownames(tmp_df)[tmp_df$in_tissue == FALSE], 
                        shifted = shifted, edge_threshold = edge_threshold, min_cluster_size = min_cluster_size)
    return(result)
  })
  
  # Store edge results
  colData(spe)[[paste0(name, "_edge")]] <- FALSE
  edge_spots <- unlist(genes_edges)
  if (length(edge_spots) > 0) {
    colData(spe)[edge_spots, paste0(name, "_edge")] <- TRUE
  }
  
  message("Number of samples with edges detected: ", 
          sum(sapply(genes_edges, length) > 0))
  samples_with_edges <- names(genes_edges)[sapply(genes_edges, length) > 0]
  if (length(samples_with_edges) > 0) {
    message("Samples with edges detected: ", paste(samples_with_edges, collapse = ", "))
  }
  
  # Problem areas detection
  message("Finding problem areas...")
  genes_probs <- lapply(sampleList, function(x) {
    tmp <- colData(spe)[colData(spe)[[samples]] == x, 
                        c("in_tissue", "array_row", "array_col", 
                          paste0(qc_metric, "_3MAD_outlier_binary"))]
    
    # CRITICAL FIX: Convert DataFrame to regular data.frame for raster operations
    tmp_df <- as.data.frame(tmp)
    
    # Ensure numeric columns are properly typed
    tmp_df$array_row <- as.numeric(tmp_df$array_row)
    tmp_df$array_col <- as.numeric(tmp_df$array_col)
    tmp_df[[paste0(qc_metric, "_3MAD_outlier_binary")]] <- as.logical(tmp_df[[paste0(qc_metric, "_3MAD_outlier_binary")]])
    
    result <- problemAreas(tmp_df[, -1], 
                          rownames(tmp_df)[tmp_df$in_tissue == FALSE], 
                          uniqueIdentifier = x, 
                          shifted = shifted,
                          min_cluster_size = min_cluster_size)
    return(result)
  })
  
  # Store problem area results
  genes_probs <- do.call(rbind, genes_probs)
  colData(spe)[[paste0(name, "_problem_id")]] <- NA
  colData(spe)[[paste0(name, "_problem_size")]] <- 0
  
  if (nrow(genes_probs) > 0) {
    colData(spe)[genes_probs$spotcode, paste0(name, "_problem_id")] <- genes_probs$clumpID
    colData(spe)[genes_probs$spotcode, paste0(name, "_problem_size")] <- genes_probs$clumpSize
  }
  
  message("Edge dryspot detection completed!")
  return(spe)
}

#' Calculate edge zones for spots based on distance from tissue boundaries
#'
#' @param spe A SpatialExperiment object
#' @param samples Character string specifying the sample ID column name
#' @return Character vector of edge zone classifications ("1", "2", "3", "interior")
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom stats quantile
#' @keywords internal

calculateEdgeZones <- function(spe, samples) {
  edge_zones <- rep("interior", nrow(colData(spe)))

  # Process each sample separately
  sample_list <- unique(colData(spe)[[samples]])

  for (sample_id in sample_list) {
    sample_spots <- which(colData(spe)[[samples]] == sample_id & colData(spe)$in_tissue)

    if (length(sample_spots) > 0) {
      # Get coordinates for this sample
      coords <- spatialCoords(spe)[sample_spots, ]

      # Calculate distance from boundary for each spot
      # Simple approach: distance from min/max coordinates
      x_dist_from_edge <- pmin(coords[,1] - min(coords[,1]), max(coords[,1]) - coords[,1])
      y_dist_from_edge <- pmin(coords[,2] - min(coords[,2]), max(coords[,2]) - coords[,2])
      dist_from_edge <- pmin(x_dist_from_edge, y_dist_from_edge)

      # Classify by percentiles
      percentiles <- quantile(dist_from_edge, c(0.33, 0.66))

      edge_zones[sample_spots[dist_from_edge <= percentiles[1]]] <- "1"
      edge_zones[sample_spots[dist_from_edge <= percentiles[2] & dist_from_edge > percentiles[1]]] <- "2"
      edge_zones[sample_spots[dist_from_edge > percentiles[2]]] <- "3"
    }
  }

  return(edge_zones)
}

#' Classify sample position relative to array boundaries
#'
#' @param spe A SpatialExperiment object
#' @param samples Character string specifying the sample ID column name
#' @return Character vector of sample position classifications ("corner", "edge", "interior")
#' @importFrom SpatialExperiment spatialCoords
#' @keywords internal

calculateSamplePosition <- function(spe, samples) {
  sample_positions <- rep("interior", nrow(colData(spe)))

  # Process each sample separately
  sample_list <- unique(colData(spe)[[samples]])

  for (sample_id in sample_list) {
    sample_spots <- which(colData(spe)[[samples]] == sample_id & colData(spe)$in_tissue)

    if (length(sample_spots) > 0) {
      coords <- spatialCoords(spe)[sample_spots, ]

      # Check if sample touches array boundaries
      # Assuming standard Visium array coordinates
      touches_left <- min(coords[,1]) <= 1  # Close to minimum array position
      touches_right <- max(coords[,1]) >= 77  # Close to maximum array position
      touches_top <- min(coords[,2]) <= 1
      touches_bottom <- max(coords[,2]) >= 127

      open_sides <- sum(c(touches_left, touches_right, touches_top, touches_bottom))

      # Classify entire sample based on boundary contact
      if (open_sides >= 2) {
        sample_positions[sample_spots] <- "corner"
      } else if (open_sides == 1) {
        sample_positions[sample_spots] <- "edge"
      } else {
        sample_positions[sample_spots] <- "interior"
      }
    }
  }

  return(sample_positions)
}

#' Classify edge dryspots into categories for quality control (Enhanced Version)
#'
#' This function classifies detected edge dryspots and problem areas into 
#' different categories based on size and quality metrics for downstream 
#' filtering decisions. Enhanced with backward compatibility.
#'
#' @param spe A SpatialExperiment object that has been processed with detectEdgeDryspots()
#' @param samples Character string specifying the sample ID column name (default: "sample_id")
#' @param min_spots Minimum number of spots required for a problem area to be considered (default: 20)
#' @param low_umi_threshold UMI count threshold below which spots are considered low quality (default: 100)
#' @param removal_threshold Proportion threshold for marking problem areas for removal (default: 0.5)
#' @param name Character string matching the name used in detectEdgeDryspots() (default: "edge_dryspot")
#' @param exclude_slides Character vector of slide IDs to exclude from edge detection (default: NULL)
#'
#' @return A SpatialExperiment object with additional classification columns in colData:
#'   \item{[name]_true_edges}{Logical indicating edges after slide exclusion}
#'   \item{[name]_binary}{Character classification: "edge", "problem area", or "none"}
#'   \item{[name]_classification}{Character detailed classification: "edge", "remove", "flag", "small", or "none"}
#'   \item{[name]_enhanced}{Character enhanced classification: "edge_dryspot", "interior_outlier", etc.}
#'   \item{[name]_edge_zone}{Numeric/Character indicating edge zone: 1, 2, 3, or "interior"}
#'   \item{[name]_sample_type}{Character sample type: "corner", "edge", or "interior"}
#'
#' @examples
#' # Basic usage (after running detectEdgeDryspots)
#' \dontrun{
#' spe <- classifyEdgeDryspots(spe)
#' 
#' # View enhanced classifications
#' table(spe$edge_dryspot_enhanced)
#' table(spe$edge_dryspot_edge_zone)
#' }
#'
#' @export
classifyEdgeDryspots <- function(
    spe,
    samples = "sample_id",
    min_spots = 20,
    low_umi_threshold = 100,
    removal_threshold = 0.5,
    name = "edge_dryspot",
    exclude_slides = NULL) {

  # Input validation
  if (!inherits(spe, "SpatialExperiment")) {
    stop("Input data must be a SpatialExperiment or inherit from SpatialExperiment.")
  }
  
  required_cols <- c("sum_umi", paste0(name, "_edge"), paste0(name, "_problem_id"))
  missing_cols <- required_cols[!required_cols %in% colnames(colData(spe))]
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", "),
               "\nRun detectEdgeDryspots() first."))
  }
  
  message("Classifying edge dryspots with enhanced categories...")
  
  # Process edge classification
  colData(spe)[[paste0(name, "_true_edges")]] <- colData(spe)[[paste0(name, "_edge")]]
  
  # Exclude specified slides if provided
  if (!is.null(exclude_slides) && "slide" %in% colnames(colData(spe))) {
    colData(spe)[[paste0(name, "_true_edges")]] <- 
      ifelse(colData(spe)$slide %in% exclude_slides, FALSE, 
             colData(spe)[[paste0(name, "_edge")]])
    message("Excluding edges from slides: ", paste(exclude_slides, collapse = ", "))
  }

  # Identify problem areas for removal
  cdata <- as.data.frame(colData(spe))
  if (requireNamespace("dplyr", quietly = TRUE)) {
    tmp <- cdata %>%
      dplyr::filter(in_tissue == TRUE, 
                    !!sym(paste0(name, "_true_edges")) == FALSE) %>%
      dplyr::mutate(lowumi = sum_umi <= low_umi_threshold) %>%
      dplyr::group_by(!!sym(paste0(name, "_problem_id"))) %>%
      dplyr::summarise(n_lowumi = sum(lowumi), 
                       n_spots = dplyr::n(), 
                       prop_lowumi = n_lowumi/n_spots) %>%
      dplyr::filter(n_spots > min_spots, !is.na(!!sym(paste0(name, "_problem_id"))))
    
    remove.areas <- unique(dplyr::filter(tmp, prop_lowumi >= removal_threshold)[[paste0(name, "_problem_id")]])
  } else {
    # Fallback without dplyr
    cdata_filtered <- cdata[cdata$in_tissue == TRUE & 
                              cdata[[paste0(name, "_true_edges")]] == FALSE, ]
    cdata_filtered$lowumi <- cdata_filtered$sum_umi <= low_umi_threshold
    problem_ids <- unique(cdata_filtered[[paste0(name, "_problem_id")]])
    problem_ids <- problem_ids[!is.na(problem_ids)]
    
    remove.areas <- c()
    for (pid in problem_ids) {
      subset_data <- cdata_filtered[cdata_filtered[[paste0(name, "_problem_id")]] == pid, ]
      n_spots <- nrow(subset_data)
      n_lowumi <- sum(subset_data$lowumi)
      prop_lowumi <- n_lowumi / n_spots
      
      if (n_spots > min_spots && prop_lowumi >= removal_threshold) {
        remove.areas <- c(remove.areas, pid)
      }
    }
  }
  
  message("Number of problem areas marked for removal: ", length(remove.areas))
  
  # ============================================================================
  # BACKWARD COMPATIBLE CLASSIFICATIONS
  # ============================================================================
  
  # Create binary classification
  colData(spe)[[paste0(name, "_binary")]] <- 
    ifelse(!is.na(colData(spe)[[paste0(name, "_problem_id")]]), "problem area", "none")
  colData(spe)[colData(spe)[[paste0(name, "_edge")]], paste0(name, "_binary")] <- "edge"

  # Create detailed classification
  problem_sizes <- colData(spe)[[paste0(name, "_problem_size")]]
  colData(spe)[[paste0(name, "_classification")]] <- 
    ifelse(problem_sizes <= min_spots, "small", "flag")
  colData(spe)[[paste0(name, "_classification")]] <- 
    ifelse(problem_sizes == 0, "none", 
           colData(spe)[[paste0(name, "_classification")]])

  # Mark areas for removal
  if (length(remove.areas) > 0) {
    colData(spe)[colData(spe)[[paste0(name, "_problem_id")]] %in% remove.areas, 
                 paste0(name, "_classification")] <- "remove"
  }

  # Mark edges
  colData(spe)[colData(spe)[[paste0(name, "_true_edges")]], 
               paste0(name, "_classification")] <- "edge"
  
  # ============================================================================
  # ENHANCED CLASSIFICATIONS
  # ============================================================================
  
  # Initialize enhanced classification
  colData(spe)[[paste0(name, "_enhanced")]] <- "none"
 
  colData(spe)[[paste0(name, "_edge_zone")]] <- calculateEdgeZones(spe, samples)
 
  colData(spe)[[paste0(name, "_sample_type")]] <- calculateSamplePosition(spe, samples)
  
  # Enhanced categories with proper priority order
  # 1. HIGHEST PRIORITY: Edge dryspots (these are the actual edges detected)
  edge_spots <- colData(spe)[[paste0(name, "_true_edges")]]
  colData(spe)[edge_spots, paste0(name, "_enhanced")] <- "edge_dryspot"
  
  # 2. SECOND PRIORITY: Problem areas that should be removed (but not edges)
  if (length(remove.areas) > 0) {
    remove_spots <- colData(spe)[[paste0(name, "_problem_id")]] %in% remove.areas & 
                    !colData(spe)[[paste0(name, "_true_edges")]]  # Exclude edge spots
    colData(spe)[remove_spots, paste0(name, "_enhanced")] <- "edge_sample_artifact"
  }
  
  # 3. THIRD PRIORITY: Small clusters (but not edges or removal areas)
  small_clusters <- colData(spe)[[paste0(name, "_problem_size")]] > 0 & 
                    colData(spe)[[paste0(name, "_problem_size")]] <= min_spots &
                    !colData(spe)[[paste0(name, "_true_edges")]] &  # Exclude edge spots
                    !(colData(spe)[[paste0(name, "_problem_id")]] %in% remove.areas)  # Exclude removal areas
  colData(spe)[small_clusters, paste0(name, "_enhanced")] <- "small_cluster"
  
  message("Enhanced edge dryspot classification completed!")
  message("Classifications added:")
  message("  - Backward compatible: ", paste0(name, "_binary"), ", ", paste0(name, "_classification"))
  message("  - Enhanced: ", paste0(name, "_enhanced"), ", ", paste0(name, "_edge_zone"), ", ", paste0(name, "_sample_type"))
  message("  - Edge zones and sample types will be implemented in subsequent phases")
  
  return(spe)
}
