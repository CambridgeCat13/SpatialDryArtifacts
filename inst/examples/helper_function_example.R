library(raster)
#source("R/edge_helper_functions.R")

demo_focal_functions <- function() {
  cat("=== Focal Functions Demo ===\n")
  cat("my_fill function: fills 3x3 neighborhood\n")
  test_pattern <- c(1, 1, 0, 1, 0, 1, 0, 1, 1)  # center=0, surrounding 8 cells=1
  result <- my_fill(test_pattern)
  cat("Input pattern (center=0): [1,1,0; 1,0,1; 0,1,1] -> Result:", result, "\n\n")
  cat("my_fill_star function: star pattern filling\n")
  star_pattern <- c(0, 1, 0, 1, 0, 1, 0, 1, 0)  # cross pattern
  result <- my_fill_star(star_pattern)
  cat("Star pattern: [0,1,0; 1,0,1; 0,1,0] -> Result:", result, "\n\n")
}

demo_edge_detection <- function() {
  cat("=== Edge Detection Demo ===\n")
  coords <- expand.grid(array_row = 1:5, array_col = 1:5)
  rownames(coords) <- paste0("spot_", 1:25)
  edge_data <- cbind(coords, outlier = FALSE)
  edge_data[coords$array_col == 1, "outlier"] <- TRUE  # left edge is outliers
  cat("Input pattern (left edge as outliers):\n")
  pattern_matrix <- matrix(edge_data$outlier, nrow = 5, ncol = 5, byrow = FALSE)
  print(pattern_matrix)
  # Detect edges
  edges <- clumpEdges(edge_data, character(0), 
                     edge_threshold = 0.6, min_cluster_size = 5)
  cat("\nNumber of detected edge spots:", length(edges), "\n")
  if(length(edges) > 0) {
    cat("Edge spots:", head(edges, 5), "...\n")
  }
}

demo_problem_areas <- function() {
  cat("\n=== Problem Areas Detection Demo ===\n")
  coords <- expand.grid(array_row = 1:6, array_col = 1:6)
  rownames(coords) <- paste0("spot_", 1:36)
  problem_data <- cbind(coords, outlier = FALSE)
  center_spots <- which(coords$array_row %in% 2:3 & coords$array_col %in% 2:3)
  problem_data[center_spots, "outlier"] <- TRUE
  
  cat("Input pattern (central area as outliers):\n")
  pattern_matrix <- matrix(problem_data$outlier, nrow = 6, ncol = 6, byrow = FALSE)
  print(pattern_matrix)
  problems <- problemAreas(problem_data, character(0), 
                          uniqueIdentifier = "sample1", min_cluster_size = 2)
  cat("\nDetected problem areas:\n")
  if(nrow(problems) > 0) {
    print(head(problems))
  } else {
    cat("No problem areas detected\n")
  }
}

demo_coordinate_handling <- function() {
  cat("\n=== Coordinate Handling Demo ===\n")
  simple_coords <- data.frame(
    array_row = c(1, 1, 2, 2),
    array_col = c(1, 2, 1, 2)
  )
  rownames(simple_coords) <- c("A", "B", "C", "D")
  cat("Original coordinates:\n")
  print(simple_coords)
  test_data <- cbind(simple_coords, z = 1:4)
  test_raster <- raster::rasterFromXYZ(test_data[,c(2,1,3)])  
  cat("\nRaster matrix:\n")
  print(as.matrix(test_raster))
  lookup_positions <- data.frame(row = c(1, 2), col = c(1, 2))
  found_spots <- lookupKey(simple_coords, lookup_positions)
  cat("\nLooking up positions (1,1) and (2,2), found spots:", found_spots, "\n")
}


demo_parameter_effects <- function() {
  cat("\n=== Parameter Effects Demo ===\n")
  coords <- expand.grid(array_row = 1:8, array_col = 1:8)
  rownames(coords) <- paste0("spot_", 1:64)
  l_data <- cbind(coords, outlier = FALSE)
  l_spots <- which(coords$array_row == 1 | coords$array_col == 1)
  l_data[l_spots, "outlier"] <- TRUE
  cat("L-shaped outlier pattern:\n")
  l_matrix <- matrix(l_data$outlier, nrow = 8, ncol = 8, byrow = FALSE)
  print(l_matrix)
  cat("\nEffect of different edge_threshold values:\n")
  thresholds <- c(0.3, 0.6, 0.9)
  for(thresh in thresholds) {
    edges <- clumpEdges(l_data, character(0), 
                       edge_threshold = thresh, min_cluster_size = 5)
    cat(sprintf("threshold=%.1f: %d edges\n", thresh, length(edges)))
  }
  cat("\nEffect of different min_cluster_size values:\n")
  cluster_sizes <- c(5, 15, 30)
  for(size in cluster_sizes) {
    edges <- clumpEdges(l_data, character(0), 
                       edge_threshold = 0.6, min_cluster_size = size)
    cat(sprintf("min_size=%d: %d edges\n", size, length(edges)))
  }
}


demo_focal_functions()
demo_edge_detection()
demo_problem_areas()
demo_coordinate_handling()
demo_parameter_effects()

cat("\nNote: All thresholds return the same result because the L-shaped pattern creates a single large cluster that touches both the top (100% coverage) 
  and left (100% coverage) edges. Since all edges have 100% coverage, any threshold from 0.3 to 0.9 will detect this cluster as an edge.\n")

cat("\nNote: All cluster sizes return the same result because after focal transformations, the pattern creates a cluster larger than 30 spots,
      so none of the tested thresholds trigger the small cluster removal. To see parameter effects, try patterns with more subtle edge coverage or smaller clusters.\n")

cat("These helper functions are core components of edge dryspot detection:
    - focal functions: handle spatial morphological transformations
    - clumpEdges: detect edge dryspots 
    - problemAreas: identify internal problem regions
    - lookupKey: coordinate mapping and lookup
    - parameters are adjustable to meet different detection sensitivity requirements\n")