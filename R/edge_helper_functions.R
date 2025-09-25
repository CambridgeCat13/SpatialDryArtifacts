
#' Fill isolated spots in 3x3 focal window
#'
#' Internal focal function that fills a center pixel if it is completely 
#' surrounded by outlier spots in a 3x3 window.
#'
#' @param x Numeric vector of length 9 representing a 3x3 window where x[5] is the center
#'
#' @return Integer value: 1 if center should be filled (marked as outlier), 
#'   0 if it should remain unchanged, NA_integer_ if center is NA
#'
#' @details This function is used internally by focal_transformations to fill
#'   gaps in outlier regions. It implements the first morphological operation
#'   to connect nearby outlier spots.
#'
#' @keywords internal
#' @noRd
my_fill <- function(x) {
  if(is.na(x[5])) return(NA_integer_)
  if(x[5]==1) return(1)
  x[is.na(x)] = 1
  if(sum(x[-5])==8) return(1)
  else return(0)
}

#' Fill spots with outliers in cardinal directions
#'
#' Internal focal function that fills a center pixel if all four cardinal 
#' neighbors (North, South, East, West) are outliers.
#'
#' @param x Numeric vector of length 9 representing a 3x3 window in star pattern
#'
#' @return Integer value: 1 if center should be filled, 0 otherwise, 
#'   NA_integer_ if center is NA
#'
#' @details Implements a star-shaped focal pattern where only cardinal directions
#'   are considered. This helps connect outlier regions along major axes while
#'   preserving diagonal boundaries.
#'
#' @keywords internal
#' @noRd
my_fill_star <- function(x) {
  if(is.na(x[5])) return(NA_integer_)
  if(x[5]==1) return(1)
  if(sum(x[-5], na.rm=T)==4) return(1)
  else return(0)
}

#' Fill spots outlined by outliers in 5x5 window
#'
#' Internal focal function that fills a center pixel if it is completely 
#' outlined by outlier spots in a 5x5 window perimeter.
#'
#' @param x Numeric vector of length 25 representing a 5x5 window where x[13] is the center
#'
#' @return Integer value: 1 if center should be filled, 0 otherwise,
#'   NA_integer_ if center is NA
#'
#' @details Checks if all 16 border pixels of a 5x5 window are outliers.
#'   This captures larger-scale patterns and fills interior gaps within
#'   outlier regions.
#'
#' @keywords internal
#' @noRd
my_outline <- function(x) {
  if(is.na(x[13])) return(NA_integer_)
  if(x[13]==1) return(1)
  x[is.na(x)] = 1
  sum.edges = sum(x[1:5])+sum(x[21:25])+sum(x[c(6,11,16)])+sum(x[c(10,15,20)])
  if(sum.edges==16) return(1)
  else(return(0))
}

#' Apply morphological transformations to identify connected outlier regions
#'
#' Performs a series of focal operations to clean and connect outlier regions
#' in spatial transcriptomics data through morphological operations.
#'
#' @param raster_object A raster object with binary values where 1 indicates 
#'   outlier spots and 0/NA indicates normal spots
#' @param min_cluster_size Numeric value specifying the minimum size for 
#'   isolated clusters. Clusters smaller than this threshold will be filled 
#'   (default: 40)
#'
#' @return A processed raster object with cleaned and connected outlier regions
#'
#' @details The function applies four sequential morphological operations:
#'   \enumerate{
#'     \item 3x3 fill: Fills spots completely surrounded by outliers
#'     \item 5x5 outline: Fills spots outlined by outliers in larger window
#'     \item Star pattern: Fills spots with outliers in cardinal directions
#'     \item Small cluster removal: Removes isolated normal regions below threshold
#'   }
#'
#' @examples
#' \dontrun{
#' # Create a raster from outlier data
#' r <- raster::rasterFromXYZ(outlier_coords)
#' # Apply morphological cleaning
#' r_cleaned <- focal_transformations(r, min_cluster_size = 30)
#' }
#'
#' @importFrom raster focal clump values values<-
#' @export
focal_transformations <- function(raster_object, min_cluster_size = 40) {
  r2 <- raster::focal(raster_object, w=matrix(1,3,3), fun=my_fill, pad=TRUE, padValue=1)
  r3 <- raster::focal(r2, w=matrix(1,5,5), fun=my_outline, pad=TRUE, padValue=1)
  
  star.m = matrix(rep(c(0,1), length.out=9), 3, 3)
  star.m[5] = 1
  r3_s <- raster::focal(r3, star.m, fun=my_fill_star, pad=TRUE, padValue=1)
  
  rev_r3 = r3_s
  rev_r3[r3_s==0] = 1
  rev_r3[r3_s==1] = 0
  rev_c3 = raster::clump(rev_r3)
  
  rev_c3_values <- raster::values(rev_c3)
  if (is.null(rev_c3_values)) {
    return(r3_s)
  }
  
  tbl = table(rev_c3_values[!is.na(rev_c3_values)])
  
  if (length(tbl) > 0) {
    flip_clump = as.numeric(names(tbl)[tbl < min_cluster_size])
    r4 = r3_s
    
    if (length(flip_clump) > 0) {
      r4_values <- raster::values(r4)
      r4_values[rev_c3_values %in% flip_clump] <- 1
      raster::values(r4) <- r4_values
    }
    
    return(r4)
  } else {
    return(r3_s)
  }
}

#' Map raster matrix indices to spot identifiers
#'
#' Internal function that maps matrix row/column indices from raster operations
#' back to original spot identifiers.
#'
#' @param .xy Data frame containing original spot coordinates with row names 
#'   as spot identifiers. Must have columns for array coordinates.
#' @param results_df Data frame with matrix indices from raster operations,
#'   typically containing row and column indices of detected spots
#'
#' @return Character vector of spot identifiers corresponding to the matrix indices
#'
#' @details This function creates a mapping between the raster matrix representation
#'   and the original spot identifiers, enabling translation of raster analysis
#'   results back to the original coordinate system.
#'
#' @keywords internal
#' @noRd
lookupKey <- function(.xy, results_df) {
  if (nrow(results_df) == 0) {
    return(character(0))
  }
  
  key1 = cbind(.xy, "numeric_key"=seq(nrow(.xy)))
  t1 = raster::rasterFromXYZ(key1)
  rast_values <- raster::values(t1)
  rast_dims <- c(nrow(t1), ncol(t1))
  rast <- matrix(rast_values, nrow = rast_dims[1], ncol = rast_dims[2], byrow = TRUE)
  
  num_key = sapply(1:nrow(results_df), function(x) {
    v1 = as.numeric(results_df[x,])
    if (v1[1] <= nrow(rast) && v1[2] <= ncol(rast)) {
      rast[v1[1], v1[2]]
    } else {
      NA
    }
  })
  
  num_key <- num_key[!is.na(num_key)]
  edge_spots = rownames(key1)[key1$numeric_key %in% num_key]
  
  return(edge_spots)
}

#' Map raster clusters to spot identifiers with metadata
#'
#' Internal function that maps raster cluster information including cluster IDs
#' and sizes back to original spot identifiers.
#'
#' @param .xy Data frame containing original spot coordinates with row names
#'   as spot identifiers
#' @param results_df Data frame with columns for row index, column index,
#'   cluster ID, and cluster size from raster clustering operations
#'
#' @return Data frame with columns:
#'   \itemize{
#'     \item spotcode: Original spot identifier
#'     \item clumpID: Cluster identifier
#'     \item clumpSize: Number of spots in the cluster
#'   }
#'
#' @details Extended version of lookupKey that preserves cluster metadata
#'   for problem area identification.
#'
#' @keywords internal
#' @noRd
lookupKeyDF <- function(.xy, results_df) {
  if (nrow(results_df) == 0) {
    return(data.frame(spotcode = character(0), 
                     clumpID = character(0), 
                     clumpSize = numeric(0)))
  }
  
  key1 = cbind(.xy, "numeric_key"=seq(nrow(.xy)))
  t1 = raster::rasterFromXYZ(key1)
  rast_values <- raster::values(t1)
  rast_dims <- c(nrow(t1), ncol(t1))
  rast <- matrix(rast_values, nrow = rast_dims[1], ncol = rast_dims[2], byrow = TRUE)
  
  num_key = do.call(rbind, lapply(1:nrow(results_df), function(x) {
    if (results_df[x,1] <= nrow(rast) && results_df[x,2] <= ncol(rast)) {
      cbind.data.frame(
        "numeric_key"=rast[results_df[x,1], results_df[x,2]], 
        "clump_id"=results_df[x,3], 
        "size"=results_df[x,4]
      )
    } else {
      NULL
    }
  }))
  
  if (is.null(num_key) || nrow(num_key) == 0) {
    return(data.frame(spotcode = character(0), 
                     clumpID = character(0), 
                     clumpSize = numeric(0)))
  }
  
  matching = match(num_key$numeric_key, key1$numeric_key)
  matching <- matching[!is.na(matching)]
  
  if (length(matching) == 0) {
    return(data.frame(spotcode = character(0), 
                     clumpID = character(0), 
                     clumpSize = numeric(0)))
  }
  
  return(cbind.data.frame(
    "spotcode"=rownames(key1)[matching],
    "clumpID"=num_key$clump_id[!is.na(matching)], 
    "clumpSize"=num_key$size[!is.na(matching)]
  ))
}


#' Detect edge dryspots using spatial clustering
#'
#' Identifies clusters of low-quality spots that touch tissue boundaries,
#' indicating potential edge dryspot artifacts from incomplete reagent coverage.
#'
#' @param .xyz Data frame with spot coordinates and QC metrics. Must contain:
#'   \itemize{
#'     \item Column 1: array_row coordinates
#'     \item Column 2: array_col coordinates
#'     \item Column 3: Binary outlier indicator (TRUE/FALSE or 1/0)
#'     \item Row names: Spot identifiers
#'   }
#' @param offTissue Character vector of spot identifiers that are off-tissue
#' @param shifted Logical indicating whether to apply coordinate adjustment 
#'   for hexagonal arrays (default: FALSE)
#' @param edge_threshold Numeric value between 0 and 1 specifying the minimum
#'   proportion of border coverage required for edge detection (default: 0.75)
#' @param min_cluster_size Numeric value for minimum cluster size in 
#'   morphological operations (default: 40)
#'
#' @return Character vector of spot identifiers classified as edge dryspots
#'
#' @details The function performs the following steps:
#'   \enumerate{
#'     \item Converts spot coordinates to raster format
#'     \item Applies morphological transformations to connect outlier regions
#'     \item Identifies connected components (clusters)
#'     \item Checks each tissue border (N, S, E, W) for cluster coverage
#'     \item Returns spots from clusters that exceed edge_threshold on any border
#'   }
#'
#' @examples
#' \dontrun{
#' # Prepare coordinate data with outliers
#' spot_data <- data.frame(
#'   array_row = spe$array_row,
#'   array_col = spe$array_col,
#'   outlier = spe$low_quality
#' )
#' rownames(spot_data) <- colnames(spe)
#' 
#' # Detect edge dryspots
#' edge_spots <- clumpEdges(
#'   spot_data,
#'   offTissue = rownames(spot_data)[!spe$in_tissue],
#'   edge_threshold = 0.75
#' )
#' }
#'
#' @importFrom raster rasterFromXYZ focal clump values
#' @export
clumpEdges <- function(.xyz, offTissue, shifted=FALSE, edge_threshold=0.75, min_cluster_size=40) {
  # Handle NA values
  if(ncol(.xyz) >= 3) {
    .xyz[, 3][is.na(.xyz[, 3])] <- FALSE
  }
  
  if(sum(.xyz[,3])==0) {
    return(c())
  }
  
  if(shifted==TRUE) {
    odds = seq(1, max(.xyz[,"array_col"]), by=2)
    .xyz[.xyz[,"array_col"] %in% odds, "array_col"] = 
      .xyz[.xyz[,"array_col"] %in% odds, "array_col"]-1
  }
  
  t1 = raster::rasterFromXYZ(.xyz)
  t2 <- focal_transformations(t1, min_cluster_size = min_cluster_size)
  rast_values = raster::values(t2)
  rast_dims = c(nrow(t2), ncol(t2))
  rast = matrix(rast_values, nrow = rast_dims[1], ncol = rast_dims[2], byrow = TRUE)
  
  c1 = raster::clump(t2, directions=8)
  
  clump_values = raster::values(c1)
  clump_dims = c(nrow(c1), ncol(c1))
  clumps = matrix(clump_values, nrow = clump_dims[1], ncol = clump_dims[2], byrow = TRUE)
  
  edgeClumps = c()

  north_total = sum(!is.na(rast[1,]))
  if (north_total > 0) {
    north = sum(!is.na(clumps[1,]))/north_total
    if(!is.na(north) && north >= edge_threshold) {
      edgeClumps = c(edgeClumps, unique(clumps[1,]))
    }
  }
  
  east_total = sum(!is.na(rast[,ncol(rast)]))
  if (east_total > 0) {
    east = sum(!is.na(clumps[,ncol(clumps)]))/east_total
    if(!is.na(east) && east >= edge_threshold) {
      edgeClumps = c(edgeClumps, unique(clumps[,ncol(clumps)]))
    }
  }
  
  south_total = sum(!is.na(rast[nrow(rast),]))
  if (south_total > 0) {
    south = sum(!is.na(clumps[nrow(clumps),]))/south_total
    if(!is.na(south) && south >= edge_threshold) {
      edgeClumps = c(edgeClumps, unique(clumps[nrow(clumps),]))
    }
  }
  
  west_total = sum(!is.na(rast[,1]))
  if (west_total > 0) {
    west = sum(!is.na(clumps[,1]))/west_total
    if(!is.na(west) && west >= edge_threshold) {
      edgeClumps = c(edgeClumps, unique(clumps[,1]))
    }
  }
  
  edgeClumps = unique(edgeClumps[!is.na(edgeClumps)])
  
  if(length(edgeClumps)==0) {
    return(c())
  }
  
  res <- vector("list", length(edgeClumps))
  names(res) <- as.character(edgeClumps)
  
  for (i in edgeClumps){
    cluster_indices <- which(clumps == i, arr.ind = TRUE)
    if (nrow(cluster_indices) > 0) {
      res[[as.character(i)]] <- as.data.frame(cluster_indices)
    }
  }
  
  res <- res[!sapply(res, is.null)]
  
  if (length(res) == 0) {
    return(c())
  }
  
  res.df = do.call(rbind, res)
  edgeSpots = lookupKey(.xyz[,1:2], res.df)
  result <- setdiff(edgeSpots, offTissue)
  
  return(result)
}


#' Identify all problem areas in spatial transcriptomics data
#'
#' Detects and characterizes all clusters of low-quality spots (problem areas)
#' in tissue sections, including both edge and interior artifacts.
#'
#' @param .xyz Data frame with spot coordinates and QC metrics. Must contain:
#'   \itemize{
#'     \item Column 1: array_row coordinates
#'     \item Column 2: array_col coordinates
#'     \item Column 3: Binary outlier indicator (TRUE/FALSE or 1/0)
#'     \item Row names: Spot identifiers
#'   }
#' @param offTissue Character vector of spot identifiers that are off-tissue
#' @param uniqueIdentifier Character string used as prefix for cluster IDs.
#'   If NA, "X" will be used (default: NA)
#' @param shifted Logical indicating whether to apply coordinate adjustment
#'   for hexagonal arrays (default: FALSE)
#' @param min_cluster_size Numeric value for minimum cluster size in
#'   morphological operations (default: 40)
#'
#' @return Data frame with the following columns:
#'   \itemize{
#'     \item spotcode: Spot identifier
#'     \item clumpID: Unique cluster identifier (format: "prefix_number")
#'     \item clumpSize: Number of spots in the cluster
#'   }
#'   Returns empty data frame if no problem areas are found.
#'
#' @details This function identifies ALL connected components of outlier spots,
#'   not just those touching edges. Each cluster is assigned a unique ID and
#'   its size is calculated. This enables downstream filtering based on
#'   cluster characteristics.
#'
#' @examples
#' \dontrun{
#' # Identify all problem areas in a sample
#' problem_df <- problemAreas(
#'   spot_data,
#'   offTissue = rownames(spot_data)[!spe$in_tissue],
#'   uniqueIdentifier = "Sample1"
#' )
#' 
#' # Filter for large problem areas
#' large_problems <- problem_df[problem_df$clumpSize > 50, ]
#' 
#' # Get spots in specific cluster
#' cluster1_spots <- problem_df[problem_df$clumpID == "Sample1_1", "spotcode"]
#' }
#'
#' @importFrom raster rasterFromXYZ clump values
#' @export
problemAreas <- function(.xyz, offTissue, uniqueIdentifier=NA, shifted=FALSE, min_cluster_size=40) {
  # Handle NA values
  if(ncol(.xyz) >= 3) {
    .xyz[, 3][is.na(.xyz[, 3])] <- FALSE
  }
  
  if(sum(.xyz[,3])==0) {
    return(data.frame(spotcode = character(0), 
                     clumpID = character(0), 
                     clumpSize = numeric(0)))
  }

  if(shifted==TRUE) {
    odds = seq(1, max(.xyz[,"array_col"]), by=2)
    .xyz[.xyz[,"array_col"] %in% odds, "array_col"] = 
      .xyz[.xyz[,"array_col"] %in% odds, "array_col"]-1
  }
  
  t1 = raster::rasterFromXYZ(.xyz)
  t2 = focal_transformations(t1, min_cluster_size = min_cluster_size)
  
  c1 = raster::clump(t2, directions=8)
  
  clump_values = raster::values(c1)
  clump_values <- clump_values[!is.na(clump_values)]
  
  if (length(clump_values) == 0) {
    return(data.frame(spotcode = character(0), 
                     clumpID = character(0), 
                     clumpSize = numeric(0)))
  }
  
  clump_dims = c(nrow(c1), ncol(c1))
  clumps = matrix(raster::values(c1), nrow = clump_dims[1], ncol = clump_dims[2], byrow = TRUE)
  
  tot <- max(clumps, na.rm=TRUE)
  
  if (is.na(tot) || tot == 0) {
    return(data.frame(spotcode = character(0), 
                     clumpID = character(0), 
                     clumpSize = numeric(0)))
  }
  
  res <- vector("list", tot)
  if(is.na(uniqueIdentifier)) uniqueIdentifier = "X"
  
  for (i in 1:tot){
    cluster_indices <- which(clumps == i, arr.ind = TRUE)
    if (nrow(cluster_indices) > 0) {
      res[[i]] <- cbind.data.frame(
        cluster_indices,
        "clump_id"=paste(uniqueIdentifier, i, sep="_"),
        "size"=nrow(cluster_indices)
      )
    }
  }
  
  res <- res[!sapply(res, is.null)]
  
  if (length(res) == 0) {
    return(data.frame(spotcode = character(0), 
                     clumpID = character(0), 
                     clumpSize = numeric(0)))
  }
  
  res.df = do.call(rbind.data.frame, res)
  pAreas = lookupKeyDF(.xyz[,1:2], res.df)
  
  if (nrow(pAreas) == 0) {
    return(data.frame(spotcode = character(0), 
                     clumpID = character(0), 
                     clumpSize = numeric(0)))
  }
  
  result <- pAreas[!pAreas$spotcode %in% offTissue,]
  
  return(result)
}
