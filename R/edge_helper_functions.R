# ==============================================================================
# Enhanced Edge Detection Helper Functions
# Phase 1: Parameterized version of Jacqui's functions WITH COORDINATE FIX
# ==============================================================================

my_fill <- function(x) {
  if(is.na(x[5])) return(NA_integer_)
  if(x[5]==1) return(1)
  x[is.na(x)] = 1
  if(sum(x[-5])==8) return(1)
  else return(0)
}

my_fill_star <- function(x) {
  if(is.na(x[5])) return(NA_integer_)
  if(x[5]==1) return(1)
  if(sum(x[-5], na.rm=T)==4) return(1)
  else return(0)
}

my_outline <- function(x) {
  if(is.na(x[13])) return(NA_integer_)
  if(x[13]==1) return(1)
  x[is.na(x)] = 1
  sum.edges = sum(x[1:5])+sum(x[21:25])+sum(x[c(6,11,16)])+sum(x[c(10,15,20)])
  if(sum.edges==16) return(1)
  else(return(0))
}

#' Enhanced focal transformations with parameterized cluster size
#' 
#' @param raster_object Raster object for processing
#' @param min_cluster_size Minimum cluster size threshold (default: 40)
#' @return Processed raster object
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
  tbl = table(rev_c3_values, useNA = "no")
  
  flip_clump = as.numeric(names(tbl)[tbl < min_cluster_size])
  r4 = r3_s
  
  if (length(flip_clump) > 0) {
    r4_values <- raster::values(r4)
    r4_values[rev_c3_values %in% flip_clump] <- 1
    raster::values(r4) <- r4_values
  }
  
  return(r4)
}

lookupKey <- function(.xy, results_df) {
  if (nrow(results_df) == 0) {
    return(character(0))
  }
  
  result_spots <- c()
  
  unique_rows <- sort(unique(.xy[,1]))  # Original array_row values
  unique_cols <- sort(unique(.xy[,2]))  # Original array_col values
  
  for(i in 1:nrow(results_df)) {
    mat_row <- results_df[i,1]  # Matrix row index
    mat_col <- results_df[i,2]  # Matrix col index
    
    if (mat_col <= length(unique_rows) && mat_row <= length(unique_cols)) {
      orig_row_coord <- unique_rows[mat_col]
      orig_col_coord <- unique_cols[length(unique_cols) - mat_row + 1]
      
      match_idx <- which(abs(.xy[,1] - orig_row_coord) < 0.1 & 
                         abs(.xy[,2] - orig_col_coord) < 0.1)
      if(length(match_idx) > 0) {
        result_spots <- c(result_spots, rownames(.xy)[match_idx])
      }
    }
  }
  
  return(unique(result_spots))
}


lookupKeyDF <- function(.xy, results_df) {
  if (nrow(results_df) == 0) {
    return(data.frame(spotcode = character(0), 
                     clumpID = character(0), 
                     clumpSize = numeric(0)))
  }
  
  result_list <- list()

  unique_rows <- sort(unique(.xy[,1]))  # Original array_row values
  unique_cols <- sort(unique(.xy[,2]))  # Original array_col values
  
  for(i in 1:nrow(results_df)) {
    mat_row <- results_df[i,1]  # Matrix row index
    mat_col <- results_df[i,2]  # Matrix col index
    
    if (mat_col <= length(unique_rows) && mat_row <= length(unique_cols)) {
      orig_row_coord <- unique_rows[mat_col]
      orig_col_coord <- unique_cols[length(unique_cols) - mat_row + 1]
      
      match_idx <- which(abs(.xy[,1] - orig_row_coord) < 0.1 & 
                         abs(.xy[,2] - orig_col_coord) < 0.1)
      if(length(match_idx) > 0) {
        result_list[[i]] <- data.frame(
          spotcode = rownames(.xy)[match_idx],
          clumpID = results_df[i,3],
          clumpSize = results_df[i,4],
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  if (length(result_list) == 0) {
    return(data.frame(spotcode = character(0), 
                     clumpID = character(0), 
                     clumpSize = numeric(0)))
  }
  
  return(do.call(rbind, result_list))
}

#' Enhanced clump edges detection with parameterized thresholds AND COORDINATE FIX
#' 
#' @param .xyz Coordinate matrix with QC data
#' @param offTissue Vector of off-tissue spot names
#' @param shifted Whether to apply coordinate adjustment (default: FALSE)
#' @param edge_threshold Edge coverage threshold (default: 0.6, was 0.75)
#' @param min_cluster_size Minimum cluster size for morphological cleaning (default: 40)
#' @return Vector of edge spot names
clumpEdges <- function(.xyz, offTissue, shifted=FALSE, edge_threshold=0.6, min_cluster_size=40) {
  if(ncol(.xyz) >= 3) {
    .xyz[, 3][is.na(.xyz[, 3])] <- FALSE
  }
  
  if(sum(.xyz[,3])==0) {
    return(c())
  }

  if(shifted==TRUE) {
    odds = seq(1,max(.xyz[,"array_col"]), by=2)
    .xyz[.xyz[,"array_col"] %in% odds, "array_col"] = 
      .xyz[.xyz[,"array_col"] %in% odds, "array_col"]-1
  }
  
  xyz_corrected <- .xyz
  xyz_corrected[, c(1, 2)] <- .xyz[, c(2, 1)]
  
  t1 = raster::rasterFromXYZ(xyz_corrected)
  t2 <- focal_transformations(t1, min_cluster_size = min_cluster_size)
  rast_values = raster::values(t2)
  rast_dims = c(nrow(t2), ncol(t2))
  rast = matrix(rast_values, nrow = rast_dims[1], ncol = rast_dims[2])
  c1 = raster::clump(t2, directions=8)
  
  clump_values = raster::values(c1)
  clump_dims = c(nrow(c1), ncol(c1))
  clumps = matrix(clump_values, nrow = clump_dims[1], ncol = clump_dims[2])
  edgeClumps = c()
  
  north = sum(!is.na(clumps[1,]))/sum(!is.na(rast[1,]))
  if(!is.na(north) && north>=edge_threshold) edgeClumps = c(edgeClumps, unique(clumps[1,]))
  
  east = sum(!is.na(clumps[,ncol(clumps)]))/sum(!is.na(rast[,ncol(rast)]))
  if(!is.na(east) && east>=edge_threshold) edgeClumps = c(edgeClumps, unique(clumps[,ncol(clumps)]))
  
  south = sum(!is.na(clumps[nrow(clumps),]))/sum(!is.na(rast[nrow(rast),]))
  if(!is.na(south) && south>=edge_threshold) edgeClumps = c(edgeClumps, unique(clumps[nrow(clumps),]))
  
  west = sum(!is.na(clumps[,1]))/sum(!is.na(rast[,1]))
  if(!is.na(west) && west>=edge_threshold) edgeClumps = c(edgeClumps, unique(clumps[,1]))
  
  edgeClumps = edgeClumps[!is.na(edgeClumps)]
  if(length(edgeClumps)==0) {
    return(c())
  }
  
  res <- vector("list", length(edgeClumps))
  names(res) <- as.character(edgeClumps)
  for (i in edgeClumps){
    res[[as.character(i)]] <- as.data.frame(which(clumps == i, arr.ind = TRUE))
  }
  res.df = do.call(rbind, res) 
  
  edgeSpots = lookupKey(.xyz[,1:2], res.df)
  result <- setdiff(edgeSpots, offTissue)
  
  return(result)
}

#' Enhanced problem areas detection with parameterized cluster size AND COORDINATE FIX
#' 
#' @param .xyz Coordinate matrix with QC data
#' @param offTissue Vector of off-tissue spot names
#' @param uniqueIdentifier Sample identifier (default: NA)
#' @param shifted Whether to apply coordinate adjustment (default: FALSE)
#' @param min_cluster_size Minimum cluster size for morphological cleaning (default: 40)
#' @return Data frame with problem area information
problemAreas <- function(.xyz, offTissue, uniqueIdentifier=NA, shifted=FALSE, min_cluster_size=40) { 
  if(ncol(.xyz) >= 3) {
    .xyz[, 3][is.na(.xyz[, 3])] <- FALSE
  }
  
  if(sum(.xyz[,3])==0) {
    return(data.frame())
  }
  
  if(shifted==TRUE) {
    odds = seq(1,max(.xyz[,"array_col"]), by=2)
    .xyz[.xyz[,"array_col"] %in% odds, "array_col"] = 
      .xyz[.xyz[,"array_col"] %in% odds, "array_col"]-1
  }
  
  # CRITICAL FIX: Swap coordinates for rasterFromXYZ
  xyz_corrected <- .xyz
  xyz_corrected[, c(1, 2)] <- .xyz[, c(2, 1)]  # Swap row and col
  
  t1 = raster::rasterFromXYZ(xyz_corrected)
  t2 = focal_transformations(t1, min_cluster_size = min_cluster_size)
  rast_values = raster::values(t2)
  rast_dims = c(nrow(t2), ncol(t2))
  rast = matrix(rast_values, nrow = rast_dims[1], ncol = rast_dims[2])
  c1 = raster::clump(t2, directions=8)
  
  clump_values = raster::values(c1)
  clump_dims = c(nrow(c1), ncol(c1))
  clumps = matrix(clump_values, nrow = clump_dims[1], ncol = clump_dims[2])
  tot <- max(clumps, na.rm=TRUE)
  res <- vector("list",tot)
  if(is.na(uniqueIdentifier)) uniqueIdentifier = "X"
  
  for (i in 1:tot){
    res[i] <- list(which(clumps == i, arr.ind = TRUE))
    res[[i]] <- cbind.data.frame(
      res[[i]],
      "clump_id"=paste(uniqueIdentifier, i, sep="_"),
      "size"=nrow(res[[i]])
    )
  }
  res.df = do.call(rbind.data.frame, res) 
  pAreas = lookupKeyDF(.xyz[,1:2], res.df)
  result <- pAreas[!pAreas$spotcode %in% offTissue,]
  
  return(result)
}
