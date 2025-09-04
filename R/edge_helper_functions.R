# ==============================================================================
# Enhanced Edge Detection Helper Functions
# Phase 1: Parameterized version of Jacqui's functions
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
  r2 <- raster::focal(raster_object, w=matrix(1,3,3), fun=my_fill, pad=T, padValue=1)
  r3 <- raster::focal(r2, w=matrix(1,5,5), fun=my_outline, pad=T, padValue=1)
  
  star.m = matrix(rep(c(0,1), length.out=9), 3, 3)
  star.m[5] = 1
  r3_s <- raster::focal(r3, star.m, fun=my_fill_star, pad=T, padValue=1)
  
  rev_r3 = r3_s
  rev_r3[r3_s==0] = 1
  rev_r3[r3_s==1] = 0
  rev_c3 = raster::clump(rev_r3)
  
  message("[focal_transformations] Step 6: Extracting values...")
  vals <- raster::getValues(rev_c3)     # change
  tbl  <- table(vals, useNA = "no")     # change
  # FIXED: Now uses parameter instead of hard-coded 40
  flip_clump = as.numeric(names(tbl)[tbl < min_cluster_size])
  r4 = r3_s
  r4[rev_c3 %in% flip_clump] = 1
  
  return(r4)
}

lookupKey <- function(.xy, results_df) {
  key1 = cbind(.xy, "numeric_key"=seq(nrow(.xy)))
  t1 = raster::rasterFromXYZ(key1)
  rast = as.matrix(t1)
  
  num_key = sapply(1:nrow(results_df), function(x) {
    v1 = as.numeric(results_df[x,])
    rast[v1[1],v1[2]]
  })
  edge_spots = rownames(key1)[key1$numeric_key %in% num_key] 
  
  return(edge_spots)
}

lookupKeyDF <- function(.xy, results_df) {
  key1 = cbind(.xy, "numeric_key"=seq(nrow(.xy)))
  t1 = raster::rasterFromXYZ(key1)
  rast = as.matrix(t1)
  
  num_key = do.call(rbind, lapply(1:nrow(results_df), function(x) 
    cbind.data.frame("numeric_key"=rast[results_df[x,1],results_df[x,2]], "clump_id"=results_df[x,3], "size"=results_df[x,4])
  ))
  
  matching = match(num_key$numeric_key, key1$numeric_key)
  result <- cbind.data.frame("spotcode"=rownames(key1)[matching],"clumpID"=num_key$clump_id, "clumpSize"=num_key$size)
  
  return(result)
}

#' Enhanced clump edges detection with parameterized thresholds
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
  
  t1 = raster::rasterFromXYZ(.xyz)
  # ENHANCED: Pass min_cluster_size to focal_transformations
  t2 <- focal_transformations(t1, min_cluster_size = min_cluster_size)
  rast = as.matrix(t2)
  c1 = raster::clump(t2, directions=8)
  
  clumps = as.matrix(c1)
  edgeClumps = c()
  
  # ENHANCED: Use parameterized edge_threshold instead of hard-coded 0.75
  north = sum(!is.na(clumps[1,]))/sum(!is.na(rast[1,]))
  if(north>=edge_threshold) edgeClumps = c(edgeClumps, unique(clumps[1,]))
  
  east = sum(!is.na(clumps[,ncol(clumps)]))/sum(!is.na(rast[,ncol(rast)]))
  if(east>=edge_threshold) edgeClumps = c(edgeClumps, unique(clumps[,ncol(clumps)]))
  
  south = sum(!is.na(clumps[nrow(clumps),]))/sum(!is.na(rast[nrow(rast),]))
  if(south>=edge_threshold) edgeClumps = c(edgeClumps, unique(clumps[nrow(clumps),]))
  
  west = sum(!is.na(clumps[,1]))/sum(!is.na(rast[,1]))
  if(west>=edge_threshold) edgeClumps = c(edgeClumps, unique(clumps[,1]))
  
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

#' Enhanced problem areas detection with parameterized cluster size
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
  
  t1 = raster::rasterFromXYZ(.xyz)
  # ENHANCED: Pass min_cluster_size to focal_transformations
  t2 = focal_transformations(t1, min_cluster_size = min_cluster_size)
  rast = as.matrix(t2)
  c1 = raster::clump(t2, directions=8)
  
  clumps = as.matrix(c1)
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
