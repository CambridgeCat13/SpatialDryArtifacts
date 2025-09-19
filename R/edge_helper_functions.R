# ==============================================================================
# Fixed Edge Detection Helper Functions
# Solves both dimnames array error AND maintains correct spatial detection
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
  
  # Safe matrix conversion using values()
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

lookupKey <- function(.xy, results_df) {
  if (nrow(results_df) == 0) {
    return(character(0))
  }
  
  # Create a key mapping with original coordinates (no swap)
  key1 = cbind(.xy, "numeric_key"=seq(nrow(.xy)))
  
  # Create raster for mapping - using original coordinates
  t1 = raster::rasterFromXYZ(key1)
  
  # Safely convert to matrix
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

lookupKeyDF <- function(.xy, results_df) {
  if (nrow(results_df) == 0) {
    return(data.frame(spotcode = character(0), 
                     clumpID = character(0), 
                     clumpSize = numeric(0)))
  }
  
  # Create a key mapping with original coordinates (no swap)
  key1 = cbind(.xy, "numeric_key"=seq(nrow(.xy)))
  
  # Create raster for mapping
  t1 = raster::rasterFromXYZ(key1)
  
  # Safely convert to matrix
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

clumpEdges <- function(.xyz, offTissue, shifted=FALSE, edge_threshold=0.75, min_cluster_size=40) {
  # Handle NA values
  if(ncol(.xyz) >= 3) {
    .xyz[, 3][is.na(.xyz[, 3])] <- FALSE
  }
  
  if(sum(.xyz[,3])==0) {
    return(c())
  }
  
  # Apply shifting if needed
  if(shifted==TRUE) {
    odds = seq(1, max(.xyz[,"array_col"]), by=2)
    .xyz[.xyz[,"array_col"] %in% odds, "array_col"] = 
      .xyz[.xyz[,"array_col"] %in% odds, "array_col"]-1
  }
  
  # IMPORTANT: Use original coordinates, no swap needed
  t1 = raster::rasterFromXYZ(.xyz)
  t2 <- focal_transformations(t1, min_cluster_size = min_cluster_size)
  
  # Safely extract values
  rast_values = raster::values(t2)
  rast_dims = c(nrow(t2), ncol(t2))
  rast = matrix(rast_values, nrow = rast_dims[1], ncol = rast_dims[2], byrow = TRUE)
  
  c1 = raster::clump(t2, directions=8)
  
  clump_values = raster::values(c1)
  clump_dims = c(nrow(c1), ncol(c1))
  clumps = matrix(clump_values, nrow = clump_dims[1], ncol = clump_dims[2], byrow = TRUE)
  
  edgeClumps = c()
  
  # Check each border
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
  
  # Apply shifting if needed
  if(shifted==TRUE) {
    odds = seq(1, max(.xyz[,"array_col"]), by=2)
    .xyz[.xyz[,"array_col"] %in% odds, "array_col"] = 
      .xyz[.xyz[,"array_col"] %in% odds, "array_col"]-1
  }
  
  # IMPORTANT: Use original coordinates, no swap needed
  t1 = raster::rasterFromXYZ(.xyz)
  t2 = focal_transformations(t1, min_cluster_size = min_cluster_size)
  
  c1 = raster::clump(t2, directions=8)
  
  # Safely extract values
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