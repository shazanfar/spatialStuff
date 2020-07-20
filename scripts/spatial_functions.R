# last updated: 15/06/2020
# some utility functions for spatial SeqFISH analysis


# function to rotate the spatial coords
rotateDF = function(DF, 
                    xname = "segmentation_vertices_x_global_affine", 
                    yname = "segmentation_vertices_y_global_affine", 
                    ang = 0) {
  # ang is a numeric vector named three values corresponding to embryo1, embryo2, and embryo3
  
  ang_long = ang[as.character(DF$embryo)]
  ang_rad = ang_long/180
  
  x = DF[,xname]
  y = DF[,yname]
  
  x_turn = x*cos(ang_rad) - y*sin(ang_rad)
  y_turn = x*sin(ang_rad) + y*cos(ang_rad)
  
  # reset the columns and then return the DF
  
  DF[,xname] <- x_turn
  DF[,yname] <- y_turn
  
  return(DF)
  
}

filterAreas = function(areas) {
  # areas is a numeric vector, output is a logical vector giving TRUE (keep)
  # or FALSE (remove)
  
  sqrt_centre = median(sqrt(areas))
  sqrt_scale = mad(sqrt(areas))
  
  pnorm.fit = pnorm(sqrt(areas), mean = sqrt_centre, 
                    sd = sqrt_scale, lower.tail = F)
  area_thresh.sqrt = min(sqrt(areas)[which(p.adjust(pnorm.fit, method = "BH") < 0.01)])
  area_thresh = area_thresh.sqrt^2
  
  return(areas < area_thresh)
}


getSegmentationVerticesDF = function(DF,
                                     xname = "segmentation_vertices_x_global",
                                     yname = "segmentation_vertices_y_global",
                                     othercols = c("uniqueID","z")) {
  # DF is a DataFrame object
  # othercols is the others to keep
  stopifnot( sum(colnames(DF) %in% c(xname, yname)) == 2 ) 
  stopifnot( sum(colnames(DF) %in% othercols) == length(othercols) ) 
  
  long_x = unlist(DF[,xname])
  long_y = unlist(DF[,yname])
  
  if (length(long_x) != length(long_y)) stop("x and y need to be same length")
  
  long_xy = data.frame(
    long_x,
    long_y
  )
  
  long_DF = cbind(
    rep(DF[,othercols], times = unlist(lapply(DF[,xname], length))),
    long_xy
  )
  colnames(long_DF) <- c(othercols, xname, yname)
  return(as.data.frame(long_DF))
}
