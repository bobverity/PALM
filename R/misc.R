
#------------------------------------------------
# The following commands ensure that package dependencies are listed in the NAMESPACE file.

#' @useDynLib PALM
#' @import sp
#' @import rgeos
#' @import raster
#' @import viridis
#' @importFrom Rcpp evalCpp
NULL

# -----------------------------------
# mat_to_Rcpp
# takes matrix as input, converts to list format for use within Rcpp code
# (not exported)

mat_to_Rcpp <- function(x) {
    return(split(x,f=1:nrow(x)))
}

# -----------------------------------
# Rcpp_to_mat
# Takes list format returned from Rcpp and converts to matrix.
# (not exported)

Rcpp_to_mat <- function(x) {
    ret <- matrix(unlist(x), nrow=length(x), byrow=TRUE)
    return(ret)
}

# -----------------------------------
#' Convert polygons to raster
#' 
#' Takes an object of class "SpatialPolygons" and converts to raster image. The argument \code{z} gives the value of each polygon that will be used when generating the raster.
#' 
#' @parameter p object of class "SpatialPolygon"
#' 
#' @export

poly_to_raster <- function(p, z, rows=1e2, cols=1e2) {
	ex <- extent(p)
	r <- raster(ncol=cols, nrow=rows, xmn=ex[1], xmx=ex[2], ymn=ex[3], ymx=ex[4])
	ret <- rasterize(p, r, field=z)
	return(ret)
}

# -----------------------------------
#' Convert polygons to segments
#' 
#' Takes an object of class "SpatialPolygons" and extracts the coordinates of all polygons. Output is in the form of a matrix with columns "x0", "y0", "x1", "y1", corresponding to the start and end coordinates of each edge. These can be plotted using the \code{segments()} function, which is often faster than plotting raw polygons.
#' 
#' @export

poly_to_segment <- function(p) {
	ret <- lapply(p@polygons, function(x){
			coords <- t(x@Polygons[[1]]@coords)
			rbind(coords[,-ncol(coords)], coords[,-1])
		})
	ret <- matrix(unlist(ret), ncol=4, byrow=TRUE)
    colnames(ret) <- c("x0", "y0", "x1", "y1")
    return(ret)
}


