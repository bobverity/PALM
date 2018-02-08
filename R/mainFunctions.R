
#------------------------------------------------
# fit the model
#' @export
fitModel <- function(dat, stat, reps=10, cellSize=1, buffer=2*cellSize) {
	
	# get convex hull of data in SpatialPolygons format
	ch <- chull(dat)
	coords <- dat[c(ch,ch[1]),]
	sp_poly <- SpatialPolygons(list(Polygons(list(Polygon(coords)), ID=1)))
	sp_poly <- gBuffer(sp_poly, width=buffer)
	
	# get hex centre points and polygons
	hex_pts <- spsample(sp_poly, type="hexagonal", cellsize=cellSize)
	hex_polys <- HexPoints2SpatialPolygons(hex_pts)
	nhex <- length(hex_polys)
	
	# get list of neighbors to each hex
	d <- as.matrix(dist(cbind(hex_pts$x, hex_pts$y)))
	hex_neighbors <- apply(d, 1, function(x){ as.vector(which(round(x/cellSize)==1))-1 })
	names(hex_neighbors) <- NULL
	
	# get which hex each data point falls in
	hex_data <- over(SpatialPoints(dat), hex_polys) - 1
	
	# check all data within map
	if (any(is.na(hex_data))) {
		stop("data falling outside hex map. Increase size of buffer.")
	}
	
	# initialise friction values
	friction <- 1 + 0*(abs(hex_pts$x) < 0.5)
	
	# process output
	output <- list()
	output$hex_pts <- hex_pts
	output$hex_polys <- hex_polys
	output$hex_neighbors <- hex_neighbors
	output$hex_data <- hex_data
	
	#return(output)
	
	# run efficient C++ code
	args_data <- list(hex_data=hex_data, hex_neighbors=hex_neighbors, stat=mat_to_Rcpp(stat))
	args_model <- list(reps=reps, friction=friction)
	output_raw <- fitModel_cpp(args_data, args_model)
	
	# add to output
	output <- c(output, output_raw)

	return(output)
}
