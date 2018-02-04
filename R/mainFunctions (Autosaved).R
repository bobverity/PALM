
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
	
	# run efficient C++ code
	args_data <- list(hex_data=hex_data, hex_neighbors=hex_neighbors, stat=mat_to_Rcpp(stat))
	args_model <- list(reps=reps)
	output_raw <- fitModel_cpp(args_data, args_model)
	
	# process output
	output <- list()
	output$hex_pts <- hex_pts
	output$hex_polys <- hex_polys
	output$hex_neighbors <- hex_neighbors
	output$hex_data <- hex_data
	output$path_mat <- output_raw$path_mat
	output$s <- Rcpp_to_mat(output_raw$s)
	output$hex_npath <- 	output_raw$hex_npath
	output$friction <- output_raw$friction
	output$alpha <- output_raw$alpha
	output$beta <- output_raw$beta
	return(output)
	
	
	#output <- list()
	#output$grid <- raster(Rcpp_to_mat(output_raw$grid)[21:1,], xmn=-10.5, xmx=10.5, ymn=-10.5, ymx=10.5)
	#output$g <- Rcpp_to_mat(output_raw$grid)
	#output$s <- Rcpp_to_mat(output_raw$s)
	
	#return(output)
}
