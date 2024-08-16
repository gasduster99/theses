library(dismo)

#
U = cbind(runif(100), runif(100))
r = raster(matrix(1:1000000, 1000, 1000))
#
v=voronoi(U)
rast = rasterize(v, r, U[,2])
#
plot(rast, col='red')

