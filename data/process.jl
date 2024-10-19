using MAT
using GLMakie, MakieLayers
using Ipaper, Ipaper.sf, ArchGDAL

cellsize = 0.1
b = bbox(-180, -90, 180, 90)
lon, lat = bbox2dims(b; cellsize)

soil = matread("./data/inpara/Soilraster.mat")["Soilraster"]' |> collect
write_gdal(rast(soil, b), "data/param_Soil_G010.tif")
# images(lon, lat, soil)

Topt = matread("./data/inpara/Topt.mat")["Topt_new"]' |> collect
write_gdal(rast(Topt, b), "data/param_Topt_G010.tif")