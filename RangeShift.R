sp_models <- list.files(path = "acutiacuti",
                        pattern = "*min_L_0.5.tif$",full.names = TRUE)
sp_models

presente<-raster(sp_models[1]) 
plot(presente)
fc45bin<-raster(sp_models[2])
plot(fc45bin)
fc85bin<-raster(sp_models[3])
plot(fc85bin)
fl45bin<-raster(sp_models[4])
plot(fl45bin)
fl85bin<-raster(sp_models[5])
plot(fl85bin)

m <- c(0, 0, 0,  0,1, 10)
rclmat <- matrix(m, ncol=2, byrow=TRUE)


rc_fc45bin <- reclassify(fc45bin, rclmat)
rc_fc45bin
plot(rc_fc45bin)
rangeShift <- presente+rc_fc45bin
plot(rangeShift)
raster::writeRaster(rangeShift, file.path("C:/CONABIO/ModelosDarwinUICN/acutiacuti/_rsfc45.tif"))