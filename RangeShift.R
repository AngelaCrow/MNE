#' Report statistics about region gain and loss for future environment
#' 
#' @param raster_presente_bin Raster object with values 0, 1 that shows the
#'                            current species distribution
#' @param raster_futuro_bin Raster object with values 0, 1 that show the future
#'                          ditribution
#'
#' @return A list with the frequency stats for habitat gain or loss coded as 
#'         follows:
#'         \itemize{
#'           \item{0 no habitat in neither two times}
#'           \item{1 habitat loss}
#'           \item{10 habitat gain}
#'           \item{11 preserved habitat} 
#'         }
#'         and a raster with the same coding.
#' 
#' @examples
#' \dontrun{
#' r_pres <- raster('present.tif')
#' r_fut <- raster('future.tif')
#' habitat_change <- reportaFuturo(r_pres, r_fut)
#' plot(habitat_change$raster_cambio)
#' print(habitat_change$stats)
#' }
#' 

setwd("C:/CONABIO/ModelosDarwinUICN/")
######
reportaFuturo <- function(raster_presente_bin, raster_futuro_bin){
  reclass_matrix <- matrix(c(0, 1, 0, 10), ncol=2)
  offset_raster <- raster::reclassify(raster_futuro_bin, reclass_matrix)
  
  result <- raster_presente_bin + offset_raster
  
  return(list(stats=raster::freq(result, useNA='no'),
              raster_cambio=result))
}

sp_models <- list.files(path = "vulgasilvestre/Con_clamp_con_extra/",
                        pattern = "*thresholded.asc$",full.names = TRUE)
sp_models
presente<-raster(sp_models[5]) 
plot(presente)

fc45bin<-raster(sp_models[1])
plot(fc45bin)
proyeccion <- reportaFuturo(presente, fc45bin)
statsfc45bin<-data.frame(proyeccion$stats)
write.csv(statsfc45bin, "vulgasilvestre/Con_clamp_con_extra/statsfc45bin.csv")
plot(proyeccion$raster_cambio)
raster::writeRaster(proyeccion$raster_cambio, 
                    file.path("C:/CONABIO/ModelosDarwinUICN/vulgasilvestre/Con_clamp_con_extra/vulgasilvestre_rsfc45.tif"),
                    overwrite = T)

fc85bin<-raster(sp_models[3])
plot(fc85bin)
proyeccion <- reportaFuturo(presente, fc85bin)
statsfc85bin<-data.frame(proyeccion$stats)
write.csv(statsfc85bin, "vulgasilvestre/Con_clamp_con_extra/statsfc85bin.csv")
plot(proyeccion$raster_cambio)
raster::writeRaster(proyeccion$raster_cambio, 
                    file.path("C:/CONABIO/ModelosDarwinUICN/vulgasilvestre/Con_clamp_con_extra/vulgasilvestre_rsfc85.tif"),
                    overwrite = T)
sp_models

fl45bin<-raster(sp_models[2])
plot(fl45bin)
proyeccion <- reportaFuturo(presente, fl45bin)
statsfl45bin<-data.frame(proyeccion$stats)
write.csv(statsfl45bin, "vulgasilvestre/Con_clamp_con_extra/statsfl45bin.csv")
plot(proyeccion$raster_cambio)
raster::writeRaster(proyeccion$raster_cambio, 
                    file.path("C:/CONABIO/ModelosDarwinUICN/vulgasilvestre/Con_clamp_con_extra/vulgasilvestre_rsfl45.tif"),
                    overwrite = T)


fl85bin<-raster(sp_models[4])
plot(fl85bin)
proyeccion <- reportaFuturo(presente, fl85bin)
statsfl85bin<-data.frame(proyeccion$stats)
write.csv(statsfl85bin, "vulgasilvestre/Con_clamp_con_extra/Con_clamp_con_extrastatsfl85bin.csv")
plot(proyeccion$raster_cambio)
raster::writeRaster(proyeccion$raster_cambio, 
                    file.path("C:/CONABIO/ModelosDarwinUICN/vulgasilvestre/Con_clamp_con_extra/vulgasilvestre_rsfl85.tif"),
                    overwrite = T)


