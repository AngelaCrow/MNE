# Imports ---------------------------------------------------------------------
library(raster)
library(rgdal)
library(rgeos)


# Functions -------------------------------------------------------------------

save_difference <- function(raster1, raster2, output_name, output_folder) {
  output_raster <- raster::overlay(raster1,
                                   raster2,
                                   fun = function(x, y) {
                                     x + y
                                   })
  
  output_polygon <- raster::rasterToPolygons(output_raster,
                                             fun = function(x) {
                                               x == 1
                                             },
                                             dissolve = TRUE)
  
  rgdal::writeOGR(output_polygon,
                  dsn = output_folder,
                  layer = output_name,
                  driver = "ESRI Shapefile",
                  overwrite_layer = TRUE)
}

# Main ------------------------------------------------------------------------

## Definitions ----------------------------------------------------------------
BASE_DIRECTORY <- ("C:/CONABIO/ModelosDarwinUICN/") #"Poner la direccion del folder con las carpetas de las especies"
NAME_FOLDERS1 <- ("Con_clamp_con_extra")
NAME_FOLDERS2 <- ("Sin_clamp_sin_extra")


## Load species names ---------------------------------------------------------
species_folders <- list.dirs(BASE_DIRECTORY, recursive = FALSE)
species_folders <- species_folders[grep("_pvalue$", species_folders, invert = TRUE)]

## Read and write rasters  ----------------------------------------------------
for (folder in species_folders) {
  print(file.path(folder,
                  NAME_FOLDERS1))
  species_rasters_fn <- list.files(path = file.path(folder,
                                                    NAME_FOLDERS1),
                                   pattern = "thresholded.asc$",
                                   full.names = FALSE,
                                   recursive = FALSE,
                                   include.dirs = FALSE)
  # print(species_rasters_fn)
  for (raster_file in species_rasters_fn) {
    print(raster_file)
    raster1 <- raster::raster(file.path(folder,
                                        NAME_FOLDERS1,
                                        raster_file))
    
    raster2 <- raster::raster(file.path(folder,
                                        NAME_FOLDERS2,
                                        raster_file))
    model_name <- tools::file_path_sans_ext(raster_file)
    
    save_difference(raster1, 
                    raster2, 
                    paste0('uncertainty_', model_name),
                    file.path(folder))
  }
}


