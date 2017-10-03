
#### Input data ####
# Regionalization shapefile folder
shapePath <- '../data/shapes/'
shapeLayer <- "wwf_terr_ecos_a" 
regionalizacion <- rgdal::readOGR(shapePath, shapeLayer)

# Present raster covariables folder
covarDataFolder <- '../data/covar_presente' 
# IMPORTANT: The raster files on Present covarDataFolder and Future covarDataFolder
# must have the same name in order to the model can be evaluated.

# Future climate (2015-2039) and two rcp´s
covarDataFolder_fc45 <- '../data/covar_fc45'
covarDataFolder_fc85 <- '../data/covar_fc85'
# Future climate (2079-2099) and two rcp´s
covarDataFolder_fl45 <- '../data/covar_fl45'
covarDataFolder_fl85 <- '../data/covar_fl85'

# Esta parte necesita ser explicada
args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Please enter a single parameter (input file).\n", call. = FALSE)
} else if (length(args) == 1) {
  print(paste("Processing model for file ", args[1]))
} else {
  stop("Single parameter is needed (input file).\n", call. = FALSE)
}

inputDataFile <- args[1]
outputFolder <- inputDataFile %>%
  basename %>%
  file_path_sans_ext

if (!dir.exists(outputFolder)) {
  dir.create(outputFolder, recursive = TRUE)
}

####Cleaning duplicate records on a cell####
crs.wgs84 <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
occsData <- readr::read_csv(inputDataFile)
sp::coordinates(occsData) <- c("Dec_Long", "Dec_Lat")
sp::proj4string(occsData) <- crs.wgs84

occsData <- sp::remove.duplicates(occsData, zero=0.00833333333)

write.csv(cbind(occsData@data, coordinates(occsData)),
          file = file.path(outputFolder, "data_wo_duplicates.csv"),
          row.names = FALSE)


#### ENVIROMENTAL VARIABLES####
#Present
covarFileList <- list_files_with_exts(covarDataFolder, "tif")
enviromentalVariables <- raster::stack(covarFileList)

# Extract envorimental varibales with species occurrences
covarData <- raster::extract(enviromentalVariables, occsData)
covarData <- cbind(occsData, covarData)

completeDataCases <- covarData@data %>% 
  dplyr::select_(.dots=names(enviromentalVariables)) %>%
  complete.cases
covarData <- covarData[completeDataCases, ]

####Variables selection####
speciesCol <- match("Presence", names(occsData))
varCols <- ncol(occsData) + 1

correlacion <- corSelect(
  data = covarData@data,
  sp.cols = speciesCol,
  var.cols = varCols:ncol(covarData),
  cor.thresh = 0.8,
  use = "pairwise.complete.obs"
)

select_var <- correlacion$selected.vars
write(select_var, file = file.path(outputFolder, "selected_variables.txt"))

# Raster covariables selected for model calibration
selectedVariables <- enviromentalVariables[[select_var]]

# Selects the M of the species, base on Olson´s ecoregions
# Download: https://www.worldwildlife.org/publications/terrestrial-ecoregions-of-the-world
# Intersects the occurrence data with polygons
ecoregionsOfInterest <- sp::over(occsData, regionalizacion) %>%
  filter(!is.na(ECO_ID))

idsEcoRegions <- unique(ecoregionsOfInterest$ECO_ID)
polygonsOfInterest <- regionalizacion[regionalizacion$ECO_ID %in% idsEcoRegions, ]
writeOGR(polygonsOfInterest, layer = 'ecoregionsOI', outputFolder, driver="ESRI Shapefile")

# Mask present rasters with ecoregions of interest
selectedVariablesCrop <- raster::crop(selectedVariables, polygonsOfInterest)
env <- raster::mask(selectedVariablesCrop, 
                    polygonsOfInterest) #Species variables delimited by M
writeRaster(env,
            file.path(outputFolder, "covars.tif"), 
            bylayer = T, suffix='names',
            overwrite = TRUE)


# Mask future raster with ecoregions of interest
# fc_45
covarFileList_fc45 <- list_files_with_exts(covarDataFolder_fc45, "tif")
enviromentalVariables_fc45 <- raster::stack(covarFileList_fc45)
selectedVariables_fc45 <- enviromentalVariables_fc45[[select_var]]
selectedVariablesCrop_fc45 <- raster::crop(selectedVariables_fc45, polygonsOfInterest)
env_fc45 <- raster::mask(selectedVariablesCrop_fc45, 
                    polygonsOfInterest) 

# fc_85
covarFileList_fc85 <- list_files_with_exts(covarDataFolder_fc85, "tif")
enviromentalVariables_fc85 <- raster::stack(covarFileList_fc85)
selectedVariables_fc85 <- enviromentalVariables_fc85[[select_var]]
selectedVariablesCrop_fc85 <- raster::crop(selectedVariables_fc85, polygonsOfInterest)
env_fc85 <- raster::mask(selectedVariablesCrop_fc85, 
                         polygonsOfInterest) 

# fl_45
covarFileList_fl45 <- list_files_with_exts(covarDataFolder_fl45, "tif")
enviromentalVariables_fl45 <- raster::stack(covarFileList_fl45)
selectedVariables_fl45 <- enviromentalVariables_fl45[[select_var]]
selectedVariablesCrop_fl45 <- raster::crop(selectedVariables_fl45, polygonsOfInterest)
env_fl45 <- raster::mask(selectedVariablesCrop_fl45, 
                         polygonsOfInterest) 

# fl_85
covarFileList_fl85 <- list_files_with_exts(covarDataFolder_fl85, "tif")
enviromentalVariables_fl85 <- raster::stack(covarFileList_fl85)
selectedVariables_fl85 <- enviromentalVariables_fl85[[select_var]]
selectedVariablesCrop_fl85 <- raster::crop(selectedVariables_fl85, polygonsOfInterest)
env_fl85 <- raster::mask(selectedVariablesCrop_fl85, 
                         polygonsOfInterest) 


