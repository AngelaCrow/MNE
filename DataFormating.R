
# Regionalization shapefile folder
shapePath <- '../data/shapes/'
shapeLayer <- "wwf_terr_ecos_a" 
regionalizacion <- rgdal::readOGR(shapePath, shapeLayer)

# Raster covariables folder
covarDataFolder <- '../data/covar_presente' #present climate variables

# Raster covariables folder where the model will be projected
# IMPORTANT: The raster files on `covarDataFolder` and `covarAOIDataFolder` 
# must have the same name in order to the model can be evaluated.

# Future climate (2015-2039) and two rcp´s
covarAOIDataFolder_fc45 <- '../data/covar_fc45'
covarAOIDataFolder_fc85 <- '../data/covar_fc85'
# Future climate (2079-2099) and two rcp´s
covarAOIDataFolder_fl45 <- '../data/covar_fl45'
covarAOIDataFolder_fl85 <- '../data/covar_fl85'

# species records
args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Please enter a single parameter (input file).\n", call. = FALSE)
} else if (length(args) == 1) {
  print(paste("Processing model for file ", args[1]))
} else {
  stop("Single parameter is needed (input file).\n", call. = FALSE)
}

# inputDataFile <- args[1]
inputDataFile <- '../data/Physalis_subrepens.csv'
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

covarAOIFileList <- list_files_with_exts(covarAOIDataFolder, "tif")
enviromentalVariablesAOI <- raster::stack(covarAOIFileList)

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

selectedVariables <- enviromentalVariables[[select_var]]
selectedVariablesAOI <- enviromentalVariablesAOI[[select_var]]

#Future
covarAOIFileList_fc45 <- list_files_with_exts(covarAOIDataFolder_fc45, "tif")
enviromentalVariablesAOI_fc45 <- raster::stack(covarAOIFileList_fc45)
selectedVariablesAOI_fc45 <- enviromentalVariablesAOI_fc45[[select_var]]

covarAOIFileList_fc85 <- list_files_with_exts(covarAOIDataFolder_fc85, "tif")
enviromentalVariablesAOI_fc85 <- raster::stack(covarAOIFileList_fc85)
selectedVariablesAOI_fc85 <- enviromentalVariablesAOI_fc85[[select_var]]

covarAOIFileList_fl45 <- list_files_with_exts(covarAOIDataFolder_fl45, "tif")
enviromentalVariablesAOI_fl45 <- raster::stack(covarAOIFileList_fl45)
selectedVariablesAOI_fl45 <- enviromentalVariablesAOI_fl45[[select_var]]

covarAOIFileList_fl85 <- list_files_with_exts(covarAOIDataFolder_fl85, "tif")
enviromentalVariablesAOI_fl85 <- raster::stack(covarAOIFileList_fl85)
selectedVariablesAOI_fl85 <- enviromentalVariablesAOI_fl85[[select_var]]





