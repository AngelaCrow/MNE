# 
# This code was developed by
# - Angela P. Cuervo-Robayo ancuervo@gmail.com
# - Juan M. Barrrios j.m.barrios@gmail.com
# 
# The clean_dup function was develop by Luis Osorio as part of the NicheToolbox
# project [https://github.com/luismurao/nichetoolbox]
#

library("rgdal")
library("fuzzySim")
library("ENMeval")
library("ROCR")
library("magrittr")
library("dplyr")
library("tools")
library("data.table")
library("raster")

source("utils/clean_dup.R")

set.seed(1)

args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Al menos un parametro es requerido (input file).\n", call. = FALSE)
} else if (length(args) == 1) {
  print(paste("Se calcula el modelo para", args[1]))
} else {
  stop("Solo se acepta un parametro.\n", call. = FALSE)
}

# Regionalization shapefile folder
shapePath <- file.path('.', 'IUCN_data', 'shapes')
shapeLayer <- "wwf_terr_ecos_a"
regionalizacion <- readOGR(shapePath, shapeLayer)

# Raster covariables folder
covarDataFolder <- file.path('.', 'IUCN_data', "covar_rasters")

# Raster covariables folder where the model will be projected
# IMPORTANT: The raster files on `covarDataFolder` and `covarAOIDataFolder` 
# must have the same name in order to the model can be evaluated.
covarAOIDataFolder <- file.path('.', 'IUCN_data', "covar_raster_PSC")

inputDataFile <- args[1]
outputFolder <- inputDataFile %>%
  basename %>%
  file_path_sans_ext

if (!dir.exists(outputFolder)) {
  dir.create(outputFolder, recursive = TRUE)
}

####Cleaning duplicate records on a cell####
occsData <- read.csv(inputDataFile, header = TRUE, stringsAsFactors = FALSE) %>%
  clean_dup("Dec_Long", "Dec_Lat", 0.00833333333)

write.csv(occsData,
          file = file.path(outputFolder, "clean_data.csv"),
          row.names = FALSE)

#### ENVIROMENTAL VARIABLES####
covarFileList <- list_files_with_exts(covarDataFolder, "tif")
enviromentalVariables <- raster::stack(covarFileList)

covarAOIFileList <- list_files_with_exts(covarAOIDataFolder, "tif")
enviromentalVariablesAOI <- raster::stack(covarAOIFileList)

#### VARIABLES + PRESENCIAS####
sp_coor <- occsData[c("Dec_Long", "Dec_Lat")]

covarData <- raster::extract(enviromentalVariables, sp_coor)
covarData <- cbind(occsData, covarData)

covarData <- covarData[!is.na(covarData$bio_1), ]

####SELECCION DE VARIABLES####
speciesCol <- match("Presence", names(occsData))
varCols <- ncol(occsData) + 1

correlacion <- corSelect(
  data = covarData,
  sp.cols = speciesCol,
  var.cols = varCols:ncol(covarData),
  cor.thresh = 0.8,
  use = "pairwise.complete.obs"
)

select_var <- correlacion$selected.vars
write(select_var, file = file.path(outputFolder, "selected_variables.txt"))
selectedVariables <- enviromentalVariables[[select_var]]

selectedVariablesAOI <- enviromentalVariablesAOI[[select_var]]

####TRAINNING###
# Divides your data into trainining and test data sets. 70/30 %
sampleDataPoints <- sample.int(
  nrow(covarData),
  size = floor(0.7*nrow(covarData))
)

selectedValues <- rep(0, nrow(covarData)) %>% inset(sampleDataPoints, 1)

covarData$isTrain <- selectedValues
write.csv(covarData, file.path(outputFolder, "finaldatabase.csv"),
          row.names = FALSE)

# Selects the M of the species, base on Olson´s ecoregions
# Download: https://www.worldwildlife.org/publications/terrestrial-ecoregions-of-the-world
coordinates(sp_coor) <- ~Dec_Long+Dec_Lat
proj4string(sp_coor) <- proj4string(regionalizacion)

# Intersects the occurrence data with polygons
dataenpoly <- over(sp_coor, regionalizacion, fn = NULL)
enpolyindex <- which(!is.na(dataenpoly$ECO_NAME))
polydatadf <- dataenpoly[enpolyindex, ]
id_polys <- unique(polydatadf$ECO_NAME)
poligonofilter <- regionalizacion[regionalizacion$ECO_NAME %in% id_polys, ]
# extract by mask
selectedVariablesCrop <- crop(selectedVariables, poligonofilter)
env <- mask(selectedVariablesCrop, poligonofilter) #Species variables delimited by M

# MAXENT
# We used ENMeval packeage to estimate optimal model complexity (Muscarrella et al. 2014)
# Modeling process, first separate the calibration and validation data
occsCalibracion <- covarData %>%
  dplyr::filter(isTrain == 1) %>%
  dplyr::select(Dec_Long, Dec_Lat)

occsValidacion <- covarData %>%
  dplyr::filter(isTrain == 0) %>%
  dplyr::select(Dec_Long, Dec_Lat) %>%
  as.data.frame

# Background
bg <- randomPoints(env[[1]], n = 10000)
# bg <- randomPoints(env[[1]], n = 100)
bg.df <- as.data.frame(bg)
write.csv(bg.df, file = file.path(outputFolder, "background_data.csv"),
          row.names = FALSE)

# ENMeval
sp <- ENMevaluate(occsCalibracion, env, bg.df, RMvalues = seq(0.5, 4, 0.5),
                 fc = c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"),
                 method = "randomkfold", kfolds = 5, bin.output = TRUE,
                 parallel = TRUE, numCores = parallel::detectCores())

resultados_enmeval <- sp@results
write.csv(resultados_enmeval,
          file = file.path(outputFolder, "enmeval_results.csv"),
          row.names = FALSE)

# delta_aic <- which(resultados_enmeval$delta.AICc == 0)
modelsAIC0 <- resultados_enmeval %>%
  mutate(index = rownames(resultados_enmeval)) %>%
  filter(delta.AICc == 0) %>%
  select(index, settings) %>%
  mutate(index = as.numeric(index), settings = as.character(settings))



saveRasterWithSettings <- function(models, predictions, prefix) {
  raster::writeRaster(predictions[[models["settings"]]],
              file.path(outputFolder, paste0(prefix,
                                             models["settings"],
                                             ".tif")),
              overwrite = TRUE)
}

apply(modelsAIC0, 1, saveRasterWithSettings,
      predictions = sp@predictions, prefix = "ENM_prediction_M_raw_")

predictAndSave <- function(model, models, data, prefix, occs) {
  choicedModel <- models[[as.integer(model["index"])]]
  predictions <- dismo::predict(choicedModel, data)
  raster::writeRaster(predictions,
                      file.path(outputFolder, paste0(prefix,
                                                     "log_",
                                                     model["settings"],
                                                     ".tif")),
                      overwrite = TRUE)

  occsValues <- raster::extract(predictions, occs)
  minValOcc <- min(occsValues, na.rm = TRUE)
  raster::writeRaster(reclassify(predictions,
                                 c(-Inf, minValOcc, 0, minValOcc, Inf, 1)),
                      file.path(outputFolder, paste0(prefix,
                                                     "bin_min_",
                                                     model["settings"],
                                                     ".tif")),
                      overwrite = TRUE)

  q10ValOcc <- quantile(occsValues, 0.1, na.rm = TRUE)
  raster::writeRaster(reclassify(predictions,
                                 c(-Inf, q10ValOcc, 0, q10ValOcc, Inf, 1)),
                      file.path(outputFolder, paste0(prefix,
                                                     "bin_q10_",
                                                     model["settings"],
                                                     ".tif")),
                      overwrite = TRUE)
}

apply(modelsAIC0, 1, predictAndSave,
      models = sp@models, data = env, prefix = "ENM_prediction_M_",
      occs = occsCalibracion)

apply(modelsAIC0, 1, predictAndSave,
      models = sp@models, data = selectedVariablesAOI, prefix = "ENM_",
      occs = occsCalibracion)

####VALIDACION####
#Independiente de umbral
#AUC
aucCalculator <- function(prediction, occs, bgPoints) {
  data <- rbind(occs, setNames(bgPoints, names(occs)))
  labels <- c(rep(1, nrow(occs)),
              rep(0, nrow(bgPoints)))
  scores <- raster::extract(prediction, data)
  pred <- prediction(scores, labels)
  # perf <- performance(pred, "tpr", "fpr")
  auc <- performance(pred, "auc")@y.values[[1]]
  return(auc)
}

aucStatistcs <- function(model, models, env, occs, bgPoints) {
  result <- apply(model, 1, function(x, models, env, occs, bgPoints){
    choicedModel <- models[[as.integer(model["index"])]]
    prediction <- dismo::predict(choicedModel, env)
    auc <- aucCalculator(prediction, occs, bgPoints)
    return(c(model["settings"], auc))
  },
  models = models,
  env = env,
  occs = occs,
  bgPoints = bgPoints)

  result <- data.frame(
    matrix(unlist(result), nrow = nrow(model), byrow = TRUE),
    stringsAsFactors = FALSE
    )

  names(result) <- c("settings", "AUC")

  result <- result %>% mutate(AUC = as.numeric(AUC))

  return(result)
}

resultsAUC <- aucStatistcs(modelsAIC0, sp@models, env, occsValidacion, bg.df)
write.csv(resultsAUC,
          file = file.path(outputFolder, "data_auc.csv"),
          row.names = FALSE)

