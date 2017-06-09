# This code was developed by
# - Angela P. Cuervo-Robayo ancuervo@gmail.com y
# - Juan Martin Barrrios j.m.barrios@gmail.com.
# It uses NicheToolBox package functions
# of Luis Osorio: https://github.com/luismurao/nichetoolbox

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

# Data folder configuration
baseDataFolder <- "../IUCN_data"

# Directory of the shapePolygon to select species'  M area
shapePath <- file.path(baseDataFolder, 'shapes')
shapeLayer <- "wwf_eco_mesoa"
regionalizacion <- readOGR(shapePath, shapeLayer)

# Directory for covariables
covarDataFolder <- file.path(baseDataFolder, "covar_rasters")

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
pdf(file = file.path(outputFolder, "env_plot.pdf"))
  plot(env[[1]], legend = FALSE)
  points(bg.df, col = 'red')
dev.off()
write.csv(bg.df, file = file.path(outputFolder, "background_data.csv"),
          row.names = FALSE)

# ENMeval
sp <- ENMevaluate(occsCalibracion, env, bg.df, RMvalues = seq(0.5, 4, 0.5),
                 fc = c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"),
                 method = "randomkfold", kfolds = 4, bin.output = TRUE,
                 parallel = TRUE, numCores = parallel::detectCores())

# identificar el nombre del modelo más parsimonioso
results_enmeval <- sp@results
write.csv(results_enmeval,
          file = file.path(outputFolder, "results_enmeval.csv")
          ,row.names = FALSE)
delta_aic <- which(results_enmeval$delta.AICc == 0)

# ENM EN RASTER
# seleccionar raster del modelo más parsomonioso
predictions <- sp@predictions
model.delta.aic <- predictions[[delta_aic]]
writeRaster(model.delta.aic,
            file.path(outputFolder, "ENMeval_prediction.tif"),
            overwrite = TRUE)

# elegir los parametros de maxent del modelo más parsimonioso y
# proyectar a otras variables
model.maxent <- sp@models[[delta_aic]]

# Raster del modelo en la M, en escala logistica
mode.in.m <- predict(model.maxent, env)
#plot(mode.in.m)
writeRaster(mode.in.m,
            file.path(outputFolder, "ENM_prediction_M_log.tif"),
            overwrite = TRUE)

##tranferirlo a otro tiempo o area mas grande
mapa.area.grande <- predict(model.maxent, selectedVariables)
#plot(mapa.area.grande)
writeRaster(mapa.area.grande,
            file.path(outputFolder, "ENM_log.tif"),
            overwrite = TRUE)


####VALIDACION####
#Independiente de umbral
#AUC
testpp <- raster::extract(mode.in.m, occsValidacion)
abs <- raster::extract(mode.in.m, bg.df)
combined <- c(testpp, abs)
label <- c(rep(1,length(testpp)),rep(0,length(abs)))
pred <- prediction(combined, label)
perf <- performance(pred, "tpr", "fpr")
auc <- performance(pred, "auc")@y.values[[1]]
auc
write.csv(auc,
          file = file.path(outputFolder, "data_auc.csv"),
          row.names = FALSE)

#Dependiente de umbral
# source("funciones_LAE.R")
#reclasificar mapa de la calibracion
rcl.m <- na.omit(raster::extract(mode.in.m, occsCalibracion))

#usando el valor de minimo de idoneidad que tienen los puntos de occurencia
rcl.min <- min(rcl.m) # extraer el minimo valor de presencia

mode.in.m.bin <- reclassify(mode.in.m, c(-Inf, rcl.min, 0, rcl.min, Inf, 1)) # reclasificar - cambie su valor en donde esta el valor decimal
writeRaster(mode.in.m.bin,
            file.path(outputFolder, "ENM_bin_min.tif"),
            overwrite = TRUE)

mode.in.mg.MTPbin <- reclassify(mapa.area.grande, c(-Inf, rcl.min, 0, rcl.min, Inf, 1)) # reclasificar - cambie su valor en donde esta el valor decimal
writeRaster(mode.in.mg.MTPbin,
            file.path(outputFolder, "ENM_binG_MTP.tif"),
            overwrite = TRUE)

#10 percentil
rcl.10 <- quantile(na.omit(rcl.m), .10) # extraer valor por percentil

mode.in.m.bin10 <- reclassify(mode.in.m, c(-Inf, rcl.10, 0, rcl.10, Inf, 1)) # reclasificar - cambie su valor en donde esta el valor decimal
writeRaster(mode.in.m.bin10,
            file.path(outputFolder, "ENM_bin_10.tif"),
            overwrite = TRUE)

mode.in.mg.bin10 <- reclassify(mapa.area.grande, c(-Inf, rcl.10, 0, rcl.10, Inf, 1)) # reclasificar - cambie su valor en donde esta el valor decimal
writeRaster(mode.in.mg.bin10,
            file.path(outputFolder, "ENM_binG_10.tif"),
            overwrite = TRUE)
#
# ##Para validar los modelos binarios usar el codigo que se llama AllMetrics.R.
