#Este codigo fue desarrollado por Angela P. Cuervo-Robayo ancuervo@gmail.com y Juan Martin Barrrios j.m.barrios@gmail.com . 

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
source("utils/utilidades.R")

set.seed(1)

args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Al menos un parametro es requerido (input file).\n", call. = FALSE)
} else if (length(args) == 1) {
  print(paste("Se calcula el modelo para", args[1]))
} else {
  stop("Solo se acepta un parametro.\n", call. = FALSE)
}

# Configuracion de folders da datos
baseDataFolder <- "../IUCN_data"

# Ruta en la que se encuentra el poligono para hacer M de cada especie
shapePath <- file.path(baseDataFolder, 'shapes/wwf_eco_mesoa.shp')
shapeLayer <- "wwf_eco_mesoa"
regionalizacion <- readOGR(shapePath, shapeLayer)

# Directorio de covariables
covarDataFolder <- file.path(baseDataFolder, "covar_rasters")

inputDataFile <- args[1]
outputFolder <- inputDataFile %>%
  basename %>%
  file_path_sans_ext

if (!dir.exists(outputFolder)) {
  dir.create(outputFolder, recursive = TRUE)
}

####LIMPIEZA DE DUPLICADOS####
occsData <- read.csv(inputDataFile, header = TRUE, stringsAsFactors = FALSE) %>%
  clean_dup("Dec_Long", "Dec_Lat", 0.00833333333)

write.csv(occsData,
          file = file.path(outputFolder, "clean_data.csv"),
          row.names = FALSE)

#### VARIABLES AMBIENTALES####
covarFileList <- list_files_with_exts(covarDataFolder, "tif")
clima <- raster::stack(covarFileList)

#### VARIABLES + PRESENCIAS####
sp_coor <- occsData[c("Dec_Long", "Dec_Lat")]

covarData <- raster::extract(clima, sp_coor)
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
climaImportantes <- clima[[select_var]]

####CALIBRACION####
# Seleccionar el 70 de los datos para calibrar y el resto para validar
sampleDataPoints <- sample.int(
  nrow(covarData),
  size = floor(0.7*nrow(covarData))
)

selectedValues <- rep(0, nrow(covarData)) %>% inset(sampleDataPoints, 1)

covarData$isTrain <- selectedValues

# Ahora cortar los raster con las ecoregiones donde exisan puntos de la especie
coordinates(sp_coor) <- ~Dec_Long+Dec_Lat
proj4string(sp_coor) <- proj4string(regionalizacion)

# obtener los datos de los poligonos que tienen sitios de colecta
dataenpoly <- over(sp_coor, regionalizacion, fn = NULL)
enpolyindex <- which(!is.na(dataenpoly$ECO_NAME))
polydatadf <- dataenpoly[enpolyindex, ]
id_polys <- unique(polydatadf$ECO_NAME)
poligonofilter <- regionalizacion[regionalizacion$ECO_NAME %in% id_polys, ]
# recortamos el grid
# climaImportantesRecortadas <- crop(climaImportantes, poligonofilter)
env <- mask(climaImportantes, poligonofilter) #M de la especie

# MAXENT
# Proceso de modelacion, primero separar los datos de calibracion y validacion
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
resultados_enmeval <- sp@results
write.csv(resultados_enmeval,
          file = file.path(outputFolder, "resultados_enmeval.csv")
          ,row.names = FALSE)
delta_aic <- which(resultados_enmeval$delta.AICc == 0)
#
# ENM EN RASTER
# seleccionar raster del modelo más parsomonioso
predictions <- sp@predictions
modelo.delta.aic <- predictions[[delta_aic]]
writeRaster(modelo.delta.aic,
            file.path(outputFolder, "ENM_in_m.tif"),
            overwrite = TRUE)

# elegir los parametros de maxent del modelo más parsimonioso y
# proyectar a otras variables
modelo.maxent <- sp@models[[delta_aic]]

# Raster del modelo en la M, en escala logistica
mapa.en.m <- predict(modelo.maxent, env)
#plot(mapa.en.m)
writeRaster(mapa.en.m,
            file.path(outputFolder, "ENM_in_mlog.tif"),
            overwrite = TRUE)

##tranferirlo a otro tiempo o area mas grande
mapa.area.grande <- predict(modelo.maxent, climaImportantes)
#plot(mapa.area.grande)
writeRaster(mapa.area.grande,
            file.path(outputFolder, "ENM_log.tif"),
            overwrite = TRUE)


####VALIDACION####
#Independiente de umbral
#AUC
testpp <- extract(mapa.en.m, occs_val)
abs <- extract(mapa.en.m, bg.df)
combined <- c(testpp, abs)       
label <- c(rep(1,length(testpp)),rep(0,length(abs)))  
pred <- prediction(combined, label)   
perf <- performance(pred, "tpr", "fpr")               
auc <- performance(pred, "auc")@y.values[[1]]       
auc                      
write.csv(auc,file = paste0("outputs/",datos$Taxon[1],"_auc.csv"),row.names = FALSE)

#Dependiente de umbral
source("funciones_LAE.R")
#reclasificar mapa de la calibracion
rcl.m <- na.omit(extract(mapa.en.m, occs_cal)) 

mapa.area.grande

#usando el valor de minimo de idoneidad que tienen los puntos de occurencia
rcl.min<-min(rcl.m) # extraer el minimo valor de presencia
mapa.en.m.bin <- reclassify(mapa.en.m, c(-Inf,rcl.min,0,rcl.min,Inf,1)) # reclasificar - cambie su valor en donde esta el valor decimal
writeRaster(mapa.en.m.bin, paste0("outputs/",datos$Taxon[1],"_bin_min.tif"),overwrite=TRUE)
mapa.en.mg.MTPbin <- reclassify(mapa.area.grande, c(-Inf,rcl.min,0,rcl.min,Inf,1)) # reclasificar - cambie su valor en donde esta el valor decimal
writeRaster(mapa.en.mg.MTPbin, paste0("outputs/",datos$Taxon[1],"_binG_MTP.tif"),overwrite=TRUE)
#10 percentil
rcl.10 <- quantile(na.omit(rcl.m),.10)
mapa.en.m.bin10 <- reclassify(mapa.en.m, c(-Inf,rcl.10,0,rcl.10,Inf,1)) # reclasificar - cambie su valor en donde esta el valor decimal
writeRaster(mapa.en.m.bin10, paste0("outputs/",datos$Taxon[1],"_bin_10.tif"),overwrite=TRUE)
mapa.en.mg.bin10 <- reclassify(mapa.area.grande, c(-Inf,rcl.10,0,rcl.10,Inf,1)) # reclasificar - cambie su valor en donde esta el valor decimal
writeRaster(mapa.en.mg.bin10, paste0("outputs/",datos$Taxon[1],"_binG_10.tif"),overwrite=TRUE)
#plot(mapa.en.m.bin)
#
# ##Para validar los modelos binarios usar el codigo que se llama AllMetrics.R.
