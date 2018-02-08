#### Calibration ####
# Divides your data into trainining and test data sets. 70/30 %
sampleDataPoints <- sample.int(
  nrow(covarData),
  size = floor(0.7*nrow(covarData))
)

selectedValues <- rep(0, nrow(covarData)) %>% inset(sampleDataPoints, 1)

covarData$isTrain <- selectedValues
write.csv(cbind(covarData@data, coordinates(covarData)), file.path(outputFolder, "speciesCovarDB.csv"),
          row.names = FALSE)

# MAXENT calibration
# We used ENMeval package to estimate optimal model complexity (Muscarrella et al. 2014)
# Modeling process, first separate the calibration and validation data
occsCalibracion <- covarData %>%
  as.data.frame() %>%
  dplyr::filter(isTrain == 1) %>%
  dplyr::select(Dec_Long, Dec_Lat)

occsValidacion <- covarData %>%
  as.data.frame() %>%
  dplyr::filter(isTrain == 0) %>%
  dplyr::select(Dec_Long, Dec_Lat) 

# Background
bg.df <- dismo::randomPoints(env[[1]], n = 10000) %>% as.data.frame()

#Divide backgeound into train and test 
sample.bg <- sample.int(
  nrow(bg.df),
  size = floor(0.7*nrow(bg.df))
)
selectedValues.bg <- rep(0, nrow(bg.df)) %>% inset(sample.bg, 1)

bg.df$isTrain <- selectedValues.bg

sp::coordinates(bg.df) <- c("x", "y")
sp::proj4string(bg.df) <- crs.wgs84
bg.dfbio <- raster::extract(enviromentalVariables, bg.df) 

bg.df<-as.data.frame(bg.df)
bg.dfbio <- cbind(bg.df, bg.dfbio) %>% as.data.frame()

write.csv(bg.dfbio, file = file.path(outputFolder, "background_points.csv"),
          row.names = FALSE)

#training background
bg.df.cal <- bg.df %>%
  dplyr::filter(isTrain == 1) %>%
  dplyr::select(x, y)


# ENMeval
sp.models <- ENMevaluate(occsCalibracion, env, bg.df.cal, RMvalues = seq(0.5, 4, 0.5),
                         fc = c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"),
                         method = "randomkfold", kfolds = 4, bin.output = TRUE,
                         parallel = TRUE, numCores = parallel::detectCores()-1, 
                         updateProgress = TRUE)

resultados_enmeval <- sp.models@results

write.csv(resultados_enmeval,
          file = file.path(outputFolder, "enmeval_results.csv"),
          row.names = FALSE)

# delta_aic <- which(resultados_enmeval$delta.AICc == 0)
modelsAIC0 <- resultados_enmeval %>%
  mutate(index = rownames(resultados_enmeval)) %>%
  filter(delta.AICc == 0) %>%
  select(index, settings) %>%
  mutate(index = as.numeric(index), settings = as.character(settings))

aic.opt <- sp.models@models[[which(sp.models@results$delta.AICc==0)]]
importa <- var.importance(aic.opt)
write.csv( importa, 
           file = file.path(outputFolder, "varImportance.csv"),
           row.names = FALSE)

# save species niche (raw output) model over raster 
saveRasterWithSettings <- function(models, predictions, prefix) {
  raster::writeRaster(predictions[[models["settings"]]],
                      file.path(outputFolder, paste0(prefix,
                                                     models["settings"],
                                                     ".tif")),
                      overwrite = TRUE)
}

apply(modelsAIC0, 1, saveRasterWithSettings,
      predictions = sp.models@predictions, prefix = "ENM_prediction_M_raw_")

rm(list=ls())
