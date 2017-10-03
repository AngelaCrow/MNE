#### Projections ####
# predict choicemodel over current climate variables
predictAndSave <- function(model, models, data, prefix, occs) {
  choicedModel <- models[[as.integer(model["index"])]]
  predictions <- dismo::predict(choicedModel, data)
  raster::writeRaster(predictions,
                      file.path(outputFolder, paste0(prefix,
                                                     "log_",
                                                     model["settings"],
                                                     ".tif")),
                      overwrite = TRUE)
  
#Threshold prection using minimum traning (min) and 10 percentil (q10) values  
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
      models = sp@models, data = env, prefix = "ENM_",
      occs = occsCalibracion)

apply(modelsAIC0, 1, predictAndSave,
      models = sp@models, data = env_fc45, prefix = "ENM_fc45",
      occs = occsCalibracion)

apply(modelsAIC0, 1, predictAndSave,
      models = sp@models, data = env_fc85, prefix = "ENM_fc85",
      occs = occsCalibracion)

apply(modelsAIC0, 1, predictAndSave,
      models = sp@models, data = env_fl45, prefix = "ENM_fl45",
      occs = occsCalibracion)

apply(modelsAIC0, 1, predictAndSave,
      models = sp@models, data = env_fl85, prefix = "ENM_fl85",
      occs = occsCalibracion)