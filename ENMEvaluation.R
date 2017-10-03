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
