####ENMTest####
#source("funciones_LAE.R")
#Threslhold independent

#AUC
aucCalculator <- function(prediction, occs, bgPoints) {
  data <- rbind(occs, setNames(bgPoints, names(occs)))
  labels <- c(rep(1, nrow(occs)),
              rep(0, nrow(bgPoints)))
  scores <- raster::extract(prediction, data)
  pred <- ROCR::prediction(scores, labels)
  # perf <- performance(pred, "tpr", "fpr")
  auc <- performance(pred, "auc")@y.values[[1]]
  return(auc)
}

aucStatistcs <- function(model, models, env, occs, bgPoints) {
  result <- apply(model, 1, function(x, models, env, occs, bgPoints){
    choicedModel <- models[[as.integer(x["index"])]]
    prediction <- dismo::predict(choicedModel, env)
    auc <- aucCalculator(prediction, occs, bgPoints)
    return(c(x["settings"], auc))
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


# Testing background
bg.df.test <- bg.df %>%
  dplyr::filter(isTrain == 0) %>%
  dplyr::select(x, y)

resultsAUC <- aucStatistcs(modelsAIC0, sp.models@models, env, occsValidacion, bg.df.test)
write.csv(resultsAUC,
          file = file.path(outputFolder, "data_auc.csv"),
          row.names = FALSE)

