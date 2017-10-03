# MAXENT calibration
# We used ENMeval package to estimate optimal model complexity (Muscarrella et al. 2014)
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
bg.df <- as.data.frame(bg)

#Divide backgeound into train and test 
sample.bg <- sample.int(
  nrow(bg.df),
  size = floor(0.7*nrow(bg.df))
)
selectedValues.bg <- rep(0, nrow(bg.df)) %>% inset(sample.bg, 1)

bg.df$isTrain <- selectedValues.bg
write.csv(bg.df, file = file.path(outputFolder, "background_data.csv"),
          row.names = FALSE)

#training background
bg.df.cal <- bg.df %>%
  dplyr::filter(isTrain == 1) %>%
  dplyr::select(x, y)


# ENMeval
sp <- ENMevaluate(occsCalibracion, env, bg.df.cal, RMvalues = seq(0.5, 4, 0.5),
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