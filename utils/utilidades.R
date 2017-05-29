#' Funcion para crear un archivo sin registros duplicados
#'
#' @param csvFilePath Ruta del archivo de registros en formato csv
#' @param outCsvFilePath Ruta del archivo donde se guardaran los resultados
#' @param longColumnName Nombre de la columna donde se encuentra la longitud
#' @param latColumnName Nombre de la columna donde se encuentra la latitud
#' @param threshold Máxima distancia entre puntos para que se consideren
#'                  duplicados
#' @return Ruta del archivo con los registros limpios
CreateCleanData <- function(csvFilePath, outCsvFilePath = "",
                            longColumnName, latColumnName,
                            threshold = 0.00833333333) {
  data <- read.csv(csvFilePath, header = TRUE)
  data_clean <- clean_dup(data, longColumnName, latColumnName, threshold)

  if (outCsvFilePath == "") {
    outCsvFilePath <- paste0(tools::file_path_sans_ext(csvFilePath),
                             "_clean.csv")
  }

  write.csv(data_clean, file = outCsvFilePath, row.names = FALSE)

  return(outCsvFilePath)
}