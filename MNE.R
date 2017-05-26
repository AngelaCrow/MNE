library("rgdal")
library("raster")
library("fuzzySim")
library("ENMeval")
library("ROCR")

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


setwd("J:/USUARIOS/DarwinUICN/reuniones_talleres/taller lista roja/Resultados_TLR/Mapas_validos/MapasValidosPorGenero")
source("clean_dup.R")

set.seed(1)
#Ruta en la que se encuentra el poligono para hacer M de cada especie
shapePath <- "H:/CoberturasRestringidas/DarwinCLIMA"
shapeFile <- "wwf_eco_mesoa"
regionalizacion <- readOGR(shapePath, shapeFile)

#En esta primera parte se eliminan los registros en una celda de aproximadamene 1km2 (0.00833333333 arcos de seg)
setwd("J:/USUARIOS/DarwinUICN/reuniones_talleres/taller lista roja/Resultados_TLR/Mapas_validos/MapasValidosPorGenero/Capsicum")
dir.create("CleanUp_1km")

####LIMPIEZA DE DUPLICADOS####
sp_files <- list.files("FormatoIUCN",pattern = "*.csv$",full.names = TRUE)
x <- sp_files[[2]]
occs <- read.csv(x,header = T)
print(occs$Binomial[1])
clean_sp <- clean_dup(occs,"Dec_Long","Dec_Lat",threshold = 0.00833333333)
write.csv(clean_sp,file = paste0("CleanUp_1km/",clean_sp$Binomial[1],".csv"),row.names = FALSE)
rm(sp_files,occs)

#### VARIABLES AMBIENTALES####
setwd("H:/CoberturasRestringidas/DarwinCLIMA")
clima<- stack(list.files("Clima",pattern="*.tif$",full.names = TRUE))


#### VARIABLES + PRESENCIAS####
setwd("J:/USUARIOS/DarwinUICN/reuniones_talleres/taller lista roja/Resultados_TLR/Mapas_validos/MapasValidosPorGenero/Capsicum/CleanUp_1km")
dir.create("extract")
#sp_files <- list.files(pattern = "*.csv$",full.names = TRUE)
coor<-c("Dec_Long","Dec_Lat")
sp_coor<-clean_sp[coor]
extract_cl <- extract(clima,sp_coor)
extract_cl<-cbind(clean_sp,extract_cl)
na_varname_index <- which(is.na(extract_cl$bio_1))
if(length(na_varname_index)>0L) extract_cl <- extract_cl[-na_varname_index,]
datos<-extract_cl
rm(na_varname_index, sp_coor)

####SELECCION DE VARIABLES####
dir.create("variables")

print(dim(datos))
# analizar las correlaciones de dos en dos, y combinadas con la significacion sobre la especie:
correlacion<-corSelect(data = datos, sp.cols = 3, var.cols = 23:ncol(datos), cor.thresh = 0.8, use = "pairwise.complete.obs")
#Seleccionar solo las variables no correlacionadas del stack de las variables
select_var<-as.data.frame(correlacion$selected.vars)
vars <- levels(select_var$`correlacion$selected.vars`)
write.csv(vars,file = paste0("variables/",datos$Binomial[1],"_var.csv"),row.names = FALSE)

clima_sp=clima[[vars]] #Variables de la especie sin recortar
rm(select_var, vars)

####CALIBRACION####
# Seleccionar el 70 de los datos para calibrar y el resto para validar
datos$cal_val <- NA
ncalibracion <- floor(dim(datos)[1]*0.7)
datos_cal <- sample(1:dim(datos)[1],size = ncalibracion)
datos$cal_val[datos_cal] <- 1
datos$cal_val[is.na(datos$cal_val)] <- 0

##Ahora cortar los raster con las ecoregiones donde exisan puntos de la especie
coordinates(extract_cl) <- ~Dec_Long+Dec_Lat
proj4string(extract_cl) <- proj4string(regionalizacion)

#obtener los datos de los poligonos que tienen sitios de colecta
dataenpoly<-over(extract_cl, regionalizacion, fn=NULL)
enpolyindex<-which(!is.na(dataenpoly$ECO_NAME))
polydatadf<-dataenpoly[enpolyindex,]
id_polys<-unique(polydatadf$ECO_NAME)
poligonofilter<-regionalizacion[regionalizacion$ECO_NAME %in% id_polys,]
#recortamos el grid
clima_cr<-crop(clima_sp,poligonofilter)
env<-mask(clima_sp, poligonofilter) #M de la especie
rm(clima_cr)

#MAXENT#
#Proceso de modelacion, primero separar los datos de calibracion y validacion
occs_cal <- datos[datos$cal_val==1,]
occs_cal<-occs_cal[,9:10]
occs_cal<-occs_cal[c(2,1)]
occs_val <- datos[datos$cal_val==0,]
occs_val<-occs_val[,9:10]
occs_val<-occs_val[c(2,1)]

#Background
dir.create("outputs")
bg <- randomPoints(env[[1]], n=10000)
bg.df <- as.data.frame(bg)
#plot(env[[1]], legend=FALSE)
#points(bg.df, col='red')
write.csv(bg.df,file = paste0("outputs/",datos$Binomial[1],"_bg.csv"),row.names = FALSE)

#ENMeval
sp <-ENMevaluate(occs_cal, env, bg.df, RMvalues = seq(0.5, 4, 0.5),
                 fc = c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"),
                 method = "randomkfold", kfolds = 4, bin.output = TRUE,
                 parallel = TRUE, numCores = 8)


##identificar el nombre del modelo más parsimonioso
resultados_enmeval<-sp@results
write.csv(resultados_enmeval,file = paste0("outputs/",datos$Binomial[1],"_enmeval.csv"),row.names = FALSE)
delta_aic<-which (resultados_enmeval$delta.AICc == 0)

####ENM EN RASTER####
##seleccionar raster del modelo más parsomonioso
predictions<-sp@predictions
modelo.delta.aic<-predictions[[delta_aic]]
writeRaster(modelo.delta.aic,paste0("outputs/", datos$Binomial[1],"_in_m.tif"),
            overwrite=TRUE)

#elegir los parametros de maxent del modelo más parsimonioso y
#proyectar a otras variables
modelo.maxent<-sp@models[[which (sp@results$delta.AICc == 0) ]]

##Raster del modelo en la M, en escala logistica
mapa.en.m <- predict(modelo.maxent, env)
#plot(mapa.en.m)
writeRaster(mapa.en.m,paste0("outputs/", datos$Binomial[1],"_in_mlog.tif"),overwrite=TRUE)

##tranferirlo a otro tiempo o area mas grande
mapa.area.grande <- predict(modelo.maxent, clima_sp)
#plot(mapa.area.grande)
writeRaster(mapa.area.grande, paste0("outputs/",datos$Binomial[1],"_log.tif"),overwrite=TRUE)


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
write.csv(auc,file = paste0("outputs/",datos$Binomial[1],"_auc.csv"),row.names = FALSE)

#Dependiente de umbral
source("funciones_LAE.R")
#reclasificar mapa de la calibracion
rcl.m <- na.omit(extract(mapa.en.m, occs_cal)) 

#usando el valor de minimo de idoneidad que tienen los puntos de occurencia
rcl.min<-min(rcl.m) # extraer el minimo valor de presencia
mapa.en.m.bin <- reclassify(mapa.en.m, c(-Inf,rcl.min,0,rcl.min,Inf,1)) # reclasificar - cambie su valor en donde esta el valor decimal
writeRaster(mapa.en.m.bin, paste0("outputs/",datos$Binomial[1],"_bin_min.tif"),overwrite=TRUE)
#10 percentil
rcl.10 <- quantile(na.omit(rcl.m),.10)
mapa.en.m.bin10 <- reclassify(mapa.en.m, c(-Inf,rcl.10,0,rcl.10,Inf,1)) # reclasificar - cambie su valor en donde esta el valor decimal
writeRaster(mapa.en.m.bin10, paste0("outputs/",datos$Binomial[1],"_bin_10.tif"),overwrite=TRUE)
#plot(mapa.en.m.bin)

##Para validar los modelos binarios usar el codigo que se llama AllMetrics.R. 
