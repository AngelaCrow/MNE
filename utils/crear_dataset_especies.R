dataFilePath <- '../IUCN_data/data/Master_ValidIUCN.csv'
outputDataFolder <- '../IUCN_data/data_especies'

especieData <- data.table::fread(dataFilePath, encoding = "Latin-1")
especieData[,
            data.table::fwrite(
              .SD,
              file.path(outputDataFolder,
                paste0("data_",
                       gsub('([[:punct:]])|\\s+','_', Binomial),
                       ".csv"))),
            by = .(Binomial)]