library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(vegan)


#Cargar solo el que nos interesa

  cult_code <- "_all_"

  #16S semillas
  OTUs <- read.csv2("OTUs/OTUs_16S_semillas.csv", header=T, sep=",")
  codigo <- "sem_16S" 
  
  #16S suelos
  OTUs <- read.csv2("OTUs/OTUs_16S_suelo.csv", header=T, sep=",")
  codigo <- "suel_16S" 
  
  #ITS semilla
  OTUs <- read.csv2("OTUs/OTUs_ITS_semillas.csv", header=T, sep=",")
  codigo <- "sem_ITS" 
  
  #ITS suelos
  OTUs <- read.csv2("OTUs/OTUs_ITS_suelo.csv", header=T, sep=",")
  codigo <- "suel_ITS" 


#METADATA SEMILLAS
  meta.data <- as.data.frame(readxl::read_excel("tax_otus_metadatos_semillas_17_01_25.xlsx", sheet = 3))
  meta.data <- meta.data[, c("Samples", "Poblacion", "ELC")]
  colnames(meta.data)[colnames(meta.data) == "Samples"] <- "Sample"
  
  
  

#METADATA SUELOS
  meta.data <- as.data.frame(readxl::read_excel("tax_otus_metadatos_suelos_20_01_25.xlsx", sheet = 1))
  meta.data <- meta.data[, c("Codigo Illumina", "Poblacion", "ELC")]
  colnames(meta.data)[colnames(meta.data) == "Codigo Illumina"] <- "Sample"
  


  
  
clim <- readxl::read_excel("C:/Users/carlos.fceldran/OneDrive - Universidad Rey Juan Carlos/2024/Data/PoblacionesLupinus.xlsx", 
                           sheet = "PobsTotal")
clim <- clim[, c("Acronimo", "Longitud", "Latitud", "cmi06")]
clim$Latitud <- as.numeric(clim$Latitud)
clim$Longitud <- as.numeric(clim$Longitud)

metadata <- left_join(meta.data, clim, by = c("Poblacion" = "Acronimo"))
metadata$cmi06 <- as.numeric(metadata$cmi06)


metadata$ELC <- metadata$ELC %>% 
  tidyr::replace_na(0) #cultivares ELC=0


###CUIDADO###
###ELIMINAR CULTIVARES###
  #para solo utilizar Wild Accession Data: solo hay en semillas
  metadata <- metadata %>% filter(ELC != 0) #metadatos poblaciones naturales
  cult_code <- ""
###

  

#cambiamos código por Sample
OTUs$Sample <- stringr::str_extract(OTUs$X, "(?<=_)\\d+") #extrae num de samples
OTUs <- OTUs %>% 
  select(Sample, everything()) %>%    #mueve samples al principio, facilita redibilidad
  select(-X)                           #elimina columna código


#data frame solo información taxonómica
otu_columns <- OTUs %>% 
  select(-Sample)
otu_columns <- as.data.frame(otu_columns)
rownames(otu_columns) <- OTUs$Sample #añadimos sample a rownames para no perder la información

#hacemos coincidir OTUs con metadata (i.e.: hemos eliminado cultivares)
otu_columns <- otu_columns %>% 
  filter(rownames(otu_columns) %in% metadata$Sample)


#PARA eliminar OTU 2 y 4 (muy mayoritarios) otu_columns <- otu_columns[, -c(2,4)]



#otu_column presencia/ausencia
otu_columns_presence <- otu_columns
otu_columns_presence[otu_columns_presence > 0] <- 1

#calcular abundancia relativa
otu_columns_relative <- otu_columns / rowSums(otu_columns)
otu_columns_log_relative <- log1p(otu_columns_relative)  # log(1 + x)


#en algunas salen filas de todo 0s, por lo que salen como NAs en otu_columns_relative; recomiendo quitarlas de todas
  na_rows <- apply(otu_columns_relative, 1, function(x) any(!is.na(x)))
  # Limpiarlas de otu_columns y volver a cargar el resto
  otu_columns <- otu_columns[na_rows, ]
  rm(na_rows)
  # eliminamos datos de metadata también
  metadata <- metadata[match(rownames(otu_columns), metadata$Sample), ]


#otu_column presencia/ausencia
otu_columns_presence <- otu_columns
otu_columns_presence[otu_columns_presence > 0] <- 1

#calcular abundancia relativa
otu_columns_relative <- otu_columns / rowSums(otu_columns)
otu_columns_log_relative <- log1p(otu_columns_relative)  # log(1 + x)  
  
  
#calcular la frecuencia de cada OTU
otu_frequencies <- colSums(otu_columns_relative > 0) / nrow(otu_columns_relative) # >0 devuelve TRUE si hay presencia, colSums cuenta las instancias de TRUE entre numero de OTUS que hay



#Filtrar OTUs que aparecen en menos del X% de las muestras
otu_filtered <- otu_columns_relative[, otu_frequencies > 0.05]






#No lo uso, pero aquí esta por si acaso:
  #join info in the same data.frame
  all <- left_join(metadata, OTUs, by = "Samples")
  all$ELC <- as.factor(all$ELC)
  all$Poblacion <- as.factor(all$Poblacion)
  
  