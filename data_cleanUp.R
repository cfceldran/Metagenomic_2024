#This is the first piece of code that should be run. This code is a template for proper data load and organization in the enviroment, please change it for the 
# specific needs of the user, but maintain object names, as they will be use throught the rest of the codes.

#Requirements
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(vegan)

#input files should be OTU abundance matrixes. Only one should be charge per run of the code.
  #in my example I have 4 distinct files of seed and soil samples

  #codigo or _code are use for dinamic naming of the files, and is useful for saving csv and graphs latter on

  #this acts as a file name code to in case you want to perform the filtering by group, i.e. excluding some of them 
  cult_code <- "_all_" #all refers to all the data

  #16S semillas
  OTUs <- read.csv2("file1.csv", header=T, sep=",")
  codigo <- "file_1" #change them as needed
  
  #16S suelos
  OTUs <- read.csv2("file2.csv", header=T, sep=",")
  codigo <- "file_2" 
  
  #ITS semilla
  OTUs <- read.csv2("file3.csv", header=T, sep=",")
  codigo <- "file_3" 
  
  #ITS suelos
  OTUs <- read.csv2("file4.csv", header=T, sep=",")
  codigo <- "file_4" 


#this is for metadata loading, in my case they were in .xslx files, but the command can be changed for you data.
#Please keep the Sample code in the metadata, its gonna be used in the _cleanUp and posterior analysis.

#METADATA
  meta.data <- as.data.frame(readxl::read_excel("METADATA.xlsx", sheet = 3))
  meta.data <- meta.data[, c("Samples", "Population", "ELC")]                       #this is just to filter the metadata of interest
  colnames(meta.data)[colnames(meta.data) == "Samples"] <- "Sample"                 #the renaming is important for later, this is the code
  

#my analysis was centeredin climatic data, so here I am adding it to metadata. This step is not necesarry if you metadata has all the information you need, otherwise, you can use this to add any data to your metadata
clim <- readxl::read_excel("ClimFile.xlsx", 
                           sheet = "SheetName")
clim <- clim[, c("Acronimo", "Longitud", "Latitud", "cmi06")]                       #filter the climData I want             
metadata <- left_join(meta.data, clim, by = c("Population" = "Acronimo"))           # join with metadata. In my case was population.



### WARNING ###

###this is used to filter the data in case you want to eliminate certain groups for certain analysis###
  
  metadata <- metadata %>% filter(ELC != 0) #metadatos poblaciones naturales

cult_code <- ""                                                                   # as mention before, add the code so you can track the group of data you are using

###

  

#This is to clean the ID-code that the original matrix had, provided by the sequencing company. It is important to rename the ID-code to "Sample".
  #this is very niche of my data, cause metadata sample-code was done only with the numeric part of the ID-code. This step might be useless for you, but just in case.

OTUs$Sample <- stringr::str_extract(OTUs$X, "(?<=_)\\d+")                         #the regular expresion should be addapted. This one extracts the number after '_'
OTUs <- OTUs %>% 
  select(Sample, everything()) %>%                                                #moves Sample col to the begining, eases readability, but not necesary for the code
  select(-X)                                                                      #deletes original code


#data frame only with the abundance information
otu_columns <- OTUs %>% 
  select(-Sample)
otu_columns <- as.data.frame(otu_columns)
rownames(otu_columns) <- OTUs$Sample                                               #add Sample to rownames not to loose the information

#we make otu_columns and metadata coincide (i.e.: we have previously filtered metadata)
otu_columns <- otu_columns %>% 
  filter(rownames(otu_columns) %in% metadata$Sample)





#otu_column presence/absence
otu_columns_presence <- otu_columns
otu_columns_presence[otu_columns_presence > 0] <- 1

#calculate relative abundance 
otu_columns_relative <- otu_columns / rowSums(otu_columns)
otu_columns_log_relative <- log1p(otu_columns_relative)  # log(1 + x)


#in case you have rows that do not have any OTUs present (empty rows), they will have NAs in otu_columns_relative. This is to eliminate those, without having to manually check
  na_rows <- apply(otu_columns_relative, 1, function(x) any(!is.na(x)))
  # Clean otu_columns and reload
  otu_columns <- otu_columns[na_rows, ]
  rm(na_rows)
  # match it with Metadata
  metadata <- metadata[match(rownames(otu_columns), metadata$Sample), ]

#recalculate everything 
                   
  #otu_column presence/absence
  otu_columns_presence <- otu_columns
  otu_columns_presence[otu_columns_presence > 0] <- 1
  
  #calculate relative abundance 
  otu_columns_relative <- otu_columns / rowSums(otu_columns)
  otu_columns_log_relative <- log1p(otu_columns_relative)  # log(1 + x)
    
    
  #calculate frecuency of each OTU
  otu_frequencies <- colSums(otu_columns_relative > 0) / nrow(otu_columns_relative) # >0 devuelve TRUE si hay presencia, colSums cuenta las instancias de TRUE entre numero de OTUS que hay
  


#Filter OTUs that have less than X% de las muestras
otu_filtered <- otu_columns_relative[, otu_frequencies > 0.05] # 5% in this example

  

  
