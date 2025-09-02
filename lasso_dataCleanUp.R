###LASSO DATA CLEAN-UP###

##este script carga los datos necesarios para que los script de lasso funcionen. 
## Se asume que se ha corrido data_cleanUp antes que este. Sino hay que definir el codigo: codigo <- "XXX_###" #e.g.: sem_16S




##Datos asignación taxonómica

  #fichero de asignación taxonómica
  asv_to_otu <- paste0("asv_to_otu/asv_to_otu_", codigo, ".xlsx")
  
  #datos taxonómicos
  otu_genus <- readxl::read_excel(asv_to_otu) #CAMBIAR NOMBRE
  #con esto seleccionamos la clasificación taxonómica más completa para cada OTU
  otu_genus <- otu_genus %>% select(c(OTU, Kingdom, Phylum, Class, Order, Family, Genus, Species)) %>%
    mutate(num_taxa_filled = rowSums(!is.na(across(Kingdom:Genus)) & across(Kingdom:Genus) != "")) %>%
    group_by(OTU) %>%
    slice_max(num_taxa_filled, with_ties = TRUE) %>% #with_ties: si hay un empate guarda las dos
    unique() %>% #eliminamos si el empate da resultados que son iguales (misma clasificación)
    mutate(
      Species = ifelse(
        length(unique(Species)) > 1,  # Si hay más de una especie, incluyendo si está vacía
        paste0("(", paste(unique(Species[!is.na(Species) & Species != ""]), collapse = " / "), ")"), #si
        unique(Species[!is.na(Species) & Species != ""])  # else: si solo hay una especie, la mantenemos tal cual
      )
    ) %>%
    unique() %>% #resolvemos los empates que se han resuelto uniendo nombres
    ungroup() %>%
    select(-num_taxa_filled) #elimina la columna que hemos creado
  
  #añadimos organización taxonómica de forma jerárquica
  otu_taxa <- data.frame(
    OTU_id = trimws(paste0("OTU", otu_genus$OTU)),
    Taxonomy = ifelse(
      !is.na(otu_genus$Species), paste(otu_genus$Genus, otu_genus$Species),
      ifelse(!is.na(otu_genus$Genus), paste(otu_genus$Genus, "spp."),
             ifelse(!is.na(otu_genus$Family), otu_genus$Family,
                    ifelse(!is.na(otu_genus$Order), otu_genus$Order,
                           ifelse(!is.na(otu_genus$Class), otu_genus$Class,
                                  ifelse(!is.na(otu_genus$Phylum), paste("Unknown", otu_genus$Phylum),
                                         "Unknown bacteria")))))),
    stringsAsFactors = FALSE
  )



#datos fenotípicos
  
  #datos fenotípicos:
  rawData <- readxl::read_excel("C:/Users/carlos.fceldran/OneDrive - Universidad Rey Juan Carlos/2024/Data/MedidasCultivo/MedidasLupinus_2023_2024.xlsx")
  growth <- read.csv("C:/Users/carlos.fceldran/OneDrive - Universidad Rey Juan Carlos/2024/Data/MedidasCultivo/growthRates_lm.csv")
  
  allData <- left_join(rawData, growth, by = c("pop_code", "sample"))
  
  medidas <- allData[, c("pop_code", "sowing","day_till_flower", "RWC", "dry_weight", 
                         "root_length_cm",	"root_collar_mm",	"N_second_roots", 
                         "N_branch",	"N_fruits",	"N_seeds",	"seed_mg", 
                         "shoot_mass",	"root_mass",	"root_shoot", "slope_day")]
  
  medidas <- medidas %>%
    filter(sowing =="si") %>%
    mutate(across(3:ncol(.), as.numeric))
  
  #agrupamos haciendo la media: un solo dato por poblacion
  medidas <- medidas %>% group_by(pop_code) %>% 
    summarise(across(where(is.numeric), mean, na.rm = TRUE))
  
  #seleccionamos solo datos fenotípicos
  vars <- medidas %>% select_if(is.numeric)%>% colnames()
  
  
  
  