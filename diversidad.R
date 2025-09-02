### data_cleanUp Script must be run before to save data in the 
### enviroment and have it in the proper format


#nos aseguramos que metadata y otu_columns sigan el mismo orden
metadata <- metadata[match(rownames(otu_columns), metadata$Sample), ]


###DIVERSIDAD ALPHA###
  
richness <- rowSums(otu_columns > 0)

shannon_index <- vegan::diversity(otu_columns, index = "shannon")
simpson_index <- vegan::diversity(otu_columns, index = "simpson")


riqueza <- data_frame(
    ELC = as.factor(metadata$ELC),
    Population = metadata$Poblacion,
    Riqueza = richness,
    Shannon = shannon_index,
    Simpson = simpson_index)


anova(lm(formula = Riqueza ~ ELC, data= riqueza))
anova(lm(formula = Shannon ~ ELC, data= riqueza))
anova(lm(formula = Simpson ~ ELC, data= riqueza))


#Riqueza vs ELC
ggsave(
  filename = paste0("output/", codigo, "/riqueza_", cult_code, codigo, ".png"),
  plot =  ggplot(riqueza, aes(x = ELC, y = Riqueza, fill = ELC)) +
  geom_boxplot() +
  scale_fill_manual(values = c("0" = "black", "1" = "red3", "2" = "yellow3", "3" = "green3"),
                     name = "ELCr") +
  labs(y = "Richness") +
  theme_minimal(),
  width = 7, height = 5, dpi = 300
)

#Riqueza vs Poblacion
ggsave(
  filename = paste0("output/", codigo, "/riqueza_pop_", cult_code, codigo, ".png"),
  plot =  ggplot(riqueza, aes(x = Population, y = Riqueza, color = as.factor(ELC))) +
  geom_point(size = 3) + 
  scale_color_manual(values = c("0" = "black", "1" = "red3", "2" = "yellow3", "3" = "green3"),
                    name = "ELCr") +
  labs(x = "Population", y = "Richness") +
  theme_minimal(),
  width = 7, height = 5, dpi = 300
)

#Shannon vs ELC
ggsave(
  filename = paste0("output/", codigo, "/shannon_", cult_code, codigo, ".png"),
  plot =  ggplot(riqueza, aes(x = ELC, y = Shannon, fill = ELC)) +
  geom_boxplot() +
  scale_fill_manual(values = c("0" = "black", "1" = "red3", "2" = "yellow3", "3" = "green3"),
                    name = "ELCr") +
  labs(x = "ELCr", y = "Shannon Index") +
  theme_minimal(),
  width = 7, height = 5, dpi = 300
)

#Simpson vs ELC
ggsave(
  filename = paste0("output/", codigo, "/simpson_", cult_code, codigo, ".png"),
  plot =  ggplot(riqueza, aes(x = ELC, y = Simpson, fill = ELC)) +
  geom_boxplot() +
  scale_fill_manual(values = c("0" = "black", "1" = "red3", "2" = "yellow3", "3" = "green3"),
                    name = "ELCr") +
  labs(x = "ELCr", y = "Simpson Index") +
  theme_minimal(),
  width = 7, height = 5, dpi = 300
)

#Shannon vs riqueza
ggsave(
  filename = paste0("output/", codigo, "/shannon_vs_riqueza_", cult_code, codigo, ".png"),
  plot = ggplot(riqueza, aes(x = Riqueza, y = Shannon, color = ELC)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("0" = "black", "1" = "red3", "2" = "yellow3", "3" = "green3"),
                     name = "ELCr") +
  labs(x = "Species Richness", y = "Shannon Index") +
  theme_minimal(),
  width = 7, height = 5, dpi = 300
)

#Simpson vs riqueza
ggsave(
  filename = paste0("output/", codigo, "/simpson_vs_riqueza_", cult_code, codigo, ".png"),
  plot = ggplot(riqueza, aes(x = Riqueza, y = Simpson, color = ELC)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("0" = "black", "1" = "red3", "2" = "yellow3", "3" = "green3"),
                     name = "ELCr") +
  labs(x = "Species Richness", y = "Simpson Index") +
  theme_minimal(),
  width = 7, height = 5, dpi = 300
)






###DIVERSIDAD BETA###

###ABUNDANCIAS

###Distance matrix, PCoA, RDA y NMDS

bray_matrix <- vegan::vegdist(otu_columns, method = "bray")


#PCoA
  pcoa <- cmdscale(bray_matrix, eig = TRUE, k = nrow(otu_columns) - 1)
  pcoa_df <- data.frame(Sample = rownames(otu_columns), PCoA1 = pcoa$points[,1], PCoA2 = pcoa$points[,2])
  pcoa_df <- left_join(pcoa_df, metadata, by = "Sample")
  
  ggsave(
    filename = paste0("output/", codigo, "/PCoA_", cult_code, codigo, ".png"),
    plot = ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = as.factor(ELC))) +
    geom_point(size = 3) +
    scale_color_manual(values = c("0" = "black", "1" = "red3", "2" = "yellow3", "3" = "green3"),
                       name = "ELCr") +
    labs(title = "PCoA (Bray-Curtis)", x = "PCoA 1", y = "PCoA 2") +
    theme_minimal(),
    width = 7, height = 5, dpi = 300
  )


#NMDS
  nmds <- vegan::metaMDS(bray_matrix, k = 2)
  nmds$stress
  
  nmds_coords <- as.data.frame(vegan::scores(nmds))  # get NMDS coordinates
  nmds_coords$Sample <- rownames(nmds_coords)  # add sample names as a column
  nmds_coords <- left_join(nmds_coords, metadata, by = "Sample")  # merge with metadata
  
  ggsave(
    filename = paste0("output/", codigo, "/NMDS_", cult_code, codigo, ".png"),
    plot = ggplot(nmds_coords, aes(x = NMDS1, y = NMDS2, color = as.factor(ELC))) +
    geom_point(size = 3) +
    scale_color_manual(values = c("0" = "black", "1" = "red3", "2" = "yellow3", "3" = "green3"),
                       name = "ELCr") +
    labs(title = "NMDS Bray",
         x = "NMDS 1", y = "NMDS 2", color = "ELC") +
    theme_minimal(),
    width = 7, height = 5, dpi = 300
  )


##RDA
  #rda necesita coordenadas en espacios euclidianos, las matrices de distancias
  #podemos transformalas con PCoA

  #con datos crudos
  rda <- list() #limpiamos la lista
  rda <- vegan::rda(otu_columns ~ ELC + Longitud, data = metadata)

  anova(rda, by = "margin", permutations = 9999)
  
  
  rda_df <- as.data.frame(vegan::scores(rda, display = "sites"))
  rda_df <- rda_df %>% mutate(Sample = rownames(rda_df))
  rda_df <- left_join(rda_df, metadata, by = "Sample")
  
  
  betadisp_test <- vegan::betadisper(dist(rda_df[, 1:2]), metadata$ELC)
  anova(betadisp_test)
  
  
  #el circo este de .data[[]], es porque ggplot lo recomienda, 
  #pero con .data no puedo usar indice numerico
  
  ggsave(
    filename = paste0("output/", codigo, "/RDA_ELClat_", cult_code, codigo, ".png"),
    plot = ggplot(rda_df, aes(x = .data[[colnames(rda_df)[1]]], y = .data[[colnames(rda_df)[2]]], 
                              color = as.factor(ELC))) +
      geom_point(size = 3) +
      scale_color_manual(values = c("0" = "black", "1" = "red3", "2" = "yellow3", "3" = "green3"), 
                         name = "ELCr") +
      labs(title = "RDA ~ ELC + Latitud", 
           x = colnames(rda_df)[1], y = colnames(rda_df)[2]) +
      theme_minimal(),
    width = 7, height = 5, dpi = 300
  )
  
  
  
  
  #con matriz bray-PCoA
  db_rda <- list()
  db_rda <- capscale(bray_matrix ~ ELC + Longitud, data = metadata)
  anova(db_rda, by = "terms", permutations = 9999)
  
  
  db_rda_df <- as.data.frame(vegan::scores(db_rda, display = "sites"))
  db_rda_df <- db_rda_df %>% mutate(Sample = rownames(db_rda_df))
  db_rda_df <- left_join(db_rda_df, metadata, by = "Sample")
  
  
  betadisp_test <- vegan::betadisper(dist(db_rda_df[, 1:2]), metadata$ELC)
  anova(betadisp_test)

  
  ggsave(
    filename = paste0("output/", codigo, "/dbRDA_ELC.Long_", cult_code, codigo, ".png"),
    plot = ggplot(db_rda_df, aes(x = .data[[colnames(db_rda_df)[1]]], y = .data[[colnames(db_rda_df)[2]]], 
                     color = as.factor(ELC))) +
    geom_point(size = 3) +
    scale_color_manual(values = c("0" = "black", "1" = "red3", "2" = "yellow3", "3" = "green3"), 
                       name = "ELCr") +
    labs(title = "RDA (Bray-Curtis) ~ ELCr + Longitude", 
         x = colnames(db_rda_df)[1], y = colnames(db_rda_df)[2]) +
    theme_minimal(),
    width = 7, height = 5, dpi = 300
  )



#PERMANOVA y beta dispersion
  
  vegan::adonis2(bray_matrix ~ ELC, 
                 data = metadata, 
                 permutations = 9999)
  
  betadisp_test <- vegan::betadisper(bray_matrix, metadata$ELC)
  anova(betadisp_test)
  
  
  vegan::adonis2(bray_matrix ~ cmi06, 
                 data = metadata, 
                 permutations = 9999)
  
  vegan::adonis2(bray_matrix ~ Poblacion, 
                 data = metadata, 
                 permutations = 9999)
  
  vegan::adonis2(bray_matrix ~ Poblacion + cmi06, 
                 data = metadata, 
                 permutations = 9999) 
  
  vegan::adonis2(bray_matrix ~ cmi06 + Poblacion, 
                 data = metadata, 
                 permutations = 9999)
  
  betadisp_test <- vegan::betadisper(bray_matrix, metadata$cmi06)
  anova(betadisp_test) #lo que estamos mirando es si la dispersión dentro del grupo es sign distinta a la de los otros




### GDM para matriz de distancias y cmi06 (https://cran.r-project.org/web//packages//gdm/gdm.pdf)

#convertir el objeto dist en un objeto matriz
bray_matrix <- as.matrix(bray_matrix)
#añadir columna de Samples (gdm lo necesita)
bray_matrix <- cbind(Sample = as.numeric(rownames(bray_matrix)), bray_matrix) 


#tengo que transformar longitud y latitud en coordenadas (para distancia espacial del gdm)
sf_points <- sf::st_as_sf(metadata, coords = c("Longitud", "Latitud"), crs = 4326) #sistema de coordenadas de GPS
sf_projected <- sf::st_transform(sf_points, crs = 25830) #sistema coordenadas españa peninsular

metadata <- cbind(metadata, sf::st_coordinates(sf_projected))

#estructuramos metadata para que lo pueda usar gdm; necesita que todo sea numérico
metadata_gdm <- data.frame(
  Sample = as.numeric(metadata$Sample),
  cmi06 = as.numeric(metadata$cmi06),
  ELC = as.numeric(metadata$ELC),
  X = as.numeric(metadata$X),
  Y = as.numeric(metadata$Y)
)

#le damos un formato amigable para gdm
site_pair_data <- gdm::formatsitepair(bray_matrix, bioFormat = 3, predData=metadata_gdm, siteColumn = "Sample", XColumn = "X", YColumn = "Y")

gdm_model <- gdm::gdm(site_pair_data, geo = TRUE)

summary(gdm_model)

plot(gdm_model, plot.layout = c(1,2))

gdm_pred <- predict(gdm_model, site_pair_data)  # Predicciones del modelo

plot(site_pair_data$distance, gdm_pred,
     xlab = "Distancia observada",
     ylab = "Distancia predicha",
     main = "Disimilitud Observada vs. Predicha",
     col = "blue", pch = 16)
abline(0, 1, col = "red", lwd = 2, lty = 2)  # Línea de ajuste perfecto








###PRESENCIA/AUSENCIA


#Distancia Jaccard
  jaccard_matrix <- vegan::vegdist(otu_columns_presence, method = "jaccard")

#PCoA  
  pcoa <- cmdscale(jaccard_matrix, eig = TRUE, k = 2)
  pcoa_df <- data.frame(Sample = rownames(otu_columns_presence), PCoA1 = pcoa$points[,1], PCoA2 = pcoa$points[,2])
  pcoa_df <- left_join(pcoa_df, metadata, by = "Sample")
  
  ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = as.factor(ELC))) +
    geom_point(size = 3) +
    scale_color_manual(values = c("0" = "black", "1" = "red3", "2" = "yellow3", "3" = "green3"),
                       name = "ELCr") +
    labs(title = "PCoA (Jaccard)", x = "PCoA 1", y = "PCoA 2") +
    theme_minimal()
  
  
#NMDS  
  nmds <- vegan::metaMDS(jaccard_matrix, k = 2)
  nmds$stress
  
  nmds_coords <- as.data.frame(vegan::scores(nmds))  # get NMDS coordinates
  nmds_coords$Sample <- rownames(nmds_coords)  # add sample names as a column
  nmds_coords <- left_join(nmds_coords, metadata, by = "Sample")  # merge with metadata
  
  ggplot(nmds_coords, aes(x = NMDS1, y = NMDS2, color = as.factor(ELC))) +
    geom_point(size = 3) +
    scale_color_manual(values = c("0" = "black", "1" = "red3", "2" = "yellow3", "3" = "green3"),
                       name = "ELCr") +
    labs(title = "NMDS Jaccard",
         x = "NMDS 1", y = "NMDS 2", color = "ELC") +
    theme_minimal()


##RDA
  
  #con datos crudos
  rda <- vegan::rda(otu_columns_presence ~ ELC, data = metadata)
  
  #con matriz Jaccard-PCoA
  pcoa <- cmdscale(jaccard_matrix, eig = TRUE, k = nrow(otu_columns) - 1)
  pcoa_points <- as.data.frame(pcoa$points)
  rda <- vegan::rda(pcoa_points ~ ELC, data = metadata)
  
  
  #analisis
  anova(rda, by = "terms", permutations = 9999)
  summary(rda)
  
  
  rda_df <- as.data.frame(vegan::scores(rda, display = "sites"))
  rda_df <- rda_df %>% mutate(Sample = rownames(rda_df))
  rda_df <- left_join(rda_df, metadata, by = "Sample")
  
  
  betadisp_test <- vegan::betadisper(dist(rda_df[, 1:2]), metadata$ELC)
  anova(betadisp_test)
  
  
  #el circo este de .data[[]], es porque ggplot lo recomienda, 
  #pero con .data no puedo usar indice numerico
  
  ggplot(rda_df, aes(x = .data[[colnames(rda_df)[1]]], y = .data[[colnames(rda_df)[2]]], 
                     color = as.factor(ELC))) +
    geom_point(size = 3) +
    scale_color_manual(values = c("0" = "black", "1" = "red3", "2" = "yellow3", "3" = "green3"), 
                       name = "ELCr") +
    labs(title = "RDA (Jaccard) ~ ELC", 
         x = colnames(rda_df)[1], y = colnames(rda_df)[2]) +
    theme_minimal()


#PERMANOVA y beta dispersion

  vegan::adonis2(jaccard_matrix ~ ELC, data = metadata, permutations = 9999)
  
  betadisp_test <- vegan::betadisper(jaccard_matrix, metadata$ELC)
  anova(betadisp_test)
  
  vegan::adonis2(jaccard_matrix ~ Poblacion, data = metadata, permutations = 9999)
  
  betadisp_test <- vegan::betadisper(jaccard_matrix, metadata_wild$Poblacion)
  anova(betadisp_test)

  

### GDM para matriz de distancias y cmi06 (https://cran.r-project.org/web//packages//gdm/gdm.pdf)
  
  #convertir el objeto dist en un objeto matriz
  jaccard_matrix <- as.matrix(jaccard_matrix)
  #añadir columna de Samples (gdm lo necesita)
  jaccard_matrix <- cbind(Sample = as.numeric(rownames(jaccard_matrix)), jaccard_matrix) 
  
  
  #tengo que transformar longitud y latitud en coordenadas (para distancia espacial del gdm)
  sf_points <- sf::st_as_sf(metadata, coords = c("Longitud", "Latitud"), crs = 4326) #sistema de coordenadas de GPS
  sf_projected <- sf::st_transform(sf_points, crs = 25830) #sistema coordenadas españa peninsular
  
  metadata <- cbind(metadata, sf::st_coordinates(sf_projected))
  
  #estructuramos metadata para que lo pueda usar gdm; necesita que todo sea numérico
  metadata_gdm <- data.frame(
    Sample = as.numeric(metadata$Sample),
    cmi06 = as.numeric(metadata$cmi06),
    ELC = as.numeric(metadata$ELC),
    X = as.numeric(metadata$X),
    Y = as.numeric(metadata$Y)
  )
  
  #le damos un formato amigable para gdm
  site_pair_data <- gdm::formatsitepair(jaccard_matrix, bioFormat = 3, predData=metadata_gdm, siteColumn = "Sample", XColumn = "X", YColumn = "Y")
  
  gdm_model <- gdm::gdm(site_pair_data, geo = TRUE)
  
  summary(gdm_model)
  
  plot(gdm_model, plot.layout = c(1,2))
  
  gdm_pred <- predict(gdm_model, site_pair_data)  # Predicciones del modelo
  
  plot(site_pair_data$distance, gdm_pred,
       xlab = "Distancia observada",
       ylab = "Distancia predicha",
       main = "Disimilitud Observada vs. Predicha",
       col = "blue", pch = 16)
  abline(0, 1, col = "red", lwd = 2, lty = 2)  # Línea de ajuste perfecto 
  
  
  
  

### OTUs core, shared y exclusive por ELC###

otu_ELC_wild_presence <- aggregate(otu_columns_wild_presence, 
                                   by = list(metadata_wild$ELC), 
                                   FUN = function(x) ifelse(sum(x) > 0, 1, 0))
colnames(otu_ELC_wild_presence)[1] <- "ELC"

otu_ELC_matrix <- as.matrix(otu_ELC_wild_presence[, -1])
rownames(otu_ELC_matrix) <- otu_ELC_wild_presence$ELC

core_OTU <- data.frame(OTU = character(), 
                       stringsAsFactors = FALSE)

shared_OTU <- data.frame(OTU = character(), 
                            ELC = character(), 
                            stringsAsFactors = FALSE)

exclusive_OTU <- data.frame(OTU = character(), 
                               ELC = character(), 
                               stringsAsFactors = FALSE)

for (i in 1:ncol(otu_ELC_matrix)) {
  num <- colSums(otu_ELC_matrix)[i]  # en cuantos ELCs el OTU está presente
  otu_name <- colnames(otu_ELC_matrix)[i]  #OTU name
  elc <- rownames(otu_ELC_matrix)[which(otu_ELC_matrix[, i] == 1)]

  
  if (num == 3) {
    core_OTU <- rbind(core_OTU, 
                      data.frame(OTU = otu_name)) #concatena el nuevo OTU a core_OTU
    }
    else if (num == 2) {
      shared_OTU <- rbind(shared_OTU,
                          data.frame(OTU = otu_name,
                                     ELC = paste(elc, collapse = ", ")))
      }
    else if (num == 1) {
      exclusive_OTU <- rbind(exclusive_OTU,
                          data.frame(OTU = otu_name,
                                     ELC = elc))
      }
}



