library(vegan)

chelsa <- readxl::read_excel("C:/Users/carlos.fceldran/OneDrive - Universidad Rey Juan Carlos/2024/Data/PoblacionesLupinus.xlsx", 
                            sheet = "CHELSA")
chelsa <- chelsa %>% select(-Poblacion, -Region, -Latitud, -Longitud)

allMetadata <- left_join(metadata, chelsa, by = c("Poblacion" = "Acronimo")) %>% na.omit(.)
#nos aseguramos que todas son numéricas
allMetadata$ELC <- as.factor(allMetadata$ELC)
allMetadata[, 4:ncol(allMetadata)] <- lapply(allMetadata[, 4:ncol(allMetadata)], as.numeric)



##RDA

#esto es necesario porque no tenemos los datos de chelsa de ALG
rda_otu_columns <- otu_columns[match(allMetadata$Sample, rownames(otu_columns)), ]

#con datos crudos
rda <- list() #limpiamos la lista

vars <- paste(colnames(allMetadata)[4:ncol(allMetadata)], collapse = " + ") 
formula <- as.formula(paste("rda_otu_columns ~ ", vars))
                      
rda <- vegan::rda(formula, data = allMetadata)
rda_null <- vegan::rda(rda_otu_columns ~ 1, data = allMetadata)



#reescribo la formula por la óptima que me devuelve ordiR2step
formula <- formula(vegan::ordiR2step(rda_null, scope = formula(rda), direction = "both", permutations = 9999))
rda <- vegan::rda(formula, data = allMetadata)

anova(rda, by = "terms", permutations = 9999)
anova(rda, by = "margin", permutations = 9999)


rda_df <- as.data.frame(vegan::scores(rda, display = "sites"))
rda_df <- rda_df %>% mutate(Sample = rownames(rda_df))
rda_df <- left_join(rda_df, metadata, by = "Sample")


betadisp_test <- vegan::betadisper(dist(rda_df[, 1:2]), allMetadata$ELC)
anova(betadisp_test)


#el circo este de .data[[]], es porque ggplot lo recomienda, 
#pero con .data no puedo usar indice numerico

ggsave(
  filename = paste0("output/", codigo, "/RDA_", cult_code, codigo, ".png"),
  plot = ggplot(rda_df, aes(x = .data[[colnames(rda_df)[1]]], y = .data[[colnames(rda_df)[2]]], 
                            color = as.factor(ELC))) +
    geom_point(size = 3) +
    scale_color_manual(values = c("0" = "black", "1" = "red3", "2" = "yellow3", "3" = "green3"), 
                       name = "ELCr") +
    labs(title = "RDA ~ ELC", 
         x = colnames(rda_df)[1], y = colnames(rda_df)[2]) +
    theme_minimal(),
  width = 7, height = 5, dpi = 300
)




#con matriz bray-PCoA
bray_matrix <- vegan::vegdist(rda_otu_columns, method = "bray")
pcoa <- cmdscale(bray_matrix, eig = TRUE, k = nrow(rda_otu_columns) - 1)
pcoa_points <- as.data.frame(pcoa$points)

formula <- as.formula(paste("pcoa_points ~ ", vars))
rda <- vegan::rda(formula, data = allMetadata)

#analisis
anova(rda, by = "terms", permutations = 9999)


rda_df <- as.data.frame(vegan::scores(rda, display = "sites"))
rda_df <- rda_df %>% mutate(Sample = rownames(rda_df))
rda_df <- left_join(rda_df, metadata, by = "Sample")


betadisp_test <- vegan::betadisper(dist(rda_df[, 1:2]), metadata$ELC)
anova(betadisp_test)


ggsave(
  filename = paste0("output/", codigo, "/RDA_PCoA_", cult_code, codigo, ".png"),
  plot = ggplot(rda_df, aes(x = .data[[colnames(rda_df)[1]]], y = .data[[colnames(rda_df)[2]]], 
                            color = as.factor(ELC))) +
    geom_point(size = 3) +
    scale_color_manual(values = c("0" = "black", "1" = "red3", "2" = "yellow3", "3" = "green3"), 
                       name = "ELCr") +
    labs(title = "RDA (Bray-Curtis) ~ ELC", 
         x = colnames(rda_df)[1], y = colnames(rda_df)[2]) +
    theme_minimal(),
  width = 7, height = 5, dpi = 300
)
