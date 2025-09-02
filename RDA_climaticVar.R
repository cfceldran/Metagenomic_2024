library(dplyr)

ecoData <- readxl::read_excel("C:/Users/carlos.fceldran/OneDrive - Universidad Rey Juan Carlos/2024/Data/PoblacionesLupinus.xlsx", sheet = "CHELSA")
metadata_full <- left_join(meta.data, ecoData, by = c("Poblacion" = "Acronimo"))

clima <- ecoData %>% select_if(is.numeric) 
clim_names <- clima %>% colnames() #tiene que haber 118
clim_bio <- clim_names[grep("Che_bio", clim_names)]

bray_matrix <- vegan::vegdist(otu_columns, method = "bray")
pcoa <- cmdscale(bray_matrix, eig = TRUE, k = nrow(otu_columns) - 1)
pcoa_points <- as.data.frame(pcoa$points)

MuMIN_data <- pcoa_points %>%
  tibble::rownames_to_column(var = "Sample") %>%
  left_join(metadata_full, by = "Sample")

for (var in clima_names){
  
  formula <- as.formula(paste(colnames(pcoa_points[1]), "~", var, " + Poblacion"))
  global_model <- lme4::lmer(formula, data = MuMIN_data, na.action = na.fail)
  
  model_selection <- MuMIn::dredge(global_model, rank = "AICc")
  best_model <- MuMIn::get.models(model_selection, 1)[[1]]
  best_aicc <- model_selection$AICc[1]
  
  rda_raw <- vegan::rda(otu_columns ~ var, data = metadata_full)
  rda_matrix <- vegan::rda(pcoa_points ~ var, data = metadata_full)

  }


#para hacer las combinaciones de todas las variables
all_combinations <- unlist(lapply(1:length(clim_names), function(x) combn(clim_names, x, simplify = FALSE)), recursive = FALSE)

#esto es una forma poco elegante de elimnar las de ALG, porque no tienen datos en Chelsa idk why
metadata_full <- metadata_full %>% filter(!Sample %in% 110:112)
otu_columns <- otu_columns[!rownames(otu_columns) %in% as.character(110:112), ]



# Evalúa cada combinación con RDA
rda_results <- lapply(clim_bio, function(predictors) {
  formula <- as.formula(paste("otu_columns ~", predictors))
  model <- vegan::rda(formula, data = metadata_full)
  aicc <- AIC(model, k = log(nrow(metadata_full)))  # AICc manual (suponiendo loglik disponible)
  data.frame(model = paste(predictors, collapse = " + "), AICc = aicc)
})
