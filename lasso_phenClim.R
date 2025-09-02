### data_cleanUp and lasso_dataCleanUp Scripts must be run before to save data 
### in the enviroment and have it in the proper format
## si se ha corrido el de diversidad probablemente hay que refrescar metadata, porque el GDM mete más columnas


library(glmnet)

#Output es una tabla de los rasgos y los OTUs asociados
output <- paste0("output/", codigo, "/", codigo, "_otu_phenClim.xlsx")

#output_long es la misma tabla, desglosada y con información taxonómica
output_long <- paste0("output/", codigo, "/", codigo, "_otu_phenClim_long.xlsx")


#arreglamos variable dependiente (OTUs*cmi)
  X <- as.matrix(otu_columns_relative)   # Predictors
  
  X <- scale(X)                          #estandarizamos datos
  X <- X[, colSums(is.na(X)) == 0]
  
  cmi <- scale(metadata$cmi06)
  
  X_interaction <- sweep(X, 1, cmi, `*`)  # cada columna de X multiplicada por CMI
  
  X_interaction <- scale(X_interaction)
  

#arreglamos variable independiente (datos fenotipicos)
  
  #con esto ordenamos medidas por Sample
  phen <- left_join(metadata, medidas, by = c("Poblacion" = "pop_code"))
  
  #estandarizamos datos
  for (col in names(phen)){
    if (col %in% vars){
      phen[[col]] <- scale(phen[[col]])
      phen[[col]] <- as.vector(phen[[col]]) #scale devuelve una matriz, lo reconvertimos en un vector
    }}

  
# Fit Lasso with cross-validation, get best lambda
#se hacen iteraciones porque los resultados no son estables

  min.lambdas <- list() #limpiamos lista
  
  for (var in vars){
    lambda_temp <- numeric(100)
    for(i in 1:100){
      lambda_temp[i] <- cv.glmnet(y= phen[[var]], 
                                  x= X_interaction,
                                  type.measure = c("deviance"), 
                                  nfolds = 7,
                                  alpha = 1,
                                  family="gaussian")$lambda.min
      #gaussian para cmi(continua), multinomial para ELC (factor)
      
    }
    print(var)
    min.lambdas[[var]] <- lambda_temp
    print(mean(min.lambdas[[var]]))
  }



# Fit Lasso model
  lasso_phen <- list()
  
  otu_phen <- data.frame(
    Trait = character(),
    OTU = c()
  )
  
  otu_phen_long <- data.frame(
    Trait = character(),
    OTU = c(),
    Taxonomy = character(),
    Coefficient_CMI = numeric(),
    Coefficient = numeric()
  )
  
  for (var in vars){
    lasso_phen[[var]] <- glmnet(y = phen[[var]], 
                                x = X_interaction,
                                family = "gaussian",
                                nlambda = 1000, #Default
                                alpha = 1,
                                lambda = mean(min.lambdas[[var]])) #opción conservadora
    
    selected_otus <- as.data.frame(as.matrix(coef(lasso_phen[[var]]))[-1, , drop = FALSE]) #eliminamos el intercepto
    selected_otus <- selected_otus %>% 
      filter(s0 != 0) 
    selected_otus_names <- selected_otus %>% rownames()
    
    # Si no hay OTUs seleccionados, pasar a la siguiente variable
    if (nrow(selected_otus) == 0) next
    
    coef_cmi <- selected_otus$s0
    
    coef <- c()
    r2 <- c()
      
      for (otu_name in selected_otus_names) {
        
        # Crear la variable de interacción
        X_single <- X[, otu_name, drop = FALSE]
        
        # Ajustar el modelo de interacción para ese OTU y el trait
        fit <- lm(phen[[var]] ~ X_single)
        coef_temp <- coef(fit)[-1]  # Eliminamos el intercepto
        r2_temp <- summary(fit)$r.squared
        
        coef <- c(coef, coef_temp)
        r2 <- c(r2, r2_temp)
      }
    
    
    #rellenamos data.frames
    otu_phen <- rbind(otu_phen, data.frame(
      Trait = var,
      OTU = paste(selected_otus_names, collapse = ", "),
      stringsAsFactors = FALSE))
    
    
    taxonomy <- otu_taxa$Taxonomy[match(selected_otus_names, otu_taxa$OTU_id)]   #encuentra el taxa del OTU
    
    
    otu_phen_long <- rbind(otu_phen_long, data.frame(
      Trait = var,
      OTU = selected_otus_names,
      Taxonomy = taxonomy,
      Coefficient_CMI = coef_cmi,
      Coefficient = coef,
      R2 = r2,
      stringsAsFactors = FALSE))
  }


openxlsx::write.xlsx(otu_phen, file = output, rowNames = FALSE)
openxlsx::write.xlsx(otu_phen_long, file = output_long, rowNames = FALSE)
