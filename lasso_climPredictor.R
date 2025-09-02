### data_cleanUp Script must be run before to save data in the 
### enviroment and have it in the proper format
## si se ha corrido el de diversidad probablemente hay que refrescar metadata, porque el GDM mete más columnas


library(glmnet)


#Output es una tabla de los OTU asociados a cada ELC y al cmi en general
output <- paste0("output/", codigo, "/", codigo, "_otu_clim.xlsx")

#output_long es la misma tabla, desglosada y con información taxonómica
output_long <- paste0("output/", codigo, "/", codigo, "_otu_clim_long.xlsx")


#datos dependientes (OTUs)
  X <- as.matrix(otu_columns_relative)   # Predictors
  #estandarizamos datos
  X <- scale(X)
  X <- X[, colSums(is.na(X)) == 0]

  
#cmi06
  Y <- scale(metadata$cmi06)           # Response
  

  # Fit Lasso with cross-validation, get best lambda
  #se hacen iteraciones porque los resultados no son establess
  
  min.lambdas <- c()
  
  for(i in 1:100){
    min.lambdas[i] <- cv.glmnet(y=Y, 
                                x=X,
                                type.measure = c("deviance"), 
                                nfolds = 7,
                                alpha = 1,
                                family="gaussian")$lambda.min
                                #gaussian para cmi(continua), multinomial para ELC (factor)
    print(i)
  }
  
  print(mean(min.lambdas))
  
  # Fit Lasso model
  lasso_cmi <- glmnet(y=Y, 
                        x=X,
                        family = "gaussian",
                        nlambda = 1000, #Default
                        alpha = 1,
                        lambda = mean(min.lambdas)) #opción conservadora
  
  # Check selected OTUs
  selected_cmi <- as.data.frame(as.matrix(coef(lasso_cmi)))
  selected_cmi_names <- selected_cmi %>% 
    filter(selected_cmi$s0 != 0) %>% rownames()
  selected_cmi_names <- setdiff(selected_cmi_names, "(Intercept)")

  
   
  
  
#ELC
  Y <- metadata$ELC           # Response
  
  
  # Fit Lasso with cross-validation, get best lambda
  #esto se hace porque los resultados no son estables
  
  min.lambdas <- c()
  
  for(i in 1:100){
    min.lambdas[i] <- cv.glmnet(y=Y, 
                                x=X,
                                type.measure = c("deviance"), 
                                nfolds = 7,
                                alpha = 1,
                                family="multinomial")$lambda.min
                                #gaussian para cmi(continua), multinomial para ELC (factor)
    
    print(i)
  }
  
  print(mean(min.lambdas))
  
  # Fit Lasso model
  lasso_elc <- glmnet(y=Y, 
                      x=X,
                      family = "multinomial",
                      nlambda = 1000, #Default
                      alpha = 1,
                      lambda = mean(min.lambdas)) #opción conservadora
  
  # Check selected OTUs
  selected_elc <- data.frame(
    ELC1 = lasso_elc$beta[1]$'1'[,1],
    ELC2 = lasso_elc$beta[2]$'2'[,1],
    ELC3 = lasso_elc$beta[3]$'3'[,1])

  elc_otus <- data.frame(
    ELC = character(),
    OTU = character()
  )
  
  elc_otus_long <- data.frame(
    ELC = character(),
    OTU = character(),
    Taxonomy = character(),
    Coefficient = numeric()
  )
  
  for (elc in colnames(selected_elc)){
    
    selected_elc_names <- rownames(selected_elc)[selected_elc[[elc]] != 0]
    
    elc_otus <- rbind(elc_otus, data.frame(
      ELC = elc,
      OTU = paste(selected_elc_names, collapse = ", "),
      stringsAsFactors = FALSE
    ))
    
    for (otu in selected_elc_names){
      
      taxonomy <- otu_taxa$Taxonomy[match(otu, otu_taxa$OTU_id)]   #encuentra el taxa del OTU
      
      elc_otus_long <- rbind(elc_otus_long, data.frame(
        ELC = elc,
        OTU = otu,
        Taxonomy = taxonomy,
        Coefficient = selected_elc[otu, elc],
        stringsAsFactors = FALSE
    ))}
  }
  
#para contruir una sola tabla que tenga también el cmi
  elc_otus <- rbind(elc_otus, data.frame(
    ELC = "cmi",
    OTU = paste(selected_cmi_names, collapse = ", "),
    stringsAsFactors = FALSE
  ))
  
  for (otu in selected_cmi_names){
    taxonomy <- otu_taxa$Taxonomy[match(otu, otu_taxa$OTU_id)]   #encuentra el taxa del OTU
    
    elc_otus_long <- rbind(elc_otus_long, data.frame(
      ELC = "cmi",
      OTU = otu,
      Taxonomy = taxonomy,
      Coefficient = selected_cmi[otu, 's0'],
      stringsAsFactors = FALSE
    ))
  }


openxlsx::write.xlsx(elc_otus, file = output, rowNames = FALSE)
openxlsx::write.xlsx(elc_otus_long, file = output_long, rowNames = FALSE)

