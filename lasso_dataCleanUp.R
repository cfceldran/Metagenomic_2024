###LASSO DATA CLEAN-UP###

## this Scrpt loads pertinent data for LASSO analysis. 
## data_cleanUp should had been run befor this one. Otherwise code should be defined: codigo <- "XXX_###" #e.g.: sem_16S




##Taxonomic assignation

  #load ASV taxonomic and OTU assignation
  asv_to_otu <- paste0("asv_to_otu/asv_to_otu_", codigo, ".xlsx") #in my case I have the taxonomic assignation of ASV samed with the 'codigo'
  otu_genus <- readxl::read_excel(asv_to_otu) 
  
  #This code will asign each OTU with a taxon extracted from the ASVs that form the OTU. It will assign the most complete taxonomic asignation.
  # as species classification can be variable (in case there is one), we are maintainig the specie assignation but mark it between parentesis to clarify uncertanty.
  
otu_genus <- otu_genus %>% select(c(OTU, Kingdom, Phylum, Class, Order, Family, Genus, Species)) %>%                
    mutate(num_taxa_filled = rowSums(!is.na(across(Kingdom:Genus)) & across(Kingdom:Genus) != "")) %>%               # num_taxa_filled counts in what dept tax assign is. e.g.: if Genus is assigned then num_taxa_filled = 6
    group_by(OTU) %>%
    slice_max(num_taxa_filled, with_ties = TRUE) %>%                                                                 # with_ties= TRUE: will save both if there is a tie
    unique() %>%                                                                                                     # delete if the tie is the same taxonomic assignation
    mutate(
      Species = ifelse(
        length(unique(Species)) > 1,                                                                                  # if more than 1 spp class, includying if there is empty ones
        paste0("(", paste(unique(Species[!is.na(Species) & Species != ""]), collapse = " / "), ")"), 
        unique(Species[!is.na(Species) & Species != ""])                                                              # else: if spp class is consistent, its maintained.
      )
    ) %>%
    unique() %>%                                                                                                      # resolve new ties
    ungroup() %>%
    select(-num_taxa_filled)                                                                                          # delete temporal column

  #add OTU assigment jeraquically
  otu_taxa <- data.frame(
    OTU_id = trimws(paste0("OTU", otu_genus$OTU)),
    Taxonomy = ifelse(
      !is.na(otu_genus$Species), paste(otu_genus$Genus, otu_genus$Species),
      ifelse(!is.na(otu_genus$Genus), paste(otu_genus$Genus, "spp."),
             ifelse(!is.na(otu_genus$Family), otu_genus$Family,
                    ifelse(!is.na(otu_genus$Order), otu_genus$Order,
                           ifelse(!is.na(otu_genus$Class), otu_genus$Class,
                                  ifelse(!is.na(otu_genus$Phylum), otu_genus$Phylum,
                                         "Unknown bacteria")))))),
    stringsAsFactors = FALSE
  )




#Response Variable

  #In my case I am using OTU abundances as the predictive variable.
    
    #Load Response variable. In my case: phenotypic data
    medidas <- readxl::read_excel("phenotipicData.xslx")
    
    #You need to have one response variable per sample. In my case each sample had several individuals, thus several measurements of the same trait, so I made the average of the trait.
    medidas <- medidas %>% group_by(pop_code) %>% 
      summarise(across(where(is.numeric), mean, na.rm = TRUE))
    
    #here I make an object that holds the colnames of the phenotypic traits
    vars <- medidas %>% select_if(is.numeric)%>% colnames()
  
  
  

  
