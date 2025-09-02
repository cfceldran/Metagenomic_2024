Code original of Carlos Celdrán at cfceldran@gmail.com or carlos.fceldran@urjc.es

NOTE 1: this pipeline has been made as general as possible so it can be usniversally used, but was built for metagenomics of plant populations and analysis with climatic data, thus it might have certain specificities explain by its original use.

NOTE 2: main explanations of the codes are included in annotations within the code. But ongoing editing

_cleanUp files are code for uploading and organazing the data in the enviroment in a way that the rest of the code will run smoothly. These are the files that should be changed to fit your data, just note that object names should be maintained for the code to run smoothly. The rest of the files do not need be edited (except obviously to adapt to user necessities).
  These files are unrunable, but serve as an example for the user.

diversidad.R is code for general alpha and beta diversity of the populations.
  it includes: richness, Simpson and Shannon indixes; distant matrix Analysis, PCoA, NMDS, RDA, dbRDA and GDM.

lasso_ files convay code for Lasso Regression analysis of the communities as predictors and several response variables.
