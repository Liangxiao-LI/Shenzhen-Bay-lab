##Read Package

install.packages("readxl")

library(readxl)

setwd("/Users/ryan/Documents/GitHub/Shenzhen Bay lab/spike_protein_3structures_markerd_with_variants_from_cncb")

data <- read_excel("6vxx_variants.xls")

library(ggplot2)






##Apply Mourn's index

install.packages("gdal")
install.packages("spdep")
install.packages("rgdal")

library(spdep)
library(rgdal)

install.packages("sp")

library(sp)

install.packages("spatialreg")

library(spatialreg)
