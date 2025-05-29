# This script selects the tree species 

species.vect <- readRDS("data/DE_BWI3_regeneration_species.rds")

spec.omit <- c("DE_BWI.other.abies",
               "DE_BWI.other.coniferous",
               "DE_BWI.other.deciduous",
               "DE_BWI.other.picea",
               "DE_BWI.other.pinus",
               "Sorbus.domestica" # no regeneration occurence
               )

species.vect <- species.vect[!species.vect %in% spec.omit] 

