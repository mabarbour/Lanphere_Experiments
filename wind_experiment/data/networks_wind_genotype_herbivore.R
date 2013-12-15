# set working directory
setwd("~/Documents/Lanphere_Experiments/data")
library("betalink")
source('wind_arthropod_analysis_script.R')
source('~/Documents/GenotypeNetworks/Network 28 betalink functions.r')

# run for loops to extract genotype-herbivore networks for each block in unexposed and exposed sites
herbivores_only <- subset(all_arthropods, select = c("leaf_hoppers_blue_brown","Leaf.Hoppers.Camo","Light.Green.Aphids", "stem_gall", "volcano_gall", "mite_gall", "yellow_knife_leaf_hopper", "twisty_gall", "giant_willow_aphid", "psyllidae", "sawflies", "LH.nymph.other", "aphid_species_used_in_ant_aphid_experiment", "all_Caloptilia", "all_tortricid",  "Genotype", "plant_code", "Block", "Wind.Exposure", "Plant.Position", "Dead."))
herbivores_only$Genotype <- as.character(herbivores_only$Genotype)

exposed_genotype_herbivore <- list()
for(i in 1:10) {
  exposed_genotype_herbivore[[i]] <- as.matrix(subset(herbivores_only, Block == i & Wind.Exposure == "Exposed", select= leaf_hoppers_blue_brown:Genotype))
  dimnames(exposed_genotype_herbivore[[i]])[[1]] <- exposed_genotype_herbivore[[i]][ ,16]
  exposed_genotype_herbivore[[i]] <- exposed_genotype_herbivore[[i]][ ,1:15]
  class(exposed_genotype_herbivore[[i]]) <- "numeric"
}

unexposed_genotype_herbivore <- list()
for(i in 1:10) {
  unexposed_genotype_herbivore[[i]] <- as.matrix(subset(herbivores_only, Block == i & Wind.Exposure == "Unexposed", select= leaf_hoppers_blue_brown:Genotype))
  dimnames(unexposed_genotype_herbivore[[i]])[[1]] <- unexposed_genotype_herbivore[[i]][ ,16]
  unexposed_genotype_herbivore[[i]] <- unexposed_genotype_herbivore[[i]][ ,1:15]
  class(unexposed_genotype_herbivore[[i]]) <- "numeric"
}

# combine all of the different networks together
realizations_genotype_herbivore_wind <- c(exposed_genotype_herbivore, unexposed_genotype_herbivore)
e1 <- exposed_genotype_herbivore[[1]]
class(e1) <- "array"
u1 <- unexposed_genotype_herbivore[[1]]
class(u1) <- "array"

# exploratory analysis, done with "Network 28 betalink functions" 
wind_metaweb <- aggregate.metaweb(realizations_genotype_herbivore_wind, F)
#beta.os_prime(realizations_genotype_herbivore_wind, bf=) # apparently there is a warning when I using the first beta diversity measure, but this measure appears to have the most variation in dissimilarity...
wind_dist <- betalink.dist(realizations_genotype_herbivore_wind, triangular=T, bf="canberra") # warnings that I'm using a binary matrix, but I don't know how to code it otherwise.

exposed <- rep("Exposed",10)
unexposed <- rep("Unexposed",10)
wind_exposure <- as.factor(c(exposed,unexposed))

anosim(wind_dist$WN, wind_exposure)

betalink.plot(realizations_genotype_herbivore_wind[[1]], realizations_genotype_herbivore_wind[[13]])




# older analysis comparing exposed and unexposed communities
exposed_agg_herb <- subset(agg_herb, Wind.Exposure=="Exposed")
mat_exposed_agg_herb <- as.matrix(exposed_agg_herb[ ,3:17])
dimnames(mat_exposed_agg_herb)[[1]] <- as.character(exposed_agg_herb[ ,1])

unexposed_agg_herb <- subset(agg_herb, Wind.Exposure == "Unexposed")
mat_unexposed_agg_herb <- as.matrix(unexposed_agg_herb[ ,3:17])
dimnames(mat_unexposed_agg_herb)[[1]] <- as.character(unexposed_agg_herb[ ,1])

visweb(mat_exposed_agg_herb, type="diagonal")
visweb(mat_unexposed_agg_herb, type="diagonal")
