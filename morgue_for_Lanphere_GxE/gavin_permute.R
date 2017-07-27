## load required libraries ----
library(vegan)
library(dplyr)
library(lme4)
library(car)

####################################################
## My main question for you is whether I'm using the permutation package correctly to test my hypotheses of interest.
####################################################

###################################################
## I conducted an experiment to examine how the composition of mycorrhizal communities varies in response to plant genotype, the environment (wind exposure), and the interaction between the two (GxE).

## My experiment consists of 10 blocks. Nested within each block are two wind exposure treatments (E = 'exposed' and U = 'unexposed'). Within each wind treatment, I randomly planted the same 10 willow genotypes (F, G, I, J, L, S, T, U, W, X), although there was some mortality resulting in an unbalanced design
###################################################

## upload sample data (500 OTUs) and visualize my split-plot experimental design ----
gav <- read.csv("gavin.csv") %>%
  tbl_df() %>%
  mutate(Block = as.factor(Block)) 

# I feel this table illustrates my experimental design
with(gav, table(Wind.Exposure, Genotype, Block)) # appears to be some dups...# note that the data is unbalanced, which hopefully isn't too much of a problem...

## Univariate Response ----
# I would use a mixed-effect model to test for an effect of genotype, environment, and GxE interaction on a single OTU. 
# I've included Wind.Exposure nested within Block to model my split-splot design.
uni <- lmer(OTU_1 ~ Wind.Exposure*Genotype + (1|Block/Wind.Exposure), data = gav, contrasts = list(Wind.Exposure = "contr.sum", Genotype = "contr.sum"))
summary(uni)
Anova(uni, test = "F", type = 3) # note that the residual degrees of freedom for Wind.Exposure is close to 9, which is what I would expect based on my experimental design.

#########################################################
## I want to use your permute package to mimic the hypothesis test above by accounting for the effects of block and wind exposure nested within block. Note that although I used a type-3 sum of squares above, I don't think this is an option in rda() or adonis() so I'm happy with testing my hypotheses based on type-2 sum of squares for those analyses (i.e. marginal effects).
## Below, I show you what I've come up with so far.
#########################################################

## Multivariate Response: Method 1 ----
f.comm <- decostand(select(gav, -(X:Genotype)), method = "hellinger") # get and transform community data for RDA

# create permutation design
samp <- nrow(f.comm)
block.perm <- how(block = gav$Block)
shuffle(samp, block.perm) # Permuting at the block level makes sense to me for evaluating at least the GxE effect. This is because the effect of wind exposure and genotype are decoupled from each other within each block of the experiment. I'm concerned though whether this is appropriate for testing the main effects of Genotype and Wind.Exposure.

# test GxE and evaluate homogeneity of variance assumption
mult1.GxE <- rda(f.comm ~ Wind.Exposure*Genotype, data = gav)
anova(mult.GxE, by = "margin", permutations = block.perm) 

disper.GxE <- betadisper(d = vegdist(f.comm, method = "euclidean"), group = interaction(gav$Genotype, gav$Wind.Exposure), bias.adjust = TRUE)
boxplot(disper.GxE)
anova(disper.GxE) # seems to satisfy assumption of homogeneity of variance
with(gav, table(Wind.Exposure, Genotype)) # although clearly some sample sizes are pretty low.

# test main effects of G and E and evaluate homogeneity of variance assumption
mult1.main <- rda(f.comm ~ Genotype + Wind.Exposure, data = gav)
anova(mult1.main, by = "margin", permutations = block.perm)

disper.G <- betadisper(d = vegdist(f.comm, method = "euclidean"), group = gav$Genotype, bias.adjust = TRUE)
boxplot(disper.G)
permutest(disper.G, permutations = block.perm) # seems to satisfy assumption of homogeneity of variance

disper.E <- betadisper(d = vegdist(f.comm, method = "euclidean"), group = gav$Wind.Exposure, bias.adjust = TRUE)
boxplot(disper.E)
permutest(disper.E, permutations = block.perm) # seems to satisfy assumption of homogeneity of variance

# I think permuting at the block level is appropriate for testing the GxE effect, and it may even be okay for testing the effect of Genotype since I've included Wind.Exposure in the model. However, I don't think it is okay for testing the effect of Wind.Exposure. This is because I believe the F statistic for Wind.Exposure is calculated based on a residual degrees of freedom of 117, rather than 9 as was used in the mixed-effects model. Therefore, this the currest test of the wind-effect has an inflated Type-1 error due to pseudoreplication (I think).

## Testing for Wind Effect: Method 2 ----
# summarize data at the plot level. This should control for the pseudoreplication I'm worried about in Method 1.
f.wind <- gav %>%
  select(-X, -Genotype) %>%
  group_by(Block, Wind.Exposure, Plot_code) %>%
  summarise_each(funs(mean)) %>%
  ungroup()
f.wind.comm <- decostand(select(f.wind, -(Block:Plot_code)), method = "hellinger")

# create permutation design
samp2 <- nrow(f.wind.comm)
block2.perm <- how(block = f.wind$Block)
shuffle(samp2, block2.perm) # Permuting at the block level makes sense to me for evaluating the effect of wind exposure in the data summarized at the plot level.

# test main effect of wind exposure and evaluate homogeneity of variance assumption.
mult2.main <- rda(f.wind.comm ~ Wind.Exposure, data = f.wind)
anova(mult2.main, by = "margin", permutations = block2.perm)

disper2.E <- betadisper(d = vegdist(f.wind.comm, method = "euclidean"), group = f.wind$Wind.Exposure, bias.adjust = TRUE)
boxplot(disper2.E)
permutest(disper2.E, permutations = block2.perm) # seems to satisfy assumption of homogeneity of variance


############################################################
## Any advice???
############################################################