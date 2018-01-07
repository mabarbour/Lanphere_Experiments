## LOAD LIBRARIES ----
library(tidyverse)
library(lme4)


## TRAIT DATA ----
w.trait.df <- read.csv('final_data/wind_trait_df.csv') %>% tbl_df() %>% mutate(Block = as.factor(Block), Plot_code=paste(Block, Wind.Exposure, sep="_"))
w.trait.2012 <- filter(w.trait.df, Year=="2012")

aa.trait.df <- read.csv('final_data/ant_aphid_trait_df.csv') %>% tbl_df() %>% mutate(Block=as.factor(Block), Ant.mound.dist=as.factor(Ant.mound.dist), Plot_code=paste(Block, Ant.mound.dist, sep="_"))
aa.trait.2012 <- filter(aa.trait.df, Year=="2012") %>% mutate(num.Ant.mound.dist = as.numeric(as.character(Ant.mound.dist)),
                                                              c.num.Ant.mound.dist = num.Ant.mound.dist - mean(num.Ant.mound.dist))
aa.trait.2012$Aphid.treatment <- relevel(aa.trait.2012$Aphid.treatment, "none")
aa.trait.2012$Aphid.treatment <- C(aa.trait.2012$Aphid.treatment, "contr.sum")

## VISUALIZE ----
interaction.plot(x.factor = w.trait.2012$Wind.Exposure, 
                 trace.factor = w.trait.2012$Genotype,
                 response = w.trait.2012$trait.PC1)

interaction.plot(x.factor = aa.trait.2012$Aphid.treatment, 
                 trace.factor = aa.trait.2012$Genotype,
                 response = aa.trait.2012$trait.PC1)

interaction.plot(x.factor = aa.trait.2012$Ant.mound.dist, 
                 trace.factor = aa.trait.2012$Genotype,
                 response = aa.trait.2012$trait.PC1)

interaction.plot(x.factor = aa.trait.2012$num.Ant.mound.dist, 
                 trace.factor = aa.trait.2012$Genotype,
                 response = aa.trait.2012$trait.PC1)

interaction.plot(x.factor = filter(aa.trait.2012, Aphid.treatment=="aphid")$num.Ant.mound.dist, 
                 trace.factor = filter(aa.trait.2012, Aphid.treatment=="aphid")$Genotype,
                 response = filter(aa.trait.2012, Aphid.treatment=="aphid")$trait.PC1)
interaction.plot(x.factor = filter(aa.trait.2012, Aphid.treatment=="none")$num.Ant.mound.dist, 
                 trace.factor = filter(aa.trait.2012, Aphid.treatment=="none")$Genotype,
                 response = filter(aa.trait.2012, Aphid.treatment=="none")$trait.PC1)

## RANDOM SLOPE MODEL ----
rs.mod <- lmer(trait.PC1 ~ Wind.Exposure + (1+Wind.Exposure|Genotype) + (1|Block) + (1|Plot_code), data=w.trait.2012)
summary(rs.mod)
FE.var <- var(as.numeric(fixef(rs.mod) %*% t(model.matrix(rs.mod))))

as.numeric(FE.var/(FE.var + VarCorr(rs.mod)$Plot_code + VarCorr(rs.mod)$Block + VarCorr(rs.mod)$Genotype[1] + VarCorr(rs.mod)$Genotype[4] + sigma(rs.mod)))
as.numeric(VarCorr(rs.mod)$Genotype[1]/(FE.var + VarCorr(rs.mod)$Plot_code + VarCorr(rs.mod)$Block + VarCorr(rs.mod)$Genotype[1] + VarCorr(rs.mod)$Genotype[4] + sigma(rs.mod)))
as.numeric(VarCorr(rs.mod)$Genotype[4]/(FE.var + VarCorr(rs.mod)$Plot_code + VarCorr(rs.mod)$Block + VarCorr(rs.mod)$Genotype[1] + VarCorr(rs.mod)$Genotype[4] + sigma(rs.mod)))

rs.aa.mod <- lmer(trait.PC1 ~ Aphid.treatment*num.Ant.mound.dist + (1+Aphid.treatment*num.Ant.mound.dist|Genotype) + (1|Block) + (1|Plot_code), data=aa.trait.2012)
summary(rs.aa.mod)
VarCorr(rs.aa.mod)$Genotype

## CHARACTER STATE MODEL ----
cs.mod <- lmer(trait.PC1 ~ Wind.Exposure + (0+Wind.Exposure|Genotype) + (1|Block) + (1|Plot_code), data=w.trait.2012)
summary(cs.mod)
FE.var <- var(as.numeric(fixef(cs.mod) %*% t(model.matrix(cs.mod))))

# should I be including the variance in the other environment into these estimates?
as.numeric(FE.var/(FE.var + VarCorr(cs.mod)$Plot_code + VarCorr(cs.mod)$Block + VarCorr(cs.mod)$Genotype[1] + VarCorr(cs.mod)$Genotype[4] + sigma(cs.mod)))
as.numeric(VarCorr(cs.mod)$Genotype[1]/(FE.var + VarCorr(cs.mod)$Plot_code + VarCorr(cs.mod)$Block + VarCorr(cs.mod)$Genotype[1] + VarCorr(cs.mod)$Genotype[4] + sigma(cs.mod)))
as.numeric(VarCorr(cs.mod)$Genotype[4]/(FE.var + VarCorr(cs.mod)$Plot_code + VarCorr(cs.mod)$Block + VarCorr(cs.mod)$Genotype[1] + VarCorr(cs.mod)$Genotype[4] + sigma(cs.mod)))


cs.aa.mod <- lmer(trait.PC1 ~ Aphid.treatment*c.num.Ant.mound.dist + (0+Aphid.treatment*c.num.Ant.mound.dist|Genotype) + (1|Block) + (1|Plot_code), data=aa.trait.2012)
summary(cs.aa.mod)

cs.aa.mod.alt <- lmer(trait.PC1 ~ Aphid.treatment*Ant.mound.dist + (0+Aphid.treatment*Ant.mound.dist|Genotype) + (1|Block) + (1|Plot_code), data=aa.trait.2012)
summary(cs.aa.mod.alt)


aa.Fixed.aphid <- fixef(cs.aa.mod)[2] * t(model.matrix(cs.aa.mod))[2, ]
aa.Fixed.ant <- fixef(cs.aa.mod)[3] * t(model.matrix(cs.aa.mod))[3, ]
aa.Fixed.AxA <- fixef(cs.aa.mod)[4] * t(model.matrix(cs.aa.mod))[4, ]
VarF <- var(aa.Fixed.aphid + aa.Fixed.ant + aa.Fixed.AxA); VarF
varFover <- var(aa.Fixed.aphid) + var(aa.Fixed.ant) + var(aa.Fixed.AxA)
varFover-VarF # perfect partition since it's a balanced experiment...I thought there was some mortality though...

VarF <- var(as.vector(fixef(cs.aa.mod) %*% t(model.matrix(cs.aa.mod)))); VarF # An alternative way for getting the same result

## NESTED MODEL ----
n.mod <- lmer(trait.PC1 ~ Wind.Exposure + (1|Genotype) + (1|Genotype:Wind.Exposure) + (1|Block) + (1|Plot_code), data=w.trait.2012)
summary(n.mod)
FE.var <- var(as.numeric(fixef(n.mod) %*% t(model.matrix(n.mod))))

as.numeric(FE.var/(FE.var + VarCorr(n.mod)$Plot_code + VarCorr(n.mod)$Block + VarCorr(n.mod)$Genotype + VarCorr(n.mod)$'Genotype:Wind.Exposure' + sigma(n.mod)))
as.numeric(VarCorr(n.mod)$Genotype/(FE.var + VarCorr(n.mod)$Plot_code + VarCorr(n.mod)$Block + VarCorr(n.mod)$Genotype + VarCorr(n.mod)$'Genotype:Wind.Exposure' + sigma(n.mod)))
as.numeric(VarCorr(n.mod)$'Genotype:Wind.Exposure'/(FE.var + VarCorr(n.mod)$Plot_code + VarCorr(n.mod)$Block + VarCorr(n.mod)$Genotype + VarCorr(n.mod)$'Genotype:Wind.Exposure' + sigma(n.mod)))


## COMPARE MODELS ----
# pretty much all give the same predictions
plot(predict(rs.mod) ~ predict(n.mod))
plot(predict(rs.mod) ~ predict(cs.mod))
cor.test(predict(rs.mod),predict(cs.mod))
plot(predict(cs.mod) ~ predict(n.mod))

## After further exploration, I think the character state model is the best way to go.
## I can explicitly measure the variance in genotype in two different environments, as well
## as in the slopes between different environments (e.g. c.Ant.mound.dist). With the bayesian
## approach, I don't have to worry about partitioning, but just reporting SD +/- intervals
## to illustrate the confidence in different estimates and the magnitude of these different effects.
## Also, with the character state approach, I can explicitly calculate the correlations in 
## genotypic effects across the different environments (with estimates of error), which will 
## further give insight to the nature of GxE effects (i.e. genotype ranks preserved across environments
## or do they change?). If I want to report variance explained in traits though, I need to account
## for the covariance (I think). I should look at the Bommer R code again for the example of the
## repeatability calculation. I should also look at Nakagawa and Schielzeth to see if they give
## more advice on adjusted vs. unadjusted repeatability.
## Or, if I use a random slope model, then I can still just report the SD +/- intervals. If I
## do this though, then I should use contr.sum for all of my factors and center the distance from ant mound 
## effect. That way, I can calculate the average effect of genotype across environments and the GxE effect
## will be indicated by the variance. Actually, since the predictions are exactly the same
## (at least for GLMER), then I can just use the transformation (like recommended by Paul Johnson)
## to get insight to the different processes at hand! So they're the same! Each has their advantages.
## I think the random slope model is better for the macro picture, but I can use the character
## state approach if I want to calculate broad-sense heritabilities in each environment and the correlations
## amongst environments. And I can partition most of the variance in my fixed effects (at least with the balanced data)
