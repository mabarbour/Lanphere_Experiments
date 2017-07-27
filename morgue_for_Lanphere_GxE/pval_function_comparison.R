
dat <- read.csv("example_data.csv") %>%
  mutate(X = as.factor(X),
         year = as.factor(year),
         block = as.factor(block))
glimpse(dat) # everything is specified correctly

mod.mixed <- mixed(response ~ treatment*year*genotype + 
                        (1|block) + 
                        (1|block:treatment) +
                        (1|plant.id),
                      data = dat,
                      contrasts = list(treatment = "contr.sum",
                                       genotype = "contr.sum",
                                       year = "contr.sum"),
                      method = "KR")

mod.lmer <- lme4::lmer(response ~ treatment*year*genotype + 
                    (1|block) + 
                    (1|block:treatment) +
                    (1|plant.id),
                  data = dat,
                  contrasts = list(treatment = "contr.sum",
                                   genotype = "contr.sum",
                                   year = "contr.sum"))#,
                  #method = "KR")

mod.lmerTest <- lmerTest::lmer(response ~ treatment*year*genotype + 
                         (1|block) + 
                         (1|block:treatment) +
                         (1|plant.id),
                       data = dat,
                       contrasts = list(treatment = "contr.sum",
                                        genotype = "contr.sum",
                                        year = "contr.sum"))#,
#method = "KR")

mod.mixed
Anova(mod.lmer, test = "F", type = 3)
anova(mod.lmerTest, ddf = "Kenward-Roger", type = 3)
