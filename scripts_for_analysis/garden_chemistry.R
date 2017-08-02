
## LOAD REQUIRED LIBRARIES ----
library(tidyverse)

garden.data <- read_csv('final_data/herb.trait.data.csv') # downloaded from dryad
filter.genos <- c("U","I","T","G","J","X","W","S","F","L")

base <- ggplot(garden.data %>% filter(Genotype %in% filter.genos), aes(x = Genotype))

base + geom_boxplot(aes(y = sal_tannin.PC1))

base + geom_boxplot(aes(y = cinn.PC1))

base + geom_boxplot(aes(y = cinn.PC2))

base + geom_boxplot(aes(y = flavonOLES.PC1))

base + geom_boxplot(aes(y = flavonOLES.PC2))

base + geom_boxplot(aes(y = flavanonOLES.PC1))

chem.data <- garden.data %>% 
  filter(Genotype %in% filter.genos) %>% 
  select(Genotype, C_N_imputed, sal_tannin.PC1:flavanonOLES.PC1) %>%
  group_by(Genotype) %>%
  summarise_all(funs(mean(., na.rm = TRUE)))

# test wind effect
w.f.geno.hell <- betadisper(vegdist(f.hell, method = "euclidean"),  interaction(fungal.df$Genotype, fungal.df$Block), bias.adjust = TRUE)

w.f.geno.centr.hell <- data.frame(w.f.geno.hell$centroids, 
                                   id = rownames(w.f.geno.hell$centroids)) %>%
  separate(col = id, into = c("Genotype","Block")) %>%
  mutate(Block = as.factor(Block),
         Genotype = as.factor(Genotype)) %>%
  left_join(., chem.data)

adonis(w.f.geno.hell$centroids ~ Block + sal_tannin.PC1, 
       data = w.f.geno.centr.hell, 
       method = "euclidean") 

adonis(w.f.geno.hell$centroids ~ cinn.PC1, 
       data = w.f.plots.centr.hell, 
       method = "euclidean") 

adonis(w.f.geno.hell$centroids ~ cinn.PC2, 
       data = w.f.geno.centr.hell, 
       method = "euclidean") 

adonis(w.f.geno.hell$centroids ~ flavonOLES.PC1, 
       data = w.f.geno.centr.hell, 
       method = "euclidean") 

adonis(w.f.geno.hell$centroids ~ flavonOLES.PC2, 
       data = w.f.geno.centr.hell, 
       method = "euclidean") 

adonis(w.f.geno.hell$centroids ~ flavanonOLES.PC1, 
       data = w.f.geno.centr.hell, 
       method = "euclidean") 

fungal.df.chem <- fungal.df %>% left_join(., chem.data)

adonis(f.hell ~ Block + Wind.Exposure + cinn.PC1 + flavonOLES.PC1 + flavanonOLES.PC1 + Genotype, data = fungal.df.chem, method = "euclidean", permutations = how(block = fungal.df$Block, nperm = 999))

adonis(w.f.geno.hell$centroids ~ sal_tannin.PC1, 
       data = w.f.plots.centr.hell, 
       method = "euclidean") 

adonis(w.f.geno.hell$centroids ~ cinn.PC1, 
       data = w.f.plots.centr.hell, 
       method = "euclidean") 

adonis(w.f.geno.hell$centroids ~ cinn.PC2, 
       data = w.f.plots.centr.hell, 
       method = "euclidean") 

adonis(w.f.geno.hell$centroids ~ flavonOLES.PC1, 
       data = w.f.plots.centr.hell, 
       method = "euclidean") 

adonis(w.f.geno.hell$centroids ~ flavonOLES.PC2, 
       data = w.f.plots.centr.hell, 
       method = "euclidean") 

adonis(w.f.geno.hell$centroids ~ flavanonOLES.PC1, 
       data = w.f.plots.centr.hell, 
       method = "euclidean") 