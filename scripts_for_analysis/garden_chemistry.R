

## LOAD & MANAGE DATA ----
garden.data <- read_csv('final_data/herb.trait.data.csv') # downloaded from dryad
filter.genos <- c("U","I","T","G","J","X","W","S","F","L") # only use genotypes from Lanphere Experiment

chem.data <- garden.data %>% 
  filter(Genotype %in% filter.genos) %>% 
  select(Genotype, C_N_imputed, sal_tannin.PC1:flavanonOLES.PC1) %>%
  group_by(Genotype) %>%
  summarise_all(funs(mean(., na.rm = TRUE)))
