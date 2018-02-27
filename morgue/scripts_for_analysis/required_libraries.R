
## LIBRARIES ----

library(merTools) # must be loaded before dplyr, otherwise 'select()' function is masked

library(tidyverse) # easy dataframe manipulation

library(lme4) # analyzing mixed-effect models

library(effects) # calculating mean and confidence intervals of treatment and genotype effects.

library(psych) # correlation tests

library(vegan) # analyzing community data

library(cowplot) # souped up version of ggplot

library(RCurl) # for pulling code from other github repos

devtools::install_github("gavinsimpson/ggvegan") # run if 
library(ggvegan) # for using fortify on rda objects

## FUNCTIONS ----

eval(parse(text = getURL("https://raw.githubusercontent.com/mabarbour/miscellaneous_R/master/model_diagnostic_functions.R", ssl.verifypeer = FALSE)))

eval(parse(text = getURL("https://raw.githubusercontent.com/mabarbour/miscellaneous_R/master/autoplot.custom.R", ssl.verifypeer = FALSE)))

eval(parse(text = getURL("https://raw.githubusercontent.com/mabarbour/miscellaneous_R/master/veganCovEllipse.R", ssl.verifypeer = FALSE)))


