

#source('~/Documents/miscellaneous_R/model_diagnostic_functions.R')
library(merTools) # must be loaded before dplyr
library(dplyr)
library(tidyr)
#library(reshape2)
#library(ggplot2)
library(lme4)
library(effects) # calculating mean and confidence intervals of treatment and genotype effects.
library(psych) # for correlation tests
#library(car) # for Anova function
#library(broom) # for tidying up model outputs
library(RCurl)
script <- getURL("https://raw.githubusercontent.com/mabarbour/miscellaneous_R/master/model_diagnostic_functions.R", ssl.verifypeer = FALSE)
eval(parse(text = script))

#source('~/Documents/miscellaneous_R/model_diagnostic_functions.R')
library(merTools) # must load before dplyr since this package requires 'MASS' which requires 'plyr'
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(psych)
library(lme4)
library(car)
library(vegan)
library(effects)
library(ggplot2)
library(cowplot)

library(RCurl) # for loading github files directly.

script <- getURL("https://raw.githubusercontent.com/mabarbour/miscellaneous_R/master/autoplot.custom.R", ssl.verifypeer = FALSE)
eval(parse(text = script))

script <- getURL("https://raw.githubusercontent.com/mabarbour/miscellaneous_R/master/veganCovEllipse.R", ssl.verifypeer = FALSE)
eval(parse(text = script))

library(RCurl)
script <- getURL("https://raw.githubusercontent.com/mabarbour/miscellaneous_R/master/model_diagnostic_functions.R", ssl.verifypeer = FALSE)
eval(parse(text = script))
script <- getURL("https://raw.githubusercontent.com/mabarbour/miscellaneous_R/master/veganCovEllipse.R", ssl.verifypeer = FALSE)
eval(parse(text = script))
#script <- getURL("https://raw.githubusercontent.com/mabarbour/miscellaneous_R/master/beta.div.R", ssl.verifypeer = FALSE)
#eval(parse(text = script))
#source('~/Documents/miscellaneous_R/model_diagnostic_functions.R')
#source('~/Documents/miscellaneous_R/veganCovEllipse.R')
#source('~/Documents/miscellaneous_R/beta.div.R')
library(merTools) # must load before dplyr since this package requires 'MASS' which requires 'plyr'
library(effects)
library(dplyr)
library(tidyr)
library(vegan)
library(lme4)
library(ggplot2)
library(cowplot) # throwing an error for unknown reason


