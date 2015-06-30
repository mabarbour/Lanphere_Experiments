# Some metadata and notes.
# For wind measurements taken on Sep 10, 2012 between 3 and 4pm, the wind in Arcata was projected to be about 14mph around that time (if I remember correctly). Either way, it was a a very windy day and represents near the maximum wind gusting experienced by these plants. Wnd measurements were made with an anenometer (unknown make and model) in kilometers per hour at a height of about 37 cm aboveground.  Max and min measurement readings were taken over a 30 sec. period and haphazardly started on either an exposed or unexposed site first.

## upload data
speed <- read.csv('~/Documents/Lanphere_Experiments/wind_experiment/data/wind_speed_data_2013.csv') %>%
  tbl_df() %>%
  mutate(block = as.factor(block))

speed1 <- filter(speed, date == "September 10, 2012")
speed2 <- filter(speed, date == "before July 7, 2013")
speed3 <- filter(speed, date == "July 28, 2013")

speed1.lmer <- lmer(max.speed ~ treatment + (1|block), speed1)
summary(speed1.lmer)
Anova(speed1.lmer, test.statistic = "F")

speed2.lmer <- lmer(max.speed ~ treatment + (1|block), speed2)
summary(speed2.lmer)
Anova(speed2.lmer, test.statistic = "F")

# marginal effect. paired t-test gives the exact same results. May just through this one out, especially since it was done by a different person (Brendan).
speed3.lmer <- lmer(max.speed ~ treatment + (1|block), speed3)
summary(speed3.lmer)
Anova(speed3.lmer, test.statistic = "F") # marginal effect is surprsing...
