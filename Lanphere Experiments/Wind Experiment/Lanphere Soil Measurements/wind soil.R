### load packages
library("reshape")

#### upload data

# soil moisture
wind.sep.18.2012 <- read.csv("~/Dropbox/Willow Project for sharing/Wind Experiment/Lanphere Dunes - Soil Moisture-Temp-EC - Sep 18 2012.csv")

# total organic matter
windTOM <- read.csv("~/Dropbox/Willow Project for sharing/Wind Experiment/Wind - Soil Total Organic Matter.csv")

# manage soil moisture data
agg1 <- aggregate(cbind(moisture,temp,EC)~Site+Treatment,FUN=mean,wind.sep.18.2012)

agg2 <- aggregate(cbind(moisture.1,temp.1,EC.1)~Site+Treatment,FUN=mean,wind.sep.18.2012)

agg3 <- aggregate(cbind(moisture.2,temp.2,EC.2)~Site+Treatment,FUN=mean,wind.sep.18.2012)

total <- cbind.data.frame(agg1,agg2,agg3)
total$avgMoist <- with(total,(moisture+moisture.1+moisture.2)/3)

total$avgTemp <- with(total, (temp+temp.1+temp.2)/3)

total$avgEC <- with(total, (EC+EC.1+EC.2)/3)

plot(avgMoist~Treatment,total)
plot(avgTemp~Treatment,total)
plot(avgEC~Treatment,total)

with(subtotal, interaction.plot(Treatment,Site,avgMoist, fun=mean, col=1:10))

subtotal <- subset(total, select = c("Site","Treatment","avgMoist","avgTemp","avgEC"))
msubtotal <- melt(subtotal, id=c("Site","Treatment"))
treat <- cast(msubtotal, Site~variable+Treatment,mean)

### data suggests that there is no difference in soil moisture, temperature, or EC between wind exposed and unexposed sites.
t.test(treat$avgMoist_exposed,treat$avgMoist_unexposed,paired=T)
t.test(treat$avgTemp_exposed,treat$avgTemp_unexposed,paired=T)
t.test(treat$avgEC_exposed,treat$avgEC_unexposed,paired=T)

### total organic matter analysis
windTOM$OvenWt <- windTOM$Cruc....Oven.Dry.Wt. - windTOM$Crucible.Wt.
windTOM$IgniteWt <- windTOM$Cruc....Ignited.Wt. - windTOM$Crucible.Wt.

windTOM$PercTOM <- with(windTOM, (OvenWt - IgniteWt)/OvenWt)

TOMsub <- subset(windTOM, BAD.Sample == "n")

plot(PercTOM~Wind.Treatment, TOMsub)

TOMmelt <- melt(TOMsub, id.vars=c("Wind.Site.No.","Wind.Treatment"), measure.vars="PercTOM")
TOMcast <- cast(TOMmelt, Wind.Site.No.~Wind.Treatment+variable)

with(TOMcast, t.test(exposed_PercTOM, unexposed_PercTOM, Tpaired=T))

interaction.plot(TOMsub$Wind.Treatment,TOMsub$Wind.Site.No., TOMsub$PercTOM, col=1:10) # note that sites seven and 1 don't show up because I dropped samples from them.  Wind does not appear to influence percent organic matter in soil.