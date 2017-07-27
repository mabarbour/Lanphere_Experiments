
# upload data
wind_plant_traits_field_summer_2012 <- read.csv("~/Documents/Lanphere_Experiments/data/wind_plant_traits_field_summer_2012.csv")

# change NAs to zeros, because I'm virtually positive they reflect real zeros. Need to triple check though.
wind_plant_traits_field_summer_2012[ ,8:9][is.na(wind_plant_traits_field_summer_2012[ ,8:9])] <- 0

# Change NAs to zeros, for mature and immature shoots
wind_plant_traits_field_summer_2012[ ,12:28][is.na(wind_plant_traits_field_summer_2012[ ,12:28])] <- 0

# create new variables to sum up mature and immature shoot length
wind_plant_traits_field_summer_2012$mature_shoot_total <- with(wind_plant_traits_field_summer_2012, X1.shoot.M + X2.shoot.M + X3.shoot.M + X4.shoot.M + X5.shoot.M + X6.shoot.M + X7.shoot.M + X8.shoot.M + X9.shoot.M + X10.shoot.M)

wind_plant_traits_field_summer_2012$immature_shoot_total <- with(wind_plant_traits_field_summer_2012, X1.shoot.I + X2.shoot.I + X3.shoot.I + X4.shoot.I + X5.shoot.I + X6.shoot.I)

# create unique plant code to match other data sets
wind_plant_traits_field_summer_2012$plant_code <- with(wind_plant_traits_field_summer_2012, paste(Wind.Exposure, Block, Plant.Position, sep="_"))

# focal plant data set
sub_plant_traits <- subset(wind_plant_traits_field_summer_2012, select=c("Mature.Green.Leaves","Shoots.sprouting.at.base", "plant_code", "mature_shoot_total", "immature_shoot_total"))

