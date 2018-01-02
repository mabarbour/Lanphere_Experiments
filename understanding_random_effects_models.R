
## COPIED FROM: http://anythingbutrbitrary.blogspot.ch/2012/06/random-regression-coefficients-using.html


rm(list = ls())
set.seed(5432)

J <- 15
N <- 30

test.df <- data.frame( unit = sort(rep(c(1:N),J)), 
                       J = rep(c(1:J),N) , x = rnorm(n = J*N) )

beta <- 3 + .2*rnorm(N)
test.df$beta <- beta[test.df$unit]
test.df$y <- 1 + test.df$x * test.df$beta + .75*rnorm(n = J*N)
head(test.df, 18)

beta.hat <- list()
for(i in 1:N){
  unit.lm <- lm(y ~ x, data = subset(test.df, unit == i) )
  beta.hat[i] <- coef(unit.lm)[2]
}
beta.hat <- as.numeric(beta.hat)

par(mfrow = c(2, 1))
hist(beta, main = "Histogram of Actual Slopes", col = "blue",
     xlab = expression(beta[i]), cex.axis = 1.5, cex.lab = 1.5,
     breaks = seq(from = 2.4, to = 3.6, by = .1) )
hist(as.numeric(beta.hat), main = "Histogram of Estimated Slopes",
     xlab = expression(hat(beta)[i]), col = "blue", cex.axis = 1.5,
     cex.lab = 1.5, breaks = seq(from = 2.4, to = 3.6, by = .1))

library(lme4)
re.lm <- lmer(y ~ x + (1+x|unit), data = test.df) 
summary(re.lm)

re.lm <- lmer(y ~ x + (0+x|unit), data = test.df) 
summary(re.lm)

coef(re.lm)

par(mfrow = c(1,1))
plot(re.beta ~ beta.hat, xlim = c(2.4, 3.6), ylim = c(2.4,3.6), 
     pch = 20, col = "blue",  cex = 1.6, cex.axis = 1.2,
     cex.lab = 1.5, xlab = expression(hat(beta)[i]), 
     ylab =  "", 
     main = "Posterior modes versus least squares estimates")  
title(ylab =  expression(paste("mode of ", hat(italic(p)), 
                               "(", beta[i], "|", bold(y)[i], ",", bold(x)[i],")"),
                         col.main = "red" ) , cex.lab = 1.5, mgp = c(2.4,0,0))
abline(h = 3.0, v = 3.0, col = "red", lty = 2)
abline(h = 3.05594, col = "grey", lty = 3)

re.mse <- mean( (re.beta - beta)^2 )
lm.mse <- mean( (beta.hat - beta)^2 )

cat("Mean Squared Error for first pass linear models:", lm.mse, 
    "\n                 versus random effects modeling:", re.mse)