# This R code examplifies how one can run a Bayesian Random Regression (RR) analysis on data collected by taking 
# repeated measures of individuals (subjects) in different environmental conditions.
# The posteriors are used to calculate the posterior distribution of the between-individual variance in each environment
# and the cross-environmental corrrelation from the RR estimates of variance in elevation and slope and their covariance.
# The example considers the simple case of two environments, but the approach is general
# Also included is a character state implementation of the same problem as well as an example on 
# how one can plot the data with the expected RR distribution on top in order to evaluate the model output.
# The code can be extended to consider what power these approaches have to detect hypothesised effect sizes.

#load the packages needed : install these if not installed
library("MCMCglmm")    #Bayesian mixed models
library(MSBVAR)

############################################3
## SIMULATION of example data set
##########################################33333
#this uses MSBVAR library, which has function
#rmultnorm(1, matrix(c(1,2),2,1), vmat=matrix(c(1,1,0,1),2,2))
#rmultnorm(n, mu, vmat, tol = 1e-10)
#n Number of variates to draw. 
#mu m column matrix of multivariate means 
#vmat m x m covariance matrix  
#tol Tolerance level used for SVD of the covariance. Default is 1e-10 
##############
# below simulates a dataset for Nind individuals which are assumed to have been measured two times in each of two environments
# sim.v1: between-individual variance in environment 1
# sim.v2: between-individual variance in environment 2
# sim.r12: individual correlation in values in E1 and E2 
Nind=200
sim.v1=.2
sim.v2=.3
sim.r12=-0.3
sim.cov12=sim.r12*sqrt(sim.v1*sim.v2)
# residual variance assumed to be 0.8
res.var=0.8
#simulate the individual specific trait values: means assumed to be 0 and 1 and E1 and E2
sim<-rmultnorm(Nind, matrix(c(0,1),2,1), vmat=matrix(c(sim.v1,sim.cov12,sim.cov12,sim.v2),2,2))
#store simulated data in a data.frame
sim.data=list()
sim.data$ind=c(1:Nind,1:Nind,1:Nind,1:Nind) #ID of individuals
sim.data$E=c(rep(0,2*Nind),rep(1,2*Nind)) # the environmental value (0 or 1)
sim.data$Efac=c(rep('a',2*Nind),rep('b',2*Nind)) #environemnt as a factor
#the observed value 'z' is the individual specific value plus a randomly drawn value based on the residual variance
sim.data$z=c(sim[1:Nind,1],sim[1:Nind,1],sim[1:Nind,2],sim[1:Nind,2])+rnorm(4*Nind,sd=sqrt(res.var)) 
sim.data=data.frame(sim.data)


#########1###########################################################################
#####          random regression approach with univariate errors ##################
#####################################################################################

#specify the prior  : non informative for the residual variance ; inverse gamma for the RR estimates
prior.m1.RR <- list(R = list(V = 1, n = 1), G = list(G1 = list(V = diag(2)*0.08, n=2)))
#below specifies a paramter expanded prior as an alternative
prior.m1.RR <- list(R = list(V = 1, n = 1), G = list(G1 = list(V = diag(2)*0.1, nu=diag(2),alpha.mu=c(0,0),alpha.V=diag(2)*25^2)))
#specify the model, first order RR: the covariate is 'E' which is 0 or 1. Use 'poly' to specify higher orders
#this model specification takes some time to run because it used lots of iterations.
#This is an overkill but tends to produce decent uncorrelated chain. Should be checked and possibly adjusted to fit individual needs
mRR.sim <- MCMCglmm(z ~ E, random = ~us(1 + E):ind, data = sim.data, verbose = FALSE, pr = TRUE, prior = prior.m1.RR,  nitt=550000,thin=500,burnin=50000)

#check the output
summary(mRR.sim)
plot(mRR.sim$VCV) #plots the posteriors of the random regression covariances
autocorr(mRR.sim$VCV) #the autocorrelation between the covariances which should be low

########2############################
####### CROSS ENVIRONMENT ###########
#####################################
#derive the environment-specific variances and the covariance between these
phi=matrix(c(1,1,0,1),2,2)  #the second column of this matrix denotes the points where you want to evaluate the RR
#below declares the output we want to get, V(ariance in environment)1, R(epeatability in environment)1, etc...
var.r<-data.frame(V1=numeric(),V2=numeric(),R1=numeric(),R2=numeric(),r12=numeric())
for (i in 1:length(mRR.sim$VCV[,1])) {
  K=matrix(mRR.sim$VCV[i,1:4],2,2) #K matrix for the i'th posterior
  P <- phi%*%K%*%t(phi)  #apply the relevant equations
  temp<-data.frame(V1=P[1,1], V2=P[2,2], R1=P[1,1]/(P[1,1]+as.vector(mRR.sim$VCV[i,5])), R2=P[2,2]/(P[2,2]+as.vector(mRR.sim$VCV[i,5])), r12=P[1,2]/sqrt(P[1,1]*P[2,2]))
  var.r<-rbind(var.r,temp)     #store the data
} #for (i in

#########3################################################################
##### code below implements the character state approach ##################
###########################################################################

#lines below fit the Character State model with a univariate residual error (as in the RR approach above)
priorb = list(R = list(V = diag(1), nu = 0.002), G = list(G1 = list(V = diag(2) * 0.2, nu = 4)))
m.cs <- MCMCglmm(z ~ E, random = ~us(Efac):ind, data = sim.data, verbose = FALSE, prior = priorb,nitt=150000,thin=100,burnin=50000)
summary(m.cs)

#### MAB CHECKING FOR UNDERSTANDING
m.cs$VCV
table(sim.data$Efac, sim.data$ind)

library(rstan)
library(brms)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

m.cs.brm <- brm(z ~ E + (0 + Efac|ind),data = sim.data, prior = c(prior(normal(0,1), class='b'), prior(normal(0,1), class='sd')), control=list(adapt_delta=0.9, max_treedepth=20))
summary(m.cs.brm)
VarCorr(m.cs.brm)

####

######4##########################################
## Organise the output in a data.frame 
################################################
output=list()
output$approach=c("BayesianRR","BayesianCS")
output$V1=c( posterior.mode(as.mcmc(var.r$V1)),posterior.mode(m.cs$VCV[,1]) )
output$V1.low=c(HPDinterval(as.mcmc(var.r$V1))[1],HPDinterval(m.cs$VCV[,1])[1] )
output$V1.high=c(HPDinterval(as.mcmc(var.r$V1))[2],HPDinterval(m.cs$VCV[,1])[2] )
output$V2=c( posterior.mode(as.mcmc(var.r$V2)),posterior.mode(m.cs$VCV[,4]) )
output$V2.low=c(HPDinterval(as.mcmc(var.r$V2))[1],HPDinterval(m.cs$VCV[,4])[1] )
output$V2.high=c(HPDinterval(as.mcmc(var.r$V2))[2],HPDinterval(m.cs$VCV[,4])[2] )
output$Res=c( posterior.mode(mRR.sim$VCV[,5]),posterior.mode(m.cs$VCV[,5]) )
output$Res.low=c(HPDinterval(mRR.sim$VCV[,5])[1],HPDinterval(m.cs$VCV[,5])[1] )
output$Res.high=c(HPDinterval(mRR.sim$VCV[,5])[2],HPDinterval(m.cs$VCV[,5])[2] )
#repeatability in each environment
output$R1=c(posterior.mode( as.mcmc(var.r$V1/(var.r$V1+mRR.sim$VCV[,5]) ) ), posterior.mode(m.cs$VCV[,1]/(m.cs$VCV[,1]+m.cs$VCV[,5]) ) ) 
output$R1.low=c(HPDinterval( as.mcmc(var.r$V1/(var.r$V1+mRR.sim$VCV[,5]) ) )[1], HPDinterval(m.cs$VCV[,1]/(m.cs$VCV[,1]+m.cs$VCV[,5]) )[1] )
output$R1.high=c(HPDinterval( as.mcmc(var.r$V1/(var.r$V1+mRR.sim$VCV[,5]) ) )[2], HPDinterval(m.cs$VCV[,1]/(m.cs$VCV[,1]+m.cs$VCV[,5]) )[2] )
output$R2=c(posterior.mode( as.mcmc(var.r$V2/(var.r$V2+mRR.sim$VCV[,5]) ) ), posterior.mode(m.cs$VCV[,4]/(m.cs$VCV[,4]+m.cs$VCV[,5]) ) ) 
output$R2.low=c(HPDinterval( as.mcmc(var.r$V2/(var.r$V2+mRR.sim$VCV[,5]) ) )[1], HPDinterval(m.cs$VCV[,4]/(m.cs$VCV[,4]+m.cs$VCV[,5]) )[1] )
output$R2.high=c(HPDinterval( as.mcmc(var.r$V2/(var.r$V2+mRR.sim$VCV[,5]) ) )[2], HPDinterval(m.cs$VCV[,4]/(m.cs$VCV[,4]+m.cs$VCV[,5]) )[2] )
#cross-environmental correlation
output$r12=c(posterior.mode(as.mcmc(var.r$r12)),posterior.mode(m.cs$VCV[,2]/sqrt(m.cs$VCV[,1]*m.cs$VCV[,4]) ) )
output$r12.low=c(HPDinterval(as.mcmc(var.r$r12))[1],HPDinterval(m.cs$VCV[,2]/sqrt(m.cs$VCV[,1]*m.cs$VCV[,4]) )[1] )
output$r12.high=c(HPDinterval(as.mcmc(var.r$r12))[2],HPDinterval(m.cs$VCV[,2]/sqrt(m.cs$VCV[,1]*m.cs$VCV[,4]) )[2] )

data.frame(output) #present a summary of the outputs produced by these approaches

#######5#######################################################
## raw data correlation between E's (for comparison)
###############################################################

ind.data<-aggregate(sim.data$z,list(sim.data$ind,sim.data$E),mean) #mean z per E per individual
cor(ind.data[1:Nind,3],ind.data[(Nind+1):(2*Nind),3]) #raw data individual-level correlation between E1 and E2


#####6###################################################
## plotting the output of the RR model
########################################################

### following gives an example of plotting the data and the random regression expectations on top of it
plot(z ~ jitter(E), data = sim.data, cex.lab = 1.5)   #plots the data, a random jitter is added along the X-axis
mus=c(posterior.mode(mRR.sim$Sol[,1]),posterior.mode(mRR.sim$Sol[,1])+posterior.mode(mRR.sim$Sol[,2])) #means in season 0 and 1
#below draws lines in the graph
lines( mus ~ c(0,1), lwd = 2)
#below calculates the expected SD per environment as the sum of the residual and the environment-specific
#between-individual variance
sds=c(sqrt(posterior.mode(as.mcmc(var.r$V1+mRR.sim$VCV[,5]) ) ), sqrt(posterior.mode(as.mcmc(var.r$V2+mRR.sim$VCV[,5]) ) ) )
#draw the expected 95% confidence around the intercepts, assuming the variances are normally distributed
lines(I(mus + 1.96 * sds) ~c(0,1), lty = 2, lwd = 2)
lines(I(mus - 1.96 * sds) ~c(0,1), lty = 2, lwd = 2)
#clearly, a satisfactory fit of the RR should produce an expectation of the means (solid lines) and the variance(dashed
# lines) which reflects the data.
# See Hadfield's MCMCglmm Course notes for another example (type: vignette("CourseNotes", "MCMCglmm") in the R Console)
