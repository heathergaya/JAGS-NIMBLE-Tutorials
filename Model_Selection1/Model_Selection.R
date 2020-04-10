###### Model Selection - Indicator Variable Selection by Heather Gaya #####
###### Last updated April 10, 2020.                          #######
###### Email heather.e.gaya@gmail.com for questions          #######
###### find me on Twitter @doofgradstudent                   ####### 

#set your working directory and get your libraries like a good bub
setwd("~/Desktop/U Georgia/R_Workshop/JAGS:NIMBLE/Model_Selection")
library(runjags)
library(nimble)
library(coda)

######Fake Scenario ######
#There are 40 sites, each one visited 6 times. We recorded present/absent data for Canada Warblers (CAWA) as well as noise levels at the time of detection, and elevation and average rainfall at the site.
#We want to know if bird precense is best predicted by elevation, rainfall, both variables or neither. We will use a model selection techniques to try to estimate which model is the best predictor and get a model averaged estimate of total occupancy.

Birds <- read.csv("Bird_Obs.csv")
Sites <- read.csv("Site_Info.csv")


######### Indicator Variable Selection in JAGS #######
# We'll add two "switches" to our occupancy model that turn the variables on and off between iterations. We'll then monitor the switch parameters to see which model was chosen most often. 
#We'll use the same model we used in the "simple occupancy" tutorial I wrote last week - see my website for more info 

modelstring.occ = "
model {
  ## Loop over sites
  for (i in 1:n.sites) { #40 different sites
  
  # Linear model for true occupancy; we'll call it psi to follow convention
    logit(psi[i]) <- psi.b0 + switch[1]*psi.b1*elevation[i] + switch[2]*psi.b2*rainfall[i] #add in switches
    occ[i] ~ dbern(psi[i]) 
    

#Now we need to deal with the detection process     
    
    ## Loop over replicates within site
    for (t in 1:6) {
      logit(p[i,t]) <- p.b0 + p.b1*noise[i,t]
      
      #The actually observed data is obs[i,t]
      obs[i,t] ~ dbern(p[i,t] * occ[i])
   }
  }
  
  totalocc <- sum(occ[]) #this just says sum across all the sites in occ 
  
  ## Priors
  psi.b0 ~ dnorm(0, 0.37)
  #note that psi.b1 and psi.b2 are not here - we will deal with those in a bit
 
  p.b0 ~ dnorm(0, 0.37) 
  p.b1 ~ dnorm(0, 0.37)
  
  ## Switch Priors
  #Prior for each predictor variable
  #should be uninformative when on but very informative when off 
  # because if it's off, then we want the value to be 0
  psi.b1 ~ dnorm(0, prec[1])
  psi.b2 ~ dnorm(0, prec[2])
  
  for (k in 1:2){ #one for each variable that controls psi 
    # Create 'slab & spike' precision
    prec[k] <- (1-switch[k])*(1000) + switch[k]*(.37)
    #if the switch is on, then switch [k] =1 and the precision = .37 
    #if the siwtch is off, switch[k] = 0 and precision = 1000 (so psi.b# will equal 0)

    # Sampling and prior for indicator variable (switch)
    switch[k] ~ dbern(p.w[k])  #is the switch on or off this iteration? 
    p.w[k] ~ dunif(0,1) #probability the switch is on
  }

  # Converts the indicator variable pair into a model number
  for (b in 1:4){
    m[b] <- (b == 1 + switch[1] + 2*switch[2])
     #In NIMBLE (and JAGS, and R), TRUE = 1 and FALSE = 0. Very convenient. 
  }
  #mod 1 is when both switches are off (when 1 + switch[1] + 2*switch[2] = 1). Mod 2 is when switch 1 (elevation) is on (2 = 1 + 1 +2*0). Mod 3 = rainfall on, elevation off, and mod 4 is both covariates are used.
}
"

#so what's actually happening in this model? 
#basically, at each iteration, JAGS is "trying out" possible values for all the parameters - sometimes that means elevation is used to model occupancy, sometimes rainfall, sometimes both (and sometimes neither!). Each time it runs an iteration, it looks at the likelihood that the data came from a distribution that relies on those predictors. *IF* one model really is a better model, we would expect the data to have a much higher likelihood under that model - a much better "fit".
# As a better model is "found", that combo of variables will be accepted more often than other combinations of our predictor variables. Then, after the model runs, we can ask "how many times was this model "accepted"?" 
#and in theory, the best model should be selected the most often. Let's see what happens. 

# Bundle data
jd.Foo <- list(n.sites = nrow(Birds), elevation = Sites$Elevation, rainfall = Sites$Rainfall,
               noise = as.matrix(Birds[,8:13]), obs = as.matrix(Birds[,2:7]))

ji.Foo <- function(){list(psi.b0 = rnorm(1, 0, 0.1), psi.b1 = rnorm(1, 0, 0.1), 
                          psi.b2 = rnorm(1, 0, 0.1), switch = rep(0,2), #start with both off
                          p.w = rnorm(2, 0.5, 0.1), p.b0 = runif(1,0,1), occ = rep(1,nrow(Birds)))}

# Parameters to estimate
params <- c("psi.b0", "psi.b1", "psi.b2", "p.b0", "p.b1", "switch", "m", "p.w", "totalocc")

#Run the model
Foo <- run.jags(model = modelstring.occ, monitor = params, data = jd.Foo, n.chains = 3, inits = ji.Foo, burnin = 5000, sample = 12000, adapt = 1000, method = "parallel")
summary(Foo)
plot(Foo)
#You'll notice some of the plots look weird - almost like they didn't converge. This is why this method is sometimes called the "slab and spike" method. When the beta values (psi.b1, psi.b2) are turned off, they are forced to 0, causing strange jumps and horizontal lines in the MCMC chains. This is to be expected from this method - though when in doubt more iterations don't hurt

#We can look at the mean values of the "m" parameters and the posterior distribution of the p.w values to tell us our model "weights" 
#In this case, mod 1 (no covariates for psi) was chosen about 9% of the time, mod2 (elevation) about 76%, mod 3 (rainfall), 1.2% and mod 4 (elevation + rainfall) about 14%. From this, we can see that model 2 is likely the best model! 
#And fun fact, it is also the model I simulated this data from so yay! 

#p.w tells us how often that predictor was used in any model. About 63% of the time elevation was used (the switch was "on") and about 38% of the time rainfall was "on"

#The COOLEST part about all this (I think) is that it's doing model selection AND model averaging all in the same step!
#Pretty nifty, eh? 


######## Indicator Variable Selection in NIMBLE #######
## A few changes
#first telling NIMBLE it's a NIMBLE code
# most importantly, NIMBLE thinks the word SWITCH means something else, so we'll change it to "sw"

nimble.birds <- nimbleCode({ 
  
  ## Loop over sites
  for (i in 1:n.sites) { #40 different sites
  
  # Linear model for true occupancy; we'll call it psi to follow convention
    logit(psi[i]) <- psi.b0 + sw[1]*psi.b1*elevation[i] + sw[2]*psi.b2*rainfall[i] #add in switches
    occ[i] ~ dbern(psi[i]) 
    

#Now we need to deal with the detection process     
    
    ## Loop over replicates within site
    for (t in 1:6) {
      logit(p[i,t]) <- p.b0 + p.b1*noise[i,t]
      
      #The actually observed data is obs[i,t]
      obs[i,t] ~ dbern(p[i,t] * occ[i])
   }
  }
  
  totalocc <- sum(occ[1:n.sites]) #have to change indexing for NIMBLE
  
  ## Priors
  psi.b0 ~ dnorm(0, 0.37)
  #note that psi.b1 and psi.b2 are not here - we will deal with those in a bit
 
  p.b0 ~ dnorm(0, 0.37) 
  p.b1 ~ dnorm(0, 0.37)
  
  ## Switch Priors
  #Prior for each predictor variable
  #should be uninformative when on but very informative when off 
  # because if it's off, then we want the value to be 0
  psi.b1 ~ dnorm(0, prec[1])
  psi.b2 ~ dnorm(0, prec[2])
  
  for (k in 1:2){ #one for each variable that controls psi 
    # Create 'slab & spike' precision
    prec[k] <- (1-sw[k])*(1000) + sw[k]*(.37)
    #if the switch is on, then switch [k] =1 and the precision = .37 
    #if the siwtch is off, switch[k] = 0 and precision = 1000 (so psi.b# will equal 0)

    # Sampling and prior for indicator variable (switch)
    sw[k] ~ dbern(p.w[k])  #is the switch on or off this iteration? 
    p.w[k] ~ dunif(0,1) #probability the switch is on
  }

  # Converts the indicator variable pair into a model number
  for (b in 1:4){
    m[b] <- ((1 + sw[1] + 2*sw[2]) == b)
    
    #In NIMBLE (and JAGS, and R), TRUE = 1 and FALSE = 0. Very convenient. 
      }
  
  #mod 1 is when both switches are off (when 1 + sw[1] + 2*sw[2] = 1). Mod 2 is when switch 1 (elevation) is on (2 = 1 + 1 +2*0). Mod 3 = rainfall on, elevation off, and mod 4 is both covariates are used.

})


# Next we need to give NIMBLE data and constants. Our only "data" is the observations; everything else is a constant

data <- list(obs = as.matrix(Birds[,2:7]))
constants <- list(n.sites = nrow(Birds), elevation = Sites$Elevation, rainfall = Sites$Rainfall,
               noise = as.matrix(Birds[,8:13]))
inits <- list(psi.b0 = rnorm(1, 0, 0.1), psi.b1 = rnorm(1, 0, 0.1),p.b1 = rnorm(1, 0, 0.1),
              psi.b2 = rnorm(1, 0, 0.1), sw = rep(0,2), #start with both off
              p.w = rnorm(2, 0.5, 0.1), p.b0 = runif(1,0,1), occ = rep(1,nrow(Birds)))
#make sure we rename "switch" to "sw"
params <- c("psi.b0", "psi.b1", "psi.b2", "p.b0", "p.b1", "sw", "m", "p.w", "totalocc")

#Fun fact - before writing this code, I had never encountered the "switch" problem before. It threw a very exciting error - "unable to create shared library" which was 0% helpful when I googled it. Turns out that error can mean many things! Ah, debugging. A daily struggle. 

#run it the speedy way
b1 <- nimbleMCMC(code = nimble.birds, constants = constants, 
                 data = data, inits = inits, monitors = params, niter = 30000, thin = 1, nchains =3, nburnin = 20000, samplesAsCodaMCMC = TRUE)

summary(b1)
gelman.diag(b1, multivariate = F)
plot(b1)
#The graphs are... lovely. 


###### A table of results, for comparison ######

jags <- summary(Foo)
nimb <- summary(b1)$statistics
data.frame(Mod = 1:4, 
           Weights.jags = round(jags[8:11,4], digits =2), 
           Weights.nimble = round(nimb[1:4,1], digits =2))


