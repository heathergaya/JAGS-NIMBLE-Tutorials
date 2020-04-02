### ################ Occupancy Modeling by Heather Gaya.  ########
###### Last updated April 2, 2020.                          #######
###### Email heather.e.gaya@gmail.com for questions          #######
###### find me on Twitter @doofgradstudent                   ##### 

## Occupancy modeling is awesome! And much more applicable to what I actually do in my real life than linear regression... though of course it includes linear regression! Wheee!
# Anyway. Occupancy modeling is more than just "does this thing exist in this space?" It's also about estimating "why does this thing exist in this space but NOT that space?" 

#We will start with single species occupancy modeling under a very simple scenario. 

#### A Fake Scenario ###### 
## Let's say we have some birds. They're cool, we like them. Let's say they're Canada Warblers (or CAWA), because those are the things I actually get paid to study, despite my love of all things herpetology. We go out and do point counts - this involves standing in one spot and listening/watching for birds and recording what you see/hear. In addition to counting birds, we also record things about the environment on the day we go there - namely environmental noise. We do repeated observations at each point. We have previously recorded the elevation of the sites where we do point counts. 

# For this example, let's say we are only interested in presence/absence. The actual number of birds does not matter. We also assume the presence/absence of birds does not change with time. Obviously these are big assumptions, but let's start simple. 

##Let's set our working directory and grab our data
setwd("~/Desktop/U Georgia/R_Workshop/JAGS:NIMBLE/Occupancy_Models")
Birds <- read.csv("Bird_Obs.csv")
Sites <- read.csv("Site_Info.csv")

#There are 40 sites, each one visited 6 times. Sites have elevation that's the same for all time periods (we hope!). Each bird observation period has covariates of noise associated with them. 

#Let's see how many sites we even saw birds at:
sum(rowSums(Birds[,2:7]) > 0) #looks like 17 sites for sure had birds; This is our raw occupancy 
#Let's also do a quick check of our variables, just to give ourselves an idea of what we're working with. Let's plot the site elevation vs. our raw occupancy estimates

plot(Sites$Elevation, rowSums(Birds[,2:7]) > 0) #this turns the row sums into T/F, which plots as 1s and 0's! 
#So maybe a positive trend with elevation 
#Remember that just becuase we didn't see birds, doesn't mean birds aren't at that site. This is why we go to all the trouble of modeling this stuff at all. 

##### Conceptualizing the Model #######
# Let's think about the data we have - 1's and 0's. What process creates a 1? That's pretty easy - that means a bird was there and we detected it. But what causes a 0? 
# 0's can be from one of two processes - either the bird was there (the site was suitable for birds) and we missed it, or the bird wasn't there and so we saw nothing. So that means we have two separate processes here that are mixing together - yay mixture models! 
#Let's think about the process that causes birds to be present or absent. For CAWA, we know they seek out higher elevation - possibly b/c of lower temps or competition - but sometimes the birds aren't there even when the site seems (to us) as pretty suitable. So this gives us the idea that there's some probability involved, but generally speaking CAWA are more likely to be found at higher elevation. 
#There are multiple ways to model this, but I like to think of it as 
#prob(suit) <- intercept + beta1*elevation 
#and then:
# actually suitable ~ bernouli(prob(suit))
# What about the process of seeing birds? 
# If the bird is not there, the probability of seeing the bird is 0 (we hope). If the bird is there... well, we know that noise in the environment probably lowers detection probability, so we can once again model this is a linear regression:
# prob(detect) <- intercept + beta1*noise
# detected ~ bernouli(prob(detect))

#Of course, detecting a bird also requires that the species is present, so p(detect) is conditional on the animal existing at the site

### Let's set this up in JAGS #####
library(runjags)
modelstring.occ = "
model {
  ## Loop over sites
  for (i in 1:n.sites) { #40 different sites
  
  # Linear model for true occupancy; we'll call it psi to follow convention
    logit(psi[i]) <- psi.b0 + psi.b1*elevation[i]
    
    #note that this needs to be a probability, so we'll put it on the logit scale to force it to stay between 0 and 1 
    
    # Occupancy is then drawn with probability psi; This comes out as a 1 or a 0. 
    occ[i] ~ dbern(psi[i]) 
    

#Now we need to deal with the detection process     
    
    ## Loop over replicates within site
    for (t in 1:6) {
      # each i is a site, each t is a time period 
      # Link function for detection model
      logit(p[i,t]) <- p.b0 + p.b1*noise[i,t]
      #Notice that have to use nested indexing to account for the 5 possible noise levels
      #again, the probability needs to be between 0 and 1 so we put it on the logit scale
      
      #The actually observed data is obs[i,t]
      obs[i,t] ~ dbern(p[i,t] * occ[i])
      #if occ is 0, we cannot see birds, no matter how high our detection prob is. 
      #if occ is 1, we possibly can see birds, but depends on detection. 
      # The result (obs[i,t]) is our data. 
    }
  }
  
  
  ## Priors
  # Because our model is on the LOGIT scale, these will turn into approximately  uniform priors 
  psi.b1 ~ dnorm(0, 0.37)
  psi.b0 ~ dnorm(0, 0.37)
  p.b0 ~ dnorm(0, 0.37) 
  p.b1 ~ dnorm(0, 0.37)
  
  # How many sites are truly occupied?		
  totalocc <- sum(occ[]) #this just says sum across all the sites in occ 
  
  #Let's also graph the relationship between occupancy and elevation. 
  #We'll give R a sequence of values from .2km to 1.5km (reasonable, given our data) and ask for the probability of occupancy. n.graph and the fake.elev values will come in as data.
    for (k in 1:n.graph){
    logit(graph.me[k]) <- psi.b0 + psi.b1*fake.elev[k]
    }
  
}
"

#nothing special about the data
jd.Foo <- list(n.sites = nrow(Birds),elevation = Sites$Elevation, noise = as.matrix(Birds[,8:13]), obs = as.matrix(Birds[,2:7]), fake.elev =seq(.2, 1.5, by = .025), n.graph = 53)

#Sometimes JAGS has trouble if it initializes the occupancy of site as 0 but then also has conflicting information with the betas about the occupancy of the site. To avoid issues starting the model, it's easier to just initialize all sites as occupied and let the MCMC chains drop the number down to whatever the true value may be.                                               
ji.Foo <- function(){list(psi.b0 = runif(1,0,1), psi.b1 = runif(1,0,1), p.b0 = runif(1,0,1), p.b1 = runif(1,0,1), occ = rep(1,nrow(Sites)))}

params <- c("psi.b0", "psi.b1", "p.b0", "p.b1", "totalocc", "graph.me")

Foo <- run.jags(model = modelstring.occ, monitor = params, data = jd.Foo, n.chains = 3, inits = ji.Foo, burnin = 4000, sample = 10000, adapt = 1000, method = "parallel")
mod <- summary(Foo)
mod
plot(Foo) # I reccomend running the code with the "graph.me" commented out in the params line above, plotting and checking for convergence then rerunning it with the graph.me added back in. Saves a lot of time when looking through the plots.

#Looks like we found a positive trend in elevation for these birds (as expected) and a negative relationship between noise and detection probability (again, expected). Our model predicts 18-33 of the 40 sites are actually occupied by CAWA. 

#real occupancy in the simulated data = 25 sites so we didn't do too badly! 

##### Time for plotting! ######
birdplot <- data.frame(Elevation =seq(.2, 1.5, by = .025), Mean = mod[6:58,4], Low = mod[6:58,1], High = mod[6:58,3])

#Base R soluation 
plot(birdplot$Elevation, birdplot$Mean, type = "l", main = "Site Suitability", xlab = "Elevation (km)", ylab = "Probability Suitable", ylim = c(0,1))
lines(birdplot$Elevation, birdplot$Low, lty = 2)
lines(birdplot$Elevation, birdplot$High, lty = 2)

#ggplot solution
library(ggplot2)
G2 <- ggplot(birdplot, aes(x = Elevation))+
  geom_smooth(aes(y = Mean, ymin = Low, ymax = High), fill = "lightblue", colour = "black", stat = "Identity")+
  labs(x = "Elevation (km)", y = "Probability Suitable")+
  ylim(0,1)+
  theme_classic(base_size = 18)+
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background =element_blank())
G2

########### The NIMBLE Version ########
library(nimble)
library(coda)
nimble.birds <- nimbleCode({ #this is change #1
  ## Loop over sites
  for (i in 1:n.sites) { 
    logit(psi[i]) <- psi.b0 + psi.b1*elevation[i]
    occ[i] ~ dbern(psi[i]) 
    
    ## Loop over replicates within site
    for (t in 1:6) {
      logit(p[i,t]) <- p.b0 + p.b1*noise[i,t]
    
      #The actually observed data is obs[i,t]
      obs[i,t] ~ dbern(p[i,t] * occ[i])
    }
  }

  ## Priors
  psi.b1 ~ dnorm(0, 0.37)
  psi.b0 ~ dnorm(0, 0.37)
  p.b0 ~ dnorm(0, 0.37) 
  p.b1 ~ dnorm(0, 0.37)
  
  # How many sites are truly occupied?		
  totalocc <- sum(occ[1:n.sites]) #this is difference #2. NIMBLE makes you specify indexing. 
  
  #Let's also graph the relationship between occupancy and elevation. 
    for (k in 1:n.graph){
    logit(graph.me[k]) <- psi.b0 + psi.b1*fake.elev[k]
    }
  
})

## Just two changes and we've got NIMBLE code. Whooo. Onward to running it. 
params <- c("psi.b0", "psi.b1", "p.b0", "p.b1", "totalocc", "graph.me")

#Next we need to give NIMBLE data and constants. Our only "data" is the observations; everything else is a constant

data <- list(obs = as.matrix(Birds[,2:7]))
constants <- list(n.sites = nrow(Birds),elevation = Sites$Elevation, noise = as.matrix(Birds[,8:13]), fake.elev =seq(.2, 1.5, by = .025), n.graph = 53)

#Initial values
inits <- list(psi.b0 = runif(1,0,1), psi.b1 = runif(1,0,1), p.b0 = runif(1,0,1), p.b1 = runif(1,0,1), occ = rep(1,nrow(Sites)))

#Run the model - fast way
b1 <- nimbleMCMC(code = nimble.birds, constants = constants, 
data = data, inits = inits, monitors = params, niter = 30000, thin = 1, nchains =3, nburnin = 10000, samplesAsCodaMCMC = TRUE)

# The long way (but useful to learn!!!)
prepbirds <- nimbleModel(code = nimble.birds, constants = constants, 
                         data = data, inits = inits) 
prepbirds$initializeInfo() #everything is good to go! 
mcmcbirds <- configureMCMC(prepbirds, monitors = params, print = T )
# this tells us that psi and p are going to be using RW samplers (continous) and occupancy will be determined with binary samplers (1s or 0s)
birdsMCMC <- buildMCMC(mcmcbirds, enableWAIC = TRUE) # a change here to allow us to use WAIC. We don't have anything to compare here, but if we did, WAIC could be useful. 
Cmodel <- compileNimble(prepbirds) #compiling the model itself in C++; it says "may take a minute" but on big models this can take HOURS. 
Compbirds <- compileNimble(birdsMCMC, project = prepbirds) # compile the samplers next 
bird.mod.nimble <- runMCMC(Compbirds, niter = 30000, thin = 1, nchains =3, nburnin = 10000, samplesAsCodaMCMC = TRUE)

#Look at convergence for the non-graphing nodes:
gelman.diag(mcmc.list(bird.mod.nimble$chain1[,54:58],bird.mod.nimble$chain2[,54:58],bird.mod.nimble$chain3[,54:58]) ) 

#summary stats:
summary(mcmc.list(bird.mod.nimble$chain1[,54:58],bird.mod.nimble$chain2[,54:58],bird.mod.nimble$chain3[,54:58]))

#plot of chains:
plot(mcmc.list(bird.mod.nimble$chain1[,54:58],bird.mod.nimble$chain2[,54:58],bird.mod.nimble$chain3[,54:58]))

####Graphing is same idea as JAGS
mod.b <- summary(mcmc.list(bird.mod.nimble$chain1[,1:53],bird.mod.nimble$chain2[,1:53],bird.mod.nimble$chain3[,1:53]))
birdplot <- data.frame(Elevation =seq(.2, 1.5, by = .025), Mean = mod.b$statistics[,1], Low = mod.b$quantiles[,1], High = mod.b$quantiles[,5])

#Base R soluation 
par(mfrow = c(1,1))
plot(birdplot$Elevation, birdplot$Mean, type = "l", main = "Site Suitability", xlab = "Elevation (km)", ylab = "Probability Suitable", ylim = c(0,1))
lines(birdplot$Elevation, birdplot$Low, lty = 2)
lines(birdplot$Elevation, birdplot$High, lty = 2)

#ggplot solution
library(ggplot2)
G2 <- ggplot(birdplot, aes(x = Elevation))+
  geom_smooth(aes(y = Mean, ymin = Low, ymax = High), fill = "lightblue", colour = "black", stat = "Identity")+
  labs(x = "Elevation (km)", y = "Probability Suitable")+
  ylim(0,1)+
  theme_classic(base_size = 18)+
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background =element_blank())
G2

