################ Linear Regression Part 2 and 3 by Heather Gaya.  ########
###### Last updated March 28, 2020.                          #######
###### Email heather.e.gaya@gmail.com for questions          #######
###### find me on Twitter @doofgradstudent                   ##### 

#Topics covered: Categorical variables in linear regression, DIC, and graphing output with credible intervals!

# This code assumes you already have JAGS and/or NIMBLE installed on your computer (they are standalone programs)

# Quick reminders about R - 
# the "#" symbol turns things into comments (instead of code)
# set your working directory! Session -> set working directory --> (to where your files are)

### A fake scenario, same as in Part 1 #####
# Joe, our favorite fake researcher, goes out into the world to study some fake frogs
# He collects 50 frogs and records some measurements on them 
# specifically, we records:
# age (in days)
# weight (g)
# back leg length (cm)
# distance from road (m)
# and species (A or B)

#Previously, Joe was only interested in how weight was related to age, leg length and distance from the road, but now he wants to model the two species separately. As before, we are kind souls and love stats, so of course we will help him. 

#First, set your working directory like always
setwd("~/Desktop/U Georgia/R_Workshop/JAGS:NIMBLE/Linear_Regression_2") #you will need to put your own file path in here

# fetch the data again
Frogs <- read.csv("Fakefrogs.csv")
head(Frogs)

# R packages we will be using to analyze 
# You will need to install these first if you haven't 
# install.packages("name-of-your-package")

library(runjags) #my preferred JAGS package
library(nimble) # only need this if you want to use NIMBLE 
library(parallel) # for faster computing
library(coda) #for inspecting chains

#Previously, we were only working with continous variables. Our equation was roughly:
#Expected weight of frog <- intercept + age*something + leg*something2 + distance*something3
#  and the actual weight of frog is from normal distribution with mean = expected weight and standard deviation = some sd

#But now we want to model the frog species sepeately.
# There are many options here for the actual equation, depending on what we think the relationship is with species - do we think the intercept is different? Maybe one frog grows faster than the other, so we would expect a difference in the relationship with age? We will test out both of these models and attempt to find the best model using DIC. 

#Let's start with just the intercepts being different. Maybe the relationship with the variables are the same, but one frog just weighs more in general. Conceptually, we might think of our expected weight equation as:
#E(weight) <- interceptA*(1 if species A, 0 if species B) + interceptB*(1 if species B, 0 if species A) + age*something + leg*something2 + distance*something3
# and actual weight = Normal(mean = expected weight and standard deviation = some sd) 

#But we can't use those "1 if A, 0 if B" type statements in JAGS or NIMBLE. Instead, we use something called "nested indexing". 
#To make these easy on ourselves, we can create a new "variable", which we'll call "s" for species that turns the species into a number instead of a letter.
#So if s =1, that means species A, if s = 2, that means species B. There's a quick way to grab this from our data exploiting the fact that species is already a FACTOR in R. Factors have levels, and levels have a number associated with them (the 1st level = 1, the 2nd level =2, etc.)
levels(Frogs$Species)
s <- as.numeric(Frogs$Species) #and just like that, the A's are all 1's and the B's are all 2's! 
s
#Okay, so now we have a vector of 1's and 2's... yay? Let's see how this helps us. 

#Instead of just one value, we can turn the intercept (beta0) into a vector of two values, the first representing A and the second representing species B. We then use our vector s to indicate if we want the first or second value. Here's a quick example in R:
fake.intercepts <- c(2, 5)
fake.intercepts[s[1]] #s[1] = the first value in our species vector, which happens to equal 1. So this selects fake.intercepts[1] 
fake.intercepts[s[30]] #s[30] = 2, so this statement equals fake.intercepts[2] 
fake.intercepts[s] #this is what the resulting intercepts would be for all the indiviudals in our data, if these intercepts weren't just fake numbers

#Now that we know that, let's put that into JAGS: 

### Writing JAGS model # 1 ###### 
modelstring.Frogs1 = "
  model 
{
for (i in 1:n.frogs){ 
  meanweight[i] <- beta0[s[i]] + age[i]*beta1 + leg[i]*beta2 + distance[i]*beta3
  #note that the only thing that changes is the indexing of beta0
  weight[i] ~ dnorm(meanweight[i], prec)
}

#but now beta0 can be two different values, we need to have two priors
beta0[1] ~ dunif(-50,50) 
beta0[2] ~ dunif(-50,50) 
beta1 ~ dunif(-50,50)
beta2 ~ dunif(-50,50)
beta3 ~ dunif(-50,50)


#I like to convert precision to standard deviation because it confuses me otherwise. 
prec <- 1/(sd *sd)
sd ~ dunif(0.0001, 100)
}
" 

#That's all there is to continous variables. If species had more levels, we would need to give a prior for beta0 for each level. 

#Let's run the model now
# As before we need parameters to model. Note that because beta0 is INDEXED, we can just say we want "beta0" and jags will give us back beta0[1] and beta0[2]
params <- c("beta0", "beta1", "beta2", "beta3","sd") 

#Next we need to give JAGS data. This is the same as in part 1, with the addition of the new "s" variable. 
data <- list(weight = Frogs$Weight, distance = Frogs$Dist, 
             age = Frogs$Age, leg = Frogs$Leg,
             n.frogs = 50, s = s)

#Initial values are similar to before, but we need to give beta0 2 values instead of 1. 
#reminder that runif arguments are runif(how many, minimum, maximum)
inits <- function(){list(beta0 = runif(2, -10, 10), beta1 = runif(1,-10,10), beta2 = runif(1,-10,10), beta3 = runif(1,-10, 10))}

#Run the model
Frog.mod1 <- run.jags(model = modelstring.Frogs1, 
                     monitor = params, inits = inits,
                     data = data, n.chains = 3, 
                     sample = 10000, method = "parallel")

#Look at the summary
Frog.mod1 

#looks like the model converged, but we should visually inspect the chains as well 

plot(Frog.mod1)

#okay, not bad! Things seem to have converged. 
# But what do the results mean? Well, when we look at the output for the two beta0 values, they seem to be fairly similar. The credible intervals almost completely overlap.  So how do we know if this model is any better than the other versions? 
#We can use DIC to help us estimate if the addition of more parameters is "worth it" 
# DIC does have restrictions on when it can be used - it is often inappropriate to use it in mixture models - but for simple cases like linear regression (or later, CMR type models), DIC can be used fairly safely.
#Like AIC, we want DIC to be as low as possible but the number itself isn't particularly meaningful 

Dic_1 <- extract(Frog.mod1, what = "dic")
# Mean deviance is the value of interest
# But since we haven't run any other models, we can't compare it to anything! So let's try a few more models below:

### JAGS Models # 2 and # 3 #######

#Maybe one frog grows faster than the other, so we would expect a difference in the relationship with age, but no difference in the intercept. We will inspect this in model 2:

modelstring.Frogs2 = "
  model 
{
for (i in 1:n.frogs){ 
  meanweight[i] <- beta0 + age[i]*beta1 + leg[i]*beta2 + distance[i]*beta3[s[i]]
  #note that the only thing that changes is the indexing of beta0 and beta3
  weight[i] ~ dnorm(meanweight[i], prec)
}

#but now beta0 is only one value and  beta3 can be two different values
beta0 ~ dunif(-50,50) 
beta1 ~ dunif(-50,50)
beta2 ~ dunif(-50,50)
beta3[1] ~ dunif(-50,50)
beta3[2] ~ dunif(-50,50)


#I like to convert precision to standard deviation because it confuses me otherwise. 
prec <- 1/(sd *sd)
sd ~ dunif(0.0001, 100)
}
" 

#we haven't added any new things to monitor (remember, we were already monitoring beta3), so we can just run everything without having to re-write all the data/params/etc. lines. We do need to rewrite the inits line though, to make sure the dimensions of beta0 and beta3 match our new model. 

inits2 <- function(){list(beta0 = runif(1, -10, 10), beta1 = runif(1,-10,10), beta2 = runif(1,-10,10), beta3 = runif(2,-10, 10))}

Frog.mod2 <- run.jags(model = modelstring.Frogs2, 
                      monitor = params, inits = inits2,
                      data = data, n.chains = 3, 
                      sample = 15000, method = "parallel")

Frog.mod2
plot(Frog.mod2)
Dic_2 <- extract(Frog.mod2, what = "dic")
Dic_2

#We can see this might be a slightly better model, but not *overwhelmingly* so. Let's also compare to third model for funsies. Maybe Joe thinks about his frogs and again and thinks that maybe weight is only related to age, but he still suspects there's a different equation (both intercept and slope) for both species. Let's see what that model looks like. 
modelstring.Frogs3 = "
  model 
{
for (i in 1:n.frogs){ 
  meanweight[i] <- beta0[s[i]] + age[i]*beta1[s[i]]
  #got to reindex everything 
  weight[i] ~ dnorm(meanweight[i], prec)
}

#adjust priors to match
beta0[1] ~ dunif(-50,50) 
beta0[2] ~ dunif(-50,50) 
beta1[1] ~ dunif(-50,50)
beta1[2] ~ dunif(-50,50)

prec <- 1/(sd *sd)
sd ~ dunif(0.0001, 100)
}
" 

#Of course, now we actually have to adjust all the data and such, since the params have changed. 
params3 <- c("beta0", "beta1","sd") 
data3 <- list(weight = Frogs$Weight, 
             age = Frogs$Age, n.frogs = 50, s = s)

#Initial values are similar to before, but we need to give beta0 2 values instead of 1. 
#reminder that runif arguments are runif(how many, minimum, maximum)
inits3 <- function(){list(beta0 = runif(2, -10, 10), beta1 = runif(2,-10,10))}

#Run the model
Frog.mod3 <- run.jags(model = modelstring.Frogs3, 
                      monitor = params3, inits = inits3,
                      data = data3, n.chains = 3, 
                      sample = 10000, method = "parallel")
Frog.mod3
#appears to have converged, double check with plot
plot(Frog.mod3)
Dic_3 <- extract(Frog.mod3, what = "dic")
Dic_3

#Let's look at the 3 DICs together in a table
#DICs are stored as deviances for each data point (in this case 50) separated, but we need to sum them together to get the number we want. 

data.frame(Mods = 1:3, DICS = c(sum(Dic_1$deviance), sum(Dic_2$deviance), sum(Dic_3$deviance)))

#Based on DIC, we can see that model 2 remains the best model *OF THE ONES WE TESTED* 
# This doesn't mean model 2 is the *true* model, it just means it is likely the best model of the ones available. This is an important distinction and one that can trip you up if you're not careful. 

#But, okay, awesome. We have our results! But now we want to graph these results - maybe a nice line graph showing the relationship between expected weight and distance to road. 
#This isn't impossible to calculate ourselves, but it's a lot easier to just let JAGS do it. 

###### Graphing the CI #######
#Inside our JAGS model we have the model make point estimates for various weights at various distances for the road - for both species. This will give us the CIs at those point estimates and we'll let ggplot do the rest. 

modelstring.Frogs2graph = "
  model 
{
for (i in 1:n.frogs){ 
  meanweight[i] <- beta0 + age[i]*beta1 + leg[i]*beta2 + distance[i]*beta3[s[i]]
  weight[i] ~ dnorm(meanweight[i], prec)
}

beta0 ~ dunif(-50,50) 
beta1 ~ dunif(-50,50)
beta2 ~ dunif(-50,50)
beta3[1] ~ dunif(-50,50)
beta3[2] ~ dunif(-50,50)


#I like to convert precision to standard deviation because it confuses me otherwise. 
prec <- 1/(sd *sd)
sd ~ dunif(0.0001, 100)

#somewhere in the model, wherever we want, we add in our little graphing loop. Let's graph distances from 0 to 200 m from the road (reasonable, given our data ranges from 1.7 to 198)
#we'll hold age and leg constant, at the mean values we have in our data

for (k in 1:201 ){ #can't iterate starting at 0, so we'll just use k-1 in our equation and start at 1 instead 

  graph_w1[k] <- beta0 + 392*beta1 + 4*beta2 + (k-1)*beta3[1] #species 1; only thing changing is the distance, all other variables held constant 
  graph_w2[k] <- beta0 + 392*beta1 + 4*beta2 + (k-1)*beta3[2] #species 2
}

}
" 
#now we'll want to monitor these new variables graph_w1 and graph_w2

params2 <- c("beta0", "beta1", "beta2", "beta3","sd", "graph_w1", "graph_w2") 

Frog.mod2graph <- run.jags(model = modelstring.Frogs2graph, 
                      monitor = params2, inits = inits2,
                      data = data, n.chains = 3, 
                      sample = 15000, method = "parallel")

#now we've produced 402 more nodes we so don't really want to view all of them.
# instead we'll just save them into an object

res <- summary(Frog.mod2graph)
head(res, n = 10)

#now put into a dataframe for ggplot or base R 

graphme <- data.frame(low = res[7:408,1], upper  = res[7:408,3], 
                      mean = res[7:408,4], Species = rep(c("A", "B"), each = 201), 
                      dist = c(0:200, 0:200))

#in base R:
plot(0:200, graphme$mean[1:201], type = "l", ylim = c(0, 70), xlab = "Distance to Road (m)", ylab = "Predicted Frog Weight (g)", main = "Predicted Frog Weight by Species") #mean
polygon(c(0:200, 200:0, 0), c(graphme$low[1:201], graphme$upper[201:1], graphme$low[1]), lty= "dashed", col = "lightblue", border = "blue")
polygon(c(0:200, 200:0, 0), c(graphme$low[202:402], graphme$upper[402:202], graphme$low[202]), lty= "dashed", col = "azure3", border = "grey20")
lines(0:200, graphme$mean[202:402], col = "grey20") #mean sp 2
lines(0:200, graphme$mean[1:201], col = "blue") #mean sp 1
legend("topleft", c("Species A", "Species B"), col = c("blue", "grey20"), lty = 1)

#obviously you could make this plot prettier. But that's the general idea. 
#ggplot is much easier (if you know ggplot, that is). 
library(ggplot2)
ggplot(data = graphme, aes(x = dist, group = Species, col = Species, fill = Species))+
  geom_smooth(aes(y = mean, ymin = low, ymax = upper), stat = "Identity")+
  labs(x = "Distance to Road (m)", y = "Predicted Frog Weight (g)")+
  ylim(0,70)+
  theme_minimal(base_size = 18)+
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background =element_blank())

### Yay, graphing! Though sadly not a particularly interesting graph all things considered. 

## That's it for Linear Regression Part 2! Part 3 (Random error) coming soon.


############ Model #1 in NIMBLE #########
#Reminder that at this point, we might think of our expected weight equation as:
#E(weight) <- interceptA*(1 if species A, 0 if species B) + interceptB*(1 if species B, 0 if species A) + age*something + leg*something2 + distance*something3
# and actual weight = Normal(mean = expected weight and standard deviation = some sd) 

nimbleFrogs1 <-
  nimbleCode({ #don't forget this part 
  
for (i in 1:n.frogs){ 
  meanweight[i] <- beta0[s[i]] + age[i]*beta1 + leg[i]*beta2 + distance[i]*beta3
  weight[i] ~ dnorm(meanweight[i], sd = sd )
}

beta0[1] ~ dunif(-50,50) 
beta0[2] ~ dunif(-50,50) 
beta1 ~ dunif(-50,50)
beta2 ~ dunif(-50,50)
beta3 ~ dunif(-50,50)

sd ~ dunif(0.0001, 100)
})


#Let's run the model now
# As before we need parameters to model. Note that because beta0 is INDEXED, we can just say we want "beta0" and Nimble will give us back beta0[1] and beta0[2]
params <- c("beta0", "beta1", "beta2", "beta3","sd") 

#Next we need to give NIMBLE data and constants. This is the same as in part 1, with the addition of the new "s" variable. 
data <- list(weight = Frogs$Weight, distance = Frogs$Dist, 
             age = Frogs$Age, leg = Frogs$Leg)
             
constants <- list(n.frogs = 50, s = s)

#Initial values are similar to before, but we need to give beta0 2 values instead of 1. 
#reminder that runif arguments are runif(how many, minimum, maximum)
inits <- list(beta0 = runif(2, -10, 10), beta1 = runif(1,-10,10), beta2 = runif(1,-10,10), beta3 = runif(1,-10, 10), sd = runif(1, .001, 100))

#Run the model
prepfrogs <- nimbleModel(code = nimbleFrogs1, constants = constants, 
                         data = data, inits = inits) 
prepfrogs$initializeInfo() #everything is good to go! 
mcmcfrogs <- configureMCMC(prepfrogs, monitors = params, print = T )
# this tells us that all our params of interest are going to be modeled with RW samplers 
# RW samplers = MetropolisHastings adaptive random-walk sampler
frogsMCMC <- buildMCMC(mcmcfrogs, enableWAIC = TRUE) # a change here to allow us to use WAIC (similar to DIC)
Cmodel <- compileNimble(prepfrogs) #compiling the model itself in C++; it says "may take a minute" but on big models this can take HOURS. 
Compfrogs <- compileNimble(frogsMCMC, project = prepfrogs) # compile the samplers next 
Frog.mod.nimble <- runMCMC(Compfrogs, niter = 30000, thin = 1, nchains =3, nburnin = 10000, samplesAsCodaMCMC = TRUE, WAIC = TRUE)

#Look at the summary and check convergence 
summary(Frog.mod.nimble$samples)
gelman.diag(Frog.mod.nimble$samples)
#looks like the model converged, but we should visually inspect the chains as well 

plot(mcmc.list(Frog.mod.nimble$samples))

#looking good! 

## In NIMBLE, there's no automatic function for DIC. You can calculate it yourself if you choose, but you can instead use the built in WAIC criteria instead. WAIC is very similar to DIC, but is actually more robust - from the literature it seems people tend to like DIC for it's ease of use but WAIC performs better overall. There's no one "right" way to do it, but since WAIC is built in, we'll go with it :) 

Frog.mod.nimble$WAIC #here's our mean deviance, but again, it's not useful if we don't have other models to compare it to. Let's run the other two models that we ran in JAGS. 

### NIMBLE Models # 2 and # 3 #######
#Maybe one frog grows faster than the other, so we would expect a difference in the relationship with age, but no difference in the intercept. We will inspect this in model 2:

nimbleFrogs2 <-
  nimbleCode({ #don't forget this part 
for (i in 1:n.frogs){ 
  meanweight[i] <- beta0 + age[i]*beta1 + leg[i]*beta2 + distance[i]*beta3[s[i]]
  #note that the only thing that changes is the indexing of beta0 and beta3
  weight[i] ~ dnorm(meanweight[i], sd = sd)
}

#but now beta0 is only one value and  beta3 can be two different values
beta0 ~ dunif(-50,50) 
beta1 ~ dunif(-50,50)
beta2 ~ dunif(-50,50)
beta3[1] ~ dunif(-50,50)
beta3[2] ~ dunif(-50,50)

sd ~ dunif(0.0001, 100)
})

#we haven't added any new things to monitor (remember, we were already monitoring beta3), so we can just run everything without having to re-write all the data/params/etc. lines. We do need to rewrite the inits line though, to make sure the dimensions of beta0 and beta3 match our new model. 

inits2 <- list(beta0 = runif(1, -10, 10), beta1 = runif(1,-10,10), beta2 = runif(1,-10,10), beta3 = runif(2,-10, 10), sd = runif(1,.1,100))

prepfrogs <- nimbleModel(code = nimbleFrogs2, constants = constants, 
                           data = data, inits = inits2) 
prepfrogs$initializeInfo()
mcmcfrogs <- configureMCMC(prepfrogs, monitors = params, print = T )
frogsMCMC <- buildMCMC(mcmcfrogs, enableWAIC = TRUE) #actually build the code for those samplers
Cmodel <- compileNimble(prepfrogs) #compiling the model itself in C++; 
Compfrogs <- compileNimble(frogsMCMC, project = prepfrogs) # compile the samplers next 
Frog.mod.nimble2 <- runMCMC(Compfrogs, niter = 30000, thin = 1, nburnin = 15000, nchains = 3, samplesAsCodaMCMC = TRUE, WAIC = TRUE)
summary(Frog.mod.nimble2$samples)
gelman.diag(Frog.mod.nimble2$samples)
#looks like the model converged, but we should visually inspect the chains as well 
plot(Frog.mod.nimble2$samples) 

#Grab WAIC and compare with model # 1
Frog.mod.nimble2$WAIC #this model
Frog.mod.nimble$WAIC #model #1 

#We can see this might be a slightly better model, but not *overwhelmingly* so

nimble.Frogs3 <- 
  nimbleCode({ #don't forget this part 
for (i in 1:n.frogs){ 
  meanweight[i] <- beta0[s[i]] + age[i]*beta1[s[i]]
  #got to reindex everything 
  weight[i] ~ dnorm(meanweight[i], sd = sd)
}

#adjust priors to match
beta0[1] ~ dunif(-50,50) 
beta0[2] ~ dunif(-50,50) 
beta1[1] ~ dunif(-50,50)
beta1[2] ~ dunif(-50,50)

sd ~ dunif(0.0001, 100)
})

#adjust all the data, since the params have changed. constants have not changed
params3 <- c("beta0", "beta1","sd") 
data3 <- list(weight = Frogs$Weight, 
              age = Frogs$Age)

#Initial values are similar to before, but we need to give beta0 2 values instead of 1. 
#reminder that runif arguments are runif(how many, minimum, maximum)
inits3 <- list(beta0 = runif(2, -10, 10), beta1 = runif(2,-10,10), sd = runif(1,.01, 10))

prepfrogs <- nimbleModel(code = nimble.Frogs3, constants = constants, 
                         data = data3, inits = inits3) 
prepfrogs$initializeInfo()
mcmcfrogs <- configureMCMC(prepfrogs, monitors = params3, print = T )
frogsMCMC <- buildMCMC(mcmcfrogs, enableWAIC = TRUE) #actually build the code for those samplers
Cmodel <- compileNimble(prepfrogs) #compiling the model itself in C++; 
Compfrogs <- compileNimble(frogsMCMC, project = prepfrogs) # compile the samplers next 
Frog.mod.nimble3 <- runMCMC(Compfrogs, niter = 30000, thin = 1, nburnin = 15000, nchains = 3, samplesAsCodaMCMC = TRUE, WAIC = TRUE)
summary(Frog.mod.nimble3$samples)
gelman.diag(Frog.mod.nimble3$samples)
#looks like the model converged, but we should visually inspect the chains as well! 
plot(Frog.mod.nimble3$samples) 
Frog.mod.nimble3$WAIC

#Model 2 remains the best model, but let's look at all 3 WAICs in a table together:
data.frame(Mods = 1:3, WAIC = c(Frog.mod.nimble$WAIC, Frog.mod.nimble2$WAIC, Frog.mod.nimble3$WAIC))

#yaaay! 


## graphing from NIMBLE will be very similar to graphing from JAGS so I won't go over it, but I can always tack something on if it's confusing! 

