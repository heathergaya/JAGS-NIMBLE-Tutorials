################ Linear Regression Part 1 by Heather Gaya.  ########
###### Last updated March 19, 2020.                          #######
###### Email heather.e.gaya@gmail.com for questions          #######
###### find me on Twitter @doofgradstudent                   ##### 


#Linear Regression Part 1! Whoo, everyone's favorite. 
# This code assumes you already have JAGS and/or NIMBLE installed on your computer (they are standalone programs)


# Quick reminders about R - 
# the "#" symbol turns things into comments (instead of code)
# set your working directory! Session -> set working directory --> (to where your files are)


### A fake scenario as an example #####
# A researcher, let's call him Joe, goes out into the world to study frogs
# He collects 50 frogs and records some measurements on them 
# specifically, we records:
# age (in days, shhhh don't worry about how he knows)
# weight (in g)
# back leg length (cm)
# distance from road (m)
# and species (A or B)

# Now he is back at home, sitting at his computer, trying to see if the weight of frogs is related to any of the other variables. We shall help him! 

#First, I set my working directory like a good programmer
setwd("~/Desktop/U Georgia/R_Workshop/JAGS:NIMBLE/Linear_Regression_1") #you will need to put your own file path in here

# fetch the data 
Frogs <- read.csv("Fakefrogs.csv")
head(Frogs)

#R packages we will be using to analyze 
# you will need to install these first if you haven't 
# install.packages("name-of-your-package") will do the trick :) 

library(runjags) #my preferred JAGS package
library(nimble) # only need this if you want to use NIMBLE 
library(parallel) # for faster computing
library(coda) #for inspecting chains


# Joe suspects weight is related to age, back leg length and distance, so he chooses those three variables to model. 
# Essentially, joe thinks:
# Expected weight of frog <- intercept + age*something + leg*something2 + distance*something3
#  and the actual weight of frog is from normal distribution with mean = expected weight and standard deviation = some sd
#calling it "something" is not really accepted in scientific journals, so we use "beta" instead. The intercept is often referred to as beta0 


### Writing the JAGS model; continuous variables only ###### 
#everything inside these quotes is now written in the BUGS language. NOT in R. 
#this means the order of things DOES NOT MATTER
# with a few v. small caveats...

modelstring.Frogs = "
  model 
{
for (i in 1:n.frogs){ 

#we want to iterate through every frog in our data. We do this using a for loop.
  meanweight[i] <- beta0 + age[i]*beta1 + leg[i]*beta2 + distance[i]*beta3 
  weight[i] ~ dnorm(meanweight[i], prec)
  #in JAGS, d = distribution and norm = the normal distribution
  #the indexing is to let JAGS know what value to use for each frog
  #note that dnorm requires two arguments, the mean (in this case meanweight[i]) and the PRECISION which is = 1/variance = 1/(standard_deviation ^2).
  
}
  
#Yay, all well and good... but what are the betas? And what is prec? 
# In bayes, we like draw unknowns from distributions. 
# In this case, we have no idea what the betas are, so we just give them really vague priors. 
# the ~ means it's coming from a distribution. 
#in JAGS, d = distribution and unif = uniform, so we're saying the value can be anything from -50 to 50, we have no idea 

beta0 ~ dunif(-50,50) 
beta1 ~ dunif(-50,50)
beta2 ~ dunif(-50,50)
beta3 ~ dunif(-50,50)

#I like to convert precision to standard deviation because it confuses me otherwise. 
prec <- 1/(sd *sd)
sd ~ dunif(0.0001, 100) #can't be exactly 0 or precision = 1/0 and it gets angry.

#awesome! we have now defined everything. Excpet the value n.frogs, which we can give to JAGS outside the model string.  
  
}
" 


#note we haven't run anything yet, we just made a model that can be run later. 
# time to do the running part. 


############# Running the model (JAGS) ############ 
# In JAGS, we only receive information back on the parameters that we monitor. So in this case, we want estimates for all the betas and the sd, so we monitor those. We don't need to monitor the data variables (distance, age, etc.) because we already know those values! 
params <- c("beta0", "beta1", "beta2", "beta3","sd")

# we also have to tell JAGS our data. This is given as variable-name-in-model = the-data-for-that-variable. You have to give it to jags as a list. 
#reminder that in R, you can use "$" and the column name to select a column from a data frame. So Frogs$Weight takes all the weight values from our data and puts it in a vector
#we also need to tell JAGS what n.frogs is (the total # of frogs)
data <- list(weight = Frogs$Weight, distance = Frogs$Dist, 
          age = Frogs$Age, leg = Frogs$Leg,
          n.frogs = 50)
#note that we do not give it anything for "expectedweight" because that is something it will already solve for. It will yell if you try to provide it values and ask it to solve for it. 

#the other thing we want to give JAGS are initial values. If you have some idea of the expected values, you can specify them here too. But if you don't know anything, you can choose random values within the bounds that you set up above.
#so for instance let's give it values for the betas that are between -10 and 10. This won't affect the answer, just gives it a starting place 

inits <- function(){list(beta0 = runif(1, -10, 10), beta1 = runif(1,-10,10), beta2 = runif(1,-10,10), beta3 = runif(1,-10, 10))}
# The "function" means it will generate a new starting value for each chain that you run
# In general, chains are just there to help you with processing. Your comptuer has (at least) 4 cores, so using 3 of them for parallel processing (aka 3 chains) will help you find your answer with maximum efficiency! There are better explanations of this out there, but basically use at least 2 or 3 if you want to be able to estimate if your chains have converged, aka, have reached an answer. 

#Okay, time to run! 
Frog.mod <- run.jags(model = modelstring.Frogs, 
                     monitor = params, inits = inits,
                     data = data, n.chains = 3, 
                     sample = 10000, method = "parallel")
#you may see a warning about You attempted to start parallel chains without setting different PRNG... etc etc. Ignore this warning

#now we look at what our output says:
Frog.mod

#lower95 and upper95 are the bounds of your credible interval 
#median and mean are the median and the mean
# psrf is a convergence diagnositic. If this number is larger than ~1.1, you need to run the model for more time (more samples)

#let's plot the results as well

plot(Frog.mod)

#the final plot (the first you'll see) is about correlation between parameters. 
# the next plot you scroll to is SD
# the top left panel is the chains. You want this to look like a blur of 3 colors (from your 3 chains). This means the algorithm is exploring the mathematical landscape correctly and the chains have converged! Yee! 
# the ECDF shows a cumulative density plot of all the values selected for sd. you can see most of the values are less than 5, since the cdf starts to hit almost 100% around 5. 
#the left bottom plot shows you a histogram of the possible values drawn for sd. In this case, you can see the mean is about 4.2 and most of the values are between ~ 4 and ~6.5. 
# The last plot is autocorrelation. This tells you how far apart samples have to be before they are "independent" draws. Essentially you want this to be near 0 for most of the lag values. A lag of 5 indicates 5 iterations apart, a lag of 20 = 20 iterations apart, etc. etc. 

# we will discuss the results after we try out the NIMBLE method! 

########## The NIMBLE model; continuous variables only ##########
#notice the different structure of the modelstring. instead of "model {" at the beginning, we use "nimbleCode({" 

nimbleFrogs <- 
  nimbleCode({ #so this is a new change 
for (i in 1:n.frogs){ 

  meanweight[i] <- beta0 + age[i]*beta1 + leg[i]*beta2 + distance[i]*beta3 
  weight[i] ~ dnorm(meanweight[i], sd = sd)
  # note that NIMBLE can use sd, var OR precision for dnorm! You can choose what you like best. tau = precision, sd = standard deviation, var = variance. Whoo flexible!
  
}

beta0 ~ dunif(-50,10) 
beta1 ~ dunif(-50,10)
beta2 ~ dunif(-50,10)
beta3 ~ dunif(-50,10)

sd ~ dunif(0.0001, 100) 

#we don't need precision but otherwise not much changes when converting to NIMBLE! 

}) 

############# Running the model (NIMBLE) ############ 
#Pros of NIMBLE - fast!!! Very fast! same language as JAGS. googlegroup for NIMBLE is small and people respond really fast. 
# Cons - compiles in C, so errors can be weird or hard to understand; less documentation on the internet; more code to type out to get it to run. 
# In NIMBLE, we only receive information back on the parameters that we monitor. So in this case, we want estimates for all the betas and the sd, so we monitor those. We don't need to monitor the data variables (distance, age, etc.) because we already know those values! 
params <- c("beta0", "beta1", "beta2", "beta3","sd")

#NIMBLE wants constants (things that don't have the ~ sign associated with them) in a separate list. In this case we only have one constant - n.frogs
constants <- list(n.frogs = 50)

#now we tell NIMBLE the data as well 
data.frogs <- list(weight = Frogs$Weight, distance = Frogs$Dist, 
             age = Frogs$Age, leg = Frogs$Leg)
#note that we do not give it anything for "expectedweight" because that is something it will already solve for.


#the other thing we want to give NIMBLE are initial values. 
#so for instance in our model code we said the betas are between -10 and 10. 
inits <- list(beta0 = runif(1, -10, 10), beta1 = runif(1,-10,10), beta2 = runif(1,-10,10), beta3 = runif(1,-10, 10), sd = runif(1,0, 100))

#Okay, time to run! 
#NIMBLE can be run all at once but I prefer to do it in steps so you have control and it's not a scary black box. 

prepfrogs <- nimbleModel(code = nimbleFrogs, constants = constants, 
                         data = data.frogs, inits = inits) 
#notice that you don't give it params to monitor yet, you will do that later. 
#you may get a lot of things showing up in the screen - that's okay! 
#let's make sure everything is initialized and ready to go
prepfrogs$initializeInfo() #everything is good to go! 
mcmcfrogs <- configureMCMC(prepfrogs, monitors = params, print = T )
# this tells us that all our params of interest are going to be modeled with RW samplers 
# RW samplers = MetropolisHastings adaptive random-walk sampler
# this is what we want. For now don't worry too much about this, I will hopefully explain MCMC samplers at a later point. 

frogsMCMC <- buildMCMC(mcmcfrogs) #actually build the code for those samplers
Cmodel <- compileNimble(prepfrogs) #compiling the model itself in C++; it says "may take a minute" but on big models this can take HOURS. 
Compfrogs <- compileNimble(frogsMCMC, project = prepfrogs) # compile the samplers next 
Frog.mod.nimble <- runMCMC(Compfrogs, niter = 15000, thin = 1, nchains =3, nburnin = 1000, samplesAsCodaMCMC = TRUE)
# in NIMBLE, you CANNOT EXTEND RUNS. This is very annoying. But alas, you cannot. Also note that unlike JAGS, the samples you get out = niter - nburnin. In JAGS they are separate things, but in NIMBLE you are throwing away the first nburnin iterations from your data. samplesAsCodaMCMC just makes the output nice for later :) 

#IF you want to run this all in parallel, instead of doing chain 1, then chain 2, then chain 3 you can send the output to your computers cores via: 

cl <- makeCluster(3) #tell the computer you're going to use 3 cores
clusterExport(cl = cl, varlist = c("constants", "data.frogs", "inits", "params", "nimbleFrogs")) #send all the information to each core to prepare 
system.time(frog.out <- clusterEvalQ(cl = cl,{
  library(nimble) #you're now in a totally different environment so have to load the package again
  prepfrogs <- nimbleModel(code = nimbleFrogs, constants = constants, 
                           data = data.frogs, inits = inits) 
  prepfrogs$initializeInfo()
  mcmcfrogs <- configureMCMC(prepfrogs, monitors = params, print = T )
  frogsMCMC <- buildMCMC(mcmcfrogs) #actually build the code for those samplers
  Cmodel <- compileNimble(prepfrogs) #compiling the model itself in C++; 
  Compfrogs <- compileNimble(frogsMCMC, project = prepfrogs) # compile the samplers next 
  Frog.mod.nimble <- runMCMC(Compfrogs, niter = 15000, thin = 1, nburnin = 1000, samplesAsCodaMCMC = TRUE)
}))
#elapsed tells you the time it took to run in seconds
# on my computer, about 30 seconds 

Frog.mod.nimble <- mcmc.list(frog.out)
stopCluster(cl) #close the parallel computing or it will take up memory/space on your computer 

#now we look at what our output says:
summary(Frog.mod.nimble)

gelman.diag(Frog.mod.nimble)
#Reminder that we want these values to all be under ~1.1 



plot(mcmc.list(Frog.mod.nimble))
#same deal with the 3 color plots, but you'll notice the output is simpler. The plots on the right are a smoothed version of a histogram of the values from all the chains.


############## Looking at the Results ###############
Frog.mod
summary(Frog.mod.nimble)

#There are slight differences in the output between the two programs, but essentially we found the same results. 

#Recall that our original model (in english) was:
#Expected weight of frog <- intercept + age*something + leg*something2 + distance*something3
#  and the actual weight of frog is from normal distribution with mean = expected weight and standard deviation = some sd

#Now we have:
#Expected weight = -8.0 + age*1.2 - leg * 2.04 + distance*.17
#and actual weight ~ normal(mean = Expected weight, sd = 5.4)

#Let's make some plots to see what this looks like compared to our data. 

age = 0:800 #range of ages (in days) we want to graph

#for 3.5 cm legged frogs at 60 m from the road, what's the expected weight based on age?
dev.off()
plot(age, -8+age*.12 - 3.5*2.04 + 60*.17, type = "l")
points(Frogs$Age, Frogs$Weight)

#Not too terrible a fit! 


