#Base DS model no covariates 
writeLines("
model
{
	for (i in 1:(nind +nz)) {
		w[i] ~ dbern(psi)					# augmentation
		x[i] ~ dunif(0,Bx)				# distances
		logp[i] <- -((x[i]*x[i])/(2*sigma*sigma))
		p[i] <- exp(logp[i])
		mu[i] <- w[i]*p[i]
		y[i] ~ dbern(mu[i])
	}

	sigma ~ dunif(0, 20)
	psi~ dunif(0,1)			#exists or not

	N <- sum(w)
	D <- N/(2*L*Bx)			#burrow density
}
", "BaseDS.txt")

#Size adjustment model 
writeLines("
model
{
	for (i in 1:(nind +nz)) {
		w[i] ~ dbern(psi)					##augmentation 
		x[i] ~ dunif(0,Bx)					#distance from line for the missed ones; Bx = max(distances) 
		
		z[i] ~ dgamma(a[i],c.beta[i])T(4,50)		                
		a[i] <- shape[clust[i]]		## a depends on which cluster the size comes from
		c.beta[i] <- betaClust[clust[i]]	#precision is allowed to vary between clusters
    clust[i] ~ dcat( pClust[1:Nclust] )	## clusters are defined as categories, assuming Nclust clusters

		sigma[i] <- exp(sigma.int+sigma.beta*z[i])	#log link
		logp[i] <- -((x[i]*x[i])/(2*sigma[i]*sigma[i]))	#This is the half normal distribution with an adjustment for covariates
		p[i] <- exp(logp[i])*xi[i]				# probabilty of seeing the thing, regardless of us actually seeing it; scaled by xi
		xi[i] <- ifelse(z[i] < b.point, m*z[i]+intercept, 1) #if less than b.point, probability is linear; if larger than b.point, perfect detection on line
		
		mu[i] <- w[i]*p[i] 					## probabilty of seeing it for all the ones we DID see (where w[i] = 1)
		y[i] ~ dbern(mu[i])
  
  o[i]~ dbin(o2[i], 1)
  logit(o2[i]) <- o.int + z.beta*z[i]
  ow[i] <- o[i]*w[i] # was this burrow both occupied and real
	}			 		

	intercept ~ dunif(.1,.8)
	b.point ~ dunif(15,25) #what size do we reach perfect detection on the line 
	p.online <- m*4+intercept #assuming min tortoise burrow size is 4 cm wide 
	m <- (1-intercept)/b.point 	#slope for detection on the line for smaller objects		
		
	
	sigma.int~ dnorm(0,s.int.tau)
		s.int.tau <- 1/(s.int.sd*s.int.sd)
		s.int.sd ~ dunif(.00001,5)	
	sigma.beta~ dnorm(0,s.beta.tau)
		s.beta.tau <- 1/(s.beta.sd*s.beta.sd)
		s.beta.sd ~ dunif(.00001,5)
	o.int~ dnorm(0,o.int.tau)
  o.int.tau <- 1/(o.int.sd*o.int.sd)
  o.int.sd ~ dunif(.00001,5)		
  z.beta~ dnorm(0,z.beta.tau)
  z.beta.tau <- 1/(z.beta.sd*z.beta.sd)
  z.beta.sd ~ dunif(.00001,5)			

	for (clustIdx in 1: Nclust) {
		shape[clustIdx] ~ dunif(1,100)
     		betaClust[clustIdx] ~ dunif(.2,2)
		}
  pClust[1:Nclust] ~ ddirch(psizes) ## probability of each cluster = the probability of that category in the ddirch distribution
	
															
	psi~ dunif(0,1)			#exists or not		
	N <- sum(w)		
	D <- N/(2*L*Bx)			#burrow density
	Nt <- sum(wo)			 
	Dt <- Nt/(2*L*Bx)			#tort density
}
", "Covmod.txt")

# Veg Model

writeLines("
model
{
  for (i in 1:(nind +nz)) {
  w[i] ~ dbern(psi)					#augmentation 
  x[i] ~ dunif(0,Bx)					#distance from line; Bx = max(distances) 
  
  z[i] ~ dgamma(a[i],c.beta[i])T(4,50)		#covariate                
  a[i] <- shape[clust[i]]		## a depends on which cluster the size comes from
  c.beta[i] <- betaClust[clust[i]]	#precision is allowed to vary between clusters
  clust[i] ~ dcat( pClust[1:Nclust] )	## clusters are defined as categories, assuming Nclust clusters
  
  v[i] ~ dbeta(d,e)T(.05,.95)				## vegetation is unknown so beta prior
  sigma[i] <- exp(sigma.int+sigma.beta*z[i]+sigma.gamma*v[i])  #sigma of detection
  logp[i] <- -((x[i]*x[i])/(2*sigma[i]*sigma[i]))	#This is the normal distribution with an adjustment for covs
  p[i] <- exp(logp[i])*xi[i]
  xi[i] <- ifelse(z[i] < b.point, m*z[i]+intercept, 1) #if less than b.point, probability is linear; if larger than b.point, perfect detection on line
  
  mu[i] <- w[i]*p[i] 					# probabilty of seeing it for all the ones we DID see (where w[i] = 1)
  
  y[i] ~ dbern(mu[i])         #found vs. missed
  
  o[i]~ dbin(o2[i], 1)        #occupancy
  logit(o2[i]) <- o.int + z.beta*z[i]
  wo[i] <- o[i]*w[i]  # real tortoises
  woz[i] <- o[i]*w[i]*z[i] #real tortoise sizes
  }
  
  for (q in 1:q) {
    v2[q] ~ dbeta (d,e)T(.05,.95)		# all veg measurements at the site
  }
  
  intercept ~ dunif(.1,.8)
	b.point ~ dunif(15,25) #what size do we reach perfect detection on the line 
	p.online <- m*4+intercept #assuming min tortoise burrow size is 4 cm wide 
	m <- (1-intercept)/b.point 	#slope for detection on the line for smaller objects		
	
  
  sigma.int~ dnorm(0,s.int.tau)
  s.int.tau <- 1/(s.int.sd*s.int.sd)
  s.int.sd ~ dunif(.00001,5)	
  sigma.beta~ dnorm(0,s.beta.tau)
  s.beta.tau <- 1/(s.beta.sd*s.beta.sd)
  s.beta.sd ~ dunif(.00001,5)
  o.int~ dnorm(0,o.int.tau)
  o.int.tau <- 1/(o.int.sd*o.int.sd)
  o.int.sd ~ dunif(.00001,5)		
  z.beta~ dnorm(0,z.beta.tau)
  z.beta.tau <- 1/(z.beta.sd*z.beta.sd)
  z.beta.sd ~ dunif(.00001,5)	
  sigma.gamma~ dnorm(0,s.gamma.tau)
  s.gamma.tau <- 1/(s.gamma.sd*s.gamma.sd)
  s.gamma.sd ~ dunif(.00001,5)
  
  d~dunif(.1,40)
  e~dunif(.1,40)
  
  for (clustIdx in 1: Nclust) {
  shape[clustIdx] ~ dunif(1,80)
  betaClust[clustIdx] ~ dunif(0.2,2)
  }
  
  pClust[1:Nclust] ~ ddirch(onesRepNclust) # probability of each cluster = the probability of that category in the ddirch distribution
  
  
  psi~ dunif(0,1)			#exists or not		
  
  N <- sum(w)		
	D <- N/(2*L*Bx)			#burrow density
	Nt <- sum(wo)			 
	Dt <- Nt/(2*L*Bx)			#tort density

  juvi1 <- sum(woz < 13)/Nt  
  juvi2 <-  (sum(woz < 21)- sum(woz <13))/Nt
  juvi3 <- sum(woz > 21)/Nt
}
", "VegMod.txt")

#Point Count Model (run in NIMBLE)

writeLines("
  nimbleCode({
    for (k in 1:n.species){
      for (i in 1:n.sites) {
        n.gof[i,k] <- n[i,k,8]
      }}
    
    
    for (i in 1:n.sites) {
      for (k in 1:n.species){
        log(psi[i,k,1]) <- psi.b2[k] +  psi.b3[k]*temp[i,1] + psi.b4[k]*precip[i,1]
        abund[i,k,1] ~dpois(psi[i,k,1])
        #state_like[i,k,1] <- dpois(abund[i,k,1], psi[i,k,1], log = T)
        for (q in 2:n.years){
          
          psi[i,k,q] <- psi[i,k,q-1]*phi[i,k,q] + psi[i,k,q-1]*gamma[i,k,q]
          abund[i,k,q] <- S[i,k,q] + R[i,k,q]
          
          S[i,k,q] ~ dbin(phi[i,k,q], abund[i,k,q-1]) #reductions from pop
          #log_S[i,k,q] <- dbinom(S[i,k,q],  abund[i,k,q-1],phi[i,k,q], log = T)
          
          R[i,k,q] ~ dpois(psi[i,k,q-1]*gamma[i,k,q])
          #log_R[i,k,q] <- dpois(R[i,k,q], psi[i,k,q-1]*gamma[i,k,q], log = T)
          
          log(gamma[i,k,q]) <- gamma.b2[k] + gamma.b3[k]*temp[i,q] + gamma.b4[k]*precip[i,q]
          
          
          logit(phi[i,k,q]) <- phi.b2[k] + phi.b3[k]*temp[i,q] + phi.b4[k]*precip[i,q]
          
          #state_like[i,k,q] <- log_S[i,k,q]+log_R[i,k,q]
          
        } #end q
      } #end k
    } #end i
    
    #detection loop
    
    
    
    for (k in 1:n.species){
      for (i in 1:n.sites) {
        for (q in 1:(n.years-1)){ #note change for GOF
          #loglike[i,k,q] <- like_obs[i,k,q]+like_time[i,k,q]+like_n[i,k,q]+state_like[i,k,q]
          
          
          
          obs.pis[i,k,q,1:bin.dist] <- p.dist.c[i,k,q,1:bin.dist]*surv[i,q]
          obs[i,k,q,1:bin.dist] ~ dmulti(obs.pis[i,k,q,1:bin.dist], n[i,k,q]) #what distance was it at when observed 
          #like_obs[i,k,q] <- dmulti(obs[i,k,q,1:bin.dist],  n[i,k,q],obs.pis[i,k,q,1:bin.dist], log = T)
          
          timeclass.pis[i,k,q,1:4] <- pi.pa.c[i,k,q,1:4]*surv[i,q]
          timeclass[i,k,q,1:4] ~ dmulti(timeclass.pis[i,k,q,1:4], n[i,k,q]) #what time interval first observed 
          #like_time[i,k,q] <- dmulti(timeclass[i,k,q,1:4],  n[i,k,q],timeclass.pis[i,k,q,1:4], log = T)
        } #end q for GOF
        
        for (q in 1:(n.years)){ #for GOF
          
          n[i,k,q] ~ dbin(p.dist.sum[i,k,q]* pi.pa.sum[i,k,q]*surv[i,q], abund[i,k,q]) #number detected
          #like_n[i,k,q] <- dbinom(n[i,k,q],  abund[i,k,q],p.dist.sum[i,k,q]* pi.pa.sum[i,k,q]*surv[i,q],log = T)
          #
          
          p.dist.sum[i,k,q] <- sum(pi[i,k,q,1:bin.dist])
          
          sigma[i,k,q] <- exp(p.b0 + p.b4[observer[i,q]]+p.b2*noise[i,q])
          for (b in 1:bin.dist){
            pbar[i,k,q,b] <- (sigma[i,k,q]^2 * (1- exp(-bb[b+1]^2/(2*sigma[i,k,q]^2))) - 
                                sigma[i,k,q]^2 * (1- exp(-bb[b]^2/(2*sigma[i,k,q]^2)))) *
              2 * 3.141593/area[b]
            pi[i,k,q,b] <- ps[b]*pbar[i,k,q,b]
            
            p.dist.c[i,k,q,b] <- (pi[i,k,q,b] / p.dist.sum[i,k,q]) #conditional probability 
            
          } #end b
          
          
          logit(p.avail[i,k,q]) <-  p.avail.b0 + p.avail.b1[k] + p.avail.b2*time[i,q] #time availability
          pi.pa.sum[i,k,q] <- sum(pi.pa[i,k,q,1:4])
          
          for (t in 1:4){
            pi.pa[i,k,q,t]  <- (p.avail[i,k,q])* pow(1-(p.avail[i,k,q]), (t-1)) 
            pi.pa.c[i,k,q,t] <- pi.pa[i,k,q,t]/pi.pa.sum[i,k,q]
          } #end t
        } #end i
      } #end k
      
    } #end q
    
    for (k in 1:n.species){
      p.avail.b1[k] ~ dnorm(0, tau.avail[group[k]]) 
      psi.b2[k] ~  dnorm(psi.bar[group[k]], tau.psib2[group[k]])
      psi.b3[k] ~ dnorm(0, 1)
      psi.b4[k] ~ dnorm(0, 1)
      gamma.b2[k] ~ dnorm(gamma.bar[group[k]], tau.gammab2[group[k]])
      gamma.b3[k] ~ dnorm(0, 1)
      gamma.b4[k] ~ dnorm(0, 1)
      gamma.b5[k] ~ dnorm(0, 1)
      phi.b2[k] ~  dnorm(phi.bar[group[k]], tau.phib2[group[k]])
      phi.b3[k] ~ dnorm(0, 1)
      phi.b4[k] ~ dnorm(0, 1)
      phi.b5[k] ~ dnorm(0, 1)
      
      
    } #end k 
    
    
    ## Group Priors
    
    for (m in 1:2){
      sd.avail[m]  ~ dexp(1)
      tau.avail[m] <- 1/(sd.avail[m]^2)
      sd.psib2[m]  ~ dexp(1) #dunif(.005, .5)
      sd.phib2[m]  ~ dexp(1) #dunif(.005, .5)
      sd.gammab2[m]  ~ dexp(1) #dunif(.005, .5)
      tau.psib2[m]  <- 1/(sd.psib2[m])^2
      tau.phib2[m]  <- 1/(sd.phib2[m])^2
      tau.gammab2[m]  <- 1/(sd.gammab2[m])^2
      
      psi.bar[m] ~ dnorm(0, 1)
      gamma.bar[m] ~ dnorm(0, 1)
      phi.bar[m] ~ dnorm(0, 1)
    }
    
    p.avail.b0 ~ dnorm(0, 1)
    p.b0 ~ dnorm(0, 1) #intercept
    p.b2 ~ dnorm(0, 1) #noise 
    p.avail.b2 ~ dnorm(0, 1) #time of survey
    
    
    for (y in 1:n.observers){
      p.b4[y] ~ dnorm(0, obs.prec)
    }
    obs.prec <- 1/(obs.sd*obs.sd)
    obs.sd ~ dunif(0,2)
    
  })
", "PointCountmod.txt")