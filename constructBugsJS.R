cat("model{
  #phi survival probablility
  #gamma removal entry probability
  #p capture probability
  
  #States: 1=not yet entered, 2=alive, 3=dead
  #Observations: 1=seen, 2=not seen
  
  #Priors and constraints

  for( s in 1:nSamples ){     
      for( r in 1:(nRivers) ){  

         phi[ s,r ] <- meanPhi
         gamma[ s,r ] ~ dunif(0,1)
         p [ s,r ] <- meanP
    
      }
    }
  meanPhi ~ dunif(0,1)
  meanP ~ dunif(0,1)

  
  #Define state-transition and observation matrices

    #State transition
  for( s in 1:nSamples ){     
      for( r in 1:(nRivers) ){ 
        ps[1,s,r,1]<-1-gamma[s,r]
        ps[1,s,r,2]<-gamma[s,r]
        ps[1,s,r,3]<-0
        ps[2,s,r,1]<-0
        ps[2,s,r,2]<-phi[s,r]
        ps[2,s,r,3]<-1-phi[s,r]
        ps[3,s,r,1]<-0
        ps[3,s,r,2]<-0
        ps[3,s,r,3]<-1
      
      #Observation
        po[1,s,r,1]<-0
        po[1,s,r,2]<-1
        po[2,s,r,1]<-p[s,r]
        po[2,s,r,2]<-1-p[s,r]
        po[3,s,r,1]<-0
        po[3,s,r,2]<-1
      }#r
    }#s
  
  #Likelihood
  #Define latent state at first occasion
#   for(i in 1:nFirstObsRows){
#     z[firstObsRows[i]]<-1 #All individuals start in unentered state
#   }

  for(i in 1:nEvalRows){
    #State process
    z[evalRows[i]]~dcat(ps[z[evalRows[i]-1],
                        sample[evalRows[i]],
                        river[evalRows[i]-1],])
    #observation process
    encDATA[evalRows[i]]~dcat(po[z[evalRows[i]],
                                 sample[evalRows[i]],
                                 river[evalRows[i]],])
  }#i
  
  #Calculate derived population parameters
    for(s in 1:nSamples){
      for( r in 1:(nRivers) ){ 
        qgamma[s,r]<-1-gamma[s,r]
      }
    }
  
  for(r in 1:nRivers){
    cprob[1,r]<-1
  }
   
  for(s in 2:nSamples){
    for( r in 1:(nRivers) ){ 
      cprob[s,r]<-gamma[s,r]*prod(qgamma[1:(s-1),r])
    }
  }

  for(r in 1:nRivers){
    psi[r]<-sum(cprob[,r]) #Inclusion probability
  }

  for(s in 2:nSamples){
    for( r in 1:(nRivers) ){
      b[s,r]<-cprob[s,r]/psi[r]
    }
  }

#   for(i in 1:M){
#     for(t in 2:n.occasions){
#       al[i,t-1]<-equals(z[i,t],2)
#     }
#     for(t in 1:(n.occasions-1)){
#       d[i,t]<-equals(z[i,t]-al[i,t],0)
#     }
#     alive[i]<-sum(al[i,])
#   }
#   for(t in 1:(n.occasions-1)){
#     N[t]<-sum(al[,t])
#     B[t]<-sum(d[,t])
#   }
#   for(i in 1:M){
#     w[i]<-1-equals(alive[i],0)
#   }
# Nsuper<-sum(w[])
}",file="bugsJS.txt")