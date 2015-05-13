cat("model{
  #phi survival probablility
  #gamma removal entry probability
  #p capture probability
  
  #States: 1=not yet entered, 2=alive, 3=dead
  #Observations: 1=seen, 2=not seen
  
  #Priors and constraints

  for( s in 1:nSamples ){     
      for( r in 1:(nRivers) ){  

         phiYoy[ s,r ] ~ dunif(0,1)
         phiAdult[ s,r ] ~ dunif(0,1)
         gamma[ s,r ] ~ dunif(0,1)
         pYoy [ s,r ] ~ dunif(0,1)
         pAdult [ s,r ] ~ dunif(0,1)
    
      }
    }

  
  #Define state-transition and observation matrices

    #State transition for non-summer (YOY stay YOY)
  for( s in 1:nNonSummerSamples ){     
      for( r in 1:(nRivers) ){ 
        ps[1,nonSummerSamples[s],r,1]<-1-gamma[nonSummerSamples[s],r]
        ps[1,nonSummerSamples[s],r,2]<-gamma[nonSummerSamples[s],r]
        ps[1,nonSummerSamples[s],r,3]<-0
        ps[1,nonSummerSamples[s],r,4]<-0
        ps[2,nonSummerSamples[s],r,1]<-0
        ps[2,nonSummerSamples[s],r,2]<-phiYoy[nonSummerSamples[s],r]
        ps[2,nonSummerSamples[s],r,3]<-0
        ps[2,nonSummerSamples[s],r,4]<-1-phiYoy[nonSummerSamples[s],r]
        ps[3,nonSummerSamples[s],r,1]<-0
        ps[3,nonSummerSamples[s],r,2]<-0
        ps[3,nonSummerSamples[s],r,3]<-phiAdult[nonSummerSamples[s],r]
        ps[3,nonSummerSamples[s],r,4]<-1-phiAdult[nonSummerSamples[s],r]
        ps[4,nonSummerSamples[s],r,1]<-0
        ps[4,nonSummerSamples[s],r,2]<-0
        ps[4,nonSummerSamples[s],r,3]<-0
        ps[4,nonSummerSamples[s],r,4]<-1
      }
  }
  
    #State transition for summer (YOY mature)
  for( s in 1:nSummerSamples ){     
      for( r in 1:(nRivers) ){ 
        ps[1,summerSamples[s],r,1]<-1-gamma[summerSamples[s],r]
        ps[1,summerSamples[s],r,2]<-gamma[summerSamples[s],r]
        ps[1,summerSamples[s],r,3]<-0
        ps[1,summerSamples[s],r,4]<-0
        ps[2,summerSamples[s],r,1]<-0
        ps[2,summerSamples[s],r,2]<-0
        ps[2,summerSamples[s],r,3]<-phiYoy[summerSamples[s],r]
        ps[2,summerSamples[s],r,4]<-1-phiYoy[summerSamples[s],r]
        ps[3,summerSamples[s],r,1]<-0
        ps[3,summerSamples[s],r,2]<-0
        ps[3,summerSamples[s],r,3]<-phiAdult[summerSamples[s],r]
        ps[3,summerSamples[s],r,4]<-1-phiAdult[summerSamples[s],r]
        ps[4,summerSamples[s],r,1]<-0
        ps[4,summerSamples[s],r,2]<-0
        ps[4,summerSamples[s],r,3]<-0
        ps[4,summerSamples[s],r,4]<-1
      }
    }
    
    for(s in 1:nSamples){
      for(r in 1:nRivers){
      #Observation
        po[1,s,r,1]<-0
        po[1,s,r,2]<-1
        po[2,s,r,1]<-pYoy[s,r]
        po[2,s,r,2]<-1-pYoy[s,r]
        po[3,s,r,1]<-pAdult[s,r]
        po[3,s,r,2]<-1-pAdult[s,r]
        po[4,s,r,1]<-0
        po[4,s,r,2]<-1
      }#r
    }#s

  #Likelihood

#Passed as data so commented out
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
   for( r in 1:(nRivers) ){     
      for(s in 1:nSamples){
        qgamma[s,r]<-1-gamma[s,r]
      }
    }
  
  for(r in 1:nRivers){
    cprob[1,r]<-gamma[1,r]

    for(s in 2:nSamples){
      cprob[s,r]<-gamma[s,r]*prod(qgamma[1:(s-1),r])
    }
  }

  for(r in 1:nRivers){
    psi[r]<-sum(cprob[2:nSamples,r]) #Inclusion probability

    for(s in 1:nSamples){
      b[s,r]<-cprob[s,r]/psi[r]
    }
  }

  for(r in 1:nRivers){
    for(s in 1:nSamples){
      for(n in 1:nSampleRows[s,r]){
        yoyAlive[n,s,r]<-equals(z[sampleRows[n,s,r]],2)
        adultAlive[n,s,r]<-equals(z[sampleRows[n,s,r]],3)
      }
      nYoy[s,r]<-sum(yoyAlive[1:nSampleRows[s,r],s,r])
      nAdult[s,r]<-sum(adultAlive[1:nSampleRows[s,r],s,r])
    }
  }
}",file="bugsJS.txt")

