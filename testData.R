#phi survival probablility
#gamma removal entry probability
#p capture probability

#States: 1=not yet entered, 2=alive, 3=dead
#Observations: 1=seen, 2=not seen

N<-200

nSamples<-20
phiYoy<-0.7
phiAdult<-0.9
gamma<-c(0.20,rep(0.11,nSamples-1))
pYoy<-0.5
pAdult<-0.8

#state transition matrix
ps<-array(NA,dim=c(4,nSamples,4))
ps[1,,1]<-1-gamma
ps[1,,2]<-gamma
ps[1,,3]<-0
ps[1,,4]<-0
ps[2,,1]<-0
ps[2,,2]<-0
ps[2,,3]<-phiYoy
ps[2,,4]<-1-phiYoy
ps[3,,1]<-0
ps[3,,2]<-0
ps[3,,3]<-phiAdult
ps[3,,4]<-1-phiAdult
ps[4,,1]<-0
ps[4,,2]<-0
ps[4,,3]<-0
ps[4,,4]<-1

#Observation state matrix
po<-matrix(NA,ncol=2,nrow=4)
po[1,1]<-0
po[1,2]<-1
po[2,1]<-pYoy
po[2,2]<-1-pYoy
po[3,1]<-pAdult
po[3,2]<-1-pAdult
po[4,1]<-0
po[4,2]<-1


data<-data.table(id=rep(1:N,each=nSamples+1),sample=rep(0:nSamples,N),state=1)
evalRows<-which(data$sample!=0)

for(i in evalRows){
  data$state[i]<-which(rmultinom(1,1,ps[data$state[i-1],data$sample[i],])==1)
}

for(i in 1:nrow(data)){
  data$enc[i]<-which(rmultinom(1,1,po[data$state[i],])==1)
}

data[,z:=as.numeric(NA)]
data[enc==1,z:=state]

knownZ<-function(z,state){
  if(any(!is.na(z))){
  lastKnown<-max(which(!is.na(z)))
  z[1:lastKnown]<-state[1:lastKnown]
  }
  return(z)
}

data[,z:=knownZ(z,state),by=id]
data[,sample:=sample+1]

nExtras<-round(data[,length(unique(id))]*0.1)
aug<-data.table(id=NA,sample=rep(1:max(data$sample),nExtras),
                state=NA,enc=2,z=NA)
data<-rbind(data,aug)

data[sample==1,z:=1]

nSampleRows<-rep(NA,nSamples+1)
sampleRows<-array(NA,dim=c(20000,nSamples+1))
for(s in 1:(nSamples+1)){
  rows<-data[,which(sample==s)]
  if(length(rows)>0){
    sampleRows[1:length(rows),s]<-rows
  }
}

nSampleRows<-apply(sampleRows,2,function(x){return(length(na.omit(x)))})
sampleRows<-sampleRows[1:max(nSampleRows),]

winData<-list(z=data$z,
              enc=data$enc,
              sample=data$sample,
              evalRows=which(data$sample!=1),
              nEvalRows=length(which(data$sample!=1)),
              nSamples=length(unique(data$sample)),
              sampleRows=sampleRows,
              nSampleRows=nSampleRows)



cat("model{
  phiAdult~dunif(0,1)
  phiYoy~dunif(0,1)
  pAdult~dunif(0,1)
  pYoy~dunif(0,1)
  
  for(s in 1:nSamples){
    gamma[s]~dunif(0,1)
  }
  
  for(s in 1:nSamples){
  #state transition probabilities
  ps[1,s,1]<-1-gamma[s]
  ps[1,s,2]<-gamma[s]
  ps[1,s,3]<-0
  ps[1,s,4]<-0
  ps[2,s,1]<-0
  ps[2,s,2]<-0
  ps[2,s,3]<-phiYoy
  ps[2,s,4]<-1-phiYoy
  ps[3,s,1]<-0
  ps[3,s,2]<-0
  ps[3,s,3]<-phiAdult
  ps[3,s,4]<-1-phiAdult
  ps[4,s,1]<-0
  ps[4,s,2]<-0
  ps[4,s,3]<-0
  ps[4,s,4]<-1
}
  
  #Observation state matrix
  po[1,1]<-0
  po[1,2]<-1
  po[2,1]<-pYoy
  po[2,2]<-1-pYoy
  po[3,1]<-pAdult
  po[3,2]<-1-pAdult
  po[4,1]<-0
  po[4,2]<-1
  
  
  for(i in 1:nEvalRows){
   z[evalRows[i]]~dcat(ps[z[evalRows[i]-1],sample[evalRows[i]],])
   enc[evalRows[i]]~dcat(po[z[evalRows[i]],]) 
  }

  #Calculate derived population parameters
  for(s in 2:nSamples){
      qgamma[s-1]<-1-gamma[s]
  }
  

    cprob[1]<-gamma[2]
    
    for(s in 3:nSamples){
      cprob[s-1]<-gamma[s]*prod(qgamma[1:(s-1)])
  }
  

    psi<-sum(cprob[1:(nSamples-1)]) #Inclusion probability
    
    for(s in 2:nSamples){
      b[s-1]<-cprob[s-1]/psi
    }



  for(s in 1:nSamples){
    for(n in 1:nSampleRows[s]){
      yoyAlive[n,s]<-equals(z[sampleRows[n,s]],2)
      adultAlive[n,s]<-equals(z[sampleRows[n,s]],3)
    }
    nYoy[s]<-sum(yoyAlive[1:nSampleRows[s],s])
    nAdult[s]<-sum(adultAlive[1:nSampleRows[s],s])
  }

  }",file="test.txt")



inits<- function(){
  list(phiAdult=runif(1,0,1),
       phiYoy=runif(1,0,1),
       gamma=runif(winData$nSamples,0,1),
       pAdult=runif(1,0,1),
       pYoy=runif(1,0,1)
  )      
}



# MCMC settings
nb <- 1000
ni <- 3000
nt <- 2
nc <- 3

varsToMonitor<-c(
  'phiYoy'
  ,'phiAdult'
  ,'pYoy'
  ,'pAdult'
  ,'nYoy'
  ,'nAdult'
  ,'psi'
  ,'gamma'
)


  out <- jags(
    data=winData,
    inits=inits,
    model = "test.txt",
    parameters.to.save = varsToMonitor,
    n.chains=nc,
    n.iter = ni,
    n.thin = nt,
    n.burnin=nb)
sims<-out$BUGSoutput$sims.list
  
getSummary<-function(parameter){
  mean<-apply(sims[[parameter]],2,mean)
  quants<-t(apply(sims[[parameter]],2,quantile,probs=c(0.025,0.5,0.975)))
  
  result<-data.table(cbind(mean,quants))
  setnames(result,c("mean","lower","median","upper"))
  result[,sample:=0:20]
  return(result)
}

nYoy<-getSummary('nYoy')
nAdult<-getSummary('nAdult')

plot(mean~sample,data=nYoy,ylim=c(0,max(nYoy)))
points(median~sample,data=nYoy,pch=19)
with(nYoy,error.bar(sample,mean,upper,lower,interval.type='asdf'))

plot(mean~sample,data=nAdult,ylim=c(0,max(nAdult)))
points(median~sample,data=nAdult,pch=19)
with(nAdult,error.bar(sample,mean,upper,lower,interval.type='asdf'))
