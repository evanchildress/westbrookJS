bugsData<-list(
  #Long format data
  encDATA=d$enc
  ,river=d$riverN
  ,sample=d$sampleNumAdj
  ,z=d$zKnown
  
  #Control pieces
  ,nSamples=evalList$nSamples
  ,nRivers=evalList$nRivers
  ,evalRows=evalList$evalRows
  ,nEvalRows=evalList$nEvalRows
  ,summerSamples=evalList$summerSamples
  ,nonSummerSamples=evalList$nonSummerSamples
  #,firstObsRows=evalList$firstObsRows
  #,nFirstObsRows=evalList$nFirstObsRows
  )

inits<- function(){
  list(meanPhi=runif(1,0,1),
       gamma=array(runif(bugsData$nSamples*bugsData$nRivers,0,1),
                   dim=c(bugsData$nSamples,bugsData$nRivers)),
       meanP=runif(1,0,1)
  )      
}



# MCMC settings
na <- 500
nb <- 1000
ni <- 2000
nt <- 3
nc <- 3

varsToMonitor<-c(
  'meanPhi'
  ,'meanP'
  ,'psi'
  )

rm(dMData)
rm(evalList)
gc()

(beforeJags<-Sys.time())

out <- jags(
  data=bugsData,
  inits=inits,
  model = "bugsJS.txt",
  parameters.to.save = varsToMonitor,
  n.chains=nc,
  n.iter = ni,
  n.thin = nt,
  n.burnin=nb)

( done <- Sys.time() ) 
print(done - beforeJags)

