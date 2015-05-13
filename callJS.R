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
  ,nSummerSamples=evalList$nSummerSamples
  ,nonSummerSamples=evalList$nonSummerSamples
  ,nNonSummerSamples=evalList$nNonSummerSamples
  ,sampleRows=evalList$sampleRows
  ,nSampleRows=evalList$nSampleRows
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
na <- 1
nb <- 1
ni <- 3
nt <- 1
nc <- 3

varsToMonitor<-c(
   'phiYoy'
  ,'phiAdult'
  ,'pYoy'
  ,'pAdult'
  ,'b'
  ,'nYoy'
  ,'nAdult'
  ,'psi'
  )

rm(dMData)
rm(evalList)
gc()

(beforeJags<-Sys.time())
if(parallel==F){
out <- jags(
  data=bugsData,
  inits=inits,
  model = "bugsJS.txt",
  parameters.to.save = varsToMonitor,
  n.chains=nc,
  n.iter = ni,
  n.thin = nt,
  n.burnin=nb)
} else {
  coda.samples.wrapper <- function(j)
  { 
    temp.model = jags.model("bugsJS.txt", 
                            inits=inits, 
                            data=bugsData,
                            n.chains=nc,
                            n.adapt=na)
    coda.samples(temp.model, varsToMonitor, n.iter=ni, thin=nt) 
  }
  
  snow.start.time = proc.time()
  cl <- makeCluster(nc, "SOCK")
  ##Make sure the rjags library is loaded in each worker
  clusterEvalQ(cl, library(rjags))
  ##Send data to workers, then fit models. One disadvantage of this
  ##parallelization is that you lose the ability to watch the progress bar.
  clusterExport(cl, list('bugsData','ni','nt','nc','na','varsToMonitor'))
  par.samples = clusterApply(cl, 1:nc, coda.samples.wrapper)
  ##Reorganize 'par.samples' so that it is recognizeable as an 'mcmc.list' object
  for(i in 1:length(par.samples)) { par.samples[[i]] <- par.samples[[i]][[1]] }
  class(par.samples) <- "mcmc.list"
  stopCluster(cl)
  snow.end.time = proc.time()
  snow.dtime = snow.end.time - snow.start.time
}

( done <- Sys.time() ) 
print(done - beforeJags)

