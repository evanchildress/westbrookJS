load("~/westbrookJS/resultsBnt/outJS.rDATA")
outArray<-as.array(out)
outArray<-outArray[1001:2000,,]



summarize<-function(x,margin=c(2)){
  getSummary<-function(x){
  meanX<-mean(x)
  quantiles<-quantile(x,c(0.025,0.975))
  return(c(meanX,quantiles))
  }
  stringSplit<-function(x,split,element){
    sapply(strsplit(x,split), "[[", element)
  }
  
  result<-data.table(dcast(melt(apply(x,MARGIN=margin,FUN=getSummary)),
                           var~Var1))
  setnames(result,c("var","mean","lower","upper"))
  result[,var:=as.character(var)]
  
  result[,parameter:=unlist(strsplit(var,"[[]"))[seq(1,nrow(result)*2-1,2)]]
  result[grep(",",var),
         sample:=stringSplit(stringSplit(var,"[[]",2),",",1)]
  result[,river:=substr(var,start=nchar(var)-1,stop=nchar(var)-1)]
  result[,sample:=as.numeric(sample)]
  
  return(result)
}

results<-summarize(outArray)
setkey(results,parameter,river,sample)
results[,season:=(sample/4-floor(sample/4))*4]
results[season==0,season:=4]
results[,year:=ceiling(sample/4)+2000]

toSummarize<-c("pYoy","pAdult","phiYoy","phiAdult","nYoy","nAdult")

rivers<-c("jimmy","mitchell","west brook")

for(i in toSummarize){
  tiff.par(paste0("~/westbrookJS/resultsBnt/",i,".tif"),mfrow=c(3,1),mar=c(3,3,1,0),mgp=c(2,0.5,0))
  for(r in 1:3){
    plot(mean~sample,data=results[river==r & parameter==i],
         pch=19,col=palette()[season],
         ylim=if(grepl('n',i)) c(0,max(upper)) else c(0,1.1),
         main=rivers[r],xlab="",ylab=i)
    with(results[river==r & parameter==i],
         error.bar(sample,mean,
                   upper.y=upper,lower.y=lower,
                   interval.type='bla')) 
  }
  #legend('bottom',c("spring","summer","fall","winter"),palette()[1:4],pch=19)
  title(xlab="Sample Number")
  dev.off()
}

tiff.par("~/westbrookJS/resultsBnt/nYoyFall.tif",mfrow=c(3,1),mar=c(3,3,1,0))
for(r in 1:3){
  plot(mean~year,data=results[river==r & season==3 & parameter=="nYoy"],
       pch=19,type='b',col=palette()[r])
}
dev.off()

pAdult<-results[parameter=='pAdult',list(mean,lower,upper,
                                         sample,river,season,year)]
pAdult[river==3,river:="4"]
obear<-data.table(mean=1,lower=1,upper=1,
                  sample=pAdult[river==1,sample],
                  river=3,
                  season=pAdult[river==1,season],
                  year=pAdult[river==1,year])
pAdult<-rbind(pAdult,obear)

pYoy<-results[parameter=='pYoy',list(mean,lower,upper,
                                     sample,river,season,year)]
pYoy[river==3,river:="4"]
obear<-data.table(mean=1,lower=1,upper=1,
                  sample=pYoy[river==1,sample],
                  river=3,
                  season=pYoy[river==1,season],
                  year=pYoy[river==1,year])
pYoy<-rbind(pYoy,obear)

nAdult<-results[parameter=='nAdult',list(mean,lower,upper,
                                     sample,river,season,year)]
nAdult[river==3,river:="4"]
obear<-data.table(mean=0,lower=0,upper=0,
                  sample=nAdult[river==1,sample],
                  river=3,
                  season=nAdult[river==1,season],
                  year=nAdult[river==1,year])
nAdult<-rbind(nAdult,obear)

densityPlot<-function(x,nchains=3,lwd=2,...){
  plot(density(x[,1]),col=palette()[1],...)
  for(i in 2:nchains){
    points(density(x[,i]),col=palette()[i],
    type='l')
  }
}

for(i in 506:672){
  densityPlot(outArray[,i,],main=dimnames(outArray)[[2]][i])
}
                  
saveRDS(pAdult,"~/westbrookJS/resultsBnt/pAdult.rds")
saveRDS(pYoy,"~/westbrookJS/resultsBnt/pYoy.rds")
saveRDS(nAdult,"~/westbrookJS/resultsBnt/nAdult.rds")
