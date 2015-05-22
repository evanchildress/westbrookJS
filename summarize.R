load("~/westbrookJS/results/outJS.rDATA")
sims<-out$BUGSoutput$sims.list

summarize<-function(x,margin=c(2,3)){
  getSummary<-function(x){
  meanX<-mean(x)
  quantiles<-quantile(x,c(0.025,0.975))
  return(c(meanX,quantiles))
  }
  result<-data.table(dcast(melt(apply(x,MARGIN=margin,FUN=getSummary)),
                           Var2+Var3~Var1))
  setnames(result,c("sample","river","mean","lower","upper"))
  result[,season:=rep(rep(1:4,each=4),14)]
  return(result)
}

toSummarize<-c("pYoy","pAdult","phiYoy","phiAdult","nYoy","nAdult")
for(i in toSummarize){
  assign(i,summarize(sims[[i]]))
}

rivers<-c("jimmy","mitchell","obear","west brook")

for(i in toSummarize){
  tiff.par(paste0("results/",i,".tif"),mfrow=c(4,1),mar=c(3,3,1,0),mgp=c(2,0.5,0))
  for(r in 1:4){
    plot(mean~sample,data=get(i)[river==r],
         pch=19,col=palette()[season],
         ylim=if(grepl('n',i)) c(0,max(upper)) else c(0,1.1),
         main=rivers[r],xlab="",ylab=i)
    with(get(i)[river==r],
         error.bar(sample,mean,
                   upper.y=upper,lower.y=lower,
                   interval.type='bla')) 
  }
  #legend('bottom',c("spring","summer","fall","winter"),palette()[1:4],pch=19)
  title(xlab="Sample Number")
  dev.off()
}

pAdult[,year:=ceiling(sample/4)+2000]
pYoy[,year:=ceiling(sample/4)+2000]
nAdult[,year:=ceiling(sample/4)+2000]

saveRDS(pAdult,"~/westbrookJS/results/pAdult.rds")
saveRDS(pYoy,"~/westbrookJS/results/pYoy.rds")
saveRDS(nAdult,"~/westbrookJS/results/nAdult.rds")
