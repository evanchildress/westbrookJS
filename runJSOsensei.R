parallel=T

simNum <- 1
rDataName <- 'dMDataOutbkt1996_2014.RData'
dataStore<-""

##########
comment<-paste(" model",simNum, " MS Sim "  )

home <- file.path("~/westbrookJS")
setwd(home)

source("constructBugsJS.r")

directory <- tempfile( pattern="output-", tmpdir ='.', fileext='-JS')
dir.create(directory)

bugsName <- paste0('./bugsJS','.txt')

file.copy(from='./callJS.R', to=paste(directory,'callJS.R',sep='/'))
file.copy(from=bugsName , to=paste(directory,bugsName ,sep='/'))
file.copy(from='./run.R', to=paste(directory,'run.R',sep='/'))
##################

fileOutName <- "outJS.RData" 

#########################################
load(file.path(dataStore,rDataName))

source('./callJS.R')

#writeLines(text=paste(date(),directory,afterAdapt - beforeAdapt,done - beforeJags), con='../latest_directory')
writeLines(text=paste(date(),directory,afterAdapt - beforeAdapt,done - beforeJags,"[", comment,"]"), con='./info.txt')
getwd()

save(d, out, file = fileOutName)
save(out,file=file.path(directory,fileOutName))
     