##
 # FILE R_readMetData
 # AUTHOR olivier crespo
 # LANGUAGE R
 ########################################################################

## all format get out as
#list(	"file"=path to the file
#	"station"=list("id","lat","lon","alt","comm")
#	"clim"=list("tav","amp","refht","wndht")
#	"period"=list("start","end","type")	# type is CSAG style {0-real,1-365,2-360 days per year}
#	"data"=list("date","yyyy","mm","dd","juld","tmin","tmax","rain","srad","wind","dewp","vprs","rhum")
# MISSING VALUES = NA

########################################################################
## create empty forma
createNULLlist <- function()
{
	file <- 	NULL;
	station <- 	list("id"=NULL,"lat"=NULL,"lon"=NULL,"alt"=NULL,"comm"=NULL)
	clim <- 	list("tav"=NULL,"amp"=NULL,"refht"=NULL,"wndht"=NULL)
	period <- 	list("start"=NULL,"end"=NULL,"type"=NULL)
	data <- 	list("date"=NULL,"yyyy"=NULL,"mm"=NULL,"dd"=NULL,"juld"=NULL,"tmin"=NULL,"tmax"=NULL,"rain"=NULL,"srad"=NULL,"wind"=NULL,"dewp"=NULL,"vprs"=NULL,"rhum"=NULL)

	l <- 	list("file"=file,"station"=station,"clim"=clim,"period"=period,"data"=data)
return(l)
}

########################################################################
## replace any kind of missing value by NA
replace_missing <- function(data)
{
	for (d in 1:length(data)){
		data[[d]][data[[d]][]<(-90)] <- NA 
		data[[d]][is.na(data[[d]][])] <- NA
# add whatever you think of
	}

return(data)
}

########################################################################
##### NEW CSAG get tmax, tmin, ppt folder names in UPPER CASE
########################################################################
#inFo<- '/home/olivier/Desktop/Wine-shared/Projects/2013-2014_FFC/EasternCape/Climate/CSAGn'
#fNa<- '0078227_9.txt'
#fNa<- '0078227A3.txt'
#fNa<- '0078617_A.txt'
# inFolder is the direct parent folder of !! CAPITALES !! TMIN, TMAX, PPT
# fName is the station file name
# e.g. the tmax file is here: inFolder/tmax/fName
read_newCSAGformat <- function(inFolder=inFo,fName=fNa)
{
	cVar <- list.files(inFolder)	# very likely tmin, tmax, ppt
	head <- array(NA,dim=c(length(cVar),8))
	colnames(head) <- c('ID','LAT','LON','ALT','sDate','eDate','TYPE','COMM')
	rownames(head) <- cVar

	csag <- createNULLlist()
	csag$file <- 	paste(inFolder,'{TMIN,TMAX,PPT,RAD}',fName,sep="/")

	# read header
	inconsistences <- FALSE
	for (v in 1:length(cVar)){
		inFile<- paste(inFolder,cVar[v],fName,sep="/")
		line <- readLines(inFile,n=66)
		line<-	strsplit(line,split=" ")
		# rm blank cells
		for(l in 1:length(line)){
			line[[l]] <- line[[l]][line[[l]][]!=""]
		}
		# fill in head
		head[v,1] <- line[[grep("ID",line)[length(grep("ID",line))]]][3]	# ID
		head[v,2] <- line[[grep("LATITUDE",line)[length(grep("LATITUDE",line))]]][3]	# LAT
		head[v,3] <- line[[grep("LONGITUDE",line)[length(grep("LONGITUDE",line))]]][3]	# LON
		head[v,4] <- line[[grep("ALTITUDE",line)[length(grep("ALTITUDE",line))]]][3]	# ALT
		head[v,5] <- as.character(as.Date(line[[grep("START_DATE",line)[length(grep("START_DATE",line))]]][3],"%Y%m%d"))	# sDate
		head[v,6] <- as.character(as.Date(line[[grep("END_DATE",line)[length(grep("END_DATE",line))]]][3],"%Y%m%d"))	# eDate
## WATCH OUT !!!!
		head[v,7] <- 0	# TYPE
		head[v,8] <- line[[grep("NAME",line)[length(grep("NAME",line))]]][3]	# COMM
	
		# data
		if(cVar[v]=="PPT")	rain <- read.csv(inFile,skip=67,as.is=TRUE,nrows=-1)[,4]
		if(cVar[v]=="TMAX")	tmax <- read.csv(inFile,skip=67,as.is=TRUE,nrows=-1)[,4]
		if(cVar[v]=="TMIN")	tmin <- read.csv(inFile,skip=67,as.is=TRUE,nrows=-1)[,4]

		if(v>1){
			if(!all(head[v,]==head[1,])) inconsistences <- TRUE
		}
	}
	if(inconsistences){
		print("### > WARNING")
		print("#     there is inconsistencies in the file headers")
		print("#     we will assume ppt (id,lat,lon,alt,type) header")
		print("#     and take the longest possible period (filling with NA)")
		print(head[,1:7])
#		browser()
	}

	# make the local format
	csag$station$id <- 	head[rownames(head)=="PPT",colnames(head)=='ID']
	csag$station$lat <- 	as.numeric(head[rownames(head)=="PPT",colnames(head)=='LAT'])
	csag$station$lon <- 	as.numeric(head[rownames(head)=="PPT",colnames(head)=='LON'])
	csag$station$alt <- 	as.numeric(head[rownames(head)=="PPT",colnames(head)=='ALT'])
	csag$station$comm <- 	head[rownames(head)=="PPT",colnames(head)=='COMM']

	csag$clim$tav <- 	NA
	csag$clim$amp <- 	NA
	csag$clim$refht <- 	2	# WMO standard, may be different in some cases...
	csag$clim$wndht <- 	NA

#	csag$period$start <- 	as.Date(max(head[,5]))
#	csag$period$end <-	as.Date(min(head[,6]))
	csag$period$start <- 	as.Date(min(head[,5]))
	csag$period$end <-	as.Date(max(head[,6]))
	csag$period$type <-	as.numeric(head[rownames(head)=="PPT",7])

	csag$data$date <- 	seq(csag$period$start,csag$period$end,1)
	csag$data$yyyy <- 	as.numeric(format(csag$data$date,"%Y"))
	csag$data$mm <- 	as.numeric(format(csag$data$date,"%m"))
	csag$data$dd <- 	as.numeric(format(csag$data$date,"%d"))
	csag$data$juld <- 	as.numeric(format(csag$data$date,"%j"))
	csag$data$date <- 	as.numeric(format(csag$data$date,"%Y%j"))

	csag$data$tmin <- 	array(NA,dim=length(csag$data$date))
	csag$data$tmax <- 	array(NA,dim=length(csag$data$date))
	csag$data$rain <- 	array(NA,dim=length(csag$data$date))

	csag$data$srad <- 	array(NA,dim=length(csag$data$date))
	csag$data$wind <- 	array(NA,dim=length(csag$data$date))
	csag$data$dewp <- 	array(NA,dim=length(csag$data$date))
	csag$data$vprs <- 	array(NA,dim=length(csag$data$date))
	csag$data$rhum <- 	array(NA,dim=length(csag$data$date))


#browser()
	# longest possible period
	for (v in 1:length(cVar)){
		#start
		startD <- as.numeric(difftime(as.Date(head[rownames(head)==cVar[v],5]),csag$period$start,units="days"))
#		if (startD>0){
			if(cVar[v]=="PPT")	csag$data$rain[(startD+1):(length(rain)+startD)] <- rain
			if(cVar[v]=="TMAX")	csag$data$tmax[(startD+1):(length(tmax)+startD)] <- tmax
			if(cVar[v]=="TMIN")	csag$data$tmin[(startD+1):(length(tmin)+startD)] <- tmin
#		}
		#end
#		endD <- difftime(csag$period$end,as.Date(head[rownames(head)==cVar[v],6]),units="days")
#		if (endD>0){
#			if(cVar[v]=="ppt")			csag$data$rain <- csag$data$rain[1:(length(csag$data$rain)-endD)]
#			if(any(cVar[v]=="tmax",cVar[v]=="tmx"))	csag$data$tmax <- csag$data$tmax[1:(length(csag$data$tmax)-endD)]
#			if(any(cVar[v]=="tmin",cVar[v]=="tmn"))	csag$data$tmin <- csag$data$tmin[1:(length(csag$data$tmin)-endD)]
#		}
	}

	# different types (real,365,360)
## fill in if necessary

	# possible different data length
	if (length(csag$data$tmin)!=length(csag$data$date)){
		print(paste("tmin data length",length(csag$data$tmin),sep=" : "))
		print(paste("date data length",length(csag$data$date),sep=" : "))
		browser()
	}
	if (length(csag$data$tmax)!=length(csag$data$date)){
		print(paste("tmax data length",length(csag$data$tmax),sep=" : "))
		print(paste("date data length",length(csag$data$date),sep=" : "))
		browser()
	}
	if (length(csag$data$rain)!=length(csag$data$date)){
		print(paste("rain data length",length(csag$data$rain),sep=" : "))
		print(paste("date data length",length(csag$data$date),sep=" : "))
		browser()
	}

	# missing values
	csag$data <- replace_missing(csag$data)

return(csag)
rm(cVar,head,line,d,o,jD,startD,endD)
}


########################################################################
##### AGMIP
########################################################################
# input:
# path to the file to read in
#
# output:
# list filled up acording to empty list above
#########################################################################
#inFi<-'/home/olivier/Desktop/Wine-shared/Projects/2013-2014_FFC/EasternCape/Climate/MERRA/ZAFB0QXX.AgMIP'
#inFi<-'/home/crespo/Desktop/12_AgMIP/2012-10-01_fastTrack/AMIP/MerraData_CM/SABA0QXX.AgMIP'
#inFi<-'/home/crespo/Desktop/Link\ to\ WinShared/12_AgMIP/2012-10-01_fastTrack/AMIP/MerraData_CM/SABA0QXX.AgMIP'
## AgMIP
read_AgMIPformat <- function(inFile=inFi)
{
	sLine <- readLines(inFile)
	sData <- grep("DATE",sLine)
	sHead <- grep("INSI",sLine)
	rm(sLine)

	# read AgMIP data
	data <- scan(inFile,what='raw',skip=sData,quiet=T)
	data <- matrix(data,ncol=(dim(read.table(inFile,skip=sData,nrows=1))[2]),byrow=T)
	colnames(data) <- scan(inFile,what='raw',skip=(sData-1),nlines=1,quiet=T)[1:dim(data)[2]]
	# read header
	head <- read.table(inFile,skip=sHead,nrows=1)

	tmp <- scan(inFile,what='raw',skip=(sHead-1),nlines=1,quiet=T)	# because there is a space in between @ and INSI ...
	colnames(head) <- tmp[2:length(tmp)]
	rm(tmp)

	agmip <- createNULLlist()
	agmip$file <- 	inFile

	agmip$station$id <- 	as.character(head[1,1])
	agmip$station$lat <- 	as.numeric(head[colnames(head)=="LAT"])
	agmip$station$lon <- 	as.numeric(head[colnames(head)=="LONG"])
	agmip$station$alt <- 	as.numeric(head[colnames(head)=="ELEV"])
	agmip$station$comm <- 	scan(inFile,what='character',sep="\n",nlines=1,quiet=T)
	agmip$station$comm <- 	strsplit(agmip$station$comm,split=":")[[1]][length(strsplit(agmip$station$comm,split=":")[[1]])]

	agmip$clim$tav <- 	as.numeric(head[colnames(head)=="TAV"])
	agmip$clim$amp <- 	as.numeric(head[colnames(head)=="AMP"])
	agmip$clim$refht <- 	as.numeric(head[colnames(head)=="REFHT"])
	agmip$clim$wndht <- 	as.numeric(head[colnames(head)=="WNDHT"])

	agmip$period$start <-	as.Date(as.character(data[1,colnames(data)=="@DATE"]),"%Y%j")
	agmip$period$end <- 	as.Date(as.character(data[dim(data)[1],colnames(data)=="@DATE"]),"%Y%j")
	agmip$period$type <-	0

	agmip$data$date <- 	seq(agmip$period$start,agmip$period$end,1)
	agmip$data$yyyy <- 	as.numeric(format(agmip$data$date,"%Y"))
	agmip$data$mm <- 	as.numeric(format(agmip$data$date,"%m"))
	agmip$data$dd <- 	as.numeric(format(agmip$data$date,"%d"))
	agmip$data$juld <- 	as.numeric(format(agmip$data$date,"%j"))
	agmip$data$date <- 	as.numeric(format(agmip$data$date,"%Y%j"))

	agmip$data$srad <- 	as.numeric(data[,colnames(data)=="SRAD"])
	agmip$data$tmax <- 	as.numeric(data[,colnames(data)=="TMAX"])
	agmip$data$tmin <- 	as.numeric(data[,colnames(data)=="TMIN"])
	agmip$data$rain <- 	as.numeric(data[,colnames(data)=="RAIN"])
	if(any(colnames(data)=="WIND"))	{	agmip$data$wind <- 	as.numeric(data[,colnames(data)=="WIND"])
	}else{					agmip$data$wind <- 	NULL
	}
	if(any(colnames(data)=="DEWP"))	{	agmip$data$dewp <- 	as.numeric(data[,colnames(data)=="DEWP"])
	}else{					agmip$data$dewp <- 	NULL
	}
	if(any(colnames(data)=="VPRS"))	{	agmip$data$vprs <- 	as.numeric(data[,colnames(data)=="VPRS"])
	}else{					agmip$data$vprs <- 	NULL
	}
	if(any(colnames(data)=="RHUM"))	{	agmip$data$rhum <- 	as.numeric(data[,colnames(data)=="RHUM"])
	}else{					agmip$data$rhum <- 	NULL
	}

	# missing values
	agmip$data <- replace_missing(agmip$data)

return(agmip)
rm(data,head)
}

########################################################################
##### DSSAT
########################################################################
# input:
# path to the file to read in
#
# output:
# list filled up acording to empty list above
#########################################################################
#inFi<-'~/Desktop/MZ_00.WTH'
read_DSSATformat <- function(inFile=inFi,hist=TRUE)	# DSSAT does says 50 for 1950 and 2050, if hist=TRUE, then it will start with 19xx 
{
	sLine <- readLines(inFile)
	sData <- grep("@DATE",sLine)
	sHead <- grep("INSI",sLine)	# usually "@INSI", wiltrud "@ INSI"
	rm(sLine)
	# read DSSAT format data
	data <-	scan(inFile,what='raw',skip=sData)
	data <- matrix(data,ncol=(dim(read.table(inFile,skip=sData,nrows=1))[2]),byrow=T)
	colnames(data) <- scan(inFile,what='raw',skip=(sData-1),nlines=1)
	# read header
	head <- read.table(inFile,skip=sHead,nrows=1)
	tmp <- scan(inFile,what='raw',skip=(sHead-1),nlines=1)
	colnames(head) <- tmp[(length(tmp)-7):length(tmp)]	# for 8 values
	dssat <- createNULLlist()
	dssat$file <- 	inFile

	dssat$station$id <- 	as.character(head[1,grep("INSI",colnames(head))])
	dssat$station$lat <- 	as.numeric(head[1,colnames(head)=="LAT"])
	dssat$station$lon <- 	as.numeric(head[1,colnames(head)=="LONG"])
	dssat$station$alt <- 	as.numeric(head[1,colnames(head)=="ELEV"])
	dssat$station$comm <- 	scan(inFile,what='character',sep="\n",nlines=1)

	dssat$clim$tav <- 	as.numeric(head[colnames(head)=="TAV"])
	dssat$clim$amp <- 	as.numeric(head[colnames(head)=="AMP"])
	dssat$clim$refht <- 	as.numeric(head[colnames(head)=="REFHT"])
	dssat$clim$wndht <- 	as.numeric(head[colnames(head)=="WNDHT"])

	xxStart <- as.numeric(substr(as.character(data[1,colnames(data)=="@DATE"]),1,2))
	xxEnd <- as.numeric(substr(as.character(data[dim(data)[1],colnames(data)=="@DATE"]),1,2))

	if((xxEnd-xxStart)<0){
		XXs <- "19"
		XXe <- "20"
	}else{
		if(hist)	XXs <- XXe <- "19"
		if(!hist)	XXs <- XXe <- "20"
	}

	XXdate <- paste(XXs,as.character(data[1,colnames(data)=="@DATE"]),sep="")
	dssat$period$start <-	as.Date(XXdate,"%Y%j")
	XXdate <- paste(XXe,as.character(data[dim(data)[1],colnames(data)=="@DATE"]),sep="")
	dssat$period$end <- 	as.Date(XXdate,"%Y%j")
	if(dssat$period$end<dssat$period$star)	 dssat$period$end <- as.Date(paste((as.numeric(format(dssat$period$end,"%Y"))+100),format(dssat$period$end,"%m"),format(dssat$period$end,"%d"),sep="-"))
	dssat$period$type <-	0

	dssat$data$date <- 	seq(dssat$period$start,dssat$period$end,1)
	dssat$data$yyyy <- 	as.numeric(format(dssat$data$date,"%Y"))
	dssat$data$mm <- 	as.numeric(format(dssat$data$date,"%m"))
	dssat$data$dd <- 	as.numeric(format(dssat$data$date,"%d"))
	dssat$data$juld <- 	as.numeric(format(dssat$data$date,"%j"))
	dssat$data$date <- 	as.numeric(format(dssat$data$date,"%Y%j"))

	dssat$data$srad <- 	as.numeric(data[,colnames(data)=="SRAD"])
	dssat$data$tmax <- 	as.numeric(data[,colnames(data)=="TMAX"])
	dssat$data$tmin <- 	as.numeric(data[,colnames(data)=="TMIN"])
	dssat$data$rain <- 	as.numeric(data[,colnames(data)=="RAIN"])
	if(any(colnames(data)=="WIND"))	{	dssat$data$wind <- 	as.numeric(data[,colnames(data)=="WIND"])
	}else{						dssat$data$wind <- 	array(NA,dim=length(dssat$data$date))
	}
	if(any(colnames(data)=="DEWP"))	{	dssat$data$dewp <- 	as.numeric(data[,colnames(data)=="DEWP"])
	}else{						dssat$data$dewp <- 	array(NA,dim=length(dssat$data$date))
	}
	if(any(colnames(data)=="VPRS"))	{	dssat$data$vprs <- 	as.numeric(data[,colnames(data)=="VPRS"])
	}else{						dssat$data$vprs <- 	array(NA,dim=length(dssat$data$date))
	}
	if(any(colnames(data)=="RHUM"))	{	dssat$data$rhum <- 	as.numeric(data[,colnames(data)=="RHUM"])
	}else{						dssat$data$rhum <- 	array(NA,dim=length(dssat$data$date))
	}

	# missing values
	dssat$data <- replace_missing(dssat$data)

return(dssat)
rm(sData,sHead,data,head)
}

########################################################################
##### old CSAG
########################################################################
#inFo<- '/home/crespo/Desktop/12_AgMIP/2012-10-01_fastTrack/AMIP/E_RCP85/BNU-ESM-rcp85'
#fNa<- '0261516.1.txt'
#fNa<- '0261516B.txt'
## old CSAG
# inFolder is the direct parent folder of tmin,tmax,ppt
# fName is the station file name
# e.g. the tmax file is here: inFolder/tmax/fName
read_oldCSAGformat <- function(inFolder=inFo,fName=fNa)
{
	cVar <- list.files(inFolder)	# very likely tmin, tmax, ppt
	head <- array(NA,dim=c(length(cVar),8))
	colnames(head) <- c('ID','LAT','LON','ALT','sDate','eDate','TYPE','COMM')
	rownames(head) <- cVar

	csag <- createNULLlist()
	csag$file <- 	paste(inFolder,'{tmin,tmax,ppt,rad}',fName,sep="/")

	# read all headers
	inconsistences <- FALSE
	for (v in 1:length(cVar)){
		inFile<- paste(inFolder,cVar[v],fName,sep="/")
		line <- readLines(inFile,n=3)
		line<-	strsplit(line,split=" ")
		# rm blank cells
		for(l in 1:length(line)){
			line[[l]] <- line[[l]][line[[l]][]!=""]
		}
		# fill in head
		head[v,1] <- line[[1]][1]	# ID
		head[v,2] <- line[[1]][2]	# LAT
		head[v,3] <- line[[1]][3]	# LON
		head[v,4] <- line[[1]][4]	# ALT
		head[v,5] <- as.character(as.Date(line[[2]][1],"%Y%m%d"))	# sDate
		head[v,6] <- as.character(as.Date(line[[2]][2],"%Y%m%d"))	# eDate
		head[v,7] <- line[[2]][3]	# TYPE
		head[v,8] <- line[[3]][1]	# COMM
	
		# data
		if(cVar[v]=="ppt")			rain <- as.array(as.numeric(scan(inFile,what="numeric",sep="\n",skip=3,nlines=-1,quiet=TRUE)));
		if(any(cVar[v]=="tmax",cVar[v]=="tmx"))	tmax <- as.array(as.numeric(scan(inFile,what="numeric",sep="\n",skip=3,nlines=-1,quiet=TRUE)));
		if(any(cVar[v]=="tmin",cVar[v]=="tmn"))	tmin <- as.array(as.numeric(scan(inFile,what="numeric",sep="\n",skip=3,nlines=-1,quiet=TRUE)));

		if(v>1){
			if(!all(head[v,]==head[1,])) inconsistences <- TRUE
		}
	}
	if(inconsistences){
		print("### > WARNING")
		print("#     there is inconsistencies in the file headers")
		print("#     we will assume ppt (id,lat,lon,alt,type) header")
		print("#     and take the longest possible period (filling with NA)")
		print(head[,1:7])
#		browser()
	}

	# make the local format
	csag$station$id <- 	head[rownames(head)=="ppt",colnames(head)=='ID']
	csag$station$lat <- 	as.numeric(head[rownames(head)=="ppt",colnames(head)=='LAT'])
	csag$station$lon <- 	as.numeric(head[rownames(head)=="ppt",colnames(head)=='LON'])
	csag$station$alt <- 	as.numeric(head[rownames(head)=="ppt",colnames(head)=='ALT'])
	csag$station$comm <- 	head[rownames(head)=="ppt",colnames(head)=='COMM']

	csag$clim$tav <- 	NA
	csag$clim$amp <- 	NA
	csag$clim$refht <- 	2	# WMO standard, may be different in some cases...
	csag$clim$wndht <- 	NA

#	csag$period$start <- 	as.Date(max(head[,5]))
#	csag$period$end <-	as.Date(min(head[,6]))
	csag$period$start <- 	as.Date(min(head[,5]))
	csag$period$end <-	as.Date(max(head[,6]))
	csag$period$type <-	as.numeric(head[rownames(head)=="ppt",7])

	csag$data$date <- 	seq(csag$period$start,csag$period$end,1)
	csag$data$yyyy <- 	as.numeric(format(csag$data$date,"%Y"))
	csag$data$mm <- 	as.numeric(format(csag$data$date,"%m"))
	csag$data$dd <- 	as.numeric(format(csag$data$date,"%d"))
	csag$data$juld <- 	as.numeric(format(csag$data$date,"%j"))
	csag$data$date <- 	as.numeric(format(csag$data$date,"%Y%j"))

	csag$data$tmin <- 	array(NA,dim=length(csag$data$date))
	csag$data$tmax <- 	array(NA,dim=length(csag$data$date))
	csag$data$rain <- 	array(NA,dim=length(csag$data$date))

	csag$data$srad <- 	array(NA,dim=length(csag$data$date))
	csag$data$wind <- 	array(NA,dim=length(csag$data$date))
	csag$data$dewp <- 	array(NA,dim=length(csag$data$date))
	csag$data$vprs <- 	array(NA,dim=length(csag$data$date))
	csag$data$rhum <- 	array(NA,dim=length(csag$data$date))

	# longest possible period
	for (v in 1:length(cVar)){
		#start
		startD <- difftime(as.Date(head[rownames(head)==cVar[v],5]),csag$period$start,units="days")
#		if (startD>0){
			if(cVar[v]=="ppt")			csag$data$rain[(startD+1):(length(rain)+startD)] <- rain
			if(cVar[v]=="tmax" || cVar[v]=="tmx")	csag$data$tmax[(startD+1):(length(tmax)+startD)] <- tmax
			if(cVar[v]=="tmin" || cVar[v]=="tmn")	csag$data$tmin[(startD+1):(length(tmin)+startD)] <- tmin
#		}
		#end
#		endD <- difftime(csag$period$end,as.Date(head[rownames(head)==cVar[v],6]),units="days")
#		if (endD>0){
#			if(cVar[v]=="ppt")			csag$data$rain <- csag$data$rain[1:(length(csag$data$rain)-endD)]
#			if(any(cVar[v]=="tmax",cVar[v]=="tmx"))	csag$data$tmax <- csag$data$tmax[1:(length(csag$data$tmax)-endD)]
#			if(any(cVar[v]=="tmin",cVar[v]=="tmn"))	csag$data$tmin <- csag$data$tmin[1:(length(csag$data$tmin)-endD)]
#		}
	}

	# different types (real,365,360)
#	source('metTransformations.r')
	if (csag$period$type==1){
		# add up NA on 29-FEB
		feb29 <- csag$data$date[csag$data$mm==2 & csag$data$dd==29]
		for (l in 1:length(feb29)){
			w <- which(csag$data$date==feb29[l])
			csag$data$tmin <- c(csag$data$tmin[1:(w-1)],NA,csag$data$tmin[w:length(csag$data$tmin)]) 
			csag$data$tmax <- c(csag$data$tmax[1:(w-1)],NA,csag$data$tmax[w:length(csag$data$tmax)]) 
			csag$data$rain <- c(csag$data$rain[1:(w-1)],NA,csag$data$rain[w:length(csag$data$rain)]) 
		}
		if(all(is.na(csag$data$tmin[(length(csag$data$date)+1):length(csag$data$tmin)])))	csag$data$tmin<-csag$data$tmin[1:length(csag$data$date)]
		if(all(is.na(csag$data$tmax[(length(csag$data$date)+1):length(csag$data$tmax)])))	csag$data$tmax<-csag$data$tmax[1:length(csag$data$date)]
		if(all(is.na(csag$data$rain[(length(csag$data$date)+1):length(csag$data$rain)])))	csag$data$rain<-csag$data$rain[1:length(csag$data$date)]
	}	

	if (csag$period$type==2){
		# on non-leap year add
		add_l <- NULL
		add_l <- c(add_l,csag$data$date[csag$data$mm==2 & csag$data$dd==29])
		add_l <- c(add_l,csag$data$date[csag$data$mm==4 & csag$data$dd==15])
		add_l <- c(add_l,csag$data$date[csag$data$mm==6 & csag$data$dd==15])
		add_l <- c(add_l,csag$data$date[csag$data$mm==8 & csag$data$dd==15])
		add_l <- c(add_l,csag$data$date[csag$data$mm==10 & csag$data$dd==15])
		add_l <- c(add_l,csag$data$date[csag$data$mm==12 & csag$data$dd==15])
		for (l in 1:length(add_l)){
			w <- which(csag$data$date==add_l[l])
			csag$data$tmin <- c(csag$data$tmin[1:(w-1)],NA,csag$data$tmin[w:length(csag$data$tmin)]) 
			csag$data$tmax <- c(csag$data$tmax[1:(w-1)],NA,csag$data$tmax[w:length(csag$data$tmax)]) 
			csag$data$rain <- c(csag$data$rain[1:(w-1)],NA,csag$data$rain[w:length(csag$data$rain)]) 
		}
		if(all(is.na(csag$data$tmin[(length(csag$data$date)+1):length(csag$data$tmin)])))	csag$data$tmin<-csag$data$tmin[1:length(csag$data$date)]
		if(all(is.na(csag$data$tmax[(length(csag$data$date)+1):length(csag$data$tmax)])))	csag$data$tmax<-csag$data$tmax[1:length(csag$data$date)]
		if(all(is.na(csag$data$rain[(length(csag$data$date)+1):length(csag$data$rain)])))	csag$data$rain<-csag$data$rain[1:length(csag$data$date)]

		print("please check that this work properly")
		browser()
	}	

	# possible different data length
	if (length(csag$data$tmin)!=length(csag$data$date)){
		print(paste("tmin data length",length(csag$data$tmin),sep=" : "))
		print(paste("date data length",length(csag$data$date),sep=" : "))
		browser()
	}
	if (length(csag$data$tmax)!=length(csag$data$date)){
		print(paste("tmax data length",length(csag$data$tmax),sep=" : "))
		print(paste("date data length",length(csag$data$date),sep=" : "))
		browser()
	}
	if (length(csag$data$rain)!=length(csag$data$date)){
		print(paste("rain data length",length(csag$data$rain),sep=" : "))
		print(paste("date data length",length(csag$data$date),sep=" : "))
		browser()
	}

	# missing values
	csag$data <- replace_missing(csag$data)

return(csag)
rm(cVar,head,line,d,o,jD,startD,endD)
}

########################################################################
##### APSIM
########################################################################
# inFi<-'~/Desktop/12_AgMIP/2012-10-01_fastTrack/CLIP/Bulawayo_51_2002.met'
# inFi<-'/home/crespo/Desktop/Link\ to\ WinShared/12_AgMIP/2012-10-01_fastTrack/AMIP/MerraData_CM/SABA0QXX.AgMIP'
## AgMIP
read_APSIMformat <- function(inFile=inFi)
{
	# read APSIM format data
	head <- readLines(inFile,n=50)
	dTop <- grep("MJ/m2",head)
	tmp <- strsplit(head[(dTop-1)],split=" ")
	cNames <- array(tmp[[1]][tmp[[1]]!=""],dim=length(tmp[[1]][tmp[[1]]!=""]))
	data <- scan(inFile,what='raw',skip=dTop)
	data <- matrix(data,ncol=length(cNames),byrow=T)
	colnames(data) <- cNames

	# read header
	head <- head[1:dTop]
	head<-strsplit(head,split=" ")
	for(i in 1:length(tmp)){
		head[[i]]<-head[[i]][tmp[[i]]!=""]
		if (length(head[[i]])==0) head[[i]] <- NA
	}
	rLat <- head[[grep("atitude",head)]]	# row with lat
	rLon <- ifelse(length(grep("ongitude",head))==0,NA,head[[grep("ongitude",head)]])	# row with lon (if exists, not mandatory in APSIM format)
	rAlt <- ifelse(length(grep("ltitude",head))==0,NA,head[[grep("ltitude",head)]])		# row with alt (if exists, not mandatory in APSIM format)
# apparently problems with getting lon and alt
	
	apsim <- createNULLlist()
	apsim$file <- inFile
	apsim$station$id <- 	data[1,colnames(data)=="site"]
	apsim$station$lat <- 	as.numeric(rLat[grep("\\.",rLat)])
	apsim$station$lon <- 	ifelse(is.na(rLon),NA,as.numeric(rLon[grep("\\.",rLon)]))
	apsim$station$alt <- 	ifelse(is.na(rAlt),NA,as.numeric(rAlt[grep("\\.",rAlt)]))
	apsim$station$comm <- 	NA

	tTav <- grep("tav",head)
	tAmp <- grep("amp",head)
	rTav <- tTav[which((tTav-tAmp)!=0)]
	rAmp <- tAmp[which((tTav-tAmp)!=0)]
	cTav <- array(head[[rTav]][head[[rTav]]!=""],dim=length(head[[rTav]][head[[rTav]]!=""]))
	cAmp <- array(head[[rAmp]][head[[rAmp]]!=""],dim=length(head[[rAmp]][head[[rAmp]]!=""]))
	apsim$clim$tav <- 	as.numeric(cTav[grep("\\.",cTav)])
	apsim$clim$amp <- 	as.numeric(cAmp[grep("\\.",cAmp)])
	apsim$clim$refht <- 	NA
	apsim$clim$wndht <- 	NA

	apsim$data$yyyy <- 	as.numeric(data[,colnames(data)=="year"])
	apsim$data$juld <- 	as.numeric(data[,colnames(data)=="day"])
	apsim$period$start <- 	as.Date(paste(apsim$data$yyyy[1],apsim$data$juld[1],sep="-"),"%Y-%j")
	apsim$period$end <-	as.Date(paste(apsim$data$yyyy[length(apsim$data$yyyy)],apsim$data$juld[length(apsim$data$juld)],sep="-"),"%Y-%j")
	apsim$period$type <-	0

	apsim$data$date <- 	seq(apsim$period$start,apsim$period$end,1)
	apsim$data$mm <- 	as.numeric(format(apsim$data$date,"%m"))
	apsim$data$dd <- 	as.numeric(format(apsim$data$date,"%d"))
	apsim$data$date <- 	as.numeric(format(apsim$data$date,"%Y%j"))
	if(length(apsim$data$yyyy)!=length(apsim$data$date)){
		print(" ### > WARNING in read_APSIMformat",quote=F)
		print("     > data length issue",quote=F)
		stop("cannot resolve")
	}

	apsim$data$srad <- 	as.numeric(data[,colnames(data)=="radn"])
	apsim$data$tmax <- 	as.numeric(data[,colnames(data)=="maxt"])
	apsim$data$tmin <- 	as.numeric(data[,colnames(data)=="mint"])
	apsim$data$rain <- 	as.numeric(data[,colnames(data)=="rain"])
	apsim$data$wind <- 	array(NA,dim=length(apsim$data$date))
	apsim$data$dewp <- 	array(NA,dim=length(apsim$data$date))
	apsim$data$vprs <- 	array(NA,dim=length(apsim$data$date))
	apsim$data$rhum <- 	array(NA,dim=length(apsim$data$date))

	# missing values
	for (d in 1:length(apsim$data)){
		apsim$data[[d]][apsim$data[[d]][]<(-99)] <- NA 
		apsim$data[[d]][apsim$data[[d]][]==NaN] <- NA
	}

return(apsim)
rm(data,head,dTop,tmp,cNames,tLat,rLon,rAlt,rTav,rAmp,tTav,tAmp,cTav,cAmp)
}

########################################################################
##### NCDF/ASCII
########################################################################
# actually not ncdf,
# I used CDO to subselect and reformat to ascii
# so I'm reading ONE column of this ascii file 
########################################################################
#inFo <- '/home/crespo/Desktop/12_AgMIP/2012-10-01_fastTrack/CLIP/D_RCP45_raw/ASCII'
#gc <- 'bnu-esm'
# inFolder is the direct parent folder of tmin,tmax,ppt
#la <- -18.75
#lo <- 28.75
#al <- 1200
#co <- 5
#id <- 'BNU-45'
#com <-'GCM watch grid -18.5,-19 28.5,29'
#sDa <- as.Date('1960-01-01')
#eDa <- as.Date('2100-12-31')
read_ASCIIformat <- function(inFolder=inFo,gcm=gc,col=co,lat=la,lon=lo,alt=al,sDate=sDa,eDate=eDa,iden=id,comm=com)
{
print(gcm)
	cVar <- list.files(inFolder)	# very likely tmin, tmax, ppt
	ascii <- createNULLlist()
	ascii$file <- 	paste(inFolder,'{tmin,tmax,ppt}',gcm,sep="/")
	# read all headers
	for (v in 1:length(cVar)){
		fiName <- list.files(paste(inFolder,cVar[v],sep="/"))		
		inFile<- paste(inFolder,cVar[v],fiName[grep(fiName,pattern=gcm)],sep="/")
print(inFile)
		if(cVar[v]=="ppt")	ascii$data$rain <- read.table(inFile)[,col]
		if(cVar[v]=="tmax")	ascii$data$tmax <- read.table(inFile)[,col]-273.15	# as far as I know it is in absolute degress celsius
		if(cVar[v]=="tmin")	ascii$data$tmin <- read.table(inFile)[,col]-273.15
	}

	# make the local format
	ascii$station$id <- 	iden
	ascii$station$lat <- 	lat
	ascii$station$lon <- 	lon
	ascii$station$alt <- 	alt
	ascii$station$comm <- 	comm

	ascii$clim$tav <- 	NA
	ascii$clim$amp <- 	NA
	ascii$clim$refht <- 	NA	# 2: WMO standard, may be different in some cases...
	ascii$clim$wndht <- 	NA

	ascii$period$start <- 	sDate
	ascii$period$end <-	eDate

	ascii$data$date <- 	seq(ascii$period$start,ascii$period$end,1)
	ascii$data$yyyy <- 	as.numeric(format(ascii$data$date,"%Y"))
	ascii$data$mm <- 	as.numeric(format(ascii$data$date,"%m"))
	ascii$data$dd <- 	as.numeric(format(ascii$data$date,"%d"))
	ascii$data$juld <- 	as.numeric(format(ascii$data$date,"%j"))
	ascii$data$date <- 	as.numeric(format(ascii$data$date,"%Y%j"))
	if(length(ascii$data$date)==length(ascii$data$rain)){	ascii$period$type <- 0
	}else{
		if(floor((ascii$data$yyyy[length(ascii$data$yyyy)]-ascii$data$yyyy[1]+1)/4)==(length(ascii$data$date)-length(ascii$data$rain))){
			ascii$period$type <- 1
		}else{
			print('### WARNING')
			print('### there is a problem here : dataRead:read_ASCIIformat:lengths')
			print(paste('length difference',length(ascii$data$date)-length(ascii$data$rain),sep=" > "))
			print(paste('(assumed) number of leap year within',floor((ascii$data$yyyy[length(ascii$data$yyyy)]-ascii$data$yyyy[1]+1)/4),sep=" > "))		
			print(paste('unexplained',(length(ascii$data$date)-length(ascii$data$rain))-(floor((ascii$data$yyyy[length(ascii$data$yyyy)]-ascii$data$yyyy[1]+1)/4)),sep=" > "))
#			browser()
# for now
		}
	}

	ascii$data$srad <- 	array(NA,dim=length(ascii$data$date))
	ascii$data$wind <- 	array(NA,dim=length(ascii$data$date))
	ascii$data$dewp <- 	array(NA,dim=length(ascii$data$date))
	ascii$data$vprs <- 	array(NA,dim=length(ascii$data$date))
	ascii$data$rhum <- 	array(NA,dim=length(ascii$data$date))

	# longest common period
	# does not deal with that, maybe should	

	# different types (real,365,360)
#	source('metTransformations.r')
	if (ascii$period$type==1){
		# add up NA on 29-FEB
		feb29 <- ascii$data$date[ascii$data$mm==2 & ascii$data$dd==29]
		for (l in 1:length(feb29)){
			w <- which(ascii$data$date==feb29[l])
			ascii$data$tmin <- c(ascii$data$tmin[1:(w-1)],NA,ascii$data$tmin[w:length(ascii$data$tmin)]) 
			ascii$data$tmax <- c(ascii$data$tmax[1:(w-1)],NA,ascii$data$tmax[w:length(ascii$data$tmax)]) 
			ascii$data$rain <- c(ascii$data$rain[1:(w-1)],NA,ascii$data$rain[w:length(ascii$data$rain)]) 
		}
	}	
	if (ascii$period$type==2){
		print("#### > WARNING : do not deal with type==2 yet")
		browser();
	}	

	# possible different data length
	if (length(ascii$data$tmin)!=length(ascii$data$date)){
		print(paste("tmin data length",length(ascii$data$tmin),sep=" : "))
		print(paste("date data length",length(ascii$data$date),sep=" : "))
		browser()
	}
	if (length(ascii$data$tmax)!=length(ascii$data$date)){
		print(paste("tmax data length",length(ascii$data$tmax),sep=" : "))
		print(paste("date data length",length(ascii$data$date),sep=" : "))
		browser()
	}
	if (length(ascii$data$rain)!=length(ascii$data$date)){
		print(paste("rain data length",length(ascii$data$rain),sep=" : "))
		print(paste("date data length",length(ascii$data$date),sep=" : "))
		browser()
	}

	# missing values
	for (d in 1:length(ascii$data)){
		ascii$data[[d]][ascii$data[[d]][]<(-99999)] <- NA 	#NAM data
		ascii$data[[d]][ascii$data[[d]][]<(-99)] <- NA
		ascii$data[[d]][ascii$data[[d]][]==NaN] <- NA 
	}

return(ascii)
rm(cVar,ascii,fiName,inFile,v,feb29,l,d)
}

######################################################################################
######################################################################################
### template for a 30-31 days per line
### down year-months
read_YMlines <- function(newD,rawD,var,init=FALSE)
{

#	rawD<-read.csv("~/Desktop/wine_shared/12_AgMIP/LES-AMIIP/Baseline/mejametalana1.csv")
#	newD <- createNULLlist()
#	newD$station$ ...
#	newD$period$ ...

	if (init){
		newD$data$date <- 	seq(newD$period$start,newD$period$end,1)
		newD$data$yyyy <- 	as.numeric(format(newD$data$date,"%Y"))
		newD$data$mm <- 	as.numeric(format(newD$data$date,"%m"))
		newD$data$dd <- 	as.numeric(format(newD$data$date,"%d"))
		newD$data$juld <- 	as.numeric(format(newD$data$date,"%j"))
		newD$data$date <- 	as.numeric(format(newD$data$date,"%Y%j"))
		newD$data$tmin <- 	array(NA,dim=length(newD$data$date))
		newD$data$tmax <- 	array(NA,dim=length(newD$data$date))
		newD$data$rain <- 	array(NA,dim=length(newD$data$date))
	}

	#first column/line/day
	day_1to31 <- 4:34
	yearmm_c <- 3
	var_c <- 2 	

	cur_r <-1
	while((cur_r<=dim(rawD)[1])&&(rawD[cur_r,var_c]!=var))	{cur_r <- cur_r+1}
	while(!is.na(rawD[cur_r,yearmm_c])){
		ym <- as.Date(paste(as.character(rawD[cur_r,yearmm_c]),"01",sep="-"))
		dayInMonth<-maxNo_days(as.numeric(format(ym,"%Y")),as.numeric(format(ym,"%m")))
		day<-0
		for(col in day_1to31){
			# depending on variable {tmin,tmax,rain} data (the rest)
			year <- as.numeric(format(ym,"%Y"))
			month <- as.numeric(format(ym,"%m"))
			day <- day+1
			if (day>dayInMonth) {break}
			switch(var,{	#1
				},{	#2:tmax
					newD$data$tmax[(newD$data$yyyy==year)&(newD$data$mm==month)&(newD$data$dd==day)] <- rawD[cur_r,col]
				},{	#3:tmin
					newD$data$tmin[(newD$data$yyyy==year)&(newD$data$mm==month)&(newD$data$dd==day)] <- rawD[cur_r,col]
				},{	#4
				},{	#5:rain
					newD$data$rain[(newD$data$yyyy==year)&(newD$data$mm==month)&(newD$data$dd==day)] <- rawD[cur_r,col]
				}
			)
		}
		cur_r <- cur_r+1
		while((cur_r<=dim(rawD)[1])&&(rawD[cur_r,var_c]!=var))	{cur_r <- cur_r+1}
	}
	# missing values
	for (d in 1:length(newD$data)){
		newD$data[[d]][newD$data[[d]][]<(-90)] <- NA 
		newD$data[[d]][is.na(newD$data[[d]][])] <- NA
	}
	

return(newD)
rm(day_1to31,yearmm_c,var_c,cur_r,ym,day,year,month)
}

######################################################################################
######################################################################################
### template for a 1 day per line
### down year-month-day
read_YMDlines <- function(newD,rawD,var,init=FALSE)
{

#	rawD<-read.csv("~/Desktop/wine_shared/12_AgMIP/LES-AMIIP/Baseline/mejametalana1.csv")
#	newD <- createNULLlist()
#	newD$station$ ...
#	newD$period$ ...

	if (init){
		newD$data$date <- 	seq(newD$period$start,newD$period$end,1)
		newD$data$yyyy <- 	as.numeric(format(newD$data$date,"%Y"))
		newD$data$mm <- 	as.numeric(format(newD$data$date,"%m"))
		newD$data$dd <- 	as.numeric(format(newD$data$date,"%d"))
		newD$data$juld <- 	as.numeric(format(newD$data$date,"%j"))
		newD$data$date <- 	as.numeric(format(newD$data$date,"%Y%j"))
		newD$data$tmin <- 	array(NA,dim=length(newD$data$date))
		newD$data$tmax <- 	array(NA,dim=length(newD$data$date))
		newD$data$rain <- 	array(NA,dim=length(newD$data$date))
	}

	#init column/line/day
	day_c <- 7
	month_c <- 6
	year_c <- 5
	var_c <- 3
	val <- 8
	cur_r <-1
	while((cur_r<=dim(rawD)[1])&&(rawD[cur_r,var_c]!=var))	{cur_r <- cur_r+1}
	while(!is.na(rawD[cur_r,yearmm_c])){
		# depending on variable {tmin,tmax,rain} data (the rest)
		year <- rawD[cur_r,year_c]
		month <- rawD[cur_r,month_c]
		day <- rawD[cur_r,day_c]
		switch(var,{	#1
			},{	#2:tmax
				newD$data$tmax[(newD$data$yyyy==year)&(newD$data$mm==month)&(newD$data$dd==day)] <- rawD[cur_r,val]
			},{	#3:tmin
				newD$data$tmin[(newD$data$yyyy==year)&(newD$data$mm==month)&(newD$data$dd==day)] <- rawD[cur_r,val]
			},{	#4
			},{	#5:rain
				newD$data$rain[(newD$data$yyyy==year)&(newD$data$mm==month)&(newD$data$dd==day)] <- rawD[cur_r,val]
			}
		)
		cur_r <- cur_r+1
	}
	# missing values
	for (d in 1:length(newD$data)){
		newD$data[[d]][newD$data[[d]][]<(-90)] <- NA 
		newD$data[[d]][is.na(newD$data[[d]][])] <- NA
	}

return(newD)
rm(day_1to31,yearmm_c,var_c,cur_r,ym,day,year,month)
}
