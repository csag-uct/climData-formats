##
 # FILE dataWrite.r
 # AUTHOR olivier crespo
 # https://r-forge.r-project.org/projects/xpos-r/
 ###############################################################################

## all format get in as
#list(	"file"=path to the file
#	"station"=list("id","lat","lon","alt","comm")
#	"clim"=list("tav","amp","refht","wndht")
#	"period"=list("start","end","type")	# type is metD style {0-real,1-365,2-360 days per year}
#	"data"=list("date","yyyy","mm","dd","juld","srad","tmax","tmin","rain","wind","dewp","vprs","rhum")
# MISSING VALUES = NA

########################################################################
##### AGMIP
########################################################################
#outFiAgmip<-'/home/crespo/Desktop/12_AgMIP/2012-10-01_fastTrack/AMIP/WorkingDir/SABA0Qtest.AgMIP'
write_AgMIPformat <- function(metD,outFile=outFiAgmip,checkMiss=TRUE)
{
	
	# template order
	# @DATE    YYYY  MM  DD  SRAD  TMAX  TMIN  RAIN  WIND  DEWP  VPRS  RHUM
	table <- array(NA,dim=c(length(metD$data[[1]]),12))
	for (v in 1:dim(table)[2]){
		switch(v,
			table[,v] <- metD$data[[which(names(metD$data)=="date")]][],
			table[,v] <- metD$data[[which(names(metD$data)=="yyyy")]][],
			table[,v] <- metD$data[[which(names(metD$data)=="mm")]][],
			table[,v] <- metD$data[[which(names(metD$data)=="dd")]][],
			table[,v] <- metD$data[[which(names(metD$data)=="srad")]][],
			table[,v] <- metD$data[[which(names(metD$data)=="tmax")]][],
			table[,v] <- metD$data[[which(names(metD$data)=="tmin")]][],
			table[,v] <- metD$data[[which(names(metD$data)=="rain")]][],
			table[,v] <- metD$data[[which(names(metD$data)=="wind")]][],
			table[,v] <- metD$data[[which(names(metD$data)=="dewp")]][],
			table[,v] <- metD$data[[which(names(metD$data)=="vprs")]][],
			table[,v] <- metD$data[[which(names(metD$data)=="rhum")]][]
		)
	}

	# missing values
	# there should be no missing values at this stage
	if(checkMiss){
	for (d in 1:length(metD$data)){
		if(any(is.na(metD$data[[d]][]))){
			print("### > WARNING missing values",quote=F)
			print("    > Replacing by -99",quote=F)
		}
	}}
	for(r in 1:dim(table)[1]){
		for(c in 1:dim(table)[2]){
			if(is.na(table[r,c]))	table[r,c]<-(-99)
		}
	}

	# file name
	if(!file.copy("./agmipTemplates/agmipTemplate.AgMIP",outFile,overwrite=TRUE)){
		print("### > FAILED TO CREATE OUTPUT FILE",quote=F)
		stop("cannot resolve")
	}
	ifelse(Sys.info()["sysname"]=="Linux",newline <- "\r\n",newline<-"\n")

## INSPIRED FROM L.ESTES
  	output <- rbind(	sprintf("%s %s", "*WEATHER DATA :", metD$station$comm),	# First row of WTH file
				c(""),
			 	c('@ INSI      LAT     LONG  ELEV   TAV   AMP REFHT WNDHT'),	# Second tier
			 	sprintf("%6s %8.2f %8.2f %5.0f %5.1f %5.1f %5.2f %5.2f", metD$station$id, round(metD$station$lat,2),round(metD$station$lon,2),round(metD$station$alt,0),round(metD$clim$tav,1),round(metD$clim$amp,1),round(metD$clim$refht,2),round(metD$clim$wndht,2)),
				c('@DATE    YYYY  MM  DD  SRAD  TMAX  TMIN  RAIN  WIND  DEWP  VPRS  RHUM'),	# Data tier
				cbind(sprintf("%5i %5i %3i %3i %5.1f %5.1f %5.1f %5.1f %5i %5i %5i %5i", table[,1], table[,2], table[,3], table[,4], table[,5], table[,6], table[,7], table[,8], table[,9], table[,10], table[,11], table[,12]))
		)
	write(output,outFile)

### OLD
#	write.table(paste("*WEATHER DATA : ",metD$station$id,sep=''),outFile,append=FALSE,quote=FALSE,col.names=FALSE,row.names=FALSE,eol=newline);
#	write.table('',outFile,append=TRUE,quote=FALSE,col.names=FALSE,row.names=FALSE,eol=newline);
#	write.table(format.AsIs(t(c('@ INSI','LAT','LONG','ELEV','TAV','AMP','REFHT','WNDHT')),justify="right",width=4),outFile,append=TRUE,quote=FALSE,col.names=FALSE,row.names=FALSE,eol=newline)
#	h1 <- array(c(paste('  ',metD$station$id,sep=''),round(metD$station$lat,2),round(metD$station$lon,2),round(metD$station$alt,0),round(metD$clim$tav,1),round(metD$clim$amp,1),round(metD$clim$refht,2),round(metD$clim$wndht,2)),dim=8)
#	write.table(format.AsIs(t(h1),justify="right",width=4),outFile,append=TRUE,quote=FALSE,col.names=FALSE,row.names=FALSE,eol=newline)

##	h2 <- array(NA,dim=12)
##	h2[1] <- format.AsIs('@DATE',justify="right",width=8)
##	h2[2:12] <- format.AsIs(c('YYYY','MM','DD','SRAD','TMAX','TMIN','RAIN','WIND','DEWP','VPRS','RHUM'),justify="right",width=4)
##	write.table(t(h2),outFile,append=TRUE,quote=FALSE,col.names=FALSE,row.names=FALSE,eol=newline)

#	tmp <- array(NA,dim=c(dim(table)[1]+1,dim(table)[2]))
#	tmp[,1] <- format.AsIs(c('@DATE',table[,1]),justify="left")
#	tmp[,2] <- format.AsIs(c('YYYY',table[,2]),justify="right")
#	tmp[,3] <- format.AsIs(c('MM',table[,3]),justify="right")
#	tmp[,4] <- format.AsIs(c('DD',table[,4]),justify="right")
#	tmp[,5] <- format.AsIs(c('SRAD',table[,5]),justify="right")
#	tmp[,6] <- format.AsIs(c('TMAX',table[,6]),justify="right")
#	tmp[,7] <- format.AsIs(c('TMIN',table[,7]),justify="right")
#	tmp[,8] <- format.AsIs(c('RAIN',table[,8]),justify="right")
#	tmp[,9] <- format.AsIs(c('WIND',table[,9]),justify="right")
#	tmp[,10] <- format.AsIs(c('DEWP',table[,10]),justify="right")
#	tmp[,11] <- format.AsIs(c('VPRS',table[,11]),justify="right")
#	tmp[,12] <- format.AsIs(c('RHUM',table[,12]),justify="right")
#	write.table(table,outFile,append=TRUE,quote=FALSE,col.names=FALSE,row.names=FALSE,eol=newline,sep=" ")

#	tmp <- format.AsIs(t(table),justify="right",width=4)
#	write.table(tmp,outName,quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE);
#	write(t(tmp),outFile,ncolumns=dim(table)[2],append=TRUE)

	## to handle the linux/windows format issue
	if(Sys.info()["sysname"]=="Linux"){
		tempFile <- readLines(outFile,n=-1,warn=FALSE);
		writeLines(tempFile,outFile,sep="\r\n");
	}

rm(table,output)
}

########################################################################
##### oldCSAG
########################################################################
#outFometD<-'/home/crespo/Desktop/TMP/CSAG'
write_oldCSAGformat <- function(metD,outFolder=outFometD)
{
	# if you built a full set, then switch type to 0
	if(metD$period$type==1)
	{	# so far in that case I have added NAs on February 29th, let's remove those
		feb29 <- metD$data$date[metD$data$mm==2 & metD$data$dd==29]
		for (l in length(feb29):1){
			w <- which(metD$data$date==feb29[l])
			metD$data$tmin <- c(metD$data$tmin[1:(w-1)],metD$data$tmin[(w+1):length(metD$data$tmin)]) 
			metD$data$tmax <- c(metD$data$tmax[1:(w-1)],metD$data$tmax[(w+1):length(metD$data$tmax)]) 
			metD$data$rain <- c(metD$data$rain[1:(w-1)],metD$data$rain[(w+1):length(metD$data$rain)]) 
		}
	}
	if(metD$period$type==2)
	{
		print("### > WARNING type == 2")
		print("    > nothing done for that yet")
		browser()
	}
	
	# missing values
	for (d in 6:8){		# only for tmin, tmax and ppt
		if(any(is.na(metD$data[[d]][]))){
			print("### > WARNING missing values",quote=F)
		}
		metD$data[[d]][is.na(metD$data[[d]][])] <- -999
	}	

	# folders
	if(!file.exists(paste(outFolder,"tmin",sep="/")))	dir.create(paste(outFolder,"tmin",sep="/"),showWarnings=TRUE,recursive=FALSE,mode="0777")
	if(!file.exists(paste(outFolder,"tmax",sep="/")))	dir.create(paste(outFolder,"tmax",sep="/"),showWarnings=TRUE,recursive=FALSE,mode="0777")
	if(!file.exists(paste(outFolder,"ppt",sep="/")))	dir.create(paste(outFolder,"ppt",sep="/"),showWarnings=TRUE,recursive=FALSE,mode="0777")
	
	# files
	outTmin <- paste(outFolder,"tmin",paste(metD$station$id,"txt",sep="."),sep="/")
	outTmax <- paste(outFolder,"tmax",paste(metD$station$id,"txt",sep="."),sep="/")
	outPpt <- paste(outFolder,"ppt",paste(metD$station$id,"txt",sep="."),sep="/")
	if(!file.copy("/home/crespo/Desktop/Optimisation/xpos-r/ClimateDataTools/ClimFormats/csagTemplates/csagTemplate.txt",outTmin,overwrite=TRUE)){
		print("### > FAILED TO CREATE OUTPUT FILE",quote=F)
		stop("cannot resolve")
	}
	if(Sys.info()["sysname"]=="Linux") newline <- "\r\n"

	# change header
	header <- readLines(outTmin)
	header <- sub("v_id",		metD$station$id,			header)
	header <- sub("v_lat",		round(metD$station$lat,2),		header)
	header <- sub("v_lon",		round(metD$station$lon,2),		header)
	header <- sub("v_alt",		round(metD$station$alt,0),		header)
	header <- sub("v_sDate",	format(metD$period$start,"%Y%m%d12"),	header)
	header <- sub("v_eDate",	format(metD$period$end,"%Y%m%d12"),	header)
	header <- sub("v_type",		metD$period$type,			header)
	header <- sub("v_description",	metD$station$comm,			header)
	write.table(header,outTmin,append=FALSE,quote=FALSE,col.names=FALSE,row.names=FALSE,eol=newline)
	# duplicate header
	if(!file.copy(outTmin,outTmax,overwrite=TRUE)){
		print("### > FAILED TO COPY TMIN HEADER INTO TMAX",quote=F)
		stop("cannot resolve")
	}
	if(!file.copy(outTmin,outPpt,overwrite=TRUE)){
		print("### > FAILED TO COPY TMIN HEADER INTO PPT",quote=F)
		stop("cannot resolve")
	}
	
	# add up body
	write.table(metD$data$tmin,outTmin,append=TRUE,quote=FALSE,col.names=FALSE,row.names=FALSE,eol=newline,sep=" ")
	write.table(metD$data$tmax,outTmax,append=TRUE,quote=FALSE,col.names=FALSE,row.names=FALSE,eol=newline,sep=" ")
	write.table(metD$data$rain,outPpt,append=TRUE,quote=FALSE,col.names=FALSE,row.names=FALSE,eol=newline,sep=" ")

rm(d,outTmin,outTmax,outPpt,header)
}


