# AUTHOR EVE
# DATE 29-DEC-2015
# FILE R_eve.r
############

source('R_readMetData.r')
source('R_writeMetData.r')
source('R_climAgro.r')

findMonthMean <- function(month,dataset,var='tmin',period=seq(1980,1999,1)){
# tmin column is 6
# tmax column is 7
return(mean(dataset$data[[grep(var,names(dataset$data))]][dataset$data$mm==month & !is.na(match(dataset$data$yyyy,period))]))
}

computeMonthMean <- function(dataset,var='tmax'){
	avVec <- array(NA,dim=12)
	for(i in 1:12){
		avVec[i] <- findMonthMean(i,dataset,var)
	}
return(avVec)
}

diffAll <- function(observed,agmerra,var='tmax'){

	differences <- array(NA,dim=12)
	differences <- computeMonthMean(agmerra,var)-computeMonthMean(observed,var)
return(differences)
}


## tests
obsFi <- '08940XXX.AgMIP'
agmFile <- 'F3430QXX.AgMIP'
obsInFo <- "/home/olivier/Desktop/Wine-shared/Projects/2015-2016_AgMIP2/RZA-allFS/Obshis_agmip"
agmInFo <- "/home/olivier/Desktop/Wine-shared/Projects/2015-2016_AgMIP2/RZA-allFS/AgMERRA"
outFo <- "/home/olivier/Desktop/Wine-shared/Projects/2015-2016_AgMIP2/RZA-allFS/RZA-allFS-Baselines"
#
# DO IT ALL IN THERE
#
# inputs:
# outputs:
biasCorrect <- function(agmInFolder=agmInFo, obsInFolder=obsInFo, outFolder=outFo){#, agmFile=agmFi, obsFile=obsFi

	obsFiles_l <- list.files(obsInFolder)
	for(obsF in obsFiles_l[391:(length(obsFiles_l))]){
print(paste(" processing observation ",obsF," (",match(obsF,obsFiles_l),"/",length(obsFiles_l),")",sep=""),quote=F)
		# read in
		observed <- read_AgMIPformat(inFile=paste(obsInFolder,obsF,sep="/"))

		# observed$station$lat
		# observed$station$lon
		# observed$station$alt
		agmFiles_l <- list.files(agmInFolder)
		for(agmF in agmFiles_l){
			agmLon <- as.numeric(scan(paste(agmInFolder,agmF,sep="/"),what="raw",skip=3,nlines=1,quiet=T)[3])
			
			# if i find the right LON
			if(observed$station$lon<=(agmLon+0.25/2) && observed$station$lon>=(agmLon-0.25/2)){
				agm_i <- match(agmF, agmFiles_l)				
				agmLat <- as.numeric(scan(paste(agmInFolder,agmFiles_l[agm_i],sep="/"),what="raw",skip=3,nlines=1,quiet=T)[2])
				# while I did not found LAT				
				while(observed$station$lat>(agmLat+0.25/2) || observed$station$lat<(agmLat-0.25/2)){
					agm_i <- agm_i +23
					if (agm_i>length(agmFiles_l)){
						print(" # >> Could not find AgMERRA files corresponding",sep="")
						browser()
					}
					agmLat <- as.numeric(scan(paste(agmInFolder,agmFiles_l[agm_i],sep="/"),what="raw",skip=3,nlines=1,quiet=T)[2])	
				}
				break # found LAT and LON
			}
		}

print(paste("with AgMERRA ",agmFiles_l[agm_i]," (",agm_i,"/",length(agmFiles_l),")",sep=""),quote=F)
		agmBiased<- agmerra <- read_AgMIPformat(inFile=paste(agmInFolder,agmFiles_l[agm_i],sep="/"))

		# set up obsBiased as period-cut observed set
		obsBiased <- pert_period(observed,as.Date("1980-01-01"),as.Date("2010-12-31"),verb=F)

		# compute AgMERRA fully biased
		minbias <- diffAll(observed,agmerra,'tmin')
		maxbias <- diffAll(observed,agmerra,'tmax')
		radbias <- diffAll(observed,agmerra,'srad')
		for (i in 1:12){
			#tmin	
			agmBiased$data$tmin[agmBiased$data$mm==i] <- round((agmerra$data$tmin[agmerra$data$mm==i])-minbias[i], digits=2)
			#tmax
			agmBiased$data$tmax[agmBiased$data$mm==i] <- round((agmerra$data$tmax[agmerra$data$mm==i])-maxbias[i], digits=2)
			#srad
			agmBiased$data$srad[agmBiased$data$mm==i] <- round((agmerra$data$srad[agmerra$data$mm==i])-radbias[i], digits=1)
		}

		# replace NAs in tmin, tmax, rain by AgMERRA biased
		for(v in 6:9){ #tmin:6 tmax:7 rain:8 srad:9
			for (r in 1:length(obsBiased$data[[v]])){
				if (!is.na(obsBiased$data[[v]][r])){
					next
				}else{
					obsBiased$data[[v]][r] <- agmBiased$data[[v]][r]
				}
			}
		}	

		# re-compute TAV and AMP
		obsBiasedWithTAVnAMP <- agro_tavamp(obsBiased)
	
		# comments
		obsBiasedWithTAVnAMP$station$comm <- paste(obsBiasedWithTAVnAMP$station$comm," amended with AgMERRA biased corrected",sep=" ")
	
	write_AgMIPformat(obsBiasedWithTAVnAMP,outFile=paste(outFolder,obsF,sep="/"),checkMiss=F)
	}
}

## CHANGE THE PERIOD
########################################################################
# metD for the major list
# sDate and eDate as.Date
pert_period <- function(metD,sDate,eDate,verb=T)
{
	newD <- metD

	# change start
	if(metD$period$start==sDate){	# nothing
	}else{
		newD$file <- paste("changed from :",metD$file,sep=" ")
		# tav and amp have to be recomputed, not refht, not wndht
		newD$clim$tav <- NA
		newD$clim$amp <- NA
		newD$period$start <- sDate

		# add up NAs
		if(metD$period$start>sDate){
			# if non existing, add NA
			if(verb)	print("### > original start is later than asked starting date",quote=F)
			if(verb)	print("    > feeling with NA",quote=F)
			dayDiff <- difftime(metD$period$start,newD$period$start,units="days")
			for (i in 1:length(newD$data)){
				newD$data[[i]] <- c(array(NA,dim=dayDiff),newD$data[[i]])
			}
			
			# can fill date related fields
			v <- sDate+(0:(dayDiff-1))
			newD$data$yyyy[1:dayDiff] <- as.numeric(format(v,"%Y"))
			newD$data$mm[1:dayDiff] <- as.numeric(format(v,"%m"))
			newD$data$dd[1:dayDiff] <- as.numeric(format(v,"%d"))
			newD$data$juld[1:dayDiff] <- as.numeric(format(v,"%j"))
			newD$data$date[1:dayDiff] <- as.numeric(format(v,"%Y%j"))

		# cut down the set
		}else{
			dayDiff <- difftime(newD$period$start,metD$period$start,units="days")
			for (i in 1:length(newD$data)){
				newD$data[[i]] <- newD$data[[i]][(dayDiff+1):length(newD$data[[i]])]
			}
		}
	}

	# change end
	if(metD$period$end==eDate){	# nothing
	}else{
		newD$file <- paste("changed from :",metD$file,sep=" ")
		# tav and amp have to be recomputed, not refht, not wndht
		newD$clim$tav <- NA
		newD$clim$amp <- NA
		newD$period$end <- eDate

		# add up NAs
		if(metD$period$end<eDate){
			# if non existing, add NA
			if(verb)	print("### > original end is earlier than asked ending date",quote=F)
			if(verb)	print("    > feeling with NA",quote=F)
			dayDiff <- difftime(newD$period$end,metD$period$end,units="days")
			for (i in 1:length(newD$data)){
				newD$data[[i]] <- c(newD$data[[i]],array(NA,dim=dayDiff))
			}
			# can fill date related fields
			v <- eDate+((-dayDiff+1):0)
			newD$data$yyyy[(length(newD$data$yyyy)-dayDiff+1):length(newD$data$yyyy)] <- as.numeric(format(v,"%Y"))
			newD$data$mm[(length(newD$data$mm)-dayDiff+1):length(newD$data$mm)] <- as.numeric(format(v,"%m"))
			newD$data$dd[(length(newD$data$dd)-dayDiff+1):length(newD$data$dd)] <- as.numeric(format(v,"%d"))
			newD$data$juld[(length(newD$data$juld)-dayDiff+1):length(newD$data$juld)] <- as.numeric(format(v,"%j"))
			newD$data$date[(length(newD$data$date)-dayDiff+1):length(newD$data$date)] <- as.numeric(format(v,"%Y%j"))

		# cut down the set
		}else{
			dayDiff <- difftime(metD$period$end,newD$period$end,units="days")
			for (i in 1:length(newD$data)){
				newD$data[[i]] <- newD$data[[i]][1:(length(newD$data[[i]])-dayDiff)]
			}
		}
	}

	# test lengths
	periodL <- difftime(newD$period$end,newD$period$start,units="days")+1
	for(i in 1:length(newD$data)){
		if(length(newD$data[[i]])!=periodL){
			print("### > data length issue (see cutPasteData.r)",quote=F)
			print("    > please remember that if you used GCM data, you may have uncomplete year length (type>0)",quote=F)
			browser()
		}				
	}

return(newD)
rm(newD,periodL,dayDiff,v)
}

