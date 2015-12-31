##
 # FILE dataPerturbe.r
 # AUTHOR olivier crespo
 # https://r-forge.r-project.org/projects/xpos-r/
 ###############################################################################

## all format get out as
## see readFormats for updates

#list(	"file"=path to the file
#	"station"=list("id","lat","lon","alt","comm")
#	"clim"=list("tav","amp","refht","wndht")
#	"period"=list("start","end","type")	# type is CSAG style {0-real,1-365,2-360 days per year}
#	"data"=list("date","yyyy","mm","dd","juld","srad","tmax","tmin","rain","wind","dewp","vprs","rhum")
# MISSING VALUES = NA

## CHANGE THE PERIOD
########################################################################
#source('/home/crespo/Desktop/Optimisation/xpos-r/ClimateDataTools/Climatology/climAgro.r')
# metD for the major list
# sDate and eDate as.Date
pert_period <- function(metD,sDate,eDate)
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
			print("### > original start is later than asked starting date",quote=F)
			print("    > feeling with NA",quote=F)
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
			print("### > original end is earlier than asked ending date",quote=F)
			print("    > feeling with NA",quote=F)
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

########################################################################
## old : fill up an incomplete vector (1 var) set by masking a complete (1 var) set
## new : all data set
## lengths have to match
pert_mask <- function(NA_V,fullV)
{
#	if(class(NA_V)!="numeric" || class(fullV)!="numeric"){ # but a data set with only NA is class:"logical"
#		print("### > one of the data set is not numeric",quote=F)
#		stop("NAfill_mask: cannot resolve")
#	}
#	if (length(NA_V)!=length(fullV)){
#		print("### > the 2 data sets are not of the same lenght",quote=F)
#		stop("NAfill_mask: cannot resolve")
#	}
#	if (length(fullV[is.na(fullV)])>0){
#		print("### > there is missing values in the supposidely full set",quote=F)
#		print("    > keep going, we'll check at the end if you're lucky",quote=F)
#	}

	newV <- NA_V
	for (r in 1:length(newV$data$date)){
		if(is.na(newV$data$tmin[r])||is.na(newV$data$tmax[r])||is.na(newV$data$rain[r])){
			newV$data$tmin[r] <- fullV$data$tmin[r]
			newV$data$tmax[r] <- fullV$data$tmax[r]
			newV$data$rain[r] <- fullV$data$rain[r]		
		}
	}

	if(any(is.na(newV$data$tmin))||any(is.na(newV$data$tmax))||any(is.na(newV$data$rain))){
		print('missing value where there should not be')
		browser()
	}
return(newV)
rm(newV)
}

## PERTURBE TEMPERATURES
########################################################################
# metD for the major list
# pTmin and pTmax are dim=12 vectors of monthly shift
# oDif is the tmax - tmin diff on which these pTmin and pTmax occured
pert_temp <- function(metD,pTmin,pTmax=pTmin,tDif=NULL,check=FALSE)
{
	newD <- metD

	#
	# It appears in some cases that a simple shift ends up creating tmin > tmax
	# So I go for a linearly proportinal shift (related to the tmax-tmin difference)
	# Same for rain given the maximum rain occurence
	# N.B. mathematically it would not prevent tmin > tmax, but in practice, it should
	#
	# this is better I believe yet not the only reason
	# sometimes only tmax or only tmin is replaced, and even though the other set may be good
	# in some cases it ends up with tmin > tmax
	# at this stage it happens so rarely I did not take care of that one
	#
	oDif <- metD$data$tmax - metD$data$tmin
	if(is.null(tDif)){
		tDif <- array(NA,dim=12)
		for(m in 1:12){
			tDif[m] <- mean(oDif[metD$data$mm==m],na.rm=T)	# assume it is a simple shift
		}
	}

	# something wrong here here MOZ case ended up with 3 (10219,10220,10964) above 100's
	# in tmin case seems that a 2 of 20 is replaced by 13 or 14

	# change start
	newD$file <- paste("changed from :",metD$file,sep=" ")
	# tav and amp have to be recomputed, not refht, not wndht
	newD$clim$tav <- NA
	newD$clim$amp <- NA

	for(l in 1:length(newD$data$tmin)){
		if (is.na(newD$data$tmin[l])) next
		newD$data$tmin[l] <- round(newD$data$tmin[l]+pTmin[newD$data$mm[l]]*oDif[l]/tDif[newD$data$mm[l]],digits=1)
	}
	for(l in 1:length(newD$data$tmax)){
		if (is.na(newD$data$tmax[l])) next
		newD$data$tmax[l] <- round(newD$data$tmax[l]+pTmax[newD$data$mm[l]]*oDif[l]/tDif[newD$data$mm[l]],digits=1)
	}
	if(any(newD$data$tmin>newD$data$tmax,na.rm=T)){
		print(" ### > WARNING Tmin > Tmax",quote=F)
		browser()
	}

	if (check){
		m <- stat_monthlyMeans(metD$data)
		n <- stat_monthlyMeans(newD$data)
		print("### > actual monthly means shift",quote=F)
		actualChanges <- array(NA,dim=c(2,12))
		for (i in 1:12){
			actualChanges[1,i] <- mean(n$tmin[n$mm==i],na.rm=T)-mean(m$tmin[n$mm==i],na.rm=T)
			actualChanges[2,i] <- mean(n$tmax[n$mm==i],na.rm=T)-mean(m$tmax[n$mm==i],na.rm=T)
		}
		rownames(actualChanges) <- c('tmin','tmax')
		colnames(actualChanges) <- c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
		print(round(actualChanges,digits=1))
		print(paste("old tmin mean",mean(m$tmin,na.rm=T),"old tmax mean",mean(m$tmax,na.rm=T),sep=" : "),quote=FALSE)
		print(paste("new tmin mean",mean(n$tmin,na.rm=T),"new tmax mean",mean(n$tmax,na.rm=T),sep=" : "),quote=FALSE)
	}

return(newD)
rm(newD)
}


## PERTURBE RAINFALL
########################################################################
# metD for the major list
# pRain is a dim=12 vector of monthly factor change
pert_facRain <- function(metD,pRain,check=FALSE)
{
	newD <- metD

	# change start
	newD$file <- paste("changed from :",metD$file,sep=" ")

	for(l in 1:length(newD$data$rain)){
		if (is.na(newD$data$rain[l])) next
		newD$data$rain[l] <- round(newD$data$rain[l]*pRain[newD$data$mm[l]],digits=1)
		if(!is.na(newD$data$rain[l])){
			if(newD$data$rain[l]<0) newD$data$rain[l]<-0
		}
	}

	if (check){
		m <- stat_monthlyTotals(metD$data)
		n <- stat_monthlyTotals(newD$data)

		print("### > actual monthly factor changes",quote=F)
		actualChanges <- array(NA,dim=12)
		for (i in 1:12){
			actualChanges[i] <- mean(n$rain[n$mm==i],na.rm=T)/mean(m$rain[m$mm==i],na.rm=T)
		}
		rownames(actualChanges) <- c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
		print(round(actualChanges,digits=1))
		print(paste("old mean",mean(m$rain,na.rm=T),sep=" : "),quote=FALSE)
		print(paste("new mean",mean(n$rain,na.rm=T),sep=" : "),quote=FALSE)
	}

return(newD)
rm(newD,l,m_vec,check,m,n,actualChanges)
}

# metD for the major list
# pRain is a shift (single value)
pert_shiRain <- function(metD,shift,check=FALSE)
{
	m <- stat_monthlyTotals(metD$data)
	pRain <- array((mean(m$rain,na.rm=T)+shift)/mean(m$rain,na.rm=T),dim=12)
	newD <- pert_facRain(metD,pRain,check)

return(newD)
rm(newD,m,check,pRain)
}


########################################################################
# merge 2 data sets
# create a set of max(length(metD1),length(metD2))
# take tmin from metTmin, tmax from metTmax, and rain from metRain
source('/home/crespo/Desktop/Optimisation/xpos-r/ClimateDataTools/ClimFormats/dataRead.r')
pert_merge <- function(metTmin,metTmax,metRain)
{
	merD<-createNULLlist()
	merD$period$start<-min(metTmin$period$start,metTmax$period$start,metRain$period$start)
	merD$period$end<-max(metTmin$period$end,metTmax$period$end,metRain$period$end)

	# create daily date info
	merD$data$date <- 	seq(merD$period$start,merD$period$end,1)
	merD$data$yyyy <- 	as.numeric(format(merD$data$date,"%Y"))
	merD$data$mm <- 	as.numeric(format(merD$data$date,"%m"))
	merD$data$dd <- 	as.numeric(format(merD$data$date,"%d"))
	merD$data$juld <- 	as.numeric(format(merD$data$date,"%j"))
	merD$data$date <- 	as.numeric(format(merD$data$date,"%Y%j"))

	# initialisa tmin, tmax, rain
	merD$data$tmin <- 	array(NA,dim=length(merD$data$date))
	merD$data$tmax <- 	array(NA,dim=length(merD$data$date))
	merD$data$rain <- 	array(NA,dim=length(merD$data$date))
	merD$data$srad <- 	array(NA,dim=length(merD$data$date))
	merD$data$wind <- 	array(NA,dim=length(merD$data$date))
	merD$data$dewp <- 	array(NA,dim=length(merD$data$date))
	merD$data$vprs <- 	array(NA,dim=length(merD$data$date))
	merD$data$rhum <- 	array(NA,dim=length(merD$data$date))

	# NA for everything else
	merD$station$lat <- 	NA
	merD$station$lon <- 	NA
	merD$station$alt <- 	NA
	merD$clim$tav <- 	NA
	merD$clim$amp <- 	NA
	merD$clim$refht <- 	NA
	merD$clim$wndht <- 	NA

	#metTmin
	sLin <- grep(metTmin$data$date[1],merD$data$date)
	merD$data$tmin[sLin:(sLin+length(metTmin$data$tmin)-1)] <- metTmin$data$tmin[1:length(metTmin$data$tmin)]
		
	#metTmax
	sLin <- grep(metTmax$data$date[1],merD$data$date)
	merD$data$tmax[sLin:(sLin+length(metTmin$data$tmax)-1)] <- metTmax$data$tmax[1:length(metTmax$data$tmax)]

	#metRain
	sLin <- grep(metRain$data$date[1],merD$data$date)
	merD$data$rain[sLin:(sLin+length(metRain$data$rain)-1)] <- metRain$data$rain[1:length(metRain$data$rain)]

return(merD)
rm(merD)
}
