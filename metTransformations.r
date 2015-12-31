##
 # FILE metTransformations.r
 # AUTHOR olivier crespo
 # https://r-forge.r-project.org/projects/xpos-r/convertBruceFormat
 ###############################################################################
 #
 # functions required to adapt, complete transform weather data
 #
 ###############################################################################


## OBSOLETE? see climAgro.r
 # COMPUTE tav AND amp APSIM CONSTANTS
 # annual average ambient temperature (TAV)
 # annual amplitude in mean monthly temperature (AMP)
 ###############################################################################
 # > ref
 # http://www.apsim.info/Wiki/public/Upload/OtherProducts/tav_amp.pdf
 # results confirmed in face of the tav_amp.exe dos application
 ###############################################################################
 # you need daily data for (tmn,tmx,year,julDay), here in the 'data' list
 ###############################################################################
compute_tavNamp <- function(data)
{

### AMP is limited to complete year only since for uncomplete ones yearlyAMP will have NA values, then mean is processed with na.rm=T
###
	# daily mean
	dMean <- (data$tmn+data$tmx)/2;

	firstYear <- as.numeric(data$year[1]);
	yearlyAMP <- array(NA,dim=(as.numeric(data$year[dim(data$year)])-as.numeric(data$year[1])+1));
	monthlyMean <- array(0,dim=c(12,5));
	year <- firstYear;
	for (line in 1:length(data$year)){
		# yearly AMP
		if (data$year[line]!=year){
			yearlyAMP[year-firstYear+1] <- max(monthlyMean[,3]/monthlyMean[,4])-min(monthlyMean[,3]/monthlyMean[,4]);
			year <- as.numeric(data$year[line]);
			monthlyMean[,3] <- 0;
			monthlyMean[,4] <- 0;
		}
		# monthly mean
		month <- as.numeric(format(as.Date(data$julDay[line]-1,origin=paste(data$year[line],"-01-01",sep="")),"%m"));
		if (!is.na(dMean[line])){ # avoid NAs
			monthlyMean[month,1] <- monthlyMean[month,1]+dMean[line];
			monthlyMean[month,2] <- monthlyMean[month,2]+1;
			monthlyMean[month,3] <- monthlyMean[month,3]+dMean[line];
			monthlyMean[month,4] <- monthlyMean[month,4]+1;
		}
	}
	
	# complete
	yearlyAMP[year-firstYear+1] <- max(monthlyMean[,3]/monthlyMean[,4])-min(monthlyMean[,3]/monthlyMean[,4]);
	monthlyMean[,5] <- monthlyMean[,1]/monthlyMean[,2];

	amp <- ifelse(all(is.na(yearlyAMP)),NA,format(mean(yearlyAMP,na.rm=T),digits=3));
	tav <- ifelse(all(is.na(monthlyMean[,5])),NA,format(mean(monthlyMean[,5]),digits=3));

	if(amp!="NA" && (as.numeric(amp)<0 || as.numeric(amp)>50))	{
		print("# WARNING: amp is not in the APSIM range allowed (0>amp>50)",quote=FALSE);
	}
	if(tav!="NA" && (as.numeric(tav)<0 || as.numeric(tav)>50))	{
		print("# WARNING: tav is not in the APSIM range allowed (0>tav>50)",quote=FALSE);
	}

	data<-list("tmn"=data$tmn,"tmx"=data$tmx,"ppt"=data$ppt,"year"=data$year,"julDay"=data$julDay,"sRad"=data$sRad,"amp"=amp,"tav"=tav);
return(data);
}
##
 # CREATE YEAR AND JULIAN DAYS
 ###############################################################################
createYearJulianDays <- function(data,fileHead)
{
	firstYear <- format(fileHead$period$start,"%Y");
	lastYear <- format(fileHead$period$end,"%Y");
	firstDay <- fileHead$period$start-as.Date(paste(firstYear,"01","01",sep="-"),"%Y-%m-%d")+1;
	lastDay <- fileHead$period$end-as.Date(paste(lastYear,"01","01",sep="-"),"%Y-%m-%d")+1;
	
	apsim_year <- array(firstYear,dim=1);
	apsim_julDay <- array(firstDay,dim=1);
	li <- 1;
	for (y in firstYear:lastYear){
		if (y == firstYear) d <- firstDay;
		repeat{
			# boundaries conditions
			if(y==lastYear && d==lastDay)	break;
			if(!is.leapYear(y) && d==365)	{d<-0; break;}
			if(d==366)	{d<-0; break;}
			d <- d+1;
			li <- li+1;
			# apsim table update
			apsim_year <- array(c(apsim_year,y),dim=dim(apsim_year)+1); 
			apsim_julDay <- array(c(apsim_julDay,d),dim=dim(apsim_julDay)+1);
		}
	}

	data <- list("tmn"=data$tmn,"tmx"=data$tmx,"ppt"=data$ppt,"year"=apsim_year,"julDay"=apsim_julDay);
return(data)
}

##
 # is that year a leap year?
 ###############################################################################
is.leapYear <- function(year)
{
	## trusting R date classe
	## otherwise there would be a conflict sooner or later anyway
	start <- as.Date(paste("01","01",year,sep="-"),"%d-%m-%Y");
	end <- as.Date(paste("31","12",year,sep="-"),"%d-%m-%Y");

	dayNo <- end-start +1;
	switch(dayNo-364,
		leap <- FALSE,
		leap <- TRUE
	);
return(leap);
}

##
 # TRANSFORM 365 DAYS A YEAR INTO 366
 ###############################################################################
 # WARNING : if you change the day added
 # go and check in "removeAddedDays" fct as well
 ###############################################################################
transform_365into366 <- function(oldYear,fillIn)
{
     ## my guess
	# 365 days is for real years without leap years
	# i.e. with no 29th of February
	# to be confirmed
	
     ## add 1 day in the middle of february
	# ppt is 0
	# tmin and tmax are copied (randomly) either from previous or following day

	# 01-01 +44 days is the 14-02
	# 01-01 +45 days is the 15-02
	ifelse(runif(1)<0.5,day<-44,day<-45);
	new1 <- oldYear[day,1];
	new2 <- oldYear[day,2];
	new3 <- 0;

	# 01-01 +44 as before
	# 15-02 is new
	# 01-01 +45 starts the 16-02 
	newYear <- oldYear[1:44,1:dim(oldYear)[2]];
	ifelse(fillIn,newYear <- rbind(newYear,c(new1,new2,new3)),newYear <- rbind(newYear,c(NA,NA,NA)));
	newYear <- rbind(newYear,oldYear[45:dim(oldYear)[1],1:dim(oldYear)[2]]);

return(newYear);
}

##
 # TRANSFORM 360 DAYS A YEAR INTO REAL DAYS No
 ###############################################################################
transform_360intoREAL <- function(oldYear,year,fillIn)
{
     ## leap year: add 1 day in the middle of jan, mar, may, jul, sep, nov
     ## non leap year: add 1 day in the middle of mar, may, jul, sep, nov
	# ppt is 0
	# tmin and tmax are copied (randomly) either from previous or following day

	ifelse (is.leapYear(year),step<-52,step<-60);
	d<-step;
	while(d<(365-step)){
		ifelse(runif(1)<0.5,day<-d,day<-d+1);
		new1 <- oldYear[day,1];
		new2 <- oldYear[day,2];
		new3 <- 0;

		newYear <- oldYear[1:d,1:dim(oldYear)[2]];
		ifelse(fillIn,newYear <- rbind(newYear,c(new1,new2,new3)),newYear <- rbind(newYear,c(NA,NA,NA)));
		newYear <- rbind(newYear,oldYear[(d+1):dim(oldYear)[1],1:dim(oldYear)[2]]);
		oldYear <- newYear;

		d <- d+step;
	}

return(newYear);
}

##
 # CHECK THAT DIM OF THE TABLE FITS THE REAL No OF DAYS
 ###############################################################################
check_dayVSdim <- function(sDate,eDate,linNo)
{
	dayNo <- eDate-sDate +1;
	if (dayNo != linNo){
		stop("*** wrong number of Days: check_dayVSdim (metTransformation.r)");
	}
}

##
 # TRANSFORM DAY No OVER A PERIOD TO REAL _ type 365 days
 ###############################################################################
transform_type1 <- function(data,head,fillIn)
{
	table <- data$tmn;
	table <- cbind(table,data$tmx);
	table <- cbind(table,data$ppt);

	for(y in format(head$period$start,"%Y"):format(head$period$end,"%Y")){
		dayNo_bef <- as.Date(paste("01","01",y,sep="-"),"%d-%m-%Y") - head$period$start;
		if (is.leapYear(y) && dayNo_bef>=0){
			# cut before
			if (dayNo_bef==0)	table_bef <- NULL
			if (dayNo_bef>0)	table_bef <- table[1:dayNo_bef,1:dim(table)[2]]
			# pull out section
			oldSection <- table[(dayNo_bef+1):(dayNo_bef+365),1:dim(table)[2]];
			# cut after
			if(y!=format(head$period$end,"%Y")){
				table_aft <- table[(dayNo_bef+365+1):dim(table)[1],1:dim(table)[2]];
			}else{ # no after
				table_aft <- NULL;
			}
			# transform
			newSection <- transform_365into366(oldSection,fillIn);
			# paste				if(data$julDay[line]==4*step) removeTheLine();

			table <- table_bef;
			table <- rbind(table,newSection);
			table <- rbind(table,table_aft);
		}
	}
	check_dayVSdim(head$period$start,head$period$end,dim(table)[1]);

	data$tmn <- as.array(table[,1],dim=dim(table)[1]);
	data$tmx <- as.array(table[,2],dim=dim(table)[1]);
	data$ppt <- as.array(table[,3],dim=dim(table)[1]);
return(data);
}

##
 # TRANSFORM DAY No OVER A PERIOD TO REAL _ type 360 days
 ###############################################################################
transform_type2 <- function(data,head,fillIn)
{
	table <- data$tmn;
	table <- cbind(table,data$tmx);
	table <- cbind(table,data$ppt);

#browser();
	for(y in format(head$period$start,"%Y"):format(head$period$end,"%Y")){
		# cut before
		dayNo_bef <- as.Date(paste("01","01",y,sep="-"),"%d-%m-%Y") - head$period$start;
		if(dayNo_bef>0){
			table_bef <- table[1:dayNo_bef,1:dim(table)[2]];
		}else{
			table_bef <- NULL;
		}
		# pull out section
		oldSection <- table[(dayNo_bef+1):(dayNo_bef+360),1:dim(table)[2]];
		# cut after
		if(y!=format(head$period$end,"%Y")){
			table_aft <- table[(dayNo_bef+360+1):dim(table)[1],1:dim(table)[2]];
		}else{ # no after
			table_aft <- NULL;
		}
		# transform
		newSection <- transform_360intoREAL(oldSection,y,fillIn);
		# paste
		table <- table_bef;
		table <- rbind(table,newSection);
		table <- rbind(table,table_aft);
	}
	check_dayVSdim(head$period$start,head$period$end,dim(table)[1]);

	data$tmn <- as.array(table[,1],dim=dim(table)[1]);
	data$tmx <- as.array(table[,2],dim=dim(table)[1]);
	data$ppt <- as.array(table[,3],dim=dim(table)[1]);
return(data);
}

##
 # REMOVE days added earlier for radn and ETo estimation
 ###############################################################################
 # julian day is required for radn and ETo computation
 # I decided to first make real years, compute those and then remove the added days
 #
 # looking for NA is not satifying as there might be others or even none of fillIn wasn't specified
 # depending on type, look for days above 360 or leap year's day=366
 ###############################################################################
removeTheLine <- function(data,line)
{
	prior <- NULL; post <- NULL;

	att <- 0;
	while(att < length(data)){
		att <- att+1;
		if (is.null(dim(data[[att]])))	next;
		# cut down
		if(line>1)
			prior <- array(data[[att]][1:(line-1)],dim=(line-1));
		if(line<(dim(data[[att]])-1))
			post <- array(data[[att]][(line+1):dim(data[[att]])],dim=(dim(data[[att]])-line));
		
		# paste
		if(is.null(prior)){
			data[[att]] <- array(post,dim=dim(post)); next;
		}
		if(is.null(post)){
			data[[att]] <- array(prior,dim=dim(prior)); next;
		}
		data[[att]] <- array(c(prior,post),dim=(dim(prior)+dim(post)));
	}

	data <- list(	"tmn"=data$tmn,
			"tmx"=data$tmx,
			"ppt"=data$ppt,
			"julDay"=data$julDay,
			"year"=data$year,
			"sRad"=data$sRad,
			"ETo"=data$ETo,
			"AI"=data$AI,
			"doItAgain"=data$doItAgain,
			"arid"=data$arid
		);

return(data);
}
removeAddedDays <- function(data,head)
{
	# so far I'm using that only for csag format, which does not include year and julDay
	# so I do not really care of the consistence of my julian Day...
	# yet you might!
	type <- head$period$type+1;
	line <- 0;
	while (line < dim(data$year)){
		line <- line+1;
		switch(type,
			{	# 0: real
			},{	# 1: 365
				if(is.leapYear(data$year[line]) && data$julDay[line]==45)
					data<-removeTheLine(data,line);
			},{	# 2: 360
##
## never checked so far
##
				ifelse (is.leapYear(data$year[line]),step<-52,step<-60);
				if(data$julDay[line]==step) data<-removeTheLine(data,line);
				if(data$julDay[line]==2*step) data<-removeTheLine(data,line);
				if(data$julDay[line]==3*step) data<-removeTheLine(data,line);
				if(data$julDay[line]==4*step) data<-removeTheLine(data,line);
				if(data$julDay[line]==5*step) data<-removeTheLine(data,line);
				if(is.leapYear(data$year[line]) && data$julDay[line]==6*step) data<-removeTheLine(data,line);
			}
		);
	}
	data <- list(	"tmn"=data$tmn,
			"tmx"=data$tmx,
			"ppt"=data$ppt,
			"julDay"=data$julDay,
			"year"=data$year,
			"sRad"=data$sRad,
			"ETo"=data$ETo,
			"AI"=data$AI,
			"doItAgain"=data$doItAgain,
			"arid"=data$arid
		);

return(data);
}

##
 # CHANGE an old known variable (characters string) into a new one
 ###############################################################################
 # not really attached to bruce data
 # but useful for all, and bruce data functions are used by all
 ###############################################################################
changeVar <- function(oldVar,newVar,oldFile,newFile)
{	
	tempFile <- readLines(oldFile,n=-1,warn=FALSE);
	tempFile <- sub(oldVar, newVar, tempFile);
	# for double occurences on the same line (e.g. sowing date at day)
	tempFile <- sub(oldVar, newVar, tempFile);
	writeLines(tempFile,newFile,sep="\n");
}

