##
 # FILE convertToApsim.r
 # AUTHOR olivier crespo
 # https://r-forge.r-project.org/projects/xpos-r/convertBruceFormat
 ###############################################################################
 #
 # functions required for APSIM format only
 #
 ###############################################################################

# NOW IN METTRANSFORMATIONS.R
##
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
#compute_tavNamp <- function(data)
#{
#	# daily mean
#	dMean <- (data$tmn+data$tmx)/2;
#
#	firstYear <- as.numeric(data$year[1]);
#	yearlyAMP <- array(NA,dim=(as.numeric(data$year[dim(data$year)])-as.numeric(data$year[1])+1));
#	monthlyMean <- array(0,dim=c(12,5));
#	year <- firstYear;
#	for (line in 1:dim(data$year)){
#		# yearly AMP
#		if (data$year[line]!=year){
#			yearlyAMP[year-firstYear+1] <- max(monthlyMean[,3]/monthlyMean[,4])-min(monthlyMean[,3]/monthlyMean[,4]);
#			year <- as.numeric(data$year[line]);
#			monthlyMean[,3] <- 0;
#			monthlyMean[,4] <- 0;
#		}
#		# monthly mean
#		month <- as.numeric(format(as.Date(data$julDay[line]-1,origin=paste(data$year[line],"-01-01",sep="")),"%m"));
#		monthlyMean[month,1] <- monthlyMean[month,1]+dMean[line];
#		monthlyMean[month,2] <- monthlyMean[month,2]+1;
#		monthlyMean[month,3] <- monthlyMean[month,3]+dMean[line];
#		monthlyMean[month,4] <- monthlyMean[month,4]+1;
#	}
#	
#	# complete
#	yearlyAMP[year-firstYear+1] <- max(monthlyMean[,3]/monthlyMean[,4])-min(monthlyMean[,3]/monthlyMean[,4]);
#	monthlyMean[,5] <- monthlyMean[,1]/monthlyMean[,2];
#
#	amp <- format(mean(yearlyAMP),digits=3);
#	tav <- format(mean(monthlyMean[,5]),digits=3);
#
#	if(amp!="NA" && (as.numeric(amp)<0 || as.numeric(amp)>50))	{
#		print("# WARNING: amp is not in the APSIM range allowed (0>amp>50)",quote=FALSE);
#	}
#	if(tav!="NA" && (as.numeric(tav)<0 || as.numeric(tav)>50))	{
#		print("# WARNING: tav is not in the APSIM range allowed (0>tav>50)",quote=FALSE);
#	}
#
#	data<-list("tmn"=data$tmn,"tmx"=data$tmx,"ppt"=data$ppt,"year"=data$year,"julDay"=data$julDay,"sRad"=data$sRad,"amp"=amp,"tav"=tav);
#return(data);
#}

##
 # FORMAT AND WRITE DATA INTO OUTPUT .met APSIM FILE
 ###############################################################################
formatToMetFile <- function(data,fileHead,path)
{
# waiting for a nicer solution within apsim ( at least necessary until apsim 7.0 )
# make 'year' a fake, and 'realY' the actual year
# still one issue: 2100 is not a leap year, and any fake within the last century (with my handmade solution) will be leap
# solution so far: do not simulate 2100 with apsim !!
	fake_year <- array(as.numeric(data$year)-ifelse(format(fileHead$period$end,"%Y")>=2065,100,0),dim=dim(data$year));

# make one table from all the data
#	apsim_table <- as.numeric(fake_year);				# fake year titled	'year'
# with RCP we'll have some problems
	apsim_table <- as.numeric(data$year);				# fake year titled	'year'
	apsim_table <- cbind(apsim_table,as.integer(data$year));	# real year titled	'realY'
	apsim_table <- cbind(apsim_table,as.integer(data$julDay));	# julian day titled	'day'
	apsim_table <- cbind(apsim_table,as.numeric(data$tmn));		# temp min titled	'mint'
	apsim_table <- cbind(apsim_table,as.numeric(data$tmx));		# temp max titled	'maxt'
	apsim_table <- cbind(apsim_table,as.numeric(data$ppt));		# precipitation titled	'rain'
	apsim_table <- cbind(apsim_table,as.numeric(data$sRad));	# radiation titled	'mint'
	apsim_table[,4:7] <- format(apsim_table[,4:7],digits=3);

# write it into a .met file
	outName <- paste(path$output,strsplit(path$file$temp,"txt")[[1]][1],"met",sep="");
	if(!file.exists(path$output)){	# create output dir if does not exists
		dir.create(path$output, showWarnings=F, recursive=T, mode = "0777");
	}else{				# remove output file if does exists
		if(file.exists(outName))	file.remove(outName);
	}
	# head
	file.copy(	"./apsimTemplates/apsimTemplate.met",	outName,overwrite=TRUE);
	changeVar(	"station_id",	fileHead$station$id,	outName,outName);
	changeVar(	"station_comm",	fileHead$comment,	outName,outName);
	changeVar(	"stat_lat",	fileHead$station$lat,	outName,outName);
	changeVar(	"stat_lon",	fileHead$station$lon,	outName,outName);
	changeVar(	"stat_alt",	fileHead$station$alt,	outName,outName);
	changeVar(	"inLand",	ifelse(path$inland,"0.16 (in land)","0.19 (coastal)"),	outName,outName);
	changeVar(	"period_tav",	data$tav,		outName,outName);
	changeVar(	"period_amp",	data$amp,		outName,outName);
	# body
	apsim_table <- format(apsim_table,justify="right",width=7);
	write.table(apsim_table,outName,quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE);

	## to handle the linux/windows format issue
	if(Sys.info()["sysname"]=="Linux"){
		tempFile <- readLines(outName,n=-1,warn=FALSE);
		writeLines(tempFile,outName,sep="\r\n");
	}
}

