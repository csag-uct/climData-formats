##
 # FILE convertToDssat.r
 # AUTHOR olivier crespo
 # https://r-forge.r-project.org/projects/xpos-r/convertBruceFormat
 ###############################################################################
 #
 # functions required for DSSAT format only
 #
 ###############################################################################

##
 # FORMAT AND WRITE DATA INTO OUTPUT .WTH DSSAT FILE
 ###############################################################################
formatToWTHfile <- function(data,fileHead,path)
{
	date_table <- array(c(data$year,data$julDay),dim=c(dim(data$year),2));
	date_table <- cbind(date_table,NA);
	for(i in 1:dim(date_table)[1]){
		## fill in with zeros where needed
		if(as.numeric(date_table[i,2])<10){
			date_table[i,3] <- paste(date_table[i,1],date_table[i,2],sep="00");
		}
		if(as.numeric(date_table[i,2])>9 && as.numeric(date_table[i,2])<100){
			date_table[i,3] <- paste(date_table[i,1],date_table[i,2],sep="0");
		}	
		if(as.numeric(date_table[i,2])>99){
			date_table[i,3] <- paste(date_table[i,1],date_table[i,2],sep="");
		}
		
		## remove 2 first digits
		tmp <- strsplit(date_table[i,3],split=NULL);
		date_table[i,3] <- paste(tmp[[1]][3],tmp[[1]][4],tmp[[1]][5],tmp[[1]][6],tmp[[1]][7],sep="");
	}

# make one table from all the data
	dssat_table <- array(date_table[,3],dim=dim(data$year));		# year and jul day titled 'DATE'
	dssat_table <- cbind(dssat_table,round(as.numeric(data$sRad)*10)/10);	# radiation titled	'SRAD'
	dssat_table <- cbind(dssat_table,round(as.numeric(data$tmx)*10)/10);				# temp max titled	'TMAX'
	dssat_table <- cbind(dssat_table,round(as.numeric(data$tmn)*10)/10);				# temp min titled	'TMIN'
	dssat_table <- cbind(dssat_table,round(as.numeric(data$ppt)*10)/10);				# precipitation titled	'RAIN'

# write it into a .WTH file
	outName <- paste(path$output,strsplit(path$file$temp,"\\.")[[1]][1],".WTH",sep="");
	if(!file.exists(path$output)){	# create output dir if does not exists
		dir.create(path$output, showWarnings = FALSE, recursive = TRUE, mode = "0777");
	}else{				# remove output file if does exists
		if(file.exists(outName))	file.remove(outName);
	}
	# head
	file.copy(	"/home/crespo/Desktop/Optimisation/xpos-r/ClimateDataTools/ClimFormats/dssatTemplates/dssatTemplate.WTH",outName,overwrite=TRUE);
	changeVar(	"station_id",	fileHead$station$id,			outName,outName);
	changeVar(	"station_comm",	fileHead$comment,			outName,outName);
	changeVar(	"stat_lat",	format(fileHead$station$lat,nsmall=3),	outName,outName);
	changeVar(	"stat_lon",	format(fileHead$station$lon,nsmall=3),	outName,outName);
	changeVar(	"stat_alt",	format(fileHead$station$alt,digits=4),			outName,outName);
	changeVar(	"inLand",	ifelse(path$inland,"0.16 (in land)","0.19 (coastal)"),	outName,outName);
	changeVar(	"period_tav",	format(data$tav,nsmall=1),		outName,outName);
	changeVar(	"period_amp",	format(data$amp,nsmall=1),		outName,outName);
	
	# body
	tmp <- format(t(dssat_table),justify="right",width=4)
#	write.table(tmp,outName,quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE);
	write(tmp,outName,ncolumns=5,append=TRUE);

	## to handle the linux/windows format issue
	if(Sys.info()["sysname"]=="Linux"){
		tempFile <- readLines(outName,n=-1,warn=FALSE);
		writeLines(tempFile,outName,sep="\r\n");
	}
}

