##
 # FILE convertToAquaCrop.r
 # AUTHOR olivier crespo
 # https://r-forge.r-project.org/projects/xpos-r/convertBruceFormat
 ###############################################################################
 #
 # functions required for AQUACROP format only
 #
 ###############################################################################

##
 # FORMAT AND WRITE DATA INTO OUTPUT .TMP AQUACROP FILE
 ###############################################################################
formatToTMPFile <- function(data,fileHead,path)
{
# I do not know yet if aquacrop is dealing with future date
# I hope so :)

# make one table from all the data
	table <- array(as.numeric(data$tmn),dim=dim(data$tmn));	# temp min titled	'Tmin'
	table <- cbind(table,as.numeric(data$tmx));		# temp max titled	'Tmax'

# write it into a .TMP file
	outName <- paste(path$output,strsplit(path$file$temp,"\\.")[[1]][1],".TMP",sep="");
	if(!file.exists(path$output)){	# create output dir if does not exists
		dir.create(path$output, showWarnings = FALSE, recursive = FALSE, mode = "0777");
	}else{				# remove output file if does exists
		if(file.exists(outName))	file.remove(outName);
	}
	# head
	file.copy(	"./aquaCropTemplates/aquaCropTemplate.TMP",	outName,overwrite=TRUE);
	changeVar(	"DESCRIPTION",	fileHead$station$id,		outName,outName);
	changeVar(	"FD",		format(fileHead$period$start,"%d"),	outName,outName);
	changeVar(	"FM",		format(fileHead$period$start,"%m"),	outName,outName);
	changeVar(	"FY",		format(fileHead$period$start,"%Y"),	outName,outName);
	# body
	table <- format(table,justify="right",width=6);
	write.table(table,outName,quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE);
}

##
 # FORMAT AND WRITE DATA INTO OUTPUT .PLU AQUACROP FILE
 ###############################################################################
formatToPLUFile <- function(data,fileHead,path)
{
# I do not know yet if aquacrop is dealing with future date
# I hope so :)

# make one table from all the data
	table <- array(as.numeric(data$ppt),dim=dim(data$ppt));	# precip titled		'Rain'

# write it into a .PLU file
	outName <- paste(path$output,strsplit(path$file$temp,"\\.")[[1]][1],".PLU",sep="");
	if(!file.exists(path$output)){	# create output dir if does not exists
		dir.create(path$output, showWarnings = FALSE, recursive = FALSE, mode = "0777");
	}else{				# remove output file if does exists
		if(file.exists(outName))	file.remove(outName);
	}
	# head
	file.copy(	"./aquaCropTemplates/aquaCropTemplate.PLU",	outName,overwrite=TRUE);
	changeVar(	"DESCRIPTION",	fileHead$station$id,		outName,outName);
	changeVar(	"FD",		format(fileHead$period$start,"%d"),	outName,outName);
	changeVar(	"FM",		format(fileHead$period$start,"%m"),	outName,outName);
	changeVar(	"FY",		format(fileHead$period$start,"%Y"),	outName,outName);
	# body
	table <- format(table,justify="right",width=6);
	write.table(table,outName,quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE);
}

##
 # FORMAT AND WRITE DATA INTO OUTPUT .TMP AQUACROP FILE
 ###############################################################################
formatToEToFile <- function(data,fileHead,path)
{
# I do not know yet if aquacrop is dealing with future date
# I hope so :)

# make one table from all the data
	table <- array(as.numeric(data$ETo),dim=dim(data$ETo));
	table <- format(table,digits=3);

# write it into a .ETo file
	outName <- paste(path$output,strsplit(path$file$temp,"\\.")[[1]][1],".ETo",sep="");
	if(!file.exists(path$output)){	# create output dir if does not exists
		dir.create(path$output, showWarnings = FALSE, recursive = FALSE, mode = "0777");
	}else{				# remove output file if does exists
		if(file.exists(outName))	file.remove(outName);
	}
	# head
	file.copy(	"./aquaCropTemplates/aquaCropTemplate.ETo",	outName,overwrite=TRUE);
	changeVar(	"DESCRIPTION",	paste(fileHead$station$id," - ",data$AI," (UNCDD): ",sep=""),outName,outName);
	changeVar(	"FD",		format(fileHead$period$start,"%d"),	outName,outName);
	changeVar(	"FM",		format(fileHead$period$start,"%m"),	outName,outName);
	changeVar(	"FY",		format(fileHead$period$start,"%Y"),	outName,outName);
	# body
	table <- format(table,justify="right",width=6);
	write.table(table,outName,quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE);

	## to handle the linux/windows format issue
	if(Sys.info()["sysname"]=="Linux"){
		tempFile <- readLines(outName,n=-1,warn=FALSE);
		writeLines(tempFile,outName,sep="\r\n");
	}
}

