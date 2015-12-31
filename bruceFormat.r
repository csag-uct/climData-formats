##
 # FILE bruceFormat.r
 # AUTHOR olivier crespo
 # https://r-forge.r-project.org/projects/xpos-r/convertBruceFormat
 ###############################################################################
 #
 # functions required to read, and work on bruce data format
 #
 ###############################################################################

## 
 # HEAD FILES FORMAT
 ###############################################################################
read_bruceHeadFile <- function(path2file)
{
	## first line
	temp <- scan(path2file,what="character",sep=" ",nlines=1,quiet=TRUE);
	temp <- temp[temp!=""];
	station <- list("id"=temp[1],"lat"=as.numeric(temp[2]),"lon"=as.numeric(temp[3]),"alt"=as.numeric(temp[4]));

	## second line
	temp <- scan(path2file,what="character",sep=" ",skip=1,nlines=1,quiet=TRUE);
	temp <- temp[temp!=""];
	ss <- strsplit(temp,"");
	period <- list("start"=as.Date(temp[1],"%Y%m%d"),"end"=as.Date(temp[2],"%Y%m%d"),"type"=as.numeric(temp[3]),"startH"=paste(ss[[1]][9],ss[[1]][10],sep=""),"endH"=paste(ss[[2]][9],ss[[2]][10],sep=""));
	## for info
	# year	=	format(date,"%Y")
	# month	=	format(date,"%m")
	# day	=	format(date,"%d")
	
	## noticed that for downscaled hindacst only so far...
	# because end date is 30-dec while it still adds one day and actually finishes on 31-dec
	if(period$type==2 && format(period$end,"%m")==12 && format(period$end,"%d")==30){
		period$end<-period$end+1;
	}
	
	## third line
	temp <- readLines(path2file,n=3);
	comm <- temp[3];

return(list("station"=station,"period"=period,"comment"=comm));
}

## 
 # IMPORT DATA according to common period
 ###############################################################################
importData <- function(path,fileHead)
{
# import all
	# tmn
	path2file <- paste(path$input,path$folder$tmn,path$file$temp,sep="");
	head_tmn <- read_bruceHeadFile(path2file);
	tmn <- as.array(as.numeric(scan(path2file,what="numeric",sep="\n",skip=3,nlines=-1,quiet=TRUE)));
	# tmx
	path2file <- paste(path$input,path$folder$tmx,path$file$temp,sep="");
	head_tmx <- read_bruceHeadFile(path2file);
	tmx <- as.array(as.numeric(scan(path2file,what="numeric",sep="\n",skip=3,nlines=-1,quiet=TRUE)));
	# ppt
	path2file <- paste(path$input,path$folder$ppt,path$file$prec,sep="");
	head_ppt <- read_bruceHeadFile(path2file);
	ppt <- as.array(as.numeric(scan(path2file,what="numeric",sep="\n",skip=3,nlines=-1,quiet=TRUE)));

# cut the tail of files that do not belong into the common period
	# tmn
	tail <- 0; tail <- head_tmn$period$end-fileHead$period$end;
	if(tail>0) tmn <- tmn[1:(dim(tmn)[1]-tail)];
	# tmx
	tail <- 0; tail <- head_tmx$period$end-fileHead$period$end;
	if(tail>0) tmx <- tmx[1:(dim(tmx)[1]-tail)];
	# ppt
	tail <- 0; tail <- head_ppt$period$end-fileHead$period$end;
	if(tail>0) ppt <- ppt[1:(dim(ppt)[1]-tail)];

# cut the head of files that do not belong into the common period
	# tmn
	head <- 0; head <- fileHead$period$start-head_tmn$period$start;
	if(head>0) tmn <- tmn[(head+1):dim(tmn)[1]];
	# tmx
	head <- 0; head <- fileHead$period$start-head_tmx$period$start;
	if(head>0) tmx <- tmx[(head+1):dim(tmx)[1]];
	# ppt
	head <- 0; head <- fileHead$period$start-head_ppt$period$start;
	if(head>0) ppt <- ppt[(head+1):dim(ppt)[1]];

# check period consistencies
	if( dim(tmn)!=dim(tmx) || dim(tmn)!=dim(ppt) || dim(tmx)!=dim(ppt)){
		print("# ERROR: not your fault, but the importData function (in bruceFormat.r) is failing producing files of same length :(",quote=FALSE);
		stop();
	}

return(list("tmn"=tmn,"tmx"=tmx,"ppt"=ppt));
}

##
 # FORMAT AND WRITE DATA INTO OUTPUT .eto CSAG LIKE FILE
 ###############################################################################
formatToEToCSAG <- function(data,head,path)
{
# make one table from all the data
	table <- array(as.numeric(data$ETo),dim=dim(data$ETo));
	table <- format(table,digits=1,nsmall=2);

# write it into a .eto file
	outName <- paste(path$output,strsplit(path$file$temp,"\\.")[[1]][1],".eto",sep="");
	if(!file.exists(path$output)){	# create output dir if does not exists
		dir.create(path$output, showWarnings = FALSE, recursive = FALSE, mode = "0777");
	}else{				# remove output file if does exists
		if(file.exists(outName))	file.remove(outName);
	}
	# head
	file.copy(	"./csagTemplates/csagTemplate",	outName,overwrite=TRUE);
	changeVar(	"ID",		head$station$id,outName,outName);
	changeVar(	"LAT",		format(head$station$lat,nsmall=2),outName,outName);
	changeVar(	"LON",		format(head$station$lon,nsmall=2),outName,outName);
	changeVar(	"ALT",		format(head$station$alt,nsmall=2),outName,outName);
	changeVar(	"START",	paste(format(head$period$start,"%Y%m%d"),head$period$startH,sep=""),outName,outName);
	changeVar(	"END",		paste(format(head$period$end,"%Y%m%d"),head$period$endH,sep=""),outName,outName);
	changeVar(	"TYPE",		head$period$type,outName,outName);
	changeVar(	"DESCRIPTION",	"ETo_ESTIMATION(TMN,TMX,PPT)",outName,outName);
	# body
	write.table(table,outName,quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE);
}

##
 # FORMAT AND WRITE DATA INTO OUTPUT .rdn CSAG LIKE FILE
 ###############################################################################
formatToRdnCSAG <- function(data,head,path)
{
# make one table from all the data
	table <- array(as.numeric(data$sRad),dim=dim(data$sRad));
	table <- format(table,digits=1,nsmall=2);

# write it into a .rdn file
	outName <- paste(path$output,strsplit(path$file$temp,"\\.")[[1]][1],".rdn",sep="");
	if(!file.exists(path$output)){	# create output dir if does not exists
		dir.create(path$output, showWarnings = FALSE, recursive = FALSE, mode = "0777");
	}else{				# remove output file if does exists
		if(file.exists(outName))	file.remove(outName);
	}
	# head
	file.copy(	"./csagTemplates/csagTemplate",	outName,overwrite=TRUE);
	changeVar(	"ID",		head$station$id,outName,outName);
	changeVar(	"LAT",		format(head$station$lat,nsmall=2),outName,outName);
	changeVar(	"LON",		format(head$station$lon,nsmall=2),outName,outName);
	changeVar(	"ALT",		format(head$station$alt,nsmall=2),outName,outName);
	changeVar(	"START",	paste(format(head$period$start,"%Y%m%d"),head$period$startH,sep=""),outName,outName);
	changeVar(	"END",		paste(format(head$period$end,"%Y%m%d"),head$period$endH,sep=""),outName,outName);
	changeVar(	"TYPE",		head$period$type,outName,outName);
	changeVar(	"DESCRIPTION",	"Rdn_ESTIMATION(TMN,TMX,PPT)",outName,outName);
	# body
	write.table(table,outName,quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE);
}

##
 # FOR CHECK ONLY
 ###############################################################################
formatToTmnCSAG <- function(data,head,path)
{
# make one table from all the data
	table <- array(as.numeric(data$tmn),dim=dim(data$tmn));
	table <- format(table,digits=1,nsmall=2);

# write it into a .eto file
	outName <- paste(path$output,strsplit(path$file$temp,"\\.")[[1]][1],".tmn.check",sep="");
	if(!file.exists(path$output)){	# create output dir if does not exists
		dir.create(path$output, showWarnings = FALSE, recursive = FALSE, mode = "0777");
	}else{				# remove output file if does exists
		if(file.exists(outName))	file.remove(outName);
	}
	# head
	file.copy(	"./csagTemplates/csagTemplate",	outName,overwrite=TRUE);
	changeVar(	"ID",		head$station$id,outName,outName);
	changeVar(	"LAT",		format(head$station$lat,nsmall=2),outName,outName);
	changeVar(	"LON",		format(head$station$lon,nsmall=2),outName,outName);
	changeVar(	"ALT",		format(head$station$alt,nsmall=2),outName,outName);
	changeVar(	"START",	paste(format(head$period$start,"%Y%m%d"),head$period$startH,sep=""),outName,outName);
	changeVar(	"END",		paste(format(head$period$end,"%Y%m%d"),head$period$endH,sep=""),outName,outName);
	changeVar(	"TYPE",		head$period$type,outName,outName);
	changeVar(	"DESCRIPTION",	"TMN_FOR_CHECK_ONLY",outName,outName);
	# body
	write.table(table,outName,quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE);
}

