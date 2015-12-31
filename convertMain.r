##
 # FILE convertMain.r
 # AUTHOR olivier crespo
 # https://r-forge.r-project.org/projects/xpos-r/convertBruceFormat
 ###############################################################################
 #
 # main algorithm, used for any fromat
 # - convertOne
 #	intends to convert one data set (one GCM - one period) for one model
 # - convert
 #	is here to process convertOne
 #	through lots of data sets and so far APSIM and AQUACROP crop models
 #
 ###############################################################################

## # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ###
##									      ##
##									      ##
##              USER JOB: FILL IN AND ADAPT ALL init_* FUNCTIONS	      ##
##									      ##
##									      ##
## # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ###

loop_on_paths <- function()
{
	for (it in 1:10){
		convert("ap",allGCM=T,iteration=it);
	}
}

##
 # INITIALISE YOUR DATA PATHS AND NAMES
 ###############################################################################
 # separated folders with "/" even for windows
 # finish all paths and folder name with "/"
 ###############################################################################
init_paths <- function(it)
{
	# in which folder to read the data
	input <- "/home/crespo/Desktop/Optimisation/xpos-r/ClimateDataTools/Input";
#	it<-1;
	# in which folder to write out the data
	if (is.null(it)){
		output <- "/home/crespo/Desktop/Optimisation/xpos-r/ClimateDataTools/Output";
	}else{
		output <- paste("/home/crespo/Desktop/11_START/ApsimMetFiles/rep",it,"/",sep="");
	}

	# what are the name of the data folders
#	folder <- list	(	"tmn"=	paste("tmin-",it,"/",sep=""),	# folder name for minimal temperatures
#				"tmx"=	paste("tmax-",it,"/",sep=""),	# folder name for maximal temperatures
#				"ppt"=	paste("ppt-",it,"/",sep="")	# folder name for precipitation
#			);	
	folder <- list	(	"tmn"=	"tmin/",	# folder name for minimal temperatures
				"tmx"=	"tmax/",	# folder name for maximal temperatures
				"ppt"=	"ppt/"	# folder name for precipitation
			);	

return(list("input"=input,"output"=output,"folder"=folder));
}

## # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ###
##									      ##
##	           PLAY FREELY WITH WHATEVER IS ABOVE ...		      ##
##       BUT								      ##
##	           PLAY CAREFULLY WITH WHATEVER IS BELOW		      ##
##									      ##
## # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ###

##
 # CONVERSION ROUTINE FOR MULTIPLE GCMs
 ###############################################################################
 # > model
 #	"br" for bruce format like
 # 	"ap" for APSIM
 #	"aq" for AQUACROP
 #	"cs" for CSAG (basically you want to compute ETo)
 #	"ds" for DSSAT
 #	"ag" for AgMIP
 #	"all" for both
 # > allGCM
 #	FALSE (default) plays only with stat list (defined in init_stat)
 #	TRUE plays also with GCMs list (defined in init_GCMs) 
 # > requires
 # 	init_paths and init_stat anyway,
 #	+ init_GCMs if allGCM is TRUE
 ###############################################################################
convert <- function(model,inLand=TRUE,seeSteps=FALSE,fillIn=TRUE,iteration=NULL)	# ,allSRES=TRUE,allGCM=TRUE
{
### crop models
	if(model=="all") model <- 1;
	if(model=="ap") model <- 2;
	if(model=="aq") model <- 3;
	if(model=="cs") model <- 4;
	if(model=="ds") model <- 5;
	if(model=="ag") model <- 6;
	if (!is.numeric(model)){
		print("### ERROR: unknown target model");
		print("### targetModel available so far: 'ag' (AgMIP), 'cs' (CSAG), 'ap' (APSIM), 'aq' (AQUACROP), 'cs' (CSAG), 'ds' (DSSAT) or 'all'");
		stop();
	}

### initialisation
	path <- init_paths(iteration)
	parentFolder <- path$input
	rcp_t <- list.files(parentFolder)
	for (r in 1:length(rcp_t)){
		print(paste("   ",rcp_t[r],sep=" > "),quote=FALSE)
		RCPFolder <- paste(parentFolder,rcp_t[r],sep="/")
		tPe_t <- list.files(RCPFolder)
		for (t in 1:length(tPe_t)){
			print(paste("   ",rcp_t[r],tPe_t[t],sep=" > "),quote=FALSE)
			pptFolder <- paste(RCPFolder,tPe_t[t],"ppt",sep="/")
			sta_t <- list.files(pptFolder)
			for (s in 1:length(sta_t)){
				print(paste("   ",rcp_t[r],tPe_t[t],sta_t[s],sep=" > "),quote=FALSE);
				pathToStation <-	list(	"input"=paste(path$input,rcp_t[r],tPe_t[t],"",sep="/"),
								"output"=paste(path$output,rcp_t[r],tPe_t[t],"",sep="/"),
								"folder"=	list(	"tmn"=path$folder$tmn,
											"tmx"=path$folder$tmx,
											"ppt"=path$folder$ppt
										),
								"file"=		list(	"temp"=sta_t[s],
											"prec"=sta_t[s]
										),
								"inland"=inLand,
								"arid"='A'
							);
				switch(model,
					{	# all
						convertOne("ap",pathToStation,seeSteps,fillIn);
						convertOne("aq",pathToStation,seeSteps,fillIn);
#						convertOne("cs",pathToStation,seeSteps,fillIn);
						convertOne("ds",pathToStation,seeSteps,fillIn);
#						convertOne("ag",pathToStation,seeSteps,fillIn);
					},{	# apsim only
						convertOne("ap",pathToStation,seeSteps,fillIn);
					},{	# aquacrop only
						convertOne("aq",pathToStation,seeSteps,fillIn);
					},{	# csag like only
						convertOne("cs",pathToStation,seeSteps,fillIn);
					},{	# dssat only
						convertOne("ds",pathToStation,seeSteps,fillIn);
					},{	# AgMIP only
						convertOne("ag",pathToStation,seeSteps,fillIn);
					}
				);

			}
		}
	}
print(" ... process completed ...");
}

##
 # CONVERT MAIN FUNCTION
 ###############################################################################
 # is the function to be called for any conversion
 # leave pathToStation NULL for one station
 # ##############################################################################
convertOne <- function(targetModel,pathToStation=NULL,seeSteps,fillIn)
{
### sources
	source("/home/crespo/Desktop/Optimisation/xpos-r/ClimateDataTools/ClimFormats/checkFunctions.r");
	source("/home/crespo/Desktop/Optimisation/xpos-r/ClimateDataTools/ClimFormats/bruceFormat.r");
	source("/home/crespo/Desktop/Optimisation/xpos-r/ClimateDataTools/ClimFormats/metTransformations.r");
	source("/home/crespo/Desktop/Optimisation/xpos-r/ClimateDataTools/ClimFormats/agriParameters.r");
	source("/home/crespo/Desktop/Optimisation/xpos-r/ClimateDataTools/ClimFormats/agriParameters_loops.r");

	if(is.null(pathToStation)){
		print("### ERROR: no station specified !!",quote=FALSE);
		stop();
	}

#?? pathToStation$arid
### target models
	if(targetModel=="ap") targetModel <- 1;
	if(targetModel=="aq") targetModel <- 2;
	if(targetModel=="cs") targetModel <- 3;
	if(targetModel=="ds") targetModel <- 4;
	if(targetModel=="ag") targetModel <- 5;
	if (!is.numeric(targetModel)){
		print("### ERROR: unknown target model",quote=FALSE);
		print("### targetModel available so far: 'ds' (DSSAT), 'cs' (CSAG), 'ap' (APSIM), 'aq' (AQUACROP) or 'all'",quote=FALSE);
		stop();
	}

## because so far, no crop model deal with missing data
 # 1: I am only interseted in the time period covered by tmn, tmx and ppt
	if(seeSteps)	print("... check unconsistencies ...",quote=FALSE);
	fileHead <- checkData(pathToStation);
	if(seeSteps)	print("... import data ...",quote=FALSE);
	data <- importData(pathToStation,fileHead);
## because so far, no crop model deal with missing data
 # 2: if there is any missing data (that I can spot) I will warn the user
	if(seeSteps)	print("... looking for missing data ...",quote=FALSE);
	data <- checkMissing(data);
	checkTmpRain(data);
## because so far, no crop model deal with missing data
 # 3: I need to make sure years are made of real number of days
	switch(fileHead$period$type+1,
		{	# 0 is for real
		},{	# 1 is for 365
			data <- transform_type1(data,fileHead,fillIn);
		},{	# 2 is for 360
			data <- transform_type2(data,fileHead,fillIn);
		}
	);

	if(seeSteps)	print("... create years and julian days ...",quote=FALSE);
	data <- createYearJulianDays(data,fileHead);

## then starts the crop model related operations
	switch(targetModel,
		{	#################### APSIM #
			source('/home/crespo/Desktop/Optimisation/xpos-r/ClimateDataTools/ClimFormats/convertToApsim.r');
			if(seeSteps)	print("... compute radiation ...",quote=FALSE);
			data <- compute_radn(data,fileHead$station,pathToStation$inland);
			if(seeSteps)	print("... compute tav and amp ...",quote=FALSE);
			data <- compute_tavNamp(data);
			if(seeSteps)	print("... format and write data into .met file ...",quote=FALSE);
			formatToMetFile(data,fileHead,pathToStation);
		},	# APSIM ####################
		{	################# AQUACROP #
			# aquacrop deals with day, 10-days and monthly records
			# so far we deal only with day records
			source('/home/crespo/Desktop/Optimisation/xpos-r/ClimateDataTools/ClimFormats/convertToAquacrop.r');
			if(seeSteps)	print("... compute ETo ...",quote=FALSE);
			if(is.numeric(pathToStation$arid)){
				data <- compute_ETo(data,fileHead,pathToStation$inland,pathToStation$arid);
			}else{
				data <- compute_ETo(data,fileHead,pathToStation$inland,3);
				counter <- 0;					# infinite loop check
				while(data$doItAgain!=0){
					counter <- counter+1;
					if( counter > 2){
						print("ERROR: infinite loop in convertOne function, my mistake",quote=FALSE);
						stop();
					}
					data <- compute_ETo(data,fileHead,pathToStation$inland,data$arid+data$doItAgain);
				}
			}
			if(seeSteps)	print("... format and write data into .TMP, .PLU and .ETo files ...",quote=FALSE);
			formatToTMPFile(data,fileHead,pathToStation);
			formatToPLUFile(data,fileHead,pathToStation);
			formatToEToFile(data,fileHead,pathToStation);
		},	# AQUACROP ####################
		{	######################## CSAG #
			if(seeSteps)	print("... compute ETo ...",quote=FALSE);
			if(is.numeric(pathToStation$arid)){
				data <- compute_ETo(data,fileHead,pathToStation$inland,pathToStation$arid);
			}else{
				data <- compute_ETo(data,fileHead,pathToStation$inland,3);
				counter <- 0;					# infinite loop check
				while(data$doItAgain!=0){
					counter <- counter+1;
					if( counter > 2){
						print("ERROR: infinite loop in convertOne function, my mistake",quote=FALSE);
						stop();
					}
					data <- compute_ETo(data,fileHead,pathToStation$inland,data$arid+data$doItAgain);
				}
			}
			if(seeSteps)	print("... remove added data ...",quote=FALSE);
			data <- removeAddedDays(data,fileHead);
			if(seeSteps)	print("... format and write data into .rdn and .eto files ...",quote=FALSE);
			formatToEToCSAG(data,fileHead,pathToStation);
			formatToRdnCSAG(data,fileHead,pathToStation);
			# formatToTmnCSAG(data,fileHead,pathToStation);
		},	# CSAG ########################
		{	####################### DSSAT #
			source('/home/crespo/Desktop/Optimisation/xpos-r/ClimateDataTools/ClimFormats/convertToDssat.r');
			if(seeSteps)	print("... compute radiation ...",quote=FALSE);
			data <- compute_radn(data,fileHead$station,pathToStation$inland);

			if(seeSteps)	print("... compute tav and amp ...",quote=FALSE);
			data <- compute_tavNamp(data);
			if(seeSteps)	print("... format and write data into .WTH file ...",quote=FALSE);
			formatToWTHfile(data,fileHead,pathToStation);
		}	# DSSAT #########################
	);

if(seeSteps)	print("... conversion completed ...",quote=FALSE);
}

