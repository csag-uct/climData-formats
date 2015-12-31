#
#
#
#

source('R_readMetData.r')
source('R_writeMetData.r')

#
# TRANSLATING for wiltrud
#
all <- read.csv('~/Desktop/Wine-shared/Projects/2015-2016_AgMIP2/RZA-allFS/Fs_Quinarys.csv')
f_v <- as.character(all[,colnames(all)=="SUB_CAT"])
inF <- '~/Desktop/Wine-shared/Projects/2015-2016_AgMIP2/RZA-allFS/Obshis_wth'
outF <- '~/Desktop/Wine-shared/Projects/2015-2016_AgMIP2/RZA-allFS/Obshis_agmip'
tran_DSSAT2AgMIP <- function(file_v=f_v,inFolder=inF,outFolder=outF)
{
	for(f in file_v){
		print(paste(which(file_v==f),length(file_v),sep=" < "))

		if(nchar(f)<4)	f <- paste("0",f,sep="")
		tmp <- read_DSSATformat(paste(inFolder,paste(f,"wth",sep="."),sep="/"))

# need some tests
		#tmin > tmax
		if(any(tmp$data$tmin>tmp$data$tmax))	browser()		

		write_AgMIPformat(tmp,paste(outFolder,paste(paste(f,"0XXX",sep=""),"AgMIP",sep="."),sep="/"),checkMiss=FALSE)
	}
}
