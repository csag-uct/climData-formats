##
 # FILE dataMain.r
 # AUTHOR olivier crespo
 # https://r-forge.r-project.org/projects/xpos-r/
 ###############################################################################

###############################################################################
# LISA, CUTTING THE 5 YEARS TAIL OF CMIP5 DOWNSCALED DATA
###############################################################################
# inFolder : input Folder in which to find GCM-RCP folders, in which are the tmin, tmax, ppt
inFo <- '/terra/data/downscaled/lcoop/dscl-pca/cmip5/output/malawi_usaid/bruce_format'
outFo <- '/terra/data/downscaled/lcoop/dscl-pca/cmip5/output/malawi_usaid/1965-2094_format'
# sDate
sD <- as.Date("1965-01-01")
eD <- as.Date("2094-12-31")
cmip5_cutPeriod <- function(inFolder=inFo, sDate=sD, eDate=eD, outFolder=outFo)
{
	rcp_l <- NULL
	rcp_t <- list.files(inFolder)

	for(r in 1:length(rcp_t)){
		# not for ncep
		if(r == grep("ncep",rcp_t)) next		
		sta_l <- NULL
		print(paste("    > ",rcp_t[r],sep=""),quote=F)
		tmpOut <- paste(outFolder,rcp_t[r],sep="/")	
		if(!file.exists(tmpOut))	dir.create(tmpOut,showWarnings=TRUE,recursive=FALSE,mode="0777")

		tmpIn <- paste(inFolder,rcp_t[r],sep="/")
		sta_t <- list.files(paste(tmpIn,"ppt",sep="/"))
		for(s in 1:length(sta_t)){
			print(paste("    >    > ",sta_t[s],sep=""),quote=F)
			# read it
			metD <- read_oldCSAGformat(tmpIn,sta_t[s])	# requires dataRead.r
			# cut it
			newD <- pert_period(metD,sDate,eDate)
			# write it
			write_oldCSAGformat(newD,tmpOut)
		}
	}
}


