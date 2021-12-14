####
#### R script to plot heterozygosity and IBD segments across the whole genome
#### (for ConGen 2018)
####
 
###########################################################
# preliminaries: fill these out before running the script
###########################################################

minLeng <- 5                                           # minimum length of IBD segments ( in Mb) you will plot
plotChroms <- "1_RaGOO"                                  # which chromosome do you want to plot heterozygosity and IBD segments for?
plotInds <-1:10                                        # which of the individuals do you want to include in the plot

ibdDat <- read.table("conGen_100window_50minSNPs_10stepSize_9September2018_ibdTracts_100",header=FALSE)    # read in the IBD tracts
hetDat <- read.table("conGen_100window_50minSNPs_10stepSize_9September2018_hetMat_100",header=TRUE)        # read in the window heterozygosity values
windows  <- read.table("conGen_100window_50minSNPs_10stepSize_9September2018_windowInformationMatrix_100",header=TRUE)  # read in the information on the sliding windows

####################################
### get the ibd tracts
####################################

lengs <- ibdDat[,4] - ibdDat[,3]
ibdDat <- cbind(ibdDat,lengs)

#####################################
### get the heterozygosity data
#####################################

IDs <- colnames(hetDat[,1:ncol(hetDat)])

####################################
# get the window information
####################################

hetDat <- cbind(windows[,c(1,4,5)],hetDat)
physLengs <- table(windows$chrom)
print(paste0("Sum of lengths: ", sum(physLengs[plotChroms])))

########################################################################################
# make plot of heterozygosity and ibd segments across the chromosomes of 10 individuals
########################################################################################

cols <- c(rep(c("gray","darkred"),100),"darkgreen")     # colors of the chromosomes you will plot  
par(mfrow=c(length(plotInds),1),mar=c(0.5,7,0.5,2) + 0.1,xpd=TRUE)

ids <- colnames(hetDat)[plotInds+3]   # ids of the individuals

for(i in 1:length(plotInds)){
  start <- 0
  	print(paste0("Sum of lengths: ", sum(physLengs[plotChroms])))
 	plot(c(0,1),c(-0.1,1),type="n",axes=TRUE,xlim=c(0,sum(physLengs[plotChroms])),xlab="",ylab="")
	axis(side=2,at=seq(0,1,1))
 	text(x=-0.15*sum(physLengs[plotChroms]),y=0.4,labels=ids[i])

 	indDat <- NULL
 	indDat <- hetDat[,c(1:3,which(colnames(hetDat) == ids[i]))]

	indIbd <- NULL
  	indIbd <- ibdDat[which(ibdDat[,1] == ids[i]),]

 	for(j in plotChroms){
  		thisIndHet <- indDat[which(indDat[,1] == j),4]
  		thesePos <- rowMeans(indDat[which(indDat[,1]  == j),2:3])/1000000
		lines(thesePos+start,thisIndHet,lwd=0.8,col=cols[j])

		##### add the coordinates of identified ibd tracts

		chromIbd <- NULL
		chromIbd <- indIbd[which(indIbd[,2] == j & indIbd[,5]/1000000 > minLeng),]
   	if(nrow(chromIbd) > 0){
			for (z in 1:nrow(chromIbd)){
				lines(c(chromIbd[z,3]/1000000,chromIbd[z,4]/1000000)+start,y=c(-0.1,-0.1),lwd=0.8)
			}
		}
		start <- start + max(thesePos)
	}
}
