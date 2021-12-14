####
#### R script implementing the LOD method of Pemberton et al., (2012, Am J Hum Genet), 
#### Wang et al., (2009, Genet Epidemiol), Kardos et al. (2017, Genetics), Kardos et al. (2018, Nature Ecology & Evolution) to identify runs of homozygosity (for ConGen 2018)
####
####

##########################################################
# preliminaries: fill these out before running the script
##########################################################

winSize <- 100                                                # number of SNPS to use in a single window for IBD LOD score
minSNP <- 50                                                 # minimum umber of non-missing SNPs to include a window LOD score 
stepSize <- 10                                                # the step size (in SNPs) you want to use for sliding the window of size winSize across the genome
e <- 0.02                                                     # assumed genotyping error/mutation rate (this is the rate of hets expected within IBD segments)

outName <- paste("conGen","_",winSize,"window_",minSNP,"minSNPs_",stepSize,"stepSize_9September2018",sep="")   # name for the output files

# data <- read.table("VCF_sorted.ragoo.012", header=TRUE)    # read in the genotypes
data <- read.table("VCF_sorted.ragoo.beagle.012", header=TRUE)    # read in the genotypes
# data <- read.table("../Wolf_scripts/genotypes_conGenROH.012", header=TRUE)

library(dplyr)
my_freqs <- as.data.frame(table(data["chrom"]))
names(my_freqs)[1] <- "Chrom"
names(my_freqs)[2] <- "Freq"
good <- my_freqs %>% filter(Freq >= winSize, !(Chrom %in% c("Chr0_RaGOO", "Un_RaGOO")))

data <- data %>% filter(chrom %in% good$Chrom)
chrs <- unique(as.character(data[,1]))                        # vector of the chromosome names

if (length(chrs) == 0) {
   print("Something went wrong, no chromosome left.")
   quit()
   } else if ("Chr0_RaGOO" %in% data$chrom) {
   print("Data not filtered.")
   quit()
}

############################################################
#  extract the allele frequencies
############################################################

freqs <- as.numeric(data[,4])       # pull out the allele frequencies
data <- data[,-4]

NumLoci <- nrow(data)

print(paste0("Number of loci:", NumLoci))    # print the number of loci in the data set
print(paste0("Number of chromosomes:", length(chrs)))

##########################################################################################################
# build vectors of genotype probabilities for all possible genotypes at each locus
# assuming non-IBD and assuming IBDStatus... equations are from Wang et al. (2009, Genetic Epidemiology)
##########################################################################################################

# note that "a" stands for the alternative allele and "r" stands for reference... Entries in the 012 file give the number of alternative alleles

aaAuto <- NULL # genotype probs assuming autozygosity
arAuto <- NULL
rrAuto <- NULL
missAuto <- NULL

aaNAuto <- NULL # genotype probs assuming non-autozygosity
arNAuto <- NULL
rrNAuto <- NULL
missNAuto <- 1

aaAuto <- (1-e)*freqs+e*(freqs^2)
arAuto <- 2*e*freqs*(1-freqs)
rrAuto <- (1-e)*(1-freqs)+e*((1-freqs)^2)
missAuto <- 1

aaNAuto <- freqs^2
arNAuto <- 2*(freqs)*(1-freqs)
rrNAuto <- (1-freqs)^2
missNAuto <- 1

#### for each individual, build two vectors of genotype probilities: one assuming autozygosity, and one assuming non-autozyosity

ids <- colnames(data)[4:ncol(data)]

autoProbs <- NULL     # object to store the autozygous genotype probs
nAutoProbs <- NULL    # object to store the non-autozygous probs

for(i in 1:length(ids)){
   thisInd <- NULL
   thisInd <- data[,i+3]
   
   thisAutoProbs <- rep(NA,length(thisInd))
   thisNAutoProbs <- rep(NA,length(thisInd))
   
   thisAutoProbs[which(thisInd == 0)] <- rrAuto[which(thisInd == 0)]
   thisAutoProbs[which(thisInd == 1)] <- arAuto[which(thisInd == 1)]
   thisAutoProbs[which(thisInd == 2)] <- aaAuto[which(thisInd == 2)]
   thisAutoProbs[which(thisInd == -1)] <- missAuto

   thisNAutoProbs[which(thisInd == 0)] <- rrNAuto[which(thisInd == 0)]
   thisNAutoProbs[which(thisInd == 1)] <- arNAuto[which(thisInd == 1)]
   thisNAutoProbs[which(thisInd == 2)] <- aaNAuto[which(thisInd == 2)]
   thisNAutoProbs[which(thisInd == -1)] <- missNAuto

   autoProbs <- cbind(autoProbs,thisAutoProbs)
   nAutoProbs <- cbind(nAutoProbs,thisNAutoProbs)
   
   # print(i)
}

print("Initialization finished.")

############################################################################################################
# split the gemone into windows of SNPs and calculate the LOD score for each window for each individual
############################################################################################################

#### get the number of SNPs on each chromosome

nSNPs <- rep(NA,length(chrs))
for (i in 1:length(nSNPs)){
	nSNPs[i] <- nrow(data[which(data[,1] == chrs[i]),])  # nrow = number of rows. Basically, assign to the vector the number of SNPs found in that chromosome.
}

snpMat <- NULL
snpMat <- cbind(chrs,nSNPs)  # cdbind <- combine the columns

winVec <- NULL            # vector to store the ID of the IBD window for each snp

chromVec <- NULL          # vector to store the chromosome ID for each window
startVec <- NULL          # vectors to store the ID of the start and end SNPs of the IBD windows
endVec <- NULL            # vector to store the number of missing genotypes for each 60 SNP window
startPosVec <- NULL       # vector to store the bp position of the start of each window
endPosVec <- NULL         # vector to store the bp position of the end of each window

missMat <- NULL           # matrix to store the number of missing genotypes in each window for each individual

winIterator <- 0          # interator for windows

chrVec <- as.character(data[,1])   # vector of chromosome names

LODMat <- NULL            # matrix to store the LOD scores for each window for each individual
HetMat <- NULL            # matrix to store window heterozygosity
nSNPMat <- NULL           # matrix to store the number of typed SNPs for each individual in each window

for(i in 1:length(ids)){
  # print(paste0("Doing individual ", i, " ", ids[i]));

  thisGenoVec <- NULL
  thisGenoVec <- data[,i+3]   # gneotypes at all loci for the ith individual
  
  thisAutoProbs<- NULL
  thisAutoProbs <- autoProbs[,i]
  
  thisNAutoProbs<- NULL
  thisNAutoProbs <- nAutoProbs[,i]

  LODVec <- NULL
  HetVec <- NULL
  nSNPVec <- NULL
  thisMissVec <- NULL
  chrom_vec <- NULL
  for(k in 1:length(chrs)){
    # print(paste0("Doing chromosome ", k, " for individual ", i))
    chrGenos <- NULL
    chrAutoProbs <- NULL
    chrNAutoProbs <- NULL
    starts <- NULL
    chrSNPPos <- NULL
    startPos <- NULL 
    endPos <- NULL
    chrGenos <- thisGenoVec[which(chrVec == chrs[k])]   # genotypes at loci on the kth chromosome
    chrAutoProbs <- thisAutoProbs[which(chrVec == chrs[k])]
    chrNAutoProbs <- thisNAutoProbs[which(chrVec == chrs[k])]
    # print(paste0(length(chrAutoProbs), chrs[k], k)
    starts <- seq(1, length(chrAutoProbs), stepSize)     # start SNP of each IBD window in the genome
    starts <- starts[-which(starts > (length(chrAutoProbs)-winSize+1))]   # trim the end of starts so that there are no windows with fewer then winSize SNPs
    ends <- starts + winSize - 1
    if(sum(ends > length(chrAutoProbs)) >0){ends[which(ends > length(chrAutoProbs))] <- length(chrAutoProbs)}
    chrSNPPos <- data[which(as.character(data[,1]) == chrs[k]),3]
    startPos <- chrSNPPos[starts]
    endPos <- chrSNPPos[ends]
          
    #-------------------------------
    # save the window information
    #-------------------------------
    if(i == 1){
		chromVec <- append(chromVec,rep(chrs[k],length(starts)))
		startVec <- append(startVec,starts)
		endVec <- append(endVec,ends)
	   	startPosVec <- append(startPosVec,startPos)
		endPosVec <- append(endPosVec,endPos)

		theseWins <- NULL
    
		for (z in 1:length(starts)){
	  	  theseWins <- append(theseWins,rep(z,winSize))
	        }

		if(length(theseWins) > length(chrAutoProbs)){theseWins <- theseWins[1:length(chrAutoProbs)]}
		theseWins <- theseWins + winIterator  
		winVec <- append(winVec,theseWins)
		winIterator <- winIterator + k
		}     
    #--------------------------------------------------------
    # get the LOD scores and heterozygosity with each window
    #--------------------------------------------------------
    chr_chrom_vec <- rep(chrs[k], length(starts))
    chrom_vec <- append(chrom_vec, chr_chrom_vec)
    chrLODVec <- rep(NA,length(starts))   # vector to store the LOD score for each window
    chrHetVec <- rep(NA,length(starts))   # vector to store proportion het snps within each window
    chrNSNPVec <- rep(NA,length(starts));

    # This causes a bug
    # if (length(starts) == 0) {
    #     next;
    # }

    for (j in 1:length(starts)){  	
            tryCatch(thisMissVec <- append(thisMissVec,sum(chrGenos[starts[j]:ends[j]] == -1)),
	             error = function(err) { print(paste0("Error for index ", j, " for chromosome ", chrs[k], " for individual ", ids[i], " (no. ", i, "; total length of starts: ", length(starts), ")")) })
	    chrNSNPVec[j] <- sum(chrGenos[starts[j]:ends[j]] != -1)
	    if(length(which(chrGenos[starts[j]:ends[j]] != -1)) >= minSNP){
			# get LOD score
			aProbs <- NULL
			aProbs <- chrAutoProbs[starts[j]:ends[j]]
      
			naProbs <- NULL
			naProbs <- chrNAutoProbs[starts[j]:ends[j]]
      
			chrLODVec [j] <- sum( log10( aProbs/naProbs),na.rm=TRUE)
	    
			### get proportion of SNPs that are het in this window
			hets <- NULL  # number of het snps
			# POTENTIAL BUG
			loci <- # number of loci not missing in this window
	   
			hets <- length(which(chrGenos[starts[j]:ends[j]] == 1))
			loci <- length(which(chrGenos[starts[j]:ends[j]] != -1))
	   
			chrHetVec[j] <- hets/loci
			}    
	  }
    LODVec <- append(LODVec,chrLODVec)
    HetVec <- append(HetVec,chrHetVec)
    nSNPVec <- append(nSNPVec,chrNSNPVec)
    # print(paste("chr",k," done",sep=""))
  }
  LODMat <- cbind(LODMat, LODVec)
  HetMat <- cbind(HetMat,HetVec)
  nSNPMat <- cbind(nSNPMat,nSNPVec)
  # print(paste("done with individual ",i,sep=""))
}

print("Finished parsing individuals")

###################################################################################################################################
# Identify the LOD score threshold above which you will call a segment IBD using gaussian kernel density model...
# This is the local mimimum density between the two modes in the distribution of LOD scores... This follows Pemberton et al., 2012
###################################################################################################################################

lodVec <- NULL
for(i in 1:ncol(LODMat)){
	lodVec <- append(lodVec,LODMat[,i])
}
	
lodVec2 <- lodVec
# Remove NAs
if(sum(is.na(lodVec)) > 0){lodVec2 <- lodVec[-which(is.na(lodVec) == TRUE)]}

dens <- density(lodVec2,kernel="gaussian",na.rm=TRUE)
plot(dens)

plot(dens[[1]][200:350],dens[[2]][200:350])
thresh <- dens[[1]][200:350][which(dens[[2]][200:350] == min(dens[[2]][200:350]))] # 

IBDMat <- LODMat >= thresh
IBDMat <- as.data.frame(IBDMat)

winMat <- cbind(chromVec,startVec,endVec,startPosVec,endPosVec)    # matrix storing information on the order and location of each IBD window in the genome
winMat <- as.data.frame(winMat)

####################################################################################################################
# concatenate cointiguous IBD chromosome segments and save the locations of IBD tracts for each individual 
####################################################################################################################

IBDTractMat <- NULL    # store the start, stop positions, and the length of each RoH for each individudal

for (i in 1:ncol(IBDMat)){
  thisDat <- NULL       # the ith individual's IBD statuses in each window across the genome
  thisDat <- IBDMat[,i]
  
  outDat <- NULL  # matrix to store the IBD tract information for this individual
  
  for(j in 1:length(chrs)){
    chrDat <- NULL
    chrDat <- thisDat[which(as.character(winMat$chromVec) == chrs[j])]
    
    if(sum(is.na(chrDat)) > 0){chrDat[is.na(chrDat)] <- FALSE}

    startDat <- NULL                                                   # start positions for windows on this chromosome
    startDat <- winMat[which(winMat$chromVec == chrs[j]),]
      
    begins <- NULL          # the first windows in each RoH
    ends <- NULL            # the last window in each RoH
    

    # which windows are TRUE (for IBD) and the previous one is FALSE (for IBD)
    begins <- tryCatch(which(chrDat[2:length(chrDat)] == TRUE & chrDat[1:(length(chrDat)-1)] == FALSE )+1,
                       error = function(err) {
		       print("Begins not assigned correctly")
		       })
    if(is.na(chrDat[1]) == FALSE) {
      if(chrDat[1] == TRUE){
        begins <- c(1,begins)
      }
    }
    
    ends <- which(chrDat[1:(length(chrDat)-1)] == TRUE & chrDat[2:length(chrDat)] == FALSE ) # which windows are true for IBD and the next one is FALSE (for IBD)

    if(is.na(chrDat[length(chrDat)]) == FALSE) {
      if(chrDat[length(chrDat)] == TRUE){
        ends <- c(ends,length(chrDat))
      }
    }
     
    tractStartPos <- NULL
    tractEndPos <- NULL
    
    tractStartPos <- as.character(startDat$startPosVec [begins])
    tractEndPos <- as.character(startDat$endPosVec[ends])
    
    #-------------------------------------------------------------------------------------
    # if there are multiple ROH, see if any of them overlap and concatenate them
    #-------------------------------------------------------------------------------------
    
    newStartPos <- NULL
    newEndPos <- NULL
	
    if(length(tractStartPos) > 1){
	    overs <- NULL   # test if each RoH window overlaps with the previous one
	    overs <- as.numeric(tractStartPos[2:length(tractStartPos)]) <= as.numeric(tractEndPos[1:(length(tractStartPos)-1)])
	    overs <- c(NA,overs)
	
	    groups <- c(1,rep(NA,length(overs)-1))       # vector of the identity of unique RoH
	
   	  for(k in 2:length(groups)){
		    if(overs[k] == FALSE){groups[k] <- 1 + groups[k-1]}else
		    if(overs[k] == TRUE) {groups[k] <- groups[k-1]}
		  }

    	uniqGroups <- NULL  # unique RoH identifiers
   	  uniqGroups <- unique(groups)
	
   	  for(k in 1:length(uniqGroups)){
		    newStartPos <- append(newStartPos, tractStartPos[which(groups == uniqGroups[k])[1]])
		    newEndPos   <- append(newEndPos,tractEndPos[which(groups == uniqGroups[k])[length(which(groups == uniqGroups[k]))]])
		  }
    }
 
    if(length(tractStartPos) == 1){
      newStartPos <- NULL
      newEndPos <- NULL
      newStartPos <- tractStartPos
      newEndPos <- tractEndPos
    }

    outDat <- rbind(outDat,cbind(rep(ids[i],length(newStartPos)),rep(chrs[j],length(newStartPos)),newStartPos,newEndPos))
    # print(j)
  }
  IBDTractMat <- rbind(IBDTractMat,outDat)
}

lengths <- as.numeric(IBDTractMat[,4]) - as.numeric(IBDTractMat[,3])

colnames(HetMat) <- ids

#########################################################
# save files
#########################################################


write.table(IBDMat,file=paste(outName,"_IBDStatusMatrix_",winSize,sep=""),quote=FALSE,row.names=FALSE)
write.table(winMat,file=paste(outName,"_windowInformationMatrix_",winSize,sep=""),quote=FALSE,row.names=FALSE)
write.table(IBDTractMat,file=paste(outName,"_ibdTracts_",winSize,sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(LODMat,file=paste(outName,"_windowLODScoreMatrix_",winSize,sep=""),quote=FALSE,row.names=FALSE)
write.table(nSNPMat,file=paste(outName,"_nSNPMatrix_",winSize,sep=""),quote=FALSE,row.names=FALSE)
write.table(HetMat,file=paste(outName,"_hetMat_",winSize,sep=""),quote=FALSE,row.names=FALSE)
