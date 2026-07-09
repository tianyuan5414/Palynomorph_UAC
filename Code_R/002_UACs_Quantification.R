#This script aims ata quantifying UAC signals from corrected FTIR spectra in the
  #`001_Spectra treatment.R`

#Load necessary packages
{
  library(here)
  
  #Moving average calculation
  library(zoo)
  
  #Peak statistics 
  library(pavo)
  
  #filtered component
  library(signal)
  
  #function to pad and filter the data
  library(vegan)
}

#Load necessary functions
{
  #Ben's peak calculation functions
  {
    # Load peak detection script
    source("https://raw.githubusercontent.com/benbell95/peak-detection/main/r/peak_detection.r")
    # Load peak area script
    source("https://raw.githubusercontent.com/benbell95/peak-detection/main/r/peak_area.r")
  }
  
  bwFilter <- function(arrayTime, arrayData, interInterval, interType, cutoffFreq, filterType, filterOrder) {
    
    #Transfer the original timespace into a new timespace
    arrayTimeNew <- seq(from = arrayTime[1], to = arrayTime[length(arrayTime)], by = interInterval)
    
    #interpolate the original data into targeted timespace
    arrayDataNew <- approx(arrayTime, arrayData, arrayTimeNew, method = interType)[[2]]
    
    #Prescribe a Butterworth filter
    bwFil <- butter(filterOrder, cutoffFreq, type = filterType)
    
    #Apply filtering function on the raw data
    arrayDataNewFiltered <- filtfilt(bwFil, arrayDataNew)
    
    #Output the results
    return(c(arrayTimeNew, arrayDataNewFiltered))
  }
  
  my.bw <- function(Ages,Data){
    #pad the dataset with 500 yr of mean value 
    Ages <- round(Ages,0)
    Pad.Age <- c(rev(seq(min(Ages)-1, min(Ages) - 500, -1)), Ages, seq(max(Ages)+1, max(Ages)+500,1))
    Pad.Data <- c(rep(mean(Data),500),Data,rep(mean(Data), 500))
    #Filter the raw data with interval of 1
    bwf <- matrix(data = bwFilter(Pad.Age, Pad.Data, 1, 'linear', 1/100, 'low', 2),
                  nrow = length(bwFilter(Pad.Age, Pad.Data, 1, 'linear', 1/100, 'low', 2)) / 2,
                  ncol = 2, byrow = FALSE)
    bwf <- bwf[-c(1:500),]
    bwf <- head(bwf,-500)
    bw.df <- data.frame(Age = bwf[,1], Filt.Data = bwf[,2])
    
    #Output the filter results
    return(bw.df)
  }
}

#Prepare chronology information for the read spectra dataset
{
  #Chronology
  chronologyHM20 <- apply(as.matrix(read.csv(here('Data_R', 'chronologyBchron.csv'))), 
                          c(1,2), as.numeric)
  colnames(chronologyHM20)[1:2] <- c('Years', 'Depth')
  
  #Read transfering functions between IntCal20 and GICC05
  chronologyHM20TransferInte <- apply(read.csv(here('Data_R', 'chronologyHM20TransferInte.csv'),
                                               row.names = 1), c(1,2), as.numeric)
  
}

#Read smoothed and baseline corrected data, and apply clean up for the dataset
{
  #Sphagnum
  {
    #Read smoothed and baseline corrected spectra matrix
    rawSmooBCSpecMatSph <- read.csv(here('Output_R',
                                         'fingerPrintBCNorSph.csv'),
                                    header = TRUE,
                                    row.names = 1)
    
    rawSmooBCSpecMatSph <- apply(rawSmooBCSpecMatSph, c(1,2),
                                 as.numeric)
    
    #Rename the column
    colnames(rawSmooBCSpecMatSph) <- c('Batch', seq(from = 1800, to = 800,by = -4))
    
    signalMultipleApproMatScaSph <- matrix(nrow = nrow(rawSmooBCSpecMatSph),
                                        ncol = 2,
                                        dimnames = list(c(),
                                                        c('Batch',
                                                          'Peak_Area')))
    
    signalMultipleApproMatScaSph[, 1] <- rawSmooBCSpecMatSph[, 1]
    
    signalMultipleApproMatScaSph <- matrix(nrow = length(unique(rawSmooBCSpecMatSph[, 1])),
                                            ncol = 5,
                                            dimnames = list(c(),
                                                            c('Batch',
                                                              'Date',
                                                              'Peak_Area',
                                                              'Peak_Area_SD',
                                                              'n')))
    
    signalMultipleApproMatScaSph[, 1] <- unique(rawSmooBCSpecMatSph[, 1])
    
    #Loop to fill the matrix
    loopI <- 1
    while (loopI <= nrow(signalMultipleApproMatScaSph)) {
      
      signalMultipleApproMatScaSph[loopI, 2] <- 
        chronologyHM20[which(chronologyHM20[, 2] == 
                               (as.numeric(signalMultipleApproMatScaSph[loopI, 1]) * 2 )), 1]
      
      signalMultipleApproMatScaSph[loopI, 5] <- 
        length(which(rawSmooBCSpecMatSph[, 1] == signalMultipleApproMatScaSph[loopI, 1]))
      
      loopI <- loopI + 1
    }
    
    rm(loopI)
    
    #Flip the wavenumber columns
    rawSmooBCSpecMatSphFlip <-
      t(apply(rawSmooBCSpecMatSph[, 2:ncol(rawSmooBCSpecMatSph)], 1, rev))
    
    
    #Temporarily store the original peak areas into an individual variable
    oriArea <- unlist(spd_area_bands(t(rawSmooBCSpecMatSphFlip),
                                     bands = c(1492, 1540), 
                                     refine = FALSE)$area)
    
    oriAreaSphMean <- unique(rawSmooBCSpecMatSph[, 1])
    
    oriAreaSphSD <- unique(rawSmooBCSpecMatSph[, 1])
    
    #Calculate the mean and sd of each sample
    for (loopI in 1:length(unique(rawSmooBCSpecMatSph[, 1]))) {
      
      oriAreaSphMean[loopI] <- mean(oriArea[which(rawSmooBCSpecMatSph[, 1] ==
                                                    unique(rawSmooBCSpecMatSph[, 1])[loopI])])
      
      oriAreaSphSD[loopI] <- sd(oriArea[which(rawSmooBCSpecMatSph[, 1] ==
                                                unique(rawSmooBCSpecMatSph[, 1])[loopI])])
      
    }
    
    #Interpolate and smooth dataset (butterworth filter of 1/100 frequency)
    UACBwPASph <- my.bw(signalMultipleApproMatScaSph[, 2], 
                     oriAreaSphMean)
  }
  
  #Alnus
  {
    #Read smoothed and baseline corrected spectra matrix
    rawSmooBCSpecMatAln <- read.csv(here('Output_R',
                                         'fingerPrintBCNorAln.csv'),
                                    header = TRUE,
                                    row.names = 1)
    
    rawSmooBCSpecMatAln <- apply(rawSmooBCSpecMatAln, c(1,2),
                                 as.numeric)
    
    #Rename the column
    colnames(rawSmooBCSpecMatAln) <- c('Batch', seq(from = 1800, to = 800,by = -4))
    
    signalMultipleApproMatScaAln <- matrix(nrow = nrow(rawSmooBCSpecMatAln),
                                           ncol = 2,
                                           dimnames = list(c(),
                                                           c('Batch',
                                                             'Peak_Area')))
    
    signalMultipleApproMatScaAln[, 1] <- rawSmooBCSpecMatAln[, 1]
    
    signalMultipleApproMatScaAln <- matrix(nrow = length(unique(rawSmooBCSpecMatAln[, 1])),
                                           ncol = 5,
                                           dimnames = list(c(),
                                                           c('Batch',
                                                             'Date',
                                                             'Peak_Area',
                                                             'Peak_Area_SD',
                                                             'n')))
    
    signalMultipleApproMatScaAln[, 1] <- unique(rawSmooBCSpecMatAln[, 1])
    
    #Loop to fill the matrix
    loopI <- 1
    while (loopI <= nrow(signalMultipleApproMatScaAln)) {
      
      signalMultipleApproMatScaAln[loopI, 2] <- 
        chronologyHM20[which(chronologyHM20[, 2] == 
                               (as.numeric(signalMultipleApproMatScaAln[loopI, 1]) * 2 )), 1]
      
      signalMultipleApproMatScaAln[loopI, 5] <- 
        length(which(rawSmooBCSpecMatAln[, 1] == signalMultipleApproMatScaAln[loopI, 1]))
      
      loopI <- loopI + 1
    }
    
    rm(loopI)
    
    #Flip the wavenumber columns
    rawSmooBCSpecMatAlnFlip <-
      t(apply(rawSmooBCSpecMatAln[, 2:ncol(rawSmooBCSpecMatAln)], 1, rev))
    
    
    #Temporarily store the original peak areas into an individual variable
    oriArea <- unlist(spd_area_bands(t(rawSmooBCSpecMatAlnFlip),
                                     bands = c(1492, 1540), 
                                     refine = FALSE)$area)
    
    oriAreaAlnMean <- unique(rawSmooBCSpecMatAln[, 1])
    
    oriAreaAlnSD <- unique(rawSmooBCSpecMatAln[, 1])
    
    #Calculate the mean and sd of each sample
    for (loopI in 1:length(unique(rawSmooBCSpecMatAln[, 1]))) {
      
      oriAreaAlnMean[loopI] <- mean(oriArea[which(rawSmooBCSpecMatAln[, 1] ==
                                                    unique(rawSmooBCSpecMatAln[, 1])[loopI])])
      
      oriAreaAlnSD[loopI] <- sd(oriArea[which(rawSmooBCSpecMatAln[, 1] ==
                                                unique(rawSmooBCSpecMatAln[, 1])[loopI])])
      
    }
    
    #Interpolate and smooth dataset (butterworth filter of 1/100 frequency)
    UACBwPAAln <- my.bw(signalMultipleApproMatScaAln[, 2], 
                        oriAreaAlnMean)
  }
  
  #Calluna
  {
    #Read smoothed and baseline corrected spectra matrix
    rawSmooBCSpecMatCal <- read.csv(here('Output_R',
                                         'fingerPrintBCNorCal.csv'),
                                    header = TRUE,
                                    row.names = 1)
    
    rawSmooBCSpecMatCal <- apply(rawSmooBCSpecMatCal, c(1,2),
                                 as.numeric)
    
    #Rename the column
    colnames(rawSmooBCSpecMatCal) <- c('Batch', seq(from = 1800, to = 800,by = -4))
    
    signalMultipleApproMatScaCal <- matrix(nrow = nrow(rawSmooBCSpecMatCal),
                                           ncol = 2,
                                           dimnames = list(c(),
                                                           c('Batch',
                                                             'Peak_Area')))
    
    signalMultipleApproMatScaCal[, 1] <- rawSmooBCSpecMatCal[, 1]
    
    signalMultipleApproMatScaCal <- matrix(nrow = length(unique(rawSmooBCSpecMatCal[, 1])),
                                           ncol = 5,
                                           dimnames = list(c(),
                                                           c('Batch',
                                                             'Date',
                                                             'Peak_Area',
                                                             'Peak_Area_SD',
                                                             'n')))
    
    signalMultipleApproMatScaCal[, 1] <- unique(rawSmooBCSpecMatCal[, 1])
    
    #Loop to fill the matrix
    loopI <- 1
    while (loopI <= nrow(signalMultipleApproMatScaCal)) {
      
      signalMultipleApproMatScaCal[loopI, 2] <- 
        chronologyHM20[which(chronologyHM20[, 2] == 
                               (as.numeric(signalMultipleApproMatScaCal[loopI, 1]) * 2 )), 1]
      
      signalMultipleApproMatScaCal[loopI, 5] <- 
        length(which(rawSmooBCSpecMatCal[, 1] == signalMultipleApproMatScaCal[loopI, 1]))
      
      loopI <- loopI + 1
    }
    
    rm(loopI)
    
    #Flip the wavenumber columns
    rawSmooBCSpecMatCalFlip <-
      t(apply(rawSmooBCSpecMatCal[, 2:ncol(rawSmooBCSpecMatCal)], 1, rev))
    
    
    #Temporarily store the original peak areas into an individual variable
    oriArea <- unlist(spd_area_bands(t(rawSmooBCSpecMatCalFlip),
                                     bands = c(1492, 1540), 
                                     refine = FALSE)$area)
    
    oriAreaCalMean <- unique(rawSmooBCSpecMatCal[, 1])
    
    oriAreaCalSD <- unique(rawSmooBCSpecMatCal[, 1])
    
    #Calculate the mean and sd of each sample
    for (loopI in 1:length(unique(rawSmooBCSpecMatCal[, 1]))) {
      
      oriAreaCalMean[loopI] <- mean(oriArea[which(rawSmooBCSpecMatCal[, 1] ==
                                                    unique(rawSmooBCSpecMatCal[, 1])[loopI])])
      
      oriAreaCalSD[loopI] <- sd(oriArea[which(rawSmooBCSpecMatCal[, 1] ==
                                                unique(rawSmooBCSpecMatCal[, 1])[loopI])])
      
    }
    
    #Interpolate and smooth dataset (butterworth filter of 1/100 frequency)
    UACBwPACal <- my.bw(signalMultipleApproMatScaCal[, 2], 
                        oriAreaCalMean)
  }
}

#Output the raw and smoothed UAC signals
{
  #Raw UAC signals
  signalMultipleApproMatScaSph[, 3] <- oriAreaSphMean
  signalMultipleApproMatScaSph[, 4] <- oriAreaSphSD
  signalMultipleApproMatScaAln[, 3] <- oriAreaAlnMean
  signalMultipleApproMatScaAln[, 4] <- oriAreaAlnSD
  signalMultipleApproMatScaCal[, 3] <- oriAreaCalMean
  signalMultipleApproMatScaCal[, 4] <- oriAreaCalSD
  
  signalMultipleApproMatSca <- rbind(cbind(Taxa = 'Sph',
                                           signalMultipleApproMatScaSph),
                                     cbind(Taxa = 'Aln',
                                           signalMultipleApproMatScaAln),
                                     cbind(Taxa = 'Cal',
                                           signalMultipleApproMatScaCal))
  
  write.csv(signalMultipleApproMatSca, file = here('Output_R',
                                                   'signalMultipleApproMatSca.csv'))
  
  #Smoothed UAC signals
  signalMultipleApproMatScaSmoo <- rbind(cbind(Taxa = 'Sph',
                                               UACBwPASph),
                                         cbind(Taxa = 'Aln',
                                               UACBwPAAln),
                                         cbind(Taxa = 'Cal',
                                               UACBwPACal))
  
  colnames(signalMultipleApproMatScaSmoo)[2:3] <- c('Date', 'Peak_Area')
  
  write.csv(signalMultipleApproMatScaSmoo, file = here('Output_R', 'signalMultipleApproMatScaSmoo.csv'))
  
  
}