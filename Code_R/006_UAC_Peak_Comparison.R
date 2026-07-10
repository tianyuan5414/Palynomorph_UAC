#This script aims at finding UAC minima and maxima in Sphagnum UAC time series,
 #and compare correlation between three taxa sectioned by UAC maxima and minima
#Load necessary packages
{
  #Moving average calculation
  library(zoo)
}

#Declare custom functions
{
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
  
  my.corr.test2 <- function (UAC, UACAGE, OtherProxyAge, OtherData){
    A <- approx(OtherProxyAge, OtherData,UACAGE)
    Al <- length(na.omit(A$y))
    return(cor.test(as.vector(na.omit(A$y)),as.vector(na.omit(UAC[1:Al])), method = "pearson"))
  }
  
  bw.test.mc <- function(Ages, Data, UACData){
    
    cor.pearson <- c()
    
    bw.df <- my.bw(Ages,Data)
    bw.temp <- data.frame(A = bw.df$Age, B = bw.df$Filt.Data)
    C <- cor.test(UACData,
                  bw.temp$B,
                  method = "pearson")
    
    cor.pearson <- append(cor.pearson, C$estimate)
    
    return(cor.pearson)
  }
  
  #generate a set of random numbers of the same length
  rednoise.test <- function (TimeSteps, UACData){
    rednoise <- colored_noise(timesteps = TimeSteps, mean = 0, sd = 1, phi = 0.99)
    my.bw(seq(0,TimeSteps - 1 ,1), rednoise)
    bw.test.mc(0:(TimeSteps - 1),rednoise, UACData)
  }
  
  #Function for identifying the sign of number, 1 is given to positive, and -1 is to negative
  signNum <- function(number) {
    tempSign <- vector(length = length(number))
    
    for (loopI in 1:length(number)) {
      if (number[loopI] < 0) {
        tempSign[loopI] <- -1
      }else if (number[loopI] > 0) {
        tempSign[loopI] <- 1
      }else {
        tempSign[loopI] <- 0
      }
    }
   
    return(tempSign)
  }
  
  #Function for identifying the effective local extremes
  findLocalExt <- function(UACData, 
                           UACDataOri, 
                           miniN,
                           UACDataOri1) {
    diffSpec <- diff(UACData[, 3])
    
    diffSpecSign <- signNum(diffSpec)
    
    peakInd <- c()
    
    for (loopI in 1:(length(diffSpecSign) - 1)) {
      if (diffSpecSign[loopI] * diffSpecSign[loopI + 1] < 0) {
        
        peakInd[length(peakInd) + 1] <- loopI
        
      }
    }
    
    peakAge <- UACData[peakInd, 2]
    
    peakValidAge <- c()
    
    loopI <- 1
    while (loopI <= (length(peakAge) - 1)) {
      
      loopJ <- 1
      while ((loopI + loopJ) <= length(peakAge)) {
        
        if (loopI == 1) {
          
          ageLeft <- min(UACDataOri[, 3])
          
          ageRight <- (peakAge[loopI + loopJ - 1] + peakAge[loopI + loopJ]) / 2
          
        }else if (loopI == (length(peakAge) - 1)) {
          
          ageLeft <- (peakAge[loopI] + peakAge[loopI - 1]) / 2
          
          ageRight <- max(UACDataOri[, 3])
          
        } else {
          
          ageLeft <- (peakAge[loopI] + peakAge[loopI - 1]) / 2
          
          ageRight <- (peakAge[loopI + loopJ - 1] + peakAge[loopI + loopJ]) / 2
          
        }
        
        tempNumN <- length(UACDataOri[which((UACDataOri[, 3] >= ageLeft) &
                                              (UACDataOri[, 3] <= ageRight)), 3])
        
        tempNumN1 <- length(UACDataOri1[which((UACDataOri1[, 3] >= ageLeft) &
                                                (UACDataOri1[, 3] <= ageRight)), 3])
        
        if ((tempNumN >= miniN) & (tempNumN1 >= miniN)) {
          
          peakValidAge[length(peakValidAge) + 1] <- ageRight
          
          break
        } else {
          
          loopJ <- loopJ + 1
          
        }
      }
      
      loopI <- loopI + loopJ + 1
      
    }
    
    return(peakValidAge)
  }
  
  #Function for calculating Pearson correlation in sectioned parts
  correLocalExt <- function(UACDataOri,
                            UACDataOri1,
                            sectionAge) {
    
    UAC3 <- rollmean(UACDataOri[,4],3)
    UAC.AGE3 <- rollmean(UACDataOri[,3],3)
    UAC31 <- rollmean(UACDataOri1[,4],3)
    UAC.AGE31 <- rollmean(UACDataOri1[,3],3)
    
    #Matrix for sectioned correlation
    correUACProxy <- matrix(nrow = length(sectionAge) + 1, ncol = 4,
                            dimnames = list(c(),
                                            c('Sections1', 
                                              'Sections2',
                                              'r','p')))
    
    for (loopI in 1:(length(sectionAge) + 1)) {
      
      if(loopI == 1) {
        
        tempInd <- which(UAC.AGE3 <= sectionAge[loopI])
        tempInd1 <- which(UAC.AGE31 <= sectionAge[loopI])
        
        tempUAC3 <- UAC3[tempInd]
        tempUAC.AGE3 <- UAC.AGE3[tempInd]
        
        tempUAC31 <- UAC31[tempInd1]
        tempUAC.AGE31 <- UAC.AGE31[tempInd1]
        
      } else if (loopI == (length(sectionAge) + 1)) {
        
        tempInd <- which(UAC.AGE3 >= sectionAge[loopI - 1])
        tempInd1 <- which(UAC.AGE31 >= sectionAge[loopI - 1])
        
        tempUAC3 <- UAC3[tempInd]
        tempUAC.AGE3 <- UAC.AGE3[tempInd]
        
        tempUAC31 <- UAC31[tempInd1]
        tempUAC.AGE31 <- UAC.AGE31[tempInd1]
        
      } else {
        
        tempInd <- which((UAC.AGE3 >= sectionAge[loopI - 1]) & 
                           (UAC.AGE3 <= sectionAge[loopI]))
        tempInd1 <- which((UAC.AGE31 >= sectionAge[loopI - 1]) &
                            (UAC.AGE31 <= sectionAge[loopI]))
        
        tempUAC3 <- UAC3[tempInd]
        tempUAC.AGE3 <- UAC.AGE3[tempInd]
        
        tempUAC31 <- UAC31[tempInd1]
        tempUAC.AGE31 <- UAC.AGE31[tempInd1]
        
      }
      
      correUACProxy[loopI, 1:2] <- range(tempUAC.AGE3) 
      
      correUACProxy[loopI, 3] <- my.corr.test2(tempUAC3, tempUAC.AGE3, tempUAC.AGE31, tempUAC31)$estimate
      
      correUACProxy[loopI, 4] <- round(my.corr.test2(tempUAC3, tempUAC.AGE3, tempUAC.AGE31, tempUAC31)$p.value, 4)
      
    }
    
    return(correUACProxy)
  }
}

#Read UAC signals
{
  UACSignals <- read.csv(file = file.path('Output_R', 'signalMultipleApproMatScaSmoo.csv'),
                                                header = TRUE, row.names = 1)
  UACSignals <- UACSignals[which(UACSignals[,2] <= 4000),]
  
  UACSignalsOri <- read.csv(file = file.path('Output_R', 'signalMultipleApproMatSca.csv'),
                                            header = TRUE, row.names = 1)
  UACSignalsOri <- UACSignalsOri[which(UACSignalsOri[,3] <= 4000),]
  
  UACSignalsSph <- UACSignals[which(UACSignals[, 1] == 'Sph'),]
  UACSignalsAln <- UACSignals[which(UACSignals[, 1] == 'Aln'),]
  UACSignalsCal <- UACSignals[which(UACSignals[, 1] == 'Cal'),]
  
  UACSignalsOriSph <- UACSignalsOri[which(UACSignalsOri[, 1] == 'Sph'),]
  UACSignalsOriAln <- UACSignalsOri[which(UACSignalsOri[, 1] == 'Aln'),]
  UACSignalsOriCal <- UACSignalsOri[which(UACSignalsOri[, 1] == 'Cal'),]
}

#Find the local minima and maxima in Sphagnum UAC series
{
  #Sphagnum & Alnus
  correUACSectionSphAln <- correLocalExt(UACSignalsOriSph,
                                   UACSignalsOriAln,
                                   findLocalExt(UACSignalsSph,
                                                UACSignalsOriSph,
                                                11,
                                                UACSignalsOriAln))
  
  #Sphagnum & Calluna
  correUACSectionSphCal <- correLocalExt(UACSignalsOriSph,
                                         UACSignalsOriCal,
                                         findLocalExt(UACSignalsSph,
                                                      UACSignalsOriSph,
                                                      11,
                                                      UACSignalsOriCal))
  
  #Alnus & Calluna
  correUACSectionAlnCal <- correLocalExt(UACSignalsOriAln,
                                         UACSignalsOriCal,
                                         findLocalExt(UACSignalsAln,
                                                      UACSignalsOriAln,
                                                      11,
                                                      UACSignalsOriCal))
}

#Output correlation results
{
  correUACSection <- rbind(cbind(Taxa = 'Sphagnum & Alnus',
                                 correUACSectionSphAln),
                           cbind(Taxa = 'Sphagnum & Calluna',
                                 correUACSectionSphCal),
                           cbind(Taxa = 'Alnus & Calluna',
                                 correUACSectionAlnCal))
  
  write.csv(correUACSection, file = file.path('Output_R',
                                              'correUACSection.csv'))
}