#This script aims at comparing UAC signals with other proxy records (TSI, temperature, WMI, precipitation)

#Load necessary packages
{
  #Plot
  library(ggplot2)
  library(tidyverse)
  library(tidypaleo)
  library(dplyr)
  library(deeptime)
  library(patchwork)
  
  #Combine plot
  library(ggpubr)
  
  #Prevent axis from intersection
  library(ggh4x)
  
  #Load arial font into R
  library(showtext)
  font_add(family = "arial", regular = here('Fonts', 'arial.ttf'))
  showtext_auto()
  
  #Palette
  library(RColorBrewer)
  
  #Regression expression
  library(ggpmisc)
  
  #Moving average calculation
  library(zoo)
}

#Declare necessary functions
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
  
  meanProxyRecord <- function(dataFrame, proxyData) {
    tempMeanProxyRecord <- vector(length = nrow(dataFrame))
    
    for (loopI in 1:length(tempMeanProxyRecord)) {
      
      tempMeanProxyRecord[loopI] <- mean(proxyData[which((proxyData[, 1] >= dataFrame[loopI, 2]) &
                                                           (proxyData[, 1] <= dataFrame[loopI, 3])), 2])
      
    }
    
    return(tempMeanProxyRecord)
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
        
        tempNumN1 <- length(UACDataOri1[which((UACDataOri1[, 1] >= ageLeft) &
                                                (UACDataOri1[, 1] <= ageRight)), 2])
        
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
    UAC31 <- rollmean(UACDataOri1[,2],3)
    UAC.AGE31 <- rollmean(UACDataOri1[,1],3)
    
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
  
  my.corr.test2 <- function (UAC, UACAGE, OtherProxyAge, OtherData){
    A <- approx(OtherProxyAge, OtherData,UACAGE)
    Al <- length(na.omit(A$y))
    return(cor.test(as.vector(na.omit(A$y)),as.vector(na.omit(UAC[1:Al])), method = "pearson"))
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

#Read TSI records (Now using newer TSI data from Wu et al., 2018)
{
  #Read TSI data
  TSIRaw <- apply(as.matrix(read.csv(file.path('Data_R',
                                               'SATIRE-M_wu18_tsi.csv'))), c(1,2),
                  as.numeric)
  
  #Transfer tbe BC/AD ages into Cal. years BP
  TSIRaw[, 1] <- 1950 - TSIRaw[ ,1]
  
  #Reverse the order of row
  TSIRaw <- TSIRaw[order(TSIRaw[, 1]),]
  
  #Rename the columns
  colnames(TSIRaw) <- c('Years_BP', 'dTSI')
  
  #Transfer the matrix into data frame
  TSIRaw <- as.data.frame(TSIRaw)
  
  #Smooth the dTSI data at level of 80 years
  TSISmooPlot <- my.bw(TSIRaw[, 1], TSIRaw[, 2])
}

#Greenland temperature
{
  #Read raw data
  GreenTempRaw <- apply(as.matrix(read.csv(file.path('Data_R', 'Greenland_Temp.csv'))), c(1,2),
                        as.numeric)
  
  #Reverse the order of row
  GreenTempRaw <- GreenTempRaw[order(GreenTempRaw[, 1]),]
  
  #Rename the columns
  colnames(GreenTempRaw) <- c('Years_BP', 'GreenTemp',' Upper', 'Lower')
  
  #Transfer the matrix into data frame
  GreenTempRaw <- as.data.frame(GreenTempRaw)
  
  #Change data to temperatureanomaly
  GreenTempRaw[, 2] <- GreenTempRaw[, 2] - mean(GreenTempRaw[, 2])
  GreenTempRaw[, 3] <- GreenTempRaw[, 3] - mean(GreenTempRaw[, 3])
  GreenTempRaw[, 4] <- GreenTempRaw[, 4] - mean(GreenTempRaw[, 4])
  
  #Smooth the GreenTemp data at level of 80 years
  GreenTempBW <- my.bw(GreenTempRaw[, 1], GreenTempRaw[, 2])
}

#West Mediterranean index
{
  #Read raw data
  WMIRaw <- apply(as.matrix(read.csv(here('Data_R', 'WMI_AitBrahim.csv'))), c(1,2),
                  as.numeric)
  
  #Reverse the order of row
  WMIRaw <- WMIRaw[order(WMIRaw[, 1]),]
  
  #Rename the columns
  colnames(WMIRaw) <- c('Years_BP', 'WMI')
  
  #Transfer the matrix into data frame
  WMIRaw <- as.data.frame(WMIRaw)
  
  #Smooth the dWMI data at level of 80 years
  WMIBW <- my.bw(WMIRaw[, 1], WMIRaw[, 2])
}

#Scotland precipitation (stalagmite growth-rate record)
{
  #Read raw data
  SURaw <- apply(as.matrix(read.csv(here('Data_R', 'SUcomp_Baker.csv'))), c(1,2),
                 as.numeric)
  
  #Reverse the order of row
  SURaw <- SURaw[order(SURaw[, 1]),]
  
  #Rename the columns
  colnames(SURaw) <- c('Years_BP', 'SU')
  
  #Transfer the matrix into data frame
  SURaw <- as.data.frame(SURaw)
  
  #Smooth the dSU data at level of 80 years
  SUBW <- my.bw(SURaw[, 1], SURaw[, 2])
}

#Ireland bog water table
{
  #Read raw data
  IreBogRaw <- apply(as.matrix(read.csv(file.path('Data_R', 'IrishBogWT.csv'))), c(1,2),
                     as.numeric)
  
  #Reverse the order of row
  IreBogRaw <- IreBogRaw[order(IreBogRaw[, 1]),]
  
  #Rename the columns
  colnames(IreBogRaw) <- c('Years_BP', 'IreBog')
  
  #Transfer the matrix into data frame
  IreBogRaw <- as.data.frame(IreBogRaw)
  
  #Smooth the dIreBog data at level of 80 years
  IreBogBW <- my.bw(IreBogRaw[, 1], IreBogRaw[, 2])
}

#Read Walton moss data
{
  #Read WM data
  WMRaw <- read.csv(file.path('Data_R', 'Daley2010.csv'))
  
  #Rename the columns
  colnames(WMRaw) <- c('Years_BP_d18O', 'd18O', 'Years_BP_WTD', 'WTD')
  
  #Transfer the matrix into data frame
  WMRaw <- as.data.frame(WMRaw)
  
  #Smooth the dWMI data at level of 80 years
  WTDBW <- my.bw(WMRaw[, 3], WMRaw[, 4])
  d18OBW <- my.bw(WMRaw[which(!is.na(WMRaw[, 1])), 1], WMRaw[which(!is.na(WMRaw[, 1])), 2])
}

#Read VSSI data
{
  vssiRaw <- read.csv(here('Data_R', 'VolSulferMat.csv'), row.names = 1)
  
  vssiRaw[, 1] <- -vssiRaw[, 1]
  vssiRaw[, 3] <- -vssiRaw[, 3]
  
  vssiBW <- my.bw(vssiRaw[, 1], vssiRaw[, 2])
}

#Pearson correlation analysis
{
  #for 3pt mean
  UAC3Sph <- rollmean(UACSignalsOri[which(UACSignalsOri[, 1] == 'Sph'),4],3)
  UAC.AGE3Sph <- rollmean(UACSignalsOri[which(UACSignalsOri[, 1] == 'Sph'),3],3)
  
  UAC3Aln <- rollmean(UACSignalsOri[which(UACSignalsOri[, 1] == 'Aln'),4],3)
  UAC.AGE3Aln <- rollmean(UACSignalsOri[which(UACSignalsOri[, 1] == 'Aln'),3],3)
  
  UAC3Cal <- rollmean(UACSignalsOri[which(UACSignalsOri[, 1] == 'Cal'),4],3)
  UAC.AGE3Cal <- rollmean(UACSignalsOri[which(UACSignalsOri[, 1] == 'Cal'),3],3)
  
  #Pearson correlation between UAC and other proxy records (smoothed)
  correUACProxyRecord <- matrix(nrow = 15, ncol = 3,
                                dimnames = list(c(),
                                                c('Proxy_Records', 'r','p')))
  
  correUACProxyRecord[, 1] <- c('Sphagnum & TSI', 
                                'Sphagnum & GLTemp',
                                'Sphagnum & WMI',
                                'Sphagnum & WMO',
                                'Sphagnum & VSSI',
                                'Alnus & TSI', 
                                'Alnus & GLTemp',
                                'Alnus & WMI',
                                'Alnus & WMO',
                                'Alnus & VSSI',
                                'Calluna & TSI', 
                                'Calluna & GLTemp',
                                'Calluna & WMI',
                                'Calluna & WMO',
                                'Calluna & VSSI')
  
  correUACProxyRecord[1, 2] <- my.corr.test2(UAC3Sph, UAC.AGE3Sph, TSISmooPlot[,1], TSISmooPlot[,2])$estimate
  
  correUACProxyRecord[1, 3] <- round(my.corr.test2(UAC3Sph, UAC.AGE3Sph, TSISmooPlot[,1], TSISmooPlot[,2])$p.value, 4)
  
  correUACProxyRecord[2, 2] <- my.corr.test2(UAC3Sph, UAC.AGE3Sph, GreenTempBW[,1], GreenTempBW[,2])$estimate
  
  correUACProxyRecord[2, 3] <- round(my.corr.test2(UAC3Sph, UAC.AGE3Sph, GreenTempBW[,1], GreenTempBW[,2])$p.value, 4)
  
  correUACProxyRecord[3, 2] <- my.corr.test2(UAC3Sph, UAC.AGE3Sph, WMIBW[,1], WMIBW[,2])$estimate
  
  correUACProxyRecord[3, 3] <- round(my.corr.test2(UAC3Sph, UAC.AGE3Sph, WMIBW[,1], WMIBW[,2])$p.value, 4)
  
  correUACProxyRecord[4, 2] <- my.corr.test2(UAC3Sph, UAC.AGE3Sph, d18OBW[,1], d18OBW[,2])$estimate
  
  correUACProxyRecord[4, 3] <- round(my.corr.test2(UAC3Sph, UAC.AGE3Sph, d18OBW[,1], d18OBW[,2])$p.value, 4)
  
  correUACProxyRecord[5, 2] <- my.corr.test2(UAC3Sph, UAC.AGE3Sph, vssiBW[,1], vssiBW[,2])$estimate
  
  correUACProxyRecord[5, 3] <- round(my.corr.test2(UAC3Sph, UAC.AGE3Sph, vssiBW[,1], vssiBW[,2])$p.value, 4)
  
  correUACProxyRecord[6, 2] <- my.corr.test2(UAC3Aln, UAC.AGE3Aln, TSISmooPlot[,1], TSISmooPlot[,2])$estimate
  
  correUACProxyRecord[6, 3] <- round(my.corr.test2(UAC3Aln, UAC.AGE3Aln, TSISmooPlot[,1], TSISmooPlot[,2])$p.value, 4)
  
  correUACProxyRecord[7, 2] <- my.corr.test2(UAC3Aln, UAC.AGE3Aln, GreenTempBW[,1], GreenTempBW[,2])$estimate
  
  correUACProxyRecord[7, 3] <- round(my.corr.test2(UAC3Aln, UAC.AGE3Aln, GreenTempBW[,1], GreenTempBW[,2])$p.value, 4)
  
  correUACProxyRecord[8, 2] <- my.corr.test2(UAC3Aln, UAC.AGE3Aln, WMIBW[,1], WMIBW[,2])$estimate
  
  correUACProxyRecord[8, 3] <- round(my.corr.test2(UAC3Aln, UAC.AGE3Aln, WMIBW[,1], WMIBW[,2])$p.value, 4)
  
  correUACProxyRecord[9, 2] <- my.corr.test2(UAC3Aln, UAC.AGE3Aln, d18OBW[,1], d18OBW[,2])$estimate
  
  correUACProxyRecord[9, 3] <- round(my.corr.test2(UAC3Aln, UAC.AGE3Aln, d18OBW[,1], d18OBW[,2])$p.value, 4)
  
  correUACProxyRecord[10, 2] <- my.corr.test2(UAC3Aln, UAC.AGE3Aln, vssiBW[,1], vssiBW[,2])$estimate
  
  correUACProxyRecord[10, 3] <- round(my.corr.test2(UAC3Aln, UAC.AGE3Aln, vssiBW[,1], vssiBW[,2])$p.value, 4)
  
  correUACProxyRecord[11, 2] <- my.corr.test2(UAC3Cal, UAC.AGE3Cal, TSISmooPlot[,1], TSISmooPlot[,2])$estimate
  
  correUACProxyRecord[11, 3] <- round(my.corr.test2(UAC3Cal, UAC.AGE3Cal, TSISmooPlot[,1], TSISmooPlot[,2])$p.value, 4)
  
  correUACProxyRecord[12, 2] <- my.corr.test2(UAC3Cal, UAC.AGE3Cal, GreenTempBW[,1], GreenTempBW[,2])$estimate
  
  correUACProxyRecord[12, 3] <- round(my.corr.test2(UAC3Cal, UAC.AGE3Cal, GreenTempBW[,1], GreenTempBW[,2])$p.value, 4)
  
  correUACProxyRecord[13, 2] <- my.corr.test2(UAC3Cal, UAC.AGE3Cal, WMIBW[,1], WMIBW[,2])$estimate
  
  correUACProxyRecord[13, 3] <- round(my.corr.test2(UAC3Cal, UAC.AGE3Cal, WMIBW[,1], WMIBW[,2])$p.value, 4)
  
  correUACProxyRecord[14, 2] <- my.corr.test2(UAC3Cal, UAC.AGE3Cal, d18OBW[,1], d18OBW[,2])$estimate
  
  correUACProxyRecord[14, 3] <- round(my.corr.test2(UAC3Cal, UAC.AGE3Cal, d18OBW[,1], d18OBW[,2])$p.value, 4)
  
  correUACProxyRecord[15, 2] <- my.corr.test2(UAC3Cal, UAC.AGE3Cal, vssiBW[,1], vssiBW[,2])$estimate
  
  correUACProxyRecord[15, 3] <- round(my.corr.test2(UAC3Cal, UAC.AGE3Cal, vssiBW[,1], vssiBW[,2])$p.value, 4)
  
  #Output the correlation results
  write.csv(correUACProxyRecord, file = file.path('Output_R', 'correUACProxyRecord.csv'))
}

#Read sectioned correlation data
{
  correUACSection <- read.csv(file = file.path('Output_R', 'correUACSection.csv'),
                              header = TRUE, row.names = 1)
}

#Compare the strength of correlation (pearson r) between different taxa's UAC with different proxy records
{
  mulProxCom <- cbind(correUACSection, 
                      TSI = meanProxyRecord(correUACSection, TSISmooPlot),
                      GLTemp = meanProxyRecord(correUACSection, GreenTempBW),
                      WMI = meanProxyRecord(correUACSection, WMIBW),
                      ScoPre = meanProxyRecord(correUACSection, SUBW),
                      IreBog = meanProxyRecord(correUACSection, IreBogBW),
                      WMO = meanProxyRecord(correUACSection, d18OBW),
                      VSSI = meanProxyRecord(correUACSection, vssiBW))
  
  #Scatter plot with linear regression analysis
  #Sphagnum & Alnus
  {
    mulProxComSphAln <- mulProxCom[which(mulProxCom[, 1] == 'Sphagnum & Alnus'), ]
    
    scatterPlotSphAlnTSI <- ggplot(mulProxComSphAln,
                                   aes(x = r,
                                       y = TSI)) +
      geom_point() +
      geom_smooth(
                  method = "lm",
                  color = 'black') +
      stat_poly_eq(aes(
                       label = paste(after_stat(eq.label), 
                                     after_stat(rr.label), 
                                     sep = "*\", \"*")), 
                   formula = y ~ x, parse = TRUE) +
      labs(x = 'Pearson r (Sphagnum & Alnus)',
           y = 'TSI (Wm-2)') + 
      theme_classic() + 
      theme(plot.title = element_blank(), 
            axis.title.y = element_text(size = 14, family = 'arial'),
            axis.text.y = element_text(size = 12, family = 'arial'),
            axis.title.x = element_text(size = 14, family = 'arial'),
            axis.text.x = element_text(size = 12, family = 'arial'),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
            legend.position = 'none'
      )
    
    scatterPlotSphAlnGLTemp <- ggplot(mulProxComSphAln,
                                   aes(x = r,
                                       y = GLTemp)) +
      geom_point() +
      geom_smooth(
        method = "lm",
        color = 'black') +
      stat_poly_eq(aes(
        label = paste(after_stat(eq.label), 
                      after_stat(rr.label), 
                      sep = "*\", \"*")), 
        formula = y ~ x, parse = TRUE) +
      labs(x = 'Pearson r (Sphagnum & Alnus)',
           y = 'Greenland temperature anomaly') + 
      theme_classic() + 
      theme(plot.title = element_blank(), 
            axis.title.y = element_text(size = 14, family = 'arial'),
            axis.text.y = element_text(size = 12, family = 'arial'),
            axis.title.x = element_text(size = 14, family = 'arial'),
            axis.text.x = element_text(size = 12, family = 'arial'),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
            legend.position = 'none'
      )
    
    scatterPlotSphAlnWMI <- ggplot(mulProxComSphAln,
                                   aes(x = r,
                                       y = WMI)) +
      geom_point() +
      geom_smooth(
        method = "lm",
        color = 'black') +
      stat_poly_eq(aes(
        label = paste(after_stat(eq.label), 
                      after_stat(rr.label), 
                      sep = "*\", \"*")), 
        formula = y ~ x, parse = TRUE) +
      labs(x = 'Pearson r (Sphagnum & Alnus)',
           y = 'Western Mediterranean Index (WMI)') + 
      theme_classic() + 
      theme(plot.title = element_blank(), 
            axis.title.y = element_text(size = 14, family = 'arial'),
            axis.text.y = element_text(size = 12, family = 'arial'),
            axis.title.x = element_text(size = 14, family = 'arial'),
            axis.text.x = element_text(size = 12, family = 'arial'),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
            legend.position = 'none'
      )
    
    scatterPlotSphAlnScoPre <- ggplot(mulProxComSphAln,
                                   aes(x = r,
                                       y = ScoPre)) +
      geom_point() +
      geom_smooth(
        method = "lm",
        color = 'black') +
      stat_poly_eq(aes(
        label = paste(after_stat(eq.label), 
                      after_stat(rr.label), 
                      sep = "*\", \"*")), 
        formula = y ~ x, parse = TRUE) +
      labs(x = 'Pearson r (Sphagnum & Alnus)',
           y = 'Scotland stalagmite growth-rate') + 
      theme_classic() + 
      theme(plot.title = element_blank(), 
            axis.title.y = element_text(size = 14, family = 'arial'),
            axis.text.y = element_text(size = 12, family = 'arial'),
            axis.title.x = element_text(size = 14, family = 'arial'),
            axis.text.x = element_text(size = 12, family = 'arial'),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
            legend.position = 'none'
      )
    
    scatterPlotSphAlnIreBog <- ggplot(mulProxComSphAln,
                                   aes(x = r,
                                       y = IreBog)) +
      geom_point() +
      geom_smooth(
        method = "lm",
        color = 'black') +
      stat_poly_eq(aes(
        label = paste(after_stat(eq.label), 
                      after_stat(rr.label), 
                      sep = "*\", \"*")), 
        formula = y ~ x, parse = TRUE) +
      labs(x = 'Pearson r (Sphagnum & Alnus)',
           y = 'Ireland bog water table') + 
      theme_classic() + 
      theme(plot.title = element_blank(), 
            axis.title.y = element_text(size = 14, family = 'arial'),
            axis.text.y = element_text(size = 12, family = 'arial'),
            axis.title.x = element_text(size = 14, family = 'arial'),
            axis.text.x = element_text(size = 12, family = 'arial'),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
            legend.position = 'none'
      )
    
    scatterPlotSphAlnWMO <- ggplot(mulProxComSphAln,
                                   aes(x = r,
                                       y = WMO)) +
      geom_point() +
      geom_smooth(
        method = "lm",
        color = 'black') +
      stat_poly_eq(aes(
        label = paste(after_stat(eq.label), 
                      after_stat(rr.label), 
                      sep = "*\", \"*")), 
        formula = y ~ x, parse = TRUE) +
      labs(x = 'Pearson r (Sphagnum & Alnus)',
           y = 'd18O') + 
      theme_classic() + 
      theme(plot.title = element_blank(), 
            axis.title.y = element_text(size = 14, family = 'arial'),
            axis.text.y = element_text(size = 12, family = 'arial'),
            axis.title.x = element_text(size = 14, family = 'arial'),
            axis.text.x = element_text(size = 12, family = 'arial'),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
            legend.position = 'none'
      )
    
    scatterPlotSphAlnVSSI <- ggplot(mulProxComSphAln,
                                   aes(x = r,
                                       y = VSSI)) +
      geom_point() +
      geom_smooth(
        method = "lm",
        color = 'black') +
      stat_poly_eq(aes(
        label = paste(after_stat(eq.label), 
                      after_stat(rr.label), 
                      sep = "*\", \"*")), 
        formula = y ~ x, parse = TRUE) +
      labs(x = 'Pearson r (Sphagnum & Alnus)',
           y = 'VSSI') + 
      theme_classic() + 
      theme(plot.title = element_blank(), 
            axis.title.y = element_text(size = 14, family = 'arial'),
            axis.text.y = element_text(size = 12, family = 'arial'),
            axis.title.x = element_text(size = 14, family = 'arial'),
            axis.text.x = element_text(size = 12, family = 'arial'),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
            legend.position = 'none'
      )
  }
  
  #Sphagnum & Calluna
  {
    mulProxComSphCal <- mulProxCom[which(mulProxCom[, 1] == 'Sphagnum & Calluna'), ]
    
    scatterPlotSphCalTSI <- ggplot(mulProxComSphCal,
                                   aes(x = r,
                                       y = TSI)) +
      geom_point() +
      geom_smooth(
        method = "lm",
        color = 'black') +
      stat_poly_eq(aes(
        label = paste(after_stat(eq.label), 
                      after_stat(rr.label), 
                      sep = "*\", \"*")), 
        formula = y ~ x, parse = TRUE) +
      labs(x = 'Pearson r (Sphagnum & Calluna)',
           y = 'TSI (Wm-2)') + 
      theme_classic() + 
      theme(plot.title = element_blank(), 
            axis.title.y = element_text(size = 14, family = 'arial'),
            axis.text.y = element_text(size = 12, family = 'arial'),
            axis.title.x = element_text(size = 14, family = 'arial'),
            axis.text.x = element_text(size = 12, family = 'arial'),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
            legend.position = 'none'
      )
    
    scatterPlotSphCalGLTemp <- ggplot(mulProxComSphCal,
                                      aes(x = r,
                                          y = GLTemp)) +
      geom_point() +
      geom_smooth(
        method = "lm",
        color = 'black') +
      stat_poly_eq(aes(
        label = paste(after_stat(eq.label), 
                      after_stat(rr.label), 
                      sep = "*\", \"*")), 
        formula = y ~ x, parse = TRUE) +
      labs(x = 'Pearson r (Sphagnum & Calluna)',
           y = 'Greenland temperature anomaly') + 
      theme_classic() + 
      theme(plot.title = element_blank(), 
            axis.title.y = element_text(size = 14, family = 'arial'),
            axis.text.y = element_text(size = 12, family = 'arial'),
            axis.title.x = element_text(size = 14, family = 'arial'),
            axis.text.x = element_text(size = 12, family = 'arial'),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
            legend.position = 'none'
      )
    
    scatterPlotSphCalWMI <- ggplot(mulProxComSphCal,
                                   aes(x = r,
                                       y = WMI)) +
      geom_point() +
      geom_smooth(
        method = "lm",
        color = 'black') +
      stat_poly_eq(aes(
        label = paste(after_stat(eq.label), 
                      after_stat(rr.label), 
                      sep = "*\", \"*")), 
        formula = y ~ x, parse = TRUE) +
      labs(x = 'Pearson r (Sphagnum & Calluna)',
           y = 'Western Mediterranean Index (WMI)') + 
      theme_classic() + 
      theme(plot.title = element_blank(), 
            axis.title.y = element_text(size = 14, family = 'arial'),
            axis.text.y = element_text(size = 12, family = 'arial'),
            axis.title.x = element_text(size = 14, family = 'arial'),
            axis.text.x = element_text(size = 12, family = 'arial'),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
            legend.position = 'none'
      )
    
    scatterPlotSphCalScoPre <- ggplot(mulProxComSphCal,
                                      aes(x = r,
                                          y = ScoPre)) +
      geom_point() +
      geom_smooth(
        method = "lm",
        color = 'black') +
      stat_poly_eq(aes(
        label = paste(after_stat(eq.label), 
                      after_stat(rr.label), 
                      sep = "*\", \"*")), 
        formula = y ~ x, parse = TRUE) +
      labs(x = 'Pearson r (Sphagnum & Calluna)',
           y = 'Scotland stalagmite growth-rate') + 
      theme_classic() + 
      theme(plot.title = element_blank(), 
            axis.title.y = element_text(size = 14, family = 'arial'),
            axis.text.y = element_text(size = 12, family = 'arial'),
            axis.title.x = element_text(size = 14, family = 'arial'),
            axis.text.x = element_text(size = 12, family = 'arial'),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
            legend.position = 'none'
      )
    
    scatterPlotSphCalIreBog <- ggplot(mulProxComSphCal,
                                      aes(x = r,
                                          y = IreBog)) +
      geom_point() +
      geom_smooth(
        method = "lm",
        color = 'black') +
      stat_poly_eq(aes(
        label = paste(after_stat(eq.label), 
                      after_stat(rr.label), 
                      sep = "*\", \"*")), 
        formula = y ~ x, parse = TRUE) +
      labs(x = 'Pearson r (Sphagnum & Calluna)',
           y = 'Ireland bog water table') + 
      theme_classic() + 
      theme(plot.title = element_blank(), 
            axis.title.y = element_text(size = 14, family = 'arial'),
            axis.text.y = element_text(size = 12, family = 'arial'),
            axis.title.x = element_text(size = 14, family = 'arial'),
            axis.text.x = element_text(size = 12, family = 'arial'),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
            legend.position = 'none'
      )
    
    scatterPlotSphCalWMO <- ggplot(mulProxComSphCal,
                                   aes(x = r,
                                       y = WMO)) +
      geom_point() +
      geom_smooth(
        method = "lm",
        color = 'black') +
      stat_poly_eq(aes(
        label = paste(after_stat(eq.label), 
                      after_stat(rr.label), 
                      sep = "*\", \"*")), 
        formula = y ~ x, parse = TRUE) +
      labs(x = 'Pearson r (Sphagnum & Calluna)',
           y = 'd18O') + 
      theme_classic() + 
      theme(plot.title = element_blank(), 
            axis.title.y = element_text(size = 14, family = 'arial'),
            axis.text.y = element_text(size = 12, family = 'arial'),
            axis.title.x = element_text(size = 14, family = 'arial'),
            axis.text.x = element_text(size = 12, family = 'arial'),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
            legend.position = 'none'
      )
    
    scatterPlotSphCalVSSI <- ggplot(mulProxComSphCal,
                                    aes(x = r,
                                        y = VSSI)) +
      geom_point() +
      geom_smooth(
        method = "lm",
        color = 'black') +
      stat_poly_eq(aes(
        label = paste(after_stat(eq.label), 
                      after_stat(rr.label), 
                      sep = "*\", \"*")), 
        formula = y ~ x, parse = TRUE) +
      labs(x = 'Pearson r (Sphagnum & Calluna)',
           y = 'VSSI') + 
      theme_classic() + 
      theme(plot.title = element_blank(), 
            axis.title.y = element_text(size = 14, family = 'arial'),
            axis.text.y = element_text(size = 12, family = 'arial'),
            axis.title.x = element_text(size = 14, family = 'arial'),
            axis.text.x = element_text(size = 12, family = 'arial'),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
            legend.position = 'none'
      )
  }
  
  #Alnus & Calluna
  {
    mulProxComAlnCal <- mulProxCom[which(mulProxCom[, 1] == 'Alnus & Calluna'), ]
    
    scatterPlotAlnCalTSI <- ggplot(mulProxComAlnCal,
                                   aes(x = r,
                                       y = TSI)) +
      geom_point() +
      geom_smooth(
        method = "lm",
        color = 'black') +
      stat_poly_eq(aes(
        label = paste(after_stat(eq.label), 
                      after_stat(rr.label), 
                      sep = "*\", \"*")), 
        formula = y ~ x, parse = TRUE) +
      labs(x = 'Pearson r (Alnus & Calluna)',
           y = 'TSI (Wm-2)') + 
      theme_classic() + 
      theme(plot.title = element_blank(), 
            axis.title.y = element_text(size = 14, family = 'arial'),
            axis.text.y = element_text(size = 12, family = 'arial'),
            axis.title.x = element_text(size = 14, family = 'arial'),
            axis.text.x = element_text(size = 12, family = 'arial'),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
            legend.position = 'none'
      )
    
    scatterPlotAlnCalGLTemp <- ggplot(mulProxComAlnCal,
                                      aes(x = r,
                                          y = GLTemp)) +
      geom_point() +
      geom_smooth(
        method = "lm",
        color = 'black') +
      stat_poly_eq(aes(
        label = paste(after_stat(eq.label), 
                      after_stat(rr.label), 
                      sep = "*\", \"*")), 
        formula = y ~ x, parse = TRUE) +
      labs(x = 'Pearson r (Alnus & Calluna)',
           y = 'Greenland temperature anomaly') + 
      theme_classic() + 
      theme(plot.title = element_blank(), 
            axis.title.y = element_text(size = 14, family = 'arial'),
            axis.text.y = element_text(size = 12, family = 'arial'),
            axis.title.x = element_text(size = 14, family = 'arial'),
            axis.text.x = element_text(size = 12, family = 'arial'),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
            legend.position = 'none'
      )
    
    scatterPlotAlnCalWMI <- ggplot(mulProxComAlnCal,
                                   aes(x = r,
                                       y = WMI)) +
      geom_point() +
      geom_smooth(
        method = "lm",
        color = 'black') +
      stat_poly_eq(aes(
        label = paste(after_stat(eq.label), 
                      after_stat(rr.label), 
                      sep = "*\", \"*")), 
        formula = y ~ x, parse = TRUE) +
      labs(x = 'Pearson r (Alnus & Calluna)',
           y = 'Western Mediterranean Index (WMI)') + 
      theme_classic() + 
      theme(plot.title = element_blank(), 
            axis.title.y = element_text(size = 14, family = 'arial'),
            axis.text.y = element_text(size = 12, family = 'arial'),
            axis.title.x = element_text(size = 14, family = 'arial'),
            axis.text.x = element_text(size = 12, family = 'arial'),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
            legend.position = 'none'
      )
    
    scatterPlotAlnCalScoPre <- ggplot(mulProxComAlnCal,
                                      aes(x = r,
                                          y = ScoPre)) +
      geom_point() +
      geom_smooth(
        method = "lm",
        color = 'black') +
      stat_poly_eq(aes(
        label = paste(after_stat(eq.label), 
                      after_stat(rr.label), 
                      sep = "*\", \"*")), 
        formula = y ~ x, parse = TRUE) +
      labs(x = 'Pearson r (Alnus & Calluna)',
           y = 'Scotland stalagmite growth-rate') + 
      theme_classic() + 
      theme(plot.title = element_blank(), 
            axis.title.y = element_text(size = 14, family = 'arial'),
            axis.text.y = element_text(size = 12, family = 'arial'),
            axis.title.x = element_text(size = 14, family = 'arial'),
            axis.text.x = element_text(size = 12, family = 'arial'),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
            legend.position = 'none'
      )
    
    scatterPlotAlnCalIreBog <- ggplot(mulProxComAlnCal,
                                      aes(x = r,
                                          y = IreBog)) +
      geom_point() +
      geom_smooth(
        method = "lm",
        color = 'black') +
      stat_poly_eq(aes(
        label = paste(after_stat(eq.label), 
                      after_stat(rr.label), 
                      sep = "*\", \"*")), 
        formula = y ~ x, parse = TRUE) +
      labs(x = 'Pearson r (Alnus & Calluna)',
           y = 'Ireland bog water table') + 
      theme_classic() + 
      theme(plot.title = element_blank(), 
            axis.title.y = element_text(size = 14, family = 'arial'),
            axis.text.y = element_text(size = 12, family = 'arial'),
            axis.title.x = element_text(size = 14, family = 'arial'),
            axis.text.x = element_text(size = 12, family = 'arial'),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
            legend.position = 'none'
      )
    
    scatterPlotAlnCalWMO <- ggplot(mulProxComAlnCal,
                                   aes(x = r,
                                       y = WMO)) +
      geom_point() +
      geom_smooth(
        method = "lm",
        color = 'black') +
      stat_poly_eq(aes(
        label = paste(after_stat(eq.label), 
                      after_stat(rr.label), 
                      sep = "*\", \"*")), 
        formula = y ~ x, parse = TRUE) +
      labs(x = 'Pearson r (Alnus & Calluna)',
           y = 'd18O') + 
      theme_classic() + 
      theme(plot.title = element_blank(), 
            axis.title.y = element_text(size = 14, family = 'arial'),
            axis.text.y = element_text(size = 12, family = 'arial'),
            axis.title.x = element_text(size = 14, family = 'arial'),
            axis.text.x = element_text(size = 12, family = 'arial'),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
            legend.position = 'none'
      )
    
    scatterPlotAlnCalVSSI <- ggplot(mulProxComAlnCal,
                                    aes(x = r,
                                        y = VSSI)) +
      geom_point() +
      geom_smooth(
        method = "lm",
        color = 'black') +
      stat_poly_eq(aes(
        label = paste(after_stat(eq.label), 
                      after_stat(rr.label), 
                      sep = "*\", \"*")), 
        formula = y ~ x, parse = TRUE) +
      labs(x = 'Pearson r (Alnus & Calluna)',
           y = 'VSSI') + 
      theme_classic() + 
      theme(plot.title = element_blank(), 
            axis.title.y = element_text(size = 14, family = 'arial'),
            axis.text.y = element_text(size = 12, family = 'arial'),
            axis.title.x = element_text(size = 14, family = 'arial'),
            axis.text.x = element_text(size = 12, family = 'arial'),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
            legend.position = 'none'
      )
  }
  
  #Arrange on one page (AlnCalTSI has the strongest correaltion)
  {
    (scatterPlotSphAlnTSI | scatterPlotSphAlnGLTemp) /
      (scatterPlotSphAlnWMI | scatterPlotSphAlnScoPre) /
      (scatterPlotSphAlnIreBog | scatterPlotSphAlnWMO)
    
    (scatterPlotSphCalTSI | scatterPlotSphCalGLTemp) /
      (scatterPlotSphCalWMI | scatterPlotSphCalScoPre) /
      (scatterPlotSphCalIreBog | scatterPlotSphCalWMO)
    
    (scatterPlotAlnCalTSI | scatterPlotAlnCalGLTemp) /
      (scatterPlotAlnCalWMI | scatterPlotAlnCalScoPre) /
      (scatterPlotAlnCalIreBog | scatterPlotAlnCalWMO)
    
    scatterPlotSphAlnTSI /
      scatterPlotSphCalTSI /
      scatterPlotAlnCalTSI
  }
}

#Sectioned correlation analysis between UAC signals and other proxy records
{
  #Sphagnum
  {
    corrSphTSI <- correLocalExt(UACSignalsOriSph,
                                TSISmooPlot,
                                findLocalExt(UACSignalsSph,
                                             UACSignalsOriSph,
                                             7,
                                             TSISmooPlot))
    
    corrSphGreenTemp <- correLocalExt(UACSignalsOriSph,
                                GreenTempBW,
                                findLocalExt(UACSignalsSph,
                                             UACSignalsOriSph,
                                             7,
                                             GreenTempBW))
    
    corrSphGreend18O <- correLocalExt(UACSignalsOriSph,
                                      d18OBW,
                                      findLocalExt(UACSignalsSph,
                                                   UACSignalsOriSph,
                                                   7,
                                                   d18OBW))
    
    corrSphGreenVSSI <- correLocalExt(UACSignalsOriSph,
                                      vssiBW,
                                      findLocalExt(UACSignalsSph,
                                                   UACSignalsOriSph,
                                                   7,
                                                   vssiBW))
    
    corrSphGreenWMI <- correLocalExt(UACSignalsOriSph[which(UACSignalsOriSph[, 3] <= max(WMIBW[, 1])),],
                                      WMIBW,
                                      findLocalExt(UACSignalsSph[which(UACSignalsSph[, 3] <= max(WMIBW[, 1])),],
                                                   UACSignalsOriSph[which(UACSignalsOriSph[, 3] <= max(WMIBW[, 1])),],
                                                   7,
                                                   WMIBW))
  }
  
  #Alnus
  {
    corrAlnTSI <- correLocalExt(UACSignalsOriAln,
                                TSISmooPlot,
                                findLocalExt(UACSignalsAln,
                                             UACSignalsOriAln,
                                             7,
                                             TSISmooPlot))
    
    corrAlnGreenTemp <- correLocalExt(UACSignalsOriAln,
                                      GreenTempBW,
                                      findLocalExt(UACSignalsAln,
                                                   UACSignalsOriAln,
                                                   7,
                                                   GreenTempBW))
    
    corrAlnGreend18O <- correLocalExt(UACSignalsOriAln,
                                      d18OBW,
                                      findLocalExt(UACSignalsAln,
                                                   UACSignalsOriAln,
                                                   7,
                                                   d18OBW))
    
    corrAlnGreenVSSI <- correLocalExt(UACSignalsOriAln,
                                      vssiBW,
                                      findLocalExt(UACSignalsAln,
                                                   UACSignalsOriAln,
                                                   7,
                                                   vssiBW))
    
    corrAlnGreenWMI <- correLocalExt(UACSignalsOriAln[which(UACSignalsOriAln[, 3] <= max(WMIBW[, 1])),],
                                     WMIBW,
                                     findLocalExt(UACSignalsAln[which(UACSignalsAln[, 3] <= max(WMIBW[, 1])),],
                                                  UACSignalsOriAln[which(UACSignalsOriAln[, 3] <= max(WMIBW[, 1])),],
                                                  5,
                                                  WMIBW))
  }
  
  #Calluna
  {
    corrCalTSI <- correLocalExt(UACSignalsOriCal,
                                TSISmooPlot,
                                findLocalExt(UACSignalsCal,
                                             UACSignalsOriCal,
                                             10,
                                             TSISmooPlot))
    
    corrCalGreenTemp <- correLocalExt(UACSignalsOriCal,
                                      GreenTempBW,
                                      findLocalExt(UACSignalsCal,
                                                   UACSignalsOriCal,
                                                   10,
                                                   GreenTempBW))
    
    corrCalGreend18O <- correLocalExt(UACSignalsOriCal,
                                      d18OBW,
                                      findLocalExt(UACSignalsCal,
                                                   UACSignalsOriCal,
                                                   10,
                                                   d18OBW))
    
    corrCalGreenVSSI <- correLocalExt(UACSignalsOriCal,
                                      vssiBW,
                                      findLocalExt(UACSignalsCal,
                                                   UACSignalsOriCal,
                                                   10,
                                                   vssiBW))
    
    corrCalGreenWMI <- correLocalExt(UACSignalsOriCal[which(UACSignalsOriCal[, 3] <= max(WMIBW[, 1])),],
                                     WMIBW,
                                     findLocalExt(UACSignalsCal[which(UACSignalsCal[, 3] <= max(WMIBW[, 1])),],
                                                  UACSignalsOriCal[which(UACSignalsOriCal[, 3] <= max(WMIBW[, 1])),],
                                                  7,
                                                  WMIBW))
  }
  
  #Output correlation results
  {
    correProxySection <- rbind(cbind(Taxa = 'Sphagnum & TSI',
                                   corrSphTSI),
                             cbind(Taxa = 'Sphagnum & GreenTemp',
                                   corrSphGreenTemp),
                             cbind(Taxa = 'Sphagnum & d18O',
                                   corrSphGreend18O),
                             cbind(Taxa = 'Sphagnum & VSSI',
                                   corrSphGreenVSSI),
                             cbind(Taxa = 'Sphagnum & WMI',
                                   corrSphGreenWMI),
                             cbind(Taxa = 'Alnus & TSI',
                                   corrAlnTSI),
                             cbind(Taxa = 'Alnus & GreenTemp',
                                   corrAlnGreenTemp),
                             cbind(Taxa = 'Alnus & d18O',
                                   corrAlnGreend18O),
                             cbind(Taxa = 'Alnus & VSSI',
                                   corrAlnGreenVSSI),
                             cbind(Taxa = 'Alnus & WMI',
                                   corrAlnGreenWMI),
                             cbind(Taxa = 'Calluna & TSI',
                                   corrCalTSI),
                             cbind(Taxa = 'Calluna & GreenTemp',
                                   corrCalGreenTemp),
                             cbind(Taxa = 'Calluna & d18O',
                                   corrCalGreend18O),
                             cbind(Taxa = 'Calluna & VSSI',
                                   corrCalGreenVSSI),
                             cbind(Taxa = 'Calluna & WMI',
                                   corrCalGreenWMI))
    
    correProxySection <- as.data.frame(correProxySection)
    correProxySection[, -1] <- apply(correProxySection[,-1], c(1,2), as.numeric)
    
    write.csv(correProxySection, file = file.path('Output_R',
                                                'correProxySection.csv'))
  }
  
  #Compare the strength of correlation (pearson r) between UAC and different proxy records
  {
    mulProxRecordCom <- cbind(correProxySection, 
                              TSI = meanProxyRecord(correProxySection, TSISmooPlot),
                              GLTemp = meanProxyRecord(correProxySection, GreenTempBW),
                              WMI = meanProxyRecord(correProxySection, WMIBW),
                              WMO = meanProxyRecord(correProxySection, d18OBW),
                              VSSI = meanProxyRecord(correProxySection, vssiBW))
    
    #Scatter plot with linear regression analysis
    #Sphagnum & TSI
    {
      mulProxRecordComSph <- mulProxRecordCom[which(mulProxRecordCom[, 1] == 'Sphagnum & TSI'), ]
      
      scatterPlotSphTSITSI <- ggplot(mulProxRecordComSph,
                                     aes(x = r,
                                         y = TSI)) +
        geom_point() +
        geom_smooth(
          method = "lm",
          color = 'black') +
        stat_poly_eq(aes(
          label = paste(after_stat(eq.label), 
                        after_stat(rr.label), 
                        sep = "*\", \"*")), 
          formula = y ~ x, parse = TRUE) +
        labs(x = 'Pearson r (Sphagnum & Alnus)',
             y = 'TSI (Wm-2)') + 
        theme_classic() + 
        theme(plot.title = element_blank(), 
              axis.title.y = element_text(size = 14, family = 'arial'),
              axis.text.y = element_text(size = 12, family = 'arial'),
              axis.title.x = element_text(size = 14, family = 'arial'),
              axis.text.x = element_text(size = 12, family = 'arial'),
              panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
              legend.position = 'none'
        )
      
      scatterPlotSphTSIGLTemp <- ggplot(mulProxRecordComSph,
                                        aes(x = r,
                                            y = GLTemp)) +
        geom_point() +
        geom_smooth(
          method = "lm",
          color = 'black') +
        stat_poly_eq(aes(
          label = paste(after_stat(eq.label), 
                        after_stat(rr.label), 
                        sep = "*\", \"*")), 
          formula = y ~ x, parse = TRUE) +
        labs(x = 'Pearson r (Sphagnum & Alnus)',
             y = 'Greenland temperature anomaly') + 
        theme_classic() + 
        theme(plot.title = element_blank(), 
              axis.title.y = element_text(size = 14, family = 'arial'),
              axis.text.y = element_text(size = 12, family = 'arial'),
              axis.title.x = element_text(size = 14, family = 'arial'),
              axis.text.x = element_text(size = 12, family = 'arial'),
              panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
              legend.position = 'none'
        )
      
      scatterPlotSphTSIWMI <- ggplot(mulProxRecordComSph,
                                     aes(x = r,
                                         y = WMI)) +
        geom_point() +
        geom_smooth(
          method = "lm",
          color = 'black') +
        stat_poly_eq(aes(
          label = paste(after_stat(eq.label), 
                        after_stat(rr.label), 
                        sep = "*\", \"*")), 
          formula = y ~ x, parse = TRUE) +
        labs(x = 'Pearson r (Sphagnum & Alnus)',
             y = 'Western Mediterranean Index (WMI)') + 
        theme_classic() + 
        theme(plot.title = element_blank(), 
              axis.title.y = element_text(size = 14, family = 'arial'),
              axis.text.y = element_text(size = 12, family = 'arial'),
              axis.title.x = element_text(size = 14, family = 'arial'),
              axis.text.x = element_text(size = 12, family = 'arial'),
              panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
              legend.position = 'none'
        )
      
      scatterPlotSphTSIWMO <- ggplot(mulProxRecordComSph,
                                     aes(x = r,
                                         y = WMO)) +
        geom_point() +
        geom_smooth(
          method = "lm",
          color = 'black') +
        stat_poly_eq(aes(
          label = paste(after_stat(eq.label), 
                        after_stat(rr.label), 
                        sep = "*\", \"*")), 
          formula = y ~ x, parse = TRUE) +
        labs(x = 'Pearson r (Sphagnum & Alnus)',
             y = 'd18O') + 
        theme_classic() + 
        theme(plot.title = element_blank(), 
              axis.title.y = element_text(size = 14, family = 'arial'),
              axis.text.y = element_text(size = 12, family = 'arial'),
              axis.title.x = element_text(size = 14, family = 'arial'),
              axis.text.x = element_text(size = 12, family = 'arial'),
              panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
              legend.position = 'none'
        )
      
      scatterPlotSphTSIVSSI <- ggplot(mulProxRecordComSph,
                                      aes(x = r,
                                          y = VSSI)) +
        geom_point() +
        geom_smooth(
          method = "lm",
          color = 'black') +
        stat_poly_eq(aes(
          label = paste(after_stat(eq.label), 
                        after_stat(rr.label), 
                        sep = "*\", \"*")), 
          formula = y ~ x, parse = TRUE) +
        labs(x = 'Pearson r (Sphagnum & Alnus)',
             y = 'VSSI') + 
        theme_classic() + 
        theme(plot.title = element_blank(), 
              axis.title.y = element_text(size = 14, family = 'arial'),
              axis.text.y = element_text(size = 12, family = 'arial'),
              axis.title.x = element_text(size = 14, family = 'arial'),
              axis.text.x = element_text(size = 12, family = 'arial'),
              panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
              legend.position = 'none'
        )
    }
    
    #Alnus & TSI
    {
      mulProxRecordComAln <- mulProxRecordCom[which(mulProxRecordCom[, 1] == 'Alnus & TSI'), ]
      
      scatterPlotAlnTSITSI <- ggplot(mulProxRecordComAln,
                                     aes(x = r,
                                         y = TSI)) +
        geom_point() +
        geom_smooth(
          method = "lm",
          color = 'black') +
        stat_poly_eq(aes(
          label = paste(after_stat(eq.label), 
                        after_stat(rr.label), 
                        sep = "*\", \"*")), 
          formula = y ~ x, parse = TRUE) +
        labs(x = 'Pearson r (Alnus & Alnus)',
             y = 'TSI (Wm-2)') + 
        theme_classic() + 
        theme(plot.title = element_blank(), 
              axis.title.y = element_text(size = 14, family = 'arial'),
              axis.text.y = element_text(size = 12, family = 'arial'),
              axis.title.x = element_text(size = 14, family = 'arial'),
              axis.text.x = element_text(size = 12, family = 'arial'),
              panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
              legend.position = 'none'
        )
      
      scatterPlotAlnTSIGLTemp <- ggplot(mulProxRecordComAln,
                                        aes(x = r,
                                            y = GLTemp)) +
        geom_point() +
        geom_smooth(
          method = "lm",
          color = 'black') +
        stat_poly_eq(aes(
          label = paste(after_stat(eq.label), 
                        after_stat(rr.label), 
                        sep = "*\", \"*")), 
          formula = y ~ x, parse = TRUE) +
        labs(x = 'Pearson r (Alnus & Alnus)',
             y = 'Greenland temperature anomaly') + 
        theme_classic() + 
        theme(plot.title = element_blank(), 
              axis.title.y = element_text(size = 14, family = 'arial'),
              axis.text.y = element_text(size = 12, family = 'arial'),
              axis.title.x = element_text(size = 14, family = 'arial'),
              axis.text.x = element_text(size = 12, family = 'arial'),
              panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
              legend.position = 'none'
        )
      
      scatterPlotAlnTSIWMI <- ggplot(mulProxRecordComAln,
                                     aes(x = r,
                                         y = WMI)) +
        geom_point() +
        geom_smooth(
          method = "lm",
          color = 'black') +
        stat_poly_eq(aes(
          label = paste(after_stat(eq.label), 
                        after_stat(rr.label), 
                        sep = "*\", \"*")), 
          formula = y ~ x, parse = TRUE) +
        labs(x = 'Pearson r (Alnus & Alnus)',
             y = 'Western Mediterranean Index (WMI)') + 
        theme_classic() + 
        theme(plot.title = element_blank(), 
              axis.title.y = element_text(size = 14, family = 'arial'),
              axis.text.y = element_text(size = 12, family = 'arial'),
              axis.title.x = element_text(size = 14, family = 'arial'),
              axis.text.x = element_text(size = 12, family = 'arial'),
              panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
              legend.position = 'none'
        )
      
      scatterPlotAlnTSIWMO <- ggplot(mulProxRecordComAln,
                                     aes(x = r,
                                         y = WMO)) +
        geom_point() +
        geom_smooth(
          method = "lm",
          color = 'black') +
        stat_poly_eq(aes(
          label = paste(after_stat(eq.label), 
                        after_stat(rr.label), 
                        sep = "*\", \"*")), 
          formula = y ~ x, parse = TRUE) +
        labs(x = 'Pearson r (Alnus & Alnus)',
             y = 'd18O') + 
        theme_classic() + 
        theme(plot.title = element_blank(), 
              axis.title.y = element_text(size = 14, family = 'arial'),
              axis.text.y = element_text(size = 12, family = 'arial'),
              axis.title.x = element_text(size = 14, family = 'arial'),
              axis.text.x = element_text(size = 12, family = 'arial'),
              panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
              legend.position = 'none'
        )
      
      scatterPlotAlnTSIVSSI <- ggplot(mulProxRecordComAln,
                                      aes(x = r,
                                          y = VSSI)) +
        geom_point() +
        geom_smooth(
          method = "lm",
          color = 'black') +
        stat_poly_eq(aes(
          label = paste(after_stat(eq.label), 
                        after_stat(rr.label), 
                        sep = "*\", \"*")), 
          formula = y ~ x, parse = TRUE) +
        labs(x = 'Pearson r (Alnus & Alnus)',
             y = 'VSSI') + 
        theme_classic() + 
        theme(plot.title = element_blank(), 
              axis.title.y = element_text(size = 14, family = 'arial'),
              axis.text.y = element_text(size = 12, family = 'arial'),
              axis.title.x = element_text(size = 14, family = 'arial'),
              axis.text.x = element_text(size = 12, family = 'arial'),
              panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
              legend.position = 'none'
        )
    }
    
    #Calluna & TSI
    {
      mulProxRecordComCal <- mulProxRecordCom[which(mulProxRecordCom[, 1] == 'Calluna & TSI'), ]
      
      scatterPlotCalTSITSI <- ggplot(mulProxRecordComCal,
                                     aes(x = r,
                                         y = TSI)) +
        geom_point() +
        geom_smooth(
          method = "lm",
          color = 'black') +
        stat_poly_eq(aes(
          label = paste(after_stat(eq.label), 
                        after_stat(rr.label), 
                        sep = "*\", \"*")), 
          formula = y ~ x, parse = TRUE) +
        labs(x = 'Pearson r (Calluna & Alnus)',
             y = 'TSI (Wm-2)') + 
        theme_classic() + 
        theme(plot.title = element_blank(), 
              axis.title.y = element_text(size = 14, family = 'arial'),
              axis.text.y = element_text(size = 12, family = 'arial'),
              axis.title.x = element_text(size = 14, family = 'arial'),
              axis.text.x = element_text(size = 12, family = 'arial'),
              panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
              legend.position = 'none'
        )
      
      scatterPlotCalTSIGLTemp <- ggplot(mulProxRecordComCal,
                                        aes(x = r,
                                            y = GLTemp)) +
        geom_point() +
        geom_smooth(
          method = "lm",
          color = 'black') +
        stat_poly_eq(aes(
          label = paste(after_stat(eq.label), 
                        after_stat(rr.label), 
                        sep = "*\", \"*")), 
          formula = y ~ x, parse = TRUE) +
        labs(x = 'Pearson r (Calluna & Alnus)',
             y = 'Greenland temperature anomaly') + 
        theme_classic() + 
        theme(plot.title = element_blank(), 
              axis.title.y = element_text(size = 14, family = 'arial'),
              axis.text.y = element_text(size = 12, family = 'arial'),
              axis.title.x = element_text(size = 14, family = 'arial'),
              axis.text.x = element_text(size = 12, family = 'arial'),
              panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
              legend.position = 'none'
        )
      
      scatterPlotCalTSIWMI <- ggplot(mulProxRecordComCal,
                                     aes(x = r,
                                         y = WMI)) +
        geom_point() +
        geom_smooth(
          method = "lm",
          color = 'black') +
        stat_poly_eq(aes(
          label = paste(after_stat(eq.label), 
                        after_stat(rr.label), 
                        sep = "*\", \"*")), 
          formula = y ~ x, parse = TRUE) +
        labs(x = 'Pearson r (Calluna & Alnus)',
             y = 'Western Mediterranean Index (WMI)') + 
        theme_classic() + 
        theme(plot.title = element_blank(), 
              axis.title.y = element_text(size = 14, family = 'arial'),
              axis.text.y = element_text(size = 12, family = 'arial'),
              axis.title.x = element_text(size = 14, family = 'arial'),
              axis.text.x = element_text(size = 12, family = 'arial'),
              panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
              legend.position = 'none'
        )
      
      scatterPlotCalTSIWMO <- ggplot(mulProxRecordComCal,
                                     aes(x = r,
                                         y = WMO)) +
        geom_point() +
        geom_smooth(
          method = "lm",
          color = 'black') +
        stat_poly_eq(aes(
          label = paste(after_stat(eq.label), 
                        after_stat(rr.label), 
                        sep = "*\", \"*")), 
          formula = y ~ x, parse = TRUE) +
        labs(x = 'Pearson r (Calluna & Alnus)',
             y = 'd18O') + 
        theme_classic() + 
        theme(plot.title = element_blank(), 
              axis.title.y = element_text(size = 14, family = 'arial'),
              axis.text.y = element_text(size = 12, family = 'arial'),
              axis.title.x = element_text(size = 14, family = 'arial'),
              axis.text.x = element_text(size = 12, family = 'arial'),
              panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
              legend.position = 'none'
        )
      
      scatterPlotCalTSIVSSI <- ggplot(mulProxRecordComCal,
                                      aes(x = r,
                                          y = VSSI)) +
        geom_point() +
        geom_smooth(
          method = "lm",
          color = 'black') +
        stat_poly_eq(aes(
          label = paste(after_stat(eq.label), 
                        after_stat(rr.label), 
                        sep = "*\", \"*")), 
          formula = y ~ x, parse = TRUE) +
        labs(x = 'Pearson r (Calluna & Alnus)',
             y = 'VSSI') + 
        theme_classic() + 
        theme(plot.title = element_blank(), 
              axis.title.y = element_text(size = 14, family = 'arial'),
              axis.text.y = element_text(size = 12, family = 'arial'),
              axis.title.x = element_text(size = 14, family = 'arial'),
              axis.text.x = element_text(size = 12, family = 'arial'),
              panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
              legend.position = 'none'
        )
    }
    
    #Sphagnum & d18O
    {
      mulProxRecordComSph <- mulProxRecordCom[which(mulProxRecordCom[, 1] == 'Sphagnum & d18O'), ]
      
      scatterPlotSphd18OTSI <- ggplot(mulProxRecordComSph,
                                     aes(x = r,
                                         y = TSI)) +
        geom_point() +
        geom_smooth(
          method = "lm",
          color = 'black') +
        stat_poly_eq(aes(
          label = paste(after_stat(eq.label), 
                        after_stat(rr.label), 
                        sep = "*\", \"*")), 
          formula = y ~ x, parse = TRUE) +
        labs(x = 'Pearson r (Sphagnum & Alnus)',
             y = 'TSI (Wm-2)') + 
        theme_classic() + 
        theme(plot.title = element_blank(), 
              axis.title.y = element_text(size = 14, family = 'arial'),
              axis.text.y = element_text(size = 12, family = 'arial'),
              axis.title.x = element_text(size = 14, family = 'arial'),
              axis.text.x = element_text(size = 12, family = 'arial'),
              panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
              legend.position = 'none'
        )
      
      scatterPlotSphd18OGLTemp <- ggplot(mulProxRecordComSph,
                                        aes(x = r,
                                            y = GLTemp)) +
        geom_point() +
        geom_smooth(
          method = "lm",
          color = 'black') +
        stat_poly_eq(aes(
          label = paste(after_stat(eq.label), 
                        after_stat(rr.label), 
                        sep = "*\", \"*")), 
          formula = y ~ x, parse = TRUE) +
        labs(x = 'Pearson r (Sphagnum & Alnus)',
             y = 'Greenland temperature anomaly') + 
        theme_classic() + 
        theme(plot.title = element_blank(), 
              axis.title.y = element_text(size = 14, family = 'arial'),
              axis.text.y = element_text(size = 12, family = 'arial'),
              axis.title.x = element_text(size = 14, family = 'arial'),
              axis.text.x = element_text(size = 12, family = 'arial'),
              panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
              legend.position = 'none'
        )
      
      scatterPlotSphd18OWMI <- ggplot(mulProxRecordComSph,
                                     aes(x = r,
                                         y = WMI)) +
        geom_point() +
        geom_smooth(
          method = "lm",
          color = 'black') +
        stat_poly_eq(aes(
          label = paste(after_stat(eq.label), 
                        after_stat(rr.label), 
                        sep = "*\", \"*")), 
          formula = y ~ x, parse = TRUE) +
        labs(x = 'Pearson r (Sphagnum & Alnus)',
             y = 'Western Mediterranean Index (WMI)') + 
        theme_classic() + 
        theme(plot.title = element_blank(), 
              axis.title.y = element_text(size = 14, family = 'arial'),
              axis.text.y = element_text(size = 12, family = 'arial'),
              axis.title.x = element_text(size = 14, family = 'arial'),
              axis.text.x = element_text(size = 12, family = 'arial'),
              panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
              legend.position = 'none'
        )
      
      scatterPlotSphd18OWMO <- ggplot(mulProxRecordComSph,
                                     aes(x = r,
                                         y = WMO)) +
        geom_point() +
        geom_smooth(
          method = "lm",
          color = 'black') +
        stat_poly_eq(aes(
          label = paste(after_stat(eq.label), 
                        after_stat(rr.label), 
                        sep = "*\", \"*")), 
          formula = y ~ x, parse = TRUE) +
        labs(x = 'Pearson r (Sphagnum & Alnus)',
             y = 'd18O') + 
        theme_classic() + 
        theme(plot.title = element_blank(), 
              axis.title.y = element_text(size = 14, family = 'arial'),
              axis.text.y = element_text(size = 12, family = 'arial'),
              axis.title.x = element_text(size = 14, family = 'arial'),
              axis.text.x = element_text(size = 12, family = 'arial'),
              panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
              legend.position = 'none'
        )
      
      scatterPlotSphd18OVSSI <- ggplot(mulProxRecordComSph,
                                      aes(x = r,
                                          y = VSSI)) +
        geom_point() +
        geom_smooth(
          method = "lm",
          color = 'black') +
        stat_poly_eq(aes(
          label = paste(after_stat(eq.label), 
                        after_stat(rr.label), 
                        sep = "*\", \"*")), 
          formula = y ~ x, parse = TRUE) +
        labs(x = 'Pearson r (Sphagnum & Alnus)',
             y = 'VSSI') + 
        theme_classic() + 
        theme(plot.title = element_blank(), 
              axis.title.y = element_text(size = 14, family = 'arial'),
              axis.text.y = element_text(size = 12, family = 'arial'),
              axis.title.x = element_text(size = 14, family = 'arial'),
              axis.text.x = element_text(size = 12, family = 'arial'),
              panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
              legend.position = 'none'
        )
    }
    
    #Alnus & d18O
    {
      mulProxRecordComAln <- mulProxRecordCom[which(mulProxRecordCom[, 1] == 'Alnus & d18O'), ]
      
      scatterPlotAlnd18OTSI <- ggplot(mulProxRecordComAln,
                                      aes(x = r,
                                          y = TSI)) +
        geom_point() +
        geom_smooth(
          method = "lm",
          color = 'black') +
        stat_poly_eq(aes(
          label = paste(after_stat(eq.label), 
                        after_stat(rr.label), 
                        sep = "*\", \"*")), 
          formula = y ~ x, parse = TRUE) +
        labs(x = 'Pearson r (Alnus & Alnus)',
             y = 'TSI (Wm-2)') + 
        theme_classic() + 
        theme(plot.title = element_blank(), 
              axis.title.y = element_text(size = 14, family = 'arial'),
              axis.text.y = element_text(size = 12, family = 'arial'),
              axis.title.x = element_text(size = 14, family = 'arial'),
              axis.text.x = element_text(size = 12, family = 'arial'),
              panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
              legend.position = 'none'
        )
      
      scatterPlotAlnd18OGLTemp <- ggplot(mulProxRecordComAln,
                                         aes(x = r,
                                             y = GLTemp)) +
        geom_point() +
        geom_smooth(
          method = "lm",
          color = 'black') +
        stat_poly_eq(aes(
          label = paste(after_stat(eq.label), 
                        after_stat(rr.label), 
                        sep = "*\", \"*")), 
          formula = y ~ x, parse = TRUE) +
        labs(x = 'Pearson r (Alnus & Alnus)',
             y = 'Greenland temperature anomaly') + 
        theme_classic() + 
        theme(plot.title = element_blank(), 
              axis.title.y = element_text(size = 14, family = 'arial'),
              axis.text.y = element_text(size = 12, family = 'arial'),
              axis.title.x = element_text(size = 14, family = 'arial'),
              axis.text.x = element_text(size = 12, family = 'arial'),
              panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
              legend.position = 'none'
        )
      
      scatterPlotAlnd18OWMI <- ggplot(mulProxRecordComAln,
                                      aes(x = r,
                                          y = WMI)) +
        geom_point() +
        geom_smooth(
          method = "lm",
          color = 'black') +
        stat_poly_eq(aes(
          label = paste(after_stat(eq.label), 
                        after_stat(rr.label), 
                        sep = "*\", \"*")), 
          formula = y ~ x, parse = TRUE) +
        labs(x = 'Pearson r (Alnus & Alnus)',
             y = 'Western Mediterranean Index (WMI)') + 
        theme_classic() + 
        theme(plot.title = element_blank(), 
              axis.title.y = element_text(size = 14, family = 'arial'),
              axis.text.y = element_text(size = 12, family = 'arial'),
              axis.title.x = element_text(size = 14, family = 'arial'),
              axis.text.x = element_text(size = 12, family = 'arial'),
              panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
              legend.position = 'none'
        )
      
      scatterPlotAlnd18OWMO <- ggplot(mulProxRecordComAln,
                                      aes(x = r,
                                          y = WMO)) +
        geom_point() +
        geom_smooth(
          method = "lm",
          color = 'black') +
        stat_poly_eq(aes(
          label = paste(after_stat(eq.label), 
                        after_stat(rr.label), 
                        sep = "*\", \"*")), 
          formula = y ~ x, parse = TRUE) +
        labs(x = 'Pearson r (Alnus & Alnus)',
             y = 'd18O') + 
        theme_classic() + 
        theme(plot.title = element_blank(), 
              axis.title.y = element_text(size = 14, family = 'arial'),
              axis.text.y = element_text(size = 12, family = 'arial'),
              axis.title.x = element_text(size = 14, family = 'arial'),
              axis.text.x = element_text(size = 12, family = 'arial'),
              panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
              legend.position = 'none'
        )
      
      scatterPlotAlnd18OVSSI <- ggplot(mulProxRecordComAln,
                                       aes(x = r,
                                           y = VSSI)) +
        geom_point() +
        geom_smooth(
          method = "lm",
          color = 'black') +
        stat_poly_eq(aes(
          label = paste(after_stat(eq.label), 
                        after_stat(rr.label), 
                        sep = "*\", \"*")), 
          formula = y ~ x, parse = TRUE) +
        labs(x = 'Pearson r (Alnus & Alnus)',
             y = 'VSSI') + 
        theme_classic() + 
        theme(plot.title = element_blank(), 
              axis.title.y = element_text(size = 14, family = 'arial'),
              axis.text.y = element_text(size = 12, family = 'arial'),
              axis.title.x = element_text(size = 14, family = 'arial'),
              axis.text.x = element_text(size = 12, family = 'arial'),
              panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
              legend.position = 'none'
        )
    }
    
    #Calluna & d18O
    {
      mulProxRecordComCal <- mulProxRecordCom[which(mulProxRecordCom[, 1] == 'Calluna & d18O'), ]
      
      scatterPlotCald18OTSI <- ggplot(mulProxRecordComCal,
                                      aes(x = r,
                                          y = TSI)) +
        geom_point() +
        geom_smooth(
          method = "lm",
          color = 'black') +
        stat_poly_eq(aes(
          label = paste(after_stat(eq.label), 
                        after_stat(rr.label), 
                        sep = "*\", \"*")), 
          formula = y ~ x, parse = TRUE) +
        labs(x = 'Pearson r (Calluna & Alnus)',
             y = 'TSI (Wm-2)') + 
        theme_classic() + 
        theme(plot.title = element_blank(), 
              axis.title.y = element_text(size = 14, family = 'arial'),
              axis.text.y = element_text(size = 12, family = 'arial'),
              axis.title.x = element_text(size = 14, family = 'arial'),
              axis.text.x = element_text(size = 12, family = 'arial'),
              panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
              legend.position = 'none'
        )
      
      scatterPlotCald18OGLTemp <- ggplot(mulProxRecordComCal,
                                         aes(x = r,
                                             y = GLTemp)) +
        geom_point() +
        geom_smooth(
          method = "lm",
          color = 'black') +
        stat_poly_eq(aes(
          label = paste(after_stat(eq.label), 
                        after_stat(rr.label), 
                        sep = "*\", \"*")), 
          formula = y ~ x, parse = TRUE) +
        labs(x = 'Pearson r (Calluna & Alnus)',
             y = 'Greenland temperature anomaly') + 
        theme_classic() + 
        theme(plot.title = element_blank(), 
              axis.title.y = element_text(size = 14, family = 'arial'),
              axis.text.y = element_text(size = 12, family = 'arial'),
              axis.title.x = element_text(size = 14, family = 'arial'),
              axis.text.x = element_text(size = 12, family = 'arial'),
              panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
              legend.position = 'none'
        )
      
      scatterPlotCald18OWMI <- ggplot(mulProxRecordComCal,
                                      aes(x = r,
                                          y = WMI)) +
        geom_point() +
        geom_smooth(
          method = "lm",
          color = 'black') +
        stat_poly_eq(aes(
          label = paste(after_stat(eq.label), 
                        after_stat(rr.label), 
                        sep = "*\", \"*")), 
          formula = y ~ x, parse = TRUE) +
        labs(x = 'Pearson r (Calluna & Alnus)',
             y = 'Western Mediterranean Index (WMI)') + 
        theme_classic() + 
        theme(plot.title = element_blank(), 
              axis.title.y = element_text(size = 14, family = 'arial'),
              axis.text.y = element_text(size = 12, family = 'arial'),
              axis.title.x = element_text(size = 14, family = 'arial'),
              axis.text.x = element_text(size = 12, family = 'arial'),
              panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
              legend.position = 'none'
        )
      
      scatterPlotCald18OWMO <- ggplot(mulProxRecordComCal,
                                      aes(x = r,
                                          y = WMO)) +
        geom_point() +
        geom_smooth(
          method = "lm",
          color = 'black') +
        stat_poly_eq(aes(
          label = paste(after_stat(eq.label), 
                        after_stat(rr.label), 
                        sep = "*\", \"*")), 
          formula = y ~ x, parse = TRUE) +
        labs(x = 'Pearson r (Calluna & Alnus)',
             y = 'd18O') + 
        theme_classic() + 
        theme(plot.title = element_blank(), 
              axis.title.y = element_text(size = 14, family = 'arial'),
              axis.text.y = element_text(size = 12, family = 'arial'),
              axis.title.x = element_text(size = 14, family = 'arial'),
              axis.text.x = element_text(size = 12, family = 'arial'),
              panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
              legend.position = 'none'
        )
      
      scatterPlotCald18OVSSI <- ggplot(mulProxRecordComCal,
                                       aes(x = r,
                                           y = VSSI)) +
        geom_point() +
        geom_smooth(
          method = "lm",
          color = 'black') +
        stat_poly_eq(aes(
          label = paste(after_stat(eq.label), 
                        after_stat(rr.label), 
                        sep = "*\", \"*")), 
          formula = y ~ x, parse = TRUE) +
        labs(x = 'Pearson r (Calluna & Alnus)',
             y = 'VSSI') + 
        theme_classic() + 
        theme(plot.title = element_blank(), 
              axis.title.y = element_text(size = 14, family = 'arial'),
              axis.text.y = element_text(size = 12, family = 'arial'),
              axis.title.x = element_text(size = 14, family = 'arial'),
              axis.text.x = element_text(size = 12, family = 'arial'),
              panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
              legend.position = 'none'
        )
    }
    
    #Arrange on one page (Calluna & TSI vs d18O)
    {
      (scatterPlotSphTSITSI | scatterPlotSphTSIGLTemp) /
        (scatterPlotSphTSIWMI | scatterPlotSphTSIWMI) /
        (scatterPlotSphTSIWMO | scatterPlotSphTSIVSSI)
      
      (scatterPlotAlnTSITSI | scatterPlotAlnTSIGLTemp) /
        (scatterPlotAlnTSIWMI | scatterPlotAlnTSIWMI) /
        (scatterPlotAlnTSIWMO | scatterPlotAlnTSIVSSI)
      
      (scatterPlotCalTSITSI | scatterPlotCalTSIGLTemp) /
        (scatterPlotCalTSIWMI | scatterPlotCalTSIWMI) /
        (scatterPlotCalTSIWMO | scatterPlotCalTSIVSSI)
      
      (scatterPlotSphd18OTSI | scatterPlotSphd18OGLTemp) /
        (scatterPlotSphd18OWMI | scatterPlotSphd18OWMI) /
        (scatterPlotSphd18OWMO | scatterPlotSphd18OVSSI)
      
      (scatterPlotAlnTSITSI | scatterPlotAlnTSIGLTemp) /
        (scatterPlotAlnd18OWMI | scatterPlotAlnd18OWMI) /
        (scatterPlotAlnd18OWMO | scatterPlotAlnd18OVSSI)
      
      (scatterPlotCald18OTSI | scatterPlotCald18OGLTemp) /
        (scatterPlotCald18OWMI | scatterPlotCald18OWMI) /
        (scatterPlotCald18OWMO | scatterPlotCald18OVSSI)
    }
  }
}



