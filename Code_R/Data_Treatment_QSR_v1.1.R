#This script aims at reproducing all figures already made for NC manuscript,
  #including those that are going to be further added,
  #while two more taxa (Alnus and Calluna) are considered

#Load packages and declare functions
{
  library(here)
  
  #Plot
  library(ggplot2)
  library(tidyverse)
  library(tidypaleo)
  library(dplyr)
  
  #Combine plot
  library(ggpubr)
  
  #Moving average calculation
  library(zoo)
  
  #Peak statistics 
  library(pavo)
  
  #Prevent axis from intersection
  library(ggh4x)
  
  #Correlation
  library(PerformanceAnalytics)
  library(Hmisc)
  
  #Load arial font into R
  library(showtext)
  font_add(family = "arial", regular = here('Fonts', 'arial.ttf'))
  showtext_auto()
  
  #Guassian smoothing
  library(smoother)
  
  library(patchwork)
  
  #Correlation analysis from Reschke et al., 2019
  library('corit')
  
  library(reshape2)
  
  #load package for Wavelet analysis
  library(WaveletComp)
  
  #load package for Batlow scientific colour rainbow (Crameri, 2020, Nature Communication)
  library(scico)
  
  #geoChronR packages
  library(geoChronR)
  library(magrittr)
  
  #Age uncertainty
  library(lipdR)
  
  #Paralell processing
  library(doSNOW)
  
  #PCA
  library(ordr)
  library(ggallin)
  
  #Bayesian regression
  library(mlbench)
  library(rstanarm)
  library(bayestestR)
  library(bayesplot)
  library(insight)
  library(broom)
  
  #Triangle diagram
  library(ggtern)
  library(plotly)
  library(readr)
  library(dplyr)
  library(tidyr)
}

{
  #Declare a function to list all peak statistics for given peak
  peakAreaTestStatis <- function(specMat, peakRangeMat, plotDec, PlotIndex) {
    #Get the wavenumber resolution of current spectra matrix
    waveReso <- as.numeric(colnames(specMat)[1]) - as.numeric(colnames(specMat)[2])
    
    #Check the number of peak ranges input from the matrix
    if(is.null(nrow(peakRangeMat))) {
      #Get the index of column for the upper and lower limit of peak range
      indexPAU <- (as.numeric(colnames(specMat)[1]) - peakRangeMat[1]) / waveReso + 1
      
      indexPAL <- (as.numeric(colnames(specMat)[1]) - peakRangeMat[2]) / waveReso + 1
      
      #Clip the peak area according to the range input
      tempMatClip <- specMat[, indexPAU:indexPAL]
      
      #Create matrix to store peak area under different definitions
      peakARActSquYellow <- matrix(nrow = ncol(tempMatClip), 
                                   ncol = nrow(tempMatClip) + 1,
                                   dimnames = list(c(), c('wl', 
                                                          1:nrow(tempMatClip))))
      
      peakARActSquYellow[, 1] <- as.numeric(colnames(tempMatClip))
      
      #Loop to calculate the peak area, under both definitions.
      loopI <- 1
      while(loopI <= nrow(specMat)){
        #Create an artificial baseline for the peak individually, by using linear
        #interpolation between the ridges of peak range. Baseline values no higher
        #than original absorbance intensities are kept, while replaced by corresponding
        #absorbance.
        tempBaselineArti <- approx(as.numeric(colnames(tempMatClip)[c(1, ncol(tempMatClip))]), 
                                   tempMatClip[loopI, c(1, ncol(tempMatClip))],
                                   xout = as.numeric(colnames(tempMatClip)))[[2]]
        tempBaselineMat <- matrix(data = NA, nrow = 2, ncol = length(tempBaselineArti))
        tempBaselineMat[1, ] <- tempBaselineArti
        tempBaselineMat[2, ] <- tempMatClip[loopI, ]
        
        #Get peak area under yellow region
        peakARActSquYellow[, loopI + 1] <-
          tempBaselineMat[2, ] - apply(tempBaselineMat, 2, min)
        
        loopI <- loopI + 1
      }
      
      #Get the statistics for current peak
      peakARActSquYellow <- as.data.frame(peakARActSquYellow)
      peakStatis <- peakshape(peakARActSquYellow, 
                              plot = plotDec,
                              select = PlotIndex)
      
      #Output the peak area following the same sequence of raw spectra matrix,
      #under two definitions
      return(peakStatis)
    }
  }
  
  #Declare a function to calculate the derivatives of a given spectrum
  derivaSpectrum <- function(specVec, degree, step) {
    
    tempvec <- specVec
    
    #Loop to calculate the derivative
    loopI <- 1
    while (loopI <= degree) {
      
      tempVecLast <- tempvec[-1]
      
      tempVecFirst <- tempvec[-length(tempvec)]
      
      tempvec <- (tempVecLast - tempVecFirst)
      
      loopI <- loopI + 1
    }
    
    return(tempvec)
  }
  
  #Ben's peak calculation functions
  {
    # Load peak detection script
    source("https://raw.githubusercontent.com/benbell95/peak-detection/main/r/peak_detection.r")
    # Load peak area script
    source("https://raw.githubusercontent.com/benbell95/peak-detection/main/r/peak_area.r")
  }
}

#Read UAC signals of different approaches, for three taxa
{
  #Prepare chronology information for the read spectra dataset
  {
    #Chronology
    chronologyHM20 <- apply(as.matrix(read.csv(here('chronologyBchron.csv'))), 
                            c(1,2), as.numeric)
    colnames(chronologyHM20)[1:2] <- c('Years', 'Depth')
    
    #Read transfering functions between IntCal20 and GICC05
    chronologyHM20TransferInte <- apply(read.csv(here('chronologyHM20TransferInte.csv'),
                                                 row.names = 1), c(1,2), as.numeric)
  }
  
  #Sphagnum
  {
    #UAC signals
    UACSph <- read.csv(here('Output_QSR', 'Data_Publishment', 'UAC_Signals_Sph.csv'),
                       header = TRUE, row.names = 1)
    
    #Mean and SD of approaches
    UACSphMean <- apply(UACSph[, c(3,7,9,11)], 1, mean)
    UACSphSD <- apply(UACSph[, c(3,7,9,11)], 1, sd)
    
    #Chronology correction between raidocarbon and ice chronology
    UACAreaMatMeanSphCorr <- c()
    loopI <- 1
    while (loopI <= length(UACSph[, 2])) {
      
      UACAreaMatMeanSphCorr[length(UACAreaMatMeanSphCorr) + 1] <-
        chronologyHM20TransferInte[which(chronologyHM20TransferInte[, 1] ==
                                           UACSph[loopI, 2]), 2]
      
      loopI <- loopI + 1
    }
    rm(loopI)
    
    #Gaussian filter on mean vaules
    UACSphGau <- smth.gaussian(UACSphMean, window = 5)
    
    UACSphGauTime <- UACSph[, 2] - UACAreaMatMeanSphCorr
    
    #Get the mean and SD of smoothed values
    UACSphGauMean <- mean(UACSphGau[which(!is.na(UACSphGau))])
    
    UACSphGauSD <- sd(UACSphGau[which(!is.na(UACSphGau))])
  }

  #Alnus
  {
    #UAC signals
    UACAln <- read.csv(here('Output_QSR', 'Data_Publishment', 'UAC_Signals_Aln.csv'),
                       header = TRUE, row.names = 1)
    
    #Mean and SD of approaches
    UACAlnMean <- apply(UACAln[, c(3,7,9,11)], 1, mean)
    UACAlnSD <- apply(UACAln[, c(3,7,9,11)], 1, sd)
    
    #Chronology correction between raidocarbon and ice chronology
    UACAreaMatMeanAlnCorr <- c()
    loopI <- 1
    while (loopI <= length(UACAln[, 2])) {
      
      UACAreaMatMeanAlnCorr[length(UACAreaMatMeanAlnCorr) + 1] <-
        chronologyHM20TransferInte[which(chronologyHM20TransferInte[, 1] ==
                                           UACAln[loopI, 2]), 2]
      
      loopI <- loopI + 1
    }
    rm(loopI)
    
    #Gaussian filter on mean vaules
    UACAlnGau <- smth.gaussian(UACAlnMean, window = 5)
    
    UACAlnGauTime <- UACAln[, 2] - UACAreaMatMeanAlnCorr
    
    #Get the mean and SD of smoothed values
    UACAlnGauMean <- mean(UACAlnGau[which(!is.na(UACAlnGau))])
    
    UACAlnGauSD <- sd(UACAlnGau[which(!is.na(UACAlnGau))])
  }
  
  #Calluna
  {
    #UAC signals
    UACCal <- read.csv(here('Output_QSR', 'Data_Publishment', 'UAC_Signals_Cal.csv'),
                       header = TRUE, row.names = 1)
    
    #Mean and SD of approaches
    UACCalMean <- apply(UACCal[, c(3,7,9,11)], 1, mean)
    UACCalSD <- apply(UACCal[, c(3,7,9,11)], 1, sd)
    
    #Chronology correction between raidocarbon and ice chronology
    UACAreaMatMeanCalCorr <- c()
    loopI <- 1
    while (loopI <= length(UACCal[, 2])) {
      
      UACAreaMatMeanCalCorr[length(UACAreaMatMeanCalCorr) + 1] <-
        chronologyHM20TransferInte[which(chronologyHM20TransferInte[, 1] ==
                                           UACCal[loopI, 2]), 2]
      
      loopI <- loopI + 1
    }
    rm(loopI)
    
    #Gaussian filter on mean vaules
    UACCalGau <- smth.gaussian(UACCalMean, window = 5)
    
    UACCalGauTime <- UACCal[, 2] - UACAreaMatMeanCalCorr
    
    #Get the mean and SD of smoothed values
    UACCalGauMean <- mean(UACCalGau[which(!is.na(UACCalGau))])
    
    UACCalGauSD <- sd(UACCalGau[which(!is.na(UACCalGau))])
  }
  
  #Direct plot of UAC against depth
  {
    #Sphagnum
    {
      #Scaled dataset
      {
        localPeakHeightPlotSphSca <- ggplot() +
          geom_point( 
            aes(x = UACSph[, 3],
                y = UACSph[, 1] * 2), 
            size = 2, shape = 21, fill = 'red', alpha = 0.2) + 
          geom_lineh(aes(x = smth.gaussian(UACSph[, 3], window = 5),
                         y = UACSph[, 1] * 2), linewidth = 1) +
          geom_errorbarh(
            aes(y = UACSph[, 1] * 2,
                xmin= UACSph[, 3] - UACSph[, 4], 
                xmax= UACSph[, 3] + UACSph[, 4]), 
            size = 1, alpha = 0.2
          ) +
          scale_x_continuous(breaks = seq(from = -10, to = 10, by = 1),
                             limits = range(UACSph[, 3] - UACSph[, 4],
                                            UACSph[, 3] + UACSph[, 4]),
                             position = 'top') +
          scale_y_reverse(limits = c(370, 0), 
                          breaks = seq(from = 370, to = 0, by = -10),
                          labels = c(rep('', 2),
                                     350,rep('', 4),
                                     300, rep('', 4),
                                     250,rep('', 4),
                                     200, rep('', 4),
                                     150,rep('', 4),
                                     100, rep('', 4),
                                     50,rep('', 4),
                                     0)
          )  +
          guides(x = "axis_truncated",y = "axis_truncated") +
          labs(y = 'Depth (cm)', x = 'Peak - Min', title = '') + 
          theme_classic() + 
          theme(plot.title = element_blank(), 
                axis.text.y = element_text(size = 10, family = 'arial'),
                axis.title.y = element_text(size = 12, family = 'arial'),
                axis.text.x = element_text(size = 10, family = 'arial'),
                axis.title.x = element_text(size = 12, family = 'arial')
          )
        
        ratioUACOHPlotSphSca <- ggplot() +
          geom_point( 
            aes(x = UACSph[, 5],
                y = UACSph[, 1] * 2), 
            size = 2, shape = 21, fill = 'red', alpha = 0.2) + 
          geom_lineh(aes(x = smth.gaussian(UACSph[, 5], window = 5),
                         y = UACSph[, 1] * 2), linewidth = 1) +
          geom_errorbarh(
            aes(y = UACSph[, 1] * 2,
                xmin= UACSph[, 5] - UACSph[, 6], 
                xmax= UACSph[, 5] + UACSph[, 6]), 
            size = 1, alpha = 0.2
          ) +
          scale_x_continuous(breaks = seq(from = -10, to = 10, by = 1),
                             limits = range(UACSph[, 5] - UACSph[, 6],
                                            UACSph[, 5] + UACSph[, 6]),
                             position = 'top') +
          scale_y_reverse(limits = c(370, 0), 
                          breaks = seq(from = 370, to = 0, by = -10),
                          labels = c(rep('', 2),
                                     350,rep('', 4),
                                     300, rep('', 4),
                                     250,rep('', 4),
                                     200, rep('', 4),
                                     150,rep('', 4),
                                     100, rep('', 4),
                                     50,rep('', 4),
                                     0)
          )  +
          guides(x = "axis_truncated", y = "axis_truncated") +
          labs(y = 'Depth (cm)', x = 'Aro-OH Ratio', title = '') + 
          theme_classic() + 
          theme(plot.title = element_blank(), 
                axis.text.y = element_text(size = 10, family = 'arial'),
                axis.title.y = element_text(size = 12, family = 'arial'),
                axis.text.x = element_text(size = 10, family = 'arial'),
                axis.title.x = element_text(size = 12, family = 'arial')
          )
        
        ratio2ndPlotSphSca <- ggplot() +
          geom_point( 
            aes(x = UACSph[, 7],
                y = UACSph[, 1] * 2), 
            size = 2, shape = 21, fill = 'red', alpha = 0.2) + 
          geom_lineh(aes(x = smth.gaussian(UACSph[, 7], window = 5),
                         y = UACSph[, 1] * 2), linewidth = 1) +
          geom_errorbarh(
            aes(y = UACSph[, 1] * 2,
                xmin= UACSph[, 7] - UACSph[, 8], 
                xmax= UACSph[, 7] + UACSph[, 8]), 
            size = 1, alpha = 0.2
          ) +
          scale_x_continuous(breaks = seq(from = -10, to = 10, by = 1),
                             limits = range(UACSph[, 7] - UACSph[, 8],
                                            UACSph[, 7] + UACSph[, 8]),
                             position = 'top') +
          scale_y_reverse(limits = c(370, 0), 
                          breaks = seq(from = 370, to = 0, by = -10),
                          labels = c(rep('', 2),
                                     350,rep('', 4),
                                     300, rep('', 4),
                                     250,rep('', 4),
                                     200, rep('', 4),
                                     150,rep('', 4),
                                     100, rep('', 4),
                                     50,rep('', 4),
                                     0)
          )  +
          guides(x = "axis_truncated",y = "axis_truncated") +
          labs(x = '2nd derivative peak ratio', y = 'Cal years BP', title = '') + 
          theme_classic() + 
          theme(plot.title = element_blank(), 
                axis.text.y = element_text(size = 10, family = 'arial'),
                axis.title.y = element_text(size = 12, family = 'arial'),
                axis.text.x = element_text(size = 10, family = 'arial'),
                axis.title.x = element_text(size = 12, family = 'arial')
          )
        
        local2ndPlotSphSca <- ggplot() +
          geom_point( 
            aes(x = UACSph[, 9],
                y = UACSph[, 1] * 2), 
            size = 2, shape = 21, fill = 'red', alpha = 0.2) + 
          geom_lineh(aes(x = smth.gaussian(UACSph[, 9], window = 5),
                         y = UACSph[, 1] * 2), linewidth = 1) +
          geom_errorbarh(
            aes(y = UACSph[, 1] * 2,
                xmin= UACSph[, 9] - UACSph[, 10], 
                xmax= UACSph[, 9] + UACSph[, 10]), 
            size = 1, alpha = 0.2
          ) +
          scale_x_continuous(breaks = seq(from = -10, to = 10, by = 1),
                             limits = range(UACSph[, 9] - UACSph[, 10],
                                            UACSph[, 9] + UACSph[, 10]),
                             position = 'top') +
          scale_y_reverse(limits = c(370, 0), 
                          breaks = seq(from = 370, to = 0, by = -10),
                          labels = c(rep('', 2),
                                     350,rep('', 4),
                                     300, rep('', 4),
                                     250,rep('', 4),
                                     200, rep('', 4),
                                     150,rep('', 4),
                                     100, rep('', 4),
                                     50,rep('', 4),
                                     0)
          )  +
          guides(x = "axis_truncated",y = "axis_truncated") +
          labs(x = '2nd derivative peak height', y = 'Cal years BP', title = '') + 
          theme_classic() + 
          theme(plot.title = element_blank(), 
                axis.text.y = element_text(size = 10, family = 'arial'),
                axis.title.y = element_text(size = 12, family = 'arial'),
                axis.text.x = element_text(size = 10, family = 'arial'),
                axis.title.x = element_text(size = 12, family = 'arial')
          )
        
        peakAreaPlotSphSca <- ggplot() +
          geom_point( 
            aes(x = UACSph[, 11],
                y = UACSph[, 1] * 2), 
            size = 2, shape = 21, fill = 'red', alpha = 0.2) + 
          geom_lineh(aes(x = smth.gaussian(UACSph[, 11], window = 5),
                         y = UACSph[, 1] * 2), linewidth = 1) +
          geom_errorbarh(
            aes(y = UACSph[, 1] * 2,
                xmin= UACSph[, 11] - UACSph[, 12], 
                xmax= UACSph[, 11] + UACSph[, 12]), 
            size = 1, alpha = 0.2
          ) +
          scale_x_continuous(breaks = seq(from = -10, to = 10, by = 1),
                             limits = range(UACSph[, 11] - UACSph[, 12],
                                            UACSph[, 11] + UACSph[, 12]),
                             position = 'top') +
          scale_y_reverse(limits = c(370, 0), 
                          breaks = seq(from = 370, to = 0, by = -10),
                          labels = c(rep('', 2),
                                     350,rep('', 4),
                                     300, rep('', 4),
                                     250,rep('', 4),
                                     200, rep('', 4),
                                     150,rep('', 4),
                                     100, rep('', 4),
                                     50,rep('', 4),
                                     0)
          )  +
          guides(x = "axis_truncated",y = "axis_truncated") +
          labs(x = 'Peak Area', y = 'Cal years BP', title = '') + 
          theme_classic() + 
          theme(plot.title = element_blank(), 
                axis.text.y = element_text(size = 10, family = 'arial'),
                axis.title.y = element_text(size = 12, family = 'arial'),
                axis.text.x = element_text(size = 10, family = 'arial'),
                axis.title.x = element_text(size = 12, family = 'arial')
          )
        
        meanPlotSphSca <- ggplot() +
          geom_point(
            aes(y = UACSph[, 1] * 2,
                x = UACSphMean), 
            size = 2, shape = 21, fill = 'red', alpha = 0.4) + 
          geom_errorbarh(aes(y = UACSph[, 1] * 2,
                             xmin = UACSphMean - UACSphSD,
                             xmax = UACSphMean + UACSphSD),
                         alpha = 0.4) +
          geom_lineh(aes(x = UACSphGau,
                         y = UACSph[, 1] * 2), size = 1) +
          scale_x_continuous(breaks = seq(from = -10, to = 10, by = 1),
                             limits = range(UACSphMean - UACSphSD,
                                            UACSphMean + UACSphSD),
                             position = 'top') +
          scale_y_reverse(limits = c(370, 0), 
                          breaks = seq(from = 370, to = 0, by = -10),
                          labels = c(rep('', 2),
                                     350,rep('', 4),
                                     300, rep('', 4),
                                     250,rep('', 4),
                                     200, rep('', 4),
                                     150,rep('', 4),
                                     100, rep('', 4),
                                     50,rep('', 4),
                                     0),
                          position = 'right'
          )  +
          guides(x = "axis_truncated", y = "axis_truncated") +
          labs(x = 'Mean of different approaches', y = 'Cal years BP', title = '') + 
          theme_classic() + 
          theme(plot.title = element_blank(), 
                axis.text.y = element_text(size = 10, family = 'arial'),
                axis.title.y = element_text(size = 12, family = 'arial'),
                axis.text.x = element_text(size = 10, family = 'arial'),
                axis.title.x = element_text(size = 12, family = 'arial')
          )
      }
      
      #Combine plot
      {
        ggarrange(ratioUACOHPlotSphSca,
                  localPeakHeightPlotSphSca +
                    theme(
                      axis.line.y = element_blank(),
                      axis.ticks.y = element_blank(),
                      axis.title.y = element_blank(),
                      axis.text.y = element_blank(),
                    ), 
                  peakAreaPlotSphSca +
                    theme(
                      axis.line.y = element_blank(),
                      axis.ticks.y = element_blank(),
                      axis.title.y = element_blank(),
                      axis.text.y = element_blank(),
                    ), 
                  ratio2ndPlotSphSca +
                    
                    theme(
                      axis.line.y = element_blank(),
                      axis.ticks.y = element_blank(),
                      axis.title.y = element_blank(),
                      axis.text.y = element_blank(),
                    ),
                  local2ndPlotSphSca +
                    theme(
                      axis.line.y = element_blank(),
                      axis.ticks.y = element_blank(),
                      axis.title.y = element_blank(),
                      axis.text.y = element_blank(),
                    ),
                  meanPlotSphSca,
                  heights = c(1,1,1,1,1,1),
                  widths = c(1.2,1,1,1,1,1.2),
                  ncol = 6, nrow = 1,
                  align = "h")
      }
    }
    
    #Alnus
    {
      #Scaled dataset
      {
        localPeakHeightPlotAlnSca <- ggplot() +
          geom_point( 
            aes(x = UACAln[, 3],
                y = UACAln[, 1] * 2), 
            size = 2, shape = 21, fill = 'red', alpha = 0.2) + 
          geom_lineh(aes(x = smth.gaussian(UACAln[, 3], window = 5),
                         y = UACAln[, 1] * 2), linewidth = 1) +
          geom_errorbarh(
            aes(y = UACAln[, 1] * 2,
                xmin= UACAln[, 3] - UACAln[, 4], 
                xmax= UACAln[, 3] + UACAln[, 4]), 
            size = 1, alpha = 0.2
          ) +
          scale_x_continuous(breaks = seq(from = -10, to = 10, by = 1),
                             limits = range(UACAln[, 3] - UACAln[, 4],
                                            UACAln[, 3] + UACAln[, 4]),
                             position = 'top') +
          scale_y_reverse(limits = c(370, 0), 
                          breaks = seq(from = 370, to = 0, by = -10),
                          labels = c(rep('', 2),
                                     350,rep('', 4),
                                     300, rep('', 4),
                                     250,rep('', 4),
                                     200, rep('', 4),
                                     150,rep('', 4),
                                     100, rep('', 4),
                                     50,rep('', 4),
                                     0)
          )  +
          guides(x = "axis_truncated",y = "axis_truncated") +
          labs(y = 'Depth (cm)', x = 'Peak - Min', title = '') + 
          theme_classic() + 
          theme(plot.title = element_blank(), 
                axis.text.y = element_text(size = 10, family = 'arial'),
                axis.title.y = element_text(size = 12, family = 'arial'),
                axis.text.x = element_text(size = 10, family = 'arial'),
                axis.title.x = element_text(size = 12, family = 'arial')
          )
        
        ratioUACOHPlotAlnSca <- ggplot() +
          geom_point( 
            aes(x = UACAln[, 5],
                y = UACAln[, 1] * 2), 
            size = 2, shape = 21, fill = 'red', alpha = 0.2) + 
          geom_lineh(aes(x = smth.gaussian(UACAln[, 5], window = 5),
                         y = UACAln[, 1] * 2), linewidth = 1) +
          geom_errorbarh(
            aes(y = UACAln[, 1] * 2,
                xmin= UACAln[, 5] - UACAln[, 6], 
                xmax= UACAln[, 5] + UACAln[, 6]), 
            size = 1, alpha = 0.2
          ) +
          scale_x_continuous(breaks = seq(from = -10, to = 10, by = 1),
                             limits = range(UACAln[, 5] - UACAln[, 6],
                                            UACAln[, 5] + UACAln[, 6]),
                             position = 'top') +
          scale_y_reverse(limits = c(370, 0), 
                          breaks = seq(from = 370, to = 0, by = -10),
                          labels = c(rep('', 2),
                                     350,rep('', 4),
                                     300, rep('', 4),
                                     250,rep('', 4),
                                     200, rep('', 4),
                                     150,rep('', 4),
                                     100, rep('', 4),
                                     50,rep('', 4),
                                     0)
          )  +
          guides(x = "axis_truncated", y = "axis_truncated") +
          labs(y = 'Depth (cm)', x = 'Aro-OH Ratio', title = '') + 
          theme_classic() + 
          theme(plot.title = element_blank(), 
                axis.text.y = element_text(size = 10, family = 'arial'),
                axis.title.y = element_text(size = 12, family = 'arial'),
                axis.text.x = element_text(size = 10, family = 'arial'),
                axis.title.x = element_text(size = 12, family = 'arial')
          )
        
        ratio2ndPlotAlnSca <- ggplot() +
          geom_point( 
            aes(x = UACAln[, 7],
                y = UACAln[, 1] * 2), 
            size = 2, shape = 21, fill = 'red', alpha = 0.2) + 
          geom_lineh(aes(x = smth.gaussian(UACAln[, 7], window = 5),
                         y = UACAln[, 1] * 2), linewidth = 1) +
          geom_errorbarh(
            aes(y = UACAln[, 1] * 2,
                xmin= UACAln[, 7] - UACAln[, 8], 
                xmax= UACAln[, 7] + UACAln[, 8]), 
            size = 1, alpha = 0.2
          ) +
          scale_x_continuous(breaks = seq(from = -10, to = 10, by = 1),
                             limits = range(UACAln[, 7] - UACAln[, 8],
                                            UACAln[, 7] + UACAln[, 8]),
                             position = 'top') +
          scale_y_reverse(limits = c(370, 0), 
                          breaks = seq(from = 370, to = 0, by = -10),
                          labels = c(rep('', 2),
                                     350,rep('', 4),
                                     300, rep('', 4),
                                     250,rep('', 4),
                                     200, rep('', 4),
                                     150,rep('', 4),
                                     100, rep('', 4),
                                     50,rep('', 4),
                                     0)
          )  +
          guides(x = "axis_truncated",y = "axis_truncated") +
          labs(x = '2nd derivative peak ratio', y = 'Cal years BP', title = '') + 
          theme_classic() + 
          theme(plot.title = element_blank(), 
                axis.text.y = element_text(size = 10, family = 'arial'),
                axis.title.y = element_text(size = 12, family = 'arial'),
                axis.text.x = element_text(size = 10, family = 'arial'),
                axis.title.x = element_text(size = 12, family = 'arial')
          )
        
        local2ndPlotAlnSca <- ggplot() +
          geom_point( 
            aes(x = UACAln[, 9],
                y = UACAln[, 1] * 2), 
            size = 2, shape = 21, fill = 'red', alpha = 0.2) + 
          geom_lineh(aes(x = smth.gaussian(UACAln[, 9], window = 5),
                         y = UACAln[, 1] * 2), linewidth = 1) +
          geom_errorbarh(
            aes(y = UACAln[, 1] * 2,
                xmin= UACAln[, 9] - UACAln[, 10], 
                xmax= UACAln[, 9] + UACAln[, 10]), 
            size = 1, alpha = 0.2
          ) +
          scale_x_continuous(breaks = seq(from = -10, to = 10, by = 1),
                             limits = range(UACAln[, 9] - UACAln[, 10],
                                            UACAln[, 9] + UACAln[, 10]),
                             position = 'top') +
          scale_y_reverse(limits = c(370, 0), 
                          breaks = seq(from = 370, to = 0, by = -10),
                          labels = c(rep('', 2),
                                     350,rep('', 4),
                                     300, rep('', 4),
                                     250,rep('', 4),
                                     200, rep('', 4),
                                     150,rep('', 4),
                                     100, rep('', 4),
                                     50,rep('', 4),
                                     0)
          )  +
          guides(x = "axis_truncated",y = "axis_truncated") +
          labs(x = '2nd derivative peak height', y = 'Cal years BP', title = '') + 
          theme_classic() + 
          theme(plot.title = element_blank(), 
                axis.text.y = element_text(size = 10, family = 'arial'),
                axis.title.y = element_text(size = 12, family = 'arial'),
                axis.text.x = element_text(size = 10, family = 'arial'),
                axis.title.x = element_text(size = 12, family = 'arial')
          )
        
        peakAreaPlotAlnSca <- ggplot() +
          geom_point( 
            aes(x = UACAln[, 11],
                y = UACAln[, 1] * 2), 
            size = 2, shape = 21, fill = 'red', alpha = 0.2) + 
          geom_lineh(aes(x = smth.gaussian(UACAln[, 11], window = 5),
                         y = UACAln[, 1] * 2), linewidth = 1) +
          geom_errorbarh(
            aes(y = UACAln[, 1] * 2,
                xmin= UACAln[, 11] - UACAln[, 12], 
                xmax= UACAln[, 11] + UACAln[, 12]), 
            size = 1, alpha = 0.2
          ) +
          scale_x_continuous(breaks = seq(from = -10, to = 10, by = 1),
                             limits = range(UACAln[, 11] - UACAln[, 12],
                                            UACAln[, 11] + UACAln[, 12]),
                             position = 'top') +
          scale_y_reverse(limits = c(370, 0), 
                          breaks = seq(from = 370, to = 0, by = -10),
                          labels = c(rep('', 2),
                                     350,rep('', 4),
                                     300, rep('', 4),
                                     250,rep('', 4),
                                     200, rep('', 4),
                                     150,rep('', 4),
                                     100, rep('', 4),
                                     50,rep('', 4),
                                     0)
          )  +
          guides(x = "axis_truncated",y = "axis_truncated") +
          labs(x = 'Peak Area', y = 'Cal years BP', title = '') + 
          theme_classic() + 
          theme(plot.title = element_blank(), 
                axis.text.y = element_text(size = 10, family = 'arial'),
                axis.title.y = element_text(size = 12, family = 'arial'),
                axis.text.x = element_text(size = 10, family = 'arial'),
                axis.title.x = element_text(size = 12, family = 'arial')
          )
        
        meanPlotAlnSca <- ggplot() +
          geom_point(
            aes(y = UACAln[, 1] * 2,
                x = UACAlnMean), 
            size = 2, shape = 21, fill = 'red', alpha = 0.4) + 
          geom_errorbarh(aes(y = UACAln[, 1] * 2,
                             xmin = UACAlnMean - UACAlnSD,
                             xmax = UACAlnMean + UACAlnSD),
                         alpha = 0.4) +
          geom_lineh(aes(x = UACAlnGau,
                         y = UACAln[, 1] * 2), size = 1) +
          scale_x_continuous(breaks = seq(from = -10, to = 10, by = 1),
                             limits = range(UACAlnMean - UACAlnSD,
                                            UACAlnMean + UACAlnSD),
                             position = 'top') +
          scale_y_reverse(limits = c(370, 0), 
                          breaks = seq(from = 370, to = 0, by = -10),
                          labels = c(rep('', 2),
                                     350,rep('', 4),
                                     300, rep('', 4),
                                     250,rep('', 4),
                                     200, rep('', 4),
                                     150,rep('', 4),
                                     100, rep('', 4),
                                     50,rep('', 4),
                                     0),
                          position = 'right'
          )  +
          guides(x = "axis_truncated", y = "axis_truncated") +
          labs(x = 'Mean of different approaches', y = 'Cal years BP', title = '') + 
          theme_classic() + 
          theme(plot.title = element_blank(), 
                axis.text.y = element_text(size = 10, family = 'arial'),
                axis.title.y = element_text(size = 12, family = 'arial'),
                axis.text.x = element_text(size = 10, family = 'arial'),
                axis.title.x = element_text(size = 12, family = 'arial')
          )
      }
      
      #Combine plot
      {
        ggarrange(ratioUACOHPlotAlnSca,
                  localPeakHeightPlotAlnSca +
                    theme(
                      axis.line.y = element_blank(),
                      axis.ticks.y = element_blank(),
                      axis.title.y = element_blank(),
                      axis.text.y = element_blank(),
                    ), 
                  peakAreaPlotAlnSca +
                    theme(
                      axis.line.y = element_blank(),
                      axis.ticks.y = element_blank(),
                      axis.title.y = element_blank(),
                      axis.text.y = element_blank(),
                    ), 
                  ratio2ndPlotAlnSca +
                    
                    theme(
                      axis.line.y = element_blank(),
                      axis.ticks.y = element_blank(),
                      axis.title.y = element_blank(),
                      axis.text.y = element_blank(),
                    ),
                  local2ndPlotAlnSca +
                    theme(
                      axis.line.y = element_blank(),
                      axis.ticks.y = element_blank(),
                      axis.title.y = element_blank(),
                      axis.text.y = element_blank(),
                    ),
                  meanPlotAlnSca,
                  heights = c(1,1,1,1,1,1),
                  widths = c(1.2,1,1,1,1,1.2),
                  ncol = 6, nrow = 1,
                  align = "h")
      }
    }

    #Calluna
    {
      #Scaled dataset
      {
        localPeakHeightPlotCalSca <- ggplot() +
          geom_point( 
            aes(x = UACCal[, 3],
                y = UACCal[, 1] * 2), 
            size = 2, shape = 21, fill = 'red', alpha = 0.2) + 
          geom_lineh(aes(x = smth.gaussian(UACCal[, 3], window = 5),
                         y = UACCal[, 1] * 2), linewidth = 1) +
          geom_errorbarh(
            aes(y = UACCal[, 1] * 2,
                xmin= UACCal[, 3] - UACCal[, 4], 
                xmax= UACCal[, 3] + UACCal[, 4]), 
            size = 1, alpha = 0.2
          ) +
          scale_x_continuous(breaks = seq(from = -10, to = 10, by = 1),
                             limits = range(UACCal[, 3] - UACCal[, 4],
                                            UACCal[, 3] + UACCal[, 4]),
                             position = 'top') +
          scale_y_reverse(limits = c(370, 0), 
                          breaks = seq(from = 370, to = 0, by = -10),
                          labels = c(rep('', 2),
                                     350,rep('', 4),
                                     300, rep('', 4),
                                     250,rep('', 4),
                                     200, rep('', 4),
                                     150,rep('', 4),
                                     100, rep('', 4),
                                     50,rep('', 4),
                                     0)
          )  +
          guides(x = "axis_truncated",y = "axis_truncated") +
          labs(y = 'Depth (cm)', x = 'Peak - Min', title = '') + 
          theme_classic() + 
          theme(plot.title = element_blank(), 
                axis.text.y = element_text(size = 10, family = 'arial'),
                axis.title.y = element_text(size = 12, family = 'arial'),
                axis.text.x = element_text(size = 10, family = 'arial'),
                axis.title.x = element_text(size = 12, family = 'arial')
          )
        
        ratioUACOHPlotCalSca <- ggplot() +
          geom_point( 
            aes(x = UACCal[, 5],
                y = UACCal[, 1] * 2), 
            size = 2, shape = 21, fill = 'red', alpha = 0.2) + 
          geom_lineh(aes(x = smth.gaussian(UACCal[, 5], window = 5),
                         y = UACCal[, 1] * 2), linewidth = 1) +
          geom_errorbarh(
            aes(y = UACCal[, 1] * 2,
                xmin= UACCal[, 5] - UACCal[, 6], 
                xmax= UACCal[, 5] + UACCal[, 6]), 
            size = 1, alpha = 0.2
          ) +
          scale_x_continuous(breaks = seq(from = -10, to = 10, by = 1),
                             limits = range(UACCal[, 5] - UACCal[, 6],
                                            UACCal[, 5] + UACCal[, 6]),
                             position = 'top') +
          scale_y_reverse(limits = c(370, 0), 
                          breaks = seq(from = 370, to = 0, by = -10),
                          labels = c(rep('', 2),
                                     350,rep('', 4),
                                     300, rep('', 4),
                                     250,rep('', 4),
                                     200, rep('', 4),
                                     150,rep('', 4),
                                     100, rep('', 4),
                                     50,rep('', 4),
                                     0)
          )  +
          guides(x = "axis_truncated", y = "axis_truncated") +
          labs(y = 'Depth (cm)', x = 'Aro-OH Ratio', title = '') + 
          theme_classic() + 
          theme(plot.title = element_blank(), 
                axis.text.y = element_text(size = 10, family = 'arial'),
                axis.title.y = element_text(size = 12, family = 'arial'),
                axis.text.x = element_text(size = 10, family = 'arial'),
                axis.title.x = element_text(size = 12, family = 'arial')
          )
        
        ratio2ndPlotCalSca <- ggplot() +
          geom_point( 
            aes(x = UACCal[, 7],
                y = UACCal[, 1] * 2), 
            size = 2, shape = 21, fill = 'red', alpha = 0.2) + 
          geom_lineh(aes(x = smth.gaussian(UACCal[, 7], window = 5),
                         y = UACCal[, 1] * 2), linewidth = 1) +
          geom_errorbarh(
            aes(y = UACCal[, 1] * 2,
                xmin= UACCal[, 7] - UACCal[, 8], 
                xmax= UACCal[, 7] + UACCal[, 8]), 
            size = 1, alpha = 0.2
          ) +
          scale_x_continuous(breaks = seq(from = -10, to = 10, by = 1),
                             limits = range(UACCal[, 7] - UACCal[, 8],
                                            UACCal[, 7] + UACCal[, 8]),
                             position = 'top') +
          scale_y_reverse(limits = c(370, 0), 
                          breaks = seq(from = 370, to = 0, by = -10),
                          labels = c(rep('', 2),
                                     350,rep('', 4),
                                     300, rep('', 4),
                                     250,rep('', 4),
                                     200, rep('', 4),
                                     150,rep('', 4),
                                     100, rep('', 4),
                                     50,rep('', 4),
                                     0)
          )  +
          guides(x = "axis_truncated",y = "axis_truncated") +
          labs(x = '2nd derivative peak ratio', y = 'Cal years BP', title = '') + 
          theme_classic() + 
          theme(plot.title = element_blank(), 
                axis.text.y = element_text(size = 10, family = 'arial'),
                axis.title.y = element_text(size = 12, family = 'arial'),
                axis.text.x = element_text(size = 10, family = 'arial'),
                axis.title.x = element_text(size = 12, family = 'arial')
          )
        
        local2ndPlotCalSca <- ggplot() +
          geom_point( 
            aes(x = UACCal[, 9],
                y = UACCal[, 1] * 2), 
            size = 2, shape = 21, fill = 'red', alpha = 0.2) + 
          geom_lineh(aes(x = smth.gaussian(UACCal[, 9], window = 5),
                         y = UACCal[, 1] * 2), linewidth = 1) +
          geom_errorbarh(
            aes(y = UACCal[, 1] * 2,
                xmin= UACCal[, 9] - UACCal[, 10], 
                xmax= UACCal[, 9] + UACCal[, 10]), 
            size = 1, alpha = 0.2
          ) +
          scale_x_continuous(breaks = seq(from = -10, to = 10, by = 1),
                             limits = range(UACCal[, 9] - UACCal[, 10],
                                            UACCal[, 9] + UACCal[, 10]),
                             position = 'top') +
          scale_y_reverse(limits = c(370, 0), 
                          breaks = seq(from = 370, to = 0, by = -10),
                          labels = c(rep('', 2),
                                     350,rep('', 4),
                                     300, rep('', 4),
                                     250,rep('', 4),
                                     200, rep('', 4),
                                     150,rep('', 4),
                                     100, rep('', 4),
                                     50,rep('', 4),
                                     0)
          )  +
          guides(x = "axis_truncated",y = "axis_truncated") +
          labs(x = '2nd derivative peak height', y = 'Cal years BP', title = '') + 
          theme_classic() + 
          theme(plot.title = element_blank(), 
                axis.text.y = element_text(size = 10, family = 'arial'),
                axis.title.y = element_text(size = 12, family = 'arial'),
                axis.text.x = element_text(size = 10, family = 'arial'),
                axis.title.x = element_text(size = 12, family = 'arial')
          )
        
        peakAreaPlotCalSca <- ggplot() +
          geom_point( 
            aes(x = UACCal[, 11],
                y = UACCal[, 1] * 2), 
            size = 2, shape = 21, fill = 'red', alpha = 0.2) + 
          geom_lineh(aes(x = smth.gaussian(UACCal[, 11], window = 5),
                         y = UACCal[, 1] * 2), linewidth = 1) +
          geom_errorbarh(
            aes(y = UACCal[, 1] * 2,
                xmin= UACCal[, 11] - UACCal[, 12], 
                xmax= UACCal[, 11] + UACCal[, 12]), 
            size = 1, alpha = 0.2
          ) +
          scale_x_continuous(breaks = seq(from = -10, to = 10, by = 1),
                             limits = range(UACCal[, 11] - UACCal[, 12],
                                            UACCal[, 11] + UACCal[, 12]),
                             position = 'top') +
          scale_y_reverse(limits = c(370, 0), 
                          breaks = seq(from = 370, to = 0, by = -10),
                          labels = c(rep('', 2),
                                     350,rep('', 4),
                                     300, rep('', 4),
                                     250,rep('', 4),
                                     200, rep('', 4),
                                     150,rep('', 4),
                                     100, rep('', 4),
                                     50,rep('', 4),
                                     0)
          )  +
          guides(x = "axis_truncated",y = "axis_truncated") +
          labs(x = 'Peak Area', y = 'Cal years BP', title = '') + 
          theme_classic() + 
          theme(plot.title = element_blank(), 
                axis.text.y = element_text(size = 10, family = 'arial'),
                axis.title.y = element_text(size = 12, family = 'arial'),
                axis.text.x = element_text(size = 10, family = 'arial'),
                axis.title.x = element_text(size = 12, family = 'arial')
          )
        
        meanPlotCalSca <- ggplot() +
          geom_point(
            aes(y = UACCal[, 1] * 2,
                x = UACCalMean), 
            size = 2, shape = 21, fill = 'red', alpha = 0.4) + 
          geom_errorbarh(aes(y = UACCal[, 1] * 2,
                             xmin = UACCalMean - UACCalSD,
                             xmax = UACCalMean + UACCalSD),
                         alpha = 0.4) +
          geom_lineh(aes(x = UACCalGau,
                         y = UACCal[, 1] * 2), size = 1) +
          scale_x_continuous(breaks = seq(from = -10, to = 10, by = 1),
                             limits = range(UACCalMean - UACCalSD,
                                            UACCalMean + UACCalSD),
                             position = 'top') +
          scale_y_reverse(limits = c(370, 0), 
                          breaks = seq(from = 370, to = 0, by = -10),
                          labels = c(rep('', 2),
                                     350,rep('', 4),
                                     300, rep('', 4),
                                     250,rep('', 4),
                                     200, rep('', 4),
                                     150,rep('', 4),
                                     100, rep('', 4),
                                     50,rep('', 4),
                                     0),
                          position = 'right'
          )  +
          guides(x = "axis_truncated", y = "axis_truncated") +
          labs(x = 'Mean of different approaches', y = 'Cal years BP', title = '') + 
          theme_classic() + 
          theme(plot.title = element_blank(), 
                axis.text.y = element_text(size = 10, family = 'arial'),
                axis.title.y = element_text(size = 12, family = 'arial'),
                axis.text.x = element_text(size = 10, family = 'arial'),
                axis.title.x = element_text(size = 12, family = 'arial')
          )
      }
      
      #Combine plot
      {
        ggarrange(ratioUACOHPlotCalSca,
                  localPeakHeightPlotCalSca +
                    theme(
                      axis.line.y = element_blank(),
                      axis.ticks.y = element_blank(),
                      axis.title.y = element_blank(),
                      axis.text.y = element_blank(),
                    ), 
                  peakAreaPlotCalSca +
                    theme(
                      axis.line.y = element_blank(),
                      axis.ticks.y = element_blank(),
                      axis.title.y = element_blank(),
                      axis.text.y = element_blank(),
                    ), 
                  ratio2ndPlotCalSca +
                    
                    theme(
                      axis.line.y = element_blank(),
                      axis.ticks.y = element_blank(),
                      axis.title.y = element_blank(),
                      axis.text.y = element_blank(),
                    ),
                  local2ndPlotCalSca +
                    theme(
                      axis.line.y = element_blank(),
                      axis.ticks.y = element_blank(),
                      axis.title.y = element_blank(),
                      axis.text.y = element_blank(),
                    ),
                  meanPlotCalSca,
                  heights = c(1,1,1,1,1,1),
                  widths = c(1.2,1,1,1,1,1.2),
                  ncol = 6, nrow = 1,
                  align = "h")
      }
    }
  }
  
  #Distribution of UAC values for each approach
  {
    #Sphagnum
    {
      UACSphAroOHDenHistPlot <-
        ggplot() + 
        geom_histogram(aes(x = UACSph$Aro_OH_Ratio), 
                       binwidth=0.1,
                       color = 'black',
                       fill = 'grey')+
        scale_x_continuous(breaks= seq(from = -2, to = 2, by = 0.5)) +
        ylab("Count") + xlab("Aro/OH") + ggtitle("") +
        theme_bw() + theme(plot.title=element_text(size=20),
                           axis.title.y=element_text(size = 16),
                           axis.title.x=element_text(size = 16),
                           axis.text.y=element_text(size = 14),
                           axis.text.x=element_text(size = 14),
                           axis.line = element_line(colour = "black"),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           panel.border = element_blank(),
                           panel.background = element_blank())
      
      UACSphPeakMinDenHistPlot <-
        ggplot() + 
        geom_histogram(aes(x = UACSph$Peak_Min), 
                       binwidth=0.1,
                       color = 'black',
                       fill = 'grey')+
        scale_x_continuous(breaks= seq(from = -2, to = 2, by = 0.5)) +
        ylab("Count") + xlab("Local maximum") + ggtitle("") +
        theme_bw() + theme(plot.title=element_text(size=20),
                           axis.title.y=element_text(size = 16),
                           axis.title.x=element_text(size = 16),
                           axis.text.y=element_text(size = 14),
                           axis.text.x=element_text(size = 14),
                           axis.line = element_line(colour = "black"),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           panel.border = element_blank(),
                           panel.background = element_blank())
      
      UACSphPeakAreaDenHistPlot <-
        ggplot() + 
        geom_histogram(aes(x = UACSph$Peak_Area), 
                       binwidth=0.1,
                       color = 'black',
                       fill = 'grey')+
        scale_x_continuous(breaks= seq(from = -2, to = 2, by = 0.5)) +
        ylab("Count") + xlab("Peak area") + ggtitle("") +
        theme_bw() + theme(plot.title=element_text(size=20),
                           axis.title.y=element_text(size = 16),
                           axis.title.x=element_text(size = 16),
                           axis.text.y=element_text(size = 14),
                           axis.text.x=element_text(size = 14),
                           axis.line = element_line(colour = "black"),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           panel.border = element_blank(),
                           panel.background = element_blank())
      
      UACSph2ndAroOHDenHistPlot <-
        ggplot() + 
        geom_histogram(aes(x = UACSph$De_Ratio), 
                       binwidth=0.1,
                       color = 'black',
                       fill = 'grey')+
        scale_x_continuous(breaks= seq(from = -2, to = 2, by = 0.5)) +
        ylab("Count") + xlab("2nd Derivative Aro/OH") + ggtitle("") +
        theme_bw() + theme(plot.title=element_text(size=20),
                           axis.title.y=element_text(size = 16),
                           axis.title.x=element_text(size = 16),
                           axis.text.y=element_text(size = 14),
                           axis.text.x=element_text(size = 14),
                           axis.line = element_line(colour = "black"),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           panel.border = element_blank(),
                           panel.background = element_blank())
      
      UACSph2ndPeakDenHistPlot <-
        ggplot() + 
        geom_histogram(aes(x = UACSph$De_Aro_Peak), 
                       binwidth=0.1,
                       color = 'black',
                       fill = 'grey')+
        scale_x_continuous(breaks= seq(from = -2, to = 2, by = 0.5)) +
        ylab("Count") + xlab("2nd Derivative Aro peak") + ggtitle("") +
        theme_bw() + theme(plot.title=element_text(size=20),
                           axis.title.y=element_text(size = 16),
                           axis.title.x=element_text(size = 16),
                           axis.text.y=element_text(size = 14),
                           axis.text.x=element_text(size = 14),
                           axis.line = element_line(colour = "black"),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           panel.border = element_blank(),
                           panel.background = element_blank())
      
      UACSphMeanDenHistPlot <-
        ggplot() + 
        geom_histogram(aes(x = UACSphMean), 
                       binwidth=0.1,
                       color = 'black',
                       fill = 'grey')+
        scale_x_continuous(breaks= seq(from = -2, to = 2, by = 0.5)) +
        ylab("Count") + xlab("Mean UACs") + ggtitle("") +
        theme_bw() + theme(plot.title=element_text(size=20),
                           axis.title.y=element_text(size = 16),
                           axis.title.x=element_text(size = 16),
                           axis.text.y=element_text(size = 14),
                           axis.text.x=element_text(size = 14),
                           axis.line = element_line(colour = "black"),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           panel.border = element_blank(),
                           panel.background = element_blank())
      
      #Combine plot
      (UACSphAroOHDenHistPlot | UACSphPeakMinDenHistPlot) /
        (UACSphPeakAreaDenHistPlot | UACSph2ndAroOHDenHistPlot) /
        (UACSph2ndPeakDenHistPlot | UACSphMeanDenHistPlot)
    }
    
    #Alnus
    {
      UACAlnAroOHDenHistPlot <-
        ggplot() + 
        geom_histogram(aes(x = UACAln$Aro_OH_Ratio), 
                       binwidth=0.1,
                       color = 'black',
                       fill = 'grey')+
        scale_x_continuous(breaks= seq(from = -2, to = 2, by = 0.5)) +
        ylab("Count") + xlab("Aro/OH") + ggtitle("") +
        theme_bw() + theme(plot.title=element_text(size=20),
                           axis.title.y=element_text(size = 16),
                           axis.title.x=element_text(size = 16),
                           axis.text.y=element_text(size = 14),
                           axis.text.x=element_text(size = 14),
                           axis.line = element_line(colour = "black"),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           panel.border = element_blank(),
                           panel.background = element_blank())
      
      UACAlnPeakMinDenHistPlot <-
        ggplot() + 
        geom_histogram(aes(x = UACAln$Peak_Min), 
                       binwidth=0.1,
                       color = 'black',
                       fill = 'grey')+
        scale_x_continuous(breaks= seq(from = -2, to = 2, by = 0.5)) +
        ylab("Count") + xlab("Local maximum") + ggtitle("") +
        theme_bw() + theme(plot.title=element_text(size=20),
                           axis.title.y=element_text(size = 16),
                           axis.title.x=element_text(size = 16),
                           axis.text.y=element_text(size = 14),
                           axis.text.x=element_text(size = 14),
                           axis.line = element_line(colour = "black"),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           panel.border = element_blank(),
                           panel.background = element_blank())
      
      UACAlnPeakAreaDenHistPlot <-
        ggplot() + 
        geom_histogram(aes(x = UACAln$Peak_Area), 
                       binwidth=0.1,
                       color = 'black',
                       fill = 'grey')+
        scale_x_continuous(breaks= seq(from = -2, to = 2, by = 0.5)) +
        ylab("Count") + xlab("Peak area") + ggtitle("") +
        theme_bw() + theme(plot.title=element_text(size=20),
                           axis.title.y=element_text(size = 16),
                           axis.title.x=element_text(size = 16),
                           axis.text.y=element_text(size = 14),
                           axis.text.x=element_text(size = 14),
                           axis.line = element_line(colour = "black"),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           panel.border = element_blank(),
                           panel.background = element_blank())
      
      UACAln2ndAroOHDenHistPlot <-
        ggplot() + 
        geom_histogram(aes(x = UACAln$De_Ratio), 
                       binwidth=0.1,
                       color = 'black',
                       fill = 'grey')+
        scale_x_continuous(breaks= seq(from = -2, to = 2, by = 0.5)) +
        ylab("Count") + xlab("2nd Derivative Aro/OH") + ggtitle("") +
        theme_bw() + theme(plot.title=element_text(size=20),
                           axis.title.y=element_text(size = 16),
                           axis.title.x=element_text(size = 16),
                           axis.text.y=element_text(size = 14),
                           axis.text.x=element_text(size = 14),
                           axis.line = element_line(colour = "black"),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           panel.border = element_blank(),
                           panel.background = element_blank())
      
      UACAln2ndPeakDenHistPlot <-
        ggplot() + 
        geom_histogram(aes(x = UACAln$De_Aro_Peak), 
                       binwidth=0.1,
                       color = 'black',
                       fill = 'grey')+
        scale_x_continuous(breaks= seq(from = -2, to = 2, by = 0.5)) +
        ylab("Count") + xlab("2nd Derivative Aro peak") + ggtitle("") +
        theme_bw() + theme(plot.title=element_text(size=20),
                           axis.title.y=element_text(size = 16),
                           axis.title.x=element_text(size = 16),
                           axis.text.y=element_text(size = 14),
                           axis.text.x=element_text(size = 14),
                           axis.line = element_line(colour = "black"),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           panel.border = element_blank(),
                           panel.background = element_blank())
      
      UACAlnMeanDenHistPlot <-
        ggplot() + 
        geom_histogram(aes(x = UACAlnMean), 
                       binwidth=0.1,
                       color = 'black',
                       fill = 'grey')+
        scale_x_continuous(breaks= seq(from = -2, to = 2, by = 0.5)) +
        ylab("Count") + xlab("Mean UACs") + ggtitle("") +
        theme_bw() + theme(plot.title=element_text(size=20),
                           axis.title.y=element_text(size = 16),
                           axis.title.x=element_text(size = 16),
                           axis.text.y=element_text(size = 14),
                           axis.text.x=element_text(size = 14),
                           axis.line = element_line(colour = "black"),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           panel.border = element_blank(),
                           panel.background = element_blank())
      
      #Combine plot
      (UACAlnAroOHDenHistPlot | UACAlnPeakMinDenHistPlot) /
        (UACAlnPeakAreaDenHistPlot | UACAln2ndAroOHDenHistPlot) /
        (UACAln2ndPeakDenHistPlot | UACAlnMeanDenHistPlot)
    }
    
    #Calluna
    {
      UACCalAroOHDenHistPlot <-
        ggplot() + 
        geom_histogram(aes(x = UACCal$Aro_OH_Ratio), 
                       binwidth=0.1,
                       color = 'black',
                       fill = 'grey')+
        scale_x_continuous(breaks= seq(from = -2, to = 2, by = 0.5)) +
        ylab("Count") + xlab("Aro/OH") + ggtitle("") +
        theme_bw() + theme(plot.title=element_text(size=20),
                           axis.title.y=element_text(size = 16),
                           axis.title.x=element_text(size = 16),
                           axis.text.y=element_text(size = 14),
                           axis.text.x=element_text(size = 14),
                           axis.line = element_line(colour = "black"),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           panel.border = element_blank(),
                           panel.background = element_blank())
      
      UACCalPeakMinDenHistPlot <-
        ggplot() + 
        geom_histogram(aes(x = UACCal$Peak_Min), 
                       binwidth=0.1,
                       color = 'black',
                       fill = 'grey')+
        scale_x_continuous(breaks= seq(from = -2, to = 2, by = 0.5)) +
        ylab("Count") + xlab("Local maximum") + ggtitle("") +
        theme_bw() + theme(plot.title=element_text(size=20),
                           axis.title.y=element_text(size = 16),
                           axis.title.x=element_text(size = 16),
                           axis.text.y=element_text(size = 14),
                           axis.text.x=element_text(size = 14),
                           axis.line = element_line(colour = "black"),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           panel.border = element_blank(),
                           panel.background = element_blank())
      
      UACCalPeakAreaDenHistPlot <-
        ggplot() + 
        geom_histogram(aes(x = UACCal$Peak_Area), 
                       binwidth=0.1,
                       color = 'black',
                       fill = 'grey')+
        scale_x_continuous(breaks= seq(from = -2, to = 2, by = 0.5)) +
        ylab("Count") + xlab("Peak area") + ggtitle("") +
        theme_bw() + theme(plot.title=element_text(size=20),
                           axis.title.y=element_text(size = 16),
                           axis.title.x=element_text(size = 16),
                           axis.text.y=element_text(size = 14),
                           axis.text.x=element_text(size = 14),
                           axis.line = element_line(colour = "black"),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           panel.border = element_blank(),
                           panel.background = element_blank())
      
      UACCal2ndAroOHDenHistPlot <-
        ggplot() + 
        geom_histogram(aes(x = UACCal$De_Ratio), 
                       binwidth=0.1,
                       color = 'black',
                       fill = 'grey')+
        scale_x_continuous(breaks= seq(from = -2, to = 2, by = 0.5)) +
        ylab("Count") + xlab("2nd Derivative Aro/OH") + ggtitle("") +
        theme_bw() + theme(plot.title=element_text(size=20),
                           axis.title.y=element_text(size = 16),
                           axis.title.x=element_text(size = 16),
                           axis.text.y=element_text(size = 14),
                           axis.text.x=element_text(size = 14),
                           axis.line = element_line(colour = "black"),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           panel.border = element_blank(),
                           panel.background = element_blank())
      
      UACCal2ndPeakDenHistPlot <-
        ggplot() + 
        geom_histogram(aes(x = UACCal$De_Aro_Peak), 
                       binwidth=0.1,
                       color = 'black',
                       fill = 'grey')+
        scale_x_continuous(breaks= seq(from = -2, to = 2, by = 0.5)) +
        ylab("Count") + xlab("2nd Derivative Aro peak") + ggtitle("") +
        theme_bw() + theme(plot.title=element_text(size=20),
                           axis.title.y=element_text(size = 16),
                           axis.title.x=element_text(size = 16),
                           axis.text.y=element_text(size = 14),
                           axis.text.x=element_text(size = 14),
                           axis.line = element_line(colour = "black"),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           panel.border = element_blank(),
                           panel.background = element_blank())
      
      UACCalMeanDenHistPlot <-
        ggplot() + 
        geom_histogram(aes(x = UACCalMean), 
                       binwidth=0.1,
                       color = 'black',
                       fill = 'grey')+
        scale_x_continuous(breaks= seq(from = -2, to = 2, by = 0.5)) +
        ylab("Count") + xlab("Mean UACs") + ggtitle("") +
        theme_bw() + theme(plot.title=element_text(size=20),
                           axis.title.y=element_text(size = 16),
                           axis.title.x=element_text(size = 16),
                           axis.text.y=element_text(size = 14),
                           axis.text.x=element_text(size = 14),
                           axis.line = element_line(colour = "black"),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           panel.border = element_blank(),
                           panel.background = element_blank())
      
      #Combine plot
      (UACCalAroOHDenHistPlot | UACCalPeakMinDenHistPlot) /
        (UACCalPeakAreaDenHistPlot | UACCal2ndAroOHDenHistPlot) /
        (UACCal2ndPeakDenHistPlot | UACCalMeanDenHistPlot)
    }
  }
  
  #Boxplot of number of grains for each taxa
  {
    #Sphagnum
    UACSphNumPlot <-
      ggplot() +
      geom_boxplot(aes(y = UACSph$n, x = 0)) +
      ylab("Count") + xlab("Sphagnum") + ggtitle("") +
      theme_bw() + theme(plot.title=element_text(size=20),
                         axis.title.y=element_text(size = 16),
                         axis.title.x=element_text(size = 16),
                         axis.text.y=element_text(size = 14),
                         axis.text.x=element_blank(),
                         axis.ticks.x = element_blank(),
                         axis.line = element_line(colour = "black"),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.border = element_blank(),
                         panel.background = element_blank())
    
    #Alnus
    UACAlnNumPlot <-
      ggplot() +
      geom_boxplot(aes(y = UACAln$n, x = 0)) +
      ylab("Count") + xlab("Alnus") + ggtitle("") +
      theme_bw() + theme(plot.title=element_text(size=20),
                         axis.title.y=element_text(size = 16),
                         axis.title.x=element_text(size = 16),
                         axis.text.y=element_text(size = 14),
                         axis.text.x=element_blank(),
                         axis.ticks.x = element_blank(),
                         axis.line = element_line(colour = "black"),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.border = element_blank(),
                         panel.background = element_blank())
    
    #Calluna
    UACCalNumPlot <-
      ggplot() +
      geom_boxplot(aes(y = UACCal$n, x = 0)) +
      ylab("Count") + xlab("Calluna") + ggtitle("") +
      theme_bw() + theme(plot.title=element_text(size=20),
                         axis.title.y=element_text(size = 16),
                         axis.title.x=element_text(size = 16),
                         axis.text.y=element_text(size = 14),
                         axis.text.x=element_blank(),
                         axis.ticks.x = element_blank(),
                         axis.line = element_line(colour = "black"),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.border = element_blank(),
                         panel.background = element_blank())
    
    #Combine
    UACSphNumPlot | UACAlnNumPlot | UACCalNumPlot
    
      
  }
  
  #Correlation matrix between different approaches
  {
    #Sphagnum
    {
      tempMat <- cbind(UACSph[, c(5,3,11,7,9)], UACSphMean)
      colnames(tempMat)[6] <- 'Mean_UAC_Signals'
      
      #Overall correlation matrix
      chart.Correlation(tempMat, histogram=TRUE, pch=19)
      
    }
    
    #Alnus
    {
      tempMat <- cbind(UACAln[, c(5,3,11,7,9)], UACAlnMean)
      colnames(tempMat)[6] <- 'Mean_UAC_Signals'
      
      #Overall correlation matrix
      chart.Correlation(tempMat, histogram=TRUE, pch=19)
      
    }
    
    #Calluna
    {
      tempMat <- cbind(UACCal[, c(5,3,11,7,9)], UACCalMean)
      colnames(tempMat)[6] <- 'Mean_UAC_Signals'
      
      #Overall correlation matrix
      chart.Correlation(tempMat, histogram=TRUE, pch=19)
      
    }
  }
}

#Figure 1: visual comparison
{
  #Plot historical climate events
  {
    histoEventsName <- c( 'RWP', 'DACP', 'LALIA',
                          'MVP/MCA', 'LIA')
    
    histoEventsDown <- c(2300, 1550, 1414, 1150, 700)
    
    histoEventsUp <- c(1600, 1185, 1290, 700, 90)
    
    histoEventsCol <- c('#cecece','#a559aa','#59a89c','#f0c571','#e02b35')
    
    histoEventsHoriMin <- c( 0, 0, 0.25, 0, 0)
    histoEventsHoriMax <- c( 1, 1, 1.25, 1, 1)
    
    histoEventsTextPosX <- c(rep(-0.25, 2), 1.5, rep(-0.25, 2))
    histoEventsTextPosY <- (histoEventsDown + histoEventsUp) / 2
    
    histoEvents <- ggplot() + 
      geom_rect(aes(ymin=histoEventsHoriMin, ymax=histoEventsHoriMax,
                    xmin=histoEventsUp, xmax=histoEventsDown
      ), fill = histoEventsCol, col = 'black') +
      annotate('text', y = histoEventsTextPosX, x = histoEventsTextPosY,
               label = histoEventsName, size = 4, angle = 0) +
      scale_y_continuous(breaks = seq(from = -5, to = 5, by = 1),
                         limits = c(-2, 3),
                         position = 'left') +
      scale_x_reverse(limits = c(2750, -100), 
                      breaks = seq(from = 4000, to = -100, by = -100),
                      labels = c(4000, rep('', 9),
                                 3000, rep('', 9),
                                 2000, rep('', 9),
                                 1000, rep('', 9),
                                 0, '')) +
      labs(x = 'Cal years BP', y = 'Historical climate intervals', title = '') + 
      theme_classic() + theme(plot.title = element_blank(), 
                              axis.title.x = element_blank(),
                              axis.text.y = element_blank(),
                              axis.text.x = element_blank(),
                              axis.ticks.x = element_blank(),
                              axis.line.x = element_blank(),
                              axis.line.y = element_blank(),
                              axis.ticks.y = element_blank(),
                              axis.title.y = element_blank()
      )
  } 
  
  #Read radiocarbon dates
  {
    radioDatesMat <- read.csv(here('Archives_GRL', 'Radiocarbon_Dates.csv'))
    
    #Attach meadian dates to radiocarbon date matrix
    radioMedianDate <- c()
    loopI <- 1
    while (loopI <= nrow(radioDatesMat)) {
      
      radioMedianDate[length(radioMedianDate) + 1] <-
        chronologyHM20[which(as.numeric(chronologyHM20[, 2]) ==
                               as.numeric(radioDatesMat[loopI, 3])), 1]
      
      loopI <- loopI + 1
    }
    
    radioDatesMat <- cbind(radioDatesMat, radioMedianDate)
    
    radioDatesMat[, 3] <- c(1, 1, 2, 1, 2, rep(1, 8))
  }
  
  #Plot radiocarbon dates and date ranges at each date point
  {
    radioDatesPlot <- ggplot() +
      geom_point(data = radioDatesMat, aes(y = depth, x = radioMedianDate),
                 size = 2, shape = 21, fill = 'black') +
      geom_linerange(data = radioDatesMat, aes(y = depth, 
                                               xmin = date_range_low,
                                               xmax = date_range_up),
                     size = 2, alpha = 0.2) +
      scale_y_continuous(breaks = c(0, 1, 2, 3),
                         limits = c(0, 3),
                         position = 'left') +
      scale_x_reverse(limits = c(2750, -100), 
                      breaks = seq(from = 4000, to = -100, by = -100),
                      labels = c(4000, rep('', 4),
                                 3500,rep('', 4),
                                 3000, rep('', 4),
                                 2500,rep('', 4),
                                 2000, rep('', 4),
                                 1500,rep('', 4),
                                 1000, rep('', 4),
                                 500,rep('', 4),
                                 0, ''),
                      position = 'top',
                      sec.axis=sec_axis(~., 
                                        breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                        labels=c(2000, rep('', 4),
                                                 1500,rep('', 4),
                                                 1000, rep('', 4),
                                                 500,rep('', 4),
                                                 0, rep('', 4),
                                                 500,rep('', 4),
                                                 1000, rep('', 4),
                                                 1500,rep('', 4),
                                                 2000)
                      )
      ) +
      labs(y = 'Depth (cm)', x = 'Cal years BP', title = '') + 
      guides(x = "axis_truncated") +
      theme_classic() + theme(plot.title = element_blank(), 
                              axis.title.y = element_blank(),
                              axis.text.x = element_text(size = 10, family = 'arial'),
                              axis.text.y = element_blank(),
                              # axis.ticks.y = element_blank(),
                              # axis.line.y = element_blank(),
                              axis.title.x = element_text(size = 12, family = 'arial')
      )
  }
  
  #Plot Average UAC records, standalone data, mean and SD dash lines
  {
    averageUACStandaloneSph <- ggplot() +
      geom_point(
        aes(x = UACSph[, 2] - UACAreaMatMeanSphCorr,
            y = UACSphMean), 
        size = 2, shape = 21, fill = 'red', alpha = 0.4) + 
      geom_errorbar(aes(x = UACSph[, 2] - UACAreaMatMeanSphCorr,
                        ymin = UACSphMean - UACSphSD,
                        ymax = UACSphMean + UACSphSD),
                    size = 1, alpha = 0.4, color = 'grey') +
      geom_line(aes(y = UACSphGau,
                    x = UACSphGauTime), 
                color = 'black',
                alpha = 1,
                linewidth = 1.5) +
      scale_y_continuous(breaks = seq(from = -10, to = 10, by = 0.5),
                         limits = range(UACSphMean - UACSphSD,
                                        UACSphMean + UACSphSD),
                         position = 'left') +
      scale_x_reverse(limits = c(2750, -100), 
                      breaks = seq(from = 4000, to = -100, by = -100),
                      labels = c(4000, rep('', 4),
                                 3500,rep('', 4),
                                 3000, rep('', 4),
                                 2500,rep('', 4),
                                 2000, rep('', 4),
                                 1500,rep('', 4),
                                 1000, rep('', 4),
                                 500,rep('', 4),
                                 0, ''),
                      sec.axis=sec_axis(~., 
                                        breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                        labels=c(2000, rep('', 4),
                                                 1500,rep('', 4),
                                                 1000, rep('', 4),
                                                 500,rep('', 4),
                                                 0, rep('', 4),
                                                 500,rep('', 4),
                                                 1000, rep('', 4),
                                                 1500,rep('', 4),
                                                 2000)
                      )) +
      guides(y = "axis_truncated") +
      labs(y = 'Sphagnum', x = 'Cal years BP', title = '') + 
      theme_classic() + 
      theme(plot.title = element_blank(), 
            axis.text.y = element_text(size = 10, family = 'arial'),
            axis.title.y = element_text(size = 12, family = 'arial'),
            axis.text.x = element_text(size = 10, family = 'arial'),
            axis.title.x = element_text(size = 12, family = 'arial')
      )
    
    averageUACStandaloneAln <- ggplot() +
      geom_point(
        aes(x = UACAln[, 2] - UACAreaMatMeanAlnCorr,
            y = UACAlnMean), 
        size = 2, shape = 21, fill = 'red', alpha = 0.4) + 
      geom_errorbar(aes(x = UACAln[, 2] - UACAreaMatMeanAlnCorr,
                        ymin = UACAlnMean - UACAlnSD,
                        ymax = UACAlnMean + UACAlnSD),
                    size = 1, alpha = 0.4, color = 'grey') +
      geom_line(aes(y = UACAlnGau,
                    x = UACAlnGauTime), 
                color = 'black',
                alpha = 1,
                linewidth = 1.5) +
      scale_y_continuous(breaks = seq(from = -10, to = 10, by = 0.5),
                         limits = range(UACAlnMean - UACAlnSD,
                                        UACAlnMean + UACAlnSD),
                         position = 'left') +
      scale_x_reverse(limits = c(2750, -100), 
                      breaks = seq(from = 4000, to = -100, by = -100),
                      labels = c(4000, rep('', 4),
                                 3500,rep('', 4),
                                 3000, rep('', 4),
                                 2500,rep('', 4),
                                 2000, rep('', 4),
                                 1500,rep('', 4),
                                 1000, rep('', 4),
                                 500,rep('', 4),
                                 0, ''),
                      sec.axis=sec_axis(~., 
                                        breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                        labels=c(2000, rep('', 4),
                                                 1500,rep('', 4),
                                                 1000, rep('', 4),
                                                 500,rep('', 4),
                                                 0, rep('', 4),
                                                 500,rep('', 4),
                                                 1000, rep('', 4),
                                                 1500,rep('', 4),
                                                 2000)
                      )) +
      guides(y = "axis_truncated") +
      labs(y = 'Alnagnum', x = 'Cal years BP', title = '') + 
      theme_classic() + 
      theme(plot.title = element_blank(), 
            axis.text.y = element_text(size = 10, family = 'arial'),
            axis.title.y = element_text(size = 12, family = 'arial'),
            axis.text.x = element_text(size = 10, family = 'arial'),
            axis.title.x = element_text(size = 12, family = 'arial')
      )
    
    averageUACStandaloneCal <- ggplot() +
      geom_point(
        aes(x = UACCal[, 2] - UACAreaMatMeanCalCorr,
            y = UACCalMean), 
        size = 2, shape = 21, fill = 'red', alpha = 0.4) + 
      geom_errorbar(aes(x = UACCal[, 2] - UACAreaMatMeanCalCorr,
                        ymin = UACCalMean - UACCalSD,
                        ymax = UACCalMean + UACCalSD),
                    size = 1, alpha = 0.4, color = 'grey') +
      geom_line(aes(y = UACCalGau,
                    x = UACCalGauTime), 
                color = 'black',
                alpha = 1,
                linewidth = 1.5) +
      scale_y_continuous(breaks = seq(from = -10, to = 10, by = 0.5),
                         limits = c(-1.5,1.5),
                         position = 'left') +
      scale_x_reverse(limits = c(2750, -100), 
                      breaks = seq(from = 4000, to = -100, by = -100),
                      labels = c(4000, rep('', 4),
                                 3500,rep('', 4),
                                 3000, rep('', 4),
                                 2500,rep('', 4),
                                 2000, rep('', 4),
                                 1500,rep('', 4),
                                 1000, rep('', 4),
                                 500,rep('', 4),
                                 0, ''),
                      sec.axis=sec_axis(~., 
                                        breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                        labels=c(2000, rep('', 4),
                                                 1500,rep('', 4),
                                                 1000, rep('', 4),
                                                 500,rep('', 4),
                                                 0, rep('', 4),
                                                 500,rep('', 4),
                                                 1000, rep('', 4),
                                                 1500,rep('', 4),
                                                 2000)
                      )) +
      guides(y = "axis_truncated") +
      labs(y = 'Calagnum', x = 'Cal years BP', title = '') + 
      theme_classic() + 
      theme(plot.title = element_blank(), 
            axis.text.y = element_text(size = 10, family = 'arial'),
            axis.title.y = element_text(size = 12, family = 'arial'),
            axis.text.x = element_text(size = 10, family = 'arial'),
            axis.title.x = element_text(size = 12, family = 'arial')
      )
  }
  
  #Plot Figure
  {
    #Combined plot
    {
      ggarrange(radioDatesPlot +
                  scale_x_reverse(limits = c(2750, -100), 
                                  breaks = seq(from = 4000, to = 0, by = -100),
                                  labels = c(4000, rep('', 4),
                                             3500,rep('', 4),
                                             3000, rep('', 4),
                                             2500,rep('', 4),
                                             2000, rep('', 4),
                                             1500,rep('', 4),
                                             1000, rep('', 4),
                                             500,rep('', 4),
                                             0),
                                  position = 'top'), 
                histoEvents +
                  scale_x_reverse(limits = c(2750, -100), 
                                  breaks = seq(from = 4000, to = 0, by = -100),
                                  labels = c(4000, rep('', 4),
                                             3500,rep('', 4),
                                             3000, rep('', 4),
                                             2500,rep('', 4),
                                             2000, rep('', 4),
                                             1500,rep('', 4),
                                             1000, rep('', 4),
                                             500,rep('', 4),
                                             0)) +
                  theme(
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                  ),
                averageUACStandaloneSph,
                averageUACStandaloneAln,
                averageUACStandaloneCal,
                heights = c(0.3,0.5,1,1,1),
                widths = c(1,1,1,1,1),
                ncol = 1, nrow = 5,
                align = "v")
    }
  }
}

#Figure 2: comparison with TSI
{
  #Read TSI records (Now using newer TSI data from Wu et al., 2018)
  {
    #Read TSI data
    TSIRaw <- apply(as.matrix(read.csv(here('SATIRE-M_wu18_tsi.csv'))), c(1,2),
                    as.numeric)
    
    #Transfer tbe BC/AD ages into Cal. years BP
    TSIRaw[, 1] <- 1950 - TSIRaw[ ,1]
    
    #Rename the columns
    colnames(TSIRaw) <- c('Years_BP', 'dTSI')
    
    #Transfer the matrix into data frame
    TSIRaw <- as.data.frame(TSIRaw)
    
    #Intepolate the dTSI data according to time sequence determined by UAC signals
    # TSIInterPola <- approx(x = TSIRaw[, 1], y = TSIRaw[, 2],
    #                        xout = interTimeSeq)[[2]]
    
    #Smooth the dTSI data at level of 80 years
    TSISmooPlot <- smth.gaussian(TSIRaw[,2], window = 9)
  }
  
  #Plot TSI data
  {
    #TSI data
    WholeTSI <- ggplot() +
      geom_point(data = as.data.frame(TSIRaw), aes(x = Years_BP, y = dTSI), 
                 size = 2, shape = 21, fill = '#00b0be', alpha = 0.4) + 
      geom_line(aes(x = TSIRaw[,1], y = TSISmooPlot), color = '#00b0be', size = 1) +
      scale_y_continuous(breaks = seq(from = min(floor(TSIRaw$dTSI)), 
                                      to = max(ceiling(TSIRaw$dTSI)), by = 0.4),
                         limits = range(TSIRaw[, 2]),
                         position = 'right') +
      scale_x_reverse(limits = c(2750, -100), 
                      breaks = seq(from = 4000, to = -100, by = -100),
                      labels = c(4000, rep('', 4),
                                 3500,rep('', 4),
                                 3000, rep('', 4),
                                 2500,rep('', 4),
                                 2000, rep('', 4),
                                 1500,rep('', 4),
                                 1000, rep('', 4),
                                 500,rep('', 4),
                                 0, ''))  +
      labs(y = 'TSI (Wm-2)', x = 'Cal years BP', title = '') + 
      guides(y = "axis_truncated", x ="axis_truncated") +
      theme_classic() + theme(plot.title = element_blank(), 
                              axis.title.y = element_text(size = 12, family = 'arial'),
                              axis.text.x = element_text(size = 10, family = 'arial'),
                              axis.text.y = element_text(size = 10, family = 'arial'),
                              axis.title.x = element_text(size = 12, family = 'arial')
      )
  }
  
  #Moving average
  {
    #Setup environment for parallel processing
    c16 <- makeCluster(16, type = 'SOCK')
    registerDoSNOW(c16)
    
    #Sphagnum
    {
      #Setup moving correlation window size and window moving step
      movCorreWindowSize <- 700
      movCorreWindowStep <- 40
      
      #Calculate the time sequence for moving correlation
      movCorreTimeSeq <- seq(from = min(UACSph[, 2]) 
                             + movCorreWindowSize / 2, to = max(UACSph[, 2]) 
                             - movCorreWindowSize / 2, by = movCorreWindowStep)
      
      movCorreTimeSeq <- movCorreTimeSeq[which(movCorreTimeSeq<= 2400)]
      
      chronoCorrMeanMatSca <- c()
      #Loop to fill the matrix
      loopI <- 1
      while (loopI <= nrow(UACSph)) {
        
        chronoCorrMeanMatSca[length(chronoCorrMeanMatSca) + 1] <- chronologyHM20TransferInte[which(chronologyHM20TransferInte[, 1] == 
                                                                                                     (as.numeric(UACSph[loopI, 1]) * 2 )), 2]
        
        loopI <- loopI + 1
      }
      
      rm(loopI)
      
      #create a matrix for storing correlation coef. and significance
      movCorreMatSphTSI <- matrix(nrow = length(movCorreTimeSeq),
                                  ncol = 4,
                                  dimnames = list(c(),
                                                  c('Time', 'Coef', 'Sig1', 'Sig2')))
      
      
      movCorreMatSphTSI[, 1] <- movCorreTimeSeq
      
      movCorreMatSphTSI[, 2:4] <-
        try(foreach(paraInd = c(1:length(movCorreTimeSeq)), .packages = c("zoo"),
                    .combine = 'rbind') %dopar% 
              {
                try(
                  {
                    tempTime1 <- movCorreTimeSeq[paraInd] - movCorreWindowSize / 2
                    tempTime2 <- movCorreTimeSeq[paraInd] + movCorreWindowSize / 2
                    tempIndUAC <- which((UACSph[,2] >= tempTime1) & 
                                          (UACSph[,2] <= tempTime2))
                    
                    tempIndTSI <- which((TSIRaw[,1] >= tempTime1) & 
                                          (TSIRaw[,1] <= tempTime2))
                    
                    #Correlation between Sphagnum UAC and TSI using linear interpolation and Gaussian filtering
                    time.series1 <- zoo::zoo(UACSphMean[tempIndUAC], 
                                             order.by = UACSph[tempIndUAC, 2] - 
                                               chronoCorrMeanMatSca[tempIndUAC])	#create a zoo-object
                    time.series2 <- zoo::zoo(TSIRaw[tempIndTSI, 2], order.by = TSIRaw[tempIndTSI, 1])
                    Cor <- corit::CorIrregTimser(
                      timser1 = time.series1,
                      timser2 = time.series2,
                      detr = FALSE,	#remove linear trend time series
                      method = "InterpolationMethod",
                      appliedFilter = "gauss",
                      fc = 1/80,	#cut-off frequency
                      dt = 20,	#time step for the interpolation
                      int.method = "linear",	#kind of interpolation
                      filt.output = FALSE)	#return filtered time series 
                    
                    #(2) applying a significance test for the correlation estimate based on the correlation of independent noise
                    slopes <- corit::estimateTimserSlopes(	#estimate spectral slopes of the time series
                      timeseries1 = time.series1,
                      timeseries2 = time.series2,
                      int.step = 1)	#time step of the interpolated time series
                    Quant <- corit::CorQuantilesNullHyp(	#quantiles estimated based on surrogate correlations
                      timser1 = time.series1,
                      timser2 = time.series2,
                      beta.noise1 = slopes$s1,
                      beta.noise2 = slopes$s2,
                      detr = FALSE,
                      rep = 1000,	#repetition during Monte Carlo procedure
                      quant = c(0.05, 0.95),	#quantiles to be estimated
                      method = "InterpolationMethod",
                      appliedFilter = "gauss",
                      fc = 1/80,
                      dt = 20,
                      int.method = "linear")
                    
                    return(c(Cor, Quant[[2]][[1]][1], Quant[[2]][[1]][2]))
                  },
                  silent = TRUE
                )
              },
            silent = TRUE
        )
      movCorreMatSphTSI <- apply(movCorreMatSphTSI, c(1,2), as.numeric)
    }
    
    #Alnus
    {
      #Setup moving correlation window size and window moving step
      movCorreWindowSize <- 700
      movCorreWindowStep <- 40
      
      #Calculate the time sequence for moving correlation
      movCorreTimeSeq <- seq(from = min(UACAln[, 2]) 
                             + movCorreWindowSize / 2, to = max(UACAln[, 2]) 
                             - movCorreWindowSize / 2, by = movCorreWindowStep)
      
      movCorreTimeSeq <- movCorreTimeSeq[which(movCorreTimeSeq<= 2400)]
      
      chronoCorrMeanMatSca <- c()
      #Loop to fill the matrix
      loopI <- 1
      while (loopI <= nrow(UACAln)) {
        
        chronoCorrMeanMatSca[length(chronoCorrMeanMatSca) + 1] <- chronologyHM20TransferInte[which(chronologyHM20TransferInte[, 1] == 
                                                                                                     (as.numeric(UACAln[loopI, 1]) * 2 )), 2]
        
        loopI <- loopI + 1
      }
      
      rm(loopI)
      
      #create a matrix for storing correlation coef. and significance
      movCorreMatAlnTSI <- matrix(nrow = length(movCorreTimeSeq),
                                  ncol = 4,
                                  dimnames = list(c(),
                                                  c('Time', 'Coef', 'Sig1', 'Sig2')))
      
      
      movCorreMatAlnTSI[, 1] <- movCorreTimeSeq
      
      movCorreMatAlnTSI[, 2:4] <-
        try(foreach(paraInd = c(1:length(movCorreTimeSeq)), .packages = c("zoo"),
                    .combine = 'rbind') %dopar% 
              {
                try(
                  {
                    tempTime1 <- movCorreTimeSeq[paraInd] - movCorreWindowSize / 2
                    tempTime2 <- movCorreTimeSeq[paraInd] + movCorreWindowSize / 2
                    tempIndUAC <- which((UACAln[,2] >= tempTime1) & 
                                          (UACAln[,2] <= tempTime2))
                    
                    tempIndTSI <- which((TSIRaw[,1] >= tempTime1) & 
                                          (TSIRaw[,1] <= tempTime2))
                    
                    #Correlation between Alnagnum UAC and TSI using linear interpolation and Gaussian filtering
                    time.series1 <- zoo::zoo(UACAlnMean[tempIndUAC], 
                                             order.by = UACAln[tempIndUAC, 2] - 
                                               chronoCorrMeanMatSca[tempIndUAC])	#create a zoo-object
                    time.series2 <- zoo::zoo(TSIRaw[tempIndTSI, 2], order.by = TSIRaw[tempIndTSI, 1])
                    Cor <- corit::CorIrregTimser(
                      timser1 = time.series1,
                      timser2 = time.series2,
                      detr = FALSE,	#remove linear trend time series
                      method = "InterpolationMethod",
                      appliedFilter = "gauss",
                      fc = 1/80,	#cut-off frequency
                      dt = 20,	#time step for the interpolation
                      int.method = "linear",	#kind of interpolation
                      filt.output = FALSE)	#return filtered time series 
                    
                    #(2) applying a significance test for the correlation estimate based on the correlation of independent noise
                    slopes <- corit::estimateTimserSlopes(	#estimate spectral slopes of the time series
                      timeseries1 = time.series1,
                      timeseries2 = time.series2,
                      int.step = 1)	#time step of the interpolated time series
                    Quant <- corit::CorQuantilesNullHyp(	#quantiles estimated based on surrogate correlations
                      timser1 = time.series1,
                      timser2 = time.series2,
                      beta.noise1 = slopes$s1,
                      beta.noise2 = slopes$s2,
                      detr = FALSE,
                      rep = 1000,	#repetition during Monte Carlo procedure
                      quant = c(0.05, 0.95),	#quantiles to be estimated
                      method = "InterpolationMethod",
                      appliedFilter = "gauss",
                      fc = 1/80,
                      dt = 20,
                      int.method = "linear")
                    
                    return(c(Cor, Quant[[2]][[1]][1], Quant[[2]][[1]][2]))
                  },
                  silent = TRUE
                )
              },
            silent = TRUE
        )
      
      movCorreMatAlnTSI <- apply(movCorreMatAlnTSI, c(1,2), as.numeric)
    }
    
    #Calluna
    {
      #Setup moving correlation window size and window moving step
      movCorreWindowSize <- 700
      movCorreWindowStep <- 40
      
      #Calculate the time sequence for moving correlation
      movCorreTimeSeq <- seq(from = min(UACCal[, 2]) 
                             + movCorreWindowSize / 2, to = max(UACCal[, 2]) 
                             - movCorreWindowSize / 2, by = movCorreWindowStep)
      
      movCorreTimeSeq <- movCorreTimeSeq[which(movCorreTimeSeq<= 2400)]
      
      chronoCorrMeanMatSca <- c()
      #Loop to fill the matrix
      loopI <- 1
      while (loopI <= nrow(UACCal)) {
        
        chronoCorrMeanMatSca[length(chronoCorrMeanMatSca) + 1] <- chronologyHM20TransferInte[which(chronologyHM20TransferInte[, 1] == 
                                                                                                     (as.numeric(UACCal[loopI, 1]) * 2 )), 2]
        
        loopI <- loopI + 1
      }
      
      rm(loopI)
      
      #create a matrix for storing correlation coef. and significance
      movCorreMatCalTSI <- matrix(nrow = length(movCorreTimeSeq),
                                  ncol = 4,
                                  dimnames = list(c(),
                                                  c('Time', 'Coef', 'Sig1', 'Sig2')))
      
      
      movCorreMatCalTSI[, 1] <- movCorreTimeSeq
      
      movCorreMatCalTSI[, 2:4] <-
        try(foreach(paraInd = c(1:length(movCorreTimeSeq)), .packages = c("zoo"),
                    .combine = 'rbind') %dopar% 
              {
                try(
                  {
                    tempTime1 <- movCorreTimeSeq[paraInd] - movCorreWindowSize / 2
                    tempTime2 <- movCorreTimeSeq[paraInd] + movCorreWindowSize / 2
                    tempIndUAC <- which((UACCal[,2] >= tempTime1) & 
                                          (UACCal[,2] <= tempTime2))
                    
                    tempIndTSI <- which((TSIRaw[,1] >= tempTime1) & 
                                          (TSIRaw[,1] <= tempTime2))
                    
                    #Correlation between Calagnum UAC and TSI using linear interpolation and Gaussian filtering
                    time.series1 <- zoo::zoo(UACCalMean[tempIndUAC], 
                                             order.by = UACCal[tempIndUAC, 2] - 
                                               chronoCorrMeanMatSca[tempIndUAC])	#create a zoo-object
                    time.series2 <- zoo::zoo(TSIRaw[tempIndTSI, 2], order.by = TSIRaw[tempIndTSI, 1])
                    Cor <- corit::CorIrregTimser(
                      timser1 = time.series1,
                      timser2 = time.series2,
                      detr = FALSE,	#remove linear trend time series
                      method = "InterpolationMethod",
                      appliedFilter = "gauss",
                      fc = 1/80,	#cut-off frequency
                      dt = 20,	#time step for the interpolation
                      int.method = "linear",	#kind of interpolation
                      filt.output = FALSE)	#return filtered time series 
                    
                    #(2) applying a significance test for the correlation estimate based on the correlation of independent noise
                    slopes <- corit::estimateTimserSlopes(	#estimate spectral slopes of the time series
                      timeseries1 = time.series1,
                      timeseries2 = time.series2,
                      int.step = 1)	#time step of the interpolated time series
                    Quant <- corit::CorQuantilesNullHyp(	#quantiles estimated based on surrogate correlations
                      timser1 = time.series1,
                      timser2 = time.series2,
                      beta.noise1 = slopes$s1,
                      beta.noise2 = slopes$s2,
                      detr = FALSE,
                      rep = 1000,	#repetition during Monte Carlo procedure
                      quant = c(0.05, 0.95),	#quantiles to be estimated
                      method = "InterpolationMethod",
                      appliedFilter = "gauss",
                      fc = 1/80,
                      dt = 20,
                      int.method = "linear")
                    
                    return(c(Cor, Quant[[2]][[1]][1], Quant[[2]][[1]][2]))
                  },
                  silent = TRUE
                )
              },
            silent = TRUE
        )
      
      movCorreMatCalTSI <- apply(movCorreMatCalTSI, c(1,2), as.numeric)
    }
    
    stopCluster(c16)
  }
  
  #Plot Moving average
  {
    #Sph TSI
    movSphTSIPlot <- ggplot() +
      geom_line(aes(y = movCorreMatSphTSI[, 2],
                    x = movCorreMatSphTSI[, 1]), color = 'black') +
      geom_line(aes(y = movCorreMatSphTSI[, 3],
                    x = movCorreMatSphTSI[, 1]), color = 'red', linetype = 2, alpha = 0.6) +
      geom_line(aes(y = movCorreMatSphTSI[, 4],
                    x = movCorreMatSphTSI[, 1]), color = 'red', linetype = 2, alpha = 0.6) +
      geom_hline(yintercept = 0, color = 'grey', alpha = 0.6) +
      scale_y_continuous(breaks = seq(from = -10, to = 10, by = 0.5),
                         limits = c(-1,1),
                         position = 'left') +
      scale_x_reverse(limits = c(2750, -100), 
                      breaks = seq(from = 4000, to = -100, by = -100),
                      labels = c(4000, rep('', 4),
                                 3500,rep('', 4),
                                 3000, rep('', 4),
                                 2500,rep('', 4),
                                 2000, rep('', 4),
                                 1500,rep('', 4),
                                 1000, rep('', 4),
                                 500,rep('', 4),
                                 0, ''),
                      sec.axis=sec_axis(~., 
                                        breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                        labels=c(2000, rep('', 4),
                                                 1500,rep('', 4),
                                                 1000, rep('', 4),
                                                 500,rep('', 4),
                                                 0, rep('', 4),
                                                 500,rep('', 4),
                                                 1000, rep('', 4),
                                                 1500,rep('', 4),
                                                 2000)
                      )
      )  +
      guides(y = "axis_truncated", x = "axis_truncated") +
      labs(y = 'Sph&TSI', x = 'Cal years BP', title = '') + 
      theme_classic() + 
      theme(plot.title = element_blank(), 
            axis.text.x = element_text(size = 10, family = 'arial'),
            axis.title.x = element_text(size = 12, family = 'arial'),
            axis.text.y = element_text(size = 10, family = 'arial'),
            axis.title.y = element_text(size = 12, family = 'arial')
      )
    
    #Aln TSI
    movAlnTSIPlot <- ggplot() +
      geom_line(aes(y = movCorreMatAlnTSI[, 2],
                    x = movCorreMatAlnTSI[, 1]), color = 'black') +
      geom_line(aes(y = movCorreMatAlnTSI[, 3],
                    x = movCorreMatAlnTSI[, 1]), color = 'red', linetype = 2, alpha = 0.6) +
      geom_line(aes(y = movCorreMatAlnTSI[, 4],
                    x = movCorreMatAlnTSI[, 1]), color = 'red', linetype = 2, alpha = 0.6) +
      geom_hline(yintercept = 0, color = 'grey', alpha = 0.6) +
      scale_y_continuous(breaks = seq(from = -10, to = 10, by = 0.5),
                         limits = c(-1,1),
                         position = 'left') +
      scale_x_reverse(limits = c(2750, -100), 
                      breaks = seq(from = 4000, to = -100, by = -100),
                      labels = c(4000, rep('', 4),
                                 3500,rep('', 4),
                                 3000, rep('', 4),
                                 2500,rep('', 4),
                                 2000, rep('', 4),
                                 1500,rep('', 4),
                                 1000, rep('', 4),
                                 500,rep('', 4),
                                 0, ''),
                      sec.axis=sec_axis(~., 
                                        breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                        labels=c(2000, rep('', 4),
                                                 1500,rep('', 4),
                                                 1000, rep('', 4),
                                                 500,rep('', 4),
                                                 0, rep('', 4),
                                                 500,rep('', 4),
                                                 1000, rep('', 4),
                                                 1500,rep('', 4),
                                                 2000)
                      )
      )  +
      guides(y = "axis_truncated", x = "axis_truncated") +
      labs(y = 'Aln&TSI', x = 'Cal years BP', title = '') + 
      theme_classic() + 
      theme(plot.title = element_blank(), 
            axis.text.x = element_text(size = 10, family = 'arial'),
            axis.title.x = element_text(size = 12, family = 'arial'),
            axis.text.y = element_text(size = 10, family = 'arial'),
            axis.title.y = element_text(size = 12, family = 'arial')
      )
    
    #Cal TSI
    movCalTSIPlot <- ggplot() +
      geom_line(aes(y = movCorreMatCalTSI[, 2],
                    x = movCorreMatCalTSI[, 1]), color = 'black') +
      geom_line(aes(y = movCorreMatCalTSI[, 3],
                    x = movCorreMatCalTSI[, 1]), color = 'red', linetype = 2, alpha = 0.6) +
      geom_line(aes(y = movCorreMatCalTSI[, 4],
                    x = movCorreMatCalTSI[, 1]), color = 'red', linetype = 2, alpha = 0.6) +
      geom_hline(yintercept = 0, color = 'grey', alpha = 0.6) +
      scale_y_continuous(breaks = seq(from = -10, to = 10, by = 0.5),
                         limits = c(-1,1),
                         position = 'left') +
      scale_x_reverse(limits = c(2750, -100), 
                      breaks = seq(from = 4000, to = -100, by = -100),
                      labels = c(4000, rep('', 4),
                                 3500,rep('', 4),
                                 3000, rep('', 4),
                                 2500,rep('', 4),
                                 2000, rep('', 4),
                                 1500,rep('', 4),
                                 1000, rep('', 4),
                                 500,rep('', 4),
                                 0, ''),
                      sec.axis=sec_axis(~., 
                                        breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                        labels=c(2000, rep('', 4),
                                                 1500,rep('', 4),
                                                 1000, rep('', 4),
                                                 500,rep('', 4),
                                                 0, rep('', 4),
                                                 500,rep('', 4),
                                                 1000, rep('', 4),
                                                 1500,rep('', 4),
                                                 2000)
                      )
      )  +
      guides(y = "axis_truncated", x = "axis_truncated") +
      labs(y = 'Cal&TSI', x = 'Cal years BP', title = '') + 
      theme_classic() + 
      theme(plot.title = element_blank(), 
            axis.text.x = element_text(size = 10, family = 'arial'),
            axis.title.x = element_text(size = 12, family = 'arial'),
            axis.text.y = element_text(size = 10, family = 'arial'),
            axis.title.y = element_text(size = 12, family = 'arial')
      )
  }
  
  #Plot Figure
  {
    #Combined plot
    {
      ggarrange(radioDatesPlot, 
                averageUACStandaloneSph +
                  geom_hline(yintercept = UACSphGauMean, 
                             linetype = 1, 
                             color = 'blue', alpha = 0.5) +
                  geom_hline(yintercept = UACSphGauMean - UACSphGauSD, 
                             linetype = 2, 
                             color = 'blue', alpha = 0.5) +
                  geom_hline(yintercept = UACSphGauMean - 2 * UACSphGauSD, 
                             linetype = 2, 
                             color = 'blue', alpha = 0.5) +
                  geom_hline(yintercept = UACSphGauMean + UACSphGauSD, 
                             linetype = 2, 
                             color = 'blue', alpha = 0.5) +
                  geom_hline(yintercept = UACSphGauMean + 2* UACSphGauSD, 
                             linetype = 2, 
                             color = 'blue', alpha = 0.5) +
                  theme(
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                  ),
                WholeTSI +
                  theme(
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                  ),
                movSphTSIPlot,
                heights = c(0.3,1,1,1),
                widths = c(1,1,1,1),
                ncol = 1, nrow = 4,
                align = "v")
      
      ggarrange(radioDatesPlot, 
                averageUACStandaloneAln +
                  geom_hline(yintercept = UACAlnGauMean, 
                             linetype = 1, 
                             color = 'blue', alpha = 0.5) +
                  geom_hline(yintercept = UACAlnGauMean - UACAlnGauSD, 
                             linetype = 2, 
                             color = 'blue', alpha = 0.5) +
                  geom_hline(yintercept = UACAlnGauMean - 2 * UACAlnGauSD, 
                             linetype = 2, 
                             color = 'blue', alpha = 0.5) +
                  geom_hline(yintercept = UACAlnGauMean + UACAlnGauSD, 
                             linetype = 2, 
                             color = 'blue', alpha = 0.5) +
                  geom_hline(yintercept = UACAlnGauMean + 2* UACAlnGauSD, 
                             linetype = 2, 
                             color = 'blue', alpha = 0.5) +
                  theme(
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                  ),
                WholeTSI +
                  theme(
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                  ),
                movAlnTSIPlot,
                heights = c(0.3,1,1,1),
                widths = c(1,1,1,1),
                ncol = 1, nrow = 4,
                align = "v")
      
      ggarrange(radioDatesPlot, 
                averageUACStandaloneCal +
                  geom_hline(yintercept = UACCalGauMean, 
                             linetype = 1, 
                             color = 'blue', alpha = 0.5) +
                  geom_hline(yintercept = UACCalGauMean - UACCalGauSD, 
                             linetype = 2, 
                             color = 'blue', alpha = 0.5) +
                  geom_hline(yintercept = UACCalGauMean - 2 * UACCalGauSD, 
                             linetype = 2, 
                             color = 'blue', alpha = 0.5) +
                  geom_hline(yintercept = UACCalGauMean + UACCalGauSD, 
                             linetype = 2, 
                             color = 'blue', alpha = 0.5) +
                  geom_hline(yintercept = UACCalGauMean + 2* UACCalGauSD, 
                             linetype = 2, 
                             color = 'blue', alpha = 0.5) +
                  theme(
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                  ),
                WholeTSI +
                  theme(
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                  ),
                movCalTSIPlot,
                heights = c(0.3,1,1,1),
                widths = c(1,1,1,1),
                ncol = 1, nrow = 4,
                align = "v")
      
    }
  }
}

#Figure 3: comparison with differnt proxy
{
  #Read NH temperature anomalies from van Dijk et al., 2024
  {
    tempAnomaRaw <- read.csv(file = here('Archives_GRL', 'NHTemp.csv'))
    
    colnames(tempAnomaRaw) <- c('Date', 'NHTemp')
    
    tempAnomaRaw <- apply(tempAnomaRaw, c(1,2), as.numeric)
    
    tempAnomaRaw <- as.data.frame(tempAnomaRaw)
    
    #Smooth the dTSI data at level of 100 years
    tempAnomaSmooPlot <- smth.gaussian(tempAnomaRaw[,2], window = 81)
    
    tempAnomaSmooPlotTime <- tempAnomaRaw[,1]
  }
  
  #Plot NH temperature
  {
    #tempAnoma data
    WholetempAnoma <- ggplot() +
      geom_line(data = as.data.frame(tempAnomaRaw), aes(x = Date, y = NHTemp), 
                size = 1, color = '#ffcd8e', alpha = 0.4) +
      geom_line(aes(x = tempAnomaSmooPlotTime, y = tempAnomaSmooPlot), 
                size = 1, color = '#ffcd8e', alpha = 1) +
      geom_hline(yintercept = 0, color = 'grey', alpha = 0.6) +
      scale_y_continuous(
        limits = range(tempAnomaRaw[, 2]),
        position = 'left') +
      scale_x_reverse(limits = c(2750, -100), 
                      breaks = seq(from = 4000, to = -100, by = -100),
                      labels = c(4000, rep('', 4),
                                 3500,rep('', 4),
                                 3000, rep('', 4),
                                 2500,rep('', 4),
                                 2000, rep('', 4),
                                 1500,rep('', 4),
                                 1000, rep('', 4),
                                 500,rep('', 4),
                                 0, ''),
                      sec.axis=sec_axis(~., 
                                        breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                        labels=c(2000, rep('', 4),
                                                 1500,rep('', 4),
                                                 1000, rep('', 4),
                                                 500,rep('', 4),
                                                 0, rep('', 4),
                                                 500,rep('', 4),
                                                 1000, rep('', 4),
                                                 1500,rep('', 4),
                                                 2000)
                      )
      )  +
      labs(y = 'NH Temperature', x = 'Cal years BP', title = '') + 
      guides(y = "axis_truncated",x = "axis_truncated") +
      theme_classic() + theme(plot.title = element_blank(), 
                              axis.title.y = element_text(size = 12, family = 'arial'),
                              axis.text.x = element_text(size = 10, family = 'arial'),
                              axis.text.y = element_text(size = 10, family = 'arial'),
                              axis.title.x = element_text(size = 12, family = 'arial')
      )
  }
  
  #Read SSI from van Sicre et al., 2008
  {
    SSIRaw <- read.csv(file = here('SSI.csv'))
    
    colnames(SSIRaw) <- c('Date', 'SSI')
    
    SSIRaw <- apply(SSIRaw, c(1,2), as.numeric)
    
    SSIRaw <- as.data.frame(SSIRaw)
    
    #Gaussian filter
    SSISmooPlot <- smth.gaussian(SSIRaw[,2], window = 17)
    
    SSISmooPlotTime <- SSIRaw[,1]
    
    SSISmooPlotTime <- SSISmooPlotTime[which(!is.na(SSISmooPlot))]
    SSISmooPlot <- SSISmooPlot[which(!is.na(SSISmooPlot))]
    
    # #Smooth the dTSI data at level of 100 years
    # SSISmooPlot <- rollmean(SSIInte, 11, align = 'center')
    # 
    # SSISmooPlotTime <- SSIInteTime[c(-1:-5, -907:-911)]
  }
  
  #Plot SSI
  {
    #SSI data
    WholeSSI <- ggplot() +
      geom_line(data = as.data.frame(SSIRaw), aes(x = Date, y = SSI), 
                size = 1, color = '#8fd7d7', alpha = 0.4) +
      geom_line(aes(x = SSISmooPlotTime, y = SSISmooPlot), 
                size = 1, color = '#8fd7d7', alpha = 1) +
      scale_y_continuous(
        limits = range(SSIRaw[, 2]),
        position = 'right') +
      scale_x_reverse(limits = c(2750, -100), 
                      breaks = seq(from = 4000, to = -100, by = -100),
                      labels = c(4000, rep('', 4),
                                 3500,rep('', 4),
                                 3000, rep('', 4),
                                 2500,rep('', 4),
                                 2000, rep('', 4),
                                 1500,rep('', 4),
                                 1000, rep('', 4),
                                 500,rep('', 4),
                                 0, ''))  +
      labs(y = 'Sea surface temperature', x = 'Cal years BP', title = '') + 
      guides(y = "axis_truncated",x = "axis_truncated") +
      theme_classic() + theme(plot.title = element_blank(), 
                              axis.title.y = element_text(size = 12, family = 'arial'),
                              axis.text.x = element_text(size = 10, family = 'arial'),
                              axis.text.y = element_text(size = 10, family = 'arial'),
                              axis.title.x = element_text(size = 12, family = 'arial')
      )
  }
  
  #Read Atlantic ice-rafting from Bond et. al., 2001
  {
    IRDRaw <- read.csv(file = here('BondIRD.csv'))
    
    colnames(IRDRaw) <- c('Date', 'IRD')
    
    IRDRaw <- apply(IRDRaw, c(1,2), as.numeric)
    
    IRDRaw <- as.data.frame(IRDRaw)
  }
  
  #Plot IRD 
  {
    #IRD data
    WholeIRD <- ggplot() +
      geom_line(data = as.data.frame(IRDRaw), aes(x = Date, y = IRD), 
                size = 1, color = '#bdd373', alpha = 1) +
      scale_y_reverse(
        limits = c(15,0),
        breaks = seq(from = 20, to = 0, by = -5),
        position = 'right') +
      scale_x_reverse(limits = c(2750, -100), 
                      breaks = seq(from = 4000, to = -100, by = -100),
                      labels = c(4000, rep('', 4),
                                 3500,rep('', 4),
                                 3000, rep('', 4),
                                 2500,rep('', 4),
                                 2000, rep('', 4),
                                 1500,rep('', 4),
                                 1000, rep('', 4),
                                 500,rep('', 4),
                                 0, ''))  +
      labs(y = 'Atlantic IRD', x = 'Cal years BP', title = '') + 
      guides(y = "axis_truncated",x = "axis_truncated") +
      theme_classic() + theme(plot.title = element_blank(), 
                              axis.title.y = element_text(size = 12, family = 'arial'),
                              axis.text.x = element_text(size = 10, family = 'arial'),
                              axis.text.y = element_text(size = 10, family = 'arial'),
                              axis.title.x = element_text(size = 12, family = 'arial')
      )
  }
  
  #Read sea ice cover (IP25) data from Cabedo et. al., 2016
  {
    SIRaw <- read.csv(file = here('CabedoSanz2016.csv'))[,c(-1,-4)]
    
    colnames(SIRaw) <- c('Date', 'SI')
    
    SIRaw <- apply(SIRaw, c(1,2), as.numeric)
    
    SIRaw <- as.data.frame(SIRaw)
  }
  
  #Plot sea ice cover data
  {
    #Sea ice cover data
    WholeSI <- ggplot() +
      geom_line(data = as.data.frame(SIRaw), aes(x = Date, y = SI), 
                size = 1, color = '#ff8ca1', alpha = 1) +
      scale_y_reverse(
        limits = rev(range(SIRaw[, 2])),
        position = 'left') +
      scale_x_reverse(limits = c(2750, -100), 
                      breaks = seq(from = 4000, to = -100, by = -100),
                      labels = c(4000, rep('', 4),
                                 3500,rep('', 4),
                                 3000, rep('', 4),
                                 2500,rep('', 4),
                                 2000, rep('', 4),
                                 1500,rep('', 4),
                                 1000, rep('', 4),
                                 500,rep('', 4),
                                 0, ''))  +
      labs(y = 'Sea ice cover', x = 'Cal years BP', title = '') + 
      guides(y = "axis_truncated",x = "axis_truncated") +
      theme_classic() + theme(plot.title = element_blank(), 
                              axis.title.y = element_text(size = 12, family = 'arial'),
                              axis.text.x = element_text(size = 10, family = 'arial'),
                              axis.text.y = element_text(size = 10, family = 'arial'),
                              axis.title.x = element_text(size = 12, family = 'arial')
      )
  }
  
  #Plot Figure
  {
    #Combined plot
    {
      ggarrange(radioDatesPlot, 
                histoEvents +
                  theme(
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                  ),
                averageUACStandaloneSph +
                  geom_hline(yintercept = 0, 
                             linetype = 1, 
                             color = 'grey', alpha = 0.6) +
                  theme(
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                  ),
                WholeSSI +
                  theme(
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                  ),
                WholeSI +
                  theme(
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                  ),
                WholeIRD +
                  theme(
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                  ),
                WholetempAnoma,
                heights = c(0.5,0.3,1,1,1,1,1.2),
                widths = c(1,1,1,1,1,1,1),
                ncol = 1, nrow = 7,
                align = "v")
      
      ggarrange(radioDatesPlot, 
                histoEvents +
                  theme(
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                  ),
                averageUACStandaloneAln +
                  geom_hline(yintercept = 0, 
                             linetype = 1, 
                             color = 'grey', alpha = 0.6) +
                  theme(
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                  ),
                WholeSSI +
                  theme(
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                  ),
                WholeSI +
                  theme(
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                  ),
                WholeIRD +
                  theme(
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                  ),
                WholetempAnoma,
                heights = c(0.5,0.3,1,1,1,1,1.2),
                widths = c(1,1,1,1,1,1,1),
                ncol = 1, nrow = 7,
                align = "v")
      
      ggarrange(radioDatesPlot, 
                histoEvents +
                  theme(
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                  ),
                averageUACStandaloneCal +
                  geom_hline(yintercept = 0, 
                             linetype = 1, 
                             color = 'grey', alpha = 0.6) +
                  theme(
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                  ),
                WholeSSI +
                  theme(
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                  ),
                WholeSI +
                  theme(
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                  ),
                WholeIRD +
                  theme(
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                  ),
                WholetempAnoma,
                heights = c(0.5,0.3,1,1,1,1,1.2),
                widths = c(1,1,1,1,1,1,1),
                ncol = 1, nrow = 7,
                align = "v")
    }
  }
}

#Figure 4: comparison with cloud cover data
{
  #Read CC data
  {
    #Read CC data
    CCRaw <- apply(as.matrix(read.csv(here('Cloud.csv'))), c(1,2),
                   as.numeric)
    
    #Transfer tbe BC/AD ages into Cal. years BP
    CCRaw[, 1] <- 1950 - CCRaw[ ,1]
    
    #Rename the columns
    colnames(CCRaw) <- c('Years_BP', 'CC', 'CC20')
    
    #Transfer the matrix into data frame
    CCRaw <- as.data.frame(CCRaw)
    
    #Convert to anomaly
    CCRaw[,3] <- CCRaw[,3] - mean(CCRaw[,3])
  }
  
  #Plot CC data
  {
    #CC data
    WholeCC <- ggplot() +
      geom_line(aes(x = CCRaw[,1], y = CCRaw[,3]), color = 'black', size = 1) +
      scale_y_continuous(breaks = seq(from = min(floor(CCRaw$CC20)), 
                                      to = max(ceil(CCRaw$CC20)), by = 4),
                         limits = range(CCRaw[, 3]),
                         position = 'right') +
      scale_x_reverse(limits = c(1000, -100), 
                      breaks = seq(from = 4000, to = -100, by = -100),
                      labels = c(4000, rep('', 4),
                                 3500,rep('', 4),
                                 3000, rep('', 4),
                                 2500,rep('', 4),
                                 2000, rep('', 4),
                                 1500,rep('', 4),
                                 1000, rep('', 4),
                                 500,rep('', 4),
                                 0, ''),
                      sec.axis=sec_axis(~., 
                                        breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                        labels=c(2000, rep('', 4),
                                                 1500,rep('', 4),
                                                 1000, rep('', 4),
                                                 500,rep('', 4),
                                                 0, rep('', 4),
                                                 500,rep('', 4),
                                                 1000, rep('', 4),
                                                 1500,rep('', 4),
                                                 2000)
                      )
                      
      )  +
      labs(y = 'Cloud cover anomaly (%)', x = 'Cal years BP', title = '') + 
      guides(y = "axis_truncated", x ="axis_truncated") +
      theme_classic() + theme(plot.title = element_blank(), 
                              axis.title.y = element_text(size = 12, family = 'arial'),
                              axis.text.x = element_text(size = 10, family = 'arial'),
                              axis.text.y = element_text(size = 10, family = 'arial'),
                              axis.title.x = element_text(size = 12, family = 'arial')
      )
  }
  
  #Plot Figure
  {
    ggarrange(radioDatesPlot, 
              averageUACStandaloneSph +
                scale_x_reverse(limits = c(1000, -100), 
                                breaks = seq(from = 4000, to = -100, by = -100),
                                labels = c(4000, rep('', 4),
                                           3500,rep('', 4),
                                           3000, rep('', 4),
                                           2500,rep('', 4),
                                           2000, rep('', 4),
                                           1500,rep('', 4),
                                           1000, rep('', 4),
                                           500,rep('', 4),
                                           0, '')) +
                theme(
                  axis.line.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                ),
              WholeCC +
                geom_hline(yintercept = 0, 
                           linetype = 1, 
                           color = 'grey', alpha = 0.6),
              heights = c(0.3,1,1.2),
              widths = c(1,1,1),
              ncol = 1, nrow = 3,
              align = "v")
    
    ggarrange(radioDatesPlot, 
              averageUACStandaloneAln +
                scale_x_reverse(limits = c(1000, -100), 
                                breaks = seq(from = 4000, to = -100, by = -100),
                                labels = c(4000, rep('', 4),
                                           3500,rep('', 4),
                                           3000, rep('', 4),
                                           2500,rep('', 4),
                                           2000, rep('', 4),
                                           1500,rep('', 4),
                                           1000, rep('', 4),
                                           500,rep('', 4),
                                           0, '')) +
                theme(
                  axis.line.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                ),
              WholeCC +
                geom_hline(yintercept = 0, 
                           linetype = 1, 
                           color = 'grey', alpha = 0.6),
              heights = c(0.3,1,1.2),
              widths = c(1,1,1),
              ncol = 1, nrow = 3,
              align = "v")
    
    ggarrange(radioDatesPlot, 
              averageUACStandaloneCal +
                scale_x_reverse(limits = c(1000, -100), 
                                breaks = seq(from = 4000, to = -100, by = -100),
                                labels = c(4000, rep('', 4),
                                           3500,rep('', 4),
                                           3000, rep('', 4),
                                           2500,rep('', 4),
                                           2000, rep('', 4),
                                           1500,rep('', 4),
                                           1000, rep('', 4),
                                           500,rep('', 4),
                                           0, '')) +
                theme(
                  axis.line.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                ),
              WholeCC +
                geom_hline(yintercept = 0, 
                           linetype = 1, 
                           color = 'grey', alpha = 0.6),
              heights = c(0.3,1,1.2),
              widths = c(1,1,1),
              ncol = 1, nrow = 3,
              align = "v")
  }
}

#Figure 5: wavelet analysis
{
  #Sphagnum
  {
    DataSph <- cbind(UACSph[,2] - chronoCorrMeanMatSca, UACSphMean)
    
    #Remove date older than 5377 years BP
    DataSph <- DataSph[which(DataSph[, 1] <= 2600),]
    
    DataSph <- as.data.frame(DataSph)
    
    colnames(DataSph) <- c('Cal_years_BP', 'UACArea')
    
    Seq <- seq(from = min(DataSph$Cal_years_BP), 
               to = max(DataSph$Cal_years_BP), 
               by = 20) #define a regular series of ages as target for interpolation
    Age <- DataSph$Cal_years_BP
    DA <- DataSph$UACArea
    Interp <- approx(Age, DA, Seq, method = "linear") # approximate the values of the proxy at target ages by linear interpolation
    Age.int <- Interp$x #the new ages, same as "Seq"
    DA.int <- Interp$y #the interpolated DA values
    
    #prepare a data frame with the data to be analysed
    DA.frame <- data.frame(Age = Age.int , Data = DA.int )
    
    #run the wavelet calculations
    DA.wavSph <- analyze.wavelet(DA.frame, #specify the data frame
                                 "Data", #specify the name of the data
                                 loess.span = 0, #performs a detrending by loess, default is span of 0.75., zero = no detrending
                                 dt = 20, #sampling resolution in time domain, i.e. years between samples
                                 dj = 1/20, #sampling resolution in frequency domain
                                 make.pval = T, #compute p-values?
                                 method = "AR", #method for generating null hypothesis, here AR(1) red noise
                                 n.sim = 999) # how many simulations for p-value
  }
  
  #Alnus
  {
    DataAln <- cbind(UACAln[,2] - chronoCorrMeanMatSca, UACAlnMean)
    
    #Remove date older than 5377 years BP
    DataAln <- DataAln[which(DataAln[, 1] <= 2600),]
    
    DataAln <- as.data.frame(DataAln)
    
    colnames(DataAln) <- c('Cal_years_BP', 'UACArea')
    
    Seq <- seq(from = min(DataAln$Cal_years_BP), 
               to = max(DataAln$Cal_years_BP), 
               by = 20) #define a regular series of ages as target for interpolation
    Age <- DataAln$Cal_years_BP
    DA <- DataAln$UACArea
    Interp <- approx(Age, DA, Seq, method = "linear") # approximate the values of the proxy at target ages by linear interpolation
    Age.int <- Interp$x #the new ages, same as "Seq"
    DA.int <- Interp$y #the interpolated DA values
    
    #prepare a data frame with the data to be analysed
    DA.frame <- data.frame(Age = Age.int , Data = DA.int )
    
    #run the wavelet calculations
    DA.wavAln <- analyze.wavelet(DA.frame, #specify the data frame
                                 "Data", #specify the name of the data
                                 loess.span = 0, #performs a detrending by loess, default is span of 0.75., zero = no detrending
                                 dt = 20, #sampling resolution in time domain, i.e. years between samples
                                 dj = 1/20, #sampling resolution in frequency domain
                                 make.pval = T, #compute p-values?
                                 method = "AR", #method for generating null hypothesis, here AR(1) red noise
                                 n.sim = 999) # how many simulations for p-value
  }
  
  #Calluna
  {
    DataCal <- cbind(UACCal[,2] - chronoCorrMeanMatSca, UACCalMean)
    
    #Remove date older than 5377 years BP
    DataCal <- DataCal[which(DataCal[, 1] <= 2600),]
    
    DataCal <- as.data.frame(DataCal)
    
    colnames(DataCal) <- c('Cal_years_BP', 'UACArea')
    
    Seq <- seq(from = min(DataCal$Cal_years_BP), 
               to = max(DataCal$Cal_years_BP), 
               by = 20) #define a regular series of ages as target for interpolation
    Age <- DataCal$Cal_years_BP
    DA <- DataCal$UACArea
    Interp <- approx(Age, DA, Seq, method = "linear") # approximate the values of the proxy at target ages by linear interpolation
    Age.int <- Interp$x #the new ages, same as "Seq"
    DA.int <- Interp$y #the interpolated DA values
    
    #prepare a data frame with the data to be analysed
    DA.frame <- data.frame(Age = Age.int , Data = DA.int )
    
    #run the wavelet calculations
    DA.wavCal <- analyze.wavelet(DA.frame, #specify the data frame
                                 "Data", #specify the name of the data
                                 loess.span = 0, #performs a detrending by loess, default is span of 0.75., zero = no detrending
                                 dt = 20, #sampling resolution in time domain, i.e. years between samples
                                 dj = 1/20, #sampling resolution in frequency domain
                                 make.pval = T, #compute p-values?
                                 method = "AR", #method for generating null hypothesis, here AR(1) red noise
                                 n.sim = 999) # how many simulations for p-value
  }
  
  #dTSI
  {
    DataTSI <- TSIRaw
    
    #Remove date older than 3600 years BP
    DataTSI <- DataTSI[which(DataTSI[, 1] <= 2600),]
    
    colnames(DataTSI) <- c('year_BP', 'dTSI')
    
    Seq <- seq(from = min(DataTSI$year_BP), to = max(DataTSI$year_BP), by = 20) #define a regular series of ages as target for interpolation
    Age <- DataTSI$year_BP[which(DataTSI$year_BP <= 2600)]
    DA <- DataTSI$dTSI[which(DataTSI$year_BP <= 2600)]
    Interp <- approx(Age, DA, Seq, method = "linear") # approximate the values of the proxy at target ages by linear interpolation
    Age.int <- Interp$x #the new ages, same as "Seq"
    DA.int <- Interp$y #the interpolated DA values
    
    #prepare a data frame with the data to be analysed
    DA.frame <- data.frame(Age = Age.int , Data = DA.int)
    
    #run the wavelet calculations
    DA.wavTSI <- analyze.wavelet(DA.frame, #specify the data frame
                                 "Data", #specify the name of the data
                                 loess.span = 0, #performs a detrending by loess, default is span of 0.75., zero = no detrending
                                 dt = 20, #sampling resolution in time domain, i.e. years between samples
                                 dj = 1/20, #sampling resolution in frequency domain
                                 make.pval = T, #compute p-values?
                                 method = "AR", #method for generating null hypothesis, here AR(1) red noise
                                 n.sim = 999) # how many simulations for p-value
  }
  
  #Wavelet plot
  {
    #Sphagnum
    WASph <-
      wt.image(DA.wavSph, 
               main="Sphagnum", 
               color.key = "quantile",
               color.palette = "rev(scico(n.levels, palette = 'batlow'))",
               siglvl = 0.1, #i.e. 95% confidence level
               n.levels = 100, 
               timelab = "Cal yr BP", 
               legend.params = list(lab = "wavelet power levels", mar=4.7, label.digits = 4), 
               plot.ridge = F,
               spec.time.axis = list(at = (seq(from = -100, to = 2750, by = 100) -
                                             min(DataSph$Cal_years_BP))/20 + 1, 
                                     labels = c('','0', rep(' ', 9),
                                                '1000', rep(' ', 9),
                                                '2000', rep(' ', 7)),
                                     las = 2, hadj = NA, padj = NA),
               # spec.period.axis = list(at = c(64, 128, 256, 512, 1024), 
               #                         las = 1, hadj = NA, padj = NA),
               graphics.reset = FALSE
      ) 
    abline(h = log2(800), lwd = 2, lty = 2, col = '#00000050')
    text(20,log2(800),"800", cex = 1)
    abline(h = log2(700), lwd = 2, lty = 2, col = '#00000050')
    text(20,log2(700),"700", cex = 1)
    abline(h = log2(350), lwd = 2, lty = 2, col = '#00000050')
    text(20,log2(350),"350", cex = 1)
    abline(h = log2(204), lwd = 2, lty = 2, col = '#00000050')
    text(20,log2(204),"204", cex = 1)
    abline(h = log2(103), lwd = 2, lty = 2, col = '#00000050')
    text(20,log2(103),"103", cex = 1)
    
    #Sphagnum
    WAAln <-
      wt.image(DA.wavAln, 
               main="Alnus", 
               color.key = "quantile",
               color.palette = "rev(scico(n.levels, palette = 'batlow'))",
               siglvl = 0.1, #i.e. 95% confidence level
               n.levels = 100, 
               timelab = "Cal yr BP", 
               legend.params = list(lab = "wavelet power levels", mar=4.7, label.digits = 4), 
               plot.ridge = F,
               spec.time.axis = list(at = (seq(from = -100, to = 2750, by = 100) -
                                             min(DataAln$Cal_years_BP))/20 + 1, 
                                     labels = c('','0', rep(' ', 9),
                                                '1000', rep(' ', 9),
                                                '2000', rep(' ', 7)),
                                     las = 2, hadj = NA, padj = NA),
               # spec.period.axis = list(at = c(64, 128, 256, 512, 1024), 
               #                         las = 1, hadj = NA, padj = NA),
               graphics.reset = FALSE
      ) 
    abline(h = log2(800), lwd = 2, lty = 2, col = '#00000050')
    text(20,log2(800),"800", cex = 1)
    abline(h = log2(700), lwd = 2, lty = 2, col = '#00000050')
    text(20,log2(700),"700", cex = 1)
    abline(h = log2(350), lwd = 2, lty = 2, col = '#00000050')
    text(20,log2(350),"350", cex = 1)
    abline(h = log2(204), lwd = 2, lty = 2, col = '#00000050')
    text(20,log2(204),"204", cex = 1)
    abline(h = log2(103), lwd = 2, lty = 2, col = '#00000050')
    text(20,log2(103),"103", cex = 1)
    
    #Calluna
    WACal <-
      wt.image(DA.wavCal, 
               main="Calluna", 
               color.key = "quantile",
               color.palette = "rev(scico(n.levels, palette = 'batlow'))",
               siglvl = 0.1, #i.e. 95% confidence level
               n.levels = 100, 
               timelab = "Cal yr BP", 
               legend.params = list(lab = "wavelet power levels", mar=4.7, label.digits = 4), 
               plot.ridge = F,
               spec.time.axis = list(at = (seq(from = -100, to = 2750, by = 100) -
                                             min(DataCal$Cal_years_BP))/20 + 1, 
                                     labels = c('','0', rep(' ', 9),
                                                '1000', rep(' ', 9),
                                                '2000', rep(' ', 7)),
                                     las = 2, hadj = NA, padj = NA),
               # spec.period.axis = list(at = c(64, 128, 256, 512, 1024), 
               #                         las = 1, hadj = NA, padj = NA),
               graphics.reset = FALSE
      ) 
    abline(h = log2(800), lwd = 2, lty = 2, col = '#00000050')
    text(20,log2(800),"800", cex = 1)
    abline(h = log2(700), lwd = 2, lty = 2, col = '#00000050')
    text(20,log2(700),"700", cex = 1)
    abline(h = log2(350), lwd = 2, lty = 2, col = '#00000050')
    text(20,log2(350),"350", cex = 1)
    abline(h = log2(204), lwd = 2, lty = 2, col = '#00000050')
    text(20,log2(204),"204", cex = 1)
    abline(h = log2(103), lwd = 2, lty = 2, col = '#00000050')
    text(20,log2(103),"103", cex = 1)
    
    #write the wavelet diagram to pdf
    WATSI <-
      wt.image(DA.wavTSI, 
               main="dTSI", 
               color.key = "quantile",
               color.palette = "rev(scico(n.levels, palette = 'batlow'))",
               siglvl = 0.1, #i.e. 95% confidence level
               n.levels = 100, 
               timelab = "Cal yr BP", 
               legend.params = list(lab = "wavelet power levels", mar=4.7, label.digits = 4), 
               plot.ridge = F,
               spec.time.axis = list(at = (seq(from = 0, to = 4000, by = 100) -
                                             min(DataSph$Cal_years_BP))/20 + 1, 
                                     labels = c(0, rep('', 9),
                                                1000, rep('', 9),
                                                2000, rep('', 9),
                                                3000, rep('', 9),
                                                4000),
                                     las = 2, hadj = NA, padj = NA),
               # spec.period.axis = list(at = c(64, 128, 256, 512, 1024), 
               #                         las = 1, hadj = NA, padj = NA),
               graphics.reset = FALSE
      )
    abline(h = log2(800), lwd = 2, lty = 2, col = '#00000050')
    text(20,log2(800),"800", cex = 1)
    abline(h = log2(700), lwd = 2, lty = 2, col = '#00000050')
    text(20,log2(700),"700", cex = 1)
    abline(h = log2(350), lwd = 2, lty = 2, col = '#00000050')
    text(20,log2(350),"350", cex = 1)
    abline(h = log2(204), lwd = 2, lty = 2, col = '#00000050')
    text(20,log2(204),"204", cex = 1)
    abline(h = log2(103), lwd = 2, lty = 2, col = '#00000050')
    text(20,log2(103),"103", cex = 1)
  }
}

#Figure 6: UAC with age uncertainty
{
  #Loading lipd data
  {
    UACLipdSph <- lipdR::readLipd(here('Archives_QSR', 'UAC_Sph.lpd'))
    
    UACLipdAln <- lipdR::readLipd(here('Archives_QSR', 'UAC_Aln.lpd'))
    
    UACLipdCal <- lipdR::readLipd(here('Archives_QSR', 'UAC_Cal.lpd'))
  }
  
  #Age model creation
  {
    #Sphagnum
    {
      UACLipdSphWithChron <- runBacon(UACLipdSph, lab.id.var = 'AMS_Number',
                                      bacon.dir = here(),
                                      age.14c.uncertainty.var = 'age14Cuncertainty', 
                                      age.var = 'Calibrated_Ages', 
                                      age.uncertainty.var = 'Calibrated_Ages_Var',
                                      reservoir.age.14c.var = NULL, 
                                      reservoir.age.14c.uncertainty.var = NULL, 
                                      rejected.ages.var = NULL,
                                      bacon.acc.mean = 20,
                                      cc = 1,
                                      bacon.thick = 10,
                                      suggest = FALSE,
                                      max.ens = 1000,
                                      accept.suggestions = TRUE)
    }
    
    #Alnus
    {
      UACLipdAlnWithChron <- runBacon(UACLipdAln, lab.id.var = 'AMS_Number',
                                      bacon.dir = here(),
                                      age.14c.uncertainty.var = 'age14Cuncertainty', 
                                      age.var = 'Calibrated_Ages', 
                                      age.uncertainty.var = 'Calibrated_Ages_Var',
                                      reservoir.age.14c.var = NULL, 
                                      reservoir.age.14c.uncertainty.var = NULL, 
                                      rejected.ages.var = NULL,
                                      bacon.acc.mean = 20,
                                      cc = 1,
                                      bacon.thick = 10,
                                      suggest = FALSE,
                                      max.ens = 1000,
                                      accept.suggestions = TRUE)
    }

    #Calluna
    {
      UACLipdCalWithChron <- runBacon(UACLipdCal, lab.id.var = 'AMS_Number',
                                      bacon.dir = here(),
                                      age.14c.uncertainty.var = 'age14Cuncertainty', 
                                      age.var = 'Calibrated_Ages', 
                                      age.uncertainty.var = 'Calibrated_Ages_Var',
                                      reservoir.age.14c.var = NULL, 
                                      reservoir.age.14c.uncertainty.var = NULL, 
                                      rejected.ages.var = NULL,
                                      bacon.acc.mean = 20,
                                      cc = 1,
                                      bacon.thick = 10,
                                      suggest = FALSE,
                                      max.ens = 1000,
                                      accept.suggestions = TRUE)
    }
  }
  
  #Plot age-depth model
  {
    #Sphagnum
    {
      chronPlotSph <- plotChronEns(UACLipdSphWithChron,truncate.dist = 1e-6)+ggtitle(NULL)+coord_cartesian(xlim = c(6000,-100))
      print(chronPlotSph)
    }
  }
  
  #Matching age-depth model to the proxy dataset, and basic plot
  {
    #Sphagnum
    {
      UACMatchingSph <- mapAgeEnsembleToPaleoData(UACLipdSphWithChron,age.var = "ageEnsemble",
                                                  paleo.depth.var = "Depth",)
      
      UACMatchingSph.ae <- selectData(UACMatchingSph,var.name = "ageEnsemble")
      
      
      UACMatchingSph.PA <- selectData(UACMatchingSph,var.name = "Peak_Area")
      
      
      UACMatchingSph.PA.ts.plot <- plotTimeseriesEnsRibbons(X = UACMatchingSph.ae,
                                                            Y = UACMatchingSph.PA,
                                                            n.bins = 1000) +
        geom_point(
          aes(x = UACMatchingSph[["paleoData"]][[1]][["measurementTable"]][[1]][["ageMedian"]][["values"]],
              y = UACMatchingSph[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]]), 
          size = 2, shape = 21, fill = 'red', alpha = 0.4) + 
        geom_errorbar(aes(x = UACMatchingSph[["paleoData"]][[1]][["measurementTable"]][[1]][["ageMedian"]][["values"]],
                          ymin = UACMatchingSph[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] - 
                            UACMatchingSph[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]],
                          ymax = UACMatchingSph[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] + 
                            UACMatchingSph[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]]),
                      size = 0.5, alpha = 0.4, color = 'black', width = 20) +
        scale_y_continuous(position = 'left', limits = range(UACMatchingSph.PA$values,
                                                             UACMatchingSph[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] - 
                                                               UACMatchingSph[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]],
                                                             UACMatchingSph[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] + 
                                                               UACMatchingSph[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]]),
                           breaks = seq(from = -10, to = 10, by = 0.5)) +
        scale_x_reverse(limits = c(2750, -100), 
                        breaks = seq(from = 4000, to = -100, by = -100),
                        labels = c(4000, rep('', 4),
                                   3500,rep('', 4),
                                   3000, rep('', 4),
                                   2500,rep('', 4),
                                   2000, rep('', 4),
                                   1500,rep('', 4),
                                   1000, rep('', 4),
                                   500,rep('', 4),
                                   0, ''),
                        position = 'top',
                        sec.axis=sec_axis(~., 
                                          breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                          labels=c(2000, rep('', 4),
                                                   1500,rep('', 4),
                                                   1000, rep('', 4),
                                                   500,rep('', 4),
                                                   0, rep('', 4),
                                                   500,rep('', 4),
                                                   1000, rep('', 4),
                                                   1500,rep('', 4),
                                                   2000)
                        )
        ) +
        guides(y = "axis_truncated") +
        labs(y = 'Average UACs', x = 'Cal years BP', title = '') + 
        theme_classic() + 
        theme(plot.title = element_blank(), 
              axis.text.y = element_text(size = 10, family = 'arial'),
              axis.title.y = element_text(size = 12, family = 'arial'),
              axis.text.x = element_text(size = 10, family = 'arial'),
              axis.title.x = element_text(size = 12, family = 'arial')
        )
      
      
      UACSphChronUnMean <- mean(UACMatchingSph.PA.ts.plot[["layers"]][[3]][["data"]][["y"]], 
                                na.rm = TRUE)
      
      UACSphChronUnSD <- sd(UACMatchingSph.PA.ts.plot[["layers"]][[3]][["data"]][["y"]], 
                                na.rm = TRUE)
    }
    
    #Alnus
    {
      UACMatchingAln <- mapAgeEnsembleToPaleoData(UACLipdAlnWithChron,age.var = "ageEnsemble")
      
      UACMatchingAln.ae <- selectData(UACMatchingAln,var.name = "ageEnsemble")
      
      UACMatchingAln.PA <- selectData(UACMatchingAln,var.name = "Peak_Area")
      
      UACMatchingAln.PA.ts.plot <- plotTimeseriesEnsRibbons(X = UACMatchingAln.ae,
                                                            Y = UACMatchingAln.PA,
                                                            n.bins = 1000) +
        geom_point(
          aes(x = UACMatchingAln[["paleoData"]][[1]][["measurementTable"]][[1]][["ageMedian"]][["values"]],
              y = UACMatchingAln[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]]), 
          size = 2, shape = 21, fill = 'red', alpha = 0.4) + 
        geom_errorbar(aes(x = UACMatchingAln[["paleoData"]][[1]][["measurementTable"]][[1]][["ageMedian"]][["values"]],
                          ymin = UACMatchingAln[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] - 
                            UACMatchingAln[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]],
                          ymax = UACMatchingAln[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] + 
                            UACMatchingAln[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]]),
                      size = 0.5, alpha = 0.4, color = 'black', width = 20) +
        scale_y_continuous(position = 'left', limits = range(UACMatchingAln.PA$values,
                                                             UACMatchingAln[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] - 
                                                               UACMatchingAln[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]],
                                                             UACMatchingAln[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] + 
                                                               UACMatchingAln[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]]),
                           breaks = seq(from = -10, to = 10, by = 0.5)) +
        scale_x_reverse(limits = c(2750, -100), 
                        breaks = seq(from = 4000, to = -100, by = -100),
                        labels = c(4000, rep('', 4),
                                   3500,rep('', 4),
                                   3000, rep('', 4),
                                   2500,rep('', 4),
                                   2000, rep('', 4),
                                   1500,rep('', 4),
                                   1000, rep('', 4),
                                   500,rep('', 4),
                                   0, ''),
                        position = 'top',
                        sec.axis=sec_axis(~., 
                                          breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                          labels=c(2000, rep('', 4),
                                                   1500,rep('', 4),
                                                   1000, rep('', 4),
                                                   500,rep('', 4),
                                                   0, rep('', 4),
                                                   500,rep('', 4),
                                                   1000, rep('', 4),
                                                   1500,rep('', 4),
                                                   2000)
                        )
        ) +
        guides(y = "axis_truncated") +
        labs(y = 'Average UACs', x = 'Cal years BP', title = '') + 
        theme_classic() + 
        theme(plot.title = element_blank(), 
              axis.text.y = element_text(size = 10, family = 'arial'),
              axis.title.y = element_text(size = 12, family = 'arial'),
              axis.text.x = element_text(size = 10, family = 'arial'),
              axis.title.x = element_text(size = 12, family = 'arial')
        )
      
      UACAlnChronUnMean <- mean(UACMatchingAln.PA.ts.plot$layers[[3]]$data$y, 
                                na.rm = TRUE)
      
      UACAlnChronUnSD <- sd(UACMatchingAln.PA.ts.plot$layers[[3]]$data$y, 
                            na.rm = TRUE)
    }
    
    #Calluna
    {
      UACMatchingCal <- mapAgeEnsembleToPaleoData(UACLipdCalWithChron,age.var = "ageEnsemble")
      
      UACMatchingCal.ae <- selectData(UACMatchingCal,var.name = "ageEnsemble")
      
      UACMatchingCal.PA <- selectData(UACMatchingCal,var.name = "Peak_Area")
      
      UACMatchingCal.PA.ts.plot <- plotTimeseriesEnsRibbons(X = UACMatchingCal.ae,
                                                            Y = UACMatchingCal.PA,
                                                            n.bins = 1000) +
        geom_point(
          aes(x = UACMatchingCal[["paleoData"]][[1]][["measurementTable"]][[1]][["ageMedian"]][["values"]],
              y = UACMatchingCal[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]]), 
          size = 2, shape = 21, fill = 'red', alpha = 0.4) + 
        geom_errorbar(aes(x = UACMatchingCal[["paleoData"]][[1]][["measurementTable"]][[1]][["ageMedian"]][["values"]],
                          ymin = UACMatchingCal[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] - 
                            UACMatchingCal[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]],
                          ymax = UACMatchingCal[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] + 
                            UACMatchingCal[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]]),
                      size = 0.5, alpha = 0.4, color = 'black', width = 20) +
        scale_y_continuous(position = 'left', limits = range(UACMatchingCal.PA$values,
                                                             UACMatchingCal[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] - 
                                                               UACMatchingCal[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]],
                                                             UACMatchingCal[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] + 
                                                               UACMatchingCal[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]]),
                           breaks = seq(from = -10, to = 10, by = 0.5)) +
        scale_x_reverse(limits = c(2750, -100), 
                        breaks = seq(from = 4000, to = -100, by = -100),
                        labels = c(4000, rep('', 4),
                                   3500,rep('', 4),
                                   3000, rep('', 4),
                                   2500,rep('', 4),
                                   2000, rep('', 4),
                                   1500,rep('', 4),
                                   1000, rep('', 4),
                                   500,rep('', 4),
                                   0, ''),
                        position = 'top',
                        sec.axis=sec_axis(~., 
                                          breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                          labels=c(2000, rep('', 4),
                                                   1500,rep('', 4),
                                                   1000, rep('', 4),
                                                   500,rep('', 4),
                                                   0, rep('', 4),
                                                   500,rep('', 4),
                                                   1000, rep('', 4),
                                                   1500,rep('', 4),
                                                   2000)
                        )
        ) +
        guides(y = "axis_truncated") +
        labs(y = 'Average UACs', x = 'Cal years BP', title = '') + 
        theme_classic() + 
        theme(plot.title = element_blank(), 
              axis.text.y = element_text(size = 10, family = 'arial'),
              axis.title.y = element_text(size = 12, family = 'arial'),
              axis.text.x = element_text(size = 10, family = 'arial'),
              axis.title.x = element_text(size = 12, family = 'arial')
        )
      
      UACCalChronUnMean <- mean(UACMatchingCal.PA.ts.plot$layers[[3]]$data$y, 
                                na.rm = TRUE)
      
      UACCalChronUnSD <- sd(UACMatchingCal.PA.ts.plot$layers[[3]]$data$y, 
                            na.rm = TRUE)
    }
    
    #Combine all plots
    {
      ggarrange(UACMatchingSph.PA.ts.plot+
                  geom_hline(yintercept = UACSphChronUnMean, 
                             linetype = 1, 
                             color = 'blue', alpha = 0.5) +
                  geom_hline(yintercept = UACSphChronUnMean - UACSphChronUnSD, 
                             linetype = 2, 
                             color = 'blue', alpha = 0.5) +
                  geom_hline(yintercept = UACSphChronUnMean - 2 * UACSphChronUnSD, 
                             linetype = 2, 
                             color = 'blue', alpha = 0.5) +
                  geom_hline(yintercept = UACSphChronUnMean + UACSphChronUnSD, 
                             linetype = 2, 
                             color = 'blue', alpha = 0.5) +
                  geom_hline(yintercept = UACSphChronUnMean + 2* UACSphChronUnSD, 
                             linetype = 2, 
                             color = 'blue', alpha = 0.5) ,
                UACMatchingAln.PA.ts.plot+
                  geom_hline(yintercept = UACAlnChronUnMean, 
                             linetype = 1, 
                             color = 'blue', alpha = 0.5) +
                  geom_hline(yintercept = UACAlnChronUnMean - UACAlnChronUnSD, 
                             linetype = 2, 
                             color = 'blue', alpha = 0.5) +
                  geom_hline(yintercept = UACAlnChronUnMean - 2 * UACAlnChronUnSD, 
                             linetype = 2, 
                             color = 'blue', alpha = 0.5) +
                  geom_hline(yintercept = UACAlnChronUnMean + UACAlnChronUnSD, 
                             linetype = 2, 
                             color = 'blue', alpha = 0.5) +
                  geom_hline(yintercept = UACAlnChronUnMean + 2* UACAlnChronUnSD, 
                             linetype = 2, 
                             color = 'blue', alpha = 0.5) +
                  theme(
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                  ),
                UACMatchingCal.PA.ts.plot+
                  geom_hline(yintercept = UACCalChronUnMean, 
                             linetype = 1, 
                             color = 'blue', alpha = 0.5) +
                  geom_hline(yintercept = UACCalChronUnMean - UACCalChronUnSD, 
                             linetype = 2, 
                             color = 'blue', alpha = 0.5) +
                  geom_hline(yintercept = UACCalChronUnMean - 2 * UACCalChronUnSD, 
                             linetype = 2, 
                             color = 'blue', alpha = 0.5) +
                  geom_hline(yintercept = UACCalChronUnMean + UACCalChronUnSD, 
                             linetype = 2, 
                             color = 'blue', alpha = 0.5) +
                  geom_hline(yintercept = UACCalChronUnMean + 2* UACCalChronUnSD, 
                             linetype = 2, 
                             color = 'blue', alpha = 0.5) +
                  theme(
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                  ),
                WholeTSI,
                heights = c(1.2,1,1,1.2),
                ncol = 1,
                nrow = 4,
                align = 'v')
      
      ggarrange(UACMatchingSph.PA.ts.plot +
                  geom_hline(yintercept = UACSphChronUnMean, 
                             linetype = 1, 
                             color = 'blue', alpha = 0.5) +
                  geom_hline(yintercept = UACSphChronUnMean - UACSphChronUnSD, 
                             linetype = 2, 
                             color = 'blue', alpha = 0.5) +
                  geom_hline(yintercept = UACSphChronUnMean - 2 * UACSphChronUnSD, 
                             linetype = 2, 
                             color = 'blue', alpha = 0.5) +
                  geom_hline(yintercept = UACSphChronUnMean + UACSphChronUnSD, 
                             linetype = 2, 
                             color = 'blue', alpha = 0.5) +
                  geom_hline(yintercept = UACSphChronUnMean + 2* UACSphChronUnSD, 
                             linetype = 2, 
                             color = 'blue', alpha = 0.5) +
                  scale_x_reverse(limits = c(1000, -100), 
                                  breaks = seq(from = 4000, to = -100, by = -100),
                                  labels = c(4000, rep('', 4),
                                             3500,rep('', 4),
                                             3000, rep('', 4),
                                             2500,rep('', 4),
                                             2000, rep('', 4),
                                             1500,rep('', 4),
                                             1000, rep('', 4),
                                             500,rep('', 4),
                                             0, ''),
                                  sec.axis=sec_axis(~., 
                                                    breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                                    labels=c(2000, rep('', 4),
                                                             1500,rep('', 4),
                                                             1000, rep('', 4),
                                                             500,rep('', 4),
                                                             0, rep('', 4),
                                                             500,rep('', 4),
                                                             1000, rep('', 4),
                                                             1500,rep('', 4),
                                                             2000)
                                    ) 
                                  ),
                UACMatchingAln.PA.ts.plot +
                  geom_hline(yintercept = UACAlnChronUnMean, 
                             linetype = 1, 
                             color = 'blue', alpha = 0.5) +
                  geom_hline(yintercept = UACAlnChronUnMean - UACAlnChronUnSD, 
                             linetype = 2, 
                             color = 'blue', alpha = 0.5) +
                  geom_hline(yintercept = UACAlnChronUnMean - 2 * UACAlnChronUnSD, 
                             linetype = 2, 
                             color = 'blue', alpha = 0.5) +
                  geom_hline(yintercept = UACAlnChronUnMean + UACAlnChronUnSD, 
                             linetype = 2, 
                             color = 'blue', alpha = 0.5) +
                  geom_hline(yintercept = UACAlnChronUnMean + 2* UACAlnChronUnSD, 
                             linetype = 2, 
                             color = 'blue', alpha = 0.5) +
                  scale_x_reverse(limits = c(1000, -100), 
                                  breaks = seq(from = 4000, to = -100, by = -100),
                                  labels = c(4000, rep('', 4),
                                             3500,rep('', 4),
                                             3000, rep('', 4),
                                             2500,rep('', 4),
                                             2000, rep('', 4),
                                             1500,rep('', 4),
                                             1000, rep('', 4),
                                             500,rep('', 4),
                                             0, ''),
                                  sec.axis=sec_axis(~., 
                                                    breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                                    labels=c(2000, rep('', 4),
                                                             1500,rep('', 4),
                                                             1000, rep('', 4),
                                                             500,rep('', 4),
                                                             0, rep('', 4),
                                                             500,rep('', 4),
                                                             1000, rep('', 4),
                                                             1500,rep('', 4),
                                                             2000)
                                    )
                                  )+
                  theme(
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                  ),
                UACMatchingCal.PA.ts.plot+
                  geom_hline(yintercept = UACCalChronUnMean, 
                             linetype = 1, 
                             color = 'blue', alpha = 0.5) +
                  geom_hline(yintercept = UACCalChronUnMean - UACCalChronUnSD, 
                             linetype = 2, 
                             color = 'blue', alpha = 0.5) +
                  geom_hline(yintercept = UACCalChronUnMean - 2 * UACCalChronUnSD, 
                             linetype = 2, 
                             color = 'blue', alpha = 0.5) +
                  geom_hline(yintercept = UACCalChronUnMean + UACCalChronUnSD, 
                             linetype = 2, 
                             color = 'blue', alpha = 0.5) +
                  geom_hline(yintercept = UACCalChronUnMean + 2* UACCalChronUnSD, 
                             linetype = 2, 
                             color = 'blue', alpha = 0.5) +
                  scale_x_reverse(limits = c(1000, -100), 
                                  breaks = seq(from = 4000, to = -100, by = -100),
                                  labels = c(4000, rep('', 4),
                                             3500,rep('', 4),
                                             3000, rep('', 4),
                                             2500,rep('', 4),
                                             2000, rep('', 4),
                                             1500,rep('', 4),
                                             1000, rep('', 4),
                                             500,rep('', 4),
                                             0, ''),
                                  sec.axis=sec_axis(~., 
                                                    breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                                    labels=c(2000, rep('', 4),
                                                             1500,rep('', 4),
                                                             1000, rep('', 4),
                                                             500,rep('', 4),
                                                             0, rep('', 4),
                                                             500,rep('', 4),
                                                             1000, rep('', 4),
                                                             1500,rep('', 4),
                                                             2000)
                                    )
                                  )+
                  theme(
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                  ),
                WholeCC,
                heights = c(1.2,1,1,1.2),
                ncol = 1,
                nrow = 4,
                align = 'v')
    }
    
    #Subgroup
    {
      #Sphagnum
      {
        #1
        {
          UAC_1_MatchingSph <- mapAgeEnsembleToPaleoData(UACLipdSphWithChron,age.var = "ageEnsemble",
                                                         paleo.depth.var = "Depth",
                                                         paleo.depth.range = c(202,278))
          
          UAC_1_MatchingSph.ae <- selectData(UAC_1_MatchingSph,var.name = "ageEnsemble")
          
          
          UAC_1_MatchingSph.PA <- selectData(UAC_1_MatchingSph,var.name = "Peak_Area")
          
          
          UAC_1_MatchingSph.PA.ts.plot <- plotTimeseriesEnsRibbons(X = UAC_1_MatchingSph.ae,
                                                                   Y = UAC_1_MatchingSph.PA,
                                                                   n.bins = 1000) +
            geom_point(
              aes(x = UAC_1_MatchingSph[["paleoData"]][[1]][["measurementTable"]][[1]][["ageMedian"]][["values"]],
                  y = UAC_1_MatchingSph[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]]), 
              size = 2, shape = 21, fill = 'red', alpha = 0.4) + 
            geom_errorbar(aes(x = UAC_1_MatchingSph[["paleoData"]][[1]][["measurementTable"]][[1]][["ageMedian"]][["values"]],
                              ymin = UAC_1_MatchingSph[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] - 
                                UAC_1_MatchingSph[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]],
                              ymax = UAC_1_MatchingSph[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] + 
                                UAC_1_MatchingSph[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]]),
                          size = 0.5, alpha = 0.4, color = 'black', width = 20) +
            scale_y_continuous(position = 'left', limits = range(UAC_1_MatchingSph.PA$values,
                                                                 UAC_1_MatchingSph[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] - 
                                                                   UAC_1_MatchingSph[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]],
                                                                 UAC_1_MatchingSph[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] + 
                                                                   UAC_1_MatchingSph[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]]),
                               breaks = seq(from = -10, to = 10, by = 0.5)) +
            scale_x_reverse(limits = c(2750, -100), 
                            breaks = seq(from = 4000, to = -100, by = -100),
                            labels = c(4000, rep('', 4),
                                       3500,rep('', 4),
                                       3000, rep('', 4),
                                       2500,rep('', 4),
                                       2000, rep('', 4),
                                       1500,rep('', 4),
                                       1000, rep('', 4),
                                       500,rep('', 4),
                                       0, ''),
                            position = 'top',
                            sec.axis=sec_axis(~., 
                                              breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                              labels=c(2000, rep('', 4),
                                                       1500,rep('', 4),
                                                       1000, rep('', 4),
                                                       500,rep('', 4),
                                                       0, rep('', 4),
                                                       500,rep('', 4),
                                                       1000, rep('', 4),
                                                       1500,rep('', 4),
                                                       2000)
                            )
            ) +
            guides(y = "axis_truncated") +
            labs(y = 'Average UAC_1_s', x = 'Cal years BP', title = '') + 
            theme_classic() + 
            theme(plot.title = element_blank(), 
                  axis.text.y = element_text(size = 10, family = 'arial'),
                  axis.title.y = element_text(size = 12, family = 'arial'),
                  axis.text.x = element_text(size = 10, family = 'arial'),
                  axis.title.x = element_text(size = 12, family = 'arial')
            )
          
          UAC_1_SphChronUnMean <- mean(UAC_1_MatchingSph.PA.ts.plot[["layers"]][[3]][["data"]][["y"]], 
                                       na.rm = TRUE)
          
          UAC_1_SphChronUnSD <- sd(UAC_1_MatchingSph.PA.ts.plot[["layers"]][[3]][["data"]][["y"]], 
                                   na.rm = TRUE)
        }
        
        #2
        {
          UAC_2_MatchingSph <- mapAgeEnsembleToPaleoData(UACLipdSphWithChron,age.var = "ageEnsemble",
                                                         paleo.depth.var = "Depth",
                                                         paleo.depth.range = c(96,200))
          
          UAC_2_MatchingSph.ae <- selectData(UAC_2_MatchingSph,var.name = "ageEnsemble")
          
          
          UAC_2_MatchingSph.PA <- selectData(UAC_2_MatchingSph,var.name = "Peak_Area")
          
          
          UAC_2_MatchingSph.PA.ts.plot <- plotTimeseriesEnsRibbons(X = UAC_2_MatchingSph.ae,
                                                                   Y = UAC_2_MatchingSph.PA,
                                                                   n.bins = 1000) +
            geom_point(
              aes(x = UAC_2_MatchingSph[["paleoData"]][[1]][["measurementTable"]][[1]][["ageMedian"]][["values"]],
                  y = UAC_2_MatchingSph[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]]), 
              size = 2, shape = 21, fill = 'red', alpha = 0.4) + 
            geom_errorbar(aes(x = UAC_2_MatchingSph[["paleoData"]][[1]][["measurementTable"]][[1]][["ageMedian"]][["values"]],
                              ymin = UAC_2_MatchingSph[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] - 
                                UAC_2_MatchingSph[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]],
                              ymax = UAC_2_MatchingSph[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] + 
                                UAC_2_MatchingSph[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]]),
                          size = 0.5, alpha = 0.4, color = 'black', width = 20) +
            scale_y_continuous(position = 'left', limits = range(UAC_2_MatchingSph.PA$values,
                                                                 UAC_2_MatchingSph[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] - 
                                                                   UAC_2_MatchingSph[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]],
                                                                 UAC_2_MatchingSph[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] + 
                                                                   UAC_2_MatchingSph[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]]),
                               breaks = seq(from = -10, to = 10, by = 0.5)) +
            scale_x_reverse(limits = c(2750, -100), 
                            breaks = seq(from = 4000, to = -100, by = -100),
                            labels = c(4000, rep('', 4),
                                       3500,rep('', 4),
                                       3000, rep('', 4),
                                       2500,rep('', 4),
                                       2000, rep('', 4),
                                       1500,rep('', 4),
                                       1000, rep('', 4),
                                       500,rep('', 4),
                                       0, ''),
                            position = 'top',
                            sec.axis=sec_axis(~., 
                                              breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                              labels=c(2000, rep('', 4),
                                                       1500,rep('', 4),
                                                       1000, rep('', 4),
                                                       500,rep('', 4),
                                                       0, rep('', 4),
                                                       500,rep('', 4),
                                                       1000, rep('', 4),
                                                       1500,rep('', 4),
                                                       2000)
                            )
            ) +
            guides(y = "axis_truncated") +
            labs(y = 'Average UAC_2_s', x = 'Cal years BP', title = '') + 
            theme_classic() + 
            theme(plot.title = element_blank(), 
                  axis.text.y = element_text(size = 10, family = 'arial'),
                  axis.title.y = element_text(size = 12, family = 'arial'),
                  axis.text.x = element_text(size = 10, family = 'arial'),
                  axis.title.x = element_text(size = 12, family = 'arial')
            )
          
          UAC_2_SphChronUnMean <- mean(UAC_2_MatchingSph.PA.ts.plot[["layers"]][[3]][["data"]][["y"]], 
                                       na.rm = TRUE)
          
          UAC_2_SphChronUnSD <- sd(UAC_2_MatchingSph.PA.ts.plot[["layers"]][[3]][["data"]][["y"]], 
                                   na.rm = TRUE)
        }
        
        #3
        {
          UAC_3_MatchingSph <- mapAgeEnsembleToPaleoData(UACLipdSphWithChron,age.var = "ageEnsemble",
                                                         paleo.depth.var = "Depth",
                                                         paleo.depth.range = c(2,94))
          
          UAC_3_MatchingSph.ae <- selectData(UAC_3_MatchingSph,var.name = "ageEnsemble")
          
          
          UAC_3_MatchingSph.PA <- selectData(UAC_3_MatchingSph,var.name = "Peak_Area")
          
          
          UAC_3_MatchingSph.PA.ts.plot <- plotTimeseriesEnsRibbons(X = UAC_3_MatchingSph.ae,
                                                                   Y = UAC_3_MatchingSph.PA,
                                                                   n.bins = 1000) +
            geom_point(
              aes(x = UAC_3_MatchingSph[["paleoData"]][[1]][["measurementTable"]][[1]][["ageMedian"]][["values"]],
                  y = UAC_3_MatchingSph[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]]), 
              size = 2, shape = 21, fill = 'red', alpha = 0.4) + 
            geom_errorbar(aes(x = UAC_3_MatchingSph[["paleoData"]][[1]][["measurementTable"]][[1]][["ageMedian"]][["values"]],
                              ymin = UAC_3_MatchingSph[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] - 
                                UAC_3_MatchingSph[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]],
                              ymax = UAC_3_MatchingSph[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] + 
                                UAC_3_MatchingSph[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]]),
                          size = 0.5, alpha = 0.4, color = 'black', width = 20) +
            scale_y_continuous(position = 'left', limits = range(UAC_3_MatchingSph.PA$values,
                                                                 UAC_3_MatchingSph[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] - 
                                                                   UAC_3_MatchingSph[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]],
                                                                 UAC_3_MatchingSph[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] + 
                                                                   UAC_3_MatchingSph[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]]),
                               breaks = seq(from = -10, to = 10, by = 0.5)) +
            scale_x_reverse(limits = c(2750, -100), 
                            breaks = seq(from = 4000, to = -100, by = -100),
                            labels = c(4000, rep('', 4),
                                       3500,rep('', 4),
                                       3000, rep('', 4),
                                       2500,rep('', 4),
                                       2000, rep('', 4),
                                       1500,rep('', 4),
                                       1000, rep('', 4),
                                       500,rep('', 4),
                                       0, ''),
                            position = 'top',
                            sec.axis=sec_axis(~., 
                                              breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                              labels=c(2000, rep('', 4),
                                                       1500,rep('', 4),
                                                       1000, rep('', 4),
                                                       500,rep('', 4),
                                                       0, rep('', 4),
                                                       500,rep('', 4),
                                                       1000, rep('', 4),
                                                       1500,rep('', 4),
                                                       2000)
                            )
            ) +
            guides(y = "axis_truncated") +
            labs(y = 'Average UAC_3_s', x = 'Cal years BP', title = '') + 
            theme_classic() + 
            theme(plot.title = element_blank(), 
                  axis.text.y = element_text(size = 10, family = 'arial'),
                  axis.title.y = element_text(size = 12, family = 'arial'),
                  axis.text.x = element_text(size = 10, family = 'arial'),
                  axis.title.x = element_text(size = 12, family = 'arial')
            )
          
          UAC_3_SphChronUnMean <- mean(UAC_3_MatchingSph.PA.ts.plot[["layers"]][[3]][["data"]][["y"]], 
                                       na.rm = TRUE)
          
          UAC_3_SphChronUnSD <- sd(UAC_3_MatchingSph.PA.ts.plot[["layers"]][[3]][["data"]][["y"]], 
                                   na.rm = TRUE)
        }
      }
      
      #Alnus
      {
        #1
        {
          UAC_1_MatchingAln <- mapAgeEnsembleToPaleoData(UACLipdAlnWithChron,age.var = "ageEnsemble",
                                                         paleo.depth.var = "Depth",
                                                         paleo.depth.range = c(202,278))
          
          UAC_1_MatchingAln.ae <- selectData(UAC_1_MatchingAln,var.name = "ageEnsemble")
          
          
          UAC_1_MatchingAln.PA <- selectData(UAC_1_MatchingAln,var.name = "Peak_Area")
          
          
          UAC_1_MatchingAln.PA.ts.plot <- plotTimeseriesEnsRibbons(X = UAC_1_MatchingAln.ae,
                                                                   Y = UAC_1_MatchingAln.PA,
                                                                   n.bins = 1000) +
            geom_point(
              aes(x = UAC_1_MatchingAln[["paleoData"]][[1]][["measurementTable"]][[1]][["ageMedian"]][["values"]],
                  y = UAC_1_MatchingAln[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]]), 
              size = 2, shape = 21, fill = 'red', alpha = 0.4) + 
            geom_errorbar(aes(x = UAC_1_MatchingAln[["paleoData"]][[1]][["measurementTable"]][[1]][["ageMedian"]][["values"]],
                              ymin = UAC_1_MatchingAln[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] - 
                                UAC_1_MatchingAln[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]],
                              ymax = UAC_1_MatchingAln[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] + 
                                UAC_1_MatchingAln[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]]),
                          size = 0.5, alpha = 0.4, color = 'black', width = 20) +
            scale_y_continuous(position = 'left', limits = range(UAC_1_MatchingAln.PA$values,
                                                                 UAC_1_MatchingAln[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] - 
                                                                   UAC_1_MatchingAln[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]],
                                                                 UAC_1_MatchingAln[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] + 
                                                                   UAC_1_MatchingAln[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]]),
                               breaks = seq(from = -10, to = 10, by = 0.5)) +
            scale_x_reverse(limits = c(2750, -100), 
                            breaks = seq(from = 4000, to = -100, by = -100),
                            labels = c(4000, rep('', 4),
                                       3500,rep('', 4),
                                       3000, rep('', 4),
                                       2500,rep('', 4),
                                       2000, rep('', 4),
                                       1500,rep('', 4),
                                       1000, rep('', 4),
                                       500,rep('', 4),
                                       0, ''),
                            position = 'top',
                            sec.axis=sec_axis(~., 
                                              breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                              labels=c(2000, rep('', 4),
                                                       1500,rep('', 4),
                                                       1000, rep('', 4),
                                                       500,rep('', 4),
                                                       0, rep('', 4),
                                                       500,rep('', 4),
                                                       1000, rep('', 4),
                                                       1500,rep('', 4),
                                                       2000)
                            )
            ) +
            guides(y = "axis_truncated") +
            labs(y = 'Average UAC_1_s', x = 'Cal years BP', title = '') + 
            theme_classic() + 
            theme(plot.title = element_blank(), 
                  axis.text.y = element_text(size = 10, family = 'arial'),
                  axis.title.y = element_text(size = 12, family = 'arial'),
                  axis.text.x = element_text(size = 10, family = 'arial'),
                  axis.title.x = element_text(size = 12, family = 'arial')
            )
          
          UAC_1_AlnChronUnMean <- mean(UAC_1_MatchingAln.PA.ts.plot[["layers"]][[3]][["data"]][["y"]], 
                                       na.rm = TRUE)
          
          UAC_1_AlnChronUnSD <- sd(UAC_1_MatchingAln.PA.ts.plot[["layers"]][[3]][["data"]][["y"]], 
                                   na.rm = TRUE)
        }
        
        #2
        {
          UAC_2_MatchingAln <- mapAgeEnsembleToPaleoData(UACLipdAlnWithChron,age.var = "ageEnsemble",
                                                         paleo.depth.var = "Depth",
                                                         paleo.depth.range = c(96,200))
          
          UAC_2_MatchingAln.ae <- selectData(UAC_2_MatchingAln,var.name = "ageEnsemble")
          
          
          UAC_2_MatchingAln.PA <- selectData(UAC_2_MatchingAln,var.name = "Peak_Area")
          
          
          UAC_2_MatchingAln.PA.ts.plot <- plotTimeseriesEnsRibbons(X = UAC_2_MatchingAln.ae,
                                                                   Y = UAC_2_MatchingAln.PA,
                                                                   n.bins = 1000) +
            geom_point(
              aes(x = UAC_2_MatchingAln[["paleoData"]][[1]][["measurementTable"]][[1]][["ageMedian"]][["values"]],
                  y = UAC_2_MatchingAln[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]]), 
              size = 2, shape = 21, fill = 'red', alpha = 0.4) + 
            geom_errorbar(aes(x = UAC_2_MatchingAln[["paleoData"]][[1]][["measurementTable"]][[1]][["ageMedian"]][["values"]],
                              ymin = UAC_2_MatchingAln[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] - 
                                UAC_2_MatchingAln[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]],
                              ymax = UAC_2_MatchingAln[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] + 
                                UAC_2_MatchingAln[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]]),
                          size = 0.5, alpha = 0.4, color = 'black', width = 20) +
            scale_y_continuous(position = 'left', limits = range(UAC_2_MatchingAln.PA$values,
                                                                 UAC_2_MatchingAln[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] - 
                                                                   UAC_2_MatchingAln[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]],
                                                                 UAC_2_MatchingAln[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] + 
                                                                   UAC_2_MatchingAln[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]]),
                               breaks = seq(from = -10, to = 10, by = 0.5)) +
            scale_x_reverse(limits = c(2750, -100), 
                            breaks = seq(from = 4000, to = -100, by = -100),
                            labels = c(4000, rep('', 4),
                                       3500,rep('', 4),
                                       3000, rep('', 4),
                                       2500,rep('', 4),
                                       2000, rep('', 4),
                                       1500,rep('', 4),
                                       1000, rep('', 4),
                                       500,rep('', 4),
                                       0, ''),
                            position = 'top',
                            sec.axis=sec_axis(~., 
                                              breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                              labels=c(2000, rep('', 4),
                                                       1500,rep('', 4),
                                                       1000, rep('', 4),
                                                       500,rep('', 4),
                                                       0, rep('', 4),
                                                       500,rep('', 4),
                                                       1000, rep('', 4),
                                                       1500,rep('', 4),
                                                       2000)
                            )
            ) +
            guides(y = "axis_truncated") +
            labs(y = 'Average UAC_2_s', x = 'Cal years BP', title = '') + 
            theme_classic() + 
            theme(plot.title = element_blank(), 
                  axis.text.y = element_text(size = 10, family = 'arial'),
                  axis.title.y = element_text(size = 12, family = 'arial'),
                  axis.text.x = element_text(size = 10, family = 'arial'),
                  axis.title.x = element_text(size = 12, family = 'arial')
            )
          
          UAC_2_AlnChronUnMean <- mean(UAC_2_MatchingAln.PA.ts.plot[["layers"]][[3]][["data"]][["y"]], 
                                       na.rm = TRUE)
          
          UAC_2_AlnChronUnSD <- sd(UAC_2_MatchingAln.PA.ts.plot[["layers"]][[3]][["data"]][["y"]], 
                                   na.rm = TRUE)
        }
        
        #3
        {
          UAC_3_MatchingAln <- mapAgeEnsembleToPaleoData(UACLipdAlnWithChron,age.var = "ageEnsemble",
                                                         paleo.depth.var = "Depth",
                                                         paleo.depth.range = c(2,94))
          
          UAC_3_MatchingAln.ae <- selectData(UAC_3_MatchingAln,var.name = "ageEnsemble")
          
          
          UAC_3_MatchingAln.PA <- selectData(UAC_3_MatchingAln,var.name = "Peak_Area")
          
          
          UAC_3_MatchingAln.PA.ts.plot <- plotTimeseriesEnsRibbons(X = UAC_3_MatchingAln.ae,
                                                                   Y = UAC_3_MatchingAln.PA,
                                                                   n.bins = 1000) +
            geom_point(
              aes(x = UAC_3_MatchingAln[["paleoData"]][[1]][["measurementTable"]][[1]][["ageMedian"]][["values"]],
                  y = UAC_3_MatchingAln[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]]), 
              size = 2, shape = 21, fill = 'red', alpha = 0.4) + 
            geom_errorbar(aes(x = UAC_3_MatchingAln[["paleoData"]][[1]][["measurementTable"]][[1]][["ageMedian"]][["values"]],
                              ymin = UAC_3_MatchingAln[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] - 
                                UAC_3_MatchingAln[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]],
                              ymax = UAC_3_MatchingAln[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] + 
                                UAC_3_MatchingAln[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]]),
                          size = 0.5, alpha = 0.4, color = 'black', width = 20) +
            scale_y_continuous(position = 'left', limits = range(UAC_3_MatchingAln.PA$values,
                                                                 UAC_3_MatchingAln[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] - 
                                                                   UAC_3_MatchingAln[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]],
                                                                 UAC_3_MatchingAln[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] + 
                                                                   UAC_3_MatchingAln[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]]),
                               breaks = seq(from = -10, to = 10, by = 0.5)) +
            scale_x_reverse(limits = c(2750, -100), 
                            breaks = seq(from = 4000, to = -100, by = -100),
                            labels = c(4000, rep('', 4),
                                       3500,rep('', 4),
                                       3000, rep('', 4),
                                       2500,rep('', 4),
                                       2000, rep('', 4),
                                       1500,rep('', 4),
                                       1000, rep('', 4),
                                       500,rep('', 4),
                                       0, ''),
                            position = 'top',
                            sec.axis=sec_axis(~., 
                                              breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                              labels=c(2000, rep('', 4),
                                                       1500,rep('', 4),
                                                       1000, rep('', 4),
                                                       500,rep('', 4),
                                                       0, rep('', 4),
                                                       500,rep('', 4),
                                                       1000, rep('', 4),
                                                       1500,rep('', 4),
                                                       2000)
                            )
            ) +
            guides(y = "axis_truncated") +
            labs(y = 'Average UAC_3_s', x = 'Cal years BP', title = '') + 
            theme_classic() + 
            theme(plot.title = element_blank(), 
                  axis.text.y = element_text(size = 10, family = 'arial'),
                  axis.title.y = element_text(size = 12, family = 'arial'),
                  axis.text.x = element_text(size = 10, family = 'arial'),
                  axis.title.x = element_text(size = 12, family = 'arial')
            )
          
          UAC_3_AlnChronUnMean <- mean(UAC_3_MatchingAln.PA.ts.plot[["layers"]][[3]][["data"]][["y"]], 
                                       na.rm = TRUE)
          
          UAC_3_AlnChronUnSD <- sd(UAC_3_MatchingAln.PA.ts.plot[["layers"]][[3]][["data"]][["y"]], 
                                   na.rm = TRUE)
        }
      }
      
      #Calluna
      {
        #1
        {
          UAC_1_MatchingCal <- mapAgeEnsembleToPaleoData(UACLipdCalWithChron,age.var = "ageEnsemble",
                                                         paleo.depth.var = "Depth",
                                                         paleo.depth.range = c(202,278))
          
          UAC_1_MatchingCal.ae <- selectData(UAC_1_MatchingCal,var.name = "ageEnsemble")
          
          
          UAC_1_MatchingCal.PA <- selectData(UAC_1_MatchingCal,var.name = "Peak_Area")
          
          
          UAC_1_MatchingCal.PA.ts.plot <- plotTimeseriesEnsRibbons(X = UAC_1_MatchingCal.ae,
                                                                   Y = UAC_1_MatchingCal.PA,
                                                                   n.bins = 1000) +
            geom_point(
              aes(x = UAC_1_MatchingCal[["paleoData"]][[1]][["measurementTable"]][[1]][["ageMedian"]][["values"]],
                  y = UAC_1_MatchingCal[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]]), 
              size = 2, shape = 21, fill = 'red', alpha = 0.4) + 
            geom_errorbar(aes(x = UAC_1_MatchingCal[["paleoData"]][[1]][["measurementTable"]][[1]][["ageMedian"]][["values"]],
                              ymin = UAC_1_MatchingCal[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] - 
                                UAC_1_MatchingCal[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]],
                              ymax = UAC_1_MatchingCal[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] + 
                                UAC_1_MatchingCal[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]]),
                          size = 0.5, alpha = 0.4, color = 'black', width = 20) +
            scale_y_continuous(position = 'left', limits = range(UAC_1_MatchingCal.PA$values,
                                                                 UAC_1_MatchingCal[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] - 
                                                                   UAC_1_MatchingCal[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]],
                                                                 UAC_1_MatchingCal[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] + 
                                                                   UAC_1_MatchingCal[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]]),
                               breaks = seq(from = -10, to = 10, by = 0.5)) +
            scale_x_reverse(limits = c(2750, -100), 
                            breaks = seq(from = 4000, to = -100, by = -100),
                            labels = c(4000, rep('', 4),
                                       3500,rep('', 4),
                                       3000, rep('', 4),
                                       2500,rep('', 4),
                                       2000, rep('', 4),
                                       1500,rep('', 4),
                                       1000, rep('', 4),
                                       500,rep('', 4),
                                       0, ''),
                            position = 'top',
                            sec.axis=sec_axis(~., 
                                              breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                              labels=c(2000, rep('', 4),
                                                       1500,rep('', 4),
                                                       1000, rep('', 4),
                                                       500,rep('', 4),
                                                       0, rep('', 4),
                                                       500,rep('', 4),
                                                       1000, rep('', 4),
                                                       1500,rep('', 4),
                                                       2000)
                            )
            ) +
            guides(y = "axis_truncated") +
            labs(y = 'Average UAC_1_s', x = 'Cal years BP', title = '') + 
            theme_classic() + 
            theme(plot.title = element_blank(), 
                  axis.text.y = element_text(size = 10, family = 'arial'),
                  axis.title.y = element_text(size = 12, family = 'arial'),
                  axis.text.x = element_text(size = 10, family = 'arial'),
                  axis.title.x = element_text(size = 12, family = 'arial')
            )
          
          UAC_1_CalChronUnMean <- mean(UAC_1_MatchingCal.PA.ts.plot[["layers"]][[3]][["data"]][["y"]], 
                                       na.rm = TRUE)
          
          UAC_1_CalChronUnSD <- sd(UAC_1_MatchingCal.PA.ts.plot[["layers"]][[3]][["data"]][["y"]], 
                                   na.rm = TRUE)
        }
        
        #2
        {
          UAC_2_MatchingCal <- mapAgeEnsembleToPaleoData(UACLipdCalWithChron,age.var = "ageEnsemble",
                                                         paleo.depth.var = "Depth",
                                                         paleo.depth.range = c(96,200))
          
          UAC_2_MatchingCal.ae <- selectData(UAC_2_MatchingCal,var.name = "ageEnsemble")
          
          
          UAC_2_MatchingCal.PA <- selectData(UAC_2_MatchingCal,var.name = "Peak_Area")
          
          
          UAC_2_MatchingCal.PA.ts.plot <- plotTimeseriesEnsRibbons(X = UAC_2_MatchingCal.ae,
                                                                   Y = UAC_2_MatchingCal.PA,
                                                                   n.bins = 1000) +
            geom_point(
              aes(x = UAC_2_MatchingCal[["paleoData"]][[1]][["measurementTable"]][[1]][["ageMedian"]][["values"]],
                  y = UAC_2_MatchingCal[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]]), 
              size = 2, shape = 21, fill = 'red', alpha = 0.4) + 
            geom_errorbar(aes(x = UAC_2_MatchingCal[["paleoData"]][[1]][["measurementTable"]][[1]][["ageMedian"]][["values"]],
                              ymin = UAC_2_MatchingCal[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] - 
                                UAC_2_MatchingCal[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]],
                              ymax = UAC_2_MatchingCal[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] + 
                                UAC_2_MatchingCal[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]]),
                          size = 0.5, alpha = 0.4, color = 'black', width = 20) +
            scale_y_continuous(position = 'left', limits = range(UAC_2_MatchingCal.PA$values,
                                                                 UAC_2_MatchingCal[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] - 
                                                                   UAC_2_MatchingCal[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]],
                                                                 UAC_2_MatchingCal[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] + 
                                                                   UAC_2_MatchingCal[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]]),
                               breaks = seq(from = -10, to = 10, by = 0.5)) +
            scale_x_reverse(limits = c(2750, -100), 
                            breaks = seq(from = 4000, to = -100, by = -100),
                            labels = c(4000, rep('', 4),
                                       3500,rep('', 4),
                                       3000, rep('', 4),
                                       2500,rep('', 4),
                                       2000, rep('', 4),
                                       1500,rep('', 4),
                                       1000, rep('', 4),
                                       500,rep('', 4),
                                       0, ''),
                            position = 'top',
                            sec.axis=sec_axis(~., 
                                              breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                              labels=c(2000, rep('', 4),
                                                       1500,rep('', 4),
                                                       1000, rep('', 4),
                                                       500,rep('', 4),
                                                       0, rep('', 4),
                                                       500,rep('', 4),
                                                       1000, rep('', 4),
                                                       1500,rep('', 4),
                                                       2000)
                            )
            ) +
            guides(y = "axis_truncated") +
            labs(y = 'Average UAC_2_s', x = 'Cal years BP', title = '') + 
            theme_classic() + 
            theme(plot.title = element_blank(), 
                  axis.text.y = element_text(size = 10, family = 'arial'),
                  axis.title.y = element_text(size = 12, family = 'arial'),
                  axis.text.x = element_text(size = 10, family = 'arial'),
                  axis.title.x = element_text(size = 12, family = 'arial')
            )
          
          UAC_2_CalChronUnMean <- mean(UAC_2_MatchingCal.PA.ts.plot[["layers"]][[3]][["data"]][["y"]], 
                                       na.rm = TRUE)
          
          UAC_2_CalChronUnSD <- sd(UAC_2_MatchingCal.PA.ts.plot[["layers"]][[3]][["data"]][["y"]], 
                                   na.rm = TRUE)
        }
        
        #3
        {
          UAC_3_MatchingCal <- mapAgeEnsembleToPaleoData(UACLipdCalWithChron,age.var = "ageEnsemble",
                                                         paleo.depth.var = "Depth",
                                                         paleo.depth.range = c(2,94))
          
          UAC_3_MatchingCal.ae <- selectData(UAC_3_MatchingCal,var.name = "ageEnsemble")
          
          
          UAC_3_MatchingCal.PA <- selectData(UAC_3_MatchingCal,var.name = "Peak_Area")
          
          
          UAC_3_MatchingCal.PA.ts.plot <- plotTimeseriesEnsRibbons(X = UAC_3_MatchingCal.ae,
                                                                   Y = UAC_3_MatchingCal.PA,
                                                                   n.bins = 1000) +
            geom_point(
              aes(x = UAC_3_MatchingCal[["paleoData"]][[1]][["measurementTable"]][[1]][["ageMedian"]][["values"]],
                  y = UAC_3_MatchingCal[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]]), 
              size = 2, shape = 21, fill = 'red', alpha = 0.4) + 
            geom_errorbar(aes(x = UAC_3_MatchingCal[["paleoData"]][[1]][["measurementTable"]][[1]][["ageMedian"]][["values"]],
                              ymin = UAC_3_MatchingCal[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] - 
                                UAC_3_MatchingCal[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]],
                              ymax = UAC_3_MatchingCal[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] + 
                                UAC_3_MatchingCal[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]]),
                          size = 0.5, alpha = 0.4, color = 'black', width = 20) +
            scale_y_continuous(position = 'left', limits = range(UAC_3_MatchingCal.PA$values,
                                                                 UAC_3_MatchingCal[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] - 
                                                                   UAC_3_MatchingCal[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]],
                                                                 UAC_3_MatchingCal[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area"]][["values"]] + 
                                                                   UAC_3_MatchingCal[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]]),
                               breaks = seq(from = -10, to = 10, by = 0.5)) +
            scale_x_reverse(limits = c(2750, -100), 
                            breaks = seq(from = 4000, to = -100, by = -100),
                            labels = c(4000, rep('', 4),
                                       3500,rep('', 4),
                                       3000, rep('', 4),
                                       2500,rep('', 4),
                                       2000, rep('', 4),
                                       1500,rep('', 4),
                                       1000, rep('', 4),
                                       500,rep('', 4),
                                       0, ''),
                            position = 'top',
                            sec.axis=sec_axis(~., 
                                              breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                              labels=c(2000, rep('', 4),
                                                       1500,rep('', 4),
                                                       1000, rep('', 4),
                                                       500,rep('', 4),
                                                       0, rep('', 4),
                                                       500,rep('', 4),
                                                       1000, rep('', 4),
                                                       1500,rep('', 4),
                                                       2000)
                            )
            ) +
            guides(y = "axis_truncated") +
            labs(y = 'Average UAC_3_s', x = 'Cal years BP', title = '') + 
            theme_classic() + 
            theme(plot.title = element_blank(), 
                  axis.text.y = element_text(size = 10, family = 'arial'),
                  axis.title.y = element_text(size = 12, family = 'arial'),
                  axis.text.x = element_text(size = 10, family = 'arial'),
                  axis.title.x = element_text(size = 12, family = 'arial')
            )
          
          UAC_3_CalChronUnMean <- mean(UAC_3_MatchingCal.PA.ts.plot[["layers"]][[3]][["data"]][["y"]], 
                                       na.rm = TRUE)
          
          UAC_3_CalChronUnSD <- sd(UAC_3_MatchingCal.PA.ts.plot[["layers"]][[3]][["data"]][["y"]], 
                                   na.rm = TRUE)
        }
      }
    }
  }
  
  #Correlation
  {
    #Overall correlation
    {
      #Correlation analysis-TSI
      {
        #Read TSI data
        {
          TSILipd <- lipdR::readLipd(here('Archives_QSR', 'TSI_TW.lpd'))
          
          TSILipd.TSI <- selectData(TSILipd,var.name = "TSI")
          
          TSILipd.TSI.Time <- selectData(TSILipd,var.name = "year")
          
          TSILipd.TSI.Time <- convertAD2BP(TSILipd.TSI.Time)
          
        }
        
        #Correlation analysis
        {
          #Sphagnum
          {
            coroutSphTSI <- corEns(time.1 = UACMatchingSph.ae,
                                values.1 = UACMatchingSph.PA,
                                time.2 = t(matrix(replicate(1,TSILipd.TSI.Time$values),nrow = 1)),
                                values.2 = TSILipd.TSI,
                                bin.step = 40,
                                max.ens = 1000000,
                                min.obs = 5,
                                isopersistent  = TRUE,
                                isospectral = TRUE)
            
            coroutSph <- corEns(UACMatchingSph.ae,UACMatchingSph.PA,UACMatchingSph.ae,UACMatchingSph.PA,bin.step = 40,max.ens = 1000)
            corPlotSph <- plotCorEns(coroutSph,
                                     legend.position = c(0.1, 0.8),
                                     significance.option = "isospectral")+ggtitle(NULL)
            
            above90 <- round(sum(coroutSph$cor.df$r > 0.9)/nrow(coroutSph$cor.df)*100,1)
            med <- round(median(coroutSph$cor.df$r),2)
            
            plotCorEnsSphTSI <- 
              plotCorEns(coroutSphTSI,
                         bins = 20,
                         legend.position =c(.85,.8),
                         f.sig.lab.position = c(.85,.6),
                         significance.option = "isospectral",
                         use.fdr = TRUE)+ggtitle("Sphagnum - TSI")
          }
          
          #Alnus
          {
            coroutAlnTSI <- corEns(time.1 = UACMatchingAln.ae,
                                   values.1 = UACMatchingAln.PA,
                                   time.2 = t(matrix(replicate(1,TSILipd.TSI.Time$values),nrow = 1)),
                                   values.2 = TSILipd.TSI,
                                   bin.step = 40,
                                   max.ens = 1000000,
                                   min.obs = 5,
                                   isopersistent  = TRUE,
                                   isospectral = TRUE)
            
            coroutAln <- corEns(UACMatchingAln.ae,UACMatchingAln.PA,UACMatchingAln.ae,UACMatchingAln.PA,bin.step = 40,max.ens = 1000)
            corPlotAln <- plotCorEns(coroutAln,
                                     legend.position = c(0.1, 0.8),
                                     significance.option = "isospectral")+ggtitle(NULL)
            
            above90 <- round(sum(coroutAln$cor.df$r > 0.9)/nrow(coroutAln$cor.df)*100,1)
            med <- round(median(coroutAln$cor.df$r),2)
            
            plotCorEnsAlnTSI <- 
              plotCorEns(coroutAlnTSI,
                         bins = 20,
                         legend.position =c(.85,.8),
                         f.sig.lab.position = c(.85,.6),
                         significance.option = "isospectral",
                         use.fdr = TRUE)+ggtitle("Alnus - TSI")
          }
          
          #Calluna
          {
            coroutCalTSI <- corEns(time.1 = UACMatchingCal.ae,
                                   values.1 = UACMatchingCal.PA,
                                   time.2 = t(matrix(replicate(1,TSILipd.TSI.Time$values),nrow = 1)),
                                   values.2 = TSILipd.TSI,
                                   bin.step = 40,
                                   max.ens = 1000000,
                                   min.obs = 5,
                                   isopersistent  = TRUE,
                                   isospectral = TRUE)
            
            coroutCal <- corEns(UACMatchingCal.ae,UACMatchingCal.PA,UACMatchingCal.ae,UACMatchingCal.PA,bin.step = 40,max.ens = 1000)
            corPlotCal <- plotCorEns(coroutCal,
                                     legend.position = c(0.1, 0.8),
                                     significance.option = "isospectral")+ggtitle(NULL)
            
            above90 <- round(sum(coroutCal$cor.df$r > 0.9)/nrow(coroutCal$cor.df)*100,1)
            med <- round(median(coroutCal$cor.df$r),2)
            
            plotCorEnsCalTSI <- 
              plotCorEns(coroutCalTSI,
                         bins = 20,
                         legend.position =c(.85,.8),
                         f.sig.lab.position = c(.85,.6),
                         significance.option = "isospectral",
                         use.fdr = TRUE)+ggtitle("Calluna - TSI")
            

          }
          
          #Combine plot
          plotCorEnsSphTSI | plotCorEnsAlnTSI | plotCorEnsCalTSI
          
          #Subgroup
          {
            #Sphagnum
            {
              #TSI
              {
                coroutTSI_1 <- corEns(time.1 = t(matrix(replicate(1,TSILipd.TSI.Time$values),nrow = 1)),
                                      values.1 = TSILipd.TSI,
                                      time.2 = UAC_1_MatchingSph.ae,
                                      values.2 = UAC_1_MatchingSph.PA,
                                      bin.step = 40,
                                      max.ens = 1000000,
                                      isopersistent  = TRUE,
                                      isospectral = TRUE,
                                      gaussianize = FALSE)
                
                plotCorEnsSphTSI_1 <- 
                  plotCorEns(coroutTSI_1,
                             bins = 20,
                             legend.position =c(.85,.8),
                             f.sig.lab.position = c(.85,.6),
                             significance.option = "isospectral",
                             use.fdr = TRUE)+ggtitle("UAC Range 1 - TSI")
                
                coroutSph_1 <- corEns(UAC_1_MatchingSph.ae,UAC_1_MatchingSph.PA,UAC_1_MatchingSph.ae,UAC_1_MatchingSph.PA,bin.step = 40,max.ens = 1000)
                corPlotSph_1 <- plotCorEns(coroutSph_1,
                                           legend.position = c(0.1, 0.8),
                                           significance.option = "isospectral")+ggtitle(NULL)
                
                above90_1 <- round(sum(coroutSph$cor.df$r > 0.9)/nrow(coroutSph$cor.df)*100,1)
                med_1 <- round(median(coroutSph$cor.df$r),2)
                
                coroutTSI_2 <- corEns(time.1 = t(matrix(replicate(1,TSILipd.TSI.Time$values),nrow = 1)),
                                      values.1 = TSILipd.TSI,
                                      time.2 = UAC_2_MatchingSph.ae,
                                      values.2 = UAC_2_MatchingSph.PA,
                                      bin.step = 40,
                                      max.ens = 1000000,
                                      isopersistent  = TRUE,
                                      isospectral = TRUE,
                                      gaussianize = FALSE)
                
                plotCorEnsSphTSI_2 <- 
                  plotCorEns(coroutTSI_2,
                             bins = 20,
                             legend.position =c(.85,.8),
                             f.sig.lab.position = c(.85,.6),
                             significance.option = "isospectral",
                             use.fdr = TRUE)+ggtitle("UAC Range 2 - TSI")
                
                coroutSph_2 <- corEns(UAC_2_MatchingSph.ae,UAC_2_MatchingSph.PA,UAC_2_MatchingSph.ae,UAC_2_MatchingSph.PA,bin.step = 40,max.ens = 1000)
                corPlotSph_2 <- plotCorEns(coroutSph_2,
                                           legend.position = c(0.1, 0.8),
                                           significance.option = "isospectral")+ggtitle(NULL)
                
                coroutTSI_3 <- corEns(time.1 = t(matrix(replicate(1,TSILipd.TSI.Time$values),nrow = 1)),
                                      values.1 = TSILipd.TSI,
                                      time.2 = UAC_3_MatchingSph.ae,
                                      values.2 = UAC_3_MatchingSph.PA,
                                      bin.step = 40,
                                      max.ens = 1000000,
                                      isopersistent  = TRUE,
                                      isospectral = TRUE,
                                      gaussianize = FALSE)
                
                plotCorEnsSphTSI_3 <- 
                  plotCorEns(coroutTSI_3,
                             bins = 20,
                             legend.position =c(.85,.8),
                             f.sig.lab.position = c(.85,.6),
                             significance.option = "isospectral",
                             use.fdr = TRUE)+ggtitle("UAC Range 3 - TSI")
                
                coroutSph_3 <- corEns(UAC_3_MatchingSph.ae,UAC_3_MatchingSph.PA,UAC_3_MatchingSph.ae,UAC_3_MatchingSph.PA,bin.step = 40,max.ens = 1000)
                corPlotSph_3 <- plotCorEns(coroutSph_3,
                                           legend.position = c(0.1, 0.8),
                                           significance.option = "isospectral")+ggtitle(NULL)
                
                #Combine plots
                plotCorEnsSphTSI_1 | plotCorEnsSphTSI_2 | plotCorEnsSphTSI_3
              }
            }
            
            #Alnus
            {
              #TSI
              {
                coroutTSI_1 <- corEns(time.1 = t(matrix(replicate(1,TSILipd.TSI.Time$values),nrow = 1)),
                                      values.1 = TSILipd.TSI,
                                      time.2 = UAC_1_MatchingAln.ae,
                                      values.2 = UAC_1_MatchingAln.PA,
                                      bin.step = 40,
                                      max.ens = 1000000,
                                      isopersistent  = TRUE,
                                      isospectral = TRUE,
                                      gaussianize = FALSE)
                
                plotCorEnsAlnTSI_1 <- 
                  plotCorEns(coroutTSI_1,
                             bins = 20,
                             legend.position =c(.85,.8),
                             f.sig.lab.position = c(.85,.6),
                             significance.option = "isospectral",
                             use.fdr = TRUE)+ggtitle("UAC Range 1 - TSI")
                
                coroutAln_1 <- corEns(UAC_1_MatchingAln.ae,UAC_1_MatchingAln.PA,UAC_1_MatchingAln.ae,UAC_1_MatchingAln.PA,bin.step = 40,max.ens = 1000)
                corPlotAln_1 <- plotCorEns(coroutAln_1,
                                           legend.position = c(0.1, 0.8),
                                           significance.option = "isospectral")+ggtitle(NULL)
                
                above90_1 <- round(sum(coroutAln$cor.df$r > 0.9)/nrow(coroutAln$cor.df)*100,1)
                med_1 <- round(median(coroutAln$cor.df$r),2)
                
                coroutTSI_2 <- corEns(time.1 = t(matrix(replicate(1,TSILipd.TSI.Time$values),nrow = 1)),
                                      values.1 = TSILipd.TSI,
                                      time.2 = UAC_2_MatchingAln.ae,
                                      values.2 = UAC_2_MatchingAln.PA,
                                      bin.step = 40,
                                      max.ens = 1000000,
                                      isopersistent  = TRUE,
                                      isospectral = TRUE,
                                      gaussianize = FALSE)
                
                plotCorEnsAlnTSI_2 <- 
                  plotCorEns(coroutTSI_2,
                             bins = 20,
                             legend.position =c(.85,.8),
                             f.sig.lab.position = c(.85,.6),
                             significance.option = "isospectral",
                             use.fdr = TRUE)+ggtitle("UAC Range 2 - TSI")
                
                coroutAln_2 <- corEns(UAC_2_MatchingAln.ae,UAC_2_MatchingAln.PA,UAC_2_MatchingAln.ae,UAC_2_MatchingAln.PA,bin.step = 40,max.ens = 1000)
                corPlotAln_2 <- plotCorEns(coroutAln_2,
                                           legend.position = c(0.1, 0.8),
                                           significance.option = "isospectral")+ggtitle(NULL)
                
                coroutTSI_3 <- corEns(time.1 = t(matrix(replicate(1,TSILipd.TSI.Time$values),nrow = 1)),
                                      values.1 = TSILipd.TSI,
                                      time.2 = UAC_3_MatchingAln.ae,
                                      values.2 = UAC_3_MatchingAln.PA,
                                      bin.step = 40,
                                      max.ens = 1000000,
                                      isopersistent  = TRUE,
                                      isospectral = TRUE,
                                      gaussianize = FALSE)
                
                plotCorEnsAlnTSI_3 <- 
                  plotCorEns(coroutTSI_3,
                             bins = 20,
                             legend.position =c(.85,.8),
                             f.sig.lab.position = c(.85,.6),
                             significance.option = "isospectral",
                             use.fdr = TRUE)+ggtitle("UAC Range 3 - TSI")
                
                coroutAln_3 <- corEns(UAC_3_MatchingAln.ae,UAC_3_MatchingAln.PA,UAC_3_MatchingAln.ae,UAC_3_MatchingAln.PA,bin.step = 40,max.ens = 1000)
                corPlotAln_3 <- plotCorEns(coroutAln_3,
                                           legend.position = c(0.1, 0.8),
                                           significance.option = "isospectral")+ggtitle(NULL)
                
                #Combine plots
                plotCorEnsAlnTSI_1 | plotCorEnsAlnTSI_2 | plotCorEnsAlnTSI_3
              }
            }
            
            #Calluna
            {
              #TSI
              {
                coroutTSI_1 <- corEns(time.1 = t(matrix(replicate(1,TSILipd.TSI.Time$values),nrow = 1)),
                                      values.1 = TSILipd.TSI,
                                      time.2 = UAC_1_MatchingCal.ae,
                                      values.2 = UAC_1_MatchingCal.PA,
                                      bin.step = 40,
                                      max.ens = 1000000,
                                      isopersistent  = TRUE,
                                      isospectral = TRUE,
                                      gaussianize = FALSE)
                
                plotCorEnsCalTSI_1 <- 
                  plotCorEns(coroutTSI_1,
                             bins = 20,
                             legend.position =c(.85,.8),
                             f.sig.lab.position = c(.85,.6),
                             significance.option = "isospectral",
                             use.fdr = TRUE)+ggtitle("UAC Range 1 - TSI")
                
                coroutCal_1 <- corEns(UAC_1_MatchingCal.ae,UAC_1_MatchingCal.PA,UAC_1_MatchingCal.ae,UAC_1_MatchingCal.PA,bin.step = 40,max.ens = 1000)
                corPlotCal_1 <- plotCorEns(coroutCal_1,
                                           legend.position = c(0.1, 0.8),
                                           significance.option = "isospectral")+ggtitle(NULL)
                
                above90_1 <- round(sum(coroutCal$cor.df$r > 0.9)/nrow(coroutCal$cor.df)*100,1)
                med_1 <- round(median(coroutCal$cor.df$r),2)
                
                coroutTSI_2 <- corEns(time.1 = t(matrix(replicate(1,TSILipd.TSI.Time$values),nrow = 1)),
                                      values.1 = TSILipd.TSI,
                                      time.2 = UAC_2_MatchingCal.ae,
                                      values.2 = UAC_2_MatchingCal.PA,
                                      bin.step = 40,
                                      max.ens = 1000000,
                                      isopersistent  = TRUE,
                                      isospectral = TRUE,
                                      gaussianize = FALSE)
                
                plotCorEnsCalTSI_2 <- 
                  plotCorEns(coroutTSI_2,
                             bins = 20,
                             legend.position =c(.85,.8),
                             f.sig.lab.position = c(.85,.6),
                             significance.option = "isospectral",
                             use.fdr = TRUE)+ggtitle("UAC Range 2 - TSI")
                
                coroutCal_2 <- corEns(UAC_2_MatchingCal.ae,UAC_2_MatchingCal.PA,UAC_2_MatchingCal.ae,UAC_2_MatchingCal.PA,bin.step = 40,max.ens = 1000)
                corPlotCal_2 <- plotCorEns(coroutCal_2,
                                           legend.position = c(0.1, 0.8),
                                           significance.option = "isospectral")+ggtitle(NULL)
                
                coroutTSI_3 <- corEns(time.1 = t(matrix(replicate(1,TSILipd.TSI.Time$values),nrow = 1)),
                                      values.1 = TSILipd.TSI,
                                      time.2 = UAC_3_MatchingCal.ae,
                                      values.2 = UAC_3_MatchingCal.PA,
                                      bin.step = 40,
                                      max.ens = 1000000,
                                      isopersistent  = TRUE,
                                      isospectral = TRUE,
                                      gaussianize = FALSE)
                
                plotCorEnsCalTSI_3 <- 
                  plotCorEns(coroutTSI_3,
                             bins = 20,
                             legend.position =c(.85,.8),
                             f.sig.lab.position = c(.85,.6),
                             significance.option = "isospectral",
                             use.fdr = TRUE)+ggtitle("UAC Range 3 - TSI")
                
                coroutCal_3 <- corEns(UAC_3_MatchingCal.ae,UAC_3_MatchingCal.PA,UAC_3_MatchingCal.ae,UAC_3_MatchingCal.PA,bin.step = 40,max.ens = 1000)
                corPlotCal_3 <- plotCorEns(coroutCal_3,
                                           legend.position = c(0.1, 0.8),
                                           significance.option = "isospectral")+ggtitle(NULL)
                
                #Combine plots
                plotCorEnsCalTSI_1 | plotCorEnsCalTSI_2 | plotCorEnsCalTSI_3
              }
            }
          }
        }
      }
      
      #Correlation analysis-CC
      {
        #Read CC data
        {
          CCLipd <- lipdR::readLipd(here('Archives_QSR', 'CC_TW.lpd'))
          
          CCLipd.CC <- selectData(CCLipd,var.name = "CC")
          
          CCLipd.CC.Time <- selectData(CCLipd,var.name = "year")
          
          CCLipd.CC.Time <- convertAD2BP(CCLipd.CC.Time)
          
        }
        
        #Correlation analysis
        {
          #Sphagnum
          {
            coroutSphCC <- corEns(time.1 = UACMatchingSph.ae,
                                   values.1 = UACMatchingSph.PA,
                                   time.2 = t(matrix(replicate(1,CCLipd.CC.Time$values),nrow = 1)),
                                   values.2 = CCLipd.CC,
                                   bin.step = 40,
                                   max.ens = 1000000,
                                   min.obs = 5,
                                   isopersistent  = TRUE,
                                   isospectral = TRUE)
            
            plotCorEnsSphCC <- 
              plotCorEns(coroutSphCC,
                         bins = 20,
                         legend.position =c(.85,.8),
                         f.sig.lab.position = c(.85,.6),
                         significance.option = "isospectral",
                         use.fdr = TRUE)+ggtitle("Sphagnum - CC")
          }
          
          #Alnus
          {
            coroutAlnCC <- corEns(time.1 = UACMatchingAln.ae,
                                   values.1 = UACMatchingAln.PA,
                                   time.2 = t(matrix(replicate(1,CCLipd.CC.Time$values),nrow = 1)),
                                   values.2 = CCLipd.CC,
                                   bin.step = 40,
                                   max.ens = 1000000,
                                   min.obs = 5,
                                   isopersistent  = TRUE,
                                   isospectral = TRUE)
            
            plotCorEnsAlnCC <- 
              plotCorEns(coroutAlnCC,
                         bins = 20,
                         legend.position =c(.85,.8),
                         f.sig.lab.position = c(.85,.6),
                         significance.option = "isospectral",
                         use.fdr = TRUE)+ggtitle("Alnus - CC")
          }
          
          #Calluna
          {
            coroutCalCC <- corEns(time.1 = UACMatchingCal.ae,
                                   values.1 = UACMatchingCal.PA,
                                   time.2 = t(matrix(replicate(1,CCLipd.CC.Time$values),nrow = 1)),
                                   values.2 = CCLipd.CC,
                                   bin.step = 40,
                                   max.ens = 1000000,
                                   min.obs = 5,
                                   isopersistent  = TRUE,
                                   isospectral = TRUE)
            
            plotCorEnsCalCC <- 
              plotCorEns(coroutCalCC,
                         bins = 20,
                         legend.position =c(.85,.8),
                         f.sig.lab.position = c(.85,.6),
                         significance.option = "isospectral",
                         use.fdr = TRUE)+ggtitle("Calluna - CC")
          }
          
          #Combine plots
          plotCorEnsSphCC | plotCorEnsAlnCC | plotCorEnsCalCC
        }
      }
      
      #Correlation analysis-PC-TSI
      {
        #Read PC data
        {
          #PC1
          {
            PC1Lipd <- lipdR::readLipd(here('Archives_QSR', 'PC1.lpd'))
            
            UACLipdPC1WithChron <- runBacon(PC1Lipd, lab.id.var = 'AMS_Number',
                                            bacon.dir = here(),
                                            age.14c.uncertainty.var = 'age14Cuncertainty', 
                                            age.var = 'Calibrated_Ages', 
                                            age.uncertainty.var = 'Calibrated_Ages_Var',
                                            reservoir.age.14c.var = NULL, 
                                            reservoir.age.14c.uncertainty.var = NULL, 
                                            rejected.ages.var = NULL,
                                            bacon.acc.mean = 20,
                                            cc = 1,
                                            bacon.thick = 10,
                                            suggest = FALSE,
                                            max.ens = 1000,
                                            accept.suggestions = TRUE)
            
            UACMatchingPC1 <- mapAgeEnsembleToPaleoData(UACLipdPC1WithChron,age.var = "ageEnsemble",
                                                        paleo.depth.var = "Depth",)
            
            UACMatchingPC1.ae <- selectData(UACMatchingPC1,var.name = "ageEnsemble")
            
            
            UACMatchingPC1.PA <- selectData(UACMatchingPC1,var.name = "Peak_Area")
            
            
            UACMatchingPC1.PA.ts.plot <- plotTimeseriesEnsRibbons(X = UACMatchingPC1.ae,
                                                                  Y = UACMatchingPC1.PA,
                                                                  n.bins = 1000) +
              scale_y_continuous(position = 'left', limits = range(UACMatchingPC1.PA$values)) +
              scale_x_reverse(limits = c(2750, -100), 
                              breaks = seq(from = 4000, to = -100, by = -100),
                              labels = c(4000, rep('', 4),
                                         3500,rep('', 4),
                                         3000, rep('', 4),
                                         2500,rep('', 4),
                                         2000, rep('', 4),
                                         1500,rep('', 4),
                                         1000, rep('', 4),
                                         500,rep('', 4),
                                         0, ''),
                              position = 'top',
                              sec.axis=sec_axis(~., 
                                                breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                                labels=c(2000, rep('', 4),
                                                         1500,rep('', 4),
                                                         1000, rep('', 4),
                                                         500,rep('', 4),
                                                         0, rep('', 4),
                                                         500,rep('', 4),
                                                         1000, rep('', 4),
                                                         1500,rep('', 4),
                                                         2000)
                              )
              ) +
              guides(y = "axis_truncated") +
              labs(y = 'Average UACs', x = 'Cal years BP', title = '') + 
              theme_classic() + 
              theme(plot.title = element_blank(), 
                    axis.text.y = element_text(size = 10, family = 'arial'),
                    axis.title.y = element_text(size = 12, family = 'arial'),
                    axis.text.x = element_text(size = 10, family = 'arial'),
                    axis.title.x = element_text(size = 12, family = 'arial')
              )
            
            UACPC1ChronUnMean <- mean(UACMatchingPC1.PA.ts.plot[["layers"]][[3]][["data"]][["y"]], 
                                      na.rm = TRUE)
            
            UACPC1ChronUnSD <- sd(UACMatchingPC1.PA.ts.plot[["layers"]][[3]][["data"]][["y"]], 
                                  na.rm = TRUE)
          }
          
          #PC2
          {
            PC2Lipd <- lipdR::readLipd(here('Archives_QSR', 'PC2.lpd'))
            
            UACLipdPC2WithChron <- runBacon(PC2Lipd, lab.id.var = 'AMS_Number',
                                            bacon.dir = here(),
                                            age.14c.uncertainty.var = 'age14Cuncertainty', 
                                            age.var = 'Calibrated_Ages', 
                                            age.uncertainty.var = 'Calibrated_Ages_Var',
                                            reservoir.age.14c.var = NULL, 
                                            reservoir.age.14c.uncertainty.var = NULL, 
                                            rejected.ages.var = NULL,
                                            bacon.acc.mean = 20,
                                            cc = 1,
                                            bacon.thick = 10,
                                            suggest = FALSE,
                                            max.ens = 1000,
                                            accept.suggestions = TRUE)
            
            UACMatchingPC2 <- mapAgeEnsembleToPaleoData(UACLipdPC2WithChron,age.var = "ageEnsemble",
                                                        paleo.depth.var = "Depth",)
            
            UACMatchingPC2.ae <- selectData(UACMatchingPC2,var.name = "ageEnsemble")
            
            
            UACMatchingPC2.PA <- selectData(UACMatchingPC2,var.name = "Peak_Area")
            
            
            UACMatchingPC2.PA.ts.plot <- plotTimeseriesEnsRibbons(X = UACMatchingPC2.ae,
                                                                  Y = UACMatchingPC2.PA,
                                                                  n.bins = 1000) +
              scale_y_continuous(position = 'left', limits = range(UACMatchingPC2.PA$values)) +
              scale_x_reverse(limits = c(2750, -100), 
                              breaks = seq(from = 4000, to = -100, by = -100),
                              labels = c(4000, rep('', 4),
                                         3500,rep('', 4),
                                         3000, rep('', 4),
                                         2500,rep('', 4),
                                         2000, rep('', 4),
                                         1500,rep('', 4),
                                         1000, rep('', 4),
                                         500,rep('', 4),
                                         0, ''),
                              position = 'top',
                              sec.axis=sec_axis(~., 
                                                breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                                labels=c(2000, rep('', 4),
                                                         1500,rep('', 4),
                                                         1000, rep('', 4),
                                                         500,rep('', 4),
                                                         0, rep('', 4),
                                                         500,rep('', 4),
                                                         1000, rep('', 4),
                                                         1500,rep('', 4),
                                                         2000)
                              )
              ) +
              guides(y = "axis_truncated") +
              labs(y = 'Average UACs', x = 'Cal years BP', title = '') + 
              theme_classic() + 
              theme(plot.title = element_blank(), 
                    axis.text.y = element_text(size = 10, family = 'arial'),
                    axis.title.y = element_text(size = 12, family = 'arial'),
                    axis.text.x = element_text(size = 10, family = 'arial'),
                    axis.title.x = element_text(size = 12, family = 'arial')
              )
            
            UACPC2ChronUnMean <- mean(UACMatchingPC2.PA.ts.plot[["layers"]][[3]][["data"]][["y"]], 
                                      na.rm = TRUE)
            
            UACPC2ChronUnSD <- sd(UACMatchingPC2.PA.ts.plot[["layers"]][[3]][["data"]][["y"]], 
                                  na.rm = TRUE)
          }
          
          #PC3
          {
            PC3Lipd <- lipdR::readLipd(here('Archives_QSR', 'PC3.lpd'))
            
            UACLipdPC3WithChron <- runBacon(PC3Lipd, lab.id.var = 'AMS_Number',
                                            bacon.dir = here(),
                                            age.14c.uncertainty.var = 'age14Cuncertainty', 
                                            age.var = 'Calibrated_Ages', 
                                            age.uncertainty.var = 'Calibrated_Ages_Var',
                                            reservoir.age.14c.var = NULL, 
                                            reservoir.age.14c.uncertainty.var = NULL, 
                                            rejected.ages.var = NULL,
                                            bacon.acc.mean = 20,
                                            cc = 1,
                                            bacon.thick = 10,
                                            suggest = FALSE,
                                            max.ens = 1000,
                                            accept.suggestions = TRUE)
            
            UACMatchingPC3 <- mapAgeEnsembleToPaleoData(UACLipdPC3WithChron,age.var = "ageEnsemble",
                                                        paleo.depth.var = "Depth",)
            
            UACMatchingPC3.ae <- selectData(UACMatchingPC3,var.name = "ageEnsemble")
            
            
            UACMatchingPC3.PA <- selectData(UACMatchingPC3,var.name = "Peak_Area")
            
            
            UACMatchingPC3.PA.ts.plot <- plotTimeseriesEnsRibbons(X = UACMatchingPC3.ae,
                                                                  Y = UACMatchingPC3.PA,
                                                                  n.bins = 1000) +
              scale_y_continuous(position = 'left', limits = range(UACMatchingPC3.PA$values)) +
              scale_x_reverse(limits = c(2750, -100), 
                              breaks = seq(from = 4000, to = -100, by = -100),
                              labels = c(4000, rep('', 4),
                                         3500,rep('', 4),
                                         3000, rep('', 4),
                                         2500,rep('', 4),
                                         2000, rep('', 4),
                                         1500,rep('', 4),
                                         1000, rep('', 4),
                                         500,rep('', 4),
                                         0, ''),
                              position = 'top',
                              sec.axis=sec_axis(~., 
                                                breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                                labels=c(2000, rep('', 4),
                                                         1500,rep('', 4),
                                                         1000, rep('', 4),
                                                         500,rep('', 4),
                                                         0, rep('', 4),
                                                         500,rep('', 4),
                                                         1000, rep('', 4),
                                                         1500,rep('', 4),
                                                         2000)
                              )
              ) +
              guides(y = "axis_truncated") +
              labs(y = 'Average UACs', x = 'Cal years BP', title = '') + 
              theme_classic() + 
              theme(plot.title = element_blank(), 
                    axis.text.y = element_text(size = 10, family = 'arial'),
                    axis.title.y = element_text(size = 12, family = 'arial'),
                    axis.text.x = element_text(size = 10, family = 'arial'),
                    axis.title.x = element_text(size = 12, family = 'arial')
              )
            
            UACPC3ChronUnMean <- mean(UACMatchingPC3.PA.ts.plot[["layers"]][[3]][["data"]][["y"]], 
                                      na.rm = TRUE)
            
            UACPC3ChronUnSD <- sd(UACMatchingPC3.PA.ts.plot[["layers"]][[3]][["data"]][["y"]], 
                                  na.rm = TRUE)
          }
          
          #Combine plot
          {
            ggarrange(UACMatchingPC1.PA.ts.plot+
                        geom_hline(yintercept = UACPC1ChronUnMean, 
                                   linetype = 1, 
                                   color = 'blue', alpha = 0.5) +
                        geom_hline(yintercept = UACPC1ChronUnMean - UACPC1ChronUnSD, 
                                   linetype = 2, 
                                   color = 'blue', alpha = 0.5) +
                        geom_hline(yintercept = UACPC1ChronUnMean - 2 * UACPC1ChronUnSD, 
                                   linetype = 2, 
                                   color = 'blue', alpha = 0.5) +
                        geom_hline(yintercept = UACPC1ChronUnMean + UACPC1ChronUnSD, 
                                   linetype = 2, 
                                   color = 'blue', alpha = 0.5) +
                        geom_hline(yintercept = UACPC1ChronUnMean + 2* UACPC1ChronUnSD, 
                                   linetype = 2, 
                                   color = 'blue', alpha = 0.5) ,
                      UACMatchingPC2.PA.ts.plot+
                        geom_hline(yintercept = UACPC2ChronUnMean, 
                                   linetype = 1, 
                                   color = 'blue', alpha = 0.5) +
                        geom_hline(yintercept = UACPC2ChronUnMean - UACPC2ChronUnSD, 
                                   linetype = 2, 
                                   color = 'blue', alpha = 0.5) +
                        geom_hline(yintercept = UACPC2ChronUnMean - 2 * UACPC2ChronUnSD, 
                                   linetype = 2, 
                                   color = 'blue', alpha = 0.5) +
                        geom_hline(yintercept = UACPC2ChronUnMean + UACPC2ChronUnSD, 
                                   linetype = 2, 
                                   color = 'blue', alpha = 0.5) +
                        geom_hline(yintercept = UACPC2ChronUnMean + 2* UACPC2ChronUnSD, 
                                   linetype = 2, 
                                   color = 'blue', alpha = 0.5) +
                        theme(
                          axis.line.x = element_blank(),
                          axis.ticks.x = element_blank(),
                          axis.title.x = element_blank(),
                          axis.text.x = element_blank(),
                        ),
                      UACMatchingPC3.PA.ts.plot+
                        geom_hline(yintercept = UACPC3ChronUnMean, 
                                   linetype = 1, 
                                   color = 'blue', alpha = 0.5) +
                        geom_hline(yintercept = UACPC3ChronUnMean - UACPC3ChronUnSD, 
                                   linetype = 2, 
                                   color = 'blue', alpha = 0.5) +
                        geom_hline(yintercept = UACPC3ChronUnMean - 2 * UACPC3ChronUnSD, 
                                   linetype = 2, 
                                   color = 'blue', alpha = 0.5) +
                        geom_hline(yintercept = UACPC3ChronUnMean + UACPC3ChronUnSD, 
                                   linetype = 2, 
                                   color = 'blue', alpha = 0.5) +
                        geom_hline(yintercept = UACPC3ChronUnMean + 2* UACPC3ChronUnSD, 
                                   linetype = 2, 
                                   color = 'blue', alpha = 0.5) +
                        theme(
                          axis.line.x = element_blank(),
                          axis.ticks.x = element_blank(),
                          axis.title.x = element_blank(),
                          axis.text.x = element_blank(),
                        ),
                      WholeTSI,
                      heights = c(1.2,1,1,1.2),
                      ncol = 1,
                      nrow = 4,
                      align = 'v')
            
            ggarrange(UACMatchingPC1.PA.ts.plot +
                        geom_hline(yintercept = UACPC1ChronUnMean, 
                                   linetype = 1, 
                                   color = 'blue', alpha = 0.5) +
                        geom_hline(yintercept = UACPC1ChronUnMean - UACPC1ChronUnSD, 
                                   linetype = 2, 
                                   color = 'blue', alpha = 0.5) +
                        geom_hline(yintercept = UACPC1ChronUnMean - 2 * UACPC1ChronUnSD, 
                                   linetype = 2, 
                                   color = 'blue', alpha = 0.5) +
                        geom_hline(yintercept = UACPC1ChronUnMean + UACPC1ChronUnSD, 
                                   linetype = 2, 
                                   color = 'blue', alpha = 0.5) +
                        geom_hline(yintercept = UACPC1ChronUnMean + 2* UACPC1ChronUnSD, 
                                   linetype = 2, 
                                   color = 'blue', alpha = 0.5) +
                        scale_x_reverse(limits = c(1000, -100), 
                                        breaks = seq(from = 4000, to = -100, by = -100),
                                        labels = c(4000, rep('', 4),
                                                   3500,rep('', 4),
                                                   3000, rep('', 4),
                                                   2500,rep('', 4),
                                                   2000, rep('', 4),
                                                   1500,rep('', 4),
                                                   1000, rep('', 4),
                                                   500,rep('', 4),
                                                   0, ''),
                                        sec.axis=sec_axis(~., 
                                                          breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                                          labels=c(2000, rep('', 4),
                                                                   1500,rep('', 4),
                                                                   1000, rep('', 4),
                                                                   500,rep('', 4),
                                                                   0, rep('', 4),
                                                                   500,rep('', 4),
                                                                   1000, rep('', 4),
                                                                   1500,rep('', 4),
                                                                   2000)
                                        ) 
                        ),
                      UACMatchingPC2.PA.ts.plot +
                        geom_hline(yintercept = UACPC2ChronUnMean, 
                                   linetype = 1, 
                                   color = 'blue', alpha = 0.5) +
                        geom_hline(yintercept = UACPC2ChronUnMean - UACPC2ChronUnSD, 
                                   linetype = 2, 
                                   color = 'blue', alpha = 0.5) +
                        geom_hline(yintercept = UACPC2ChronUnMean - 2 * UACPC2ChronUnSD, 
                                   linetype = 2, 
                                   color = 'blue', alpha = 0.5) +
                        geom_hline(yintercept = UACPC2ChronUnMean + UACPC2ChronUnSD, 
                                   linetype = 2, 
                                   color = 'blue', alpha = 0.5) +
                        geom_hline(yintercept = UACPC2ChronUnMean + 2* UACPC2ChronUnSD, 
                                   linetype = 2, 
                                   color = 'blue', alpha = 0.5) +
                        scale_x_reverse(limits = c(1000, -100), 
                                        breaks = seq(from = 4000, to = -100, by = -100),
                                        labels = c(4000, rep('', 4),
                                                   3500,rep('', 4),
                                                   3000, rep('', 4),
                                                   2500,rep('', 4),
                                                   2000, rep('', 4),
                                                   1500,rep('', 4),
                                                   1000, rep('', 4),
                                                   500,rep('', 4),
                                                   0, ''),
                                        sec.axis=sec_axis(~., 
                                                          breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                                          labels=c(2000, rep('', 4),
                                                                   1500,rep('', 4),
                                                                   1000, rep('', 4),
                                                                   500,rep('', 4),
                                                                   0, rep('', 4),
                                                                   500,rep('', 4),
                                                                   1000, rep('', 4),
                                                                   1500,rep('', 4),
                                                                   2000)
                                        )
                        )+
                        theme(
                          axis.line.x = element_blank(),
                          axis.ticks.x = element_blank(),
                          axis.title.x = element_blank(),
                          axis.text.x = element_blank(),
                        ),
                      UACMatchingPC3.PA.ts.plot+
                        geom_hline(yintercept = UACPC3ChronUnMean, 
                                   linetype = 1, 
                                   color = 'blue', alpha = 0.5) +
                        geom_hline(yintercept = UACPC3ChronUnMean - UACPC3ChronUnSD, 
                                   linetype = 2, 
                                   color = 'blue', alpha = 0.5) +
                        geom_hline(yintercept = UACPC3ChronUnMean - 2 * UACPC3ChronUnSD, 
                                   linetype = 2, 
                                   color = 'blue', alpha = 0.5) +
                        geom_hline(yintercept = UACPC3ChronUnMean + UACPC3ChronUnSD, 
                                   linetype = 2, 
                                   color = 'blue', alpha = 0.5) +
                        geom_hline(yintercept = UACPC3ChronUnMean + 2* UACPC3ChronUnSD, 
                                   linetype = 2, 
                                   color = 'blue', alpha = 0.5) +
                        scale_x_reverse(limits = c(1000, -100), 
                                        breaks = seq(from = 4000, to = -100, by = -100),
                                        labels = c(4000, rep('', 4),
                                                   3500,rep('', 4),
                                                   3000, rep('', 4),
                                                   2500,rep('', 4),
                                                   2000, rep('', 4),
                                                   1500,rep('', 4),
                                                   1000, rep('', 4),
                                                   500,rep('', 4),
                                                   0, ''),
                                        sec.axis=sec_axis(~., 
                                                          breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                                          labels=c(2000, rep('', 4),
                                                                   1500,rep('', 4),
                                                                   1000, rep('', 4),
                                                                   500,rep('', 4),
                                                                   0, rep('', 4),
                                                                   500,rep('', 4),
                                                                   1000, rep('', 4),
                                                                   1500,rep('', 4),
                                                                   2000)
                                        )
                        )+
                        theme(
                          axis.line.x = element_blank(),
                          axis.ticks.x = element_blank(),
                          axis.title.x = element_blank(),
                          axis.text.x = element_blank(),
                        ),
                      WholeCC,
                      heights = c(1.2,1,1,1.2),
                      ncol = 1,
                      nrow = 4,
                      align = 'v')
          }
          
        }
        
        #Correlation analysis
        {
          #PC1
          {
            coroutPC1TSI <- corEns(time.1 = UACMatchingPC1.ae,
                                   values.1 = UACMatchingPC1.PA,
                                   time.2 = t(matrix(replicate(1,TSILipd.TSI.Time$values),nrow = 1)),
                                   values.2 = TSILipd.TSI,
                                   bin.step = 40,
                                   max.ens = 1000000,
                                   min.obs = 5,
                                   isopersistent  = TRUE,
                                   isospectral = TRUE)
            plotCorEnsPC1TSI <- 
              plotCorEns(coroutPC1TSI,
                         bins = 20,
                         legend.position =c(.85,.8),
                         f.sig.lab.position = c(.85,.6),
                         significance.option = "isospectral",
                         use.fdr = TRUE)+ggtitle("PC1 - TSI")
            
            coroutPC1 <- corEns(UACMatchingPC1.ae,UACMatchingPC1.PA,UACMatchingPC1.ae,UACMatchingPC1.PA,bin.step = 40,max.ens = 1000)
            corPlotPC1 <- plotCorEns(coroutPC1,
                                     legend.position = c(0.1, 0.8),
                                     significance.option = "isospectral")+ggtitle(NULL)
            
            above90 <- round(sum(coroutPC1$cor.df$r > 0.9)/nrow(coroutPC1$cor.df)*100,1)
            med <- round(median(coroutPC1$cor.df$r),2)
            
            
          }
          
          #PC2
          {
            coroutPC2TSI <- corEns(time.1 = UACMatchingPC2.ae,
                                   values.1 = UACMatchingPC2.PA,
                                   time.2 = t(matrix(replicate(1,TSILipd.TSI.Time$values),nrow = 1)),
                                   values.2 = TSILipd.TSI,
                                   bin.step = 40,
                                   max.ens = 1000000,
                                   min.obs = 5,
                                   isopersistent  = TRUE,
                                   isospectral = TRUE)
            
            plotCorEnsPC2TSI <- 
              plotCorEns(coroutPC2TSI,
                         bins = 20,
                         legend.position =c(.85,.8),
                         f.sig.lab.position = c(.85,.6),
                         significance.option = "isospectral",
                         use.fdr = TRUE)+ggtitle("PC2 - TSI")
            
            coroutPC2 <- corEns(UACMatchingPC2.ae,UACMatchingPC2.PA,UACMatchingPC2.ae,UACMatchingPC2.PA,bin.step = 40,max.ens = 1000)
            corPlotPC2 <- plotCorEns(coroutPC2,
                                     legend.position = c(0.1, 0.8),
                                     significance.option = "isospectral")+ggtitle(NULL)
            
            above90 <- round(sum(coroutPC2$cor.df$r > 0.9)/nrow(coroutPC2$cor.df)*100,1)
            med <- round(median(coroutPC2$cor.df$r),2)

          }
          
          #PC3
          {
            coroutPC3TSI <- corEns(time.1 = UACMatchingPC3.ae,
                                   values.1 = UACMatchingPC3.PA,
                                   time.2 = t(matrix(replicate(1,TSILipd.TSI.Time$values),nrow = 1)),
                                   values.2 = TSILipd.TSI,
                                   bin.step = 40,
                                   max.ens = 1000000,
                                   min.obs = 5,
                                   isopersistent  = TRUE,
                                   isospectral = TRUE)
            
            plotCorEnsPC3TSI <- 
              plotCorEns(coroutPC3TSI,
                         bins = 20,
                         legend.position =c(.85,.8),
                         f.sig.lab.position = c(.85,.6),
                         significance.option = "isospectral",
                         use.fdr = TRUE)+ggtitle("PC3 - TSI")
            
            coroutPC3 <- corEns(UACMatchingPC3.ae,UACMatchingPC3.PA,UACMatchingPC3.ae,UACMatchingPC3.PA,bin.step = 40,max.ens = 1000)
            corPlotPC3 <- plotCorEns(coroutPC3,
                                     legend.position = c(0.1, 0.8),
                                     significance.option = "isospectral")+ggtitle(NULL)
            
            above90 <- round(sum(coroutPC3$cor.df$r > 0.9)/nrow(coroutPC3$cor.df)*100,1)
            med <- round(median(coroutPC3$cor.df$r),2)
            
          }
          
          #Combine plots
          plotCorEnsPC1TSI | plotCorEnsPC2TSI | plotCorEnsPC3TSI
        }
      }
    
      #Correlation analysis-PC-CC
      {
        #Correlation analysis
        {
          #PC1
          {
            coroutPC1CC <- corEns(time.1 = UACMatchingPC1.ae,
                                   values.1 = UACMatchingPC1.PA,
                                   time.2 = t(matrix(replicate(1,CCLipd.CC.Time$values),nrow = 1)),
                                   values.2 = CCLipd.CC,
                                   bin.step = 40,
                                   max.ens = 1000000,
                                   min.obs = 5,
                                   isopersistent  = TRUE,
                                   isospectral = TRUE)
            
            plotCorEnsPC1CC <- 
              plotCorEns(coroutPC1CC,
                         bins = 20,
                         legend.position =c(.85,.8),
                         f.sig.lab.position = c(.85,.6),
                         significance.option = "isospectral",
                         use.fdr = TRUE)+ggtitle("PC1 - CC")
          }
          
          #PC2
          {
            coroutPC2CC <- corEns(time.1 = UACMatchingPC2.ae,
                                   values.1 = UACMatchingPC2.PA,
                                   time.2 = t(matrix(replicate(1,CCLipd.CC.Time$values),nrow = 1)),
                                   values.2 = CCLipd.CC,
                                   bin.step = 40,
                                   max.ens = 1000000,
                                   min.obs = 5,
                                   isopersistent  = TRUE,
                                   isospectral = TRUE)
            
            plotCorEnsPC2CC <- 
              plotCorEns(coroutPC2CC,
                         bins = 20,
                         legend.position =c(.85,.8),
                         f.sig.lab.position = c(.85,.6),
                         significance.option = "isospectral",
                         use.fdr = TRUE)+ggtitle("PC2 - CC")
          }
          
          #PC3
          {
            coroutPC3CC <- corEns(time.1 = UACMatchingPC3.ae,
                                   values.1 = UACMatchingPC3.PA,
                                   time.2 = t(matrix(replicate(1,CCLipd.CC.Time$values),nrow = 1)),
                                   values.2 = CCLipd.CC,
                                   bin.step = 40,
                                   max.ens = 1000000,
                                   min.obs = 5,
                                   isopersistent  = TRUE,
                                   isospectral = TRUE)
            
            plotCorEnsPC3CC <- 
              plotCorEns(coroutPC3CC,
                         bins = 20,
                         legend.position =c(.85,.8),
                         f.sig.lab.position = c(.85,.6),
                         significance.option = "isospectral",
                         use.fdr = TRUE)+ggtitle("PC3 - CC")
          }
          
          #Combine plots
          plotCorEnsPC1CC | plotCorEnsPC2CC | plotCorEnsPC3CC
        }
      }
    }
    
    #Sectioned correlation (moving correlation)
    {
      c16 <- makeCluster(8, type = 'SOCK')
      registerDoSNOW(c16)
      
      #TSI
      {
        #Sphagnum
        {
          #Setup sub-group size and step
          movCorreWindowHeadCut <- 10
          movCorreWindowSize <- 25
          movCorreWindowStep <- 4
          
          #Setup significant level for moving correlation confidence interval check
          sigVal <- 0.05
          
          #Create a matrix to store the median, 95% range and 95% confidential interval
          #for sectioned correlation
          movCorreMatChronUnSphTSI <- matrix(nrow = floor((nrow(UACSph) - movCorreWindowSize - movCorreWindowHeadCut) / movCorreWindowStep),
                                             ncol = 10,
                                             dimnames = list(c(),
                                                             c('Ind1', 'Ind2',
                                                               'Time', 'Time1', 'Time2',
                                                               'P1','M','P2',
                                                               'Sig1', 'Sig2')))
          
          movCorreMatChronUnSphTSI[, 1] <- seq(from = 0, to = nrow(movCorreMatChronUnSphTSI) - 1,
                                               by = 1) * movCorreWindowStep + 1
          
          movCorreMatChronUnSphTSI[, 2] <- movCorreMatChronUnSphTSI[, 1] + movCorreWindowSize
          
          movCorreMatChronUnSphTSI[, 3:10] <-
            foreach(paraInd = c(1:nrow(movCorreMatChronUnSphTSI)), .packages = c("zoo", 'geoChronR'),
                    .combine = 'rbind') %dopar% 
            {
              try(
                {
                  tempTime1 <- min(UACMatchingSph.ae$values[movCorreMatChronUnSphTSI[paraInd, 1]:movCorreMatChronUnSphTSI[paraInd, 2],])
                  tempTime2 <- max(UACMatchingSph.ae$values[movCorreMatChronUnSphTSI[paraInd, 1]:movCorreMatChronUnSphTSI[paraInd, 2],])
                  
                  tempTime3 <- quantile(UACMatchingSph.ae$values[movCorreMatChronUnSphTSI[paraInd, 1]:movCorreMatChronUnSphTSI[paraInd, 2],], probs = 0.025)
                  tempTime4 <- quantile(UACMatchingSph.ae$values[movCorreMatChronUnSphTSI[paraInd, 1]:movCorreMatChronUnSphTSI[paraInd, 2],], probs = 0.975)
                  
                  tempIndTSI <- which((TSIRaw[,1] >= tempTime1) & 
                                        (TSIRaw[,1] <= tempTime2))
                  
                  coroutSphTSITemp <- corEns(time.1 = UACMatchingSph.ae$values[movCorreMatChronUnSphTSI[paraInd, 1]:movCorreMatChronUnSphTSI[paraInd, 2], ],
                                             values.1 = UACMatchingSph.PA$values[movCorreMatChronUnSphTSI[paraInd, 1]:movCorreMatChronUnSphTSI[paraInd, 2]],
                                             time.2 = t(matrix(replicate(1000,TSIRaw[tempIndTSI, 1]),nrow = 1000)),
                                             values.2 = TSIRaw[tempIndTSI, 2],
                                             bin.step = 50,
                                             max.ens = 1000,
                                             min.obs = 3,
                                             isopersistent  = TRUE,
                                             isospectral = TRUE)
                  
                  timeMedian <- median(apply(UACMatchingSph.ae$values[movCorreMatChronUnSphTSI[paraInd, 1]:movCorreMatChronUnSphTSI[paraInd, 2], ], 2, median))
                  
                  tempSigCheck <- sort(coroutSphTSITemp[[1]][[5]]) < sigVal
                  
                  tempSigp <- sort(coroutSphTSITemp[[1]][[1]])[c(max(min(which(tempSigCheck == FALSE)) - 1, 1), 
                                                                 min(max(which(tempSigCheck == FALSE)) + 1, 
                                                                     length(sort(coroutSphTSITemp[[1]][[1]]))
                                                                 )
                  )
                  ]
                  
                  return(c(timeMedian, tempTime3, tempTime4, 
                           coroutSphTSITemp[[2]][[5]][1],
                           coroutSphTSITemp[[2]][[5]][3],
                           coroutSphTSITemp[[2]][[5]][5],
                           tempSigp))
                },
                silent = TRUE
              )
              
            }
          
          movCorreMatChronUnSphTSI <- apply(movCorreMatChronUnSphTSI, c(1,2), as.numeric)
        }
        
        #Alnus
        {
          #Setup sub-group size and step
          movCorreWindowHeadCut <- 10
          movCorreWindowSize <- 25
          movCorreWindowStep <- 4
          
          #Setup significant level for moving correlation confidence interval check
          sigVal <- 0.05
          
          #Create a matrix to store the median, 95% range and 95% confidential interval
          #for sectioned correlation
          movCorreMatChronUnAlnTSI <- matrix(nrow = floor((nrow(UACAln) - movCorreWindowSize - movCorreWindowHeadCut) / movCorreWindowStep),
                                             ncol = 10,
                                             dimnames = list(c(),
                                                             c('Ind1', 'Ind2',
                                                               'Time', 'Time1', 'Time2',
                                                               'P1','M','P2',
                                                               'Sig1', 'Sig2')))
          
          movCorreMatChronUnAlnTSI[, 1] <- seq(from = 0, to = nrow(movCorreMatChronUnAlnTSI) - 1,
                                               by = 1) * movCorreWindowStep + 1
          
          movCorreMatChronUnAlnTSI[, 2] <- movCorreMatChronUnAlnTSI[, 1] + movCorreWindowSize
          
          movCorreMatChronUnAlnTSI[, 3:10] <-
            foreach(paraInd = c(1:nrow(movCorreMatChronUnAlnTSI)), .packages = c("zoo", 'geoChronR'),
                    .combine = 'rbind') %dopar% 
            {
              try(
                {
                  tempTime1 <- min(UACMatchingAln.ae$values[movCorreMatChronUnAlnTSI[paraInd, 1]:movCorreMatChronUnAlnTSI[paraInd, 2],])
                  tempTime2 <- max(UACMatchingAln.ae$values[movCorreMatChronUnAlnTSI[paraInd, 1]:movCorreMatChronUnAlnTSI[paraInd, 2],])
                  
                  tempTime3 <- quantile(UACMatchingAln.ae$values[movCorreMatChronUnAlnTSI[paraInd, 1]:movCorreMatChronUnAlnTSI[paraInd, 2],], probs = 0.025)
                  tempTime4 <- quantile(UACMatchingAln.ae$values[movCorreMatChronUnAlnTSI[paraInd, 1]:movCorreMatChronUnAlnTSI[paraInd, 2],], probs = 0.975)
                  
                  tempIndTSI <- which((TSIRaw[,1] >= tempTime1) & 
                                        (TSIRaw[,1] <= tempTime2))
                  
                  coroutAlnTSITemp <- corEns(time.1 = UACMatchingAln.ae$values[movCorreMatChronUnAlnTSI[paraInd, 1]:movCorreMatChronUnAlnTSI[paraInd, 2], ],
                                             values.1 = UACMatchingAln.PA$values[movCorreMatChronUnAlnTSI[paraInd, 1]:movCorreMatChronUnAlnTSI[paraInd, 2]],
                                             time.2 = t(matrix(replicate(1000,TSIRaw[tempIndTSI, 1]),nrow = 1000)),
                                             values.2 = TSIRaw[tempIndTSI, 2],
                                             bin.step = 50,
                                             max.ens = 1000,
                                             min.obs = 3,
                                             isopersistent  = TRUE,
                                             isospectral = TRUE)
                  
                  timeMedian <- median(apply(UACMatchingAln.ae$values[movCorreMatChronUnAlnTSI[paraInd, 1]:movCorreMatChronUnAlnTSI[paraInd, 2], ], 2, median))
                  
                  tempSigCheck <- sort(coroutAlnTSITemp[[1]][[5]]) < sigVal
                  
                  tempSigp <- sort(coroutAlnTSITemp[[1]][[1]])[c(max(min(which(tempSigCheck == FALSE)) - 1, 1), 
                                                                 min(max(which(tempSigCheck == FALSE)) + 1, 
                                                                     length(sort(coroutAlnTSITemp[[1]][[1]]))
                                                                 )
                  )
                  ]
                  
                  return(c(timeMedian, tempTime3, tempTime4, 
                           coroutAlnTSITemp[[2]][[5]][1],
                           coroutAlnTSITemp[[2]][[5]][3],
                           coroutAlnTSITemp[[2]][[5]][5],
                           tempSigp))
                },
                silent = TRUE
              )
            }
          
          movCorreMatChronUnAlnTSI <- apply(movCorreMatChronUnAlnTSI, c(1,2), as.numeric)
          
        }
        
        #Calluna
        {
          #Setup sub-group size and step
          movCorreWindowHeadCut <- 10
          movCorreWindowSize <- 25
          movCorreWindowStep <- 4
          
          #Setup significant level for moving correlation confidence interval check
          sigVal <- 0.05
          
          #Create a matrix to store the median, 95% range and 95% confidential interval
          #for sectioned correlation
          movCorreMatChronUnCalTSI <- matrix(nrow = floor((nrow(UACCal) - movCorreWindowSize - movCorreWindowHeadCut) / movCorreWindowStep),
                                             ncol = 10,
                                             dimnames = list(c(),
                                                             c('Ind1', 'Ind2',
                                                               'Time', 'Time1', 'Time2',
                                                               'P1','M','P2',
                                                               'Sig1', 'Sig2')))
          
          movCorreMatChronUnCalTSI[, 1] <- seq(from = 0, to = nrow(movCorreMatChronUnCalTSI) - 1,
                                               by = 1) * movCorreWindowStep + 1
          
          movCorreMatChronUnCalTSI[, 2] <- movCorreMatChronUnCalTSI[, 1] + movCorreWindowSize
          
          movCorreMatChronUnCalTSI[, 3:10] <-
            foreach(paraInd = c(1:nrow(movCorreMatChronUnCalTSI)), .packages = c("zoo", 'geoChronR'),
                    .combine = 'rbind') %dopar% 
            {
              try(
                {
                  tempTime1 <- min(UACMatchingCal.ae$values[movCorreMatChronUnCalTSI[paraInd, 1]:movCorreMatChronUnCalTSI[paraInd, 2],])
                  tempTime2 <- max(UACMatchingCal.ae$values[movCorreMatChronUnCalTSI[paraInd, 1]:movCorreMatChronUnCalTSI[paraInd, 2],])
                  
                  tempTime3 <- quantile(UACMatchingCal.ae$values[movCorreMatChronUnCalTSI[paraInd, 1]:movCorreMatChronUnCalTSI[paraInd, 2],], probs = 0.025)
                  tempTime4 <- quantile(UACMatchingCal.ae$values[movCorreMatChronUnCalTSI[paraInd, 1]:movCorreMatChronUnCalTSI[paraInd, 2],], probs = 0.975)
                  
                  tempIndTSI <- which((TSIRaw[,1] >= tempTime1) & 
                                        (TSIRaw[,1] <= tempTime2))
                  
                  coroutCalTSITemp <- corEns(time.1 = UACMatchingCal.ae$values[movCorreMatChronUnCalTSI[paraInd, 1]:movCorreMatChronUnCalTSI[paraInd, 2], ],
                                             values.1 = UACMatchingCal.PA$values[movCorreMatChronUnCalTSI[paraInd, 1]:movCorreMatChronUnCalTSI[paraInd, 2]],
                                             time.2 = t(matrix(replicate(1000,TSIRaw[tempIndTSI, 1]),nrow = 1000)),
                                             values.2 = TSIRaw[tempIndTSI, 2],
                                             bin.step = 50,
                                             max.ens = 1000,
                                             min.obs = 3,
                                             isopersistent  = TRUE,
                                             isospectral = TRUE)
                  
                  timeMedian <- median(apply(UACMatchingCal.ae$values[movCorreMatChronUnCalTSI[paraInd, 1]:movCorreMatChronUnCalTSI[paraInd, 2], ], 2, median))
                  
                  tempSigCheck <- sort(coroutCalTSITemp[[1]][[5]]) < sigVal
                  
                  tempSigp <- sort(coroutCalTSITemp[[1]][[1]])[c(max(min(which(tempSigCheck == FALSE)) - 1, 1), 
                                                                 min(max(which(tempSigCheck == FALSE)) + 1, 
                                                                     length(sort(coroutCalTSITemp[[1]][[1]]))
                                                                 )
                  )
                  ]
                  
                  return(c(timeMedian, tempTime3, tempTime4, 
                           coroutCalTSITemp[[2]][[5]][1],
                           coroutCalTSITemp[[2]][[5]][3],
                           coroutCalTSITemp[[2]][[5]][5],
                           tempSigp))
                },
                silent = TRUE
              )
              
            }
          
          movCorreMatChronUnCalTSI <- apply(movCorreMatChronUnCalTSI, c(1,2), as.numeric)
        }
        
        #Plot the sectioned correlation
        {
          #Sphagnum
          movChronUnSphTSIPlot <- ggplot() +
            geom_point(aes(y = movCorreMatChronUnSphTSI[, 7],
                           x = movCorreMatChronUnSphTSI[, 3]), color = 'black',
                       size = 2) +
            geom_errorbar(aes(x = movCorreMatChronUnSphTSI[, 3],
                              y = movCorreMatChronUnSphTSI[, 7],
                              ymin = movCorreMatChronUnSphTSI[, 6],
                              ymax = movCorreMatChronUnSphTSI[, 8]),
                          fill = 'black', linewidth = 1, alpha = 0.3,
                          width = 20) +
            geom_errorbarh(aes(xmin = movCorreMatChronUnSphTSI[, 4],
                               xmax = movCorreMatChronUnSphTSI[, 5],
                               y = movCorreMatChronUnSphTSI[, 7],
                               x = movCorreMatChronUnSphTSI[, 3]
            ),
            fill = 'black', linewidth = 2, 
            alpha = 0.3) +
            geom_line(aes(y = movCorreMatChronUnSphTSI[, 10],
                          x = movCorreMatChronUnSphTSI[, 3]), color = 'red', linetype = 2, alpha = 0.6) +
            geom_line(aes(y = movCorreMatChronUnSphTSI[, 9],
                          x = movCorreMatChronUnSphTSI[, 3]), color = 'red', linetype = 2, alpha = 0.6) +
            geom_hline(yintercept = 0, color = 'grey', alpha = 0.6) +
            scale_y_continuous(breaks = seq(from = -10, to = 10, by = 0.25),
                               limits = c(0,1),
                               position = 'left') +
            scale_x_reverse(limits = c(2750, -100), 
                            breaks = seq(from = 4000, to = -100, by = -100),
                            labels = c(4000, rep('', 4),
                                       3500,rep('', 4),
                                       3000, rep('', 4),
                                       2500,rep('', 4),
                                       2000, rep('', 4),
                                       1500,rep('', 4),
                                       1000, rep('', 4),
                                       500,rep('', 4),
                                       0, ''),
                            sec.axis=sec_axis(~., 
                                              breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                              labels=c(2000, rep('', 4),
                                                       1500,rep('', 4),
                                                       1000, rep('', 4),
                                                       500,rep('', 4),
                                                       0, rep('', 4),
                                                       500,rep('', 4),
                                                       1000, rep('', 4),
                                                       1500,rep('', 4),
                                                       2000)
                            )
            )  +
            guides(y = "axis_truncated", x = "axis_truncated") +
            labs(y = 'Sph&TSI', x = 'Cal years BP', title = '') + 
            theme_classic() + 
            theme(plot.title = element_blank(), 
                  axis.text.x = element_text(size = 10, family = 'arial'),
                  axis.title.x = element_text(size = 12, family = 'arial'),
                  axis.text.y = element_text(size = 10, family = 'arial'),
                  axis.title.y = element_text(size = 12, family = 'arial')
            )
          
          #Alnus
          movChronUnAlnTSIPlot <- ggplot() +
            geom_point(aes(y = movCorreMatChronUnAlnTSI[, 7],
                           x = movCorreMatChronUnAlnTSI[, 3]), color = 'black',
                       size = 2) +
            geom_errorbar(aes(x = movCorreMatChronUnAlnTSI[, 3],
                              y = movCorreMatChronUnAlnTSI[, 7],
                              ymin = movCorreMatChronUnAlnTSI[, 6],
                              ymax = movCorreMatChronUnAlnTSI[, 8]),
                          fill = 'black', linewidth = 1, alpha = 0.3,
                          width = 20) +
            geom_errorbarh(aes(xmin = movCorreMatChronUnAlnTSI[, 4],
                               xmax = movCorreMatChronUnAlnTSI[, 5],
                               y = movCorreMatChronUnAlnTSI[, 7],
                               x = movCorreMatChronUnAlnTSI[, 3]
            ),
            fill = 'black', linewidth = 2, 
            alpha = 0.3) +
            geom_line(aes(y = movCorreMatChronUnAlnTSI[, 10],
                          x = movCorreMatChronUnAlnTSI[, 3]), color = 'red', linetype = 2, alpha = 0.6) +
            geom_line(aes(y = movCorreMatChronUnAlnTSI[, 9],
                          x = movCorreMatChronUnAlnTSI[, 3]), color = 'red', linetype = 2, alpha = 0.6) +
            geom_hline(yintercept = 0, color = 'grey', alpha = 0.6) +
            scale_y_continuous(breaks = seq(from = -10, to = 10, by = 0.25),
                               limits = c(0,1),
                               position = 'left') +
            scale_x_reverse(limits = c(2750, -100), 
                            breaks = seq(from = 4000, to = -100, by = -100),
                            labels = c(4000, rep('', 4),
                                       3500,rep('', 4),
                                       3000, rep('', 4),
                                       2500,rep('', 4),
                                       2000, rep('', 4),
                                       1500,rep('', 4),
                                       1000, rep('', 4),
                                       500,rep('', 4),
                                       0, ''),
                            sec.axis=sec_axis(~., 
                                              breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                              labels=c(2000, rep('', 4),
                                                       1500,rep('', 4),
                                                       1000, rep('', 4),
                                                       500,rep('', 4),
                                                       0, rep('', 4),
                                                       500,rep('', 4),
                                                       1000, rep('', 4),
                                                       1500,rep('', 4),
                                                       2000)
                            )
            )  +
            guides(y = "axis_truncated", x = "axis_truncated") +
            labs(y = 'Aln&TSI', x = 'Cal years BP', title = '') + 
            theme_classic() + 
            theme(plot.title = element_blank(), 
                  axis.text.x = element_text(size = 10, family = 'arial'),
                  axis.title.x = element_text(size = 12, family = 'arial'),
                  axis.text.y = element_text(size = 10, family = 'arial'),
                  axis.title.y = element_text(size = 12, family = 'arial')
            )
          
          #Calluna
          movChronUnCalTSIPlot <- ggplot() +
            geom_point(aes(y = movCorreMatChronUnCalTSI[, 7],
                           x = movCorreMatChronUnCalTSI[, 3]), color = 'black',
                       size = 2) +
            geom_errorbar(aes(x = movCorreMatChronUnCalTSI[, 3],
                              y = movCorreMatChronUnCalTSI[, 7],
                              ymin = movCorreMatChronUnCalTSI[, 6],
                              ymax = movCorreMatChronUnCalTSI[, 8]),
                          fill = 'black', linewidth = 1, alpha = 0.3,
                          width = 20) +
            geom_errorbarh(aes(xmin = movCorreMatChronUnCalTSI[, 4],
                               xmax = movCorreMatChronUnCalTSI[, 5],
                               y = movCorreMatChronUnCalTSI[, 7],
                               x = movCorreMatChronUnCalTSI[, 3]
            ),
            fill = 'black', linewidth = 2, 
            alpha = 0.3) +
            geom_line(aes(y = movCorreMatChronUnCalTSI[, 10],
                          x = movCorreMatChronUnCalTSI[, 3]), color = 'red', linetype = 2, alpha = 0.6) +
            geom_line(aes(y = movCorreMatChronUnCalTSI[, 9],
                          x = movCorreMatChronUnCalTSI[, 3]), color = 'red', linetype = 2, alpha = 0.6) +
            geom_hline(yintercept = 0, color = 'grey', alpha = 0.6) +
            scale_y_continuous(breaks = seq(from = -10, to = 10, by = 0.25),
                               limits = c(0,1),
                               position = 'left') +
            scale_x_reverse(limits = c(2750, -100), 
                            breaks = seq(from = 4000, to = -100, by = -100),
                            labels = c(4000, rep('', 4),
                                       3500,rep('', 4),
                                       3000, rep('', 4),
                                       2500,rep('', 4),
                                       2000, rep('', 4),
                                       1500,rep('', 4),
                                       1000, rep('', 4),
                                       500,rep('', 4),
                                       0, ''),
                            sec.axis=sec_axis(~., 
                                              breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                              labels=c(2000, rep('', 4),
                                                       1500,rep('', 4),
                                                       1000, rep('', 4),
                                                       500,rep('', 4),
                                                       0, rep('', 4),
                                                       500,rep('', 4),
                                                       1000, rep('', 4),
                                                       1500,rep('', 4),
                                                       2000)
                            )
            )  +
            guides(y = "axis_truncated", x = "axis_truncated") +
            labs(y = 'Cal&TSI', x = 'Cal years BP', title = '') + 
            theme_classic() + 
            theme(plot.title = element_blank(), 
                  axis.text.x = element_text(size = 10, family = 'arial'),
                  axis.title.x = element_text(size = 12, family = 'arial'),
                  axis.text.y = element_text(size = 10, family = 'arial'),
                  axis.title.y = element_text(size = 12, family = 'arial')
            )
        }
        
        #Combine plot
        {
          ggarrange(UACMatchingSph.PA.ts.plot +
                      geom_hline(yintercept = UACSphChronUnMean, 
                                 linetype = 1, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACSphChronUnMean - UACSphChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACSphChronUnMean - 2 * UACSphChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACSphChronUnMean + UACSphChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACSphChronUnMean + 2* UACSphChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) ,
                    WholeTSI +
                      theme(
                        axis.line.x = element_blank(),
                        axis.ticks.x = element_blank(),
                        axis.title.x = element_blank(),
                        axis.text.x = element_blank(),
                      ),
                    movChronUnSphTSIPlot,
                    heights = c(1,1,1),
                    widths = c(1,1,1),
                    ncol = 1, nrow = 3,
                    align = "v")
          
          ggarrange(UACMatchingAln.PA.ts.plot +
                      geom_hline(yintercept = UACAlnChronUnMean, 
                                 linetype = 1, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACAlnChronUnMean - UACAlnChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACAlnChronUnMean - 2 * UACAlnChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACAlnChronUnMean + UACAlnChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACAlnChronUnMean + 2* UACAlnChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) ,
                    WholeTSI +
                      theme(
                        axis.line.x = element_blank(),
                        axis.ticks.x = element_blank(),
                        axis.title.x = element_blank(),
                        axis.text.x = element_blank(),
                      ),
                    movChronUnAlnTSIPlot,
                    heights = c(1,1,1),
                    widths = c(1,1,1),
                    ncol = 1, nrow = 3,
                    align = "v")
          
          ggarrange(UACMatchingCal.PA.ts.plot +
                      geom_hline(yintercept = UACCalChronUnMean, 
                                 linetype = 1, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACCalChronUnMean - UACCalChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACCalChronUnMean - 2 * UACCalChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACCalChronUnMean + UACCalChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACCalChronUnMean + 2* UACCalChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) ,
                    WholeTSI +
                      theme(
                        axis.line.x = element_blank(),
                        axis.ticks.x = element_blank(),
                        axis.title.x = element_blank(),
                        axis.text.x = element_blank(),
                      ),
                    movChronUnCalTSIPlot,
                    heights = c(1,1,1),
                    widths = c(1,1,1),
                    ncol = 1, nrow = 3,
                    align = "v")
        }
      }
      
      #CC
      {
        #Sphagnum
        {
          #Setup sub-group size and step
          movCorreWindowHeadCut <- 60
          movCorreWindowSize <- 25
          movCorreWindowStep <- 4
          
          #Setup significant level for moving correlation confidence interval check
          sigVal <- 0.05
          
          #Create a matrix to store the median, 95% range and 95% confidential interval
          #for sectioned correlation
          movCorreMatChronUnSphCC <- matrix(nrow = floor((nrow(UACSph) - movCorreWindowSize - movCorreWindowHeadCut) / movCorreWindowStep),
                                             ncol = 10,
                                             dimnames = list(c(),
                                                             c('Ind1', 'Ind2',
                                                               'Time', 'Time1', 'Time2',
                                                               'P1','M','P2',
                                                               'Sig1', 'Sig2')))
          
          movCorreMatChronUnSphCC[, 1] <- seq(from = 0, to = nrow(movCorreMatChronUnSphCC) - 1,
                                               by = 1) * movCorreWindowStep + 1
          
          movCorreMatChronUnSphCC[, 2] <- movCorreMatChronUnSphCC[, 1] + movCorreWindowSize
          
          movCorreMatChronUnSphCC[, 3:10] <-
            foreach(paraInd = c(1:nrow(movCorreMatChronUnSphCC)), .packages = c("zoo", 'geoChronR'),
                    .combine = 'rbind') %dopar% 
            {
              try(
                {
                  tempTime1 <- min(UACMatchingSph.ae$values[movCorreMatChronUnSphCC[paraInd, 1]:movCorreMatChronUnSphCC[paraInd, 2],])
                  tempTime2 <- max(UACMatchingSph.ae$values[movCorreMatChronUnSphCC[paraInd, 1]:movCorreMatChronUnSphCC[paraInd, 2],])
                  
                  tempTime3 <- quantile(UACMatchingSph.ae$values[movCorreMatChronUnSphCC[paraInd, 1]:movCorreMatChronUnSphCC[paraInd, 2],], probs = 0.025)
                  tempTime4 <- quantile(UACMatchingSph.ae$values[movCorreMatChronUnSphCC[paraInd, 1]:movCorreMatChronUnSphCC[paraInd, 2],], probs = 0.975)
                  
                  tempIndCC <- which((CCRaw[,1] >= tempTime1) & 
                                        (CCRaw[,1] <= tempTime2))
                  
                  coroutSphCCTemp <- corEns(time.1 = UACMatchingSph.ae$values[movCorreMatChronUnSphCC[paraInd, 1]:movCorreMatChronUnSphCC[paraInd, 2], ],
                                             values.1 = UACMatchingSph.PA$values[movCorreMatChronUnSphCC[paraInd, 1]:movCorreMatChronUnSphCC[paraInd, 2]],
                                             time.2 = t(matrix(replicate(1000,CCRaw[tempIndCC, 1]),nrow = 1000)),
                                             values.2 = CCRaw[tempIndCC, 2],
                                             bin.step = 50,
                                             max.ens = 1000,
                                             min.obs = 3,
                                             isopersistent  = TRUE,
                                             isospectral = TRUE)
                  
                  timeMedian <- median(apply(UACMatchingSph.ae$values[movCorreMatChronUnSphCC[paraInd, 1]:movCorreMatChronUnSphCC[paraInd, 2], ], 2, median))
                  
                  tempSigCheck <- sort(coroutSphCCTemp[[1]][[5]]) < sigVal
                  
                  tempSigp <- sort(coroutSphCCTemp[[1]][[1]])[c(max(min(which(tempSigCheck == FALSE)) - 1, 1), 
                                                                 min(max(which(tempSigCheck == FALSE)) + 1, 
                                                                     length(sort(coroutSphCCTemp[[1]][[1]]))
                                                                 )
                  )
                  ]
                  
                  return(c(timeMedian, tempTime3, tempTime4, 
                           coroutSphCCTemp[[2]][[5]][1],
                           coroutSphCCTemp[[2]][[5]][3],
                           coroutSphCCTemp[[2]][[5]][5],
                           tempSigp))
                },
                silent = TRUE
              )
              
            }
          
          movCorreMatChronUnSphCC <- apply(movCorreMatChronUnSphCC, c(1,2), as.numeric)
        }
        
        #Alnus
        {
          #Setup sub-group size and step
          movCorreWindowHeadCut <- 50
          movCorreWindowSize <- 25
          movCorreWindowStep <- 4
          
          #Setup significant level for moving correlation confidence interval check
          sigVal <- 0.05
          
          #Create a matrix to store the median, 95% range and 95% confidential interval
          #for sectioned correlation
          movCorreMatChronUnAlnCC <- matrix(nrow = floor((nrow(UACAln) - movCorreWindowSize - movCorreWindowHeadCut) / movCorreWindowStep),
                                             ncol = 10,
                                             dimnames = list(c(),
                                                             c('Ind1', 'Ind2',
                                                               'Time', 'Time1', 'Time2',
                                                               'P1','M','P2',
                                                               'Sig1', 'Sig2')))
          
          movCorreMatChronUnAlnCC[, 1] <- seq(from = 0, to = nrow(movCorreMatChronUnAlnCC) - 1,
                                               by = 1) * movCorreWindowStep + 1
          
          movCorreMatChronUnAlnCC[, 2] <- movCorreMatChronUnAlnCC[, 1] + movCorreWindowSize
          
          movCorreMatChronUnAlnCC[, 3:10] <-
            foreach(paraInd = c(1:nrow(movCorreMatChronUnAlnCC)), .packages = c("zoo", 'geoChronR'),
                    .combine = 'rbind') %dopar% 
            {
              try(
                {
                  tempTime1 <- min(UACMatchingAln.ae$values[movCorreMatChronUnAlnCC[paraInd, 1]:movCorreMatChronUnAlnCC[paraInd, 2],])
                  tempTime2 <- max(UACMatchingAln.ae$values[movCorreMatChronUnAlnCC[paraInd, 1]:movCorreMatChronUnAlnCC[paraInd, 2],])
                  
                  tempTime3 <- quantile(UACMatchingAln.ae$values[movCorreMatChronUnAlnCC[paraInd, 1]:movCorreMatChronUnAlnCC[paraInd, 2],], probs = 0.025)
                  tempTime4 <- quantile(UACMatchingAln.ae$values[movCorreMatChronUnAlnCC[paraInd, 1]:movCorreMatChronUnAlnCC[paraInd, 2],], probs = 0.975)
                  
                  tempIndCC <- which((CCRaw[,1] >= tempTime1) & 
                                        (CCRaw[,1] <= tempTime2))
                  
                  coroutAlnCCTemp <- corEns(time.1 = UACMatchingAln.ae$values[movCorreMatChronUnAlnCC[paraInd, 1]:movCorreMatChronUnAlnCC[paraInd, 2], ],
                                             values.1 = UACMatchingAln.PA$values[movCorreMatChronUnAlnCC[paraInd, 1]:movCorreMatChronUnAlnCC[paraInd, 2]],
                                             time.2 = t(matrix(replicate(1000,CCRaw[tempIndCC, 1]),nrow = 1000)),
                                             values.2 = CCRaw[tempIndCC, 2],
                                             bin.step = 50,
                                             max.ens = 1000,
                                             min.obs = 3,
                                             isopersistent  = TRUE,
                                             isospectral = TRUE)
                  
                  timeMedian <- median(apply(UACMatchingAln.ae$values[movCorreMatChronUnAlnCC[paraInd, 1]:movCorreMatChronUnAlnCC[paraInd, 2], ], 2, median))
                  
                  tempSigCheck <- sort(coroutAlnCCTemp[[1]][[5]]) < sigVal
                  
                  tempSigp <- sort(coroutAlnCCTemp[[1]][[1]])[c(max(min(which(tempSigCheck == FALSE)) - 1, 1), 
                                                                 min(max(which(tempSigCheck == FALSE)) + 1, 
                                                                     length(sort(coroutAlnCCTemp[[1]][[1]]))
                                                                 )
                  )
                  ]
                  
                  return(c(timeMedian, tempTime3, tempTime4, 
                           coroutAlnCCTemp[[2]][[5]][1],
                           coroutAlnCCTemp[[2]][[5]][3],
                           coroutAlnCCTemp[[2]][[5]][5],
                           tempSigp))
                },
                silent = TRUE
              )
            }
          
          movCorreMatChronUnAlnCC <- apply(movCorreMatChronUnAlnCC, c(1,2), as.numeric)
          
        }
        
        #Calluna
        {
          #Setup sub-group size and step
          movCorreWindowHeadCut <- 100
          movCorreWindowSize <- 25
          movCorreWindowStep <- 4
          
          #Setup significant level for moving correlation confidence interval check
          sigVal <- 0.05
          
          #Create a matrix to store the median, 95% range and 95% confidential interval
          #for sectioned correlation
          movCorreMatChronUnCalCC <- matrix(nrow = floor((nrow(UACCal) - movCorreWindowSize - movCorreWindowHeadCut) / movCorreWindowStep),
                                             ncol = 10,
                                             dimnames = list(c(),
                                                             c('Ind1', 'Ind2',
                                                               'Time', 'Time1', 'Time2',
                                                               'P1','M','P2',
                                                               'Sig1', 'Sig2')))
          
          movCorreMatChronUnCalCC[, 1] <- seq(from = 0, to = nrow(movCorreMatChronUnCalCC) - 1,
                                               by = 1) * movCorreWindowStep + 1
          
          movCorreMatChronUnCalCC[, 2] <- movCorreMatChronUnCalCC[, 1] + movCorreWindowSize
          
          movCorreMatChronUnCalCC[, 3:10] <-
            foreach(paraInd = c(1:nrow(movCorreMatChronUnCalCC)), .packages = c("zoo", 'geoChronR'),
                    .combine = 'rbind') %dopar% 
            {
              try(
                {
                  tempTime1 <- min(UACMatchingCal.ae$values[movCorreMatChronUnCalCC[paraInd, 1]:movCorreMatChronUnCalCC[paraInd, 2],])
                  tempTime2 <- max(UACMatchingCal.ae$values[movCorreMatChronUnCalCC[paraInd, 1]:movCorreMatChronUnCalCC[paraInd, 2],])
                  
                  tempTime3 <- quantile(UACMatchingCal.ae$values[movCorreMatChronUnCalCC[paraInd, 1]:movCorreMatChronUnCalCC[paraInd, 2],], probs = 0.025)
                  tempTime4 <- quantile(UACMatchingCal.ae$values[movCorreMatChronUnCalCC[paraInd, 1]:movCorreMatChronUnCalCC[paraInd, 2],], probs = 0.975)
                  
                  tempIndCC <- which((CCRaw[,1] >= tempTime1) & 
                                        (CCRaw[,1] <= tempTime2))
                  
                  coroutCalCCTemp <- corEns(time.1 = UACMatchingCal.ae$values[movCorreMatChronUnCalCC[paraInd, 1]:movCorreMatChronUnCalCC[paraInd, 2], ],
                                             values.1 = UACMatchingCal.PA$values[movCorreMatChronUnCalCC[paraInd, 1]:movCorreMatChronUnCalCC[paraInd, 2]],
                                             time.2 = t(matrix(replicate(1000,CCRaw[tempIndCC, 1]),nrow = 1000)),
                                             values.2 = CCRaw[tempIndCC, 2],
                                             bin.step = 50,
                                             max.ens = 1000,
                                             min.obs = 3,
                                             isopersistent  = TRUE,
                                             isospectral = TRUE)
                  
                  timeMedian <- median(apply(UACMatchingCal.ae$values[movCorreMatChronUnCalCC[paraInd, 1]:movCorreMatChronUnCalCC[paraInd, 2], ], 2, median))
                  
                  tempSigCheck <- sort(coroutCalCCTemp[[1]][[5]]) < sigVal
                  
                  tempSigp <- sort(coroutCalCCTemp[[1]][[1]])[c(max(min(which(tempSigCheck == FALSE)) - 1, 1), 
                                                                 min(max(which(tempSigCheck == FALSE)) + 1, 
                                                                     length(sort(coroutCalCCTemp[[1]][[1]]))
                                                                 )
                  )
                  ]
                  
                  return(c(timeMedian, tempTime3, tempTime4, 
                           coroutCalCCTemp[[2]][[5]][1],
                           coroutCalCCTemp[[2]][[5]][3],
                           coroutCalCCTemp[[2]][[5]][5],
                           tempSigp))
                },
                silent = TRUE
              )
              
            }
          
          movCorreMatChronUnCalCC <- apply(movCorreMatChronUnCalCC, c(1,2), as.numeric)
        }
        
        #Plot the sectioned correlation
        {
          #Sphagnum
          movChronUnSphCCPlot <- ggplot() +
            geom_point(aes(y = movCorreMatChronUnSphCC[, 7],
                           x = movCorreMatChronUnSphCC[, 3]), color = 'black',
                       size = 2) +
            geom_errorbar(aes(x = movCorreMatChronUnSphCC[, 3],
                              y = movCorreMatChronUnSphCC[, 7],
                              ymin = movCorreMatChronUnSphCC[, 6],
                              ymax = movCorreMatChronUnSphCC[, 8]),
                          fill = 'black', linewidth = 1, alpha = 0.3,
                          width = 20) +
            geom_errorbarh(aes(xmin = movCorreMatChronUnSphCC[, 4],
                               xmax = movCorreMatChronUnSphCC[, 5],
                               y = movCorreMatChronUnSphCC[, 7],
                               x = movCorreMatChronUnSphCC[, 3]
            ),
            fill = 'black', linewidth = 2, 
            alpha = 0.3) +
            geom_line(aes(y = movCorreMatChronUnSphCC[, 10],
                          x = movCorreMatChronUnSphCC[, 3]), color = 'red', linetype = 2, alpha = 0.6) +
            geom_line(aes(y = movCorreMatChronUnSphCC[, 9],
                          x = movCorreMatChronUnSphCC[, 3]), color = 'red', linetype = 2, alpha = 0.6) +
            geom_hline(yintercept = 0, color = 'grey', alpha = 0.6) +
            scale_y_continuous(breaks = seq(from = -10, to = 10, by = 0.25),
                               limits = c(0,1),
                               position = 'left') +
            scale_x_reverse(limits = c(1000, -100), 
                            breaks = seq(from = 4000, to = -100, by = -100),
                            labels = c(4000, rep('', 4),
                                       3500,rep('', 4),
                                       3000, rep('', 4),
                                       2500,rep('', 4),
                                       2000, rep('', 4),
                                       1500,rep('', 4),
                                       1000, rep('', 4),
                                       500,rep('', 4),
                                       0, ''),
                            sec.axis=sec_axis(~., 
                                              breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                              labels=c(2000, rep('', 4),
                                                       1500,rep('', 4),
                                                       1000, rep('', 4),
                                                       500,rep('', 4),
                                                       0, rep('', 4),
                                                       500,rep('', 4),
                                                       1000, rep('', 4),
                                                       1500,rep('', 4),
                                                       2000)
                            )
            )  +
            guides(y = "axis_truncated", x = "axis_truncated") +
            labs(y = 'Sph&CC', x = 'Cal years BP', title = '') + 
            theme_classic() + 
            theme(plot.title = element_blank(), 
                  axis.text.x = element_text(size = 10, family = 'arial'),
                  axis.title.x = element_text(size = 12, family = 'arial'),
                  axis.text.y = element_text(size = 10, family = 'arial'),
                  axis.title.y = element_text(size = 12, family = 'arial')
            )
          
          #Alnus
          movChronUnAlnCCPlot <- ggplot() +
            geom_point(aes(y = movCorreMatChronUnAlnCC[, 7],
                           x = movCorreMatChronUnAlnCC[, 3]), color = 'black',
                       size = 2) +
            geom_errorbar(aes(x = movCorreMatChronUnAlnCC[, 3],
                              y = movCorreMatChronUnAlnCC[, 7],
                              ymin = movCorreMatChronUnAlnCC[, 6],
                              ymax = movCorreMatChronUnAlnCC[, 8]),
                          fill = 'black', linewidth = 1, alpha = 0.3,
                          width = 20) +
            geom_errorbarh(aes(xmin = movCorreMatChronUnAlnCC[, 4],
                               xmax = movCorreMatChronUnAlnCC[, 5],
                               y = movCorreMatChronUnAlnCC[, 7],
                               x = movCorreMatChronUnAlnCC[, 3]
            ),
            fill = 'black', linewidth = 2, 
            alpha = 0.3) +
            geom_line(aes(y = movCorreMatChronUnAlnCC[, 10],
                          x = movCorreMatChronUnAlnCC[, 3]), color = 'red', linetype = 2, alpha = 0.6) +
            geom_line(aes(y = movCorreMatChronUnAlnCC[, 9],
                          x = movCorreMatChronUnAlnCC[, 3]), color = 'red', linetype = 2, alpha = 0.6) +
            geom_hline(yintercept = 0, color = 'grey', alpha = 0.6) +
            scale_y_continuous(breaks = seq(from = -10, to = 10, by = 0.25),
                               limits = c(-0.1,1),
                               position = 'left') +
            scale_x_reverse(limits = c(1000, -100), 
                            breaks = seq(from = 4000, to = -100, by = -100),
                            labels = c(4000, rep('', 4),
                                       3500,rep('', 4),
                                       3000, rep('', 4),
                                       2500,rep('', 4),
                                       2000, rep('', 4),
                                       1500,rep('', 4),
                                       1000, rep('', 4),
                                       500,rep('', 4),
                                       0, ''),
                            sec.axis=sec_axis(~., 
                                              breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                              labels=c(2000, rep('', 4),
                                                       1500,rep('', 4),
                                                       1000, rep('', 4),
                                                       500,rep('', 4),
                                                       0, rep('', 4),
                                                       500,rep('', 4),
                                                       1000, rep('', 4),
                                                       1500,rep('', 4),
                                                       2000)
                            )
            )  +
            guides(y = "axis_truncated", x = "axis_truncated") +
            labs(y = 'Aln&CC', x = 'Cal years BP', title = '') + 
            theme_classic() + 
            theme(plot.title = element_blank(), 
                  axis.text.x = element_text(size = 10, family = 'arial'),
                  axis.title.x = element_text(size = 12, family = 'arial'),
                  axis.text.y = element_text(size = 10, family = 'arial'),
                  axis.title.y = element_text(size = 12, family = 'arial')
            )
          
          #Calluna
          movChronUnCalCCPlot <- ggplot() +
            geom_point(aes(y = movCorreMatChronUnCalCC[, 7],
                           x = movCorreMatChronUnCalCC[, 3]), color = 'black',
                       size = 2) +
            geom_errorbar(aes(x = movCorreMatChronUnCalCC[, 3],
                              y = movCorreMatChronUnCalCC[, 7],
                              ymin = movCorreMatChronUnCalCC[, 6],
                              ymax = movCorreMatChronUnCalCC[, 8]),
                          fill = 'black', linewidth = 1, alpha = 0.3,
                          width = 20) +
            geom_errorbarh(aes(xmin = movCorreMatChronUnCalCC[, 4],
                               xmax = movCorreMatChronUnCalCC[, 5],
                               y = movCorreMatChronUnCalCC[, 7],
                               x = movCorreMatChronUnCalCC[, 3]
            ),
            fill = 'black', linewidth = 2, 
            alpha = 0.3) +
            geom_line(aes(y = movCorreMatChronUnCalCC[, 10],
                          x = movCorreMatChronUnCalCC[, 3]), color = 'red', linetype = 2, alpha = 0.6) +
            geom_line(aes(y = movCorreMatChronUnCalCC[, 9],
                          x = movCorreMatChronUnCalCC[, 3]), color = 'red', linetype = 2, alpha = 0.6) +
            geom_hline(yintercept = 0, color = 'grey', alpha = 0.6) +
            scale_y_continuous(breaks = seq(from = -10, to = 10, by = 0.25),
                               limits = c(0,1),
                               position = 'left') +
            scale_x_reverse(limits = c(1000, -100), 
                            breaks = seq(from = 4000, to = -100, by = -100),
                            labels = c(4000, rep('', 4),
                                       3500,rep('', 4),
                                       3000, rep('', 4),
                                       2500,rep('', 4),
                                       2000, rep('', 4),
                                       1500,rep('', 4),
                                       1000, rep('', 4),
                                       500,rep('', 4),
                                       0, ''),
                            sec.axis=sec_axis(~., 
                                              breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                              labels=c(2000, rep('', 4),
                                                       1500,rep('', 4),
                                                       1000, rep('', 4),
                                                       500,rep('', 4),
                                                       0, rep('', 4),
                                                       500,rep('', 4),
                                                       1000, rep('', 4),
                                                       1500,rep('', 4),
                                                       2000)
                            )
            )  +
            guides(y = "axis_truncated", x = "axis_truncated") +
            labs(y = 'Cal&CC', x = 'Cal years BP', title = '') + 
            theme_classic() + 
            theme(plot.title = element_blank(), 
                  axis.text.x = element_text(size = 10, family = 'arial'),
                  axis.title.x = element_text(size = 12, family = 'arial'),
                  axis.text.y = element_text(size = 10, family = 'arial'),
                  axis.title.y = element_text(size = 12, family = 'arial')
            )
        }
        
        #Combine plot
        {
          ggarrange(UACMatchingSph.PA.ts.plot +
                      geom_hline(yintercept = UACSphChronUnMean, 
                                 linetype = 1, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACSphChronUnMean - UACSphChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACSphChronUnMean - 2 * UACSphChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACSphChronUnMean + UACSphChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACSphChronUnMean + 2* UACSphChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) +
                      scale_x_reverse(limits = c(1000, -100), 
                                      breaks = seq(from = 4000, to = -100, by = -100),
                                      labels = c(4000, rep('', 4),
                                                 3500,rep('', 4),
                                                 3000, rep('', 4),
                                                 2500,rep('', 4),
                                                 2000, rep('', 4),
                                                 1500,rep('', 4),
                                                 1000, rep('', 4),
                                                 500,rep('', 4),
                                                 0, ''),
                                      sec.axis=sec_axis(~., 
                                                        breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                                        labels=c(2000, rep('', 4),
                                                                 1500,rep('', 4),
                                                                 1000, rep('', 4),
                                                                 500,rep('', 4),
                                                                 0, rep('', 4),
                                                                 500,rep('', 4),
                                                                 1000, rep('', 4),
                                                                 1500,rep('', 4),
                                                                 2000)
                                      )
                      ),
                    WholeCC +
                      theme(
                        axis.line.x = element_blank(),
                        axis.ticks.x = element_blank(),
                        axis.title.x = element_blank(),
                        axis.text.x = element_blank(),
                      ),
                    movChronUnSphCCPlot,
                    heights = c(1,1,1),
                    widths = c(1,1,1),
                    ncol = 1, nrow = 3,
                    align = "v")
          
          ggarrange(UACMatchingAln.PA.ts.plot +
                      geom_hline(yintercept = UACAlnChronUnMean, 
                                 linetype = 1, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACAlnChronUnMean - UACAlnChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACAlnChronUnMean - 2 * UACAlnChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACAlnChronUnMean + UACAlnChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACAlnChronUnMean + 2* UACAlnChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) +
                      scale_x_reverse(limits = c(1000, -100), 
                                      breaks = seq(from = 4000, to = -100, by = -100),
                                      labels = c(4000, rep('', 4),
                                                 3500,rep('', 4),
                                                 3000, rep('', 4),
                                                 2500,rep('', 4),
                                                 2000, rep('', 4),
                                                 1500,rep('', 4),
                                                 1000, rep('', 4),
                                                 500,rep('', 4),
                                                 0, ''),
                                      sec.axis=sec_axis(~., 
                                                        breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                                        labels=c(2000, rep('', 4),
                                                                 1500,rep('', 4),
                                                                 1000, rep('', 4),
                                                                 500,rep('', 4),
                                                                 0, rep('', 4),
                                                                 500,rep('', 4),
                                                                 1000, rep('', 4),
                                                                 1500,rep('', 4),
                                                                 2000)
                                      )
                      ),
                    WholeCC +
                      theme(
                        axis.line.x = element_blank(),
                        axis.ticks.x = element_blank(),
                        axis.title.x = element_blank(),
                        axis.text.x = element_blank(),
                      ),
                    movChronUnAlnCCPlot,
                    heights = c(1,1,1),
                    widths = c(1,1,1),
                    ncol = 1, nrow = 3,
                    align = "v")
          
          ggarrange(UACMatchingCal.PA.ts.plot +
                      geom_hline(yintercept = UACCalChronUnMean, 
                                 linetype = 1, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACCalChronUnMean - UACCalChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACCalChronUnMean - 2 * UACCalChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACCalChronUnMean + UACCalChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACCalChronUnMean + 2* UACCalChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) +
                      scale_x_reverse(limits = c(1000, -100), 
                                      breaks = seq(from = 4000, to = -100, by = -100),
                                      labels = c(4000, rep('', 4),
                                                 3500,rep('', 4),
                                                 3000, rep('', 4),
                                                 2500,rep('', 4),
                                                 2000, rep('', 4),
                                                 1500,rep('', 4),
                                                 1000, rep('', 4),
                                                 500,rep('', 4),
                                                 0, ''),
                                      sec.axis=sec_axis(~., 
                                                        breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                                        labels=c(2000, rep('', 4),
                                                                 1500,rep('', 4),
                                                                 1000, rep('', 4),
                                                                 500,rep('', 4),
                                                                 0, rep('', 4),
                                                                 500,rep('', 4),
                                                                 1000, rep('', 4),
                                                                 1500,rep('', 4),
                                                                 2000)
                                      )
                      ),
                    WholeCC +
                      theme(
                        axis.line.x = element_blank(),
                        axis.ticks.x = element_blank(),
                        axis.title.x = element_blank(),
                        axis.text.x = element_blank(),
                      ),
                    movChronUnCalCCPlot,
                    heights = c(1,1,1),
                    widths = c(1,1,1),
                    ncol = 1, nrow = 3,
                    align = "v")
        }
      }
      
      #PC-TSI
      {
        #PC1
        {
          #Setup sub-group size and step
          movCorreWindowHeadCut <- 28
          movCorreWindowSize <- 25
          movCorreWindowStep <- 4
          
          #Setup significant level for moving correlation confidence interval check
          sigVal <- 0.05
          
          #Create a matrix to store the median, 95% range and 95% confidential interval
          #for sectioned correlation
          movCorreMatChronUnPC1TSI <- matrix(nrow = floor((length(UACMatchingPC1.PA$values) - movCorreWindowSize - movCorreWindowHeadCut) / movCorreWindowStep),
                                             ncol = 10,
                                             dimnames = list(c(),
                                                             c('Ind1', 'Ind2',
                                                               'Time', 'Time1', 'Time2',
                                                               'P1','M','P2',
                                                               'Sig1', 'Sig2')))
          
          movCorreMatChronUnPC1TSI[, 1] <- seq(from = 0, to = nrow(movCorreMatChronUnPC1TSI) - 1,
                                               by = 1) * movCorreWindowStep + 1
          
          movCorreMatChronUnPC1TSI[, 2] <- movCorreMatChronUnPC1TSI[, 1] + movCorreWindowSize
          
          movCorreMatChronUnPC1TSI[, 3:10] <-
            foreach(paraInd = c(1:nrow(movCorreMatChronUnPC1TSI)), .packages = c("zoo", 'geoChronR'),
                    .combine = 'rbind') %dopar% 
            {
              try(
                {
                  tempTime1 <- min(UACMatchingPC1.ae$values[movCorreMatChronUnPC1TSI[paraInd, 1]:movCorreMatChronUnPC1TSI[paraInd, 2],])
                  tempTime2 <- max(UACMatchingPC1.ae$values[movCorreMatChronUnPC1TSI[paraInd, 1]:movCorreMatChronUnPC1TSI[paraInd, 2],])
                  
                  tempTime3 <- quantile(UACMatchingPC1.ae$values[movCorreMatChronUnPC1TSI[paraInd, 1]:movCorreMatChronUnPC1TSI[paraInd, 2],], probs = 0.025)
                  tempTime4 <- quantile(UACMatchingPC1.ae$values[movCorreMatChronUnPC1TSI[paraInd, 1]:movCorreMatChronUnPC1TSI[paraInd, 2],], probs = 0.975)
                  
                  tempIndTSI <- which((TSIRaw[,1] >= tempTime1) & 
                                        (TSIRaw[,1] <= tempTime2))
                  
                  coroutPC1TSITemp <- corEns(time.1 = UACMatchingPC1.ae$values[movCorreMatChronUnPC1TSI[paraInd, 1]:movCorreMatChronUnPC1TSI[paraInd, 2], ],
                                             values.1 = UACMatchingPC1.PA$values[movCorreMatChronUnPC1TSI[paraInd, 1]:movCorreMatChronUnPC1TSI[paraInd, 2]],
                                             time.2 = t(matrix(replicate(1000,TSIRaw[tempIndTSI, 1]),nrow = 1000)),
                                             values.2 = TSIRaw[tempIndTSI, 2],
                                             bin.step = 50,
                                             max.ens = 1000,
                                             min.obs = 3,
                                             isopersistent  = TRUE,
                                             isospectral = TRUE)
                  
                  timeMedian <- median(apply(UACMatchingPC1.ae$values[movCorreMatChronUnPC1TSI[paraInd, 1]:movCorreMatChronUnPC1TSI[paraInd, 2], ], 2, median))
                  
                  tempSigCheck <- sort(coroutPC1TSITemp[[1]][[5]]) < sigVal
                  
                  tempSigp <- sort(coroutPC1TSITemp[[1]][[1]])[c(max(min(which(tempSigCheck == FALSE)) - 1, 1), 
                                                                 min(max(which(tempSigCheck == FALSE)) + 1, 
                                                                     length(sort(coroutPC1TSITemp[[1]][[1]]))
                                                                 )
                  )
                  ]
                  
                  return(c(timeMedian, tempTime3, tempTime4, 
                           coroutPC1TSITemp[[2]][[5]][1],
                           coroutPC1TSITemp[[2]][[5]][3],
                           coroutPC1TSITemp[[2]][[5]][5],
                           tempSigp))
                },
                silent = TRUE
              )
              
            }
          
          movCorreMatChronUnPC1TSI <- apply(movCorreMatChronUnPC1TSI, c(1,2), as.numeric)
        }
        
        #PC2
        {
          #Setup sub-group size and step
          movCorreWindowHeadCut <- 28
          movCorreWindowSize <- 25
          movCorreWindowStep <- 4
          
          #Setup significant level for moving correlation confidence interval check
          sigVal <- 0.05
          
          #Create a matrix to store the median, 95% range and 95% confidential interval
          #for sectioned correlation
          movCorreMatChronUnPC2TSI <- matrix(nrow = floor((length(UACMatchingPC2.PA$values) - movCorreWindowSize - movCorreWindowHeadCut) / movCorreWindowStep),
                                             ncol = 10,
                                             dimnames = list(c(),
                                                             c('Ind1', 'Ind2',
                                                               'Time', 'Time1', 'Time2',
                                                               'P1','M','P2',
                                                               'Sig1', 'Sig2')))
          
          movCorreMatChronUnPC2TSI[, 1] <- seq(from = 0, to = nrow(movCorreMatChronUnPC2TSI) - 1,
                                               by = 1) * movCorreWindowStep + 1
          
          movCorreMatChronUnPC2TSI[, 2] <- movCorreMatChronUnPC2TSI[, 1] + movCorreWindowSize
          
          movCorreMatChronUnPC2TSI[, 3:10] <-
            foreach(paraInd = c(1:nrow(movCorreMatChronUnPC2TSI)), .packages = c("zoo", 'geoChronR'),
                    .combine = 'rbind') %dopar% 
            {
              try(
                {
                  tempTime1 <- min(UACMatchingPC2.ae$values[movCorreMatChronUnPC2TSI[paraInd, 1]:movCorreMatChronUnPC2TSI[paraInd, 2],])
                  tempTime2 <- max(UACMatchingPC2.ae$values[movCorreMatChronUnPC2TSI[paraInd, 1]:movCorreMatChronUnPC2TSI[paraInd, 2],])
                  
                  tempTime3 <- quantile(UACMatchingPC2.ae$values[movCorreMatChronUnPC2TSI[paraInd, 1]:movCorreMatChronUnPC2TSI[paraInd, 2],], probs = 0.025)
                  tempTime4 <- quantile(UACMatchingPC2.ae$values[movCorreMatChronUnPC2TSI[paraInd, 1]:movCorreMatChronUnPC2TSI[paraInd, 2],], probs = 0.975)
                  
                  tempIndTSI <- which((TSIRaw[,1] >= tempTime1) & 
                                        (TSIRaw[,1] <= tempTime2))
                  
                  coroutPC2TSITemp <- corEns(time.1 = UACMatchingPC2.ae$values[movCorreMatChronUnPC2TSI[paraInd, 1]:movCorreMatChronUnPC2TSI[paraInd, 2], ],
                                             values.1 = UACMatchingPC2.PA$values[movCorreMatChronUnPC2TSI[paraInd, 1]:movCorreMatChronUnPC2TSI[paraInd, 2]],
                                             time.2 = t(matrix(replicate(1000,TSIRaw[tempIndTSI, 1]),nrow = 1000)),
                                             values.2 = TSIRaw[tempIndTSI, 2],
                                             bin.step = 50,
                                             max.ens = 1000,
                                             min.obs = 3,
                                             isopersistent  = TRUE,
                                             isospectral = TRUE)
                  
                  timeMedian <- median(apply(UACMatchingPC2.ae$values[movCorreMatChronUnPC2TSI[paraInd, 1]:movCorreMatChronUnPC2TSI[paraInd, 2], ], 2, median))
                  
                  tempSigCheck <- sort(coroutPC2TSITemp[[1]][[5]]) < sigVal
                  
                  tempSigp <- sort(coroutPC2TSITemp[[1]][[1]])[c(max(min(which(tempSigCheck == FALSE)) - 1, 1), 
                                                                 min(max(which(tempSigCheck == FALSE)) + 1, 
                                                                     length(sort(coroutPC2TSITemp[[1]][[1]]))
                                                                 )
                  )
                  ]
                  
                  return(c(timeMedian, tempTime3, tempTime4, 
                           coroutPC2TSITemp[[2]][[5]][1],
                           coroutPC2TSITemp[[2]][[5]][3],
                           coroutPC2TSITemp[[2]][[5]][5],
                           tempSigp))
                },
                silent = TRUE
              )
              
            }
          
          movCorreMatChronUnPC2TSI <- apply(movCorreMatChronUnPC2TSI, c(1,2), as.numeric)
        }
        
        #PC3
        {
          #Setup sub-group size and step
          movCorreWindowHeadCut <- 28
          movCorreWindowSize <- 25
          movCorreWindowStep <- 4
          
          #Setup significant level for moving correlation confidence interval check
          sigVal <- 0.05
          
          #Create a matrix to store the median, 95% range and 95% confidential interval
          #for sectioned correlation
          movCorreMatChronUnPC3TSI <- matrix(nrow = floor((length(UACMatchingPC3.PA$values) - movCorreWindowSize - movCorreWindowHeadCut) / movCorreWindowStep),
                                             ncol = 10,
                                             dimnames = list(c(),
                                                             c('Ind1', 'Ind2',
                                                               'Time', 'Time1', 'Time2',
                                                               'P1','M','P2',
                                                               'Sig1', 'Sig2')))
          
          movCorreMatChronUnPC3TSI[, 1] <- seq(from = 0, to = nrow(movCorreMatChronUnPC3TSI) - 1,
                                               by = 1) * movCorreWindowStep + 1
          
          movCorreMatChronUnPC3TSI[, 2] <- movCorreMatChronUnPC3TSI[, 1] + movCorreWindowSize
          
          movCorreMatChronUnPC3TSI[, 3:10] <-
            foreach(paraInd = c(1:nrow(movCorreMatChronUnPC3TSI)), .packages = c("zoo", 'geoChronR'),
                    .combine = 'rbind') %dopar% 
            {
              try(
                {
                  tempTime1 <- min(UACMatchingPC3.ae$values[movCorreMatChronUnPC3TSI[paraInd, 1]:movCorreMatChronUnPC3TSI[paraInd, 2],])
                  tempTime2 <- max(UACMatchingPC3.ae$values[movCorreMatChronUnPC3TSI[paraInd, 1]:movCorreMatChronUnPC3TSI[paraInd, 2],])
                  
                  tempTime3 <- quantile(UACMatchingPC3.ae$values[movCorreMatChronUnPC3TSI[paraInd, 1]:movCorreMatChronUnPC3TSI[paraInd, 2],], probs = 0.025)
                  tempTime4 <- quantile(UACMatchingPC3.ae$values[movCorreMatChronUnPC3TSI[paraInd, 1]:movCorreMatChronUnPC3TSI[paraInd, 2],], probs = 0.975)
                  
                  tempIndTSI <- which((TSIRaw[,1] >= tempTime1) & 
                                        (TSIRaw[,1] <= tempTime2))
                  
                  coroutPC3TSITemp <- corEns(time.1 = UACMatchingPC3.ae$values[movCorreMatChronUnPC3TSI[paraInd, 1]:movCorreMatChronUnPC3TSI[paraInd, 2], ],
                                             values.1 = UACMatchingPC3.PA$values[movCorreMatChronUnPC3TSI[paraInd, 1]:movCorreMatChronUnPC3TSI[paraInd, 2]],
                                             time.2 = t(matrix(replicate(1000,TSIRaw[tempIndTSI, 1]),nrow = 1000)),
                                             values.2 = TSIRaw[tempIndTSI, 2],
                                             bin.step = 50,
                                             max.ens = 1000,
                                             min.obs = 3,
                                             isopersistent  = TRUE,
                                             isospectral = TRUE)
                  
                  timeMedian <- median(apply(UACMatchingPC3.ae$values[movCorreMatChronUnPC3TSI[paraInd, 1]:movCorreMatChronUnPC3TSI[paraInd, 2], ], 2, median))
                  
                  tempSigCheck <- sort(coroutPC3TSITemp[[1]][[5]]) < sigVal
                  
                  tempSigp <- sort(coroutPC3TSITemp[[1]][[1]])[c(max(min(which(tempSigCheck == FALSE)) - 1, 1), 
                                                                 min(max(which(tempSigCheck == FALSE)) + 1, 
                                                                     length(sort(coroutPC3TSITemp[[1]][[1]]))
                                                                 )
                  )
                  ]
                  
                  return(c(timeMedian, tempTime3, tempTime4, 
                           coroutPC3TSITemp[[2]][[5]][1],
                           coroutPC3TSITemp[[2]][[5]][3],
                           coroutPC3TSITemp[[2]][[5]][5],
                           tempSigp))
                },
                silent = TRUE
              )
              
            }
          
          movCorreMatChronUnPC3TSI <- apply(movCorreMatChronUnPC3TSI, c(1,2), as.numeric)
        }
        
        #Plot the sectioned correlation
        {
          #PC1
          movChronUnPC1TSIPlot <- ggplot() +
            geom_point(aes(y = movCorreMatChronUnPC1TSI[, 7],
                           x = movCorreMatChronUnPC1TSI[, 3]), color = 'black',
                       size = 2) +
            geom_errorbar(aes(x = movCorreMatChronUnPC1TSI[, 3],
                              y = movCorreMatChronUnPC1TSI[, 7],
                              ymin = movCorreMatChronUnPC1TSI[, 6],
                              ymax = movCorreMatChronUnPC1TSI[, 8]),
                          fill = 'black', linewidth = 1, alpha = 0.3,
                          width = 20) +
            geom_errorbarh(aes(xmin = movCorreMatChronUnPC1TSI[, 4],
                               xmax = movCorreMatChronUnPC1TSI[, 5],
                               y = movCorreMatChronUnPC1TSI[, 7],
                               x = movCorreMatChronUnPC1TSI[, 3]
            ),
            fill = 'black', linewidth = 2, 
            alpha = 0.3) +
            geom_line(aes(y = movCorreMatChronUnPC1TSI[, 10],
                          x = movCorreMatChronUnPC1TSI[, 3]), color = 'red', linetype = 2, alpha = 0.6) +
            geom_line(aes(y = movCorreMatChronUnPC1TSI[, 9],
                          x = movCorreMatChronUnPC1TSI[, 3]), color = 'red', linetype = 2, alpha = 0.6) +
            geom_hline(yintercept = 0, color = 'grey', alpha = 0.6) +
            scale_y_continuous(breaks = seq(from = -10, to = 10, by = 0.25),
                               limits = c(0,1),
                               position = 'left') +
            scale_x_reverse(limits = c(2750, -100), 
                            breaks = seq(from = 4000, to = -100, by = -100),
                            labels = c(4000, rep('', 4),
                                       3500,rep('', 4),
                                       3000, rep('', 4),
                                       2500,rep('', 4),
                                       2000, rep('', 4),
                                       1500,rep('', 4),
                                       1000, rep('', 4),
                                       500,rep('', 4),
                                       0, ''),
                            sec.axis=sec_axis(~., 
                                              breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                              labels=c(2000, rep('', 4),
                                                       1500,rep('', 4),
                                                       1000, rep('', 4),
                                                       500,rep('', 4),
                                                       0, rep('', 4),
                                                       500,rep('', 4),
                                                       1000, rep('', 4),
                                                       1500,rep('', 4),
                                                       2000)
                            )
            )  +
            guides(y = "axis_truncated", x = "axis_truncated") +
            labs(y = 'PC1&TSI', x = 'PC3 years BP', title = '') + 
            theme_classic() + 
            theme(plot.title = element_blank(), 
                  axis.text.x = element_text(size = 10, family = 'arial'),
                  axis.title.x = element_text(size = 12, family = 'arial'),
                  axis.text.y = element_text(size = 10, family = 'arial'),
                  axis.title.y = element_text(size = 12, family = 'arial')
            )
          
          #PC2
          movChronUnPC2TSIPlot <- ggplot() +
            geom_point(aes(y = movCorreMatChronUnPC2TSI[, 7],
                           x = movCorreMatChronUnPC2TSI[, 3]), color = 'black',
                       size = 2) +
            geom_errorbar(aes(x = movCorreMatChronUnPC2TSI[, 3],
                              y = movCorreMatChronUnPC2TSI[, 7],
                              ymin = movCorreMatChronUnPC2TSI[, 6],
                              ymax = movCorreMatChronUnPC2TSI[, 8]),
                          fill = 'black', linewidth = 1, alpha = 0.3,
                          width = 20) +
            geom_errorbarh(aes(xmin = movCorreMatChronUnPC2TSI[, 4],
                               xmax = movCorreMatChronUnPC2TSI[, 5],
                               y = movCorreMatChronUnPC2TSI[, 7],
                               x = movCorreMatChronUnPC2TSI[, 3]
            ),
            fill = 'black', linewidth = 2, 
            alpha = 0.3) +
            geom_line(aes(y = movCorreMatChronUnPC2TSI[, 10],
                          x = movCorreMatChronUnPC2TSI[, 3]), color = 'red', linetype = 2, alpha = 0.6) +
            geom_line(aes(y = movCorreMatChronUnPC2TSI[, 9],
                          x = movCorreMatChronUnPC2TSI[, 3]), color = 'red', linetype = 2, alpha = 0.6) +
            geom_hline(yintercept = 0, color = 'grey', alpha = 0.6) +
            scale_y_continuous(breaks = seq(from = -10, to = 10, by = 0.25),
                               limits = c(0,1),
                               position = 'left') +
            scale_x_reverse(limits = c(2750, -100), 
                            breaks = seq(from = 4000, to = -100, by = -100),
                            labels = c(4000, rep('', 4),
                                       3500,rep('', 4),
                                       3000, rep('', 4),
                                       2500,rep('', 4),
                                       2000, rep('', 4),
                                       1500,rep('', 4),
                                       1000, rep('', 4),
                                       500,rep('', 4),
                                       0, ''),
                            sec.axis=sec_axis(~., 
                                              breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                              labels=c(2000, rep('', 4),
                                                       1500,rep('', 4),
                                                       1000, rep('', 4),
                                                       500,rep('', 4),
                                                       0, rep('', 4),
                                                       500,rep('', 4),
                                                       1000, rep('', 4),
                                                       1500,rep('', 4),
                                                       2000)
                            )
            )  +
            guides(y = "axis_truncated", x = "axis_truncated") +
            labs(y = 'PC2&TSI', x = 'PC3 years BP', title = '') + 
            theme_classic() + 
            theme(plot.title = element_blank(), 
                  axis.text.x = element_text(size = 10, family = 'arial'),
                  axis.title.x = element_text(size = 12, family = 'arial'),
                  axis.text.y = element_text(size = 10, family = 'arial'),
                  axis.title.y = element_text(size = 12, family = 'arial')
            )
          
          #PC3
          movChronUnPC3TSIPlot <- ggplot() +
            geom_point(aes(y = movCorreMatChronUnPC3TSI[, 7],
                           x = movCorreMatChronUnPC3TSI[, 3]), color = 'black',
                       size = 2) +
            geom_errorbar(aes(x = movCorreMatChronUnPC3TSI[, 3],
                              y = movCorreMatChronUnPC3TSI[, 7],
                              ymin = movCorreMatChronUnPC3TSI[, 6],
                              ymax = movCorreMatChronUnPC3TSI[, 8]),
                          fill = 'black', linewidth = 1, alpha = 0.3,
                          width = 20) +
            geom_errorbarh(aes(xmin = movCorreMatChronUnPC3TSI[, 4],
                               xmax = movCorreMatChronUnPC3TSI[, 5],
                               y = movCorreMatChronUnPC3TSI[, 7],
                               x = movCorreMatChronUnPC3TSI[, 3]
            ),
            fill = 'black', linewidth = 2, 
            alpha = 0.3) +
            geom_line(aes(y = movCorreMatChronUnPC3TSI[, 10],
                          x = movCorreMatChronUnPC3TSI[, 3]), color = 'red', linetype = 2, alpha = 0.6) +
            geom_line(aes(y = movCorreMatChronUnPC3TSI[, 9],
                          x = movCorreMatChronUnPC3TSI[, 3]), color = 'red', linetype = 2, alpha = 0.6) +
            geom_hline(yintercept = 0, color = 'grey', alpha = 0.6) +
            scale_y_continuous(breaks = seq(from = -10, to = 10, by = 0.25),
                               limits = c(0,1),
                               position = 'left') +
            scale_x_reverse(limits = c(2750, -100), 
                            breaks = seq(from = 4000, to = -100, by = -100),
                            labels = c(4000, rep('', 4),
                                       3500,rep('', 4),
                                       3000, rep('', 4),
                                       2500,rep('', 4),
                                       2000, rep('', 4),
                                       1500,rep('', 4),
                                       1000, rep('', 4),
                                       500,rep('', 4),
                                       0, ''),
                            sec.axis=sec_axis(~., 
                                              breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                              labels=c(2000, rep('', 4),
                                                       1500,rep('', 4),
                                                       1000, rep('', 4),
                                                       500,rep('', 4),
                                                       0, rep('', 4),
                                                       500,rep('', 4),
                                                       1000, rep('', 4),
                                                       1500,rep('', 4),
                                                       2000)
                            )
            )  +
            guides(y = "axis_truncated", x = "axis_truncated") +
            labs(y = 'PC3&TSI', x = 'PC3 years BP', title = '') + 
            theme_classic() + 
            theme(plot.title = element_blank(), 
                  axis.text.x = element_text(size = 10, family = 'arial'),
                  axis.title.x = element_text(size = 12, family = 'arial'),
                  axis.text.y = element_text(size = 10, family = 'arial'),
                  axis.title.y = element_text(size = 12, family = 'arial')
            )
        }
        
        #Combine plot
        {
          ggarrange(UACMatchingPC1.PA.ts.plot +
                      geom_hline(yintercept = UACPC1ChronUnMean, 
                                 linetype = 1, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACPC1ChronUnMean - UACPC1ChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACPC1ChronUnMean - 2 * UACPC1ChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACPC1ChronUnMean + UACPC1ChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACPC1ChronUnMean + 2* UACPC1ChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) ,
                    WholeTSI +
                      theme(
                        axis.line.x = element_blank(),
                        axis.ticks.x = element_blank(),
                        axis.title.x = element_blank(),
                        axis.text.x = element_blank(),
                      ),
                    movChronUnPC1TSIPlot,
                    heights = c(1,1,1),
                    widths = c(1,1,1),
                    ncol = 1, nrow = 3,
                    align = "v")
          
          ggarrange(UACMatchingPC2.PA.ts.plot +
                      geom_hline(yintercept = UACPC2ChronUnMean, 
                                 linetype = 1, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACPC2ChronUnMean - UACPC2ChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACPC2ChronUnMean - 2 * UACPC2ChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACPC2ChronUnMean + UACPC2ChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACPC2ChronUnMean + 2* UACPC2ChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) ,
                    WholeTSI +
                      theme(
                        axis.line.x = element_blank(),
                        axis.ticks.x = element_blank(),
                        axis.title.x = element_blank(),
                        axis.text.x = element_blank(),
                      ),
                    movChronUnPC2TSIPlot,
                    heights = c(1,1,1),
                    widths = c(1,1,1),
                    ncol = 1, nrow = 3,
                    align = "v")
          
          ggarrange(UACMatchingPC3.PA.ts.plot +
                      geom_hline(yintercept = UACPC3ChronUnMean, 
                                 linetype = 1, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACPC3ChronUnMean - UACPC3ChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACPC3ChronUnMean - 2 * UACPC3ChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACPC3ChronUnMean + UACPC3ChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACPC3ChronUnMean + 2* UACPC3ChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) ,
                    WholeTSI +
                      theme(
                        axis.line.x = element_blank(),
                        axis.ticks.x = element_blank(),
                        axis.title.x = element_blank(),
                        axis.text.x = element_blank(),
                      ),
                    movChronUnPC3TSIPlot,
                    heights = c(1,1,1),
                    widths = c(1,1,1),
                    ncol = 1, nrow = 3,
                    align = "v")
        }
      }
      
      #PC-CC
      {
        #PC1
        {
          #Setup sub-group size and step
          movCorreWindowHeadCut <- 28
          movCorreWindowSize <- 25
          movCorreWindowStep <- 4
          
          #Setup significant level for moving correlation confidence interval check
          sigVal <- 0.05
          
          #Create a matrix to store the median, 95% range and 95% confidential interval
          #for sectioned correlation
          movCorreMatChronUnPC1CC <- matrix(nrow = floor((length(UACMatchingPC1.PA$values) - movCorreWindowSize - movCorreWindowHeadCut) / movCorreWindowStep),
                                             ncol = 10,
                                             dimnames = list(c(),
                                                             c('Ind1', 'Ind2',
                                                               'Time', 'Time1', 'Time2',
                                                               'P1','M','P2',
                                                               'Sig1', 'Sig2')))
          
          movCorreMatChronUnPC1CC[, 1] <- seq(from = 0, to = nrow(movCorreMatChronUnPC1CC) - 1,
                                               by = 1) * movCorreWindowStep + 1
          
          movCorreMatChronUnPC1CC[, 2] <- movCorreMatChronUnPC1CC[, 1] + movCorreWindowSize
          
          movCorreMatChronUnPC1CC[, 3:10] <-
            foreach(paraInd = c(1:nrow(movCorreMatChronUnPC1CC)), .packages = c("zoo", 'geoChronR'),
                    .combine = 'rbind') %dopar% 
            {
              try(
                {
                  tempTime1 <- min(UACMatchingPC1.ae$values[movCorreMatChronUnPC1CC[paraInd, 1]:movCorreMatChronUnPC1CC[paraInd, 2],])
                  tempTime2 <- max(UACMatchingPC1.ae$values[movCorreMatChronUnPC1CC[paraInd, 1]:movCorreMatChronUnPC1CC[paraInd, 2],])
                  
                  tempTime3 <- quantile(UACMatchingPC1.ae$values[movCorreMatChronUnPC1CC[paraInd, 1]:movCorreMatChronUnPC1CC[paraInd, 2],], probs = 0.025)
                  tempTime4 <- quantile(UACMatchingPC1.ae$values[movCorreMatChronUnPC1CC[paraInd, 1]:movCorreMatChronUnPC1CC[paraInd, 2],], probs = 0.975)
                  
                  tempIndCC <- which((CCRaw[,1] >= tempTime1) & 
                                        (CCRaw[,1] <= tempTime2))
                  
                  coroutPC1CCTemp <- corEns(time.1 = UACMatchingPC1.ae$values[movCorreMatChronUnPC1CC[paraInd, 1]:movCorreMatChronUnPC1CC[paraInd, 2], ],
                                             values.1 = UACMatchingPC1.PA$values[movCorreMatChronUnPC1CC[paraInd, 1]:movCorreMatChronUnPC1CC[paraInd, 2]],
                                             time.2 = t(matrix(replicate(1000,CCRaw[tempIndCC, 1]),nrow = 1000)),
                                             values.2 = CCRaw[tempIndCC, 2],
                                             bin.step = 50,
                                             max.ens = 1000,
                                             min.obs = 3,
                                             isopersistent  = TRUE,
                                             isospectral = TRUE)
                  
                  timeMedian <- median(apply(UACMatchingPC1.ae$values[movCorreMatChronUnPC1CC[paraInd, 1]:movCorreMatChronUnPC1CC[paraInd, 2], ], 2, median))
                  
                  tempSigCheck <- sort(coroutPC1CCTemp[[1]][[5]]) < sigVal
                  
                  tempSigp <- sort(coroutPC1CCTemp[[1]][[1]])[c(max(min(which(tempSigCheck == FALSE)) - 1, 1), 
                                                                 min(max(which(tempSigCheck == FALSE)) + 1, 
                                                                     length(sort(coroutPC1CCTemp[[1]][[1]]))
                                                                 )
                  )
                  ]
                  
                  return(c(timeMedian, tempTime3, tempTime4, 
                           coroutPC1CCTemp[[2]][[5]][1],
                           coroutPC1CCTemp[[2]][[5]][3],
                           coroutPC1CCTemp[[2]][[5]][5],
                           tempSigp))
                },
                silent = TRUE
              )
              
            }
          
          movCorreMatChronUnPC1CC <- apply(movCorreMatChronUnPC1CC, c(1,2), as.numeric)
        }
        
        #PC2
        {
          #Setup sub-group size and step
          movCorreWindowHeadCut <- 28
          movCorreWindowSize <- 25
          movCorreWindowStep <- 4
          
          #Setup significant level for moving correlation confidence interval check
          sigVal <- 0.05
          
          #Create a matrix to store the median, 95% range and 95% confidential interval
          #for sectioned correlation
          movCorreMatChronUnPC2CC <- matrix(nrow = floor((length(UACMatchingPC2.PA$values) - movCorreWindowSize - movCorreWindowHeadCut) / movCorreWindowStep),
                                             ncol = 10,
                                             dimnames = list(c(),
                                                             c('Ind1', 'Ind2',
                                                               'Time', 'Time1', 'Time2',
                                                               'P1','M','P2',
                                                               'Sig1', 'Sig2')))
          
          movCorreMatChronUnPC2CC[, 1] <- seq(from = 0, to = nrow(movCorreMatChronUnPC2CC) - 1,
                                               by = 1) * movCorreWindowStep + 1
          
          movCorreMatChronUnPC2CC[, 2] <- movCorreMatChronUnPC2CC[, 1] + movCorreWindowSize
          
          movCorreMatChronUnPC2CC[, 3:10] <-
            foreach(paraInd = c(1:nrow(movCorreMatChronUnPC2CC)), .packages = c("zoo", 'geoChronR'),
                    .combine = 'rbind') %dopar% 
            {
              try(
                {
                  tempTime1 <- min(UACMatchingPC2.ae$values[movCorreMatChronUnPC2CC[paraInd, 1]:movCorreMatChronUnPC2CC[paraInd, 2],])
                  tempTime2 <- max(UACMatchingPC2.ae$values[movCorreMatChronUnPC2CC[paraInd, 1]:movCorreMatChronUnPC2CC[paraInd, 2],])
                  
                  tempTime3 <- quantile(UACMatchingPC2.ae$values[movCorreMatChronUnPC2CC[paraInd, 1]:movCorreMatChronUnPC2CC[paraInd, 2],], probs = 0.025)
                  tempTime4 <- quantile(UACMatchingPC2.ae$values[movCorreMatChronUnPC2CC[paraInd, 1]:movCorreMatChronUnPC2CC[paraInd, 2],], probs = 0.975)
                  
                  tempIndCC <- which((CCRaw[,1] >= tempTime1) & 
                                        (CCRaw[,1] <= tempTime2))
                  
                  coroutPC2CCTemp <- corEns(time.1 = UACMatchingPC2.ae$values[movCorreMatChronUnPC2CC[paraInd, 1]:movCorreMatChronUnPC2CC[paraInd, 2], ],
                                             values.1 = UACMatchingPC2.PA$values[movCorreMatChronUnPC2CC[paraInd, 1]:movCorreMatChronUnPC2CC[paraInd, 2]],
                                             time.2 = t(matrix(replicate(1000,CCRaw[tempIndCC, 1]),nrow = 1000)),
                                             values.2 = CCRaw[tempIndCC, 2],
                                             bin.step = 50,
                                             max.ens = 1000,
                                             min.obs = 3,
                                             isopersistent  = TRUE,
                                             isospectral = TRUE)
                  
                  timeMedian <- median(apply(UACMatchingPC2.ae$values[movCorreMatChronUnPC2CC[paraInd, 1]:movCorreMatChronUnPC2CC[paraInd, 2], ], 2, median))
                  
                  tempSigCheck <- sort(coroutPC2CCTemp[[1]][[5]]) < sigVal
                  
                  tempSigp <- sort(coroutPC2CCTemp[[1]][[1]])[c(max(min(which(tempSigCheck == FALSE)) - 1, 1), 
                                                                 min(max(which(tempSigCheck == FALSE)) + 1, 
                                                                     length(sort(coroutPC2CCTemp[[1]][[1]]))
                                                                 )
                  )
                  ]
                  
                  return(c(timeMedian, tempTime3, tempTime4, 
                           coroutPC2CCTemp[[2]][[5]][1],
                           coroutPC2CCTemp[[2]][[5]][3],
                           coroutPC2CCTemp[[2]][[5]][5],
                           tempSigp))
                },
                silent = TRUE
              )
              
            }
          
          movCorreMatChronUnPC2CC <- apply(movCorreMatChronUnPC2CC, c(1,2), as.numeric)
        }
        
        #PC3
        {
          #Setup sub-group size and step
          movCorreWindowHeadCut <- 28
          movCorreWindowSize <- 25
          movCorreWindowStep <- 4
          
          #Setup significant level for moving correlation confidence interval check
          sigVal <- 0.05
          
          #Create a matrix to store the median, 95% range and 95% confidential interval
          #for sectioned correlation
          movCorreMatChronUnPC3CC <- matrix(nrow = floor((length(UACMatchingPC3.PA$values) - movCorreWindowSize - movCorreWindowHeadCut) / movCorreWindowStep),
                                             ncol = 10,
                                             dimnames = list(c(),
                                                             c('Ind1', 'Ind2',
                                                               'Time', 'Time1', 'Time2',
                                                               'P1','M','P2',
                                                               'Sig1', 'Sig2')))
          
          movCorreMatChronUnPC3CC[, 1] <- seq(from = 0, to = nrow(movCorreMatChronUnPC3CC) - 1,
                                               by = 1) * movCorreWindowStep + 1
          
          movCorreMatChronUnPC3CC[, 2] <- movCorreMatChronUnPC3CC[, 1] + movCorreWindowSize
          
          movCorreMatChronUnPC3CC[, 3:10] <-
            foreach(paraInd = c(1:nrow(movCorreMatChronUnPC3CC)), .packages = c("zoo", 'geoChronR'),
                    .combine = 'rbind') %dopar% 
            {
              try(
                {
                  tempTime1 <- min(UACMatchingPC3.ae$values[movCorreMatChronUnPC3CC[paraInd, 1]:movCorreMatChronUnPC3CC[paraInd, 2],])
                  tempTime2 <- max(UACMatchingPC3.ae$values[movCorreMatChronUnPC3CC[paraInd, 1]:movCorreMatChronUnPC3CC[paraInd, 2],])
                  
                  tempTime3 <- quantile(UACMatchingPC3.ae$values[movCorreMatChronUnPC3CC[paraInd, 1]:movCorreMatChronUnPC3CC[paraInd, 2],], probs = 0.025)
                  tempTime4 <- quantile(UACMatchingPC3.ae$values[movCorreMatChronUnPC3CC[paraInd, 1]:movCorreMatChronUnPC3CC[paraInd, 2],], probs = 0.975)
                  
                  tempIndCC <- which((CCRaw[,1] >= tempTime1) & 
                                        (CCRaw[,1] <= tempTime2))
                  
                  coroutPC3CCTemp <- corEns(time.1 = UACMatchingPC3.ae$values[movCorreMatChronUnPC3CC[paraInd, 1]:movCorreMatChronUnPC3CC[paraInd, 2], ],
                                             values.1 = UACMatchingPC3.PA$values[movCorreMatChronUnPC3CC[paraInd, 1]:movCorreMatChronUnPC3CC[paraInd, 2]],
                                             time.2 = t(matrix(replicate(1000,CCRaw[tempIndCC, 1]),nrow = 1000)),
                                             values.2 = CCRaw[tempIndCC, 2],
                                             bin.step = 50,
                                             max.ens = 1000,
                                             min.obs = 3,
                                             isopersistent  = TRUE,
                                             isospectral = TRUE)
                  
                  timeMedian <- median(apply(UACMatchingPC3.ae$values[movCorreMatChronUnPC3CC[paraInd, 1]:movCorreMatChronUnPC3CC[paraInd, 2], ], 2, median))
                  
                  tempSigCheck <- sort(coroutPC3CCTemp[[1]][[5]]) < sigVal
                  
                  tempSigp <- sort(coroutPC3CCTemp[[1]][[1]])[c(max(min(which(tempSigCheck == FALSE)) - 1, 1), 
                                                                 min(max(which(tempSigCheck == FALSE)) + 1, 
                                                                     length(sort(coroutPC3CCTemp[[1]][[1]]))
                                                                 )
                  )
                  ]
                  
                  return(c(timeMedian, tempTime3, tempTime4, 
                           coroutPC3CCTemp[[2]][[5]][1],
                           coroutPC3CCTemp[[2]][[5]][3],
                           coroutPC3CCTemp[[2]][[5]][5],
                           tempSigp))
                },
                silent = TRUE
              )
              
            }
          
          movCorreMatChronUnPC3CC <- apply(movCorreMatChronUnPC3CC, c(1,2), as.numeric)
        }
        
        #Plot the sectioned correlation
        {
          #PC1
          movChronUnPC1CCPlot <- ggplot() +
            geom_point(aes(y = movCorreMatChronUnPC1CC[, 7],
                           x = movCorreMatChronUnPC1CC[, 3]), color = 'black',
                       size = 2) +
            geom_errorbar(aes(x = movCorreMatChronUnPC1CC[, 3],
                              y = movCorreMatChronUnPC1CC[, 7],
                              ymin = movCorreMatChronUnPC1CC[, 6],
                              ymax = movCorreMatChronUnPC1CC[, 8]),
                          fill = 'black', linewidth = 1, alpha = 0.3,
                          width = 20) +
            geom_errorbarh(aes(xmin = movCorreMatChronUnPC1CC[, 4],
                               xmax = movCorreMatChronUnPC1CC[, 5],
                               y = movCorreMatChronUnPC1CC[, 7],
                               x = movCorreMatChronUnPC1CC[, 3]
            ),
            fill = 'black', linewidth = 2, 
            alpha = 0.3) +
            geom_line(aes(y = movCorreMatChronUnPC1CC[, 10],
                          x = movCorreMatChronUnPC1CC[, 3]), color = 'red', linetype = 2, alpha = 0.6) +
            geom_line(aes(y = movCorreMatChronUnPC1CC[, 9],
                          x = movCorreMatChronUnPC1CC[, 3]), color = 'red', linetype = 2, alpha = 0.6) +
            geom_hline(yintercept = 0, color = 'grey', alpha = 0.6) +
            scale_y_continuous(breaks = seq(from = -10, to = 10, by = 0.25),
                               limits = c(0,1),
                               position = 'left') +
            scale_x_reverse(limits = c(1000, -100), 
                            breaks = seq(from = 4000, to = -100, by = -100),
                            labels = c(4000, rep('', 4),
                                       3500,rep('', 4),
                                       3000, rep('', 4),
                                       2500,rep('', 4),
                                       2000, rep('', 4),
                                       1500,rep('', 4),
                                       1000, rep('', 4),
                                       500,rep('', 4),
                                       0, ''),
                            sec.axis=sec_axis(~., 
                                              breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                              labels=c(2000, rep('', 4),
                                                       1500,rep('', 4),
                                                       1000, rep('', 4),
                                                       500,rep('', 4),
                                                       0, rep('', 4),
                                                       500,rep('', 4),
                                                       1000, rep('', 4),
                                                       1500,rep('', 4),
                                                       2000)
                            )
            )  +
            guides(y = "axis_truncated", x = "axis_truncated") +
            labs(y = 'PC1&CC', x = 'PC3 years BP', title = '') + 
            theme_classic() + 
            theme(plot.title = element_blank(), 
                  axis.text.x = element_text(size = 10, family = 'arial'),
                  axis.title.x = element_text(size = 12, family = 'arial'),
                  axis.text.y = element_text(size = 10, family = 'arial'),
                  axis.title.y = element_text(size = 12, family = 'arial')
            )
          
          #PC2
          movChronUnPC2CCPlot <- ggplot() +
            geom_point(aes(y = movCorreMatChronUnPC2CC[, 7],
                           x = movCorreMatChronUnPC2CC[, 3]), color = 'black',
                       size = 2) +
            geom_errorbar(aes(x = movCorreMatChronUnPC2CC[, 3],
                              y = movCorreMatChronUnPC2CC[, 7],
                              ymin = movCorreMatChronUnPC2CC[, 6],
                              ymax = movCorreMatChronUnPC2CC[, 8]),
                          fill = 'black', linewidth = 1, alpha = 0.3,
                          width = 20) +
            geom_errorbarh(aes(xmin = movCorreMatChronUnPC2CC[, 4],
                               xmax = movCorreMatChronUnPC2CC[, 5],
                               y = movCorreMatChronUnPC2CC[, 7],
                               x = movCorreMatChronUnPC2CC[, 3]
            ),
            fill = 'black', linewidth = 2, 
            alpha = 0.3) +
            geom_line(aes(y = movCorreMatChronUnPC2CC[, 10],
                          x = movCorreMatChronUnPC2CC[, 3]), color = 'red', linetype = 2, alpha = 0.6) +
            geom_line(aes(y = movCorreMatChronUnPC2CC[, 9],
                          x = movCorreMatChronUnPC2CC[, 3]), color = 'red', linetype = 2, alpha = 0.6) +
            geom_hline(yintercept = 0, color = 'grey', alpha = 0.6) +
            scale_y_continuous(breaks = seq(from = -10, to = 10, by = 0.25),
                               limits = c(0,1),
                               position = 'left') +
            scale_x_reverse(limits = c(1000, -100), 
                            breaks = seq(from = 4000, to = -100, by = -100),
                            labels = c(4000, rep('', 4),
                                       3500,rep('', 4),
                                       3000, rep('', 4),
                                       2500,rep('', 4),
                                       2000, rep('', 4),
                                       1500,rep('', 4),
                                       1000, rep('', 4),
                                       500,rep('', 4),
                                       0, ''),
                            sec.axis=sec_axis(~., 
                                              breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                              labels=c(2000, rep('', 4),
                                                       1500,rep('', 4),
                                                       1000, rep('', 4),
                                                       500,rep('', 4),
                                                       0, rep('', 4),
                                                       500,rep('', 4),
                                                       1000, rep('', 4),
                                                       1500,rep('', 4),
                                                       2000)
                            )
            )  +
            guides(y = "axis_truncated", x = "axis_truncated") +
            labs(y = 'PC2&CC', x = 'PC3 years BP', title = '') + 
            theme_classic() + 
            theme(plot.title = element_blank(), 
                  axis.text.x = element_text(size = 10, family = 'arial'),
                  axis.title.x = element_text(size = 12, family = 'arial'),
                  axis.text.y = element_text(size = 10, family = 'arial'),
                  axis.title.y = element_text(size = 12, family = 'arial')
            )
          
          #PC3
          movChronUnPC3CCPlot <- ggplot() +
            geom_point(aes(y = movCorreMatChronUnPC3CC[, 7],
                           x = movCorreMatChronUnPC3CC[, 3]), color = 'black',
                       size = 2) +
            geom_errorbar(aes(x = movCorreMatChronUnPC3CC[, 3],
                              y = movCorreMatChronUnPC3CC[, 7],
                              ymin = movCorreMatChronUnPC3CC[, 6],
                              ymax = movCorreMatChronUnPC3CC[, 8]),
                          fill = 'black', linewidth = 1, alpha = 0.3,
                          width = 20) +
            geom_errorbarh(aes(xmin = movCorreMatChronUnPC3CC[, 4],
                               xmax = movCorreMatChronUnPC3CC[, 5],
                               y = movCorreMatChronUnPC3CC[, 7],
                               x = movCorreMatChronUnPC3CC[, 3]
            ),
            fill = 'black', linewidth = 2, 
            alpha = 0.3) +
            geom_line(aes(y = movCorreMatChronUnPC3CC[, 10],
                          x = movCorreMatChronUnPC3CC[, 3]), color = 'red', linetype = 2, alpha = 0.6) +
            geom_line(aes(y = movCorreMatChronUnPC3CC[, 9],
                          x = movCorreMatChronUnPC3CC[, 3]), color = 'red', linetype = 2, alpha = 0.6) +
            geom_hline(yintercept = 0, color = 'grey', alpha = 0.6) +
            scale_y_continuous(breaks = seq(from = -10, to = 10, by = 0.25),
                               limits = c(0,1),
                               position = 'left') +
            scale_x_reverse(limits = c(1000, -100), 
                            breaks = seq(from = 4000, to = -100, by = -100),
                            labels = c(4000, rep('', 4),
                                       3500,rep('', 4),
                                       3000, rep('', 4),
                                       2500,rep('', 4),
                                       2000, rep('', 4),
                                       1500,rep('', 4),
                                       1000, rep('', 4),
                                       500,rep('', 4),
                                       0, ''),
                            sec.axis=sec_axis(~., 
                                              breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                              labels=c(2000, rep('', 4),
                                                       1500,rep('', 4),
                                                       1000, rep('', 4),
                                                       500,rep('', 4),
                                                       0, rep('', 4),
                                                       500,rep('', 4),
                                                       1000, rep('', 4),
                                                       1500,rep('', 4),
                                                       2000)
                            )
            )  +
            guides(y = "axis_truncated", x = "axis_truncated") +
            labs(y = 'PC3&CC', x = 'PC3 years BP', title = '') + 
            theme_classic() + 
            theme(plot.title = element_blank(), 
                  axis.text.x = element_text(size = 10, family = 'arial'),
                  axis.title.x = element_text(size = 12, family = 'arial'),
                  axis.text.y = element_text(size = 10, family = 'arial'),
                  axis.title.y = element_text(size = 12, family = 'arial')
            )
        }
        
        #Combine plot
        {
          ggarrange(UACMatchingPC1.PA.ts.plot +
                      geom_hline(yintercept = UACPC1ChronUnMean, 
                                 linetype = 1, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACPC1ChronUnMean - UACPC1ChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACPC1ChronUnMean - 2 * UACPC1ChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACPC1ChronUnMean + UACPC1ChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACPC1ChronUnMean + 2* UACPC1ChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) +
                      scale_x_reverse(limits = c(1000, -100), 
                                      breaks = seq(from = 4000, to = -100, by = -100),
                                      labels = c(4000, rep('', 4),
                                                 3500,rep('', 4),
                                                 3000, rep('', 4),
                                                 2500,rep('', 4),
                                                 2000, rep('', 4),
                                                 1500,rep('', 4),
                                                 1000, rep('', 4),
                                                 500,rep('', 4),
                                                 0, ''),
                                      sec.axis=sec_axis(~., 
                                                        breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                                        labels=c(2000, rep('', 4),
                                                                 1500,rep('', 4),
                                                                 1000, rep('', 4),
                                                                 500,rep('', 4),
                                                                 0, rep('', 4),
                                                                 500,rep('', 4),
                                                                 1000, rep('', 4),
                                                                 1500,rep('', 4),
                                                                 2000)
                                      )
                      ),
                    WholeCC +
                      theme(
                        axis.line.x = element_blank(),
                        axis.ticks.x = element_blank(),
                        axis.title.x = element_blank(),
                        axis.text.x = element_blank(),
                      ),
                    movChronUnPC1CCPlot,
                    heights = c(1,1,1),
                    widths = c(1,1,1),
                    ncol = 1, nrow = 3,
                    align = "v")
          
          ggarrange(UACMatchingPC2.PA.ts.plot +
                      geom_hline(yintercept = UACPC2ChronUnMean, 
                                 linetype = 1, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACPC2ChronUnMean - UACPC2ChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACPC2ChronUnMean - 2 * UACPC2ChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACPC2ChronUnMean + UACPC2ChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACPC2ChronUnMean + 2* UACPC2ChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) +
                      scale_x_reverse(limits = c(1000, -100), 
                                      breaks = seq(from = 4000, to = -100, by = -100),
                                      labels = c(4000, rep('', 4),
                                                 3500,rep('', 4),
                                                 3000, rep('', 4),
                                                 2500,rep('', 4),
                                                 2000, rep('', 4),
                                                 1500,rep('', 4),
                                                 1000, rep('', 4),
                                                 500,rep('', 4),
                                                 0, ''),
                                      sec.axis=sec_axis(~., 
                                                        breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                                        labels=c(2000, rep('', 4),
                                                                 1500,rep('', 4),
                                                                 1000, rep('', 4),
                                                                 500,rep('', 4),
                                                                 0, rep('', 4),
                                                                 500,rep('', 4),
                                                                 1000, rep('', 4),
                                                                 1500,rep('', 4),
                                                                 2000)
                                      )
                      ),
                    WholeCC +
                      theme(
                        axis.line.x = element_blank(),
                        axis.ticks.x = element_blank(),
                        axis.title.x = element_blank(),
                        axis.text.x = element_blank(),
                      ),
                    movChronUnPC2CCPlot,
                    heights = c(1,1,1),
                    widths = c(1,1,1),
                    ncol = 1, nrow = 3,
                    align = "v")
          
          ggarrange(UACMatchingPC3.PA.ts.plot +
                      geom_hline(yintercept = UACPC3ChronUnMean, 
                                 linetype = 1, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACPC3ChronUnMean - UACPC3ChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACPC3ChronUnMean - 2 * UACPC3ChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACPC3ChronUnMean + UACPC3ChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) +
                      geom_hline(yintercept = UACPC3ChronUnMean + 2* UACPC3ChronUnSD, 
                                 linetype = 2, 
                                 color = 'blue', alpha = 0.5) +
                      scale_x_reverse(limits = c(1000, -100), 
                                      breaks = seq(from = 4000, to = -100, by = -100),
                                      labels = c(4000, rep('', 4),
                                                 3500,rep('', 4),
                                                 3000, rep('', 4),
                                                 2500,rep('', 4),
                                                 2000, rep('', 4),
                                                 1500,rep('', 4),
                                                 1000, rep('', 4),
                                                 500,rep('', 4),
                                                 0, ''),
                                      sec.axis=sec_axis(~., 
                                                        breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                                        labels=c(2000, rep('', 4),
                                                                 1500,rep('', 4),
                                                                 1000, rep('', 4),
                                                                 500,rep('', 4),
                                                                 0, rep('', 4),
                                                                 500,rep('', 4),
                                                                 1000, rep('', 4),
                                                                 1500,rep('', 4),
                                                                 2000)
                                      )
                      ),
                    WholeCC +
                      theme(
                        axis.line.x = element_blank(),
                        axis.ticks.x = element_blank(),
                        axis.title.x = element_blank(),
                        axis.text.x = element_blank(),
                      ),
                    movChronUnPC3CCPlot,
                    heights = c(1,1,1),
                    widths = c(1,1,1),
                    ncol = 1, nrow = 3,
                    align = "v")
        }
      }
      
      stopCluster(c16)
    }
  }
  
  #Reworked UAC and other climatic records
  {
    #Sohagnum
    ggarrange(radioDatesPlot, 
              histoEvents +
                theme(
                  axis.line.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                ),
              UACMatchingSph.PA.ts.plot +
                geom_hline(yintercept = 0, 
                           linetype = 1, 
                           color = 'grey', alpha = 0.9) +
                theme(
                  axis.line.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                ),
              WholeSSI +
                theme(
                  axis.line.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                ),
              WholeSI +
                theme(
                  axis.line.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                ),
              WholeIRD +
                theme(
                  axis.line.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                ),
              WholetempAnoma,
              heights = c(0.5,0.3,1,1,1,1,1.2),
              widths = c(1,1,1,1,1,1,1),
              ncol = 1, nrow = 7,
              align = "v")
    
    #Alnus
    ggarrange(radioDatesPlot, 
              histoEvents +
                theme(
                  axis.line.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                ),
              UACMatchingAln.PA.ts.plot +
                geom_hline(yintercept = 0, 
                           linetype = 1, 
                           color = 'grey', alpha = 0.9) +
                theme(
                  axis.line.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                ),
              WholeSSI +
                theme(
                  axis.line.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                ),
              WholeSI +
                theme(
                  axis.line.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                ),
              WholeIRD +
                theme(
                  axis.line.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                ),
              WholetempAnoma,
              heights = c(0.5,0.3,1,1,1,1,1.2),
              widths = c(1,1,1,1,1,1,1),
              ncol = 1, nrow = 7,
              align = "v")
    
    #Calluna
    ggarrange(radioDatesPlot, 
              histoEvents +
                theme(
                  axis.line.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                ),
              UACMatchingCal.PA.ts.plot +
                geom_hline(yintercept = 0, 
                           linetype = 1, 
                           color = 'grey', alpha = 0.9) +
                theme(
                  axis.line.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                ),
              WholeSSI +
                theme(
                  axis.line.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                ),
              WholeSI +
                theme(
                  axis.line.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                ),
              WholeIRD +
                theme(
                  axis.line.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                ),
              WholetempAnoma,
              heights = c(0.5,0.3,1,1,1,1,1.2),
              widths = c(1,1,1,1,1,1,1),
              ncol = 1, nrow = 7,
              align = "v")
    
    ggarrange(UACMatchingSph.PA.ts.plot +
                geom_hline(yintercept = 0, 
                           linetype = 1, 
                           color = 'grey', alpha = 0.9), 
              UACMatchingAln.PA.ts.plot +
                geom_hline(yintercept = 0, 
                           linetype = 1, 
                           color = 'grey', alpha = 0.9) +
                theme(
                  axis.line.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                ),
              UACMatchingCal.PA.ts.plot +
                geom_hline(yintercept = 0, 
                           linetype = 1, 
                           color = 'grey', alpha = 0.9) +
                theme(
                  axis.line.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                ),
              HumiMatchingSph.PA.ts.plot +
                geom_hline(yintercept = 0, 
                           linetype = 1, 
                           color = 'grey', alpha = 0.9) +
                theme(
                  axis.line.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                ),
              WholetempAnoma,
              heights = c(1.4,1,1,1,1.4),
              widths = c(1,1,1,1,1),
              ncol = 1, nrow = 5,
              align = "v")
  }
}

#Figure 7: PCA and other multivariate analysis
{
  #PCA
  {
    #Data treatment
    {
      #Trim the UAC signals of three taxa, to arrange them into the same time sequence
      commonTimeSeq <- sort(unique(c(UACSph[,2], UACAln[,2], UACCal[,2])))
      commonBatch <- sort(unique(c(UACSph[,1], UACAln[,1], UACCal[,1])))
      
      #Create a matrix for PCA analysis, with UAC signals of three taxa
      pcaMatAllTaxa <- matrix(nrow = length(commonTimeSeq),
                              ncol = 8,
                              dimnames = list(
                                c(), 
                                c('Depth', 'Cal_Years_BP',
                                  'Sphagnum', 'Sphagnum_SD',
                                  'Alnus', 'Alnus_SD',
                                  'Calluna', 'Calluna_SD')
                              )
      )
      
      #Fill in the matrix
      pcaMatAllTaxa[, 1] <- commonBatch * 2
      pcaMatAllTaxa[, 2] <- commonTimeSeq
      
      for (loopI in 1:nrow(pcaMatAllTaxa)) {
        
        
        pcaMatAllTaxa[loopI, 3] <- approx(x = UACSph[,2], y = UACSphMean,
                                          xout = pcaMatAllTaxa[loopI, 2])[[2]]
        pcaMatAllTaxa[loopI, 4] <- approx(x = UACSph[,2], y = UACSphSD,
                                          xout = pcaMatAllTaxa[loopI, 2])[[2]]
        
        pcaMatAllTaxa[loopI, 5] <- approx(x = UACAln[,2], y = UACAlnMean,
                                          xout = pcaMatAllTaxa[loopI, 2])[[2]]
        pcaMatAllTaxa[loopI, 6] <- approx(x = UACAln[,2], y = UACAlnSD,
                                          xout = pcaMatAllTaxa[loopI, 2])[[2]]
        
        pcaMatAllTaxa[loopI, 7] <- approx(x = UACCal[,2], y = UACCalMean,
                                          xout = pcaMatAllTaxa[loopI, 2])[[2]]
        pcaMatAllTaxa[loopI, 8] <- approx(x = UACCal[,2], y = UACCalSD,
                                          xout = pcaMatAllTaxa[loopI, 2])[[2]]
      }
      
      #Remove rows with NA
      pcaMatAllTaxa <- pcaMatAllTaxa[which(apply(pcaMatAllTaxa, 1, all)),-1]
      
      # Ordination
      treat_pca <- ordr::ordinate(as.data.frame(pcaMatAllTaxa), cols=colnames(pcaMatAllTaxa)[c(2,4,6)], model= ~ prcomp(.), argument=c('Cal_Years_BP'))
      
      pcaMatAllTaxa %>%
        subset(select = colnames(pcaMatAllTaxa)[c(2,4,6)]) %>%
        as.matrix() %>%
        print() -> treat_data_arrows
      
      colnames(treat_data_arrows) <- colnames(pcaMatAllTaxa)[c(2,4,6)]
      
      lm(treat_data_arrows ~ get_rows(treat_pca)) %>%
        as_tbl_ord() %>%
        augment_ord() %>%
        print() -> treat_data_arrows_plot
    }
    
    #PCA plot - three components
    {
      #1 vs 2
      ordr::ggbiplot(treat_pca, sec.axes="cols", scale.factor=1,
                     mapping = aes(x = 1, y = 2),) +
        # Add points and colour via treatment
        geom_rows_point(alpha=0.65) +
        # Add ellipses - 95% confidence
        # geom_mark_ellipse(aes(label = Unit, group = Unit, color = Unit), 
        #                   stat = 'rows_ellipse',
        #                   size=1, level=.95,
        #                   label.buffer = unit(0.5, 'mm'),
        #                   con.cap = 0) +
        # Add arrows
        geom_cols_vector(data = treat_data_arrows_plot,color="navy", size=1) +
        # Add arrow labels
        geom_cols_text_radiate(data = treat_data_arrows_plot, aes(label= name), size=6, color="navy") +
        # Fix axes/ensure both x and y are the same
        # expand_limits(y=c(-0.75, 0.75)) +
        # Minimal theme +
        theme_minimal(base_size=16) +
        # Add grid
        geom_hline(yintercept=0, linetype="dashed", size=1) +
        geom_vline(xintercept=0, linetype="dashed", size=1) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              text=element_text(family="arial"),
              legend.text = element_text(size = 20, family="arial")
        )
      
      #1 vs 3
      ordr::ggbiplot(treat_pca, sec.axes="cols", scale.factor=1,
                     mapping = aes(x = 1, y = 3),) +
        # Add points and colour via treatment
        geom_rows_point(alpha=0.65) +
        # Add ellipses - 95% confidence
        # geom_mark_ellipse(aes(label = Unit, group = Unit, color = Unit), 
        #                   stat = 'rows_ellipse',
        #                   size=1, level=.95,
        #                   label.buffer = unit(0.5, 'mm'),
        #                   con.cap = 0) +
        # Add arrows
        geom_cols_vector(data = treat_data_arrows_plot,color="navy", size=1) +
        # Add arrow labels
        geom_cols_text_radiate(data = treat_data_arrows_plot, aes(label= name), size=6, color="navy") +
        # Fix axes/ensure both x and y are the same
        # expand_limits(y=c(-0.75, 0.75)) +
        # Minimal theme +
        theme_minimal(base_size=16) +
        # Add grid
        geom_hline(yintercept=0, linetype="dashed", size=1) +
        geom_vline(xintercept=0, linetype="dashed", size=1) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              text=element_text(family="arial"),
              legend.text = element_text(size = 20, family="arial")
        )
      
      #2 vd 3
      ordr::ggbiplot(treat_pca, sec.axes="cols", scale.factor=1,
                     mapping = aes(x = 2, y = 3),) +
        # Add points and colour via treatment
        geom_rows_point(alpha=0.65) +
        # Add ellipses - 95% confidence
        # geom_mark_ellipse(aes(label = Unit, group = Unit, color = Unit), 
        #                   stat = 'rows_ellipse',
        #                   size=1, level=.95,
        #                   label.buffer = unit(0.5, 'mm'),
        #                   con.cap = 0) +
        # Add arrows
        geom_cols_vector(data = treat_data_arrows_plot,color="navy", size=1) +
        # Add arrow labels
        geom_cols_text_radiate(data = treat_data_arrows_plot, aes(label= name), size=6, color="navy") +
        # Fix axes/ensure both x and y are the same
        # expand_limits(y=c(-0.75, 0.75)) +
        # Minimal theme +
        theme_minimal(base_size=16) +
        # Add grid
        geom_hline(yintercept=0, linetype="dashed", size=1) +
        geom_vline(xintercept=0, linetype="dashed", size=1) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              text=element_text(family="arial"),
              legend.text = element_text(size = 20, family="arial")
        )
    }
    
    #PC loadings for three components
    {
      #PC1
      {
        treat_pca_rotation_PC1 <-cbind(rownames(treat_pca$rotation), 
                                       treat_pca$rotation[, 1])
        colnames(treat_pca_rotation_PC1) <- c('Element', 'Loadings')
        
        #Convert the matrix to data frame
        treat_pca_rotation_PC1 <- as.data.frame(treat_pca_rotation_PC1)
        treat_pca_rotation_PC1$Loadings <- as.numeric( treat_pca_rotation_PC1$Loadings)
        
        #Plot the colomn plot
        PCAUACPC1 <-  ggplot(data = (treat_pca_rotation_PC1), aes(x = Element, y = Loadings)) +
          geom_bar(stat="identity") +
          geom_text(aes(label=round(Loadings, 4)), size = 3.5,
                    vjust = sign(treat_pca_rotation_PC1$Loadings) * -0.5)+
          # Minimal theme +
          theme_classic(base_size=16) +
          geom_hline(yintercept=0, linetype="solid", size=1) +
          scale_y_continuous(trans= ssqrt_trans) +
          ylab('PC1 loadings') +
          # Add grid
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                text=element_text(family="arial"),
                legend.text = element_text(size = 20, family="arial"),
                legend.position="none",
                # axis.title.y = element_blank(),
                # axis.text.y  = element_blank(),
                # axis.ticks.y = element_blank(),
                # axis.ticks.x = element_blank(),
                axis.title.x  = element_blank()
          )
      }
      
      #PC2
      {
        treat_pca_rotation_PC2 <-cbind(rownames(treat_pca$rotation), 
                                       treat_pca$rotation[, 2])
        colnames(treat_pca_rotation_PC2) <- c('Element', 'Loadings')
        
        #Convert the matrix to data frame
        treat_pca_rotation_PC2 <- as.data.frame(treat_pca_rotation_PC2)
        treat_pca_rotation_PC2$Loadings <- as.numeric( treat_pca_rotation_PC2$Loadings)
        
        #Plot the colomn plot
        PCAUACPC2 <-  ggplot(data = (treat_pca_rotation_PC2), aes(x = Element, y = Loadings)) +
          geom_bar(stat="identity") +
          geom_text(aes(label=round(Loadings, 4)), size = 3.5,
                    vjust = sign(treat_pca_rotation_PC2$Loadings) * -0.5)+
          # Minimal theme +
          theme_classic(base_size=16) +
          geom_hline(yintercept=0, linetype="solid", size=1) +
          scale_y_continuous(trans= ssqrt_trans) +
          ylab('PC2 loadings') +
          # Add grid
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                text=element_text(family="arial"),
                legend.text = element_text(size = 20, family="arial"),
                legend.position="none",
                # axis.title.y = element_blank(),
                # axis.text.y  = element_blank(),
                # axis.ticks.y = element_blank(),
                # axis.ticks.x = element_blank(),
                # axis.title.x  = element_blank()
          )
      }
      
      #PC3
      {
        treat_pca_rotation_PC3 <-cbind(rownames(treat_pca$rotation), 
                                       treat_pca$rotation[, 3])
        colnames(treat_pca_rotation_PC3) <- c('Element', 'Loadings')
        
        #Convert the matrix to data frame
        treat_pca_rotation_PC3 <- as.data.frame(treat_pca_rotation_PC3)
        treat_pca_rotation_PC3$Loadings <- as.numeric( treat_pca_rotation_PC3$Loadings)
        
        #Plot the colomn plot
        PCAUACPC3 <-  ggplot(data = (treat_pca_rotation_PC3), aes(x = Element, y = Loadings)) +
          geom_bar(stat="identity") +
          geom_text(aes(label=round(Loadings, 4)), size = 3.5,
                    vjust = sign(treat_pca_rotation_PC3$Loadings) * -0.5)+
          # Minimal theme +
          theme_classic(base_size=16) +
          geom_hline(yintercept=0, linetype="solid", size=1) +
          scale_y_continuous(trans= ssqrt_trans) +
          ylab('PC3 loadings') +
          # Add grid
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                text=element_text(family="arial"),
                legend.text = element_text(size = 20, family="arial"),
                legend.position="none",
                # axis.title.y = element_blank(),
                # axis.text.y  = element_blank(),
                # axis.ticks.y = element_blank(),
                # axis.ticks.x = element_blank(),
                # axis.title.x  = element_blank()
          )
      }
      
      #Scores comparison
      {
        ggarrange(PCAUACPC1,
                  PCAUACPC2,
                  PCAUACPC3,
                  align = 'v',
                  nrow = 3,
                  ncol = 1)
      }
    }
    
    #PCA scores
    {
      #Retrieve the PC scores of three components
      pcaMatAllTaxa <- cbind(pcaMatAllTaxa, treat_pca$x)
      
      #Smooth the scores with Gaussian model
      PCAScoresPC1Gau <- smth.gaussian(pcaMatAllTaxa[, 8], window = 5)
      PCAScoresPC2Gau <- smth.gaussian(pcaMatAllTaxa[, 9], window = 5)
      PCAScoresPC3Gau <- smth.gaussian(pcaMatAllTaxa[, 10], window = 5)
      
      PCAScoresPC1GauMean <- mean(PCAScoresPC1Gau, na.rm = TRUE)
      PCAScoresPC1GauSD <- sd(PCAScoresPC1Gau, na.rm = TRUE)
      
      PCAScoresPC2GauMean <- mean(PCAScoresPC2Gau, na.rm = TRUE)
      PCAScoresPC2GauSD <- sd(PCAScoresPC2Gau, na.rm = TRUE)
      
      PCAScoresPC3GauMean <- mean(PCAScoresPC3Gau, na.rm = TRUE)
      PCAScoresPC3GauSD <- sd(PCAScoresPC3Gau, na.rm = TRUE)
      
    }
    
    #Compare PC scores with TSI
    {
      #Plot PC scores
      {
        #PC1
        PCAScoresPC1 <- ggplot() +
          geom_point(
            aes(x = pcaMatAllTaxa[, 1],
                y = pcaMatAllTaxa[, 8]), 
            size = 2, shape = 21, fill = 'red', alpha = 0.4) + 
          geom_line(aes(y = PCAScoresPC1Gau,
                        x = pcaMatAllTaxa[, 1]), 
                    color = 'black',
                    alpha = 1,
                    linewidth = 1.5) +
          scale_y_continuous(breaks = seq(from = -10, to = 10, by = 0.5),
                             limits = range(pcaMatAllTaxa[, 8]),
                             position = 'left') +
          scale_x_reverse(limits = c(2750, -100), 
                          breaks = seq(from = 4000, to = -100, by = -100),
                          labels = c(4000, rep('', 4),
                                     3500,rep('', 4),
                                     3000, rep('', 4),
                                     2500,rep('', 4),
                                     2000, rep('', 4),
                                     1500,rep('', 4),
                                     1000, rep('', 4),
                                     500,rep('', 4),
                                     0, ''),
                          sec.axis=sec_axis(~., 
                                            breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                            labels=c(2000, rep('', 4),
                                                     1500,rep('', 4),
                                                     1000, rep('', 4),
                                                     500,rep('', 4),
                                                     0, rep('', 4),
                                                     500,rep('', 4),
                                                     1000, rep('', 4),
                                                     1500,rep('', 4),
                                                     2000)
                          )) +
          guides(y = "axis_truncated") +
          labs(y = 'PC1', x = 'Cal years BP', title = '') + 
          theme_classic() + 
          theme(plot.title = element_blank(), 
                axis.text.y = element_text(size = 10, family = 'arial'),
                axis.title.y = element_text(size = 12, family = 'arial'),
                axis.text.x = element_text(size = 10, family = 'arial'),
                axis.title.x = element_text(size = 12, family = 'arial')
          )
        
        #PC2
        PCAScoresPC2 <- ggplot() +
          geom_point(
            aes(x = pcaMatAllTaxa[, 1],
                y = pcaMatAllTaxa[, 9]), 
            size = 2, shape = 21, fill = 'red', alpha = 0.4) + 
          geom_line(aes(y = PCAScoresPC2Gau,
                        x = pcaMatAllTaxa[, 1]), 
                    color = 'black',
                    alpha = 1,
                    linewidth = 1.5) +
          scale_y_continuous(breaks = seq(from = -10, to = 10, by = 0.5),
                             limits = range(pcaMatAllTaxa[, 9]),
                             position = 'left') +
          scale_x_reverse(limits = c(2750, -100), 
                          breaks = seq(from = 4000, to = -100, by = -100),
                          labels = c(4000, rep('', 4),
                                     3500,rep('', 4),
                                     3000, rep('', 4),
                                     2500,rep('', 4),
                                     2000, rep('', 4),
                                     1500,rep('', 4),
                                     1000, rep('', 4),
                                     500,rep('', 4),
                                     0, ''),
                          sec.axis=sec_axis(~., 
                                            breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                            labels=c(2000, rep('', 4),
                                                     1500,rep('', 4),
                                                     1000, rep('', 4),
                                                     500,rep('', 4),
                                                     0, rep('', 4),
                                                     500,rep('', 4),
                                                     1000, rep('', 4),
                                                     1500,rep('', 4),
                                                     2000)
                          )) +
          guides(y = "axis_truncated") +
          labs(y = 'PC2', x = 'Cal years BP', title = '') + 
          theme_classic() + 
          theme(plot.title = element_blank(), 
                axis.text.y = element_text(size = 10, family = 'arial'),
                axis.title.y = element_text(size = 12, family = 'arial'),
                axis.text.x = element_text(size = 10, family = 'arial'),
                axis.title.x = element_text(size = 12, family = 'arial')
          )
        
        #PC3
        PCAScoresPC3 <- ggplot() +
          geom_point(
            aes(x = pcaMatAllTaxa[, 1],
                y = pcaMatAllTaxa[, 10]), 
            size = 2, shape = 21, fill = 'red', alpha = 0.4) + 
          geom_line(aes(y = PCAScoresPC3Gau,
                        x = pcaMatAllTaxa[, 1]), 
                    color = 'black',
                    alpha = 1,
                    linewidth = 1.5) +
          scale_y_continuous(breaks = seq(from = -10, to = 10, by = 0.5),
                             limits = range(pcaMatAllTaxa[, 10]),
                             position = 'left') +
          scale_x_reverse(limits = c(2750, -100), 
                          breaks = seq(from = 4000, to = -100, by = -100),
                          labels = c(4000, rep('', 4),
                                     3500,rep('', 4),
                                     3000, rep('', 4),
                                     2500,rep('', 4),
                                     2000, rep('', 4),
                                     1500,rep('', 4),
                                     1000, rep('', 4),
                                     500,rep('', 4),
                                     0, ''),
                          sec.axis=sec_axis(~., 
                                            breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                            labels=c(2000, rep('', 4),
                                                     1500,rep('', 4),
                                                     1000, rep('', 4),
                                                     500,rep('', 4),
                                                     0, rep('', 4),
                                                     500,rep('', 4),
                                                     1000, rep('', 4),
                                                     1500,rep('', 4),
                                                     2000)
                          )) +
          guides(y = "axis_truncated") +
          labs(y = 'PC3', x = 'Cal years BP', title = '') + 
          theme_classic() + 
          theme(plot.title = element_blank(), 
                axis.text.y = element_text(size = 10, family = 'arial'),
                axis.title.y = element_text(size = 12, family = 'arial'),
                axis.text.x = element_text(size = 10, family = 'arial'),
                axis.title.x = element_text(size = 12, family = 'arial')
          )
        
      }
      
      #Moving correlation
      {
        #Setup environment for parallel processing
        c16 <- makeCluster(16, type = 'SOCK')
        registerDoSNOW(c16)
        
        #PC1
        {
          #Setup moving correlation window size and window moving step
          movCorreWindowSize <- 700
          movCorreWindowStep <- 40
          
          #Extract the matrix of randomed chronologies
          PCScoresPC1 <- PCAScoresPC1Gau
          
          #Calculate the time sequence for moving correlation
          movCorreTimeSeq <- seq(from = min(chronRamMatSph) 
                                 + movCorreWindowSize / 2, to = max(chronRamMatSph) 
                                 - movCorreWindowSize / 2, by = movCorreWindowStep)
          
          movCorreTimeSeq <- movCorreTimeSeq[which(movCorreTimeSeq<= 2400)]
          
          #create a matrix for storing correlation coef. and significance
          movCorreMatPC1 <- matrix(nrow = length(movCorreTimeSeq),
                                      ncol = 4,
                                      dimnames = list(c(),
                                                      c('Time', 'Coef', 'Sig1', 'Sig2')))
          
          movCorreMatPC1[, 1] <- movCorreTimeSeq
          
          chronoCorrMeanMatSca <-
            approx(x = chronologyHM20TransferInte[, 1], 
                   y = chronologyHM20TransferInte[, 2],
                   xout = pcaMatAllTaxa[, 1])[[2]]
          
          movCorreMatPC1[, 2:4] <-
            try(foreach(paraInd = c(1:length(movCorreTimeSeq)), .packages = c("zoo"),
                        .combine = 'rbind') %dopar% 
                  {
                    try(
                      {
                        tempTime1 <- movCorreTimeSeq[paraInd] - movCorreWindowSize / 2
                        tempTime2 <- movCorreTimeSeq[paraInd] + movCorreWindowSize / 2
                        tempIndUAC <- which((pcaMatAllTaxa[, 1]>= tempTime1) & 
                                              (pcaMatAllTaxa[, 1] <= tempTime2))
                        
                        tempIndTSI <- which((TSIRaw[,1] >= tempTime1) & 
                                              (TSIRaw[,1] <= tempTime2))
                        
                        #Correlation between Sphagnum UAC and TSI using linear interpolation and Gaussian filtering
                        time.series1 <- zoo::zoo(PCScoresPC1[tempIndUAC], 
                                                 order.by = pcaMatAllTaxa[, 1][tempIndUAC] - 
                                                   chronoCorrMeanMatSca[tempIndUAC])	#create a zoo-object
                        time.series2 <- zoo::zoo(TSIRaw[tempIndTSI, 2], order.by = TSIRaw[tempIndTSI, 1])
                        Cor <- corit::CorIrregTimser(
                          timser1 = time.series1,
                          timser2 = time.series2,
                          detr = FALSE,	#remove linear trend time series
                          method = "InterpolationMethod",
                          appliedFilter = "gauss",
                          fc = 1/80,	#cut-off frequency
                          dt = 20,	#time step for the interpolation
                          int.method = "linear",	#kind of interpolation
                          filt.output = FALSE)	#return filtered time series 
                        
                        #(2) applying a significance test for the correlation estimate based on the correlation of independent noise
                        slopes <- corit::estimateTimserSlopes(	#estimate spectral slopes of the time series
                          timeseries1 = time.series1,
                          timeseries2 = time.series2,
                          int.step = 1)	#time step of the interpolated time series
                        Quant <- corit::CorQuantilesNullHyp(	#quantiles estimated based on surrogate correlations
                          timser1 = time.series1,
                          timser2 = time.series2,
                          beta.noise1 = slopes$s1,
                          beta.noise2 = slopes$s2,
                          detr = FALSE,
                          rep = 1000,	#repetition during Monte Carlo procedure
                          quant = c(0.05, 0.95),	#quantiles to be estimated
                          method = "InterpolationMethod",
                          appliedFilter = "gauss",
                          fc = 1/80,
                          dt = 20,
                          int.method = "linear")
                        
                        return(c(Cor, Quant[[2]][[1]][1], Quant[[2]][[1]][2]))
                      },
                      silent = TRUE
                    )
                  },
                silent = TRUE
            )
          
          movCorreMatPC1 <- apply(movCorreMatPC1, c(1,2), as.numeric)
        }
        
        #PC2
        {
          #Setup moving correlation window size and window moving step
          movCorreWindowSize <- 700
          movCorreWindowStep <- 40
          
          #Extract the matrix of randomed chronologies
          PCScoresPC2 <- PCAScoresPC2Gau
          
          #Calculate the time sequence for moving correlation
          movCorreTimeSeq <- seq(from = min(chronRamMatSph) 
                                 + movCorreWindowSize / 2, to = max(chronRamMatSph) 
                                 - movCorreWindowSize / 2, by = movCorreWindowStep)
          
          movCorreTimeSeq <- movCorreTimeSeq[which(movCorreTimeSeq<= 2400)]
          
          #create a matrix for storing correlation coef. and significance
          movCorreMatPC2 <- matrix(nrow = length(movCorreTimeSeq),
                                   ncol = 4,
                                   dimnames = list(c(),
                                                   c('Time', 'Coef', 'Sig1', 'Sig2')))
          
          movCorreMatPC2[, 1] <- movCorreTimeSeq
          
          chronoCorrMeanMatSca <-
            approx(x = chronologyHM20TransferInte[, 1], 
                   y = chronologyHM20TransferInte[, 2],
                   xout = pcaMatAllTaxa[, 1])[[2]]
          
          movCorreMatPC2[, 2:4] <-
            try(foreach(paraInd = c(1:length(movCorreTimeSeq)), .packages = c("zoo"),
                        .combine = 'rbind') %dopar% 
                  {
                    try(
                      {
                        tempTime1 <- movCorreTimeSeq[paraInd] - movCorreWindowSize / 2
                        tempTime2 <- movCorreTimeSeq[paraInd] + movCorreWindowSize / 2
                        tempIndUAC <- which((pcaMatAllTaxa[, 1]>= tempTime1) & 
                                              (pcaMatAllTaxa[, 1] <= tempTime2))
                        
                        tempIndTSI <- which((TSIRaw[,1] >= tempTime1) & 
                                              (TSIRaw[,1] <= tempTime2))
                        
                        #Correlation between Sphagnum UAC and TSI using linear interpolation and Gaussian filtering
                        time.series1 <- zoo::zoo(PCScoresPC2[tempIndUAC], 
                                                 order.by = pcaMatAllTaxa[, 1][tempIndUAC] - 
                                                   chronoCorrMeanMatSca[tempIndUAC])	#create a zoo-object
                        time.series2 <- zoo::zoo(TSIRaw[tempIndTSI, 2], order.by = TSIRaw[tempIndTSI, 1])
                        Cor <- corit::CorIrregTimser(
                          timser1 = time.series1,
                          timser2 = time.series2,
                          detr = FALSE,	#remove linear trend time series
                          method = "InterpolationMethod",
                          appliedFilter = "gauss",
                          fc = 1/80,	#cut-off frequency
                          dt = 20,	#time step for the interpolation
                          int.method = "linear",	#kind of interpolation
                          filt.output = FALSE)	#return filtered time series 
                        
                        #(2) applying a significance test for the correlation estimate based on the correlation of independent noise
                        slopes <- corit::estimateTimserSlopes(	#estimate spectral slopes of the time series
                          timeseries1 = time.series1,
                          timeseries2 = time.series2,
                          int.step = 1)	#time step of the interpolated time series
                        Quant <- corit::CorQuantilesNullHyp(	#quantiles estimated based on surrogate correlations
                          timser1 = time.series1,
                          timser2 = time.series2,
                          beta.noise1 = slopes$s1,
                          beta.noise2 = slopes$s2,
                          detr = FALSE,
                          rep = 1000,	#repetition during Monte Carlo procedure
                          quant = c(0.05, 0.95),	#quantiles to be estimated
                          method = "InterpolationMethod",
                          appliedFilter = "gauss",
                          fc = 1/80,
                          dt = 20,
                          int.method = "linear")
                        
                        return(c(Cor, Quant[[2]][[1]][1], Quant[[2]][[1]][2]))
                      },
                      silent = TRUE
                    )
                  },
                silent = TRUE
            )
          
          movCorreMatPC2 <- apply(movCorreMatPC2, c(1,2), as.numeric)
        }
        
        #PC3
        {
          #Setup moving correlation window size and window moving step
          movCorreWindowSize <- 700
          movCorreWindowStep <- 40
          
          #Extract the matrix of randomed chronologies
          PCScoresPC3 <- PCAScoresPC3Gau
          
          #Calculate the time sequence for moving correlation
          movCorreTimeSeq <- seq(from = min(chronRamMatSph) 
                                 + movCorreWindowSize / 2, to = max(chronRamMatSph) 
                                 - movCorreWindowSize / 2, by = movCorreWindowStep)
          
          movCorreTimeSeq <- movCorreTimeSeq[which(movCorreTimeSeq<= 2400)]
          
          #create a matrix for storing correlation coef. and significance
          movCorreMatPC3 <- matrix(nrow = length(movCorreTimeSeq),
                                   ncol = 4,
                                   dimnames = list(c(),
                                                   c('Time', 'Coef', 'Sig1', 'Sig2')))
          
          movCorreMatPC3[, 1] <- movCorreTimeSeq
          
          chronoCorrMeanMatSca <-
            approx(x = chronologyHM20TransferInte[, 1], 
                   y = chronologyHM20TransferInte[, 2],
                   xout = pcaMatAllTaxa[, 1])[[2]]
          
          movCorreMatPC3[, 2:4] <-
            try(foreach(paraInd = c(1:length(movCorreTimeSeq)), .packages = c("zoo"),
                        .combine = 'rbind') %dopar% 
                  {
                    try(
                      {
                        tempTime1 <- movCorreTimeSeq[paraInd] - movCorreWindowSize / 2
                        tempTime2 <- movCorreTimeSeq[paraInd] + movCorreWindowSize / 2
                        tempIndUAC <- which((pcaMatAllTaxa[, 1]>= tempTime1) & 
                                              (pcaMatAllTaxa[, 1] <= tempTime2))
                        
                        tempIndTSI <- which((TSIRaw[,1] >= tempTime1) & 
                                              (TSIRaw[,1] <= tempTime2))
                        
                        #Correlation between Sphagnum UAC and TSI using linear interpolation and Gaussian filtering
                        time.series1 <- zoo::zoo(PCScoresPC3[tempIndUAC], 
                                                 order.by = pcaMatAllTaxa[, 1][tempIndUAC] - 
                                                   chronoCorrMeanMatSca[tempIndUAC])	#create a zoo-object
                        time.series2 <- zoo::zoo(TSIRaw[tempIndTSI, 2], order.by = TSIRaw[tempIndTSI, 1])
                        Cor <- corit::CorIrregTimser(
                          timser1 = time.series1,
                          timser2 = time.series2,
                          detr = FALSE,	#remove linear trend time series
                          method = "InterpolationMethod",
                          appliedFilter = "gauss",
                          fc = 1/80,	#cut-off frequency
                          dt = 20,	#time step for the interpolation
                          int.method = "linear",	#kind of interpolation
                          filt.output = FALSE)	#return filtered time series 
                        
                        #(2) applying a significance test for the correlation estimate based on the correlation of independent noise
                        slopes <- corit::estimateTimserSlopes(	#estimate spectral slopes of the time series
                          timeseries1 = time.series1,
                          timeseries2 = time.series2,
                          int.step = 1)	#time step of the interpolated time series
                        Quant <- corit::CorQuantilesNullHyp(	#quantiles estimated based on surrogate correlations
                          timser1 = time.series1,
                          timser2 = time.series2,
                          beta.noise1 = slopes$s1,
                          beta.noise2 = slopes$s2,
                          detr = FALSE,
                          rep = 1000,	#repetition during Monte Carlo procedure
                          quant = c(0.05, 0.95),	#quantiles to be estimated
                          method = "InterpolationMethod",
                          appliedFilter = "gauss",
                          fc = 1/80,
                          dt = 20,
                          int.method = "linear")
                        
                        return(c(Cor, Quant[[2]][[1]][1], Quant[[2]][[1]][2]))
                      },
                      silent = TRUE
                    )
                  },
                silent = TRUE
            )
          
          movCorreMatPC3 <- apply(movCorreMatPC3, c(1,2), as.numeric)
        }
        
        stopCluster(c16)
        
        #Plot moving correlation
        {
          #PC1
          PCAPC1Plot <- ggplot() +
            geom_line(aes(y = movCorreMatPC1[, 2],
                          x = movCorreMatPC1[, 1]), color = 'black') +
            geom_line(aes(y = movCorreMatPC1[, 3],
                          x = movCorreMatPC1[, 1]), color = 'red', linetype = 2, alpha = 0.6) +
            geom_line(aes(y = movCorreMatPC1[, 4],
                          x = movCorreMatPC1[, 1]), color = 'red', linetype = 2, alpha = 0.6) +
            geom_hline(yintercept = 0, color = 'grey', alpha = 0.6) +
            scale_y_continuous(breaks = seq(from = -10, to = 10, by = 0.5),
                               limits = c(-1,1),
                               position = 'left') +
            scale_x_reverse(limits = c(2750, -100), 
                            breaks = seq(from = 4000, to = -100, by = -100),
                            labels = c(4000, rep('', 4),
                                       3500,rep('', 4),
                                       3000, rep('', 4),
                                       2500,rep('', 4),
                                       2000, rep('', 4),
                                       1500,rep('', 4),
                                       1000, rep('', 4),
                                       500,rep('', 4),
                                       0, ''),
                            sec.axis=sec_axis(~., 
                                              breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                              labels=c(2000, rep('', 4),
                                                       1500,rep('', 4),
                                                       1000, rep('', 4),
                                                       500,rep('', 4),
                                                       0, rep('', 4),
                                                       500,rep('', 4),
                                                       1000, rep('', 4),
                                                       1500,rep('', 4),
                                                       2000)
                            )
            )  +
            guides(y = "axis_truncated", x = "axis_truncated") +
            labs(y = 'PC1&TSI', x = 'Cal years BP', title = '') + 
            theme_classic() + 
            theme(plot.title = element_blank(), 
                  axis.text.x = element_text(size = 10, family = 'arial'),
                  axis.title.x = element_text(size = 12, family = 'arial'),
                  axis.text.y = element_text(size = 10, family = 'arial'),
                  axis.title.y = element_text(size = 12, family = 'arial')
            )
          
          #PC2
          PCAPC2Plot <- ggplot() +
            geom_line(aes(y = movCorreMatPC2[, 2],
                          x = movCorreMatPC2[, 1]), color = 'black') +
            geom_line(aes(y = movCorreMatPC2[, 3],
                          x = movCorreMatPC2[, 1]), color = 'red', linetype = 2, alpha = 0.6) +
            geom_line(aes(y = movCorreMatPC2[, 4],
                          x = movCorreMatPC2[, 1]), color = 'red', linetype = 2, alpha = 0.6) +
            geom_hline(yintercept = 0, color = 'grey', alpha = 0.6) +
            scale_y_continuous(breaks = seq(from = -10, to = 10, by = 0.5),
                               limits = c(-1,1),
                               position = 'left') +
            scale_x_reverse(limits = c(2750, -100), 
                            breaks = seq(from = 4000, to = -100, by = -100),
                            labels = c(4000, rep('', 4),
                                       3500,rep('', 4),
                                       3000, rep('', 4),
                                       2500,rep('', 4),
                                       2000, rep('', 4),
                                       1500,rep('', 4),
                                       1000, rep('', 4),
                                       500,rep('', 4),
                                       0, ''),
                            sec.axis=sec_axis(~., 
                                              breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                              labels=c(2000, rep('', 4),
                                                       1500,rep('', 4),
                                                       1000, rep('', 4),
                                                       500,rep('', 4),
                                                       0, rep('', 4),
                                                       500,rep('', 4),
                                                       1000, rep('', 4),
                                                       1500,rep('', 4),
                                                       2000)
                            )
            )  +
            guides(y = "axis_truncated", x = "axis_truncated") +
            labs(y = 'PC2&TSI', x = 'Cal years BP', title = '') + 
            theme_classic() + 
            theme(plot.title = element_blank(), 
                  axis.text.x = element_text(size = 10, family = 'arial'),
                  axis.title.x = element_text(size = 12, family = 'arial'),
                  axis.text.y = element_text(size = 10, family = 'arial'),
                  axis.title.y = element_text(size = 12, family = 'arial')
            )
          
          #PC3
          PCAPC3Plot <- ggplot() +
            geom_line(aes(y = movCorreMatPC3[, 2],
                          x = movCorreMatPC3[, 1]), color = 'black') +
            geom_line(aes(y = movCorreMatPC3[, 3],
                          x = movCorreMatPC3[, 1]), color = 'red', linetype = 2, alpha = 0.6) +
            geom_line(aes(y = movCorreMatPC3[, 4],
                          x = movCorreMatPC3[, 1]), color = 'red', linetype = 2, alpha = 0.6) +
            geom_hline(yintercept = 0, color = 'grey', alpha = 0.6) +
            scale_y_continuous(breaks = seq(from = -10, to = 10, by = 0.5),
                               limits = c(-1,1),
                               position = 'left') +
            scale_x_reverse(limits = c(2750, -100), 
                            breaks = seq(from = 4000, to = -100, by = -100),
                            labels = c(4000, rep('', 4),
                                       3500,rep('', 4),
                                       3000, rep('', 4),
                                       2500,rep('', 4),
                                       2000, rep('', 4),
                                       1500,rep('', 4),
                                       1000, rep('', 4),
                                       500,rep('', 4),
                                       0, ''),
                            sec.axis=sec_axis(~., 
                                              breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                              labels=c(2000, rep('', 4),
                                                       1500,rep('', 4),
                                                       1000, rep('', 4),
                                                       500,rep('', 4),
                                                       0, rep('', 4),
                                                       500,rep('', 4),
                                                       1000, rep('', 4),
                                                       1500,rep('', 4),
                                                       2000)
                            )
            )  +
            guides(y = "axis_truncated", x = "axis_truncated") +
            labs(y = 'PC3&TSI', x = 'Cal years BP', title = '') + 
            theme_classic() + 
            theme(plot.title = element_blank(), 
                  axis.text.x = element_text(size = 10, family = 'arial'),
                  axis.title.x = element_text(size = 12, family = 'arial'),
                  axis.text.y = element_text(size = 10, family = 'arial'),
                  axis.title.y = element_text(size = 12, family = 'arial')
            )
        }
      }
      
      #Combine plot
      {
        ggarrange(radioDatesPlot, 
                  PCAScoresPC1 +
                    geom_hline(yintercept = PCAScoresPC1GauMean, 
                               linetype = 1, 
                               color = 'blue', alpha = 0.5) +
                    geom_hline(yintercept = PCAScoresPC1GauMean - PCAScoresPC1GauSD, 
                               linetype = 2, 
                               color = 'blue', alpha = 0.5) +
                    geom_hline(yintercept = PCAScoresPC1GauMean - 2 * PCAScoresPC1GauSD, 
                               linetype = 2, 
                               color = 'blue', alpha = 0.5) +
                    geom_hline(yintercept = PCAScoresPC1GauMean + PCAScoresPC1GauSD, 
                               linetype = 2, 
                               color = 'blue', alpha = 0.5) +
                    geom_hline(yintercept = PCAScoresPC1GauMean + 2* PCAScoresPC1GauSD, 
                               linetype = 2, 
                               color = 'blue', alpha = 0.5) +
                    theme(
                      axis.line.x = element_blank(),
                      axis.ticks.x = element_blank(),
                      axis.title.x = element_blank(),
                      axis.text.x = element_blank(),
                    ),
                  WholeTSI +
                    theme(
                      axis.line.x = element_blank(),
                      axis.ticks.x = element_blank(),
                      axis.title.x = element_blank(),
                      axis.text.x = element_blank(),
                    ),
                  PCAPC1Plot,
                  heights = c(0.3,1,1,1),
                  widths = c(1,1,1,1),
                  ncol = 1, nrow = 4,
                  align = "v")
        
        ggarrange(radioDatesPlot, 
                  PCAScoresPC2 +
                    geom_hline(yintercept = PCAScoresPC2GauMean, 
                               linetype = 1, 
                               color = 'blue', alpha = 0.5) +
                    geom_hline(yintercept = PCAScoresPC2GauMean - PCAScoresPC2GauSD, 
                               linetype = 2, 
                               color = 'blue', alpha = 0.5) +
                    geom_hline(yintercept = PCAScoresPC2GauMean - 2 * PCAScoresPC2GauSD, 
                               linetype = 2, 
                               color = 'blue', alpha = 0.5) +
                    geom_hline(yintercept = PCAScoresPC2GauMean + PCAScoresPC2GauSD, 
                               linetype = 2, 
                               color = 'blue', alpha = 0.5) +
                    geom_hline(yintercept = PCAScoresPC2GauMean + 2* PCAScoresPC2GauSD, 
                               linetype = 2, 
                               color = 'blue', alpha = 0.5) +
                    theme(
                      axis.line.x = element_blank(),
                      axis.ticks.x = element_blank(),
                      axis.title.x = element_blank(),
                      axis.text.x = element_blank(),
                    ),
                  WholeTSI +
                    theme(
                      axis.line.x = element_blank(),
                      axis.ticks.x = element_blank(),
                      axis.title.x = element_blank(),
                      axis.text.x = element_blank(),
                    ),
                  PCAPC2Plot,
                  heights = c(0.3,1,1,1),
                  widths = c(1,1,1,1),
                  ncol = 1, nrow = 4,
                  align = "v")
        
        ggarrange(radioDatesPlot, 
                  PCAScoresPC3 +
                    geom_hline(yintercept = PCAScoresPC3GauMean, 
                               linetype = 1, 
                               color = 'blue', alpha = 0.5) +
                    geom_hline(yintercept = PCAScoresPC3GauMean - PCAScoresPC3GauSD, 
                               linetype = 2, 
                               color = 'blue', alpha = 0.5) +
                    geom_hline(yintercept = PCAScoresPC3GauMean - 2 * PCAScoresPC3GauSD, 
                               linetype = 2, 
                               color = 'blue', alpha = 0.5) +
                    geom_hline(yintercept = PCAScoresPC3GauMean + PCAScoresPC3GauSD, 
                               linetype = 2, 
                               color = 'blue', alpha = 0.5) +
                    geom_hline(yintercept = PCAScoresPC3GauMean + 2* PCAScoresPC3GauSD, 
                               linetype = 2, 
                               color = 'blue', alpha = 0.5) +
                    theme(
                      axis.line.x = element_blank(),
                      axis.ticks.x = element_blank(),
                      axis.title.x = element_blank(),
                      axis.text.x = element_blank(),
                    ),
                  WholeTSI +
                    theme(
                      axis.line.x = element_blank(),
                      axis.ticks.x = element_blank(),
                      axis.title.x = element_blank(),
                      axis.text.x = element_blank(),
                    ),
                  PCAPC3Plot,
                  heights = c(0.3,1,1,1),
                  widths = c(1,1,1,1),
                  ncol = 1, nrow = 4,
                  align = "v")
        
      }
    }
    
    #Compare PC scores with other proxy records
    {
      ggarrange(radioDatesPlot, 
                histoEvents +
                  theme(
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                  ),
                PCAScoresPC1 +
                  geom_hline(yintercept = 0, 
                             linetype = 1, 
                             color = 'grey', alpha = 0.9) +
                  theme(
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                  ),
                WholeSSI +
                  theme(
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                  ),
                WholeSI +
                  theme(
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                  ),
                WholeIRD +
                  theme(
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                  ),
                WholetempAnoma,
                heights = c(0.5,0.3,1,1,1,1,1.2),
                widths = c(1,1,1,1,1,1,1),
                ncol = 1, nrow = 7,
                align = "v")
      
      ggarrange(radioDatesPlot, 
                histoEvents +
                  theme(
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                  ),
                PCAScoresPC2 +
                  geom_hline(yintercept = 0, 
                             linetype = 1, 
                             color = 'grey', alpha = 0.9) +
                  theme(
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                  ),
                WholeSSI +
                  theme(
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                  ),
                WholeSI +
                  theme(
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                  ),
                WholeIRD +
                  theme(
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                  ),
                WholetempAnoma,
                heights = c(0.5,0.3,1,1,1,1,1.2),
                widths = c(1,1,1,1,1,1,1),
                ncol = 1, nrow = 7,
                align = "v")
      
      ggarrange(radioDatesPlot, 
                histoEvents +
                  theme(
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                  ),
                PCAScoresPC3 +
                  geom_hline(yintercept = 0, 
                             linetype = 1, 
                             color = 'grey', alpha = 0.9) +
                  theme(
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                  ),
                WholeSSI +
                  theme(
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                  ),
                WholeSI +
                  theme(
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                  ),
                WholeIRD +
                  theme(
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                  ),
                WholetempAnoma,
                heights = c(0.5,0.3,1,1,1,1,1.2),
                widths = c(1,1,1,1,1,1,1),
                ncol = 1, nrow = 7,
                align = "v")
      
      ggarrange(radioDatesPlot, 
                PCAScoresPC1 +
                  scale_x_reverse(limits = c(1000, -100), 
                                  breaks = seq(from = 4000, to = -100, by = -100),
                                  labels = c(4000, rep('', 4),
                                             3500,rep('', 4),
                                             3000, rep('', 4),
                                             2500,rep('', 4),
                                             2000, rep('', 4),
                                             1500,rep('', 4),
                                             1000, rep('', 4),
                                             500,rep('', 4),
                                             0, '')) +
                  theme(
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                  ),
                WholeCC +
                  geom_hline(yintercept = 0, 
                             linetype = 1, 
                             color = 'grey', alpha = 0.6),
                heights = c(0.3,1,1.2),
                widths = c(1,1,1),
                ncol = 1, nrow = 3,
                align = "v")
      
      ggarrange(radioDatesPlot, 
                PCAScoresPC2 +
                  scale_x_reverse(limits = c(1000, -100), 
                                  breaks = seq(from = 4000, to = -100, by = -100),
                                  labels = c(4000, rep('', 4),
                                             3500,rep('', 4),
                                             3000, rep('', 4),
                                             2500,rep('', 4),
                                             2000, rep('', 4),
                                             1500,rep('', 4),
                                             1000, rep('', 4),
                                             500,rep('', 4),
                                             0, '')) +
                  theme(
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                  ),
                WholeCC +
                  geom_hline(yintercept = 0, 
                             linetype = 1, 
                             color = 'grey', alpha = 0.6),
                heights = c(0.3,1,1.2),
                widths = c(1,1,1),
                ncol = 1, nrow = 3,
                align = "v")
      
      ggarrange(radioDatesPlot, 
                PCAScoresPC3 +
                  scale_x_reverse(limits = c(1000, -100), 
                                  breaks = seq(from = 4000, to = -100, by = -100),
                                  labels = c(4000, rep('', 4),
                                             3500,rep('', 4),
                                             3000, rep('', 4),
                                             2500,rep('', 4),
                                             2000, rep('', 4),
                                             1500,rep('', 4),
                                             1000, rep('', 4),
                                             500,rep('', 4),
                                             0, '')) +
                  theme(
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                  ),
                WholeCC +
                  geom_hline(yintercept = 0, 
                             linetype = 1, 
                             color = 'grey', alpha = 0.6),
                heights = c(0.3,1,1.2),
                widths = c(1,1,1),
                ncol = 1, nrow = 3,
                align = "v")
    }
  }
  
  #(Bayesian) multivariate regression
  {
    #TSI and UAC
    {
      #Prepare matrix for multivariate regression
      bayRegUACMat <- pcaMatAllTaxa[,c(-3,-5,-7,-8:-10)]
      bayRegUACMat <- cbind(bayRegUACMat, approx(x = TSIRaw[,1],
                                                 y = TSIRaw[,2],
                                                 xout = bayRegUACMat[,1] +
                                                   chronoCorrMeanMatSca)[[2]])
      colnames(bayRegUACMat)[5] <- 'TSI'
      bayRegUACMat <- as.data.frame(bayRegUACMat)
      bayRegUACMat <- bayRegUACMat[,-1]
      
      #Ordinary liner regression model
      lmMulRegUACTSI <-lm(TSI ~., data = bayRegUACMat)
      summary(lmMulRegUACTSI)
      
      #Bayesian model
      bayMulRegUACTSI <- stan_glm(TSI~., data=bayRegUACMat, seed=111)
      print(bayMulRegUACTSI, digits = 3)
      
      #Distribution of coefficients
      mcmc_dens(bayMulRegUACTSI, pars = c("Sphagnum"))+
        vline_at(-0.040, col="red")
      mcmc_dens(bayMulRegUACTSI, pars = c("Alnus"))+
        vline_at( 0.270, col="red")
      mcmc_dens(bayMulRegUACTSI, pars = c("Calluna"))+
        vline_at( 0.099, col="red")
      
      #Describe the posterior results
      describe_posterior(bayMulRegUACTSI)
      
      #Extract the coefficient values from Bayesian results
      postUACTSI <- get_parameters(bayMulRegUACTSI)
      print(purrr::map_dbl(postUACTSI,median),digits = 3)

    }
    
    #TSI and PC
    {
      #Prepare matrix for multivariate regression
      bayRegPCMat <- pcaMatAllTaxa[,c(-2:-7)]
      bayRegPCMat <- cbind(bayRegPCMat, approx(x = TSIRaw[,1],
                                                 y = TSIRaw[,2],
                                                 xout = bayRegPCMat[,1] +
                                                   chronoCorrMeanMatSca)[[2]])
      colnames(bayRegPCMat)[5] <- 'TSI'
      bayRegPCMat <- as.data.frame(bayRegPCMat)
      bayRegPCMat <- bayRegPCMat[,-1]
      
      #Ordinary liner regression model
      lmMulRegPCTSI <-lm(TSI ~., data = bayRegPCMat)
      summary(lmMulRegPCTSI)
      
      #Bayesian model
      bayMulRegPCTSI <- stan_glm(TSI~., data=bayRegPCMat, seed=111)
      print(bayMulRegPCTSI, digits = 3)
      
      #Distribution of coefficients
      mcmc_dens(bayMulRegPCTSI, pars = c("PC1"))+
        vline_at(-0.079, col="red")
      mcmc_dens(bayMulRegPCTSI, pars = c("PC2"))+
        vline_at( -0.174, col="red")
      mcmc_dens(bayMulRegPCTSI, pars = c("PC3"))+
        vline_at( -0.221, col="red")
      
      #Describe the posterior results
      describe_posterior(bayMulRegPCTSI)
      
      #Extract the coefficient values from Bayesian results
      postPCTSI <- get_parameters(bayMulRegPCTSI)
      print(purrr::map_dbl(postPCTSI,median),digits = 3)
    }
    
    #CC and UAC
    {
      #Prepare matrix for multivariate regression
      bayRegUACMat <- pcaMatAllTaxa[,c(-3,-5,-7,-8:-10)]
      bayRegUACMat <- cbind(bayRegUACMat, approx(x = CCRaw[,1],
                                                 y = CCRaw[,3],
                                                 xout = bayRegUACMat[,1])[[2]])
      colnames(bayRegUACMat)[5] <- 'CC'
      # bayRegUACMat[, 5] <- scale(bayRegUACMat[, 5])
      bayRegUACMat <- as.data.frame(bayRegUACMat)
      bayRegUACMat <- bayRegUACMat[,-1]
      
      #Ordinary liner regression model
      lmMulRegUACCC <-lm(CC ~., data = bayRegUACMat)
      summary(lmMulRegUACCC)
      
      #Bayesian model
      bayMulRegUACCC <- stan_glm(CC~., data=bayRegUACMat, seed=111)
      print(bayMulRegUACCC, digits = 3)
      
      #Distribution of coefficients
      mcmc_dens(bayMulRegUACCC, pars = c("Sphagnum"))+
        vline_at(-0.967, col="red")
      mcmc_dens(bayMulRegUACCC, pars = c("Alnus"))+
        vline_at(2.083, col="red")
      mcmc_dens(bayMulRegUACCC, pars = c("Calluna"))+
        vline_at(0.486, col="red")
      
      #Describe the posterior results
      describe_posterior(bayMulRegUACCC)
      
      #Extract the coefficient values from Bayesian results
      postUACCC <- get_parameters(bayMulRegUACCC)
      print(purrr::map_dbl(postUACCC,median),digits = 3)
    }
    
    #CC and PC
    {
      #Prepare matrix for multivariate regression
      bayRegPCMat <- pcaMatAllTaxa[,c(-2:-7)]
      bayRegPCMat <- cbind(bayRegPCMat, approx(x = CCRaw[,1],
                                               y = CCRaw[,2],
                                               xout = bayRegPCMat[,1] +
                                                 chronoCorrMeanMatSca)[[2]])
      colnames(bayRegPCMat)[5] <- 'CC'
      bayRegPCMat <- as.data.frame(bayRegPCMat)
      bayRegPCMat <- bayRegPCMat[,-1]
      
      #Ordinary liner regression model
      lmMulRegPCCC <-lm(CC ~., data = bayRegPCMat)
      summary(lmMulRegPCCC)
      
      #Bayesian model
      bayMulRegPCCC <- stan_glm(CC~., data=bayRegPCMat, seed=111)
      print(bayMulRegPCCC, digits = 3)
      
      #Distribution of coefficients
      mcmc_dens(bayMulRegPCCC, pars = c("PC1"))+
        vline_at(-1.376, col="red")
      mcmc_dens(bayMulRegPCCC, pars = c("PC2"))+
        vline_at(-1.274, col="red")
      mcmc_dens(bayMulRegPCCC, pars = c("PC3"))+
        vline_at(1.774, col="red")
      
      #Describe the posterior results
      describe_posterior(bayMulRegPCCC)
      
      #Extract the coefficient values from Bayesian results
      postPCCC <- get_parameters(bayMulRegPCCC)
      print(purrr::map_dbl(postPCCC,median),digits = 3)
    }
  }
}

#Figure 8 triangle plot discussing seasonal patterns of three taxa in recording
  #ground UV-B
{
  #Data preparation
  {
    #Calculate the ratio of PC loadings of three taxa
    loadingsRatPC1 <- abs(treat_pca_rotation_PC1$Loadings) / sum(abs(treat_pca_rotation_PC1$Loadings))
    loadingsRatPC2 <- abs(treat_pca_rotation_PC2$Loadings) / sum(abs(treat_pca_rotation_PC2$Loadings))
    loadingsRatPC3 <- abs(treat_pca_rotation_PC3$Loadings) / sum(abs(treat_pca_rotation_PC3$Loadings))
    
    #Combine the data into one matrix
    triPlotPCMat <- rbind(c('PC1', 'PC2', 'PC3'),
                          loadingsRatPC1,
                          loadingsRatPC2,
                          loadingsRatPC3)
    
    #Rename colomns
    rownames(triPlotPCMat) <- c('PC', 'Sphagnum', 'Alnus', 'Calluna')
    
    triPlotPCMat <- as.data.frame(t(triPlotPCMat))
    
    triPlotPCMat[,2:4] <- apply(triPlotPCMat[,2:4], c(1,2), as.numeric)
    
    #Read ground UV-B monitoring data around Manchester
    groundUVB <- read.csv(here('Archives_QSR', 'Manchester_UVB.csv'),
                          header = FALSE)
    
    #Group by months
    groundUVBMon <- 
      aggregate(groundUVB[, 2], 
                by = list(Category = groundUVB[, 1]), 
                FUN=sum)
    
    #Calculate the average monthly doses
    monthResi <- groundUVBMon[, 1] %% 12
    groundUVBMonMean <- 
      aggregate(groundUVBMon[, 2], 
                by = list(Category = monthResi), 
                FUN=mean)
    
    groundUVBMonSD <- 
      aggregate(groundUVBMon[, 2], 
                by = list(Category = monthResi), 
                FUN=range)
    groundUVBMonMean[, 1] <- groundUVBMonMean[, 1] + 1
    
    groundUVBMonSD[, 1] <- groundUVBMonSD[, 1] + 1
    
    colnames(groundUVBMonMean) <- c('Month', 'Mean_UVB')
  }
  
  #Stimulate random combination of growing month of three taxa
  {
    #All one-month
    {
      #Calculate the number of rows need
      combineOneMon <- expand.grid(rep(list(1:12), 3))
      
      #Prepare the data matrix
      ratioOneMon <- cbind(combineOneMon, combineOneMon)
      
      colnames(ratioOneMon) <- c('Mon1', 'Mon2', 'Mon3',
                                 'Sphagnum', 'Alnus', 'Calluna')
      
      #Calculate the dose ratio according to the random combination
      #Setup environment for parallel processing
      c16 <- makeCluster(16, type = 'SOCK')
      registerDoSNOW(c16)
      
      ratioOneMon[, 4:6] <-
        foreach(paraInd = c(1: nrow(ratioOneMon)), .combine = 'rbind') %dopar%
        {
          temp1 <- groundUVBMonMean[
            which(groundUVBMonMean[, 1] == ratioOneMon[paraInd, 1]), 2]
          
          temp2 <- groundUVBMonMean[
            which(groundUVBMonMean[, 1] == ratioOneMon[paraInd, 2]), 2]
          
          temp3 <- groundUVBMonMean[
            which(groundUVBMonMean[, 1] == ratioOneMon[paraInd, 3]), 2]
          
          return(c(temp1, temp2, temp3) /
                   sum(temp1, temp2, temp3))
        }
      
      stopCluster(c16)
      
      #Calculate the point distance (maximum ratio difference) 
        #between random ratios and PC loading ratios
      ratioDiffOneMonPC1 <-
        unlist(lapply(apply(ratioOneMon[, 4:6], 1, function (x) (x - triPlotPCMat[1,-1])), max))
      ratioDiffOneMonPC2 <-
        unlist(lapply(apply(ratioOneMon[, 4:6], 1, function (x) (x - triPlotPCMat[2,-1])), max))
      ratioDiffOneMonPC3 <-
        unlist(lapply(apply(ratioOneMon[, 4:6], 1, function (x) (x - triPlotPCMat[3,-1])), max))
      
      #Find the closest random point
      ratioDiffOneMonPC1Min <- which.min(ratioDiffOneMonPC1)
      ratioDiffOneMonPC2Min <- which.min(ratioDiffOneMonPC2)
      ratioDiffOneMonPC3Min <- which.min(ratioDiffOneMonPC3)
      
    }
    
    #List all possible combinations, given the common pollination months for three
      #taxa
    {
      #Calculate the number of rows need
      combineOneMon <- expand.grid(c(list(5:8), 
                                     list(c(8:12, 1:3)),
                                     list(7:9)))
      
      #Prepare the data matrix
      ratioOneMon <- cbind(combineOneMon, combineOneMon)
      
      colnames(ratioOneMon) <- c('Mon1', 'Mon2', 'Mon3',
                                 'Sphagnum', 'Alnus', 'Calluna')
      
      #Calculate the dose ratio according to the random combination
      #Setup environment for parallel processing
      c16 <- makeCluster(8, type = 'SOCK')
      registerDoSNOW(c16)
      
      ratioOneMon[, 4:6] <-
        foreach(paraInd = c(1: nrow(ratioOneMon)), .combine = 'rbind') %dopar%
        {
          temp1 <- groundUVBMonMean[
            which(groundUVBMonMean[, 1] == ratioOneMon[paraInd, 1]), 2]
          
          temp2 <- groundUVBMonMean[
            which(groundUVBMonMean[, 1] == ratioOneMon[paraInd, 2]), 2]
          
          temp3 <- groundUVBMonMean[
            which(groundUVBMonMean[, 1] == ratioOneMon[paraInd, 3]), 2]
          
          return(c(temp1, temp2, temp3) /
                   sum(temp1, temp2, temp3))
        }
      
      stopCluster(c16)
    }
  }

  #Calculate the synthesised UAC signals of three taxa, with known monthly UV-B
    #ratios as weights
  {
    #Prepare matrix for all possible pollination months combinations
    seasonUACSyn <- matrix(nrow = nrow(pcaMatAllTaxa),
                           ncol = nrow(ratioOneMon))
    
    #Loop to fill the matrix
    loopI <- 1
    while (loopI <= nrow(ratioOneMon)) {
      
      seasonUACSyn[, loopI] <- pcaMatAllTaxa[, 2] * ratioOneMon[loopI, 4] +
        pcaMatAllTaxa[, 4] * ratioOneMon[loopI, 5] +
        pcaMatAllTaxa[, 6] * ratioOneMon[loopI, 6]
      
      loopI <- loopI + 1
    }
    
    #Attach depth information to this matrix
    seasonUACSyn <- cbind(commonBatch[c(-1, -190:-193)] * 2,
                          commonTimeSeq[c(-1, -190:-193)],
                          seasonUACSyn)
  }
  
  #Calculate the correlations between UAC and TSI, referring to the synthesised
    #UAC signals above
  {
    #Prepare a matrix for correlation statistics
    seasonCorMat <- matrix(data = NA,
                           nrow = nrow(ratioOneMon),
                           ncol = ncol(ratioOneMon) + 2,
                           dimnames = list(c(),
                                           c(colnames(ratioOneMon), 
                                             'TSI_Cor',
                                             'CC_Cor')))
    
    #Prepare TSI data for the time sequence analysed
    seasonTSIVec <- approx(x = TSIRaw[, 1],
                           y = TSISmooPlot,
                           xout = seasonUACSyn[, 2])[[2]]
    
    seasonCCVec <- approx(x = CCRaw[, 1],
                          y = CCRaw[, 3],
                          xout = seasonUACSyn[which(seasonUACSyn[, 2] < 1000), 2])[[2]]
    
    #Fill the matrix
    seasonCorMat[,1:6] <- as.matrix(ratioOneMon)
    loopI <- 1
    while (loopI <= nrow(seasonCorMat)) {
      
      seasonCorMat[loopI, 7] <- cor(seasonUACSyn[, loopI + 2], seasonTSIVec,
                                    method = 'pearson',
                                    use = 'pairwise.complete.obs')
      
      seasonCorMat[loopI, 8] <- cor(seasonUACSyn[which(seasonUACSyn[, 2] < 1000), loopI + 2], 
                                    seasonCCVec,
                                    method = 'pearson',
                                    use = 'pairwise.complete.obs')
      
      loopI <- loopI + 1
    }
    
    #Output the synthesised UAC signals
    write.csv(seasonUACSyn[,c(1, 2, which((seasonCorMat[, 7] == min(seasonCorMat[, 7])) |
                                    (seasonCorMat[, 7] == max(seasonCorMat[, 7]))) + 2)], file = here('Archives_QSR', 'seasonUACSyn.csv'))
    
    write.csv(seasonCorMat, file = here('Archives_QSR', 'seasonCorMat.csv'))
  }
  
  #Plot trianlge diagram
  {
    #One month
    ggtern() +
      geom_point(data=triPlotPCMat,aes(x=Sphagnum,y=Alnus, z=Calluna),
                 size = 4, color = c("blue","green",'red')) +
      geom_point(data=ratioOneMon,aes(x=Sphagnum,y=Alnus, z=Calluna),
                 size = 1, alpha = 0.1) +
      geom_point(data=ratioOneMon[c(ratioDiffOneMonPC1Min,
                                    ratioDiffOneMonPC2Min,
                                    ratioDiffOneMonPC3Min),],
                 aes(x=Sphagnum,y=Alnus, z=Calluna),
                 size = 3, alpha = 1, color = 'purple') +
      annotate(geom  = 'text',
               x     = ratioOneMon[ratioDiffOneMonPC1Min, 4],
               y     = ratioOneMon[ratioDiffOneMonPC1Min, 5],
               z     = ratioOneMon[ratioDiffOneMonPC1Min, 6],

               vjust = -0.5,
               hjust = 1.5,
               label = ratioOneMon[ratioDiffOneMonPC1Min,1],
               color = c('blue'),
               size = 6) +
      annotate(geom  = 'text',
               x     = ratioOneMon[ratioDiffOneMonPC1Min, 4],
               y     = ratioOneMon[ratioDiffOneMonPC1Min, 5],
               z     = ratioOneMon[ratioDiffOneMonPC1Min, 6],

               vjust = -0.5,
               hjust = 0.25,
               label = ratioOneMon[ratioDiffOneMonPC1Min,2],
               color = c('red'),
               size = 6) +
      annotate(geom  = 'text',
               x     = ratioOneMon[ratioDiffOneMonPC1Min, 4],
               y     = ratioOneMon[ratioDiffOneMonPC1Min, 5],
               z     = ratioOneMon[ratioDiffOneMonPC1Min, 6],

               vjust = -0.5,
               hjust = -0.5,
               label = ratioOneMon[ratioDiffOneMonPC1Min,3],
               color = c('green'),
               size = 6) +
      
      annotate(geom  = 'text',
               x     = ratioOneMon[ratioDiffOneMonPC2Min, 4],
               y     = ratioOneMon[ratioDiffOneMonPC2Min, 5],
               z     = ratioOneMon[ratioDiffOneMonPC2Min, 6],
               
               vjust = -0.5,
               hjust = 1.5,
               label = ratioOneMon[ratioDiffOneMonPC2Min,1],
               color = c('blue'),
               size = 6) +
      annotate(geom  = 'text',
               x     = ratioOneMon[ratioDiffOneMonPC2Min, 4],
               y     = ratioOneMon[ratioDiffOneMonPC2Min, 5],
               z     = ratioOneMon[ratioDiffOneMonPC2Min, 6],
               
               vjust = -0.5,
               hjust = 0.25,
               label = ratioOneMon[ratioDiffOneMonPC2Min,2],
               color = c('red'),
               size = 6) +
      annotate(geom  = 'text',
               x     = ratioOneMon[ratioDiffOneMonPC2Min, 4],
               y     = ratioOneMon[ratioDiffOneMonPC2Min, 5],
               z     = ratioOneMon[ratioDiffOneMonPC2Min, 6],
               
               vjust = -0.5,
               hjust = -0.5,
               label = ratioOneMon[ratioDiffOneMonPC1Min,3],
               color = c('green'),
               size = 6) +
      
      annotate(geom  = 'text',
               x     = ratioOneMon[ratioDiffOneMonPC3Min, 4],
               y     = ratioOneMon[ratioDiffOneMonPC3Min, 5],
               z     = ratioOneMon[ratioDiffOneMonPC3Min, 6],
               
               vjust = 1,
               hjust = 1.5,
               label = ratioOneMon[ratioDiffOneMonPC3Min,1],
               color = c('blue'),
               size = 6) +
      annotate(geom  = 'text',
               x     = ratioOneMon[ratioDiffOneMonPC3Min, 4],
               y     = ratioOneMon[ratioDiffOneMonPC3Min, 5],
               z     = ratioOneMon[ratioDiffOneMonPC3Min, 6],
               
               vjust = 1,
               hjust = 0.25,
               label = ratioOneMon[ratioDiffOneMonPC3Min,2],
               color = c('red'),
               size = 6) +
      annotate(geom  = 'text',
               x     = ratioOneMon[ratioDiffOneMonPC3Min, 4],
               y     = ratioOneMon[ratioDiffOneMonPC3Min, 5],
               z     = ratioOneMon[ratioDiffOneMonPC3Min, 6],
               
               vjust = 1,
               hjust = -0.5,
               label = ratioOneMon[ratioDiffOneMonPC1Min,3],
               color = c('green'),
               size = 6) +
      
      annotate(geom  = 'text',
               x     = triPlotPCMat[,2],
               y     = triPlotPCMat[,3],
               z     = triPlotPCMat[,4],
               angle = c(0,0,0),
               vjust = c(1,1,-0.4),
               label = triPlotPCMat[,1],
               color = c("blue","green",'red'),
               size = 6) +
      labs(title="PC loading ratios and Monthly UV-B ratio") +
      theme_rgbw()
    
    #Season weighted UACs
    ggtern() +
      geom_point(data = as.data.frame(seasonCorMat[which(seasonCorMat[, 7] >= 0.2),]),
                 aes(x = Sphagnum,y = Alnus, z = Calluna,
                     size = TSI_Cor * 20),
                 color = 'red',
                 alpha = 1) +
      geom_point(data = as.data.frame(seasonCorMat[which(seasonCorMat[, 7] <= 0.05),]),
                 aes(x = Sphagnum,y = Alnus, z = Calluna,
                     size = TSI_Cor * 20),
                 color = 'blue',
                 alpha = 1) +
      geom_point(data = as.data.frame(seasonCorMat[which((seasonCorMat[, 7] < 0.2) &
                                                           (seasonCorMat[, 7] > 0.05)),]),
                 aes(x = Sphagnum,y = Alnus, z = Calluna,
                     size = TSI_Cor * 20),
                 color = 'black',
                 alpha = 0.1) +
      scale_size(name = "Pearson correlation",
                 range = c(0,0.3) * 20,
                 breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3) * 20,
                 labels = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3)) +
      labs(title="Season weighting factors") +
      theme_rgbw()
  }
  
  #Calculate the correlation distribution of weighted UACs, for both most
    #strong and weak cases
  {
    #Loading lipd data
    {
      UACLipdSyn <- lipdR::readLipd(here('Archives_QSR', 'UAC_Syn.lpd'))
    }
    
    #Age model creation
    {
      #Syn
      {
        UACLipdSynWithChron <- runBacon(UACLipdSyn, lab.id.var = 'AMS_Number',
                                        bacon.dir = here(),
                                        age.14c.uncertainty.var = 'age14Cuncertainty', 
                                        age.var = 'Calibrated_Ages', 
                                        age.uncertainty.var = 'Calibrated_Ages_Var',
                                        reservoir.age.14c.var = NULL, 
                                        reservoir.age.14c.uncertainty.var = NULL, 
                                        rejected.ages.var = NULL,
                                        bacon.acc.mean = 20,
                                        cc = 1,
                                        bacon.thick = 10,
                                        suggest = FALSE,
                                        max.ens = 1000,
                                        accept.suggestions = TRUE)
      }
    }
    
    #Syn-Max
    {
      UACMatchingSynMax <- mapAgeEnsembleToPaleoData(UACLipdSynWithChron,age.var = "ageEnsemble",
                                                  paleo.depth.var = "Depth",)
      
      UACMatchingSynMax.ae <- selectData(UACMatchingSynMax,var.name = "ageEnsemble")
      
      
      UACMatchingSynMax.PA <- selectData(UACMatchingSynMax,var.name = "Max")
      
      
      UACMatchingSynMax.PA.ts.plot <- plotTimeseriesEnsRibbons(X = UACMatchingSynMax.ae,
                                                            Y = UACMatchingSynMax.PA,
                                                            n.bins = 1000) +
        geom_point(
          aes(x = UACMatchingSynMax[["paleoData"]][[1]][["measurementTable"]][[1]][["ageMedian"]][["values"]],
              y = UACMatchingSynMax[["paleoData"]][[1]][["measurementTable"]][[1]][["Max"]][["values"]]), 
          size = 2, shape = 21, fill = 'red', alpha = 0.4) + 
        scale_y_continuous(position = 'left', limits = range(UACMatchingSynMax.PA$values,
                                                             UACMatchingSynMax[["paleoData"]][[1]][["measurementTable"]][[1]][["Max"]][["values"]] - 
                                                               UACMatchingSynMax[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]],
                                                             UACMatchingSynMax[["paleoData"]][[1]][["measurementTable"]][[1]][["Max"]][["values"]] + 
                                                               UACMatchingSynMax[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]]),
                           breaks = seq(from = -10, to = 10, by = 0.5)) +
        scale_x_reverse(limits = c(2750, -100), 
                        breaks = seq(from = 4000, to = -100, by = -100),
                        labels = c(4000, rep('', 4),
                                   3500,rep('', 4),
                                   3000, rep('', 4),
                                   2500,rep('', 4),
                                   2000, rep('', 4),
                                   1500,rep('', 4),
                                   1000, rep('', 4),
                                   500,rep('', 4),
                                   0, ''),
                        position = 'top',
                        sec.axis=sec_axis(~., 
                                          breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                          labels=c(2000, rep('', 4),
                                                   1500,rep('', 4),
                                                   1000, rep('', 4),
                                                   500,rep('', 4),
                                                   0, rep('', 4),
                                                   500,rep('', 4),
                                                   1000, rep('', 4),
                                                   1500,rep('', 4),
                                                   2000)
                        )
        ) +
        guides(y = "axis_truncated", x = "axis_truncated") +
        labs(y = 'Average UACs', x = 'Cal years BP', title = '') + 
        theme_classic() + 
        theme(plot.title = element_blank(), 
              axis.text.y = element_text(size = 10, family = 'arial'),
              axis.title.y = element_text(size = 12, family = 'arial'),
              axis.text.x = element_text(size = 10, family = 'arial'),
              axis.title.x = element_text(size = 12, family = 'arial')
        )
      
      
      UACSynMaxChronUnMean <- mean(UACMatchingSynMax.PA.ts.plot[["layers"]][[3]][["data"]][["y"]], 
                                na.rm = TRUE)
      
      UACSynMaxChronUnSD <- sd(UACMatchingSynMax.PA.ts.plot[["layers"]][[3]][["data"]][["y"]], 
                            na.rm = TRUE)
    }
    
    #Syn-Min
    {
      UACMatchingSynMin <- mapAgeEnsembleToPaleoData(UACLipdSynWithChron,age.var = "ageEnsemble",
                                                     paleo.depth.var = "Depth",)
      
      UACMatchingSynMin.ae <- selectData(UACMatchingSynMin,var.name = "ageEnsemble")
      
      
      UACMatchingSynMin.PA <- selectData(UACMatchingSynMin,var.name = "Min")
      
      
      UACMatchingSynMin.PA.ts.plot <- plotTimeseriesEnsRibbons(X = UACMatchingSynMin.ae,
                                                               Y = UACMatchingSynMin.PA,
                                                               n.bins = 1000) +
        geom_point(
          aes(x = UACMatchingSynMin[["paleoData"]][[1]][["measurementTable"]][[1]][["ageMedian"]][["values"]],
              y = UACMatchingSynMin[["paleoData"]][[1]][["measurementTable"]][[1]][["Min"]][["values"]]), 
          size = 2, shape = 21, fill = 'red', alpha = 0.4) + 
        scale_y_continuous(position = 'left', limits = range(UACMatchingSynMin.PA$values,
                                                             UACMatchingSynMin[["paleoData"]][[1]][["measurementTable"]][[1]][["Min"]][["values"]] - 
                                                               UACMatchingSynMin[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]],
                                                             UACMatchingSynMin[["paleoData"]][[1]][["measurementTable"]][[1]][["Min"]][["values"]] + 
                                                               UACMatchingSynMin[["paleoData"]][[1]][["measurementTable"]][[1]][["Peak_Area_SD"]][["values"]]),
                           breaks = seq(from = -10, to = 10, by = 0.5)) +
        scale_x_reverse(limits = c(2750, -100), 
                        breaks = seq(from = 4000, to = -100, by = -100),
                        labels = c(4000, rep('', 4),
                                   3500,rep('', 4),
                                   3000, rep('', 4),
                                   2500,rep('', 4),
                                   2000, rep('', 4),
                                   1500,rep('', 4),
                                   1000, rep('', 4),
                                   500,rep('', 4),
                                   0, ''),
                        position = 'top',
                        sec.axis=sec_axis(~., 
                                          breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                                          labels=c(2000, rep('', 4),
                                                   1500,rep('', 4),
                                                   1000, rep('', 4),
                                                   500,rep('', 4),
                                                   0, rep('', 4),
                                                   500,rep('', 4),
                                                   1000, rep('', 4),
                                                   1500,rep('', 4),
                                                   2000)
                        )
        ) +
        guides(y = "axis_truncated", x = "axis_truncated") +
        labs(y = 'Average UACs', x = 'Cal years BP', title = '') + 
        theme_classic() + 
        theme(plot.title = element_blank(), 
              axis.text.y = element_text(size = 10, family = 'arial'),
              axis.title.y = element_text(size = 12, family = 'arial'),
              axis.text.x = element_text(size = 10, family = 'arial'),
              axis.title.x = element_text(size = 12, family = 'arial')
        )
      
      
      UACSynMinChronUnMean <- mean(UACMatchingSynMin.PA.ts.plot[["layers"]][[3]][["data"]][["y"]], 
                                   na.rm = TRUE)
      
      UACSynMinChronUnSD <- sd(UACMatchingSynMin.PA.ts.plot[["layers"]][[3]][["data"]][["y"]], 
                               na.rm = TRUE)
    }
    
    #Compined plot with TSI
    {
      ggarrange(UACMatchingSynMax.PA.ts.plot,
                UACMatchingSynMin.PA.ts.plot,
                WholeTSI ,
                heights = c(1,1,1),
                widths = c(1,1,1),
                ncol = 1, nrow = 3,
                align = "v")
    }
    
    #Correlation distribution
    {
      #Max
      coroutSynMaxTSI <- corEns(time.1 = UACMatchingSynMax.ae,
                             values.1 = UACMatchingSynMax.PA,
                             time.2 = t(matrix(replicate(1,TSILipd.TSI.Time$values),nrow = 1)),
                             values.2 = TSILipd.TSI,
                             bin.step = 40,
                             max.ens = 1000000,
                             min.obs = 5,
                             isopersistent  = TRUE,
                             isospectral = TRUE)
      
      plotCorEnsSynMaxTSI <- 
        plotCorEns(coroutSynMaxTSI,
                   bins = 20,
                   legend.position =c(.85,.8),
                   f.sig.lab.position = c(.85,.6),
                   significance.option = "isospectral",
                   use.fdr = TRUE)+ggtitle("SynMaxagnum - TSI")
      
      #Min
      coroutSynMinTSI <- corEns(time.1 = UACMatchingSynMin.ae,
                                values.1 = UACMatchingSynMin.PA,
                                time.2 = t(matrix(replicate(1,TSILipd.TSI.Time$values),nrow = 1)),
                                values.2 = TSILipd.TSI,
                                bin.step = 40,
                                max.ens = 1000000,
                                min.obs = 5,
                                isopersistent  = TRUE,
                                isospectral = TRUE)
      
      plotCorEnsSynMinTSI <- 
        plotCorEns(coroutSynMinTSI,
                   bins = 20,
                   legend.position =c(.85,.8),
                   f.sig.lab.position = c(.85,.6),
                   significance.option = "isospectral",
                   use.fdr = TRUE)+ggtitle("SynMinagnum - TSI")
      
      #Combine plot
      plotCorEnsSynMaxTSI / plotCorEnsSynMinTSI
    }
    
  }
  
  
  #Plot monthly UV-B variations
  {
    monthlyUVB <- 
      ggplot() +
      geom_col(aes(x = groundUVBMonMean[, 1],
                    y = groundUVBMonMean[, 2] / 1000),
                linewidth = 1,
                color = 'black',
               alpha = 0.4) +
      scale_y_continuous(
                         breaks = seq(from = 0, to = 18, by = 2.5)) +
      scale_x_continuous(limits = c(0,13),
                         breaks = seq(from = 1, to = 12, by = 1),
                         position = 'top') +

      labs(y = 'Monly mean UV-B (x103 Jm-2)', x = 'Month', title = '') + 
      theme_classic() +
      guides(y = "axis_truncated", x = "axis_truncated") +
      theme(plot.title = element_blank(), 
            axis.text.x = element_text(size = 10, family = 'arial'),
            axis.title.x = element_text(size = 12, family = 'arial'),
            axis.text.y = element_text(size = 10, family = 'arial'),
            axis.title.y = element_text(size = 12, family = 'arial'),
            panel.border = element_rect(color = "black", 
                                        fill = NA, 
                                        linewidth = 1)
      )
  }
  
}
