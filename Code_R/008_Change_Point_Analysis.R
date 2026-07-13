#This script aims at applying change point analysis on the quantified UAC signals, and compare them with the summer/winter
  #isolation change
#Load necessary packages
{
  library(here)
  
  #Prevent axis from intersection
  library(ggh4x)
  library(ggplot2)
  
  library(ggpubr)
  library(reshape2)
  library(showtext)
  font_add(family = "arial", regular = here('Fonts', 'arial.ttf'))
  showtext_auto()
  library(patchwork)
  
  library(prospectr)
  library(seewave)
  
  #Moving average calculation
  library(zoo)
  
  #Change point
  library(changepoint)
  library(bcp)
}

#Read UAC peak area data
{
  #Read UAC peak area series
  UACAreaMatPlot <- read.csv(here('Output_R', 'signalMultipleApproMatSca.csv'),
                                header = TRUE,
                                row.names =  1)
  UACAreaMatPlotSph <- UACAreaMatPlot[which(UACAreaMatPlot[, 1] == 'Sph'), ]
  UACAreaMatPlotAln <- UACAreaMatPlot[which(UACAreaMatPlot[, 1] == 'Aln'), ]
  UACAreaMatPlotCal <- UACAreaMatPlot[which(UACAreaMatPlot[, 1] == 'Cal'), ]
  
  #Setup the interpolation pace length
  interPaceLength <- 20
  
  #Interpolate the UAC signals to the expected time domain sequence
  UACInterPolaSph <- approx(x = UACAreaMatPlotSph[, 3], y = UACAreaMatPlotSph[, 4],
                            xout = seq(min(UACAreaMatPlotSph[, 3]), max(UACAreaMatPlotSph[, 3]),
                                           by = interPaceLength)
                            )
  UACInterPolaAln <- approx(x = UACAreaMatPlotAln[, 3], y = UACAreaMatPlotAln[, 4],
                            xout = seq(min(UACAreaMatPlotAln[, 3]), max(UACAreaMatPlotAln[, 3]),
                                                                              by = interPaceLength)
                            )
  UACInterPolaCal <- approx(x = UACAreaMatPlotCal[, 3], y = UACAreaMatPlotCal[, 4],
                            xout = seq(min(UACAreaMatPlotCal[, 3]), max(UACAreaMatPlotCal[, 3]),
                                                                              by = interPaceLength)
                            )
  
  #Transfer smoothed UAC series into into time series
  UACSmooPlotSphTS <- UACInterPolaSph[[2]]
  UACSmooPlotAlnTS <- UACInterPolaAln[[2]]
  UACSmooPlotCalTS <- UACInterPolaCal[[2]]
}

#Apply change point analysis (Bayesian model with assumption on normal distribution of residues)
{
  #Prescribe the threshold probability of normal distribution
  CPSignma <- 0.32
  
  #Sphagnum
  {
    #Change point analysis
    CPUACSph <- bcp(UACSmooPlotSphTS, 
                    return.mcmc = FALSE)
    
    bcp_sum <- as.data.frame(summary(CPUACSph))
    # Let's filter the data frame and identify the year:
    bcp_sum$id <- 1:length(UACSmooPlotSphTS)
    sel <- bcp_sum[which(CPUACSph$posterior.prob > CPSignma), ]
    # Get the year:
    CPUACSphTime <- (time(UACSmooPlotSphTS)[sel$id] - 1) * 20 + min(UACAreaMatPlotSph[, 3])
    
    plot(CPUACSph)
  }
  
  #Alnus
  {
    #Change point analysis
    CPUACAln <- bcp(UACSmooPlotAlnTS, 
                    return.mcmc = FALSE)
    
    bcp_sum <- as.data.frame(summary(CPUACAln))
    # Let's filter the data frame and identify the year:
    bcp_sum$id <- 1:length(UACSmooPlotAlnTS)
    sel <- bcp_sum[which(CPUACAln$posterior.prob > CPSignma), ]
    # Get the year:
    CPUACAlnTime <- (time(UACSmooPlotAlnTS)[sel$id] - 1) * 20+ min(UACAreaMatPlotAln[, 3])
    
    plot(CPUACAln)
  }
  
  #Calluna
  {
    #Change point analysis
    CPUACCal <- bcp(UACSmooPlotCalTS, 
                    return.mcmc = FALSE)
    
    bcp_sum <- as.data.frame(summary(CPUACCal))
    # Let's filter the data frame and identify the year:
    bcp_sum$id <- 1:length(UACSmooPlotCalTS)
    sel <- bcp_sum[which(CPUACCal$posterior.prob > CPSignma), ]
    # Get the year:
    CPUACCalTime <- (time(UACSmooPlotCalTS)[sel$id] - 1) * 20+ min(UACAreaMatPlotCal[, 3])
    
    plot(CPUACCal)
  }
}

#Plot results of change point analysis, according to posterior probabilities
{
  #Sphagnum
  {
    #Calculate mean and variance of segmented periods in change point analysis
    #Create a matrix to store the mean and variance for each block
    meanvarCPUACSphMat <- matrix(nrow = length(UACSmooPlotSphTS), 
                                 ncol = 3,
                                 dimnames = list(c(), c('Time', 'Mean', 'Var'))
    )
    
    meanvarCPUACSphMat[, 1] <- UACInterPolaSph[[1]]
    
    #Loop to calculate mean and variance
    loopI <- 1
    tempPos <- 1
    while (loopI <= length(CPUACSphTime)) {
      
      meanvarCPUACSphMat[tempPos:which(meanvarCPUACSphMat[, 1] == CPUACSphTime[loopI]), 2] <-
        mean(
          UACSmooPlotSphTS[tempPos:which(meanvarCPUACSphMat[, 1] == CPUACSphTime[loopI])]
        )
      
      meanvarCPUACSphMat[tempPos:which(meanvarCPUACSphMat[, 1] == CPUACSphTime[loopI]), 3] <-
        sd(
          UACSmooPlotSphTS[tempPos:which(meanvarCPUACSphMat[, 1] == CPUACSphTime[loopI])]
        )
      
      tempPos <- which(meanvarCPUACSphMat[, 1] == CPUACSphTime[loopI]) + 1
      
      loopI <- loopI + 1
      
    }
    
    meanvarCPUACSphMat[tempPos:nrow(meanvarCPUACSphMat), 2] <- mean(
      UACSmooPlotSphTS[tempPos:length(UACSmooPlotSphTS)]
    )
    
    meanvarCPUACSphMat[tempPos:nrow(meanvarCPUACSphMat), 3] <- sd(
      UACSmooPlotSphTS[tempPos:length(UACSmooPlotSphTS)]
    )
    
    #Normal point plot with mean and variance between change points
    CPUACSphPointPlot <- ggplot() + 
      geom_point(aes(x = UACAreaMatPlotSph[, 3],
                     y = UACAreaMatPlotSph[, 4]), size = 2, color = 'black') +
      geom_line(aes(x = meanvarCPUACSphMat[, 1],
                    y = meanvarCPUACSphMat[, 2]),
                color = 'red', size = 4) +
      geom_ribbon(aes(x = meanvarCPUACSphMat[, 1],
                      ymin = meanvarCPUACSphMat[, 2] - meanvarCPUACSphMat[, 3],
                      ymax = meanvarCPUACSphMat[, 2] + meanvarCPUACSphMat[, 3]),
                  fill = 'black',
                  alpha = 0.2
      ) +
      scale_x_reverse(limits = c(2800, -100), 
                      breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                      labels=c(2000, rep('', 4),
                               1500,rep('', 4),
                               1000, rep('', 4),
                               500,rep('', 4),
                               0, rep('', 4),
                               500,rep('', 4),
                               1000, rep('', 4),
                               1500,rep('', 4),
                               2000),
                      position = 'top',
                      guide = guide_axis(cap = TRUE)
                      
      ) +
      scale_y_continuous(breaks = seq(from = -10, to = 10, by = 1),
                         limits = range(meanvarCPUACSphMat[, 2] - meanvarCPUACSphMat[, 3],
                                        meanvarCPUACSphMat[, 2] + meanvarCPUACSphMat[, 3],
                                        na.rm = TRUE),
                         position = 'left',
                         guide = guide_axis(cap = TRUE)) +
      guides(y = "axis_truncated") +
      labs(y = 'Sphagnum', x = 'BCE/CE', title = '') + 
      theme_classic() + 
      theme(plot.title = element_blank(), 
            axis.text.y = element_text(size = 10, family = 'arial'),
            axis.title.y = element_text(size = 12, family = 'arial'),
            axis.text.x = element_text(size = 10, family = 'arial'),
            axis.title.x = element_text(size = 12, family = 'arial')
      )
  }
  
  #Alnus
  {
    #Calculate mean and variance of segmented periods in change point analysis
    #Create a matrix to store the mean and variance for each block
    meanvarCPUACAlnMat <- matrix(nrow = length(UACSmooPlotAlnTS), 
                                 ncol = 3,
                                 dimnames = list(c(), c('Time', 'Mean', 'Var'))
    )
    
    meanvarCPUACAlnMat[, 1] <- UACInterPolaAln[[1]]
    
    #Loop to calculate mean and variance
    loopI <- 1
    tempPos <- 1
    while (loopI <= length(CPUACAlnTime)) {
      
      meanvarCPUACAlnMat[tempPos:which(meanvarCPUACAlnMat[, 1] == CPUACAlnTime[loopI]), 2] <-
        mean(
          UACSmooPlotAlnTS[tempPos:which(meanvarCPUACAlnMat[, 1] == CPUACAlnTime[loopI])]
        )
      
      meanvarCPUACAlnMat[tempPos:which(meanvarCPUACAlnMat[, 1] == CPUACAlnTime[loopI]), 3] <-
        sd(
          UACSmooPlotAlnTS[tempPos:which(meanvarCPUACAlnMat[, 1] == CPUACAlnTime[loopI])]
        )
      
      tempPos <- which(meanvarCPUACAlnMat[, 1] == CPUACAlnTime[loopI]) + 1
      
      loopI <- loopI + 1
      
    }
    
    meanvarCPUACAlnMat[tempPos:nrow(meanvarCPUACAlnMat), 2] <- mean(
      UACSmooPlotAlnTS[tempPos:length(UACSmooPlotAlnTS)]
    )
    
    meanvarCPUACAlnMat[tempPos:nrow(meanvarCPUACAlnMat), 3] <- sd(
      UACSmooPlotAlnTS[tempPos:length(UACSmooPlotAlnTS)]
    )
    
    #Normal point plot with mean and variance between change points
    CPUACAlnPointPlot <- ggplot() + 
      geom_point(aes(x = UACAreaMatPlotAln[, 3],
                     y = UACAreaMatPlotAln[, 4]), size = 2, color = 'black') +
      geom_line(aes(x = meanvarCPUACAlnMat[, 1],
                    y = meanvarCPUACAlnMat[, 2]),
                color = 'red', size = 4) +
      geom_ribbon(aes(x = meanvarCPUACAlnMat[, 1],
                      ymin = meanvarCPUACAlnMat[, 2] - meanvarCPUACAlnMat[, 3],
                      ymax = meanvarCPUACAlnMat[, 2] + meanvarCPUACAlnMat[, 3]),
                  fill = 'black',
                  alpha = 0.2
      ) +
      scale_x_reverse(limits = c(2800, -100), 
                      breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                      labels=c(2000, rep('', 4),
                               1500,rep('', 4),
                               1000, rep('', 4),
                               500,rep('', 4),
                               0, rep('', 4),
                               500,rep('', 4),
                               1000, rep('', 4),
                               1500,rep('', 4),
                               2000),
                      position = 'top',
                      guide = guide_axis(cap = TRUE)
                      
      ) +
      scale_y_continuous(breaks = seq(from = -10, to = 10, by = 0.5),
                         limits = range(meanvarCPUACAlnMat[, 2] - meanvarCPUACAlnMat[, 3],
                                        meanvarCPUACAlnMat[, 2] + meanvarCPUACAlnMat[, 3]),
                         position = 'right',
                         guide = guide_axis(cap = TRUE)) +
      guides(y = "axis_truncated") +
      labs(y = 'Alnus', x = 'BCE/CE', title = '') + 
      theme_classic() + 
      theme(plot.title = element_blank(), 
            axis.text.y = element_text(size = 10, family = 'arial'),
            axis.title.y = element_text(size = 12, family = 'arial'),
            axis.text.x = element_text(size = 10, family = 'arial'),
            axis.title.x = element_text(size = 12, family = 'arial')
      )
  }
  
  #Calluna
  {
    #Calculate mean and variance of segmented periods in change point analysis
    #Create a matrix to store the mean and variance for each block
    meanvarCPUACCalMat <- matrix(nrow = length(UACSmooPlotCalTS), 
                                 ncol = 3,
                                 dimnames = list(c(), c('Time', 'Mean', 'Var'))
    )
    
    meanvarCPUACCalMat[, 1] <- UACInterPolaCal[[1]]
    
    #Loop to calculate mean and variance
    loopI <- 1
    tempPos <- 1
    while (loopI <= length(CPUACCalTime)) {
      
      meanvarCPUACCalMat[tempPos:which(meanvarCPUACCalMat[, 1] == CPUACCalTime[loopI]), 2] <-
        mean(
          UACSmooPlotCalTS[tempPos:which(meanvarCPUACCalMat[, 1] == CPUACCalTime[loopI])]
        )
      
      meanvarCPUACCalMat[tempPos:which(meanvarCPUACCalMat[, 1] == CPUACCalTime[loopI]), 3] <-
        sd(
          UACSmooPlotCalTS[tempPos:which(meanvarCPUACCalMat[, 1] == CPUACCalTime[loopI])]
        )
      
      tempPos <- which(meanvarCPUACCalMat[, 1] == CPUACCalTime[loopI]) + 1
      
      loopI <- loopI + 1
      
    }
    
    meanvarCPUACCalMat[tempPos:nrow(meanvarCPUACCalMat), 2] <- mean(
      UACSmooPlotCalTS[tempPos:length(UACSmooPlotCalTS)]
    )
    
    meanvarCPUACCalMat[tempPos:nrow(meanvarCPUACCalMat), 3] <- sd(
      UACSmooPlotCalTS[tempPos:length(UACSmooPlotCalTS)]
    )
    
    #Normal point plot with mean and variance between change points
    CPUACCalPointPlot <- ggplot() + 
      geom_point(aes(x = UACAreaMatPlotCal[, 3],
                     y = UACAreaMatPlotCal[, 4]), size = 2, color = 'black') +
      geom_line(aes(x = meanvarCPUACCalMat[, 1],
                    y = meanvarCPUACCalMat[, 2]),
                color = 'red', size = 4) +
      geom_ribbon(aes(x = meanvarCPUACCalMat[, 1],
                      ymin = meanvarCPUACCalMat[, 2] - meanvarCPUACCalMat[, 3],
                      ymax = meanvarCPUACCalMat[, 2] + meanvarCPUACCalMat[, 3]),
                  fill = 'black',
                  alpha = 0.2
      ) +
      scale_x_reverse(limits = c(2800, -100), 
                      breaks= 1950 - seq(from = -2000, to = 2000, by = 100), 
                      labels=c(2000, rep('', 4),
                               1500,rep('', 4),
                               1000, rep('', 4),
                               500,rep('', 4),
                               0, rep('', 4),
                               500,rep('', 4),
                               1000, rep('', 4),
                               1500,rep('', 4),
                               2000),
                      position = 'bottom',
                      guide = guide_axis(cap = TRUE)
                      
      ) +
      scale_y_continuous(breaks = seq(from = -10, to = 10, by = 1),
                         limits = range(meanvarCPUACCalMat[, 2] - meanvarCPUACCalMat[, 3],
                                        meanvarCPUACCalMat[, 2] + meanvarCPUACCalMat[, 3]),
                         position = 'left',
                         guide = guide_axis(cap = TRUE)) +
      guides(y = "axis_truncated") +
      labs(y = 'Calluna', x = 'BCE/CE', title = '') + 
      theme_classic() + 
      theme(plot.title = element_blank(), 
            axis.text.y = element_text(size = 10, family = 'arial'),
            axis.title.y = element_text(size = 12, family = 'arial'),
            axis.text.x = element_text(size = 10, family = 'arial'),
            axis.title.x = element_text(size = 12, family = 'arial')
      )
  }
  
  #Combine all plots
  {
    #Output as PDF
    pdf(file = file.path('Figures_R', 'UAC_CP.pdf'),   # The directory you want to save the file in
        width = 10, # The width of the plot in inches
        height = 8)
    
    ggarrange(CPUACSphPointPlot,
              CPUACAlnPointPlot +
                theme(plot.title = element_blank(), 
                      axis.text.y = element_text(size = 10, family = 'arial'),
                      axis.title.y = element_text(size = 12, family = 'arial'),
                      axis.text.x = element_blank(),
                      axis.line.x = element_blank(),
                      axis.ticks.x = element_blank(),
                      axis.title.x = element_blank()
                ),
              CPUACCalPointPlot,
              heights = c(1.2,1,1.4),
              nrow = 3,
              ncol = 1,
              align = 'v'
    )
    
    dev.off()
  }
}