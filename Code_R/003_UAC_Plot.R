#This script aims at plotting UAC signals output by '002_UACs_Quantification.r'
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
}

#Read and plot UAC signals
{
  signalMultipleApproMeanMatSca <- read.csv(file = file.path('Output_R', 'signalMultipleApproMatSca.csv'),
                                            header = TRUE, row.names = 1)
  signalMultipleApproMeanMatScaSmoo <- read.csv(file = file.path('Output_R', 'signalMultipleApproMatScaSmoo.csv'),
                                                header = TRUE, row.names = 1)
  
  #Sphagnum
  {
    oriUACMeanSph <- signalMultipleApproMeanMatSca[which(signalMultipleApproMeanMatSca[, 1] == 'Sph'),4]
    oriUACSDSph <- signalMultipleApproMeanMatSca[which(signalMultipleApproMeanMatSca[, 1] == 'Sph'),5]
    meanGSmooPlotSphScaMeanPASph <- mean(signalMultipleApproMeanMatScaSmoo[which(signalMultipleApproMeanMatScaSmoo[, 1] == 'Sph'),3])
    meanGSmooPlotSphScaSDPASph <- sd(signalMultipleApproMeanMatScaSmoo[which(signalMultipleApproMeanMatScaSmoo[, 1] == 'Sph'),3])
    
    averageUACStandalonePASph <- ggplot(data = signalMultipleApproMeanMatSca[which(signalMultipleApproMeanMatSca[, 1] == 'Sph'),],) +
      geom_point(aes(x = Date,
                     y = Peak_Area), 
                 size = 2, 
                 shape = 21, 
                 fill = 'red', 
                 alpha = 0.4) + 
      geom_errorbar(aes(x = Date,
                        ymin = Peak_Area - Peak_Area_SD,
                        ymax = Peak_Area + Peak_Area_SD),
                    linewidth = 0.5, 
                    alpha = 0.4, 
                    color = 'black', 
                    width = 20) +
      geom_line(data = signalMultipleApproMeanMatScaSmoo[which(signalMultipleApproMeanMatScaSmoo[, 1] == 'Sph'),],
                aes(y = Peak_Area,
                    x = Date), 
                color = 'black',
                alpha = 1,
                linewidth = 1) +
      geom_hline(yintercept = meanGSmooPlotSphScaMeanPASph, 
                 linetype = 1, 
                 color = 'black', alpha = 0.5) +
      geom_hline(yintercept = meanGSmooPlotSphScaMeanPASph - meanGSmooPlotSphScaSDPASph, 
                 linetype = 2, 
                 color = 'blue', alpha = 0.5) +
      geom_hline(yintercept = meanGSmooPlotSphScaMeanPASph - 2 * meanGSmooPlotSphScaSDPASph, 
                 linetype = 2, 
                 color = 'blue', alpha = 1) +
      geom_hline(yintercept = meanGSmooPlotSphScaMeanPASph + meanGSmooPlotSphScaSDPASph, 
                 linetype = 2, 
                 color = 'red', alpha = 0.5) +
      geom_hline(yintercept = meanGSmooPlotSphScaMeanPASph + 2* meanGSmooPlotSphScaSDPASph, 
                 linetype = 2, 
                 color = 'red', alpha = 1) +
      scale_y_continuous(breaks = seq(from = -10, to = 20, by = 2),
                         limits = range(oriUACMeanSph - oriUACSDSph,
                                        oriUACMeanSph + oriUACSDSph),
                         position = 'left',
                         guide = guide_axis(cap = TRUE)) +
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
      labs(y = 'UAC peak area', x = 'BCE/CE', title = 'Sphagnum') +
      theme_classic() + 
      theme(plot.title = element_text(size = 14, family = 'arial'), 
            axis.text.y = element_text(size = 10, family = 'arial'),
            axis.title.y = element_text(size = 12, family = 'arial'),
            axis.text.x = element_text(size = 10, family = 'arial'),
            axis.title.x = element_text(size = 12, family = 'arial')
      )
  }
  
  #Alnus
  {
    oriUACMeanAln <- signalMultipleApproMeanMatSca[which(signalMultipleApproMeanMatSca[, 1] == 'Aln'),4]
    oriUACSDAln <- signalMultipleApproMeanMatSca[which(signalMultipleApproMeanMatSca[, 1] == 'Aln'),5]
    meanGSmooPlotAlnScaMeanPAAln <- mean(signalMultipleApproMeanMatScaSmoo[which(signalMultipleApproMeanMatScaSmoo[, 1] == 'Aln'),3])
    meanGSmooPlotAlnScaSDPAAln <- sd(signalMultipleApproMeanMatScaSmoo[which(signalMultipleApproMeanMatScaSmoo[, 1] == 'Aln'),3])
    
    averageUACStandalonePAAln <- ggplot(data = signalMultipleApproMeanMatSca[which(signalMultipleApproMeanMatSca[, 1] == 'Aln'),],) +
      geom_point(aes(x = Date,
                     y = Peak_Area), 
                 size = 2, 
                 shape = 21, 
                 fill = 'red', 
                 alpha = 0.4) + 
      geom_errorbar(aes(x = Date,
                        ymin = Peak_Area - Peak_Area_SD,
                        ymax = Peak_Area + Peak_Area_SD),
                    linewidth = 0.5, 
                    alpha = 0.4, 
                    color = 'black', 
                    width = 20) +
      geom_line(data = signalMultipleApproMeanMatScaSmoo[which(signalMultipleApproMeanMatScaSmoo[, 1] == 'Aln'),],
                aes(y = Peak_Area,
                    x = Date), 
                color = 'black',
                alpha = 1,
                linewidth = 1) +
      geom_hline(yintercept = meanGSmooPlotAlnScaMeanPAAln, 
                 linetype = 1, 
                 color = 'black', alpha = 0.5) +
      geom_hline(yintercept = meanGSmooPlotAlnScaMeanPAAln - meanGSmooPlotAlnScaSDPAAln, 
                 linetype = 2, 
                 color = 'blue', alpha = 0.5) +
      geom_hline(yintercept = meanGSmooPlotAlnScaMeanPAAln - 2 * meanGSmooPlotAlnScaSDPAAln, 
                 linetype = 2, 
                 color = 'blue', alpha = 1) +
      geom_hline(yintercept = meanGSmooPlotAlnScaMeanPAAln + meanGSmooPlotAlnScaSDPAAln, 
                 linetype = 2, 
                 color = 'red', alpha = 0.5) +
      geom_hline(yintercept = meanGSmooPlotAlnScaMeanPAAln + 2* meanGSmooPlotAlnScaSDPAAln, 
                 linetype = 2, 
                 color = 'red', alpha = 1) +
      scale_y_continuous(breaks = seq(from = -10, to = 20, by = 2),
                         limits = range(oriUACMeanAln - oriUACSDAln,
                                        oriUACMeanAln + oriUACSDAln),
                         position = 'left',
                         guide = guide_axis(cap = TRUE)) +
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
      labs(y = 'UAC peak area', x = 'BCE/CE', title = 'Alnus') +
      theme_classic() + 
      theme(plot.title = element_text(size = 14, family = 'arial'),
            axis.text.y = element_text(size = 10, family = 'arial'),
            axis.title.y = element_text(size = 12, family = 'arial'),
            axis.text.x = element_text(size = 10, family = 'arial'),
            axis.title.x = element_text(size = 12, family = 'arial')
      )
  }
  
  #Calluna
  {
    oriUACMeanCal <- signalMultipleApproMeanMatSca[which(signalMultipleApproMeanMatSca[, 1] == 'Cal'),4]
    oriUACSDCal <- signalMultipleApproMeanMatSca[which(signalMultipleApproMeanMatSca[, 1] == 'Cal'),5]
    meanGSmooPlotCalScaMeanPACal <- mean(signalMultipleApproMeanMatScaSmoo[which(signalMultipleApproMeanMatScaSmoo[, 1] == 'Cal'),3])
    meanGSmooPlotCalScaSDPACal <- sd(signalMultipleApproMeanMatScaSmoo[which(signalMultipleApproMeanMatScaSmoo[, 1] == 'Cal'),3])
    
    averageUACStandalonePACal <- ggplot(data = signalMultipleApproMeanMatSca[which(signalMultipleApproMeanMatSca[, 1] == 'Cal'),],) +
      geom_point(aes(x = Date,
                     y = Peak_Area), 
                 size = 2, 
                 shape = 21, 
                 fill = 'red', 
                 alpha = 0.4) + 
      geom_errorbar(aes(x = Date,
                        ymin = Peak_Area - Peak_Area_SD,
                        ymax = Peak_Area + Peak_Area_SD),
                    linewidth = 0.5, 
                    alpha = 0.4, 
                    color = 'black', 
                    width = 20) +
      geom_line(data = signalMultipleApproMeanMatScaSmoo[which(signalMultipleApproMeanMatScaSmoo[, 1] == 'Cal'),],
                aes(y = Peak_Area,
                    x = Date), 
                color = 'black',
                alpha = 1,
                linewidth = 1) +
      geom_hline(yintercept = meanGSmooPlotCalScaMeanPACal, 
                 linetype = 1, 
                 color = 'black', alpha = 0.5) +
      geom_hline(yintercept = meanGSmooPlotCalScaMeanPACal - meanGSmooPlotCalScaSDPACal, 
                 linetype = 2, 
                 color = 'blue', alpha = 0.5) +
      geom_hline(yintercept = meanGSmooPlotCalScaMeanPACal - 2 * meanGSmooPlotCalScaSDPACal, 
                 linetype = 2, 
                 color = 'blue', alpha = 1) +
      geom_hline(yintercept = meanGSmooPlotCalScaMeanPACal + meanGSmooPlotCalScaSDPACal, 
                 linetype = 2, 
                 color = 'red', alpha = 0.5) +
      geom_hline(yintercept = meanGSmooPlotCalScaMeanPACal + 2* meanGSmooPlotCalScaSDPACal, 
                 linetype = 2, 
                 color = 'red', alpha = 1) +
      scale_y_continuous(breaks = seq(from = -10, to = 20, by = 2),
                         limits = range(oriUACMeanCal - oriUACSDCal,
                                        oriUACMeanCal + oriUACSDCal),
                         position = 'left',
                         guide = guide_axis(cap = TRUE)) +
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
      labs(y = 'UAC peak area', x = 'BCE/CE', title = 'Calluna') +
      theme_classic() + 
      theme(plot.title = element_text(size = 14, family = 'arial'),
            axis.text.y = element_text(size = 10, family = 'arial'),
            axis.title.y = element_text(size = 12, family = 'arial'),
            axis.text.x = element_text(size = 10, family = 'arial'),
            axis.title.x = element_text(size = 12, family = 'arial')
      )
  }
  
  #Plot the sectioned correlation output by '006_UAC_Peak_Comparison.r'
  {
    correUACSection <- read.csv(file = file.path('Output_R', 'correUACSection.csv'),
                                header = TRUE, row.names = 1)
    
    adjMid <- vector(length = nrow(correUACSection) + 3)
    adjMid[1] <- correUACSection[1, 2]
    adjMid[2] <- mean(c(correUACSection[2, 2], correUACSection[1 , 3]))
    
    loopI <- 2
    loopJ <- 2
    loopJInd <- c()
    while (loopJ <= nrow(correUACSection)) {
      if (loopJ < nrow(correUACSection)) {
        if (correUACSection[loopJ, 1] != correUACSection[loopJ + 1, 1]) {
          
          adjMid[loopI] <- mean(c(correUACSection[loopJ, 2], correUACSection[loopJ - 1, 3]))
          adjMid[loopI + 1] <- correUACSection[loopJ, 3]
          
          adjMid[loopI + 2] <- correUACSection[loopJ + 1, 2]
          
          loopI <- loopI + 3
          
          loopJInd[length(loopJInd) + 1] <- loopJ
          
          loopJ <- loopJ + 2
          
        }else {
          
          adjMid[loopI] <- mean(c(correUACSection[loopJ, 2], correUACSection[loopJ - 1, 3]))
          
          loopI <- loopI + 1
          
          loopJ <- loopJ + 1
        }
        
      }else {
        
        adjMid[loopI] <- mean(c(correUACSection[loopJ, 2], correUACSection[loopJ - 1, 3]))
        
        adjMid[loopI + 1] <- correUACSection[loopJ, 3]
        
        
        loopI <- loopI + 2
        
        loopJ <- loopJ + 2
      }
    }

    correUACSection <- rbind(correUACSection[1:loopJInd[1],],
                             correUACSection[loopJInd[1],],
                             correUACSection[(loopJInd[1] + 1):loopJInd[2],],
                             correUACSection[loopJInd[2],],
                             correUACSection[(loopJInd[2] + 1):nrow(correUACSection),],
                             correUACSection[nrow(correUACSection),])
    correUACSection <- cbind(correUACSection, mid = adjMid)

    
    sectionedCorrelationPlot <- ggplot(data = correUACSection) +
      geom_step(aes(x = mid,
                    y = r,
                    color = Taxa),
                linewidth = 1,
                alpha = 0.6,
                direction = 'vh') +
      scale_color_manual(name = 'Taxa',
                         values = c(brewer.pal(3, 'Set2')
                                    )
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
      scale_y_continuous(
                         position = 'right',
                         guide = guide_axis(cap = TRUE)) +
      labs(y = 'Pearson correlation coefficient (r)', x = 'BCE/CE', title = 'Correlation') +
      theme_classic() + 
      theme(plot.title = element_text(size = 14, family = 'arial'), 
            axis.text.y = element_text(size = 10, family = 'arial'),
            axis.title.y = element_text(size = 12, family = 'arial'),
            axis.text.x = element_text(size = 10, family = 'arial'),
            axis.title.x = element_text(size = 12, family = 'arial'),
            legend.position = 'bottom'
      )
  }
  
  #Combine all plots
  {
    #Output as PDF
    pdf(file = file.path('Figures_R', 'UAC.pdf'),   # The directory you want to save the file in
        width = 10, # The width of the plot in inches
        height = 18)
    
    ggarrange(
      averageUACStandalonePASph,
      averageUACStandalonePAAln +
        scale_y_continuous(breaks = seq(from = -10, to = 20, by = 2),
                           limits = range(oriUACMeanAln - oriUACSDAln,
                                          oriUACMeanAln + oriUACSDAln),
                           position = 'right',
                           guide = guide_axis(cap = TRUE)) +
        theme(plot.title = element_text(size = 14, family = 'arial'),
              axis.text.y = element_text(size = 10, family = 'arial'),
              axis.title.y = element_text(size = 12, family = 'arial'),
              axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.line.x = element_blank(),
              axis.ticks.x = element_blank()
        ),
      averageUACStandalonePACal +
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
        theme(plot.title = element_text(size = 14, family = 'arial'),
              axis.text.y = element_text(size = 10, family = 'arial'),
              axis.title.y = element_text(size = 12, family = 'arial'),
              axis.text.x = element_text(size = 10, family = 'arial'),
              axis.title.x = element_text(size = 12, family = 'arial'),
              # axis.line.x = element_blank(),
              # axis.ticks.x = element_blank()
        ),
      heights = c(1.2,2,2),
      widths = c(1,1,1),
      ncol = 1, nrow = 3,
      align = "v"
    )
    
    dev.off()
  }
}

