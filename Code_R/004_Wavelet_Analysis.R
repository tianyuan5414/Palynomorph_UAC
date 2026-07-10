#This script aims at applying wavelet analysis on UAC signals output by '002_UACs_Quantification.r'
#Load necessary packages
{
  #load package for Batlow scientific colour rainbow (Crameri, 2020, Nature Communication)
  library(scico)
  
  #load package for Wavelet analysis
  library(WaveletComp)
}

#Read and plot UAC signals
{
  signalMultipleApproMeanMatSca <- read.csv(file = file.path('Output_R', 'signalMultipleApproMatSca.csv'),
                                            header = TRUE, row.names = 1)
  signalMultipleApproMeanMatScaSmoo <- read.csv(file = file.path('Output_R', 'signalMultipleApproMatScaSmoo.csv'),
                                                header = TRUE, row.names = 1)
}

#Wavelet analysis on UAC signals
{
  #Sphagnum
  {
    aveTimeIntervalSph <- mean(signalMultipleApproMeanMatSca[which(signalMultipleApproMeanMatSca[, 1] == 'Sph'),3][-1] -
                                 signalMultipleApproMeanMatSca[which(signalMultipleApproMeanMatSca[, 1] == 'Sph'),3][-nrow(signalMultipleApproMeanMatSca[which(signalMultipleApproMeanMatSca[, 1] == 'Sph'),])])
    
    #Calculate the expected time domain dequence with given pace length
    interTimeSeqSph <- seq(from = min(signalMultipleApproMeanMatSca[which(signalMultipleApproMeanMatSca[, 1] == 'Sph'),][,3]),
                           to = max(signalMultipleApproMeanMatSca[which(signalMultipleApproMeanMatSca[, 1] == 'Sph'),][,3]),
                           by = aveTimeIntervalSph)
    
    #UAC data
    UACInterPolaSph <- approx(x = signalMultipleApproMeanMatSca[which(signalMultipleApproMeanMatSca[, 1] == 'Sph'),][,3], 
                              y = signalMultipleApproMeanMatSca[which(signalMultipleApproMeanMatSca[, 1] == 'Sph'),][,4],
                              xout = interTimeSeqSph)[[2]]
    
    #Create a dataframe for bi wavelet analysis
    datafraBiWaveSph <- data.frame(Age = interTimeSeqSph, 
                                   Data = UACInterPolaSph)
    #run the wavelet calculations
    DA.wavUACSph <- analyze.wavelet(datafraBiWaveSph, #specify the data frame
                                    "Data", #specify the name of the data
                                    loess.span = 0, #performs a detrending by loess, default is span of 0.75., zero = no detrending
                                    dt = aveTimeIntervalSph, #sampling resolution in time domain, i.e. years between samples
                                    dj = 1/20, #sampling resolution in frequency domain
                                    make.pval = T, #compute p-values?
                                    method = "AR", #method for generating null hypothesis, here AR(1) red noise
                                    n.sim = 999,
                                    lowerPeriod = 64,
                                    upperPeriod = 1024) # how many simulations for p-value
  }
  
  #Alnus
  {
    aveTimeIntervalAln <- mean(signalMultipleApproMeanMatSca[which(signalMultipleApproMeanMatSca[, 1] == 'Aln'),3][-1] -
                                 signalMultipleApproMeanMatSca[which(signalMultipleApproMeanMatSca[, 1] == 'Aln'),3][-nrow(signalMultipleApproMeanMatSca[which(signalMultipleApproMeanMatSca[, 1] == 'Aln'),])])
    
    #Calculate the expected time domain dequence with given pace length
    interTimeSeqAln <- seq(from = min(signalMultipleApproMeanMatSca[which(signalMultipleApproMeanMatSca[, 1] == 'Aln'),][,3]),
                           to = max(signalMultipleApproMeanMatSca[which(signalMultipleApproMeanMatSca[, 1] == 'Aln'),][,3]),
                           by = aveTimeIntervalAln)
    
    #UAC data
    UACInterPolaAln <- approx(x = signalMultipleApproMeanMatSca[which(signalMultipleApproMeanMatSca[, 1] == 'Aln'),][,3], 
                              y = signalMultipleApproMeanMatSca[which(signalMultipleApproMeanMatSca[, 1] == 'Aln'),][,4],
                              xout = interTimeSeqAln)[[2]]
    
    #Create a dataframe for bi wavelet analysis
    datafraBiWaveAln <- data.frame(Age = interTimeSeqAln, 
                                   Data = UACInterPolaAln)
    #run the wavelet calculations
    DA.wavUACAln <- analyze.wavelet(datafraBiWaveAln, #specify the data frame
                                    "Data", #specify the name of the data
                                    loess.span = 0, #performs a detrending by loess, default is span of 0.75., zero = no detrending
                                    dt = aveTimeIntervalAln, #sampling resolution in time domain, i.e. years between samples
                                    dj = 1/20, #sampling resolution in frequency domain
                                    make.pval = T, #compute p-values?
                                    method = "AR", #method for generating null hypothesis, here AR(1) red noise
                                    n.sim = 999,
                                    lowerPeriod = 64,
                                    upperPeriod = 1024) # how many simulations for p-value
  }
  
  #Calluna
  {
    aveTimeIntervalCal <- mean(signalMultipleApproMeanMatSca[which(signalMultipleApproMeanMatSca[, 1] == 'Cal'),3][-1] -
                                 signalMultipleApproMeanMatSca[which(signalMultipleApproMeanMatSca[, 1] == 'Cal'),3][-nrow(signalMultipleApproMeanMatSca[which(signalMultipleApproMeanMatSca[, 1] == 'Cal'),])])
    
    #Calculate the expected time domain dequence with given pace length
    interTimeSeqCal <- seq(from = min(signalMultipleApproMeanMatSca[which(signalMultipleApproMeanMatSca[, 1] == 'Cal'),][,3]),
                           to = max(signalMultipleApproMeanMatSca[which(signalMultipleApproMeanMatSca[, 1] == 'Cal'),][,3]),
                           by = aveTimeIntervalCal)
    
    #UAC data
    UACInterPolaCal <- approx(x = signalMultipleApproMeanMatSca[which(signalMultipleApproMeanMatSca[, 1] == 'Cal'),][,3], 
                              y = signalMultipleApproMeanMatSca[which(signalMultipleApproMeanMatSca[, 1] == 'Cal'),][,4],
                              xout = interTimeSeqCal)[[2]]
    
    #Create a dataframe for bi wavelet analysis
    datafraBiWaveCal <- data.frame(Age = interTimeSeqCal, 
                                   Data = UACInterPolaCal)
    #run the wavelet calculations
    DA.wavUACCal <- analyze.wavelet(datafraBiWaveCal, #specify the data frame
                                    "Data", #specify the name of the data
                                    loess.span = 0, #performs a detrending by loess, default is span of 0.75., zero = no detrending
                                    dt = aveTimeIntervalCal, #sampling resolution in time domain, i.e. years between samples
                                    dj = 1/20, #sampling resolution in frequency domain
                                    make.pval = T, #compute p-values?
                                    method = "AR", #method for generating null hypothesis, here AR(1) red noise
                                    n.sim = 999,
                                    lowerPeriod = 64,
                                    upperPeriod = 1024) # how many simulations for p-value
  }
  
  #Output wavelet analysis results
  {
    #Output as PDF
    pdf(file = file.path('Figures_R', 'Wavelet.pdf'),   # The directory you want to save the file in
        width = 14, # The width of the plot in inches
        height = 8)
    
    #Plot
    WAUACSph <-
      wt.image(DA.wavUACSph, 
               main="Sphagnum", 
               color.key = "quantile",
               color.palette = "rev(scico(n.levels, palette = 'batlow'))",
               siglvl = 0.1, #i.e. 90% confidence level
               n.levels = 100, 
               timelab = "Cal yr BP", 
               legend.params = list(lab = "wavelet power levels", mar=4.7, label.digits = 3), 
               plot.ridge = F,
               spec.time.axis = list(at = (seq(from = -100, to = 2700, by = 100) -
                                             min(interTimeSeqSph))/aveTimeIntervalSph + 1,
                                     labels = c('','0', rep(' ', 9),
                                                '1000', rep(' ', 9),
                                                '2000', rep(' ', 7)),
                                     las = 2, hadj = NA, padj = NA),
               spec.period.axis = list(at = c(64, 128, 256, 512, 1024),
                                       las = 1, hadj = NA, padj = NA),
               graphics.reset = FALSE
      ) 
    
    WAUACAln <-
      wt.image(DA.wavUACAln, 
               main="Alnus", 
               color.key = "quantile",
               color.palette = "rev(scico(n.levels, palette = 'batlow'))",
               siglvl = 0.1, #i.e. 90% confidence level
               n.levels = 100, 
               timelab = "Cal yr BP", 
               legend.params = list(lab = "wavelet power levels", mar=4.7, label.digits = 3), 
               plot.ridge = F,
               spec.time.axis = list(at = (seq(from = -100, to = 2700, by = 100) -
                                             min(interTimeSeqAln))/aveTimeIntervalAln + 1,
                                     labels = c('','0', rep(' ', 9),
                                                '1000', rep(' ', 9),
                                                '2000', rep(' ', 7)),
                                     las = 2, hadj = NA, padj = NA),
               spec.period.axis = list(at = c(64, 128, 256, 512, 1024),
                                       las = 1, hadj = NA, padj = NA),
               graphics.reset = FALSE
      ) 
    
    WAUACCal <-
      wt.image(DA.wavUACCal, 
               main="Calluna", 
               color.key = "quantile",
               color.palette = "rev(scico(n.levels, palette = 'batlow'))",
               siglvl = 0.1, #i.e. 90% confidence level
               n.levels = 100, 
               timelab = "Cal yr BP", 
               legend.params = list(lab = "wavelet power levels", mar=4.7, label.digits = 3), 
               plot.ridge = F,
               spec.time.axis = list(at = (seq(from = -100, to = 2700, by = 100) -
                                             min(interTimeSeqCal))/aveTimeIntervalCal + 1,
                                     labels = c('','0', rep(' ', 9),
                                                '1000', rep(' ', 9),
                                                '2000', rep(' ', 7)),
                                     las = 2, hadj = NA, padj = NA),
               spec.period.axis = list(at = c(64, 128, 256, 512, 1024),
                                       las = 1, hadj = NA, padj = NA),
               graphics.reset = FALSE
      )
    
    dev.off()
  }
}