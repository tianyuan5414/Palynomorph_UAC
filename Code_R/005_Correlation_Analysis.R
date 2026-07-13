#This script aims at applying correlation analysis on UAC signals output by '002_UACs_Quantification.r'
#Load necessary packages
{
  ###get a sense of significance by comparison with random rednoise series
  library(colorednoise)
  
  #Moving average calculation
  library(zoo)
  
  #PLSR regression
  library(pls)
  
  #Spectral treatment
  library(signal)
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
  
  #Function for linear regression between UAC signals and proxy records
  lmAnalysisTreeTaxa <- function(extData, extName, UACDF){
    
    #Interpolate external proxy records into the time sequence of UAC signals
    intoExtProxy <- approx(extData[, 1], extData[, 2],
                           xout = UACDF[, 1])[[2]]
    
    #Combine external proxy records and UAC signals
    lmDF <- cbind(UACDF, extName = intoExtProxy)
    
    #Apply pls model
    lmModel <- lm(data = lmDF,
                  lmDF[, 5] ~ Sph + Aln + Cal)
    
    #Matrix for recording necessary statistics in the linear model
    lmStaMat <- matrix(nrow = 5,
                       ncol = 3,
                       dimnames = list(c('Coefficient',
                                         'Significance',
                                         'R',
                                         'n',
                                         'F'),
                                       c('Sph',
                                         'Aln',
                                         'Cal')
                                       )
                       )
    
    lmStaMat[1:2,] <- apply(t(summary(lmModel)$coefficients)[c(1,4), -1], c(1,2), round, digits = 4)
    lmStaMat[3,1:2] <- c(summary(lmModel)$r.squared, summary(lmModel)$adj.r.squared)
    lmStaMat[4,1] <- summary(lmModel)$df[2]
    lmStaMat[5,1] <- summary(lmModel)$fstatistic[1]
    
    return(lmStaMat)
  }
}

#Read and plot UAC signals
{
  signalMultipleApproMeanMatSca <- read.csv(file = file.path('Output_R', 'signalMultipleApproMatSca.csv'),
                                            header = TRUE, row.names = 1)
  signalMultipleApproMeanMatSca <- signalMultipleApproMeanMatSca[which(signalMultipleApproMeanMatSca[,3] <= 3000),]
  
  signalMultipleApproMeanMatScaSmoo <- read.csv(file = file.path('Output_R', 'signalMultipleApproMatScaSmoo.csv'),
                                                header = TRUE, row.names = 1)
  signalMultipleApproMeanMatScaSmoo <- signalMultipleApproMeanMatScaSmoo[which(signalMultipleApproMeanMatScaSmoo[,2] <= 3000),]
}

#Pearson correlation analysis
{
  #for 3pt mean
  UAC3Sph <- rollmean(signalMultipleApproMeanMatSca[which(signalMultipleApproMeanMatSca[, 1] == 'Sph'),4],3)
  UAC.AGE3Sph <- rollmean(signalMultipleApproMeanMatSca[which(signalMultipleApproMeanMatSca[, 1] == 'Sph'),3],3)
  
  UAC3Aln <- rollmean(signalMultipleApproMeanMatSca[which(signalMultipleApproMeanMatSca[, 1] == 'Aln'),4],3)
  UAC.AGE3Aln <- rollmean(signalMultipleApproMeanMatSca[which(signalMultipleApproMeanMatSca[, 1] == 'Aln'),3],3)
  
  UAC3Cal <- rollmean(signalMultipleApproMeanMatSca[which(signalMultipleApproMeanMatSca[, 1] == 'Cal'),4],3)
  UAC.AGE3Cal <- rollmean(signalMultipleApproMeanMatSca[which(signalMultipleApproMeanMatSca[, 1] == 'Cal'),3],3)
  
  #Pearson correlation between UAC and other proxy records (smoothed)
  correUACProxy <- matrix(nrow = 3, ncol = 3,
                          dimnames = list(c(),
                                          c('Proxy_Records', 'r','p')))
  
  correUACProxy[, 1] <- c('Sphagnum & Alnus', 
                          'Sphagnum & Calluna',
                          'Alnus & Calluna')
  
  correUACProxy[1, 2] <- my.corr.test2(UAC3Sph, UAC.AGE3Sph, UAC.AGE3Aln, UAC3Aln)$estimate
  
  correUACProxy[1, 3] <- round(my.corr.test2(UAC3Sph, UAC.AGE3Sph, UAC.AGE3Aln, UAC3Aln)$p.value, 4)
  
  correUACProxy[2, 2] <- my.corr.test2(UAC3Sph, UAC.AGE3Sph, UAC.AGE3Cal, UAC3Cal)$estimate
  
  correUACProxy[2, 3] <- round(my.corr.test2(UAC3Sph, UAC.AGE3Sph, UAC.AGE3Cal, UAC3Cal)$p.value, 4)
  
  correUACProxy[3, 2] <- my.corr.test2(UAC3Aln, UAC.AGE3Aln, UAC.AGE3Cal, UAC3Cal)$estimate
  
  correUACProxy[3, 3] <- round(my.corr.test2(UAC3Aln, UAC.AGE3Aln, UAC.AGE3Cal, UAC3Cal)$p.value, 4)
}

#Random noises
{
  #Sphagnum
  {
    cor.pearsonSph <- list() 
    
    cor.pearsonSph <- replicate(n=9999, rednoise.test(max(signalMultipleApproMeanMatScaSmoo[which(signalMultipleApproMeanMatScaSmoo[, 1] == 'Sph'), 2]) -
                                                        min(signalMultipleApproMeanMatScaSmoo[which(signalMultipleApproMeanMatScaSmoo[, 1] == 'Sph'), 2]) + 1, 
                                                      signalMultipleApproMeanMatScaSmoo[which(signalMultipleApproMeanMatScaSmoo[, 1] == 'Sph'), 3]))
    cor.pearsonSph <- unlist(cor.pearsonSph, use.names=FALSE)
    
    correUACRandomSph <- matrix(nrow = 4, ncol = 3,
                                dimnames = list(c(),
                                                c('Level','Lower_r',
                                                  'Upper_r')))
    
    correUACRandomSph[1:4, 1] <- c('99%', '95%', '90%', '80%')
    
    #99% confidence lower and upper
    correUACRandomSph[1, 2:3] <- quantile(as.vector(cor.pearsonSph), probs = c(0.005,0.995))
    #95% confidence lower and upper
    correUACRandomSph[2, 2:3] <- quantile(as.vector(cor.pearsonSph), probs = c(0.025,0.975))
    #90% confidence lower and upper
    correUACRandomSph[3, 2:3] <- quantile(as.vector(cor.pearsonSph), probs = c(0.05,0.95))
    #80% confidence lower and upper
    correUACRandomSph[4, 2:3] <- quantile(as.vector(cor.pearsonSph), probs = c(0.1,0.9))
  }
  
  #Alnus
  {
    cor.pearsonAln <- list() 
    
    cor.pearsonAln <- replicate(n=9999, rednoise.test(max(signalMultipleApproMeanMatScaSmoo[which(signalMultipleApproMeanMatScaSmoo[, 1] == 'Aln'), 2]) -
                                                        min(signalMultipleApproMeanMatScaSmoo[which(signalMultipleApproMeanMatScaSmoo[, 1] == 'Aln'), 2]) + 1, 
                                                      signalMultipleApproMeanMatScaSmoo[which(signalMultipleApproMeanMatScaSmoo[, 1] == 'Aln'), 3]))
    cor.pearsonAln <- unlist(cor.pearsonAln, use.names=FALSE)
    
    correUACRandomAln <- matrix(nrow = 4, ncol = 3,
                                dimnames = list(c(),
                                                c('Level','Lower_r',
                                                  'Upper_r')))
    
    correUACRandomAln[1:4, 1] <- c('99%', '95%', '90%', '80%')
    
    #99% confidence lower and upper
    correUACRandomAln[1, 2:3] <- quantile(as.vector(cor.pearsonAln), probs = c(0.005,0.995))
    #95% confidence lower and upper
    correUACRandomAln[2, 2:3] <- quantile(as.vector(cor.pearsonAln), probs = c(0.025,0.975))
    #90% confidence lower and upper
    correUACRandomAln[3, 2:3] <- quantile(as.vector(cor.pearsonAln), probs = c(0.05,0.95))
    #80% confidence lower and upper
    correUACRandomAln[4, 2:3] <- quantile(as.vector(cor.pearsonAln), probs = c(0.1,0.9))
  }
  
  #Calluna
  {
    cor.pearsonCal <- list() 
    
    cor.pearsonCal <- replicate(n=9999, rednoise.test(max(signalMultipleApproMeanMatScaSmoo[which(signalMultipleApproMeanMatScaSmoo[, 1] == 'Cal'), 2]) -
                                                        min(signalMultipleApproMeanMatScaSmoo[which(signalMultipleApproMeanMatScaSmoo[, 1] == 'Cal'), 2]) + 1, 
                                                      signalMultipleApproMeanMatScaSmoo[which(signalMultipleApproMeanMatScaSmoo[, 1] == 'Cal'), 3]))
    cor.pearsonCal <- unlist(cor.pearsonCal, use.names=FALSE)
    
    correUACRandomCal <- matrix(nrow = 4, ncol = 3,
                                dimnames = list(c(),
                                                c('Level','Lower_r',
                                                  'Upper_r')))
    
    correUACRandomCal[1:4, 1] <- c('99%', '95%', '90%', '80%')
    
    #99% confidence lower and upper
    correUACRandomCal[1, 2:3] <- quantile(as.vector(cor.pearsonCal), probs = c(0.005,0.995))
    #95% confidence lower and upper
    correUACRandomCal[2, 2:3] <- quantile(as.vector(cor.pearsonCal), probs = c(0.025,0.975))
    #90% confidence lower and upper
    correUACRandomCal[3, 2:3] <- quantile(as.vector(cor.pearsonCal), probs = c(0.05,0.95))
    #80% confidence lower and upper
    correUACRandomCal[4, 2:3] <- quantile(as.vector(cor.pearsonCal), probs = c(0.1,0.9))
  }
}

#mlm between TSI, WMI, d18O, temperature, VSSI and three taxa's UAC signals
{
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
  
  #lm regression
  {
    #Prepare a dataframe for mlm regression
    UACAreaMatPlotSph <- signalMultipleApproMeanMatScaSmoo[which(signalMultipleApproMeanMatScaSmoo[, 1] == 'Sph'), ]
    UACAreaMatPlotAln <- signalMultipleApproMeanMatScaSmoo[which(signalMultipleApproMeanMatScaSmoo[, 1] == 'Aln'), ]
    UACAreaMatPlotCal <- signalMultipleApproMeanMatScaSmoo[which(signalMultipleApproMeanMatScaSmoo[, 1] == 'Cal'), ]
    
    #Find the common time sequence
    commonTime <- unique(signalMultipleApproMeanMatSca[, 3])
    
    #Create a matrix and fill in it with Interpolated UAC signals according to common time sequence
    lmUACThreeTaxaDF <- matrix(nrow = length(commonTime),
                              ncol = 4,
                              dimnames = list(c(),
                                              c('Time',
                                                'Sph',
                                                'Aln',
                                                'Cal')
                                              )
                              )
    
    lmUACThreeTaxaDF[, 1] <- commonTime
    lmUACThreeTaxaDF[, 2] <- approx(UACAreaMatPlotSph[, 2],
                                   UACAreaMatPlotSph[, 3],
                                   xout = commonTime)[[2]]
    lmUACThreeTaxaDF[, 3] <- approx(UACAreaMatPlotAln[, 2],
                                   UACAreaMatPlotAln[, 3],
                                   xout = commonTime)[[2]]
    lmUACThreeTaxaDF[, 4] <- approx(UACAreaMatPlotCal[, 2],
                                   UACAreaMatPlotCal[, 3],
                                   xout = commonTime)[[2]]
    lmUACThreeTaxaDF <- as.data.frame(lmUACThreeTaxaDF)
    
    #TSI
    TSIlmModel <- lmAnalysisTreeTaxa(TSISmooPlot, 'TSI', lmUACThreeTaxaDF)
    
    #Greenland temperature
    GreenTemplmModel <- lmAnalysisTreeTaxa(GreenTempBW, 'GreenTemp', lmUACThreeTaxaDF)
    
    #WMI
    WMIlmModel <- lmAnalysisTreeTaxa(WMIBW, 'WMI', lmUACThreeTaxaDF)
    
    #Scotland precipitation
    SUlmModel <- lmAnalysisTreeTaxa(SUBW, 'SU', lmUACThreeTaxaDF)
    
    #WM water table
    WMTDlmModel <- lmAnalysisTreeTaxa(WTDBW, 'WTD', lmUACThreeTaxaDF)
    
    #d18O
    d18OlmModel <- lmAnalysisTreeTaxa(d18OBW, 'd18O', lmUACThreeTaxaDF)
    
    #VSSI
    vssilmModel <- lmAnalysisTreeTaxa(vssiBW, 'VSSI', lmUACThreeTaxaDF)
    
    #Combine all statistic results
    lmAllStatistics <- matrix(nrow = 8,
                              ncol = 10,
                              dimnames = list(c('',
                                                'TSI',
                                                'Greenland temperature',
                                                'WMI',
                                                'Scotland precipitation',
                                                'WM water table',
                                                'd18O',
                                                'VSSI'),
                                              c('Sph',
                                                '',
                                                'Aln',
                                                '',
                                                'Cal',
                                                '',
                                                'R',
                                                'Adj R',
                                                'n',
                                                'F'
                                                )
                                              )
                              )
    
    lmAllStatistics[1, 1:6] <- rep(c('Coefficient', 'Significance'), 3)
    lmAllStatistics[2, ] <- c(TSIlmModel[1:2, 1], 
                              TSIlmModel[1:2, 2],
                              TSIlmModel[1:2, 3],
                              TSIlmModel[3, 1:2],
                              TSIlmModel[4:5, 1])
    lmAllStatistics[3, ] <- c(GreenTemplmModel[1:2, 1], 
                              GreenTemplmModel[1:2, 2],
                              GreenTemplmModel[1:2, 3],
                              GreenTemplmModel[3, 1:2],
                              GreenTemplmModel[4:5, 1])
    lmAllStatistics[4, ] <- c(WMIlmModel[1:2, 1], 
                               WMIlmModel[1:2, 2],
                               WMIlmModel[1:2, 3],
                               WMIlmModel[3, 1:2],
                               WMIlmModel[4:5, 1])
    lmAllStatistics[5, ] <- c(SUlmModel[1:2, 1], 
                              SUlmModel[1:2, 2],
                              SUlmModel[1:2, 3],
                              SUlmModel[3, 1:2],
                              SUlmModel[4:5, 1])
    lmAllStatistics[6, ] <- c(WMTDlmModel[1:2, 1], 
                              WMTDlmModel[1:2, 2],
                              WMTDlmModel[1:2, 3],
                              WMTDlmModel[3, 1:2],
                              WMTDlmModel[4:5, 1])
    lmAllStatistics[7, ] <- c(d18OlmModel[1:2, 1], 
                              d18OlmModel[1:2, 2],
                              d18OlmModel[1:2, 3],
                              d18OlmModel[3, 1:2],
                              d18OlmModel[4:5, 1])
    lmAllStatistics[8, ] <- c(vssilmModel[1:2, 1], 
                              vssilmModel[1:2, 2],
                              vssilmModel[1:2, 3],
                              vssilmModel[3, 1:2],
                              vssilmModel[4:5, 1])
    
    #Output lm regression results
    write.csv(lmAllStatistics, file = file.path('Output_R', 'lmAllStatistics.csv'))
  }

}

#Output correlation results
{
  write.csv(correUACProxy, file = file.path('Output_R', 'correUACProxy.csv'))
  
  correUACRandom <- rbind(cbind(Taxa = 'Sphagnum',
                                correUACRandomSph),
                          cbind(Taxa = 'Alnus',
                                correUACRandomAln),
                          cbind(Taxa = 'Calluna',
                                correUACRandomCal))
  
  write.csv(correUACRandom, file = file.path('Output_R', 'correUACRandom.csv'))
}