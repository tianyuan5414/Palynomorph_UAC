#This script aims at applying correlation analysis on UAC signals output by '002_UACs_Quantification.r'
#Load necessary packages
{
  ###get a sense of significance by comparison with random rednoise series
  library(colorednoise)
  
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