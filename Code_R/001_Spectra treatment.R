###This script aims at reading raw spectra, applying baseline correction, and
  #normalise the spectra
#Load necessary packages
{
  #Data reading
  library(here)
  
  #Smoothing
  library(prospectr)
  
  #Baseline correction
  library(baseline)
}

#Declare necessary functions
{
  #Function to normalise the spectral absorbance to mean
  zstandMax <- function(x) {
    return((x - mean(x))/ sd(x))
  }
}

#Read raw spectra
{
  #Sphagnum
  {
    #Main sequence spectra read, and smoothing applied
    {
      spectraInteMatSph <- read.csv(here('Data_R', 'spectraInteMatSph.csv'), 
                                    row.names = 1)
      
      colnames(spectraInteMatSph)[4:ncol(spectraInteMatSph)] <- seq(from = 4000,
                                                                    to = 652, 
                                                                    by = -4)
      
      spectraInteMatSph <- spectraInteMatSph[order(as.numeric(spectraInteMatSph[,3])),]
      
      #Apply Savitsky-Golay smoothing on raw spectra at first, with temporary settings
      #of m = 0, p = 2, w = 11
      specSmooSph <- savitzkyGolay(apply(spectraInteMatSph[,-1:-3], c(1,2), as.numeric),
                                   m = 0, p = 2, w = 11)
      
      #List the matrix according to batch order
      specSmooSph <- specSmooSph[order(as.numeric(spectraInteMatSph[,3])),]
      
      #Record all batch information individually
      batchRawSph <- spectraInteMatSph[order(as.numeric(spectraInteMatSph[,3])),1:3]
      
    }
    
    #Top section spectra read, and smoothing applied
    {
      spectraInteMatSphTop <- read.csv(here('Data_R', 'spectraInteMatSphTop.csv'), 
                                       row.names = 1)
      
      colnames(spectraInteMatSphTop)[4:ncol(spectraInteMatSphTop)] <- seq(from = 4000,
                                                                          to = 652, 
                                                                          by = -4)
      
      spectraInteMatSphTop <- spectraInteMatSphTop[order(as.numeric(spectraInteMatSphTop[,3])),]
      
      #Apply Savitsky-Golay smoothing on raw spectra at first, with temporary settings
      #of m = 0, p = 2, w = 11
      specSmooSphTop <- savitzkyGolay(apply(spectraInteMatSphTop[,-1:-3], c(1,2), as.numeric),
                                      m = 0, p = 2, w = 11)
      
      #List the matrix according to batch order
      specSmooSphTop <- specSmooSphTop[order(as.numeric(spectraInteMatSphTop[,3])),]
      
      #Record all batch information individually
      batchRawSphTop <- spectraInteMatSphTop[order(as.numeric(spectraInteMatSphTop[,3])),1:3]
    }
    
    #Baseline correction on finger print area, and then z-score standardise
    {
      #Baseline correction
      {
        #Main sequence
        {
          #Read unqualified spectra index
          unquaIndex <- apply(as.matrix(read.csv(here('Data_R', 'quaIndex.csv'))), 
                              c(1,2),
                              as.numeric)
          
          unquaIndexSph <- unquaIndex[, 3][!is.na(unquaIndex[, 3])]
          
          #Update the batch information
          batchRawQuaSph <- batchRawSph[-unquaIndexSph, ]
          
          #Update the raw spectral library
          spectraInteMatQuaSph <- spectraInteMatSph[-unquaIndexSph,]
          
          #Polynomial
          spectraInteBCPalySphOri <-  baseline.modpolyfit(specSmooSph[, 546:796],
                                                          degree = 3, tol = 0.00002, 
                                                          rep = 100)[[2]][-unquaIndexSph, ]
          
          #Rename the colomn
          colnames(spectraInteBCPalySphOri) <- seq(1800, 800, by = -4)
          
          #Normalise the spectra to mean
          spectraInteBCPalySph <- t(apply(spectraInteBCPalySphOri, 1, zstandMax)) 
        }
        
        #Top section
        {
          #Read unqualified spectra index
          unquaIndexTop <- apply(as.matrix(read.csv(here('Data_R', 'quaIndexTop.csv'))), 
                              c(1,2),
                              as.numeric)
          
          unquaIndexSphTop <- unquaIndex[, 3][!is.na(unquaIndex[, 3])]
          
          #Update the batch information
          batchRawQuaSphTop <- batchRawSphTop[-unquaIndexSphTop, ]
          
          #Update the raw spectral library
          spectraInteMatQuaSphTop <- spectraInteMatSphTop[-unquaIndexSphTop,]
          
          #Polynomial
          spectraInteBCPalySphTopOri <-  baseline.modpolyfit(specSmooSphTop[, 546:796],
                                                             degree = 3, tol = 0.00002, 
                                                             rep = 100)[[2]][-unquaIndexSphTop, ]
          
          #Rename the colomn
          colnames(spectraInteBCPalySphTopOri) <- seq(1800, 800, by = -4)
          
          
          #z-score standardisation
          spectraInteBCPalySphTop <- t(apply(spectraInteBCPalySphTopOri, 1, zstandMax))
        }
      }
      
      #Combine main sequence and top section
      {
        rawInteSpecQuaSph <- rbind(spectraInteMatQuaSphTop, spectraInteMatQuaSph)
        batchRawQuaSph <- rbind(batchRawQuaSphTop, batchRawQuaSph)
        rawSmooBCSpecMatSph <- rbind(spectraInteBCPalySphTopOri, spectraInteBCPalySphOri)
        rawSmooBCSpecMatSph <- cbind(as.numeric(batchRawQuaSph[, 3]), rawSmooBCSpecMatSph)
        rawSmooBCSpecMatSphNor <- rbind(spectraInteBCPalySphTop, spectraInteBCPalySph)
        rawSmooBCSpecMatSphNor <- cbind(as.numeric(batchRawQuaSph[, 3]), rawSmooBCSpecMatSphNor)
      }
      
      #Output the spectra
      {
        write.csv(rawSmooBCSpecMatSphNor, file = here('Output_R', 
                                                      
                                                      'fingerPrintBCNorSph.csv'))
        
        write.csv(batchRawQuaSph, file = here('Output_R', 
                                                      
                                                      'batchRawQuaSph.csv'))
      }
    }
  }
  
  #Alnus
  {
    #Main sequence spectra read, and smoothing applied
    {
      spectraInteMatAln <- read.csv(here('Data_R', 'spectraInteMatAln.csv'), 
                                    row.names = 1)
      
      colnames(spectraInteMatAln)[4:ncol(spectraInteMatAln)] <- seq(from = 4000,
                                                                    to = 652, 
                                                                    by = -4)
      
      spectraInteMatAln <- spectraInteMatAln[order(as.numeric(spectraInteMatAln[,3])),]
      
      #Apply Savitsky-Golay smoothing on raw spectra at first, with temporary settings
      #of m = 0, p = 2, w = 11
      specSmooAln <- savitzkyGolay(apply(spectraInteMatAln[,-1:-3], c(1,2), as.numeric),
                                   m = 0, p = 2, w = 11)
      
      #List the matrix according to batch order
      specSmooAln <- specSmooAln[order(as.numeric(spectraInteMatAln[,3])),]
      
      #Record all batch information individually
      batchRawAln <- spectraInteMatAln[order(as.numeric(spectraInteMatAln[,3])),1:3]
      
    }
    
    #Top section spectra read, and smoothing applied
    {
      spectraInteMatAlnTop <- read.csv(here('Data_R', 'spectraInteMatAlnTop.csv'), 
                                       row.names = 1)
      
      colnames(spectraInteMatAlnTop)[4:ncol(spectraInteMatAlnTop)] <- seq(from = 4000,
                                                                          to = 652, 
                                                                          by = -4)
      
      spectraInteMatAlnTop <- spectraInteMatAlnTop[order(as.numeric(spectraInteMatAlnTop[,3])),]
      
      #Apply Savitsky-Golay smoothing on raw spectra at first, with temporary settings
      #of m = 0, p = 2, w = 11
      specSmooAlnTop <- savitzkyGolay(apply(spectraInteMatAlnTop[,-1:-3], c(1,2), as.numeric),
                                      m = 0, p = 2, w = 11)
      
      #List the matrix according to batch order
      specSmooAlnTop <- specSmooAlnTop[order(as.numeric(spectraInteMatAlnTop[,3])),]
      
      #Record all batch information individually
      batchRawAlnTop <- spectraInteMatAlnTop[order(as.numeric(spectraInteMatAlnTop[,3])),1:3]
    }
    
    #Baseline correction on finger print area, and then z-score standardise
    {
      #Baseline correction
      {
        #Main sequence
        {
          #Read unqualified spectra index
          unquaIndex <- apply(as.matrix(read.csv(here('Data_R', 'quaIndex.csv'))), 
                              c(1,2),
                              as.numeric)
          
          unquaIndexAln <- unquaIndex[, 1][!is.na(unquaIndex[, 1])]
          
          #Update the batch information
          batchRawQuaAln <- batchRawAln[-unquaIndexAln, ]
          
          #Update the raw spectral library
          spectraInteMatQuaAln <- spectraInteMatAln[-unquaIndexAln,]
          
          #Polynomial
          spectraInteBCPalyAlnOri <-  baseline.modpolyfit(specSmooAln[, 546:796],
                                                          degree = 3, tol = 0.00002, 
                                                          rep = 100)[[2]][-unquaIndexAln, ]
          
          #Rename the colomn
          colnames(spectraInteBCPalyAlnOri) <- seq(1800, 800, by = -4)
          
          #Normalise the spectra to mean
          spectraInteBCPalyAln <- t(apply(spectraInteBCPalyAlnOri, 1, zstandMax)) 
        }
        
        #Top section
        {
          #Read unqualified spectra index
          unquaIndexTop <- apply(as.matrix(read.csv(here('Data_R', 'quaIndexTop.csv'))), 
                                 c(1,2),
                                 as.numeric)
          
          unquaIndexAlnTop <- unquaIndex[, 1][!is.na(unquaIndex[, 1])]
          
          #Update the batch information
          batchRawQuaAlnTop <- batchRawAlnTop[-unquaIndexAlnTop, ]
          
          #Update the raw spectral library
          spectraInteMatQuaAlnTop <- spectraInteMatAlnTop[-unquaIndexAlnTop,]
          
          #Polynomial
          spectraInteBCPalyAlnTopOri <-  baseline.modpolyfit(specSmooAlnTop[, 546:796],
                                                             degree = 3, tol = 0.00002, 
                                                             rep = 100)[[2]][-unquaIndexAlnTop, ]
          
          #Rename the colomn
          colnames(spectraInteBCPalyAlnTopOri) <- seq(1800, 800, by = -4)
          
          
          #z-score standardisation
          spectraInteBCPalyAlnTop <- t(apply(spectraInteBCPalyAlnTopOri, 1, zstandMax))
        }
      }
      
      #Combine main sequence and top section
      {
        rawInteSpecQuaAln <- rbind(spectraInteMatQuaAlnTop, spectraInteMatQuaAln)
        batchRawQuaAln <- rbind(batchRawQuaAlnTop, batchRawQuaAln)
        rawSmooBCSpecMatAln <- rbind(spectraInteBCPalyAlnTopOri, spectraInteBCPalyAlnOri)
        rawSmooBCSpecMatAln <- cbind(as.numeric(batchRawQuaAln[, 3]), rawSmooBCSpecMatAln)
        rawSmooBCSpecMatAlnNor <- rbind(spectraInteBCPalyAlnTop, spectraInteBCPalyAln)
        rawSmooBCSpecMatAlnNor <- cbind(as.numeric(batchRawQuaAln[, 3]), rawSmooBCSpecMatAlnNor)
      }
      
      #Output the spectra
      {
        write.csv(rawSmooBCSpecMatAlnNor, file = here('Output_R', 
                                                      
                                                      'fingerPrintBCNorAln.csv'))
        
        write.csv(batchRawQuaAln, file = here('Output_R', 
                                              
                                              'batchRawQuaAln.csv'))
      }
    }
  }
  
  #Calluna
  {
    #Main sequence spectra read, and smoothing applied
    {
      spectraInteMatCal <- read.csv(here('Data_R', 'spectraInteMatCal.csv'), 
                                    row.names = 1)
      
      colnames(spectraInteMatCal)[4:ncol(spectraInteMatCal)] <- seq(from = 4000,
                                                                    to = 652, 
                                                                    by = -4)
      
      spectraInteMatCal <- spectraInteMatCal[order(as.numeric(spectraInteMatCal[,3])),]
      
      #Apply Savitsky-Golay smoothing on raw spectra at first, with temporary settings
      #of m = 0, p = 2, w = 11
      specSmooCal <- savitzkyGolay(apply(spectraInteMatCal[,-1:-3], c(1,2), as.numeric),
                                   m = 0, p = 2, w = 11)
      
      #List the matrix according to batch order
      specSmooCal <- specSmooCal[order(as.numeric(spectraInteMatCal[,3])),]
      
      #Record all batch information individually
      batchRawCal <- spectraInteMatCal[order(as.numeric(spectraInteMatCal[,3])),1:3]
      
    }
    
    #Top section spectra read, and smoothing applied
    {
      spectraInteMatCalTop <- read.csv(here('Data_R', 'spectraInteMatCalTop.csv'), 
                                       row.names = 1)
      
      colnames(spectraInteMatCalTop)[4:ncol(spectraInteMatCalTop)] <- seq(from = 4000,
                                                                          to = 652, 
                                                                          by = -4)
      
      spectraInteMatCalTop <- spectraInteMatCalTop[order(as.numeric(spectraInteMatCalTop[,3])),]
      
      #Apply Savitsky-Golay smoothing on raw spectra at first, with temporary settings
      #of m = 0, p = 2, w = 11
      specSmooCalTop <- savitzkyGolay(apply(spectraInteMatCalTop[,-1:-3], c(1,2), as.numeric),
                                      m = 0, p = 2, w = 11)
      
      #List the matrix according to batch order
      specSmooCalTop <- specSmooCalTop[order(as.numeric(spectraInteMatCalTop[,3])),]
      
      #Record all batch information individually
      batchRawCalTop <- spectraInteMatCalTop[order(as.numeric(spectraInteMatCalTop[,3])),1:3]
    }
    
    #Baseline correction on finger print area, and then z-score standardise
    {
      #Baseline correction
      {
        #Main sequence
        {
          #Read unqualified spectra index
          unquaIndex <- apply(as.matrix(read.csv(here('Data_R', 'quaIndex.csv'))), 
                              c(1,2),
                              as.numeric)
          
          unquaIndexCal <- unquaIndex[, 2][!is.na(unquaIndex[, 2])]
          
          #Update the batch information
          batchRawQuaCal <- batchRawCal[-unquaIndexCal, ]
          
          #Update the raw spectral library
          spectraInteMatQuaCal <- spectraInteMatCal[-unquaIndexCal,]
          
          #Polynomial
          spectraInteBCPalyCalOri <-  baseline.modpolyfit(specSmooCal[, 546:796],
                                                          degree = 3, tol = 0.00002, 
                                                          rep = 100)[[2]][-unquaIndexCal, ]
          
          #Rename the colomn
          colnames(spectraInteBCPalyCalOri) <- seq(1800, 800, by = -4)
          
          #Normalise the spectra to mean
          spectraInteBCPalyCal <- t(apply(spectraInteBCPalyCalOri, 1, zstandMax)) 
        }
        
        #Top section
        {
          #Read unqualified spectra index
          unquaIndexTop <- apply(as.matrix(read.csv(here('Data_R', 'quaIndexTop.csv'))), 
                                 c(1,2),
                                 as.numeric)
          
          unquaIndexCalTop <- unquaIndex[, 2][!is.na(unquaIndex[, 2])]
          
          #Update the batch information
          batchRawQuaCalTop <- batchRawCalTop[-unquaIndexCalTop, ]
          
          #Update the raw spectral library
          spectraInteMatQuaCalTop <- spectraInteMatCalTop[-unquaIndexCalTop,]
          
          #Polynomial
          spectraInteBCPalyCalTopOri <-  baseline.modpolyfit(specSmooCalTop[, 546:796],
                                                             degree = 3, tol = 0.00002, 
                                                             rep = 100)[[2]][-unquaIndexCalTop, ]
          
          #Rename the colomn
          colnames(spectraInteBCPalyCalTopOri) <- seq(1800, 800, by = -4)
          
          
          #z-score standardisation
          spectraInteBCPalyCalTop <- t(apply(spectraInteBCPalyCalTopOri, 1, zstandMax))
        }
      }
      
      #Combine main sequence and top section
      {
        rawInteSpecQuaCal <- rbind(spectraInteMatQuaCalTop, spectraInteMatQuaCal)
        batchRawQuaCal <- rbind(batchRawQuaCalTop, batchRawQuaCal)
        rawSmooBCSpecMatCal <- rbind(spectraInteBCPalyCalTopOri, spectraInteBCPalyCalOri)
        rawSmooBCSpecMatCal <- cbind(as.numeric(batchRawQuaCal[, 3]), rawSmooBCSpecMatCal)
        rawSmooBCSpecMatCalNor <- rbind(spectraInteBCPalyCalTop, spectraInteBCPalyCal)
        rawSmooBCSpecMatCalNor <- cbind(as.numeric(batchRawQuaCal[, 3]), rawSmooBCSpecMatCalNor)
      }
      
      #Output the spectra
      {
        write.csv(rawSmooBCSpecMatCalNor, file = here('Output_R', 
                                                      
                                                      'fingerPrintBCNorCal.csv'))
        
        write.csv(batchRawQuaCal, file = here('Output_R', 
                                              
                                              'batchRawQuaCal.csv'))
      }
    }
  }
}