        ###############################################################
        # CALCULATION OF SELECTION COEFFICIENTS FOR A LARGE RANGE OF CUMULATIVE RISK IN THE CASE WHERE ALL SOCIO-CULTURAL FACTORS ARE INCLUDES
			
            # FUNCTIONS FOR DISTRIBUTION OF AGE AT ONSET            
				# Variance of  logistic function as a function of MAO, FAO and threshold (the cumulative percentage of case observed at FAO)
				Fvariance<-function(MAO, FAO, Threshold){
				  y<--(FAO-MAO)/(log(1/Threshold-1))
				  return(y)
				}
				
				CumRisk<-function(Age, MAO, FAO, Threshold){
				  y<-plogis(Age, location = MAO, scale = Fvariance(MAO, FAO, Threshold))
				  y[Age<FAO]<-0
				  return(y)
				}

            # FUNCTIONS FOR CALCULATING SELECTION COEFFICIENT            
            
            CoefSel <- function (vecMAO, vecFAO){
                res     <- rep(NA, length(vecMAO))
                
                the.Ld <- rep(1, length(Age))
                W.NC.Female <- Results.C.female (FALSE, FALSE, FALSE, the.Ld, lambda, Ff, Lf, Lm, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
				W.NC.Male   <- Results.C.male   (FALSE, FALSE, FALSE, the.Ld, lambda, Ff, Lf, Lm, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
				#scaling <- W.NC.Female/W.NC.Male
                W.NC <- 0.5*W.NC.Female + 0.5*W.NC.Male 
                
					for (j in 1:length(vecMAO)){
						print(vecMAO[j])
						the.Ld  <- 1- CumRisk(Age, vecMAO[j], vecFAO[j], 0.01)

						W.males      <- Results.C.male   (FALSE, FALSE, TRUE,  the.Ld, lambda, Ff, Lf, Lm, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
						W.females.1  <- Results.C.female (TRUE, FALSE, FALSE,  the.Ld, lambda, Ff, Lf, Lm, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
						W.females.2  <- Results.C.female (TRUE, TRUE, FALSE,   the.Ld, lambda, Ff, Lf, Lm, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
						
						the.W <- 0.5*W.males + 0.25*W.females.1 + 0.25*W.females.2 
						the.S <- 1- the.W/W.NC
						res[j] <- the.S
					}
                return(res)
            }

			# SCEANRIO Fertility of males is different from that of females (+ matrimonial structure), + maternal care + grandpaternal care + paternal care
				
				# Loading the table of sigma(y1,y2,y3) with no investment	
					vec.Sigma <- read.table("C:/Users/Pavard/Documents/LaNouvelleSouille/P_Selection_Old_Age_NouvelleSouille2/Results_1/vi.m.gm.f_HG.csv", sep = ";", header=FALSE)
					length(vecSigma)
					
				# Solving Population Dynamics 
					# ChosenFertMale != ChosenFertFemale
					pop <- Solver(20, Ff, Lf, Lm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, ChosenFertMale, ChosenFertFemale)
					lambda <- pop$conv.Lambda[20]
					newParAmales <- solvA (lambda, Lf, Alpham, Betam, Alphaf, Betaf, ChosenFertMale, ChosenFertFemale)
					newParFertMales <- c(newParAmales, ChosenFertMale[2:6])
					the.Fm <- Fert(Age, newParFertMales)
					matB <- Union.x1.x3 (Age, dage, lambda, Lf, Ff, the.Fm)		
			
			# Range of MAO and FAO	
					the.MAO <- seq(45,85,by=1)   ; length(the.MAO) # 40:70
					the.FAO <- seq(20,80,by=2.5) ; length(the.FAO) # 20:65
					resmat  <- matrix (NA, nrow=length(the.FAO), ncol=length(the.MAO))

			# Filling the matrix of results		
					for (j in 1:length(the.MAO)){
					  for (i in 1:length(the.FAO)){
						print(j)
						if (the.MAO[j]-the.FAO[i]>5){
							print(i)
							res <- CoefSel(the.MAO[j], the.FAO[i])
							resmat[i,j] <- res}
						else
							resmat[i,j] <- NA
					  }
					} 
      
			# Table of results 		
			write(resmat, file = "C:/Users/Pavard/Documents/LaNouvelleSouille/P_Selection_Old_Age_NouvelleSouille2/Results_2/Table_RES2.csv", sep = ";")	
			
			
			