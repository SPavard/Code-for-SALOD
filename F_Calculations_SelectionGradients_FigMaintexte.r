        ###############################################################
        # CALCULATION OF SELECTION COEFFICIENTS FOR DIFFERENT SCENARIOS
			# This code produces results of figure 4 of the main text if the chosen mortality is that of mean hunter-gatherers
        
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

			####################################################################################################################################################################
			# Scenario 1 and 6 - Fertility of males is that of females (+ no matrimonial structure), no maternal, grandpaternal and paternal investment.            
            
				# Loading the table of sigma(y1,y2,y3) with no investment - ATTENTION! ON my computer: record vecSigma over 5 columns intsead of a vector! Can mess the results!
					vec.Sigma <- read.table("C:/Users/Pavard/Documents/LaNouvelleSouille/P_Selection_Old_Age_NouvelleSouille2/Results_1/i.No.Invest_HG_fertf.csv", sep = ";", header=FALSE)

				# Solving Population Dynamics 
					# ChosenFertMale == ChosenFertFemale
					pop <- Solver(20, Ff, Lf, Lm, Alphaf, Betaf, Alphaf, Betaf, Omega, vecSigma,ChosenFertFemale, ChosenFertFemale)
					lambda <- pop$conv.Lambda[20]
					newParAmales <- solvA (lambda, Lf, Alphaf, Betaf, Alphaf, Betaf,ChosenFertFemale, ChosenFertFemale)
					newParFertMales <- c(newParAmales, ChosenFertFemale[2:6])
					the.Fm <- Fert(Age, newParFertMales)
					matB <- Union.x1.x3 (Age, dage, lambda, Lf, Lm, Ff, the.Fm)
					
					
				# Scenario 1 - No variance in age at onset (FAO=MAO)
					the.MAO <- seq(40,85,by=1) 
					the.FAO <- the.MAO - 0.5
					res <- CoefSel(the.MAO, the.FAO)
					write.table(res, 'clipboard', sep='\t', row.names=FALSE, col.names=FALSE)
					
				# Scenario 6 -Variance in age at onset (FAO=MAO-20)	
					the.MAO <- seq(40,85,by=1) 
					
					# the.MAO <- seq(40,45,by=5)
					
					the.FAO <- the.MAO - 20
					res <- CoefSel(the.MAO, the.FAO)
					write.table(res, 'clipboard', sep='\t', row.names=FALSE, col.names=FALSE)
					
			####################################################################################################################################################################		
			# Scenario 2 and 7 - Fertility of males is that of females (+ no matrimonial structure), + maternal care but no grandpaternal and paternal investment.            
					
				# Loading the table of sigma(y1,y2,y3) with no investment
					vec.Sigma <- read.table("C:/Users/Pavard/Documents/LaNouvelleSouille/P_Selection_Old_Age_NouvelleSouille2/Results_1/ii.m.only_HG_fertf.csv", sep = ";", header=FALSE)
					dim(vec.Sigma)			
					vec.Sigma <- as.matrix(vec.Sigma) ; vec.Sigma <- as.vector(vec.Sigma) ; 
					length(vec.Sigma)			
					
					
				# Solving Population Dynamics 
					# ChosenFertMale == ChosenFertFemale
					pop <- Solver(20, Ff, Lf, Lm, Alphaf, Betaf, Alphaf, Betaf, Omega, vecSigma,ChosenFertFemale, ChosenFertFemale)
					lambda <- pop$conv.Lambda[20]
					newParAmales <- solvA (lambda, Lf, Alphaf, Betaf, Alphaf, Betaf,ChosenFertFemale, ChosenFertFemale)
					newParFertMales <- c(newParAmales, ChosenFertFemale[2:6])
					the.Fm <- Fert(Age, newParFertMales)
					matB <- Union.x1.x3 (Age, dage, lambda, Lf,  Lm, Ff, the.Fm)
					
					
				# Scenario 2 - No variance in age at onset (FAO=MAO)
					the.MAO <- seq(40,85,by=1) 
					the.FAO <- the.MAO - 0.5
					res <- CoefSel(the.MAO, the.FAO)
					write.table(res, 'clipboard', sep='\t', row.names=FALSE, col.names=FALSE)
					
				# Scenario 7 -Variance in age at onset (FAO=MAO-20)	
					the.MAO <- seq(40,85,by=1) 
					the.FAO <- the.MAO - 20
					res <- CoefSel(the.MAO, the.FAO)
					write.table(res, 'clipboard', sep='\t', row.names=FALSE, col.names=FALSE)					
					
			####################################################################################################################################################################		
			# Scenario 3 and 8 - Fertility of males is that of females (+ no matrimonial structure), + maternal care + grandpaternal care but no paternal investment.            
					
				# Loading the table of sigma(y1,y2,y3) with no investment
					vec.Sigma <- read.table("C:/Users/Pavard/Documents/LaNouvelleSouille/P_Selection_Old_Age_NouvelleSouille2/Results_1/v.m.gm_HG_fertf.csv", sep = ";", header=FALSE)
					dim(vec.Sigma)			
					vec.Sigma <- as.matrix(vec.Sigma) ; vec.Sigma <- as.vector(vec.Sigma) ; 
					length(vec.Sigma)	
					
					
				# Solving Population Dynamics 
					# ChosenFertMale == ChosenFertFemale
					pop <- Solver(20, Ff, Lf, Lm, Alphaf, Betaf, Alphaf, Betaf, Omega, vecSigma,ChosenFertFemale, ChosenFertFemale)
					lambda <- pop$conv.Lambda[20]
					newParAmales <- solvA (lambda, Lf, Alphaf, Betaf, Alphaf, Betaf,ChosenFertFemale, ChosenFertFemale)
					newParFertMales <- c(newParAmales, ChosenFertFemale[2:6])
					the.Fm <- Fert(Age, newParFertMales)
					matB <- Union.x1.x3 (Age, dage, lambda, Lf,  Lm, Ff, the.Fm)
					
				# Scenario 3 - No variance in age at onset (FAO=MAO)
					the.MAO <- seq(40,85,by=1) 
					the.FAO <- the.MAO - 0.5
					res <- CoefSel(the.MAO, the.FAO)
					write.table(res, 'clipboard', sep='\t', row.names=FALSE, col.names=FALSE)
					
				# Scenario 8 -Variance in age at onset (FAO=MAO-20)	
					the.MAO <- seq(40,85,by=1) 
					the.FAO <- the.MAO - 20
					res <- CoefSel(the.MAO, the.FAO)
					write.table(res, 'clipboard', sep='\t', row.names=FALSE, col.names=FALSE)					
					
			####################################################################################################################################################################
			# Scenario 4 and 9 - Fertility of males is different from that of females (+ matrimonial structure), + maternal care + grandpaternal care but no paternal investment.            
				# Loading the table of sigma(y1,y2,y3) with no investment	
					vec.Sigma <- read.table("C:/Users/Pavard/Documents/LaNouvelleSouille/P_Selection_Old_Age_NouvelleSouille2/Results_1/v.m.gm_HG_fertf.csv", sep = ";", header=FALSE)
					dim(vec.Sigma)			
					vec.Sigma <- as.matrix(vec.Sigma) ; vec.Sigma <- as.vector(vec.Sigma) ; 
					length(vec.Sigma)	
					
				# Solving Population Dynamics 
					# ChosenFertMale != ChosenFertFemale
					pop <- Solver(20, Ff, Lf, Lm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, ChosenFertMale, ChosenFertFemale)
					lambda <- pop$conv.Lambda[20]
					newParAmales <- solvA (lambda, Lf, Alpham, Betam, Alphaf, Betaf, ChosenFertMale, ChosenFertFemale)
					newParFertMales <- c(newParAmales, ChosenFertMale[2:6])
					the.Fm <- Fert(Age, newParFertMales)
					matB <- Union.x1.x3 (Age, dage, lambda, Lf,  Lm, Ff, the.Fm)	
					
				# Scenario 4 - No variance in age at onset (FAO=MAO)
					the.MAO <- seq(40,85,by=1) 
					the.FAO <- the.MAO - 0.5
					res <- CoefSel(the.MAO, the.FAO)
					write.table(res, 'clipboard', sep='\t', row.names=FALSE, col.names=FALSE)
					
				# Scenario 9 -Variance in age at onset (FAO=MAO-20)	
					the.MAO <- seq(40,85,by=1) 
					the.FAO <- the.MAO - 20
					res <- CoefSel(the.MAO, the.FAO)
					write.table(res, 'clipboard', sep='\t', row.names=FALSE, col.names=FALSE)					
			
			####################################################################################################################################################################			
			# Scenario 5 and 10 - Fertility of males is different from that of females (+ matrimonial structure), + maternal care + grandpaternal care + paternal care.            
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
					matB <- Union.x1.x3 (Age, dage, lambda, Lf,  Lm, Ff, the.Fm)		
					
				# Scenario 5 - No variance in age at onset (FAO=MAO)
					the.MAO <- seq(40,85,by=1) 
					the.FAO <- the.MAO - 0.5
					res <- CoefSel(the.MAO, the.FAO)
					write.table(res, 'clipboard', sep='\t', row.names=FALSE, col.names=FALSE)
					
				# Scenario 10 -Variance in age at onset (FAO=MAO-20)	
					the.MAO <- seq(40,85,by=1) 
					the.FAO <- the.MAO - 20
					res <- CoefSel(the.MAO, the.FAO)
					write.table(res, 'clipboard', sep='\t', row.names=FALSE, col.names=FALSE)			
					
					
					
					
					
				