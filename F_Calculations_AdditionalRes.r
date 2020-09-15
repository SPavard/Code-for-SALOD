        ###############################################################
        # CALCULATION OF SELECTION COEFFICIENTS FOR DIFFERENT SCENARIOS
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

            # FUNCTIONS FOR CALCULATING SELECTION COEFFICIENT WITH INCOMPLETE PENETRANCE           
            CoefSel <- function (vecMAO, vecFAO, Penetrance){
                res     <- rep(NA, length(vecMAO))
                
                the.Ld <- rep(1, length(Age))
                W.NC.Female <- Results.C.female (FALSE, FALSE, FALSE, the.Ld, lambda, Ff, Lf, Lm, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
				W.NC.Male   <- Results.C.male   (FALSE, FALSE, FALSE, the.Ld, lambda, Ff, Lf, Lm, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
				#scaling <- W.NC.Female/W.NC.Male
                W.NC <- 0.5*W.NC.Female + 0.5*W.NC.Male 
                
					for (j in 1:length(vecMAO)){
						print(vecMAO[j])
						the.Ld  <- 1- (CumRisk(Age, vecMAO[j], vecFAO[j], 0.01)*Penetrance)

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
			####################################################################################################################################################################			
			# Scenario 5 and 10 - Fertility of males is different from that of females (+ matrimonial structure), + maternal care + grandpaternal care + paternal care.            
				# Loading the table of sigma(y1,y2,y3) with no investment	
					# vec.Sigma <- read.table("C:/Users/Pavard/Documents/LaNouvelleSouille/P_Selection_Old_Age_NouvelleSouille2/Results_1/vi.m.gm.f_HG.csv", sep = ";", header=FALSE)
					length(vecSigma)			
					#vec.Sigma <- as.matrix(vec.Sigma) ; vec.Sigma <- as.vector(vec.Sigma) ; 
					#length(vec.Sigma)		
					
				# Solving Population Dynamics 
					# ChosenFertMale != ChosenFertFemale
					pop <- Solver(20, Ff, Lf, Lm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, ChosenFertMale, ChosenFertFemale)
					lambda <- pop$conv.Lambda[20]
					newParAmales <- solvA (lambda, Lf, Alpham, Betam, Alphaf, Betaf, ChosenFertMale, ChosenFertFemale)
					newParFertMales <- c(newParAmales, ChosenFertMale[2:6])
					the.Fm <- Fert(Age, newParFertMales)
					matB <- Union.x1.x3 (Age, dage, lambda, Lf, Lm, Ff, the.Fm)		
					lambda

			####################################################################################################################################################################
			####################################################################################################################################################################			
			# CUMULATIVE PENETRANCE

				# Scenario 10 -Variance in age at onset (FAO=MAO-20)	
					the.MAO <- seq(40,85,by=2.5) 
					the.FAO <- the.MAO - 20
					
					res100 <- CoefSel(the.MAO, the.FAO, 1)
					res50  <- CoefSel(the.MAO, the.FAO, 0.5)
					res10  <- CoefSel(the.MAO, the.FAO, 0.1)
					res1   <- CoefSel(the.MAO, the.FAO, 0.01)	

				#plot

				  age <- seq(40,85,by=2.5) 
				  the.lwd <- 1
				  
				  plot (age, res100, col="white", bty="l", type="b",pch=1, cex.pch=2, log="y", xlab="Mean age at onset (MAO)", ylab="Selection Coefficient", ylim=c(2.5e-07, 1))
				  
				  vage <- c(38,87)
				  # 4Nes<100 => 2.5e-02
				  polygon (c(vage, rev(vage)), c(2.5e-03,2.5e-03,1,1), col="gray60", border=F)
				  # 4Nes<100 => 2.5e-03
				  polygon (c(vage, rev(vage)), c(2.5e-03,2.5e-03,1,1), col="gray60", border=F)
				  # 4Nes=1000 => 2.5e-04
				  polygon (c(vage, rev(vage)), c(2.5e-03,2.5e-03,2.5e-04,2.5e-04), col="gray70", border=F)
				  # 4Nes=10000 => 2.5e-05
				  polygon (c(vage, rev(vage)), c(2.5e-04,2.5e-04,2.5e-05,2.5e-05), col="gray80", border=F)
				  # 4Nes>10000 => 2.5e-05
				  polygon (c(vage, rev(vage)), c(2.5e-05,2.5e-05,2.5e-07,2.5e-07), col="gray90", border=F)

				abline(h=2.5e-03, col="grey", lty=3)
				abline(h=2.5e-04, col="grey", lty=3)
				abline(h=2.5e-05, col="grey", lty=3)
				  
				lines(age, res100, type="b",pch=1, lwd=2)
				lines(age, res50,  type="b",pch=6, lwd=2)
				lines(age, res10,  type="b",pch=15, lwd=2)
				lines(age, res1,   type="b",pch=8, lwd=2)
				

			####################################################################################################################################################################
			####################################################################################################################################################################			
			# RECESIVE, Mt, Y, DISEASE IN ONE SEX	
			
			# A recessive autosomal allele (the probability for the mother of a female homozygous carrier of being herself homozygous is considered as negligible and the female's mother is considered as non-carrier in equation	
			vecMAO <- the.MAO
			vecFAO <- the.FAO
			res     <- rep(NA, length(vecMAO))
            
				the.Ld <- rep(1, length(Age))
				W.NC.Female <- Results.C.female (FALSE, FALSE, FALSE, the.Ld, lambda, Ff, Lf, Lm, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
				W.NC.Male   <- Results.C.male   (FALSE, FALSE, FALSE, the.Ld, lambda, Ff, Lf, Lm, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
				W.NC <- 0.5*W.NC.Female + 0.5*W.NC.Male 
                
					for (j in 1:length(vecMAO)){
						print(vecMAO[j])
						the.Ld  <- 1- (CumRisk(Age, vecMAO[j], vecFAO[j], 0.01))

						W.males      <- Results.C.male   (FALSE, FALSE, TRUE,  the.Ld, lambda, Ff, Lf, Lm, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
						W.females.1  <- Results.C.female (TRUE, FALSE, FALSE,  the.Ld, lambda, Ff, Lf, Lm, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
						# W.females.2  <- Results.C.female (TRUE, TRUE, FALSE,   the.Ld, lambda, Ff, Lf, Lm, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
						
						the.W <- 0.5*W.males + 0.5*W.females.1
						the.S <- 1- the.W/W.NC
						res[j] <- the.S
					}
			resRec <- res
			
			# Autosomal SALOD with disease in female only
				res     <- rep(NA, length(vecMAO))
 				the.Ld <- rep(1, length(Age))
				W.NC.Female <- Results.C.female (FALSE, FALSE, FALSE, the.Ld, lambda, Ff, Lf, Lm, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
				W.NC.Male   <- Results.C.male   (FALSE, FALSE, FALSE, the.Ld, lambda, Ff, Lf, Lm, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
				W.NC <- 0.5*W.NC.Female + 0.5*W.NC.Male 
                
					for (j in 1:length(vecMAO)){
						print(vecMAO[j])
						the.Ld  <- 1- (CumRisk(Age, vecMAO[j], vecFAO[j], 0.01))

						# W.males      <- Results.C.male   (FALSE, FALSE, TRUE,  the.Ld, lambda, Ff, Lf, Lm, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
						W.females.1  <- Results.C.female (TRUE, FALSE, FALSE,  the.Ld, lambda, Ff, Lf, Lm, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
						W.females.2  <- Results.C.female (TRUE, TRUE, FALSE,   the.Ld, lambda, Ff, Lf, Lm, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
						
						the.W <- 0.5*W.NC.Male + 0.25*W.females.1 + 0.25*W.females.2  # W.Male== W.Males.NC
						the.S <- 1- the.W/W.NC
						res[j] <- the.S
					}
				resFemaleOnly <- res
			
			# Autosomal SALOD with disease in male only
				res     <- rep(NA, length(vecMAO))
 				the.Ld <- rep(1, length(Age))
				W.NC.Female <- Results.C.female (FALSE, FALSE, FALSE, the.Ld, lambda, Ff, Lf, Lm, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
				W.NC.Male   <- Results.C.male   (FALSE, FALSE, FALSE, the.Ld, lambda, Ff, Lf, Lm, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
				W.NC <- 0.5*W.NC.Female + 0.5*W.NC.Male 
                
					for (j in 1:length(vecMAO)){
						print(vecMAO[j])
						the.Ld  <- 1- (CumRisk(Age, vecMAO[j], vecFAO[j], 0.01))

						W.males      <- Results.C.male   (FALSE, FALSE, TRUE,  the.Ld, lambda, Ff, Lf, Lm, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
						# W.females.1  <- Results.C.female (TRUE, FALSE, FALSE,  the.Ld, lambda, Ff, Lf, Lm, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
						# W.females.2  <- Results.C.female (TRUE, TRUE, FALSE,   the.Ld, lambda, Ff, Lf, Lm, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
						
						the.W <- 0.5*W.males + + 0.5*W.NC.Female # W.female == W.females.NC
						the.S <- 1- the.W/W.NC
						res[j] <- the.S
					}
				resMaleOnly <- res			
			
			# A Mitochondrial SALOD
				res     <- rep(NA, length(vecMAO))
 				the.Ld <- rep(1, length(Age))
				W.NC.Female <- Results.C.female (FALSE, FALSE, FALSE, the.Ld, lambda, Ff, Lf, Lm, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
				# W.NC.Male   <- Results.C.male   (FALSE, FALSE, FALSE, the.Ld, lambda, Ff, Lf, Lm, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
				W.NC <- W.NC.Female 
                
					for (j in 1:length(vecMAO)){
						print(vecMAO[j])
						the.Ld  <- 1- (CumRisk(Age, vecMAO[j], vecFAO[j], 0.01))

						# W.males      <- Results.C.male   (FALSE, FALSE, TRUE,  the.Ld, lambda, Ff, Lf, Lm, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
						W.females.1  <- Results.C.female (TRUE, FALSE, FALSE,  the.Ld, lambda, Ff, Lf, Lm, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
						W.females.2  <- Results.C.female (TRUE, TRUE, FALSE,   the.Ld, lambda, Ff, Lf, Lm, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
						
						the.W <- 0.5*W.females.1 + 0.5*W.females.2 
						the.S <- 1- the.W/W.NC
						res[j] <- the.S
					}
				resMito <- res		

			# A Y chrom SALOD
				res     <- rep(NA, length(vecMAO))
 				the.Ld <- rep(1, length(Age))
				#W.NC.Female <- Results.C.female (FALSE, FALSE, FALSE, the.Ld, lambda, Ff, Lf, Lm, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
				W.NC.Male   <- Results.C.male   (FALSE, FALSE, FALSE, the.Ld, lambda, Ff, Lf, Lm, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
				W.NC <- W.NC.Male
                
					for (j in 1:length(vecMAO)){
						print(vecMAO[j])
						the.Ld  <- 1- (CumRisk(Age, vecMAO[j], vecFAO[j], 0.01))

						W.males      <- Results.C.male   (FALSE, FALSE, TRUE,  the.Ld, lambda, Ff, Lf, Lm, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
						# W.females.1  <- Results.C.female (TRUE, FALSE, FALSE,  the.Ld, lambda, Ff, Lf, Lm, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
						# W.females.2  <- Results.C.female (TRUE, TRUE, FALSE,   the.Ld, lambda, Ff, Lf, Lm, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
						
						the.W <- W.males 
						the.S <- 1- the.W/W.NC
						res[j] <- the.S
					}
				resY <- res		
			
			# THE PLOT

				  age <- seq(40,85,by=2.5) 
				  the.lwd <- 1
				  
				  plot (age, res100, col="white", bty="l", type="b",pch=1, cex.pch=2, log="y", xlab="Mean age at onset (MAO)", ylab="Selection Coefficient", ylim=c(2.5e-07, 1))
				  
				  vage <- c(38,87)
				  # 4Nes<100 => 2.5e-02
				  polygon (c(vage, rev(vage)), c(2.5e-03,2.5e-03,1,1), col="gray60", border=F)
				  # 4Nes<100 => 2.5e-03
				  polygon (c(vage, rev(vage)), c(2.5e-03,2.5e-03,1,1), col="gray60", border=F)
				  # 4Nes=1000 => 2.5e-04
				  polygon (c(vage, rev(vage)), c(2.5e-03,2.5e-03,2.5e-04,2.5e-04), col="gray70", border=F)
				  # 4Nes=10000 => 2.5e-05
				  polygon (c(vage, rev(vage)), c(2.5e-04,2.5e-04,2.5e-05,2.5e-05), col="gray80", border=F)
				  # 4Nes>10000 => 2.5e-05
				  polygon (c(vage, rev(vage)), c(2.5e-05,2.5e-05,2.5e-07,2.5e-07), col="gray90", border=F)

				abline(h=2.5e-03, col="grey", lty=3)
				abline(h=2.5e-04, col="grey", lty=3)
				abline(h=2.5e-05, col="grey", lty=3)
				  
				lines(age, res100		, type="b", pch=1, lwd=2)
				lines(age, resRec		, type="b", pch=3, lwd=2)
				lines(age, resFemaleOnly, type="b", pch=5, lwd=2)
				lines(age, resMaleOnly 	, type="b", pch=2, lwd=2)
				lines(age, resMito		, type="b", pch=8, lwd=2)
				lines(age, resY			, type="b", pch=4, lwd=2)
				




			
			
			