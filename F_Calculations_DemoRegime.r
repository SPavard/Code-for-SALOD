			####################################################################################################################################################################
			####################################################################################################################################################################			
			Fvariance<-function(MAO, FAO, Threshold){
				  y<--(FAO-MAO)/(log(1/Threshold-1))
				  return(y)
				}
				
			CumRisk<-function(Age, MAO, FAO, Threshold){
				  y<-plogis(Age, location = MAO, scale = Fvariance(MAO, FAO, Threshold))
				  y[Age<FAO]<-0
				  return(y)
				}
			
			CoefSel <- function (vecMAO, vecFAO, vecSurv){
                res     <- rep(NA, length(vecMAO))
                
                the.Ld <- rep(1, length(Age))
                W.NC.Female <- Results.C.female (FALSE, FALSE, FALSE, the.Ld, lambda, Ff, vecSurv, vecSurv, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
				W.NC.Male   <- Results.C.male   (FALSE, FALSE, FALSE, the.Ld, lambda, Ff, vecSurv, vecSurv, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
				#scaling <- W.NC.Female/W.NC.Male
				W.NC <- 0.5*W.NC.Female + 0.5*W.NC.Male 
								
				for (j in 1:length(vecMAO)){
				print(vecMAO[j])
				the.Ld  <- 1- CumRisk(Age, vecMAO[j], vecFAO[j], 0.01)

				W.males      <- Results.C.male   (FALSE, FALSE, TRUE,  the.Ld, lambda, Ff, vecSurv, vecSurv, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
				W.females.1  <- Results.C.female (TRUE, FALSE, FALSE,  the.Ld, lambda, Ff, vecSurv, vecSurv, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
				W.females.2  <- Results.C.female (TRUE, TRUE, FALSE,   the.Ld, lambda, Ff, vecSurv, vecSurv, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)

				the.W <- 0.5*W.males + 0.25*W.females.1 + 0.25*W.females.2 
				the.S <- 1- the.W/W.NC
				res[j] <- the.S
				}
                return(res)
				}
			
			# RESULTS HUNTER-GATHERER
			pop <- Solver(20, Ff, Lf, Lm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, ChosenFertMale, ChosenFertFemale)
					lambda <- pop$conv.Lambda[20]
					newParAmales <- solvA (lambda, Lf, Alpham, Betam, Alphaf, Betaf, ChosenFertMale, ChosenFertFemale)
					newParFertMales <- c(newParAmales, ChosenFertMale[2:6])
					the.Fm <- Fert(Age, newParFertMales)
					matB <- Union.x1.x3 (Age, dage, lambda, Lf, Lm, Ff, the.Fm)				
			lambda	

			the.MAO <- seq(40,85,by=2.5) 
			the.FAO <- the.MAO - 20
			res100 <- CoefSel(the.MAO, the.FAO, Lf)
			res100
			
			CoefSel(40, 40-0.5, Lf)
			
			# LOADING MORTALITY OF SWEDEN
			# SWEDEN SURVIVAL
			LifeT_Sweden <- read.table("C:/Users/Pavard/Documents/LaNouvelleSouille/P_Selection_Old_Age_NouvelleSouille2/Demo_Sweden/Sueden_Cohort_F3.csv", sep = ";", header=TRUE)

			L1751 <-LifeT_Sweden$lx[LifeT_Sweden$Year==1751] 
			L1751 <-L1751/L1751[Alphaf]    ;  L1751[L1751>1]<-1 

			L1800 <-LifeT_Sweden$lx[LifeT_Sweden$Year==1800] 
			L1800 <-L1800/L1800[Alphaf]    ;  L1800[L1800>1]<-1 
			
			L1850 <-LifeT_Sweden$lx[LifeT_Sweden$Year==1850] 
			L1850 <-L1850/L1850[Alphaf]    ;  L1850[L1850>1]<-1 
			
			L1900 <-LifeT_Sweden$lx[LifeT_Sweden$Year==1900] 
			L1900 <-L1900/L1900[Alphaf]    ;  L1900[L1900>1]<-1 
			
			plot (Age, Lf, type="l")
			lines(Age, L1900, type="l")
			
			# RES1900 
			pop <- Solver(20, Ff, L1900,  L1900, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, ChosenFertMale, ChosenFertFemale)
					lambda <- pop$conv.Lambda[20]
					newParAmales <- solvA (lambda, L1900, Alpham, Betam, Alphaf, Betaf, ChosenFertMale, ChosenFertFemale)
					newParFertMales <- c(newParAmales, ChosenFertMale[2:6])
					the.Fm <- Fert(Age, newParFertMales)
					matB <- Union.x1.x3 (Age, dage, lambda, L1900, L1900, Ff, the.Fm)		
					lambda
			1.006186
			
			CoefSel(85, 65, L1900)

				the.Ld <- rep(1, length(Age))
                W.NC.Female <- Results.C.female (FALSE, FALSE, FALSE, the.Ld, lambda, Ff, L1900, L1900, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
				W.NC.Male   <- Results.C.male   (FALSE, FALSE, FALSE, the.Ld, lambda, Ff, L1900, L1900, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
				W.NC        <- 0.5*W.NC.Female + 0.5*W.NC.Male # 1.137374
			
				the.Ld  <- 1- CumRisk(Age, 85, 65, 0.01)
				W.males      <- Results.C.male   (FALSE, FALSE, TRUE,  the.Ld, lambda, Ff, L1900, L1900, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB) # 1.156511
				W.females.1  <- Results.C.female (TRUE, FALSE, FALSE,  the.Ld, lambda, Ff, L1900, L1900, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
				W.females.2  <- Results.C.female (TRUE, TRUE, FALSE,   the.Ld, lambda, Ff, L1900, L1900, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
 
				W.NC.Male > W.males 
			
			res <- CoefSel(the.MAO, the.FAO, L1900)
			res1900 <-res
			
			
					x1<-Alphaf:Betaf  
					
					# The S.Alpha.Knowing.x1 
					Salphax1<-rep(NA,length(x1))
						for (i in 1:length(x1)){
						  thelist     <- P.y1.y2.y3.C.females(x1[i],FALSE, FALSE, FALSE, rep(1, length(Age)), lambda, Ff, L1900, L1900, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, matB)
						  Salphax1[i] <- S.Alpha(thelist, vecSigma)
						}
					
					# The S.Alpha.Knowing.x1 
					Salphax1<-rep(NA,length(x1))
						for (i in 1:length(x1)){
						  thelist     <- P.y1.y2.y3.C.females(x1[i],TRUE, FALSE, FALSE, the.Ld, lambda, Ff, L1900, L1900, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, matB)
						  Salphax1[i] <- S.Alpha(thelist, vecSigma)
						}
					lines(Salphax1)
					
					
					
					
					#Salpha.0100  
					Salpha.0100 <- rep(0, length(Age))
					tmp <- match((vAlphaf:vBetaf),Age)
					Salpha.0100[tmp]<- Salphax1
					
					# The Wc for females
					Wc.females <- sum(vLf*vFf*Ld*Salpha.0100)
			
			
			
			
			
			
			
			
			
			# RES1800 
			pop <- Solver(20, Ff, L1800,  L1800, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, ChosenFertMale, ChosenFertFemale)
					lambda <- pop$conv.Lambda[20]
					newParAmales <- solvA (lambda, L1800, Alpham, Betam, Alphaf, Betaf, ChosenFertMale, ChosenFertFemale)
					newParFertMales <- c(newParAmales, ChosenFertMale[2:6])
					the.Fm <- Fert(Age, newParFertMales)
					matB <- Union.x1.x3 (Age, dage, lambda, L1800,L1800, Ff, the.Fm)		
					lambda
			1.005495
			res <- CoefSel(the.MAO, the.FAO)
			res1800 <-res
			
			# RES1850 
			pop <- Solver(20, Ff, L1850,  L1850, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, ChosenFertMale, ChosenFertFemale)
					lambda <- pop$conv.Lambda[20]
					newParAmales <- solvA (lambda, L1850, Alpham, Betam, Alphaf, Betaf, ChosenFertMale, ChosenFertFemale)
					newParFertMales <- c(newParAmales, ChosenFertMale[2:6])
					the.Fm <- Fert(Age, newParFertMales)
					matB <- Union.x1.x3 (Age, dage, lambda, L1850,L1850, Ff, the.Fm)		
					lambda
			1.006445
			res <- CoefSel(the.MAO, the.FAO)
			res1850 <-res
			
			# RES1900 
			pop <- Solver(20, Ff, L1900,  L1900, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, ChosenFertMale, ChosenFertFemale)
					lambda <- pop$conv.Lambda[20]
					newParAmales <- solvA (lambda, L1900, Alpham, Betam, Alphaf, Betaf, ChosenFertMale, ChosenFertFemale)
					newParFertMales <- c(newParAmales, ChosenFertMale[2:6])
					the.Fm <- Fert(Age, newParFertMales)
					matB <- Union.x1.x3 (Age, dage, lambda, L1900, L1900, Ff, the.Fm)		
					lambda
			1.007032
			res <- CoefSel(the.MAO, the.FAO, L1900)
			res1900 <-res
			
			the.Ld <- rep(1, length(Age))
            W.NC.Female <- Results.C.female (FALSE, FALSE, FALSE, the.Ld, lambda, Ff, L1900, L1900, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
			W.NC.Male   <- Results.C.male   (FALSE, FALSE, FALSE, the.Ld, lambda, Ff, L1900, L1900, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
			W.NC <- 0.5*W.NC.Female + 0.5*W.NC.Male 
			the.Ld  <- 1- CumRisk(Age, 70, 20, 0.01)
			W.males      <- Results.C.male   (FALSE, FALSE, TRUE,  the.Ld, lambda, Ff, L1900,  L1900, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
			W.females.1  <- Results.C.female (TRUE, FALSE, FALSE,  the.Ld, lambda, Ff, L1900, L1900, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
			W.females.2  <- Results.C.female (TRUE, TRUE, FALSE,   the.Ld, lambda, Ff, L1900, L1900, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
			the.W <- 0.5*W.males + 0.25*W.females.1 + 0.25*W.females.2 
			1- the.W/W.NC
				
			
			the.Ld <- rep(1, length(Age))
            W.NC.Female <- Results.C.female (FALSE, FALSE, FALSE, the.Ld, lambda, Ff, Lf, Lf, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
			W.NC.Male   <- Results.C.male   (FALSE, FALSE, FALSE, the.Ld, lambda, Ff, Lf, Lf, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
			W.NC <- 0.5*W.NC.Female + 0.5*W.NC.Male 
			the.Ld  <- 1- CumRisk(Age, 70, 20, 0.01)
			W.males      <- Results.C.male   (FALSE, FALSE, TRUE,  the.Ld, lambda, Ff, Lf, Lf, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
			W.females.1  <- Results.C.female (TRUE, FALSE, FALSE,  the.Ld, lambda, Ff ,Lf, Lf, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
			W.females.2  <- Results.C.female (TRUE, TRUE, FALSE,   the.Ld, lambda, Ff, Lf, Lf, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
			the.W <- 0.5*W.males + 0.25*W.females.1 + 0.25*W.females.2 
			1- the.W/W.NC
				
			
			# RES1751 
			L1751[100] <- 0
			pop <- Solver(20, Ff, L1751,  L1751, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, ChosenFertMale, ChosenFertFemale)
					lambda <- pop$conv.Lambda[20]
					newParAmales <- solvA (lambda, L1751, Alpham, Betam, Alphaf, Betaf, ChosenFertMale, ChosenFertFemale)
					newParFertMales <- c(newParAmales, ChosenFertMale[2:6])
					the.Fm <- Fert(Age, newParFertMales)
					matB <- Union.x1.x3 (Age, dage, lambda, L1751, L1751, Ff, the.Fm)		
					lambda
			1.004118
			res <- CoefSel(the.MAO, the.FAO, L1751)
			res1751 <-res
			CoefSel(85, 65, L1751)
				
			
			# THE PLOT

			age <- seq(40,85,by=2.5) 
			par(mfrow=c(1,2))

			plot (Age, Lf   , type="b", pch=1 )
			lines(Age, L1751, type="b", pch=3 )
			lines(Age, L1800, type="b", pch=5 )
			lines(Age, L1850, type="b", pch=2 )
			lines(Age, L1900, type="b", pch=8 )	


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
				lines(age, res1751		, type="b", pch=3, lwd=2)
				lines(age, res1800      , type="b", pch=5, lwd=2)
				lines(age, res1850  	, type="b", pch=2, lwd=2)
				lines(age, res1900		, type="b", pch=8, lwd=2)
