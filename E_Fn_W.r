

		######################################################################################
		# Function Results Carriers : Return the reproductive success of female carrier
		# The Function Solver should have been run for finding Lambda, vFm and matB
				
				Results.C.female <- function (Mother, GrandM, Father, Ld, lambda, vFf, vLf, vLm, vFm, vAlphaf, vBetaf, vAlpham, vBetam, Omega, vecSigma, matB){
					# Possible x1
					x1<-vAlphaf:vBetaf  
					
					# The S.Alpha.Knowing.x1 
					Salphax1<-rep(NA,length(x1))
						for (i in 1:length(x1)){
						  thelist     <- P.y1.y2.y3.C.females(x1[i], Mother, GrandM, Father, Ld, lambda, vFf, vLf, vLm, vFm, vAlphaf, vBetaf, vAlpham, vBetam, Omega, matB)
						  Salphax1[i] <- S.Alpha(thelist, vecSigma)
						}
					
					#Salpha.0100  
					Salpha.0100 <- rep(0, length(Age))
					tmp <- match((vAlphaf:vBetaf),Age)
					Salpha.0100[tmp]<- Salphax1
					
					# The Wc for females
					Wc.females <- sum(vLf*vFf*Ld*Salpha.0100)
					return(Wc.females)
				} 
				
				# lambda <- 1.00432777404785
			    #newParAmales <- solvA (lambda, Lf, Alpham, Betam, Alpha, Beta, ChosenFertMale, ChosenFertFemale)
				#newParFertMales <- c(newParAmales, ChosenFertMale[2:6])
				#the.Fm <- Fert(Age, newParFertMales)
				#matB <- Union.x1.x3 (Age, dage, lambda, Lf, Ff, the.Fm)
				#Results.C.female(TRUE, TRUE, TRUE, rep(1, 100), lambda, Ff, Lf, Lm, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
				
				Results.C.male <- function (Mother, GrandM, Father, Ld, lambda, vFf, vLf, vLm, vFm, vAlphaf, vBetaf, vAlpham, vBetam, Omega, vecSigma, matB){
					# Possible x1
					x3<-vAlpham:vBetam  
					
					# The S.Alpha.Knowing.x1 
					Salphax3<-rep(NA,length(x3))
						for (i in 1:length(x3)){
						  thelist     <- P.y1.y2.y3.C.males(x3[i], Mother, GrandM, Father, Ld, lambda, vFf, vLf, vLm, vFm, vAlphaf, vBetaf, vAlpham, vBetam, Omega, matB)
						  Salphax3[i] <- S.Alpha(thelist, vecSigma)
						}
									
					#Salpha.0100 
					Salpha.0100 <- rep(0, length(Age))
					tmp <- match((vAlpham:vBetam),Age)
					Salpha.0100[tmp]<- Salphax3
					
					# The Wc for males
					Wc.males <- sum(vLm*vFm*Ld*Salpha.0100)
					return(Wc.males)
				}
				
				#Results.C.male(TRUE, TRUE, TRUE, rep(1, 100), lambda, Ff, Lf, Lm, the.Fm, Alphaf, Betaf, Alpham, Betam, Omega, vecSigma, matB)
				#plot(Ff) ; lines(the.Fm)
				#sum(Age*vLf*vFf)/sum(vLf*vFf)                  # Generation time of females
				#sum(Age*vLm*the.Fm)/sum(vLm*the.Fm)            # Generation time of males

					