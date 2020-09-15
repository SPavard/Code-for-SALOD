
	#INPUTS REQUIRED: Age, dAge, Lf, Lm, Ff, ChosenFertFemale, ChosenFertMale, vecSigma, Ld
	
	# SOLVER FUNCTION	
	
	Solver<-function(niteration, vFf, vLf, vLm, vAlphaf, vBetaf, vAlpham, vBetam, Omega, vecSigma, ChosenFertMale, ChosenFertFemale){
  
			nit <- niteration                                  # Number of iteration to slove lambda
			thelambdas <- rep(NA, nit)                         # Empty vector to fill iterated lambda
			theones    <- rep(NA, nit)                         # Empty vector to check convergence to one of the Euler sum
			theW	   <- rep(NA, nit)
			thePopS    <- rep(NA, nit)
			theTrajSf  <- rep(NA, nit)
			
			thelambdas[1] <- 1                                # Initial females lambda value
			up<-2                                             # Initial step up
			low<-0                                            # Initial step down

		for (i in 1:nit){
			print(i)
			mid    <- low+(up-low)/2
			lambda <- mid
			
			# For one lambda : Solving males fertility trajectories
				newParAmales    <- solvA (lambda, vLf, vAlpham, vBetam, vAlpha, vBeta, ChosenFertMale, ChosenFertFemale)
				newParFertMales <- c(newParAmales, ChosenFertMale[2:6])
				the.Fm <- Fert(Age, newParFertMales)
			
			# For one vFm :  Finding the matrix of unions probabilities 
				matB <- Union.x1.x3 (Age, dage, lambda, vLf, vLm, vFf, the.Fm)
			
			# From matB and Everything else : Solving the Salpha vector
				x1<-vAlphaf:vBetaf
				Salphax1<-rep(NA,length(x1))
							
					for (j in 1:length(x1)){
								thelist    <- P.y1.y2.y3.C.females(x1[j], TRUE, TRUE, TRUE, rep(1,100), lambda, vFf, vLf, vLm, the.Fm, vAlphaf, vBetaf, vAlpham, vBetam, Omega, matB)
								Salphax1[j]<- S.Alpha(thelist, vecSigma)
					}
			
			# Buiding a vector Salpha from age 0 to 99 => Salpha.0100 
				Salpha.0100 <- rep(0, length(Age))
				tmp <- match((vAlphaf:vBetaf),Age)
				Salpha.0100[tmp]<- Salphax1
			
			# It is the good lambda?
				is.one <- sum(lambda^(-dage)*vLf*vFf*Salpha.0100)
				theones[i] <- is.one
				thelambdas[i] <- lambda
				low[is.one>1] <- mid
				up [is.one<1] <- mid

			# Reproductive success
				theW[i]  <-sum(vLf*vFf*Salpha.0100)
			
			# pop mean SAlpha
				thePopS[i]   <- 1/sum(lambda^(-dage)*vLf*vFf)
			
			# traj mean SAlpha
				theTrajSf[i]  <- sum(vLf*vFf*Salpha.0100)/sum(vLf*vFf)
				
			}

		# Return Values: The convergence of lambda, the convergence of the Euler Sum, the last selective value, the S.x1 trajectories, the last mean children survival from lambda and lambda.males,
		theResults<-list(conv.Lambda=thelambdas , conv.one=theones, Repro.value=theW, Last.S.Alpha.x=Salpha.0100, PopS=thePopS, TrajS=theTrajSf, Last.Fm=the.Fm)
		return(theResults)
	}	



	
	
	
	