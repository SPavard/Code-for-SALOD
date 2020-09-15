    # FUNCTIONS FOR DISTRIBUTION PROBABILITY OF CHILDREN'S AGE AT MATERNAL, GRANDMATERNAL AND PATERNAL DEATH
	# AS A FUNCTION OF X3 - AGE OF THE FATHER AT CHILDBIRTH
	########################################################################################################

        # 1. Function pY3. males - Probability for the father to die when child is y3 while he was x3 at its birth
				pY3.males <- function(x3, y3, vLm){ 
					the.y3  <- (vLm[1+x3+y3]-vLm[1+x3+y3+1])/vLm[1+x3]
					the.y3[is.na(the.y3)==TRUE] <- 0
					return(the.y3)
				}      
				# plots and sanity checks
				# plot(pY3.males(45, 0:99, vLm)) ; lines (pY3.males(30, 0:99, vLm)) ; lines(pY3.males(20, 0:99, vLm))
				# sum (pY3.males(35, 0:99, vLm))
				# sum (pY3.males(25, 0:99, vLm))
      
        # 2. Function pY1. males - Probability for the mother to die when child is y1 while the father was x3 at its birth
				pY1.males <- function(x3, y1, lambda, vFf, vLf, vAlphaf, Omega, matB){            
				  
					# Range of possible age x1
					x1<-vAlphaf:x3
					vx1 <- x1+1
					  
					# Proportion of children born to mother at age x1  (SMean is not needed because of the ratio it cancels out)
					prop.x1 <- u.x1.knowing.x3(x3, matB)
					prop.x1 <- prop.x1[x1+1]
					
					sum(prop.x1)

					# p(y1|x1) => probability of losing the Mother at age y1 knowing the mother's age at birth x1
					y1.knowing.x1<-(vLf[1+x1+y1]-vLf[1+x1+y1+1])/vLf[1+x1]
					  
					# Setting p(y1|x1) for impossible occurence
					y1.knowing.x1[y1>(Omega-x1-1)]<-0                            
					y1.knowing.x1[is.na(y1.knowing.x1)==TRUE] <- 0
					return(sum(prop.x1*y1.knowing.x1))
				}
				# plots and sanity checks
				# tmp.y1 <- 0:99 
				# prob.y1.1 <- rep(NA,length(tmp.y1))
				# for (i in 1:length(tmp.y1)){prob.y1.1[i]<-pY1.males(49, tmp.y1[i], lambda, vFf, vLf, vAlphaf, Omega, matB)} ; sum(prob.y1.1)
				# prob.y1.2 <- rep(NA,length(tmp.y1))
				# for (i in 1:length(tmp.y1)){prob.y1.2[i]<-pY1.males(25, tmp.y1[i], lambda, vFf, vLf, vAlphaf, Omega, matB)} ;  sum(prob.y1.2)
				# plot(tmp.y1,prob.y1.1); lines (tmp.y1,prob.y1.2)
				# cbind(tmp.y1,prob.y1.2)
				

		# 3. Function P.y1.y2.y3.C.males knowing x3 - Prob y1, y2, y3 for males carriers
			# The probability for the maternal grandmother (the males mother's in law) to die at child's age y2 knowing x3 is calculated there.
			# This is because this probability depends on the mother (the males' wife age) at the childbirth x1: 
				# if the males wife is young, its mother in law is more likely being alive.
				# if the males wife is old, its mother in law is more likely being dead.
				# for each males age e have therefore to integrate over all possible wife's age.
			
				P.y1.y2.y3.C.males <-function(x3, Mother, GrandM, Father, Ld, lambda, vFf, vLf, vLm, vFm, vAlphaf, vBetaf, vAlpham, vBetam, Omega, matB){
					Y1<- 0:(Omega-vAlphaf)                                                  # possible child's age at mother's death
					Y2<- (-vBetaf):(Omega-2*vAlphaf)                                        # possible child's age at grandmother's death
					Y3<- 0:(Omega-vAlpham)
					
					# Prob. y3 as a function of x3
					tmp         <- rep(1, length(Ld))
					tmp[Father] <- Ld[Father]
					Surv.f <- vLm*tmp
					vec.probY3<-rep(NA,length(Y3))
					for (i in 1:length(Y3)){vec.probY3[i]<- pY3.males(x3,Y3[i], Surv.f)}    
				   
					# Prob. y1 as a function of x3
					tmp         <- rep(1, length(Ld))
					tmp[Mother] <- Ld[Mother]
					Surv.m <- vLf*tmp
					vec.probY1<-rep(NA,length(Y1))
					for (i in 0:length(Y1)){vec.probY1[i] <- pY1.males(x3, Y1[i], lambda, vFf, Surv.m, vAlphaf, Omega, matB)}
				
				
					# Prob. y2 as a function of x3 = p(x1|x3)p(x2|x1)p(y2|x1,x2)
					tmp         <- rep(1, length(Ld))
					tmp[GrandM] <- Ld[GrandM]
					Surv.gm <- vLf*tmp

					# Probability of x1 knowing x3
						# Range of possible age x1
						x1<-vAlphaf:x3
						vx1 <- x1+1

						# Proportion of children born to a mother at age x1 knowing x3 (SMean is not needed because of the ratio it cancels out)
						prop.x1 <- u.x1.knowing.x3(x3, matB)
						prop.x1 <- prop.x1[x1+1]
						sum(prop.x1)
						
						# Proportion of children loosing their grandmother at a given age y2 knowing x1
						tres <- matrix(NA, ncol=length(vx1), nrow = length(Y2))
						for (i in 1:length(vx1)){
						  vec.probY2<-rep(NA,length(Y2))
						  for (j in 1:length(Y2)){
							vec.probY2[j]<-pY2(x1[i], Y2[j], lambda, Surv.gm, vFf, vAlphaf, vBetaf)
							tres[,i] <- vec.probY2
						  }
						}
						# Proportion of children loosing their grandmother at a given age y2
						vec.probY2 <- as.vector(t(tres%*%prop.x1))
						
				TheListe<-list(vy1=Y1,vy2=Y2,vy3=Y3, py1=vec.probY1, py2=vec.probY2, py3=vec.probY3)
				return(TheListe)
			  }
			  
      		# Sanity check
			# thelist<- P.y1.y2.y3.C.males(70, TRUE, TRUE, TRUE, rep(1, 100), lambda, vFf, vLf, vLm, vFm, vAlphaf, vBetaf, vAlpham, vBetam, Omega, matB)
			# thelist<- P.y1.y2.y3.C.males(20, TRUE, TRUE, TRUE, Ld, lambda, vFf, vLf, vLm, vFm, vAlphaf, vBetaf, vAlpham, vBetam, Omega)
      
      
			# A Numerical calculation  for hunter gatherer where lambda = 1 (Fm  and matB for lambda=1)
			
			par(mfrow=c(2,1))
			essai <- P.y1.y2.y3.C.males(20, TRUE, TRUE, TRUE, rep(1, 100) , lambda, vFf, vLf, vLm, vFm, vAlphaf, vBetaf, vAlpham, vBetam, Omega, matB)
			plot  (essai$vy1 ,essai$py1, col="blue", xlim=c(-49, 79), ylim=c(0,0.07), type="l", lty=1, xlab="Child's Age at mat, pat and grandmat death", ylab="probability")
			polygon (c(0,15,15,0),c(0,0,0.07,0.07), col="grey", border=F)
			lines (essai$vy2[essai$py2>0] ,essai$py2[essai$py2>0], col="blue", lty=2)
			lines (essai$vy3 ,essai$py3, col="blue", lty=3)
			lines (essai$vy1 ,essai$py1, col="blue", lty=1)
			title ("x3 = 20 yrs")
			
			essai <- P.y1.y2.y3.C.males(69, TRUE, TRUE, TRUE, rep(1, 100) , lambda, vFf, vLf, vLm, vFm, vAlphaf, vBetaf, vAlpham, vBetam, Omega, matB)
			plot  (essai$vy1 ,essai$py1, col="blue", xlim=c(-49, 79), ylim=c(0,0.07), type="l", lty=1, xlab="Child's Age at mat, pat and grandmat death", ylab="probability")
			polygon (c(0,15,15,0),c(0,0,0.07,0.07), col="grey", border=F)
			lines (essai$vy2[essai$py2>0] ,essai$py2[essai$py2>0], col="blue", lty=2)
			lines (essai$vy3 ,essai$py3, col="blue", lty=3)
			lines (essai$vy1 ,essai$py1, col="blue", lty=1)
			title ("x3 = 69 yrs")
