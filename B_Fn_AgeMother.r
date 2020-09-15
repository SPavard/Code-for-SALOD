    # FUNCTIONS FOR DISTRIBUTION PROBABILITY OF CHILDREN'S AGE AT MATERNAL, GRANDMATERNAL AND PATERNAL DEATH
	# AS A FUNCTION OF X1 - AGE OF THE MOTHER AT CHILDBIRTH
	########################################################################################################

			# 1. Function pY1 knowing x1 (corresponding to SM1.2) 
				# Inputs is the adult survival vector of females 
				# y1 is defined between 0 et omega-x1
				# x1 is defined here between 15 et 49 and y1 between 0 et omega-x1-1 
				# return 0 for couple y1 and x1 such that mother dies at age > omega years
				
                pY1<-function(x1, y1, vLf){ 
                  the.y1 <- (vLf[1+x1+y1]-vLf[1+x1+y1+1])/vLf[1+x1]
                  the.y1[is.na(the.y1)==TRUE] <- 0
                  return(the.y1)
                }                                                                                                       
				# plots and sanity checks
                # plot(pY1(35, 0:99, vLf)) ; lines (pY1(25, 0:99, vLf)) ; lines(pY1(15, 0:99, vLf))	# Distribution of child's age at mother death for mother aged 35, 25 and 15 at its birth 
				# sum (pY1(35, 0:99, vLf))
                # sum (pY1(25, 0:99, vLf))
                
			# 2. Function pY2 knowing x1 and Lambda.Females (corresponding to SM2.3)  
				# y2 is defined between -vAlphaf et omega-2*vAlphaf
            
				pY2<-function(x1, y2, lambda, vLf, vFf, vAlphaf, vBetaf, Omega){
                  
					# Because GM was alive at mother's birth, y2 cannot be lower than x1-1 
					if (y2 >= (-x1)){                                    
                        # Range of possible x2
                        x2<-vAlphaf:vBetaf
                        vx2 <- x2+1
                        vx1 <- x1+1
                          
                        # Calculation of mean children survival
                        SMean<-1/(sum(lambda^(-Age+1)*vLf*vFf))
                          
                        # Proportion of mother born to grandmother at age x2 (equation 5 in SM2)
                        prop.x2<-(lambda^(-x2+1))*vLf[vx2]*vFf[vx2]*SMean
                        # sum(prop.x2) = 1
                          
                        # p(y2|x1,x2) => probability of the daughter of losing her grandmother at age y2 (equation 3 in SM2)
                        # Knowing the mother age at her birth x1 and her grandmother age at the birth of her mother x2
                        # Return NA for x1+x2+y2>99
                        y2.knowing.x1.x2  <- (vLf[1+x1+x2+y2]-vLf[1+x1+x2+y2+1])/vLf[1+x2]
                        y2.knowing.x1.x2 [is.na(y2.knowing.x1.x2)] <- 0
                        return(sum(prop.x2*y2.knowing.x1.x2))
                    }
                    else {return(0)}
                }            
				# plots and sanity checks for any value Lambda
				# tmp.y2 <- -49:69 
				# prob.y2.1<-rep(NA,length(tmp.y2))
				# for (i in 1:length(tmp.y2)){prob.y2.1[i]<-pY2(45, tmp.y2[i], lambda, vLf, vFf, vAlphaf, vBetaf, 99)} ; sum(prob.y2.1)
				# prob.y2.2<-rep(NA,length(tmp.y2))
                # for (i in 1:length(tmp.y2)){prob.y2.2[i]<-pY2(25, tmp.y2[i], lambda, vLf, vFf, vAlphaf, vBetaf, 99)} ; sum(prob.y2.2)
                # prob.y2.3<-rep(NA,length(tmp.y2))
                # for (i in 1:length(tmp.y2)){prob.y2.3[i]<-pY2(15, tmp.y2[i], lambda, vLf, vFf, vAlphaf, vBetaf, 99)} ; sum(prob.y2.3)
				# plot (tmp.y2, prob.y2.1, type="l") ; lines(tmp.y2, prob.y2.2) ; lines(tmp.y2, prob.y2.3) # Distribution of child's age at grandmother death for mother aged 45, 25 and 15 at its birth 
               
               
            # 3. Function py3 knowing x1 and Lambda.males (corresponding to SM1.4)
				# y2 is defined between -1 et omega-x1
                
				pY3<-function(x1, y3, lambda, vFm, vLm, vAlpham, vBetam, Omega, matB){            
               
                    # Range of possible age x3
                    x3 <- x1:(vBetam-1)
                    vx3 <- x3+1
                    
                    # Proportion of children born to father at age x3  knowing x1 (SMean is not needed because it cancels out in the ratio; equation SM1.10)
                    prop.x3 <- u.x3.knowing.x1(x1,matB)
                    prop.x3 <- prop.x3[x3+1]
					sum(prop.x3)
                    
                    # p(y3|x3) => probability of losing the father at age y3 knowing the father's age at birth (equation SM1.11)
                    y3.knowing.x3<-(vLm[1+x3+y3]-vLm[1+x3+y3+1])/vLm[1+x3]
                  
                    # Setting p(y3|x3) for impossible occurrence
                    y3.knowing.x3[y3>(Omega-x3-1)]<-0                             
                    y3.knowing.x3[is.na(y3.knowing.x3)==TRUE] <- 0
                    return(sum(prop.x3*y3.knowing.x3))
                }
				# matB <- Union.x1.x3 (Age, dage, Lambda, Lf, Ff, Fm)
				# plots and sanity checks for any value Lambda.males	
                # tmp.y3 <- 0:99 
                # prob.y3.1 <- rep(NA,length(tmp.y3))
                # for (i in 1:length(tmp.y3)){prob.y3.1[i]<-pY3(15, tmp.y3[i], lambdam, vFm, vLm, vAlpham, vBetam, Omega, matB)} ; sum(prob.y3.1)
                # prob.y3.2 <- rep(NA,length(tmp.y3))
                # for (i in 1:length(tmp.y3)){prob.y3.2[i]<-pY3(25, tmp.y3[i], lambdam, vFm, vLm, vAlpham, vBetam, Omega, matB)}
                # plot(tmp.y3,prob.y3.1) ; lines (tmp.y3,prob.y3.2) # Distribution of child's age at father's death for mother aged 45, 25 and 15 at its birth 

            
			# 4. Function P.y1.y2.y3 knowing x1 : Vector of probability for y1, y2, y3 for a given mother's age at child birth x1
				# Harmonizes the possible 3-tuples y1, y2 et y3 and return their probabilities
               
				P.y1.y2.y3.C.females <-function(x1, Mother, GrandM, Father, Ld, lambda, vFf, vLf, vLm, vFm, vAlphaf, vBetaf, vAlpham, vBetam, Omega, matB){
					Y1<- 0:(Omega-vAlphaf)                                                  # possible child's age at mother's death
					Y2<- (-vBetaf):(Omega-2*vAlphaf)                                         # possible child's age at grandmother's death
					Y3<- 0:(Omega-vAlpham)
			  
					tmp         <- rep(1, length(Ld))
					tmp[Mother] <- Ld[Mother]
					Surv.m <- vLf*tmp
					vec.probY1<-rep(NA,length(Y1))
					for (i in 0:length(Y1)){vec.probY1[i]<-pY1(x1, Y1[i], Surv.m)}
					  
					tmp         <- rep(1, length(Ld))
					tmp[GrandM] <- Ld[GrandM]
					Surv.gm <- vLf*tmp
					vec.probY2<-rep(NA,length(Y2))
					for (i in 1:length(Y2)){vec.probY2[i]<-pY2(x1 , Y2[i], lambda, Surv.gm, vFf, vAlphaf, vBetaf)}
				  
					tmp         <- rep(1, length(Ld))
					tmp[Father] <- Ld[Father]
					Surv.f <- vLm*tmp
					vec.probY3<-rep(NA,length(Y3))
					for (i in 1:length(Y3)){vec.probY3[i]<-pY3(x1, Y3[i], lambda, vFm, Surv.f, vAlpham, vBetam, Omega, matB)}    
			  
					TheListe<-list(vy1=Y1,vy2=Y2,vy3=Y3, py1=vec.probY1, py2=vec.probY2, py3=vec.probY3)
					return(TheListe)
				}
			
			# A Numerical calculation  for hunter gatherer where lambda = 1 (Fm  and matB for lambda=1)
			
			par(mfrow=c(2,2))
			essai <- P.y1.y2.y3.C.females(15, TRUE, TRUE, TRUE, rep(1, 100) , lambda, vFf, vLf, vLm, vFm, vAlphaf, vBetaf, vAlpham, vBetam, Omega, matB)
			plot  (essai$vy1 ,essai$py1, col="red", xlim=c(-49, 79), ylim=c(0,0.07), type="l", lty=1, xlab="Child's Age at mat, pat and grandmat death", ylab="probability")
			polygon (c(0,15,15,0),c(0,0,0.07,0.07), col="grey", border=F)
			lines (essai$vy2[essai$py2>0] ,essai$py2[essai$py2>0], col="red", lty=2)
			lines (essai$vy3 ,essai$py3, col="red", lty=3)
			lines (essai$vy1 ,essai$py1, col="red", lty=1)
			title ("x1 = 15 yrs")
			
			essai <- P.y1.y2.y3.C.females(48, TRUE, TRUE, TRUE, rep(1, 100) , lambda, vFf, vLf, vLm, vFm, vAlphaf, vBetaf, vAlpham, vBetam, Omega, matB)
			plot  (essai$vy1 ,essai$py1, col="red", xlim=c(-49, 79), ylim=c(0,0.07), type="l", lty=1, xlab="Child's Age at mat, pat and grandmat death", ylab="probability")
			polygon (c(0,15,15,0),c(0,0,0.07,0.07), col="grey", border=F)
			lines (essai$vy2[essai$py2>0] ,essai$py2[essai$py2>0], col="red", lty=2)
			lines (essai$vy3 ,essai$py3, col="red", lty=3)
			lines (essai$vy1 ,essai$py1, col="red", lty=1)
			title ("x1 = 48 yrs")

			