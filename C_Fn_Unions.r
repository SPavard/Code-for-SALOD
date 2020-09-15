
		###################################################################################################################
		# Function to scale male fertility - Find the Am parameter such thatlamb-x*lx*Fx is identical for males and females.
			solvA <- function (Lambda, Lx, Alpham, Betam, Alpha, Beta, ChosenFertMale, ChosenFertFemale){
				# Attention dage is 1:100 so lambda^x (Daniel Goodman, 1982) and Age is 0:100 	
				p1<- ChosenFertMale[1]
				p2<- ChosenFertMale[2]
				p3<- ChosenFertMale[3]
				p4<- ChosenFertMale[4]
		
				Alpham <- Alpham-1
				  
				X <- sum(Lambda^-dage * Lx * Fert(Age, ChosenFertFemale))
				  
				Z  <- Lambda^-dage * Lx * ((Age-Alpham)*((Betam-Age)^2)*((p2*Age)+(p3*(Age^2))+(p4*(Age^3))))  
				Z[Age<Alpham]<-0                                                     
				Z[Age>Betam] <-0
				Z<-sum(Z)
				  
				Y  <- Lambda^-dage * Lx *((Age-Alpham)*((Betam-Age)^2))
				Y[Age<Alpham]<-0                                                     
				Y[Age>Betam] <-0
				Y<-sum(Y)
				  
				the.y<-(X - Z)/Y
				return(the.y)
          }
		
		###################################################################################################################
		# Function UNION (x, x3) - Provide the table of the frequency of birth of children from mother at age x1 and father at age x3
			Union.x1.x3 <- function (Age, dage, Lambda, vLf, vLm, vFf, vFm){
				distf <-  Lambda^-dage * vLf * vFf
				distm <-  Lambda^-dage * vLf * vFm
				matbirth <- matrix (0, length(Age),	length(Age))	
					for (i in Age){
						restf <- distf[i]
						for (j in 1:100){
							restm         <- distm[j] - sum(matbirth[j,])
							matbirth[j,i] <- min(restm, restf)
							restf         <- restf-matbirth[j,i]
						}			
					}
				return(matbirth)	 
			}
		
		###################################################################################################################
		# Function u.x3.knowing.x1 - Provide the vector of probability that an husband is x3 if wife is x1
			u.x3.knowing.x1 <- function(x1, matB){
				y <- matB[,x1+1]/sum(matB[,x1+1])
				y[is.na(y)==TRUE]<-0
				return(y)
			}
			# Attention: Correspondence with Age (from 0 to 99)
			
		# Function u.x1.knowing.x3 - Provide the vector of probability that a wife is x1 if husband is x3
			u.x1.knowing.x3 <- function(x3, matB){
				y <- matB[x3+1,]/sum(matB[x3+1,])
				y[is.na(y)==TRUE]<-0
				return(y)
			}
			# Attention: Correspondence with Age (from 0 to 99)
			
		###################################################################################################################
		# Sanity check and graphs	
				# Sanity Check of the SolveA function - CAlculatin of Am if Lambda = 1.01
				Lambda <- 1.0
				newParAmales <- solvA(Lambda, Lf, Alpham, Betam, Alpha, Beta, ChosenFertMale, ChosenFertFemale)
				newParFertMales <- c(newParAmales, ChosenFertMale[2:6])
				Fm <- Fert(Age, newParFertMales)
				sum(Lambda^-dage * Lf * Ff)
				sum(Lambda^-dage * Lf * Fm)
				# it seems to work well...
			
		# Making a nice plot of it
			#matB <- Union.x1.x3 (Age, dage, Lambda, Lf, Ff, Fm)
			
			#layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
			#plot (colSums(matB), type="s", xlim=c(15, 55), xlab="Female age", ylab="Female birth rate") 
			#palcol <- palette(gray(seq(0.05 ,0.95,len = 60)))
			#colfem <- rgb(t(col2rgb("firebrick")), max = 255, alpha=valalpha<-150)       
			#library(colorspace)
			#palcol <- diverge_hcl(55,h=c(246,40), c=96) # of library colorspace
			
			#for (i in (Alphaf+1):(Betaf)){
			#	print(i)
		#		vec  <- matB[,i]
			#	lvec <-	length(which(vec!=0))
			#	y    <- vec[which(vec!=0)] 
			#	y1 <- c(0, cumsum(y)[-length(y)])
			#	y2 <- cumsum(y)
			#	agem <- Age [which(vec!=0)] - Alpham +1
			#	print(agem)
			#	rect(i, y1, i+1, y2, col=palcol[agem])
			#}		
		
			#lines(Lambda^-dage * Lf * Fert(Age, ChosenFertFemale), col=colfem, type="s", lwd=2) 
			#legend_image <- as.raster(matrix(rev(palcol), ncol=1)) 
			#plot(c(0,2), c(Alpham,Betam), type = 'n', axes = F, xlab = '', ylab = '', main = 'Male age')
			#text(x=1.5, y = Alpham:Betam, labels = Alpham:Betam, cex=0.5)
			#rasterImage(legend_image, 0, Alpham, 1.25,Betam) 

