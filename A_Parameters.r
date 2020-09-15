    # INPUTS AND PARAMETERS
	####################################################################
		# SURVIVAL AND FERTILITY
		  
			# PARAMETERS FOR SURVIVAL AND FERTILITY	
	
				# For averaged hunter-gatherer population 
				Par.Surv.HG 			<- c(0.422, 1.131, 0.013, 1.47E-04, 0.086) 
				Par.Fert.HG.Females 	<- c(2.32171630793466e-04, -2.10766279707403e-05,  6.67518316815165e-07 , -6.83478862092560e-09, 15, 49)
				Par.Fert.HG.Males   	<- c(-3.66757388064390e-06,  5.86546994301002e-07, -1.49909448171934e-08,  1.19628084499445e-10, 15, 70) # Parameter A fine-tuned for HG mortality 
	
	
				# For 1751 Sweden 
				Par.Surv.Swed         	<- c(0.2973080783, 0.8127970813, 0.0076763786, 0.0001888114, 0.0841270685)
				Par.Fert.Swed.Females 	<- c(4.64343220558943e-04, -4.21532510259519e-05,  1.33503644755059e-06, -1.36695749316199e-08, 15, 49)
				Par.Fert.Swed.Males   	<- c(4.927266e-05, -2.47197225583200e-06,  5.04667113890550e-08, -3.50481066261845e-10, 20, 70)         # Parameter A fine-tuned for Sweden mortality 
				
			# PARAMETERS FOR MATERNAL, GRANDMATERNAL AND PATERNAL INVESTMENT	
				# entries 1-3 correspond to the parameters of the best fitted relative risk function for maternal care (ay^b+1+c), 4 and 5 for the uniform RR function of grandmaternall and paternal care
				par.invest <- c(13.6691,  0.7086,  0.5678, 1.098, 1.763) 
				
		    # DEMOGRAPHIC FUNCTIONS FOR POPULATION SURVIVAL - SILER MODEL 
          
				# POPULATION MORTALITY HAZARD - SILER MODEL                          # Inputs: AGE + vector a1, b1, a2, a3, b3
				hx<-function(Age, par.Surv){
				a1<-par.Surv[1] ; b1<-par.Surv[2]; a2<-par.Surv[3] ; a3<-par.Surv[4] ; b3<-par.Surv[5] ;
				y<-a1*exp(-b1*Age)+a2+a3*exp(b3*Age)
				y[Age>=99]<-10000                                                    # Omega makeS everybody dying at age 99 in discrete calculation
				return(y)
				}
				  
				# POPULATION SURVIVAL PROBABILITIES - SILER MODEL                    # Inputs: AGE + vector a1, b1, a2, a3, b3
				Sx<-function(Age, par.Surv){
				a1<-par.Surv[1] ; b1<-par.Surv[2]; a2<-par.Surv[3] ; a3<-par.Surv[4] ; b3<-par.Surv[5] ;
				y<-exp((a1/b1)*(exp(-b1*Age)-1)) * exp(-a2*Age) * exp((a3/b3)*(1-exp(b3*Age)))
				y[Age>=99]<-0                                                        # Attention: Omega here to make everybody diying at age 99 in discrete calculation
				return(y)
				}

			# DEMOGRAPHIC FUNCTIONS FOR FERTILITY - BRASS POLYNOMIAL
				Fert<-function(Age, ParFert){                                        # Inputs: AGE + vector p1, p2, p3, p4, Alpha, Beta]
				p1  <-ParFert[1]                                                     # ATTENTION THE BRASS FUNCTION EQUAL 0 IN BETA AND ALPHA
				p2  <-ParFert[2]
				p3  <-ParFert[3]
				p4  <-ParFert[4]          
				ALPHA <-ParFert[5]-1                                                 # Age at first birth recorded
				BETA  <-ParFert[6]
					
				# Age at last birth recorded (menopause-1)
				y<-(Age-ALPHA)*((BETA-Age)^2)*(p1+(p2*Age)+(p3*(Age^2))+(p4*(Age^3)))
				y[Age<ALPHA]  <-0                                                    # The function is fitted for fecundity =0 before age ALPHA
				y[Age>=(BETA)]<-0                                                    # The function is fitted for fecundity =0 after age BETA
				return(y)
				}

			# CHOSEN INPUTS PARAMETERS (here we chose Hunter-Gatherers)
			  ChosenSurvMale    <- Par.Surv.HG                                  # Adult Survival - Males
				ChosenSurvFemale  <- Par.Surv.HG                                  # Adult Survival - Females
				ChosenFertFemale  <- Par.Fert.HG.Females                          # Adult Fecundity - Females
				ChosenFertMale    <- Par.Fert.HG.Males                            # Adult Fecundity - Males
	      
				Omega  <- 99
				Age    <- 0:Omega
				dage   <- 1:100
				Alphaf <- ChosenFertFemale[5]                                # First reproduction between age 15 and 16
				Betaf  <- ChosenFertFemale[6]                                # Last reproduction between  age 49 and 50
				Alpham <- ChosenFertMale[5]                                  # First reproduction between age 20 and 21
				Betam  <- ChosenFertMale[6]                                  # Last reproduction between  age 59 and 60

			    # Adults survival
				Lf <- Sx(Age, ChosenSurvFemale)     ;  Lf<-Lf/Lf[Alphaf]    ;  Lf[Lf>1]<-1   # Vector of survival probabilities for female surviving at age Alphaf
				Lm <- Sx(Age, ChosenSurvMale)       ;  Lm<-Lm/Lm[Alphaf]    ;  Lm[Lm>1]<-1   # Vector of survival probabilities for male surviving at age Alphaf

				# Fertility
				Ff <- Fert(Age, c(ChosenFertFemale))                                # Vector of fecundity - Females
                Fm <- Fert(Age, c(ChosenFertMale))                                  # Vector of fecundity - Males 

				# Same mean number of children up to the 6th decimal
				sum(Lf*Ff) ; sum(Lm*Fm)

				
		# SCENARIOS FOR MATERNAL, GRANDMATERNAL AND PATERNAL INVESTMENTS - Calculation of the sigma(y1,y2,y3)
			# Because this calculation incorporates expensive computer time integration we found it best to store all the possible sigma(y1,y2,y3) as tables that can be further re-used.
			# The following code store this tables according to the different scenarios presented in the main text: no investment, M, M+GM, M + GM + F
 				
				
			# 1. Functions for Relative risk of death RR(y) according to child's age at MATERNAL, GRANDMATERNAL AND PATERNAL deaths (see SM 4.2)
          
				# RR Mother (3 parameters exponential:  a, b and c)
				RRm<-function(t, a, b, c){y<-a*exp(-b*t) + 1 + c ; return(y)}
			  
				# RR GrandMother - 1 parameter Uniform
				RRgm <- function(t, a){y<-a ; return(y)} 
			  
				# RR Father - 1 parameter Uniform
				RRf <- function(t, a){y<-a ; return(y)}	
				
			# 2. Function Sigma - Return the survival probability at age alpha as a function of y1 y2 and y3.
				# Inputs are y1 , y2, y3
				# Inputs Alphaf=first age at repro for female , par.Surv = Surv.parameters, vecInvestment=concatenated parameters of RRm,RRgm and RRF
				# Inputs are precision=timestep for integration
				
				Sigma <- function(y1, y2, y3, vAlphaf, par.Surv, vecInvestment, precision){                 
					# Equivalent for the child
					y2[y2<0]<-0
					y3[y3<0]<-0

					# Variable for integration
					vecage<-seq(0,vAlphaf,length.out=precision)
					timestep<-vAlphaf/precision
					  
					# The Risks distribution              
					vecRRm<-rep(1, length(vecage)) ; 
					vecRRm[vecage>=y1]<-RRm(vecage[vecage>=y1], vecInvestment[1],vecInvestment[2],vecInvestment[3])    
					vecRRm[vecage>=vAlphaf] <- 1 # (useless but to be safe)
					  
					vecRRgm<-rep(1, length(vecage)) ; 
					vecRRgm[vecage>=y2]<-RRgm(vecage[vecage>=y2], vecInvestment[4])   
					vecRRgm[vecage>=vAlphaf] <- 1 
					  
					vecRRf<-rep(1, length(vecage)) ; 
					vecRRf[vecage>=y3]<-RRf(vecage[vecage>=y3], vecInvestment[5]) 
					vecRRf[vecage>=vAlphaf] <- 1 
					  
					risque<-vecRRm*vecRRgm*vecRRf              

					# Calculation of the survival (attention: here child survival is independent of child sex)              
					theS<-exp(-(sum(vecRRm*vecRRgm*vecRRf*hx(vecage, par.Surv)*timestep)))   # ATTENTION - Age 15 not included.
					return(theS)
				} 
					
				# Sanity checks
				# Sx(15, ChosenSurvFemale)
				# vecage<-seq(0,15,length.out=50000)
				# timestep<-15/50000
				# exp(-(sum(hx(vecage, ChosenSurvFemale)*timestep)))  #Integration is working
				# # Without m mg f or for y1, y2, y3>alpha
				# invest0 <-c(0,  0,  0, 1, 1)
				# Salpha(5, -10, 10, 15, ChosenSurvFemale, invest0,    500)
				# Salpha(15, 15, 15, 15, ChosenSurvFemale, par.invest, 500) # this is working at the third decimal

			          ######################################################################
			# 3. Function table.sigma built the (very) large table gathering the sigma for all the possible 3-tuple y1, y2, y3
			
				table.Sigma<-function(vAlphaf, vBetaf, vAlpham, Omega, par.Surv, vecInvestment, precision){
					Y1<- 0:(Omega-vAlphaf)                                                  # possible child's age at mother's death
					Y2<- (-vBetaf):(Omega-2*vAlphaf)                                        # possible child's age at grandmother's death
					Y3<- 0:(Omega-vAlpham)                                                  # possible child's age at father's death
				
					tmp1<-rep(Y3, length(Y2))         # Rep de Y3 pour chaque Y2
					tmp2<-rep(Y2, each=length(Y3))    # Rep de chaque valeur Y2, Y3 fois
		  
					vY1<-rep(Y1, each=length(tmp2))
					vY2<-rep(tmp2, length(Y1))
					vY3<-rep(tmp1, length(Y1))        # Carefully checked => var. y3 , then var y2 then var y1

					Corespond.Salpha<-rep(NA, length(vY1))
						for (i in 1:length(vY1)){                                                          
						Corespond.Salpha[i] <- Sigma(vY1[i], vY2[i], vY3[i], vAlphaf, par.Surv, vecInvestment, precision)
						print(i)
						}
				  return(Corespond.Salpha)    
				}
	
			# 4. Saving Table for different scenarios [par.invest <- c(13.6691,  0.7086,  0.5678, 1.098, 1.763)] 
				# (i) No investment - (+ fertility of males equals that of females => juste change the size of the table because of change in alpham and betam)
				invest0 <-c(0,  0,  0, 1, 1)
				vecSigma <- table.Sigma(Alphaf, Betaf, Alpham, Omega,  ChosenSurvFemale, invest0, 500)
				length(vecSigma)
				#write(vecSigma, file = "C:/Users/Pavard/Documents/LaNouvelleSouille/P_Selection_Old_Age_NouvelleSouille2/Results_1/i.No.Invest_HG_fertf.csv", sep = ";")

				# (ii). Mother only - (+ fertility of males equals that of females => juste change the size of the table because of change in alpham and betam)
				m.only <- c(13.6691,  0.7086,  0.5678, 1, 1) 
				vecSigma<-table.Sigma(Alphaf, Betaf, Alphaf, Omega, ChosenSurvFemale, m.only, 500)
				length(vecSigma)
				#write(vecSigma, file = "C:/Users/Pavard/Documents/LaNouvelleSouille/P_Selection_Old_Age_NouvelleSouille2/Results_1/ii.m.only_HG_fertf.csv", sep = ";")

				# (iii). Mother + Grandmother
				m.gm.only <- c(13.6691,  0.7086,  0.5678, 1.098, 1) 
				vecSigma<-table.Sigma(Alphaf, Betaf, Alphaf, Omega, ChosenSurvFemale, m.gm.only, 500) # ATTENTION: precision of 500 is needed to be correct at the 5th decimal
				#write(vecSigma, file = "C:/Users/Pavard/Documents/LaNouvelleSouille/P_Selection_Old_Age_NouvelleSouille2/Results_1/v.m.gm_HG_fertf.csv", sep = ";")
				  
				# (iv). Mother + Grandmother + Father
				vecSigma<-table.Sigma(15, 49, 15, 99, ChosenSurvFemale, par.invest, 500) # ATTENTION: precision of 500 is needed to be correct at the 5th decimal
				#write(vecSigma, file = "C:/Users/Pavard/Documents/LaNouvelleSouille/P_Selection_Old_Age_NouvelleSouille2/Results_1/vi.m.gm.f_HG.csv", sep = ";")
				
				
				
				
				
				
		
		
		