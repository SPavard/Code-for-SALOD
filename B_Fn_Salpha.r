

	# Integrate child survival over all y1, y2 and y3 tuples to gives the Survival at age alpha of a children defined by its probability ditsribution of loosing its parents and grandparents
	# Inputs are the list of function P.y1.y2.y3 and the table denoted vecSalpha which is the table regrouping all the sigma.y1.y2.y3 for any y1, y2 and y3 values.
				
	S.Alpha <- function(List.P.y1.y2.y3, vecSigma){
            
		Y1<- List.P.y1.y2.y3$py1                                    # possible child's age at mother's death
		Y2<- List.P.y1.y2.y3$py2                                    # possible child's age at grandmother's death
		Y3<- List.P.y1.y2.y3$py3                                    # possible child's age at father's death

		tmp1<-rep(Y3, length(Y2))         # Rep of Y3 for any Y2
		tmp2<-rep(Y2, each=length(Y3))    # Rep of each Y2, Y3 times

		vY1<-rep(Y1, each=length(tmp2))
		vY2<-rep(tmp2, length(Y1))
		vY3<-rep(tmp1, length(Y1))        # Carefully checked => var. y3 , then var y2 then var y1

		result<-sum(vY1*vY2*vY3*vecSigma)
		return(result)
	}
